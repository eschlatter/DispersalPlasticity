######################## Main simulation function #######################
# Inputs:
#   params: list of biological and simulation parameters
#   hab_params: list of habitat-related parameters and objects; output of f_MakeHabitat
#   output_flag: "all" (npatch x ngroup x nsteps array of abundances) or "lite" (summary stats only for each timestep)
f_RunSim <- function(params, hab_params, keep=list("abund","p","kern","sp_struct"),output_flag="all",show_plot=FALSE,output_thin=1){
  starttime <- proc.time()
  numCores <- parallelly::availableCores()
  
  # load parameters
  list2env(x=params,envir=environment())
  
  # load habitat data structures (see f_MakeHabitat for details on what's in each object)
  list2env(x=hab_params,envir=environment()) 
  if(hab_type=="points") patch_locations$K <- 1 # sometimes these end up as 0 from map resolution issues
  # save original values, in case they're modified by disturbance
  K <- patch_locations$K 
  b <- patch_locations$b
  q <- patch_locations$q
  units(patch_dists) <- NULL
  units(nav_rad) <- NULL
  units(patch_angles) <- NULL
  
  ##### Data structures to describe population ######
  
  ## 1. group_index: all unique combinations of parameters alpha, theta, and p
  group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  ngroups <- nrow(group_index)
  
  ## 2. Population objects
  previous_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  new_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  if(output_flag=="all") Pij <- array(0, dim=c(npatch,ngroups,nsteps/output_thin))  # hold everything, if required
  
  # initialize previous_pop
  if(alpha_start==0) alpha_start = 1:length(v_alphas)
  if(theta_start==0) theta_start = 1:length(v_thetas)
  if(p_start==0) p_start = 1:length(v_p)
  start_grps <- which((group_index$alpha %in% alpha_start) & (group_index$theta %in% theta_start) & (group_index$p %in% p_start))
  start_probs <- rep(0,ngroups)
  start_probs[start_grps] <- 1
  starts <- lapply(1:npatch,function(i) as.vector(rmultinom(n=1,size=patch_locations$K[i],prob=start_probs))) 
  previous_pop <- do.call(rbind,starts)
  if(output_flag=="all") Pij[,,1] <- previous_pop
  
  ## 3. mutation_destinations: for each parameter groups, what parameter groups can a single mutation reach?
  ## dimensions 1 = ngroups, 2 = number of types of mutation events (including no mutation)
  # first, list the possible mutations to each parameter
  # (each index represents a mutation event; only one parameter can change per mutation event)
  alpha_adds=c(0,1,-1,0,0,0,0)
  theta_adds=c(0,0,0,1,-1,0,0)
  p_adds=c(0,0,0,0,0,1,-1)
  # then make the matrix
  mutation_destinations <- matrix(NA, nrow=ngroups, ncol=length(alpha_adds))
  for(mut_num in 1:length(alpha_adds)) {
    for(grp in 1:ngroups) {
      dest_grp <- which(group_index$alpha == group_index$alpha[grp] + alpha_adds[mut_num] & 
                          group_index$theta == group_index$theta[grp] + theta_adds[mut_num] &
                          group_index$p == group_index$p[grp] + p_adds[mut_num])
      mutation_destinations[grp, mut_num] <- ifelse(length(dest_grp)!=0, dest_grp, grp)
    }
  }
  
  ## 4. Initialize temporary data structures
  temp_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  to_patch <- numeric(npatch) # hold numbers of immigrants to each patch during dispersal
  all_conn_mats <- vector("list", ngroups) # list to hold connectivity matrices, which will be computed as they're needed (will have to be reset if K_i changes)
  
  ## 5. Initialize output data structures
  if(output_flag=="lite"){
    # df_out = data.frame(t_i=1:nsteps,abund=NA,p_mean=NA,p_var=NA,
    #                                  fund_kernmode_mean=NA,fund_kernmode_var=NA,
    #                                  fund_kernmean_mean=NA,fund_kernmean_var=NA,
    #                                  eff_kernmode_mean=NA,eff_kernmode_var=NA,
    #                                  eff_kernmean_mean=NA,eff_kernmean_var=NA)
    df_out <- data.frame(t_i=integer(),metric=character(),mean=numeric(),var=numeric())
  }
  
  # objects to hold param values in different ways for more efficient use later
  p_by_group <- v_p[group_index$p] # value, not index
  alpha_by_group <- group_index$alpha # index, not value
  theta_by_group <- group_index$theta # index, not value
  Pij_alpha <- matrix(alpha_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_theta <- matrix(theta_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_alpha_val <- v_alphas[Pij_alpha]
  Pij_theta_val <- v_thetas[Pij_theta]
  Pij_p <- matrix(p_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # value, not index
  
  # # generate matrix of weights (inverse distance) to use in calculating Moran's I
  # # is this an okay distance function to use? I know it matters when doing the statistical test,
  # # but maybe as long as we use the same metric for all timesteps and simulations it's OK for comparing among them.
  # # I tried a few options (1/d, 1/d^2, 1/(1+d^2)) and the dynamics were quite similar.
  # moran_weights <- 1/(patch_dists) 
  # diag(moran_weights) <- 0
  
  ################ Simulate ###################
  interval_starttime <- proc.time()
  for(t_i in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    ################## Disturbance ##################
    
    if(t_i %% 10 == 0){ # disturbances last for 10 timesteps
      # first, reset all K's from any disturbance that occurred 10 timesteps previously
      if(prod(patch_locations$K == as.vector(K))==0){
        patch_locations$K = as.vector(K)
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K
        print(paste("t_i =",t_i,"reset connectivity matrices"))
      }
      
      # then make a new disturbance (maybe)
      if(rbinom(n=1,size=1,p=disturb_prob)==1){
        print(paste("t_i =",t_i,"disturbance"))
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K
        
        disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
        disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
        patch_locations$K[disturb_patches] <- 0
      }
    }
    
    ################## Reproduction and Dispersal and Mutation ##################
    pop_by_group <- colSums(previous_pop)
    
    for(g in which(pop_by_group>0)){
      # get parameter values for that parameter group
      v <- group_index[g,]
      p_penalty <- abs(v_p[v$p])*0
      
      # get population of each patch for that parameter group
      patch_pops <- previous_pop[,g]
      
      # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
      # (and accounting for the patch population x per capita output b_i from each patch)
      if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
                                                                                patch_dists=patch_dists, patch_angles=patch_angles, 
                                                                                overlap_discount=overlap_discount, 
                                                                                v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas, 
                                                                                nav_rad=nav_rad, numCores=numCores) # calculate this matrix, if it hasn't been used before
      conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
      to_patch <- (1-p_penalty)*patch_locations$b*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
      
      # Divide up to_patch among parameter groups that are the result of mutation
      temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
      
      for(mut_group in mutation_destinations[g,-1]){ # for each of the possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
        temp_pop[,mut_group] <- (mu/6)*to_patch+temp_pop[,mut_group]
      }
    }
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    # patch_abunds <- rowSums(temp_pop)
    # for(i_patch in which(patch_abunds>0)){
    #   survivors=rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
    #                       size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
    #                       prob = temp_pop[i_patch,]) # probability of each cell being chosen depends on its current abundance)
    #   new_pop[i_patch,] <- survivors
    # }
    
    patch_abunds <- rowSums(temp_pop)
    comp_results <- lapply(1:npatch,function(i) f_Competition(i_patch=i,patch_abunds=patch_abunds,patch_locations=patch_locations,
                                                              temp_pop=temp_pop,ngroups=ngroups))
    #comp_results <- mclapply(1:npatch,function(i) f_Competition(i,patch_abunds,patch_locations,temp_pop),mc.cores = numCores)
    new_pop <- do.call(rbind,comp_results)
    
    previous_pop <- new_pop
    if(output_flag=="all" & (t_i %% output_thin == 0)) {
      Pij[,,t_i/output_thin] <- new_pop
    }

    ################## Output ##################
    if(t_i %% max(1,round(nsteps/10)) == 0){
      # status updates to console every so often
      print(t_i)
      print(proc.time()-interval_starttime)
      interval_starttime <- proc.time()
      
      # plot, if requested
      if(show_plot==TRUE){
        if(hab_type=="points"){
          g_map_abund <- ggplot(reef_sf)+
            geom_sf()+
            geom_sf(data=sfc_patches,aes(color=factor(rowSums(new_pop))))+
            labs(color="Abundance",title=paste0("t_i=",t_i))+
            scale_color_manual(values=c("blue","red"),breaks=c(1,0))+
            theme_minimal()+
            annotation_scale()
          print(g_map_abund)
        }
        if(hab_type=="grid"){
          g_map_abund <- ggplot(patch_locations,aes(x=x,y=y))+
            geom_tile(aes(fill=rowSums(new_pop)))+
            labs(title=paste0("t_i=",t_i),fill="abund")+
            scale_color_continuous(limits=c(0,max(K)))
          print(g_map_abund)
        }
      }
    }
    
    ## these are for output_flag==lite. Needs work.
    if(output_flag=="lite"){
      if(t_i %% output_thin == 0){
        if(sum(previous_pop)>0){
          df_i <- data.frame(t_i=rep(t_i,4),metric=NA,mean=NA,var=NA)
          
          ###### total abundance
          #df_out$abund[t_i] <- sum(previous_pop)
          df_i[1,] <- data.frame(t_i=t_i,metric="abund",mean=sum(previous_pop),var=NA)
          
          ###### p
          # mean value of p at each timestep
          p_mean <- sum(p_by_group*colSums(previous_pop))/sum(previous_pop)
          p_var <- sum(colSums(previous_pop)*(p_by_group-p_mean)^2)/sum(previous_pop)
          df_i[2,] <- data.frame(t_i=t_i,metric="p",mean=p_mean,var=p_var)
          
          ###### kernel properties
          Pij_b <- matrix(patch_locations$b,byrow=FALSE,nrow=npatch,ncol=ngroups)
          Pij_eff <- f_plasticityb(Pij_b,Pij_p,Pij_alpha,Pij_theta,length(v_alphas),length(v_thetas)) # indices, not values
          
          ## fundamental kernel properties (i.e., based on inherited alpha and theta, but not plasticity)
          # 1. kernel mode
          fund_mode_each <- ifelse(Pij_alpha_val<1,0,(Pij_alpha_val-1)*Pij_theta_val)
          # mean and variance (over the population) of the fundamental kernel mode at each timestep
          fund_mode_mean <- sum(fund_mode_each*previous_pop)/sum(previous_pop)
          fund_kernmode_var <- sum(previous_pop*(fund_mode_each-fund_mode_mean)^2)/sum(previous_pop)
          df_i[3,] <- data.frame(t_i=t_i,metric="fund_kernmode",mean=fund_mode_mean,var=fund_kernmode_var)
          # 2. kernel mean
          fund_mean_each <- Pij_alpha_val*Pij_theta_val
          # mean and variance (over the population) of the fundamental kernel mean at each timestep
          fund_mean_mean <- sum(fund_mean_each*previous_pop)/sum(previous_pop)
          fund_kernmean_var <- sum(previous_pop*(fund_mean_each-fund_mean_mean)^2)/sum(previous_pop)
          df_i[4,] <- data.frame(t_i=t_i,metric="fund_kernmean",mean=fund_mean_mean,var=fund_kernmean_var)
          
          ## effective kernel properties (i.e., accounting for plasticity)
          #1. kernel mode -- something is weird here!!
          eff_mode_each <- ifelse(v_alphas[Pij_eff$alpha_plastic]<1,0,
                                  (v_alphas[Pij_eff$alpha_plastic]-1)*v_alphas[Pij_eff$theta_plastic])
          # mean and variance (over the population) of the effective kernel mode at each timestep
          eff_mode_mean <- sum(eff_mode_each*previous_pop)/sum(previous_pop)
          eff_mode_var <- sum(previous_pop*(eff_mode_each-eff_mode_mean)^2)/sum(previous_pop)
          df_i[5,] <- data.frame(t_i=t_i,metric="eff_kernmode",mean=eff_mode_mean,var=eff_mode_var)
          # 2. kernel mean
          eff_mean_each <- v_alphas[Pij_eff$alpha_plastic]*v_thetas[Pij_eff$theta_plastic]
          # mean and variance (over the population) of the effective kernel mean at each timestep
          eff_mean_mean <- sum(eff_mean_each*previous_pop)/sum(previous_pop)
          eff_mean_var <- sum(previous_pop*(eff_mean_each-eff_mean_mean)^2)/sum(previous_pop)
          df_i[6,] <- data.frame(t_i=t_i,metric="eff_kernmean",mean=eff_mean_mean,var=eff_mean_var)
          
          ###### spatial structure
          # # some sort of measure of the spatial autocorrelation of fundamental kernel mean
          # fund_mean_by_patch <- as.vector(rowSums(fund_mean_each*previous_pop))
          # output_list$df_fund$kernmean_moran[t_i] <- Moran.I(fund_mean_by_patch,weight=moran_weights)$observed
          # 
          # # some sort of measure of the spatial autocorrelation of effective kernel mean
          # eff_mean_by_patch <- as.vector(rowSums(eff_mean_each*previous_pop))
          # output_list$df_eff$kernmean_moran[t_i] <- Moran.I(eff_mean_by_patch,weight=moran_weights)$observed
          
          df_out <- rbind(df_out,df_i)
        }
      }
    }
    
  } #t_i
  
  time_run <- proc.time()-starttime
  output_list <- list(params=params,hab_params=hab_params,group_index=group_index,time_run=time_run)
  if(output_flag=="all") output_list <- c(output_list,list(Pij=Pij))
  if(output_flag=="lite") output_list <- c(output_list,list(df_out=df_out))
  return(output_list)
}

######################## Main simulation function #######################
# Inputs:
#   params: list of biological and simulation parameters
#   hab_params: list of habitat-related parameters and objects; output of f_MakeHabitat
#   output_flag: "all" (npatch x ngroup x nsteps array of abundances) or "lite" (summary stats only for each timestep)
f_RunSimNew <- function(params, hab_params, keep=list("abund","p","kern","sp_struct"),
                        output_flag="all",show_plot=FALSE,output_thin=1,output_file=NA){
  # run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  # if(is.na(run_i)) run_i <- 1
  # # Create a new folder in the MSI scratch directory
  # temp_dir <- as.character(paste0(round(as.numeric(Sys.time())),run_i,sample(1:1000,1)))
  # temp_path <- file.path("/scratch.global","schla103",temp_dir)
  # dir.create(temp_path)
  
  starttime <- proc.time()
  numCores <- parallelly::availableCores()
  
  # load parameters
  list2env(x=params,envir=environment())
  
  # load habitat data structures (see f_MakeHabitat for details on what's in each object)
  list2env(x=hab_params,envir=environment()) 
  if(hab_type=="points") patch_locations$K <- 1 # sometimes these end up as 0 from map resolution issues
  # save original values, in case they're modified by disturbance
  K <- patch_locations$K 
  b <- patch_locations$b
  q <- patch_locations$q
  units(patch_dists) <- NULL
  units(nav_rad) <- NULL
  units(patch_angles) <- NULL
  
  ##### Data structures to describe population ######
  
  ## 1. group_index: all unique combinations of parameters alpha, theta, and p
  group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  ngroups <- nrow(group_index)
  
  ## 2. Population objects
  previous_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  new_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  
  # initialize previous_pop
  if(alpha_start==0) alpha_start = 1:length(v_alphas)
  if(theta_start==0) theta_start = 1:length(v_thetas)
  if(p_start==0) p_start = 1:length(v_p)
  start_grps <- which((group_index$alpha %in% alpha_start) & (group_index$theta %in% theta_start) & (group_index$p %in% p_start))
  start_probs <- rep(0,ngroups)
  start_probs[start_grps] <- 1
  starts <- lapply(1:npatch,function(i) as.vector(rmultinom(n=1,size=patch_locations$K[i],prob=start_probs))) 
  previous_pop <- do.call(rbind,starts)
  
  ## 3. mutation_destinations: for each parameter groups, what parameter groups can a single mutation reach?
  ## dimensions 1 = ngroups, 2 = number of types of mutation events (including no mutation)
  # first, list the possible mutations to each parameter
  # (each index represents a mutation event; only one parameter can change per mutation event)
  alpha_adds=c(0,1,-1,0,0,0,0)
  theta_adds=c(0,0,0,1,-1,0,0)
  p_adds=c(0,0,0,0,0,1,-1)
  # then make the matrix
  mutation_destinations <- matrix(NA, nrow=ngroups, ncol=length(alpha_adds))
  for(mut_num in 1:length(alpha_adds)) {
    for(grp in 1:ngroups) {
      dest_grp <- which(group_index$alpha == group_index$alpha[grp] + alpha_adds[mut_num] & 
                          group_index$theta == group_index$theta[grp] + theta_adds[mut_num] &
                          group_index$p == group_index$p[grp] + p_adds[mut_num])
      mutation_destinations[grp, mut_num] <- ifelse(length(dest_grp)!=0, dest_grp, grp)
    }
  }
  
  ## 4. Initialize temporary data structures
  temp_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  to_patch <- numeric(npatch) # hold numbers of immigrants to each patch during dispersal
  all_conn_mats <- vector("list", ngroups) # list to hold connectivity matrices, which will be computed as they're needed (will have to be reset if K_i changes)
  
  ## 5. Initialize output data structures
  if(output_flag=="lite"){
    # df_out = data.frame(t_i=1:nsteps,abund=NA,p_mean=NA,p_var=NA,
    #                                  fund_kernmode_mean=NA,fund_kernmode_var=NA,
    #                                  fund_kernmean_mean=NA,fund_kernmean_var=NA,
    #                                  eff_kernmode_mean=NA,eff_kernmode_var=NA,
    #                                  eff_kernmean_mean=NA,eff_kernmean_var=NA)
    df_out <- data.frame(t_i=integer(),metric=character(),mean=numeric(),var=numeric())
  }
  
  # objects to hold param values in different ways for more efficient use later
  p_by_group <- v_p[group_index$p] # value, not index
  alpha_by_group <- group_index$alpha # index, not value
  theta_by_group <- group_index$theta # index, not value
  Pij_alpha <- matrix(alpha_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_theta <- matrix(theta_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_alpha_val <- v_alphas[Pij_alpha]
  Pij_theta_val <- v_thetas[Pij_theta]
  Pij_p <- matrix(p_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # value, not index
  
  # # generate matrix of weights (inverse distance) to use in calculating Moran's I
  # # is this an okay distance function to use? I know it matters when doing the statistical test,
  # # but maybe as long as we use the same metric for all timesteps and simulations it's OK for comparing among them.
  # # I tried a few options (1/d, 1/d^2, 1/(1+d^2)) and the dynamics were quite similar.
  # moran_weights <- 1/(patch_dists) 
  # diag(moran_weights) <- 0
  
  ## 6. Save input information
  metadata_list <- list(params=params,group_index=group_index,patch_locations=patch_locations,hab_file=hab_params$hab_file)
  save(metadata_list,file=paste0(output_file,"_metadata.RData"))
  
  ################ Simulate ###################
  interval_starttime <- proc.time()
  for(t_i in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    ################## Disturbance ##################
    
    if(t_i %% 10 == 0){ # disturbances last for 10 timesteps
      # first, reset all K's from any disturbance that occurred 10 timesteps previously
      if(prod(patch_locations$K == as.vector(K))==0){
        patch_locations$K = as.vector(K)
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K
        print(paste("t_i =",t_i,"reset connectivity matrices"))
      }
      
      # then make a new disturbance (maybe)
      if(rbinom(n=1,size=1,p=disturb_prob)==1){
        print(paste("t_i =",t_i,"disturbance"))
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K
        
        disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
        disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
        patch_locations$K[disturb_patches] <- 0
      }
    }
    
    ################## Reproduction and Dispersal and Mutation ##################
    pop_by_group <- colSums(previous_pop)
    
    for(g in which(pop_by_group>0)){
      # get parameter values for that parameter group
      v <- group_index[g,]
      p_penalty <- abs(v_p[v$p])*0
      
      # get population of each patch for that parameter group
      patch_pops <- previous_pop[,g]
      
      # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
      # (and accounting for the patch population x per capita output b_i from each patch)
      # if(!file.exists(paste0(temp_path,"/grp_",g))){
      #   # calculate this matrix, if it hasn't been used before
      #   conn_mat <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
      #                                   patch_dists=patch_dists, patch_angles=patch_angles, 
      #                                   overlap_discount=overlap_discount, 
      #                                   v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas, 
      #                                   nav_rad=nav_rad, numCores=numCores)
      #   # then store it
      #   write_fst(as.data.frame(conn_mat),paste0(temp_path,"/grp_",g))
      #   # otherwise, import it
      # } else conn_mat <- as.matrix(read_fst(paste0(temp_path,"/grp_",g)))
      
      if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
                                                                                patch_dists=patch_dists, patch_angles=patch_angles,
                                                                                overlap_discount=overlap_discount,
                                                                                v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas,
                                                                                nav_rad=nav_rad, numCores=numCores) # calculate this matrix, if it hasn't been used before
      conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
      to_patch <- (1-p_penalty)*patch_locations$b*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
      
      # Divide up to_patch among parameter groups that are the result of mutation
      temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
      
      for(mut_group in mutation_destinations[g,-1]){ # for each of the possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
        temp_pop[,mut_group] <- (mu/6)*to_patch+temp_pop[,mut_group]
      }
    }
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    # patch_abunds <- rowSums(temp_pop)
    # for(i_patch in which(patch_abunds>0)){
    #   survivors=rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
    #                       size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
    #                       prob = temp_pop[i_patch,]) # probability of each cell being chosen depends on its current abundance)
    #   new_pop[i_patch,] <- survivors
    # }
    
    patch_abunds <- rowSums(temp_pop)
    comp_results <- lapply(1:npatch,function(i) f_Competition(i_patch=i,patch_abunds=patch_abunds,patch_locations=patch_locations,
                                                              temp_pop=temp_pop,ngroups=ngroups))
    #comp_results <- mclapply(1:npatch,function(i) f_Competition(i,patch_abunds,patch_locations,temp_pop),mc.cores = numCores)
    new_pop <- do.call(rbind,comp_results)
    
    previous_pop <- new_pop
    
    ################## Output ##################
    if(t_i %% max(1,round(nsteps/1000)) == 0){
      # status updates to console every so often
      print(t_i)
      print(proc.time()-interval_starttime)
      interval_starttime <- proc.time()
      
      # plot, if requested
      if(show_plot==TRUE){
        if(hab_type=="points"){
          g_map_abund <- ggplot(reef_sf)+
            geom_sf()+
            geom_sf(data=sfc_patches,aes(color=factor(rowSums(new_pop))))+
            labs(color="Abundance",title=paste0("t_i=",t_i))+
            scale_color_manual(values=c("blue","red"),breaks=c(1,0))+
            theme_minimal()+
            annotation_scale()
          print(g_map_abund)
        }
        if(hab_type=="grid"){
          g_map_abund <- ggplot(patch_locations,aes(x=x,y=y))+
            geom_tile(aes(fill=rowSums(new_pop)))+
            labs(title=paste0("t_i=",t_i),fill="abund")+
            scale_color_continuous(limits=c(0,max(K)))
          print(g_map_abund)
        }
      }
    }
    
    if((output_flag=="all") & (t_i %% output_thin == 0)){
      export_mat <- cbind(rep(t_i,times=npatch),previous_pop)
      fwrite(export_mat,file=paste0(output_file,"_raw.csv"),append=TRUE)
    }
  } #t_i
  
  time_run <- proc.time()-starttime
  return(time_run)
}

######################## Main simulation function #######################
# Inputs:
#   params: list of biological and simulation parameters
#   hab_params: list of habitat-related parameters and objects; output of f_MakeHabitat
#   output_flag: "all" (npatch x ngroup x nsteps array of abundances) or "lite" (summary stats only for each timestep)
f_RunSimNewScratch <- function(params, hab_params, keep=list("abund","p","kern","sp_struct"),
                        output_flag="all",show_plot=FALSE,output_thin=1,output_file=NA){
  run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if(is.na(run_i)) run_i <- 1
  # Create a new folder in the MSI scratch directory
  temp_dir <- as.character(paste0(round(as.numeric(Sys.time())),run_i,sample(1:1000,1)))
  temp_path <- file.path("/scratch.local","schla103",temp_dir)
  dir.create(temp_path)
  
  starttime <- proc.time()
  numCores <- parallelly::availableCores()
  
  # load parameters
  list2env(x=params,envir=environment())
  
  # load habitat data structures (see f_MakeHabitat for details on what's in each object)
  list2env(x=hab_params,envir=environment()) 
  if(hab_type=="points") patch_locations$K <- 1 # sometimes these end up as 0 from map resolution issues
  # save original values, in case they're modified by disturbance
  K <- patch_locations$K 
  b <- patch_locations$b
  q <- patch_locations$q
  units(patch_dists) <- NULL
  units(nav_rad) <- NULL
  units(patch_angles) <- NULL
  
  ##### Data structures to describe population ######
  
  ## 1. group_index: all unique combinations of parameters alpha, theta, and p
  group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  ngroups <- nrow(group_index)
  
  ## 2. Population objects
  previous_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  new_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  
  # initialize previous_pop
  if(alpha_start==0) alpha_start = 1:length(v_alphas)
  if(theta_start==0) theta_start = 1:length(v_thetas)
  if(p_start==0) p_start = 1:length(v_p)
  start_grps <- which((group_index$alpha %in% alpha_start) & (group_index$theta %in% theta_start) & (group_index$p %in% p_start))
  start_probs <- rep(0,ngroups)
  start_probs[start_grps] <- 1
  starts <- lapply(1:npatch,function(i) as.vector(rmultinom(n=1,size=patch_locations$K[i],prob=start_probs))) 
  previous_pop <- do.call(rbind,starts)
  
  ## 3. mutation_destinations: for each parameter groups, what parameter groups can a single mutation reach?
  ## dimensions 1 = ngroups, 2 = number of types of mutation events (including no mutation)
  # first, list the possible mutations to each parameter
  # (each index represents a mutation event; only one parameter can change per mutation event)
  alpha_adds=c(0,1,-1,0,0,0,0)
  theta_adds=c(0,0,0,1,-1,0,0)
  p_adds=c(0,0,0,0,0,1,-1)
  # then make the matrix
  mutation_destinations <- matrix(NA, nrow=ngroups, ncol=length(alpha_adds))
  for(mut_num in 1:length(alpha_adds)) {
    for(grp in 1:ngroups) {
      dest_grp <- which(group_index$alpha == group_index$alpha[grp] + alpha_adds[mut_num] & 
                          group_index$theta == group_index$theta[grp] + theta_adds[mut_num] &
                          group_index$p == group_index$p[grp] + p_adds[mut_num])
      mutation_destinations[grp, mut_num] <- ifelse(length(dest_grp)!=0, dest_grp, grp)
    }
  }
  
  ## 4. Initialize temporary data structures
  temp_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  to_patch <- numeric(npatch) # hold numbers of immigrants to each patch during dispersal
  # all_conn_mats <- vector("list", ngroups) # list to hold connectivity matrices, which will be computed as they're needed (will have to be reset if K_i changes)
  
  ## 5. Initialize output data structures
  if(output_flag=="lite"){
    # df_out = data.frame(t_i=1:nsteps,abund=NA,p_mean=NA,p_var=NA,
    #                                  fund_kernmode_mean=NA,fund_kernmode_var=NA,
    #                                  fund_kernmean_mean=NA,fund_kernmean_var=NA,
    #                                  eff_kernmode_mean=NA,eff_kernmode_var=NA,
    #                                  eff_kernmean_mean=NA,eff_kernmean_var=NA)
    df_out <- data.frame(t_i=integer(),metric=character(),mean=numeric(),var=numeric())
  }
  
  # objects to hold param values in different ways for more efficient use later
  p_by_group <- v_p[group_index$p] # value, not index
  alpha_by_group <- group_index$alpha # index, not value
  theta_by_group <- group_index$theta # index, not value
  Pij_alpha <- matrix(alpha_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_theta <- matrix(theta_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_alpha_val <- v_alphas[Pij_alpha]
  Pij_theta_val <- v_thetas[Pij_theta]
  Pij_p <- matrix(p_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # value, not index
  
  # # generate matrix of weights (inverse distance) to use in calculating Moran's I
  # # is this an okay distance function to use? I know it matters when doing the statistical test,
  # # but maybe as long as we use the same metric for all timesteps and simulations it's OK for comparing among them.
  # # I tried a few options (1/d, 1/d^2, 1/(1+d^2)) and the dynamics were quite similar.
  # moran_weights <- 1/(patch_dists) 
  # diag(moran_weights) <- 0
  
  ## 6. Save input information
  metadata_list <- list(params=params,group_index=group_index,patch_locations=patch_locations,hab_file=hab_params$hab_file)
  save(metadata_list,file=paste0(output_file,"_metadata.RData"))
  
  ################ Simulate ###################
  interval_starttime <- proc.time()
  for(t_i in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    ################## Disturbance ##################
    
    if(t_i %% 10 == 0){ # disturbances last for 10 timesteps
      # first, reset all K's from any disturbance that occurred 10 timesteps previously
      if(prod(patch_locations$K == as.vector(K))==0){
        patch_locations$K = as.vector(K)
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K
        print(paste("t_i =",t_i,"reset connectivity matrices"))
      }
      
      # then make a new disturbance (maybe)
      if(rbinom(n=1,size=1,p=disturb_prob)==1){
        print(paste("t_i =",t_i,"disturbance"))
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K
        
        disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
        disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
        patch_locations$K[disturb_patches] <- 0
      }
    }
    
    ################## Reproduction and Dispersal and Mutation ##################
    pop_by_group <- colSums(previous_pop)
    
    for(g in which(pop_by_group>0)){
      # get parameter values for that parameter group
      v <- group_index[g,]
      p_penalty <- abs(v_p[v$p])*0
      
      # get population of each patch for that parameter group
      patch_pops <- previous_pop[,g]
      
      # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
      # (and accounting for the patch population x per capita output b_i from each patch)
      if(!file.exists(paste0(temp_path,"/grp_",g))){
        # calculate this matrix, if it hasn't been used before
        conn_mat <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
                                        patch_dists=patch_dists, patch_angles=patch_angles,
                                        overlap_discount=overlap_discount,
                                        v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas,
                                        nav_rad=nav_rad, numCores=numCores)
        # then store it
        write_fst(as.data.frame(conn_mat),paste0(temp_path,"/grp_",g),compress=50)
        # otherwise, import it
      } else conn_mat <- as.matrix(read_fst(paste0(temp_path,"/grp_",g)))
      
      # if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
      #                                                                           patch_dists=patch_dists, patch_angles=patch_angles,
      #                                                                           overlap_discount=overlap_discount,
      #                                                                           v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas,
      #                                                                           nav_rad=nav_rad, numCores=numCores) # calculate this matrix, if it hasn't been used before
      # conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
      to_patch <- (1-p_penalty)*patch_locations$b*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
      
      # Divide up to_patch among parameter groups that are the result of mutation
      temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
      
      for(mut_group in mutation_destinations[g,-1]){ # for each of the possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
        temp_pop[,mut_group] <- (mu/6)*to_patch+temp_pop[,mut_group]
      }
    }
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    # patch_abunds <- rowSums(temp_pop)
    # for(i_patch in which(patch_abunds>0)){
    #   survivors=rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
    #                       size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
    #                       prob = temp_pop[i_patch,]) # probability of each cell being chosen depends on its current abundance)
    #   new_pop[i_patch,] <- survivors
    # }
    
    patch_abunds <- rowSums(temp_pop)
    comp_results <- lapply(1:npatch,function(i) f_Competition(i_patch=i,patch_abunds=patch_abunds,patch_locations=patch_locations,
                                                              temp_pop=temp_pop,ngroups=ngroups))
    #comp_results <- mclapply(1:npatch,function(i) f_Competition(i,patch_abunds,patch_locations,temp_pop),mc.cores = numCores)
    new_pop <- do.call(rbind,comp_results)
    
    previous_pop <- new_pop
    
    ################## Output ##################
    if(t_i %% max(1,round(nsteps/1000)) == 0){
      # status updates to console every so often
      print(t_i)
      print(proc.time()-interval_starttime)
      interval_starttime <- proc.time()
      
      # plot, if requested
      if(show_plot==TRUE){
        if(hab_type=="points"){
          g_map_abund <- ggplot(reef_sf)+
            geom_sf()+
            geom_sf(data=sfc_patches,aes(color=factor(rowSums(new_pop))))+
            labs(color="Abundance",title=paste0("t_i=",t_i))+
            scale_color_manual(values=c("blue","red"),breaks=c(1,0))+
            theme_minimal()+
            annotation_scale()
          print(g_map_abund)
        }
        if(hab_type=="grid"){
          g_map_abund <- ggplot(patch_locations,aes(x=x,y=y))+
            geom_tile(aes(fill=rowSums(new_pop)))+
            labs(title=paste0("t_i=",t_i),fill="abund")+
            scale_color_continuous(limits=c(0,max(K)))
          print(g_map_abund)
        }
      }
    }
    
    if((output_flag=="all") & (t_i %% output_thin == 0)){
      export_mat <- cbind(rep(t_i,times=npatch),previous_pop)
      fwrite(export_mat,file=paste0(output_file,"_raw.csv"),append=TRUE)
    }
  } #t_i
  
  time_run <- proc.time()-starttime
  return(time_run)
}

######################## function to do competition step in parallel #######################
f_Competition <- function(i_patch,patch_abunds,patch_locations,temp_pop,ngroups){
  # decide the maximum number of larvae that can settle in each patch, based on how many arrived (patch_abunds[i_patch])
  if(patch_abunds[i_patch]>0 & patch_abunds[i_patch]<1){
    n_setts <- rbinom(1,1,prob=patch_abunds[i_patch])
  }
  else if(patch_abunds[i_patch]>=1){
    n_setts <- round(patch_abunds[i_patch])
  }
  else n_setts=0
  
  # do the competition
  if(n_setts>0){ # if the patch isn't empty
    survivors=t(rmultinom(n=1,
                          size=min(n_setts,patch_locations$K[i_patch]), # choose groups for min(n_setts, K) settlers
                          prob = temp_pop[i_patch,])) # probability of each group being chosen depends on its current abundance
  }
  else survivors <- vector(mode='integer',length=ngroups)
  return(survivors)
}

############# function to do reproduction, dispersal, and mutation steps in parallel #############
# (not using this; it actually slows things down)
f_ReprodDispMut <- function(g,group_index,previous_step,mutation_destinations,all_conn_mats,patch_locations,patch_dists,patch_angles,v_p,v_alphas,v_thetas,mu){
  temp_pop_g <- matrix(0,nrow=nrow(patch_locations),ncol=nrow(group_index)) # hold intermediate population values from this group
  
  # get parameter values for that parameter group
  v <- group_index[g,] 
  
  # get population of each patch for that parameter group
  patch_pops <- previous_step[,g]
  if(sum(patch_pops)>0){
    
    # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
    # (and accounting for the patch population x per capita output b_i from each patch)
    if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, overlap_discount, v_p, v_alphas, v_thetas) # calculate this matrix, if it hasn't been used before
    conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
    to_patch <- patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
    
    # Divide up to_patch among parameter groups that are the result of mutation
    temp_pop_g[,g] <- (1-mu)*to_patch+temp_pop_g[,g]
    
    for(mut_group in mutation_destinations[g,-1]){ # for each of the 4 possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
      temp_pop_g[,mut_group] <- (mu/4)*to_patch+temp_pop_g[,mut_group]
    }
  }
  
  return(temp_pop_g)
}

############## process Pij into a dataframe #############
f_ProcessPij <- function(Pij,patch_locations,group_index,output_thin){
  Pij <- Pij[,,seq(from=output_thin,to=dim(Pij)[3],by=output_thin)]
  npatch <- nrow(patch_locations)
  ngroups <- nrow(group_index)
  nsteps <- dim(Pij)[3]
  
  # melt into a datatable with columns patch, timestep, alpha, theta, p, popsize
  # add columns for param values at each time/patch/alpha/theta/p combo, scaled by the population size that has that combo
  sim_melt <- as.table(Pij)
  dimnames(sim_melt) <- list(1:npatch,1:ngroups,output_thin*(1:nsteps))
  sim_melt <- as.data.table(sim_melt)
  setnames(sim_melt, c("patch","group","t_i","popsize"))
  
  sim_melt[,`:=`(patch = as.numeric(patch),
                 group = as.numeric(group),
                 t_i = output_thin*as.numeric(t_i))]
  sim_melt[,`:=`(alpha = group_index$alpha[group],
                 theta = group_index$theta[group],
                 p = group_index$p[group])]
  sim_melt[,`:=`(alpha_value = v_alphas[alpha],
                 theta_value = v_thetas[theta],
                 p_value = v_p[p])]
  
  # mean param values at each timepoint
  # first take the sum across all cells of the param*popsize (numerator of the mean)
  by_t <- sim_melt[ ,.(alpha=sum(alpha_value*popsize),
                       theta=sum(theta_value*popsize),
                       p=sum(p_value*popsize),
                       popsize=sum(popsize)),
                    by=t_i]
  # then divide by total popsize (denominator of the mean)
  by_t[,`:=`(alpha = alpha/popsize,
             theta = theta/popsize,
             p = p/popsize)]
  
  return(list(sim_melt=sim_melt,
              by_t=by_t))
}

f_ProcessPopmat <- function(popmat,patch_locations,group_index,params){
  npatch <- nrow(patch_locations)
  v_alphas <- params$v_alphas
  v_thetas <- params$v_thetas
  v_p <- params$v_p
  
  # get the data into long format
  popmat <- cbind(rep(1:npatch,length.out=nrow(popmat)),popmat)
  colnames(popmat) <- c("patch","t_i",as.character(1:(ncol(popmat)-2)))
  sim_melt <- melt.data.table(popmat,measure.vars=3:ncol(popmat),variable.name="group",value.name="popsize")
  sim_melt[,`:=`(group = as.numeric(group))]
  
  # add columns for alpha, theta, and p
  sim_melt[,`:=`(alpha = group_index$alpha[group],
                 theta = group_index$theta[group],
                 p = group_index$p[group])]
  sim_melt[,`:=`(alpha_value = v_alphas[alpha],
                 theta_value = v_thetas[theta],
                 p_value = v_p[p])]
  
  # mean param values at each timepoint
  # first take the sum across all cells of the param*popsize (numerator of the mean)
  by_t <- sim_melt[ ,.(alpha=sum(alpha_value*popsize),
                       theta=sum(theta_value*popsize),
                       p=sum(p_value*popsize),
                       popsize=sum(popsize)),
                    by=t_i]
  # then divide by total popsize (denominator of the mean)
  by_t[,`:=`(alpha = alpha/popsize,
             theta = theta/popsize,
             p = p/popsize)]
  
  return(list(sim_melt=sim_melt,
              by_t=by_t))
}

# take the sim_melt dataframe, and get useful output
f_ProcessSimMelt <- function(sim_melt,patch_locations,group_index){
  sim_melt <- sim_melt[popsize>0,] # first, get rid of the rows with zero abundance
  
  # prepare other data
  n_alpha <- length(unique(group_index$alpha))
  n_theta <- length(unique(group_index$theta))
  patch_locations <- as.data.table(patch_locations)
  
  # merge b and K columns into sim_melt
  sim_melt[patch_locations, on=list(patch=id), b := i.b]
  sim_melt[patch_locations, on=list(patch=id), K := i.K]
  
  # create columns for alpha_plastic and theta_plastic (values, not indices)
  sim_melt[, `:=`(alpha_plastic = v_alphas[f_plasticityb(b,p_value,alpha,theta,n_alpha,n_theta)$alpha_plastic],
                  theta_plastic = v_thetas[f_plasticityb(b,p_value,alpha,theta,n_alpha,n_theta)$theta_plastic])
  ]
  
  # create columns for kernel mean and mode
  sim_melt[, `:=`(fund_kernmean = alpha_value*theta_value,
                  fund_kernmode = ifelse(alpha_value<1,0,(alpha_value-1)*theta_value),
                  realized_kernmean = alpha_plastic*theta_plastic,
                  realized_kernmode = ifelse(alpha_plastic<1,0,(alpha_plastic-1)*theta_plastic))]
  
  # make a by_patch, that collects all groups in each patch at each timestep
  
  
  # make a by_t, that collects all values at each timestep
  
}


############## process raw output of lite simulation function #############

f_GetOutputStats <- function(output_df,nsteps){
  equil_pt <- 0.25*nsteps # I need a better method than this.
  dat <- output_df[equil_pt:nsteps,]
  
  # a value of p (and uncertainty)
  p_mean <- mean(dat$v_pmeans)
  p_var_among_t <- var(dat$v_pmeans) # this is the among-timestep variation in the sample mean of p
  p_var_within_t <- mean(dat$v_pvars) # this is our estimate of the parametric variance of the distribution of p (i.e., sigma).
  
  # fundamental and effective kernel means (value + uncertainty)
  fund_kernmean_mean <- mean(dat$fund_mean_mean)
  fund_kernmean_var_among_t <- var(dat$fund_mean_mean) # variance among timesteps
  fund_kernmean_var_within_t <- mean(dat$fund_mean_var) # variance among patches within a timestep
  
  eff_kernmean_mean <- mean(dat$eff_mean_mean)
  eff_kernmean_var_among_t <- var(dat$eff_mean_mean) # variance among timesteps
  eff_kernmean_var_within_t <- mean(dat$eff_mean_var) # variance among patches within a timestep
  
  # spatial structure in fundamental and effective kernel means (value + uncertainty)
  fund_kernmoran_mean <- mean(dat$fund_mean_moran)
  fund_kernmoran_var_among_t <- var(dat$fund_mean_moran) # variance among timesteps
  
  eff_kernmoran_mean <- mean(dat$eff_mean_moran)
  eff_kernmoran_var_among_t <- var(dat$eff_mean_moran) # variance among timesteps
  
  df_return <- data.frame(p_mean=p_mean,p_var_among_t=p_var_among_t,p_var_within_t=p_var_within_t,
                          fund_kernmean_mean=fund_kernmean_mean,fund_kernmean_var_among_t=fund_kernmoran_var_among_t,fund_kernmean_var_within_t=fund_kernmean_var_within_t,
                          eff_kernmean_mean=eff_kernmean_mean,eff_kernmean_var_among_t=eff_kernmean_var_among_t,eff_kernmean_var_within_t=eff_kernmean_var_within_t,
                          fund_kernmoran_mean=fund_kernmoran_mean,fund_kernmoran_var_among_t=fund_kernmoran_var_among_t,
                          eff_kernmoran_mean=eff_kernmoran_mean,eff_kernmoran_var_among_t=eff_kernmoran_var_among_t)  
  return(df_return)
}
