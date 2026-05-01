####################### Main simulation function #######################
# Inputs:
#   params: list of biological and simulation parameters
#   hab_params: list of habitat-related parameters and objects; output of f_MakeHabitat
#   output_flag: "all" (npatch x ngroup x nsteps array of abundances) or "lite" (summary stats only for each timestep)
f_RunSimComboStorage_Parallel <- function(params, hab_params, keep=list("abund","p","kern","sp_struct"),
                                          output_flag="all",show_plot=FALSE,output_thin=1,output_file=NA,run_i=1,
                                          connmat_folder=NA,connmat_format="fst",patchdist_file=NULL,
                                          connmat_size_GB = 3.2, jobmem_GB = 450){
  
  # # Create a new folder in the MSI scratch directory
  # dir.create(connmat_folder)
  
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
  units(nav_rad) <- NULL
  
  # # patch_dists and patch_angles
  # load(file=patchdist_file)
  # units(patch_dists) <- NULL
  # patch_angles <- (2*nav_rad)/(2*pi*patch_dists)
  # patch_angles[is.nan(patch_angles)] <- 1
  # patch_angles[patch_angles>1] <- 1
  
  # connmat_size <- as.numeric(object.size(patch_dists)/(10^6)) # in Mb
  connmat_size <- connmat_size_GB*8000 # in Mb
  job_mem <- jobmem_GB*1000 # in MB
  
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

  all_conn_mats <- lapply(1:ngroups,
                            function(i) fn_readconnmat(g=i,connmat_folder=connmat_folder,job_mem=job_mem,connmat_format=connmat_format))
  
  all_conn_mats <- mclapply(1:ngroups,
                            function(i) fn_readconnmat(g=i,connmat_folder=connmat_folder,job_mem=job_mem,connmat_format=connmat_format),
                            mc.cores = parallelly::availableCores())
  # all_conn_mats <- vector("list", ngroups) # list to hold connectivity matrices, which will be computed as they're needed (will have to be reset if K_i changes)
  # for(g in 1:ngroups){
  #   # import it as conn_mat
  #   #conn_mat <- as.matrix(read_fst(paste0(connmat_folder,"/grp_",g)))
  #   conn_mat <- readRDS(paste0(connmat_folder,"/grp_",g,".rds"))
  #   # and save it in the list:
  #   x <- gc()
  #   current_used <- sum(x[,2]) # current memory usage in Mb
  #   if(current_used<0.75*job_mem){ # might mess with this criterion
  #     # save it in memory
  #     all_conn_mats[[g]] <- conn_mat
  #   } else{
  #     # save it on RAMdisk, and put the path in the all_conn_mats list
  #     # write_fst(as.data.frame(conn_mat),paste0("/dev/shm/grp_",g),compress=0)
  #     saveRDS(conn_mat,file=paste0("/dev/shm/grp_",g,".rds"),compress=FALSE)
  #     
  #     all_conn_mats[[g]] <- paste0("/dev/shm/grp_",g)
  #   }
  # }
  
  
  ## 5. Initialize output data structures
  # if(output_flag=="lite"){
  #   # df_out = data.frame(t_i=1:nsteps,abund=NA,p_mean=NA,p_var=NA,
  #   #                                  fund_kernmode_mean=NA,fund_kernmode_var=NA,
  #   #                                  fund_kernmean_mean=NA,fund_kernmean_var=NA,
  #   #                                  eff_kernmode_mean=NA,eff_kernmode_var=NA,
  #   #                                  eff_kernmean_mean=NA,eff_kernmean_var=NA)
  #   df_out <- data.frame(t_i=integer(),metric=character(),mean=numeric(),var=numeric())
  # }
  
  # objects to hold param values in different ways for more efficient use later
  p_by_group <- v_p[group_index$p] # value, not index
  alpha_by_group <- group_index$alpha # index, not value
  theta_by_group <- group_index$theta # index, not value
  theta_val_by_group <- v_thetas[theta_by_group]
  Pij_alpha <- matrix(alpha_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_theta <- matrix(theta_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_alpha_val <- v_alphas[Pij_alpha]
  Pij_theta_val <- v_thetas[Pij_theta]
  Pij_p <- matrix(p_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # value, not index
  
  # generate matrix of weights (inverse distance) to use in calculating Moran's I
  # is this an okay distance function to use? I know it matters when doing the statistical test,
  # but maybe as long as we use the same metric for all timesteps and simulations it's OK for comparing among them.
  # I tried a few options (1/d, 1/d^2, 1/(1+d^2)) and the dynamics were quite similar.
  load(file=patchdist_file)
  moran_weights <- 1/(drop_units(patch_dists))
  diag(moran_weights) <- 0
  
  ## 6. Save input information
  metadata_list <- list(params=params,group_index=group_index,patch_locations=patch_locations,hab_file=hab_params$hab_file)
  save(metadata_list,file=paste0(output_file,"_metadata.RData"))
  
  ## 7. Output first line of the summary stats file
  t_i=1
  df_i <- data.frame(t_i=rep(t_i,4),metric=NA,mean=NA,var=NA,MoranI=NA,corr_q=NA)
  
  ###### total abundance
  #df_out$abund[t_i] <- sum(previous_pop)
  df_i[1,] <- data.frame(t_i=t_i,metric="abund",mean=sum(previous_pop),var=NA,MoranI=NA,corr_q=NA)
  
  patch_pops <- rowSums(previous_pop)
  patch_full <- which(patch_pops!=0)
  group_pops <- colSums(previous_pop)
  group_full <- which(group_pops!=0)
  
  ###### p
  # mean value of p in each patch
  p_by_patch <- (previous_pop %*% p_by_group)/patch_pops # some of these might be NA
  # correlation between p and habitat quality
  cor_p_q <- cor(p_by_patch[patch_full],patch_locations$q[patch_full])
  # Moran's I for p
  p_Moran <- Moran.I(as.vector(p_by_patch),weight=moran_weights,na.rm=TRUE)$observed
  # overall mean and variance of p (value, not index) at the timestep
  p_mean <- sum(patch_pops * p_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
  p_var <- sum(group_pops*(p_by_group-p_mean)^2)/sum(previous_pop)
  df_i[2,] <- data.frame(t_i=t_i,metric="p",mean=p_mean,var=p_var,MoranI=p_Moran,corr_q=cor_p_q)
  
  ###### theta (fundamental)
  # mean value of theta in each patch (value, not index)
  theta_by_patch <- (previous_pop %*% theta_val_by_group)/patch_pops # some might be NA
  # correlation between theta and habitat quality
  cor_theta_q <- cor(theta_by_patch[patch_full],patch_locations$q[patch_full])
  # Moran's I for theta
  theta_Moran <- Moran.I(as.vector(theta_by_patch),weight=moran_weights,na.rm=TRUE)$observed
  # overall mean and variance of theta (value, not index) at the timestep
  theta_mean <- sum(patch_pops * theta_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
  theta_var <- sum(group_pops*(theta_by_group-theta_mean)^2)/sum(previous_pop)
  df_i[3,] <- data.frame(t_i=t_i,metric="theta",mean=theta_mean,var=theta_var,
                         MoranI=theta_Moran,corr_q=cor_theta_q)
  
  ###### theta (effective)
  # a Pij matrix with effective theta value in each patch/grp
  Pij_b <- matrix(patch_locations$b,byrow=FALSE,nrow=npatch,ncol=ngroups) # need in this format to feed the plasticity function
  Pij_theta_eff <- f_plasticityb(Pij_b,Pij_p,Pij_alpha,Pij_theta,length(v_alphas),length(v_thetas))$theta_plastic #index
  Pij_thetaval_eff <- matrix(v_thetas[Pij_theta_eff],nrow=nrow(Pij_theta_eff)) #value
  # mean effective theta value in each patch
  efftheta_by_patch <- rowSums(previous_pop*Pij_thetaval_eff)/rowSums(previous_pop) # some might be NA
  # correlation between effective theta value and habitat quality
  cor_efftheta_q <- cor(efftheta_by_patch[patch_full],patch_locations$q[patch_full])
  # Moran's I for p
  efftheta_Moran <- Moran.I(as.vector(efftheta_by_patch),weight=moran_weights,na.rm = TRUE)$observed
  # overall mean and variance of effective theta (value, not index) at the timestep
  efftheta_mean <- sum(previous_pop * Pij_thetaval_eff)/sum(previous_pop)
  efftheta_var <- sum(previous_pop*(Pij_thetaval_eff-efftheta_mean)^2)/sum(previous_pop)
  df_i[4,] <- data.frame(t_i=t_i,metric="efftheta",mean=efftheta_mean,var=efftheta_var,
                         MoranI=efftheta_Moran,corr_q=cor_efftheta_q)
  
  fwrite(df_i,file=paste0(output_file,"_raw.csv"),append=TRUE)
  
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
    
    #pops_all <- lapply(1:ngroups,f_ReprodDispMut_New)
    pops_all <- mclapply(1:ngroups,
                         function(i) f_ReprodDispMut_New(g=i,patch_locations=patch_locations,group_index=group_index,
                                                         pop_by_group=pop_by_group,v_p=v_p,previous_pop=previous_pop,
                                                         all_conn_mats=all_conn_mats,mutation_destinations=mutation_destinations,mu=mu),
                         mc.cores = parallelly::availableCores())
    pops_all <- simplify2array(pops_all)
    temp_pop <- apply(pops_all,1:2,sum)
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
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
    
    ## these are for output_flag==lite
    if(output_flag=="lite"){
      if(t_i %% output_thin == 0){
        if(sum(previous_pop)>0){
          df_i <- data.frame(t_i=rep(t_i,4),metric=NA,mean=NA,var=NA,MoranI=NA,corr_q=NA)
          
          ###### total abundance
          #df_out$abund[t_i] <- sum(previous_pop)
          df_i[1,] <- data.frame(t_i=t_i,metric="abund",mean=sum(previous_pop),var=NA,MoranI=NA,corr_q=NA)
          
          
          patch_pops <- rowSums(previous_pop)
          patch_full <- which(patch_pops!=0)
          group_pops <- colSums(previous_pop)
          group_full <- which(group_pops!=0)
          
          ###### p
          # mean value of p in each patch
          p_by_patch <- (previous_pop %*% p_by_group)/patch_pops # some of these might be NA
          # correlation between p and habitat quality
          cor_p_q <- cor(p_by_patch[patch_full],patch_locations$q[patch_full])
          # Moran's I for p
          p_Moran <- Moran.I(as.vector(p_by_patch),weight=moran_weights,na.rm=TRUE)$observed
          # overall mean and variance of p (value, not index) at the timestep
          p_mean <- sum(patch_pops * p_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
          p_var <- sum(group_pops*(p_by_group-p_mean)^2)/sum(previous_pop)
          df_i[2,] <- data.frame(t_i=t_i,metric="p",mean=p_mean,var=p_var,MoranI=p_Moran,corr_q=cor_p_q)
          
          ###### theta (fundamental)
          # mean value of theta in each patch (value, not index)
          theta_by_patch <- (previous_pop %*% theta_val_by_group)/patch_pops # some might be NA
          # correlation between theta and habitat quality
          cor_theta_q <- cor(theta_by_patch[patch_full],patch_locations$q[patch_full])
          # Moran's I for theta
          theta_Moran <- Moran.I(as.vector(theta_by_patch),weight=moran_weights,na.rm=TRUE)$observed
          # overall mean and variance of theta (value, not index) at the timestep
          theta_mean <- sum(patch_pops * theta_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
          theta_var <- sum(group_pops*(theta_by_group-theta_mean)^2)/sum(previous_pop)
          df_i[3,] <- data.frame(t_i=t_i,metric="theta",mean=theta_mean,var=theta_var,
                                 MoranI=theta_Moran,corr_q=cor_theta_q)
          
          ###### theta (effective)
          # a Pij matrix with effective theta value in each patch/grp
          Pij_b <- matrix(patch_locations$b,byrow=FALSE,nrow=npatch,ncol=ngroups) # need in this format to feed the plasticity function
          Pij_theta_eff <- f_plasticityb(Pij_b,Pij_p,Pij_alpha,Pij_theta,length(v_alphas),length(v_thetas))$theta_plastic #index
          Pij_thetaval_eff <- matrix(v_thetas[Pij_theta_eff],nrow=nrow(Pij_theta_eff)) #value
          # mean effective theta value in each patch
          efftheta_by_patch <- rowSums(previous_pop*Pij_thetaval_eff)/rowSums(previous_pop) # some might be NA
          # correlation between effective theta value and habitat quality
          cor_efftheta_q <- cor(efftheta_by_patch[patch_full],patch_locations$q[patch_full])
          # Moran's I for p
          efftheta_Moran <- Moran.I(as.vector(efftheta_by_patch),weight=moran_weights,na.rm = TRUE)$observed
          # overall mean and variance of effective theta (value, not index) at the timestep
          efftheta_mean <- sum(previous_pop * Pij_thetaval_eff)/sum(previous_pop)
          efftheta_var <- sum(previous_pop*(Pij_thetaval_eff-efftheta_mean)^2)/sum(previous_pop)
          df_i[4,] <- data.frame(t_i=t_i,metric="efftheta",mean=efftheta_mean,var=efftheta_var,
                                 MoranI=efftheta_Moran,corr_q=cor_efftheta_q)
          
          fwrite(df_i,file=paste0(output_file,"_raw.csv"),append=TRUE)
        }
      }
    } # if output_flag=="lite"
    
    
  } #t_i
  
  time_run <- proc.time()-starttime
  return(time_run)
}


############# function to do reproduction, dispersal, and mutation steps in parallel #############
# (not using this; it actually slows things down)
f_ReprodDispMut_New <- function(g,patch_locations,group_index,pop_by_group,v_p,previous_pop,
                                all_conn_mats,mutation_destinations,mu){
  # hold intermediate population values from this group
  temp_pop_g <- matrix(0,nrow=nrow(patch_locations),ncol=nrow(group_index)) 
  
  if(pop_by_group[g]>0){
    # get parameter values for that parameter group
    v <- group_index[g,]
    #system(paste0("echo '",v$p, "'"))
    p_penalty <- abs(v_p[v$p])*0
    #p_penalty <- 0
    
    # get population of each patch for that parameter group
    patch_pops <- previous_pop[,g]
    
    # get connectivity matrix from list
    # If it's already a connectivity matrix, we're done.
    conn_mat <- all_conn_mats[[g]]
    # if it's just the filepath, get the object
    if(is.character(all_conn_mats[[g]])){
      # find the appropriate file on RAMdisk, and save it as the conn_mat object
      conn_mat <- readRDS(paste0(conn_mat,".rds"))
      #conn_mat <- as.matrix(read_fst(conn_mat))
    }
    
    to_patch <- (1-p_penalty)*patch_locations$b*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
    
    # Divide up to_patch among parameter groups that are the result of mutation
    temp_pop_g[,g] <- (1-mu)*to_patch+temp_pop_g[,g]
    
    for(mut_group in mutation_destinations[g,-1]){ # for each of the possible mutations.
      temp_pop_g[,mut_group] <- (mu/6)*to_patch+temp_pop_g[,mut_group]
    }
  }
  
  return(temp_pop_g)
}

fn_readconnmat <- function(g,connmat_folder,job_mem,connmat_format){
  x <- gc()
  current_used <- sum(x[,2]) # current memory usage in Mb
  connmat_file <- ifelse(connmat_format=="rds",paste0("grp_",g,".rds"),paste0("grp_",g))
  
  if(current_used<0.75*job_mem){ # might mess with this criterion
    if(connmat_format=="rds"){
      return(readRDS(file.path(connmat_folder,connmat_file)))      
    } else{return(as.matrix(read_fst(file.path(connmat_folder,connmat_file))))}
  
  } else{
    if(connmat_format=="rds"){
      saveRDS(readRDS(connmat_file),file=file.path("/dev/shm",connmat_file),compress=FALSE)
    } else {write_fst(as.data.frame(conn_mat),file.path("/dev/shm",connmat_file),compress=0)}
    
    return(file.path("/dev/shm",connmat_file))
  }
}
