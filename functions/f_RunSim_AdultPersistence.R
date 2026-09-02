######################## Main simulation function #######################
# Inputs:
#   params: list of biological and simulation parameters
#   hab_params: list of habitat-related parameters and objects; output of f_MakeHabitat
#   output_flag: "all" (npatch x ngroup x nsteps array of abundances) or "lite" (summary stats only for each timestep)
f_RunSim_AdultPersistence <- function(params, hab_params, keep=list("abund","p","kern","sp_struct"),
                                     output_flag="all",show_plot=FALSE,output_thin=1,output_file=NA,run_i=1,
                                     connmat_folder=NA,connmat_format="fst",patchdist_file=NULL,
                                     connmat_size_GB=3.2,jobmem_GB=400,start_location="everywhere",normalize_connmats=FALSE,
                                     disturb_freq=10,disturb_radius=1){
  starttime <- proc.time()
  numCores <- parallelly::availableCores()
  n_connmats_tostore <- floor((1*jobmem_GB)/connmat_size_GB)
  
  # load parameters
  list2env(x=params,envir=environment())
  # load habitat data structures (see f_MakeHabitat for details on what's in each object)
  list2env(x=hab_params,envir=environment()) 
  load(file=patchdist_file)
  patch_dists <- collapse::qM(patch_dists) # remove units, but it's in km
  if(hab_type=="points") patch_locations$K <- 1 # sometimes these end up as 0 from map resolution issues
  # save original values, in case they're modified by disturbance
  K <- patch_locations$K 
  b <- patch_locations$b
  q <- patch_locations$q
  units(nav_rad) <- NULL
  
  
  patch_locations$vac_prob <- 0.5
  
  connmat_size <- connmat_size_GB*1000*8 # in Mb
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
  if(start_location=="local"){
    start_spot <- st_sample(reef_sf,1)
    dist_thresh <- 100
    units(dist_thresh) <- "m"
    start_patches <- which(st_distance(sfc_patches,start_spot)<dist_thresh)
  } else start_patches <- 1:npatch
  start_per_patch <- rep(0,npatch)
  start_per_patch[start_patches] <- patch_locations$K[start_patches]
  starts <- lapply(1:npatch,function(i) as.vector(rmultinom(n=1,size=start_per_patch[i],prob=start_probs))) 
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
  
  ## 6. Save input information
  metadata_list <- list(params=params,group_index=group_index,patch_locations=patch_locations,hab_file=hab_params$hab_file)
  save(metadata_list,file=paste0(output_file,"_metadata.RData"))
  
  ## 7. Output first line of the summary stats file
  t_i=1
  df_i <- data.frame(t_i=rep(t_i,4),metric=NA,mean=NA,var=NA,q05=NA,q25=NA,median=NA,q75=NA,q95=NA,corr_q=NA)
  
  ###### total abundance
  df_i[1,] <- data.frame(t_i=t_i,metric="abund",mean=NA,var=NA,q05=NA,q25=NA,median=sum(previous_pop),q75=NA,q95=NA,corr_q=NA)
  
  patch_pops <- rowSums(previous_pop)
  patch_full <- which(patch_pops!=0)
  group_pops <- colSums(previous_pop)
  group_full <- which(group_pops!=0)
  ###### p
  # mean value of p in each patch
  p_by_patch <- (previous_pop %*% p_by_group)/patch_pops # some of these might be NA
  # correlation between p and habitat quality
  cor_p_q <- cor(p_by_patch[patch_full],patch_locations$q[patch_full])
  # overall mean and variance of p (value, not index) at the timestep
  # p_median is a weighted median of p_by_group, weighted by group_pops
  p_quants <- as.numeric(Hmisc::wtd.quantile(x = p_by_group, weights = group_pops, probs = c(0.05,0.25,0.5,0.75,0.95)))
  p_mean <- sum(patch_pops * p_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
  p_var <- sum(group_pops*(p_by_group-p_mean)^2)/sum(previous_pop)
  df_i[2,] <- data.frame(t_i=t_i,metric="p",mean=p_mean,var=p_var,
                         q05=p_quants[1],q25=p_quants[2],median=p_quants[3],q75=p_quants[4],q95=p_quants[5],
                         corr_q=cor_p_q)
  
  ###### theta (fundamental)
  # mean value of theta in each patch (value, not index)
  theta_by_patch <- (previous_pop %*% theta_val_by_group)/patch_pops # some might be NA
  # correlation between theta and habitat quality
  cor_theta_q <- cor(theta_by_patch[patch_full],patch_locations$q[patch_full])
  # overall mean and variance of theta (value, not index) at the timestep
  theta_quants <- as.numeric(Hmisc::wtd.quantile(x = theta_val_by_group, weights = group_pops, probs = c(0.05,0.25,0.5,0.75,0.95)))
  theta_mean <- sum(patch_pops * theta_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
  theta_var <- sum(group_pops*(theta_by_group-theta_mean)^2)/sum(previous_pop)
  df_i[3,] <- data.frame(t_i=t_i,metric="theta",mean=theta_mean,var=theta_var,
                         q05=theta_quants[1],q25=theta_quants[2],median=theta_quants[3],q75=theta_quants[4],q95=theta_quants[5],
                         corr_q=cor_theta_q)
  
  ###### theta (effective)
  # a Pij matrix with effective theta value in each patch/grp
  Pij_b <- matrix(patch_locations$b,byrow=FALSE,nrow=npatch,ncol=ngroups) # need in this format to feed the plasticity function
  Pij_theta_eff <- f_plasticityb(Pij_b,Pij_p,Pij_alpha,Pij_theta,length(v_alphas),length(v_thetas))$theta_plastic #index
  Pij_thetaval_eff <- matrix(v_thetas[Pij_theta_eff],nrow=nrow(Pij_theta_eff)) #value
  # mean effective theta value in each patch
  efftheta_by_patch <- rowSums(previous_pop*Pij_thetaval_eff)/rowSums(previous_pop) # some might be NA
  # correlation between effective theta value and habitat quality
  cor_efftheta_q <- cor(efftheta_by_patch[patch_full],patch_locations$q[patch_full])
  # overall mean and variance of effective theta (value, not index) at the timestep
  efftheta_quants <- as.numeric(Hmisc::wtd.quantile(x = Pij_thetaval_eff, weights = previous_pop, probs = c(0.05,0.25,0.5,0.75,0.95)))
  efftheta_mean <- sum(previous_pop * Pij_thetaval_eff)/sum(previous_pop)
  efftheta_var <- sum(previous_pop*(Pij_thetaval_eff-efftheta_mean)^2)/sum(previous_pop)
  df_i[4,] <- data.frame(t_i=t_i,metric="efftheta",mean=efftheta_mean,var=efftheta_var,
                         q05=efftheta_quants[1],q25=efftheta_quants[2],median=efftheta_quants[3],q75=efftheta_quants[4],q95=efftheta_quants[5],
                         corr_q=cor_efftheta_q)
  
  fwrite(df_i,file=paste0(output_file,"_raw.csv"),append=TRUE)
  
  to_patch <- vector(length = npatch)
  
  ################ Simulate ###################
  interval_starttime <- proc.time()
  for(t_i in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    patch_abunds_adult <- rowSums(previous_pop)
    fwrite(matrix(patch_abunds_adult,nrow=1),file=paste0(output_file,"_patch_abunds_adult.csv"),append=TRUE)
    
    ################## Disturbance ##################
    
    if(t_i %% disturb_freq == 0){ # disturbances every disturb_freq timesteps
      # make a new disturbance
        # choose an origin anemone
        disturbance_center <- sample(1:npatch,1)
        # find all anemones within a certain distance from the center
        disturbed_anems <- as.numeric(which(patch_dists[,disturbance_center]<disturb_radius))
        # remove all residents of those anemones
        previous_pop[disturbed_anems,] <- 0
    }
    
    ################## Reproduction and Dispersal and Mutation ##################
    pop_by_group <- colSums(previous_pop)
    
    # decide which anemones will have a vacancy for new larvae to settle this round
    will_have_vacancy <- which(simDAG::rbernoulli(n=nrow(patch_locations),p=patch_locations$vac_prob))
    
    #print(paste0("timestep: ",t_i))
    for(g in which(pop_by_group>0)){
      # get parameter values for that parameter group
      v <- group_index[g,]
      p_penalty <- abs(v_p[v$p])*0
      
      # get population of each patch for that parameter group
      patch_pops <- previous_pop[,g]
      
      # get the connectivity matrix among patches, given the group parameter values and patch-level q's
      # (and accounting for the patch population x per capita output b_i from each patch)
      
      # if it hasn't been imported from global.scratch yet
      if(is.null(all_conn_mats[[g]])){
        #print("importing")
        # import it as conn_mat
        if(connmat_format=="fst"){conn_mat <- collapse::qM(read_fst(paste0(connmat_folder,"/grp_",g)))
        } else if(connmat_format=="rds"){conn_mat <- readRDS(paste0(connmat_folder,"/grp_",g,".rds"))}
        # if we want to normalize them, do that:
        if(normalize_connmats==TRUE) conn_mat <- sweep(conn_mat,2,colSums(conn_mat),FUN='/') 
        # and then save it in the list:
        if(g<n_connmats_tostore){ #(current_used<0.75*job_mem){ # might mess with this criterion
          # save it in memory
          all_conn_mats[[g]] <- conn_mat
          #print("saving to memory")
        } else{
          # save it on RAMdisk, and put the path in the all_conn_mats list
          write_fst(collapse::qDF(conn_mat),paste0("/dev/shm/grp_",g),compress=0)
          # saveRDS(conn_mat,file=paste0("/dev/shm/grp_",g,".rds"),compress=FALSE)
          all_conn_mats[[g]] <- paste0("/dev/shm/grp_",g)
          #print("saving to disk (/dev/shm)")
        }
        # if it has already been imported from global scratch
      } else{
        # grab what's in the list. If this is already a connectivity matrix, we're done.
        conn_mat <- all_conn_mats[[g]]
        #print("exists in memory")
        # if it's just the filepath, get the object
        if(is.character(all_conn_mats[[g]])){
          # find the appropriate file on RAMdisk, and save it as the conn_mat object
          conn_mat <- collapse::qM(read_fst(conn_mat))
          #conn_mat <- readRDS(conn_mat)
          #print("reading from disk")
        }
      }
      
      occupied_patches <- which(patch_pops>0)
      
      # vector of contribution of the population of this group to each patch
      to_patch[ ] <- NA # to_patch and temp_pop are just keeping track of larvae, so we'll keep the rows of occupied patches at NA.
      to_patch[will_have_vacancy] <- (1-p_penalty)*(conn_mat[will_have_vacancy,occupied_patches,drop=FALSE] %*% (patch_pops[occupied_patches]*patch_locations$b[occupied_patches])) 
      # Pij_larvae[,g] <- to_patch
      
      # Divide up to_patch among parameter groups that are the result of mutation
      temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
      for(mut_group in mutation_destinations[g,-1]){ # for each of the possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
        temp_pop[,mut_group] <- (mu/6)*to_patch+temp_pop[,mut_group]
      }
    } #g
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    patch_abunds <- rowSums(temp_pop)
    fwrite(matrix(patch_abunds,nrow=1),file=paste0(output_file,"_patch_abunds.csv"),append=TRUE)
    comp_results <- vapply(will_have_vacancy,function(i) f_Competition(i_patch=i,patch_abunds=patch_abunds,patch_locations=patch_locations,
                                                              temp_pop=temp_pop,ngroups=ngroups),
                           integer(ngroups))
    previous_pop[will_have_vacancy,] <- t(comp_results)
    new_pop <- previous_pop
    
    ################## Output ##################
    if(t_i %% max(1,round(nsteps/100)) == 0){
      # status updates to console every so often
      print(t_i)
      print(proc.time()-interval_starttime)
      interval_starttime <- proc.time()
      
      # plot, if requested
      if(show_plot==TRUE){
        if(hab_type=="points"){
          grp_by_patch <- sapply(asplit(new_pop,MARGIN=1),which.max)
          theta_by_patch <- theta_by_group[grp_by_patch]
          p_by_patch <- p_by_group[grp_by_patch]
          
          g_map_abund <- ggplot(reef_sf)+
            geom_sf()+
            geom_sf(data=sfc_patches,aes(color=factor(rowSums(new_pop))))+
            labs(color="Abundance",title=paste0("t_i=",t_i))+
            scale_color_manual(values=c("blue","red"),breaks=c(1,0))+
            theme_minimal()+
            annotation_scale()
          g_map_q <- ggplot(reef_sf)+
            geom_sf()+
            geom_sf(data=sfc_patches,aes(color=patch_locations$q))+
            labs(color="Hab quality")+
            theme_minimal()+
            annotation_scale()
          g_map_p <- ggplot(reef_sf)+
            geom_sf()+
            geom_sf(data=sfc_patches,aes(color=p_by_patch))+
            labs(color="Plasticity")+
            theme_minimal()+
            annotation_scale()
          g_map_theta <- ggplot(reef_sf)+
            geom_sf()+
            geom_sf(data=sfc_patches,aes(color=theta_by_patch))+
            labs(color="Kernel mean")+
            theme_minimal()+
            annotation_scale()
          g_all <- grid.arrange(g_map_q,g_map_p,g_map_theta,top=paste0("t_i=",t_i))
          
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
    
    if(t_i %% 100 == 0){
      save(previous_pop,t_i,file=paste0(output_file,"_popsnapshot.RData"))
    }
    
    ## these are for output_flag==lite
    if(output_flag=="lite"){
      if(t_i %% output_thin == 0){
        if(sum(previous_pop)>0){
          
          df_i <- data.frame(t_i=rep(t_i,5),metric=NA,mean=NA,var=NA,q05=NA,q25=NA,median=NA,q75=NA,q95=NA,corr_q=NA)
          
          ###### total abundance
          df_i[1,] <- data.frame(t_i=t_i,metric="abund",mean=NA,var=NA,q05=NA,q25=NA,median=sum(previous_pop),q75=NA,q95=NA,corr_q=NA)
          
          patch_pops <- rowSums(previous_pop)
          patch_full <- which(patch_pops!=0)
          group_pops <- colSums(previous_pop)
          group_full <- which(group_pops!=0)
          ###### p
          # mean value of p in each patch
          p_by_patch <- (previous_pop %*% p_by_group)/patch_pops # some of these might be NA
          p_val_by_patch <- (v_thetas[2]/v_thetas[1])^-p_by_patch
          # correlation between p and habitat quality
          cor_p_q <- cor(p_val_by_patch[patch_full],patch_locations$q[patch_full])
          # overall mean and variance of p (value, not index) at the timestep
          # p_median is a weighted median of p_by_group, weighted by group_pops
          p_quants <- as.numeric(Hmisc::wtd.quantile(x = p_by_group, weights = group_pops, probs = c(0.05,0.25,0.5,0.75,0.95)))
          p_mean <- sum(patch_pops * p_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
          p_var <- sum(group_pops*(p_by_group-p_mean)^2)/sum(previous_pop)
          df_i[2,] <- data.frame(t_i=t_i,metric="p",mean=p_mean,var=p_var,
                                 q05=p_quants[1],q25=p_quants[2],median=p_quants[3],q75=p_quants[4],q95=p_quants[5],
                                 corr_q=cor_p_q)
          
          ###### theta (fundamental)
          # mean value of theta in each patch (value, not index)
          theta_by_patch <- (previous_pop %*% theta_val_by_group)/patch_pops # some might be NA
          # correlation between theta and habitat quality
          cor_theta_q <- cor(theta_by_patch[patch_full],patch_locations$q[patch_full])
          # overall mean and variance of theta (value, not index) at the timestep
          theta_quants <- as.numeric(Hmisc::wtd.quantile(x = theta_val_by_group, weights = group_pops, probs = c(0.05,0.25,0.5,0.75,0.95)))
          theta_mean <- sum(patch_pops * theta_by_patch,na.rm=TRUE)/sum(patch_pops,na.rm=TRUE)
          theta_var <- sum(group_pops*(theta_by_group-theta_mean)^2)/sum(previous_pop)
          df_i[3,] <- data.frame(t_i=t_i,metric="theta",mean=theta_mean,var=theta_var,
                                 q05=theta_quants[1],q25=theta_quants[2],median=theta_quants[3],q75=theta_quants[4],q95=theta_quants[5],
                                 corr_q=cor_theta_q)
          
          ###### theta (effective)
          # a Pij matrix with effective theta value in each patch/grp
          Pij_b <- matrix(patch_locations$b,byrow=FALSE,nrow=npatch,ncol=ngroups) # need in this format to feed the plasticity function
          Pij_theta_eff <- f_plasticityb(Pij_b,Pij_p,Pij_alpha,Pij_theta,length(v_alphas),length(v_thetas))$theta_plastic #index
          Pij_thetaval_eff <- matrix(v_thetas[Pij_theta_eff],nrow=nrow(Pij_theta_eff)) #value
          # mean effective theta value in each patch
          efftheta_by_patch <- rowSums(previous_pop*Pij_thetaval_eff)/rowSums(previous_pop) # some might be NA
          # correlation between effective theta value and habitat quality
          cor_efftheta_q <- cor(efftheta_by_patch[patch_full],patch_locations$q[patch_full])
          # overall mean and variance of effective theta (value, not index) at the timestep
          efftheta_quants <- as.numeric(Hmisc::wtd.quantile(x = Pij_thetaval_eff, weights = previous_pop, probs = c(0.05,0.25,0.5,0.75,0.95)))
          efftheta_mean <- sum(previous_pop * Pij_thetaval_eff)/sum(previous_pop)
          efftheta_var <- sum(previous_pop*(Pij_thetaval_eff-efftheta_mean)^2)/sum(previous_pop)
          df_i[4,] <- data.frame(t_i=t_i,metric="efftheta",mean=efftheta_mean,var=efftheta_var,
                                 q05=efftheta_quants[1],q25=efftheta_quants[2],median=efftheta_quants[3],q75=efftheta_quants[4],q95=efftheta_quants[5],
                                 corr_q=cor_efftheta_q)
          
          df_i[5,] <- data.frame(t_i=t_i,metric="larval_abund",mean=NA,var=NA,
                                 q05=NA,q25=NA,median=sum(patch_abunds),q75=NA,q95=NA,
                                 corr_q=NA)
          
          fwrite(df_i,file=paste0(output_file,"_raw.csv"),append=TRUE)
        }
      }
    } # if output_flag=="lite"
    
    
  } #t_i
  
  time_run <- proc.time()-starttime
  return(time_run)
}
