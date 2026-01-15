####################### version that stores all timesteps #######################
f_RunMatrixLoop <- function(params){
  starttime <- proc.time()
  numCores <- parallelly::availableCores()
  
  list2env(x=params,envir=environment())
  
  # # make each run replicable, or use seed from existing run
  # if(is.null(seed)) seed = .Random.seed
  # set.seed(seed)
  
  # Data structures to describe space and dispersal
  # -------------------------------------------------------------------
  # see utility functions/f_MakeHabitat for details on what's in each object
  list2env(x=f_MakeHabitat(nx,ny,v_alphas,v_thetas,patch_locations),envir=environment())
  if(!"K_i" %in% colnames(patch_locations)) patch_locations$K_i <- as.vector(K)
  patch_locations$b_i <- as.vector(b)
  params$patch_locations <- patch_locations
  
  # Data structures to describe population
  # -------------------------------------------------------------------
  
  ## 1. group_index: all unique combinations of parameters alpha, theta, and p
  group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  ngroups <- nrow(group_index)
  
  ## 2. Pij: population object
  ## dimensions 1 = npatch, 2 = ngroups, 3 = nsteps
  Pij <- array(0, dim=c(npatch,ngroups,nsteps))
  # initialize Pij
  start_grp <- which(group_index$alpha==alpha_start & group_index$theta==theta_start & group_index$p==p_start)
  Pij[,start_grp,1] <- patch_locations$K_i
  
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
  
  # Simulate
  # -------------------------------------------------------------------
  for(t in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    ################## Disturbance ##################
    
    if(t %% 10 == 0){ # disturbances last for 10 timesteps
      # first, reset all K's from any disturbance that occurred 10 timesteps previously
      if(prod(patch_locations$K_i == as.vector(K))==0){
        patch_locations$K_i = as.vector(K)
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K_i
      }
      
      # then make a new disturbance (maybe)
      if(rbinom(n=1,size=1,p=disturb_prob)==1){
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K_i
        
        disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
        disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
        patch_locations$K_i[disturb_patches] <- 0
      }
    }
    
    ################## Reproduction and Dispersal and Mutation ##################
    pop_by_group <- colSums(Pij[,,t-1])
    
    for(g in which(pop_by_group>0)){
      # get parameter values for that parameter group
      v <- group_index[g,]
      p_penalty <- abs(v_p[v$p])*0
      
      # get population of each patch for that parameter group
      patch_pops <- Pij[,g,t-1]
      
      # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
      # (and accounting for the patch population x per capita output b_i from each patch)
      if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas, overlap_discount, numCores) # calculate this matrix, if it hasn't been used before
      conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
      to_patch <- (1-p_penalty)*patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
      
      # Divide up to_patch among parameter groups that are the result of mutation
      temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
      
      for(mut_group in mutation_destinations[g,-1]){ # for each of the possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
        temp_pop[,mut_group] <- (mu/4)*to_patch+temp_pop[,mut_group]
      }
    }
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    patch_abunds <- rowSums(temp_pop)
    for(i_patch in which(patch_abunds>0)){
      survivors=rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
                          size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                          prob = temp_pop[i_patch,]) # probability of each cell being chosen depends on its current abundance)
      Pij[i_patch,,t] <- survivors
    }
    
    if(t %% max(1,round(nsteps/10)) == 0) print(t)
  }
  
  time_run <- proc.time()-starttime
  return(list(params=params,Pij=Pij,group_index=group_index,time_run=time_run,seed=seed))
}


######################## version that only stores summary stats #######################
f_RunMatrixLoopLite <- function(params, keep=list("abund","p","kern","sp_struct"),show_plot=FALSE,makehab_output=NULL){
  starttime <- proc.time()
  numCores <- parallelly::availableCores()
  
  list2env(x=params,envir=environment())
  
  # make each run replicable, or use seed from existing run
  # if(is.null(seed)) seed = .Random.seed
  # set.seed(seed)
  
  # Data structures to describe space and dispersal
  # -------------------------------------------------------------------
  # see utility functions/f_MakeHabitat for details on what's in each object
  if(is.null(makehab_output)) makehab_output <- f_MakeHabitat(nx=nx,ny=ny,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                                                              hab_type=hab_type,nav_rad=nav_rad,numCores=numCores)
  list2env(x=makehab_output,envir=environment())

  # set patch-specific characteristics
  if(hab_type=="points") patch_locations$K_i <- 1
  K <- patch_locations$K_i # K records the initial values; the values in patch_locations$K_i could change with disturbance
  b <- patch_locations$b_i # b records the initial values; the values in patch_locations$b_i could change with disturbance
  # put habitat info into params so it's accessible in the simulation output
  params$patch_locations <- patch_locations
  params$patch_sf <- patch_sf
  
  # Data structures to describe population
  # -------------------------------------------------------------------
  
  ## 1. group_index: all unique combinations of parameters alpha, theta, and p
  group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  ngroups <- nrow(group_index)
  
  ## 2. Population objects
  ## dimensions 1 = npatch, 2 = ngroups, 3 = nsteps
  previous_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  new_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  
  # initialize previous_pop
  if(alpha_start==0) alpha_start = 1:length(v_alphas)
  if(theta_start==0) theta_start = 1:length(v_thetas)
  if(p_start==0) p_start = 1:length(v_p)
  start_grps <- which((group_index$alpha %in% alpha_start) & (group_index$theta %in% theta_start) & (group_index$p %in% p_start))
  start_probs <- rep(0,ngroups)
  start_probs[start_grps] <- 1
  starts <- lapply(1:npatch,function(i) as.vector(rmultinom(n=1,size=patch_locations$K_i[i],prob=start_probs))) 
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
  output_list=list()
  if("abund" %in% keep){
    output_list$df_abund = data.frame(t=1:nsteps,abund=NA)
  }
  if("p" %in% keep){
    output_list$df_p = data.frame(t=1:nsteps,mean=NA,var=NA)
  }
  if("kern" %in% keep){
    output_list$df_fund = data.frame(t=1:nsteps,kernmode_mean=NA,kernmode_var=NA,kernmean_mean=NA,kernmean_var=NA,kernmean_moran=NA)
    output_list$df_eff = data.frame(t=1:nsteps,kernmode_mean=NA,kernmode_var=NA,kernmean_mean=NA,kernmean_var=NA,kernmean_moran=NA)
  }
  
  p_by_group <- v_p[group_index$p] # value, not index
  alpha_by_group <- group_index$alpha # index, not value
  theta_by_group <- group_index$theta # index, not value
  Pij_alpha <- matrix(alpha_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_theta <- matrix(theta_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # index, not value
  Pij_alpha_val <- v_alphas[Pij_alpha]
  Pij_theta_val <- v_thetas[Pij_theta]
  Pij_p <- matrix(p_by_group,byrow=TRUE,nrow=npatch,ncol=ngroups) # value, not index
  # generate matrix of weights (inverse distance) to use in calculating Moran's I
  # is this an okay distance function to use? I know it matters when doing the statistical test,
  # but maybe as long as we use the same metric for all timesteps and simulations it's OK for comparing among them.
  # I tried a few options (1/d, 1/d^2, 1/(1+d^2)) and the dynamics were quite similar.
  moran_weights <- 1/(patch_dists) 
  diag(moran_weights) <- 0
  
  # Simulate
  # -------------------------------------------------------------------
  interval_starttime <- proc.time()
  for(t in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    ################## Disturbance ##################
    
    if(t %% 10 == 0){ # disturbances last for 10 timesteps
      # first, reset all K's from any disturbance that occurred 10 timesteps previously
      if(prod(patch_locations$K_i == as.vector(K))==0){
        patch_locations$K_i = as.vector(K)
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K_i
        print(paste("t =",t,"reset connectivity matrices"))
      }
      
      # then make a new disturbance (maybe)
      if(rbinom(n=1,size=1,p=disturb_prob)==1){
        print(paste("t =",t,"disturbance"))
        all_conn_mats <- vector("list", ngroups) # need to reset this if we're changing K_i
        
        disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
        disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
        patch_locations$K_i[disturb_patches] <- 0
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
      to_patch <- (1-p_penalty)*patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
      
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
    if(t %% max(1,round(nsteps/10)) == 0){
      print(t)
      print(proc.time()-interval_starttime)
      interval_starttime <- proc.time()
      
      if(show_plot==TRUE){
        if(hab_type=="points"){
          g_map_abund <- ggplot(patch_sf)+
            geom_sf(aes(color=factor(rowSums(new_pop))))+
            labs(color="Abundance",title=paste0("t=",t))+
            scale_color_manual(values=c("blue","red"),breaks=c(1,0))+
            theme_minimal()+
            annotation_scale()
          print(g_map_abund)
        }
      }
    } 
    
    if(sum(previous_pop)>0){
      if("abund" %in% keep){
        output_list$df_abund$abund[t] <- sum(previous_pop)
      }
      
      if("p" %in% keep){
        # mean value of p at each timestep
        p_mean <- sum(p_by_group*colSums(previous_pop))/sum(previous_pop)
        output_list$df_p$mean[t] <- p_mean
        # variance of p within the population at each timestep
        output_list$df_p$var[t] <- sum(colSums(previous_pop)*(p_by_group-p_mean)^2)/sum(previous_pop)
      }
      
      if("kern" %in% keep){
        Pij_b <- matrix(patch_locations$b_i,byrow=FALSE,nrow=npatch,ncol=ngroups)
        Pij_eff <- f_plasticityb(Pij_b,Pij_p,Pij_alpha,Pij_theta,length(v_alphas),length(v_thetas)) # indices, not values
        
        ## fundamental kernel properties (i.e., based on inherited alpha and theta, but not plasticity)
        ### 1. kernel mode
        fund_mode_each <- ifelse(Pij_alpha_val<1,0,(Pij_alpha_val-1)*Pij_theta_val)
        # mean value (over the population) of the fundamental kernel mode at each timestep
        fund_mode_mean <- sum(fund_mode_each*previous_pop)/sum(previous_pop)
        output_list$df_fund$kernmode_mean[t] <- fund_mode_mean
        # variance (over the population) of the fundamental kernel mode
        output_list$df_fund$kernmode_var[t] <- sum(previous_pop*(fund_mode_each-fund_mode_mean)^2)/sum(previous_pop)
        ### 2. kernel mean
        fund_mean_each <- Pij_alpha_val*Pij_theta_val
        # mean value (over the population) of the fundamental kernel mean at each timestep
        fund_mean_mean <- sum(fund_mean_each*previous_pop)/sum(previous_pop)
        output_list$df_fund$kernmean_mean[t] <- fund_mean_mean
        # variance (over the population) of the fundamental kernel mean
        output_list$df_fund$kernmean_var[t] <- sum(previous_pop*(fund_mean_each-fund_mean_mean)^2)/sum(previous_pop)
        
        ## effective kernel properties (i.e., accounting for plasticity)
        ### 1. kernel mode -- something is weird here!!
        eff_mode_each <- ifelse(v_alphas[Pij_eff$alpha_plastic]<1,0,
                                (v_alphas[Pij_eff$alpha_plastic]-1)*v_alphas[Pij_eff$theta_plastic])
        # mean value (over the population) of the effective kernel mode at each timestep
        eff_mode_mean <- sum(eff_mode_each*previous_pop)/sum(previous_pop)
        output_list$df_eff$kernmode_mean[t] <- eff_mode_mean
        # variance (over the population) of the effective kernel mode
        output_list$df_eff$kernmode_var[t] <- sum(previous_pop*(eff_mode_each-eff_mode_mean)^2)/sum(previous_pop)
        ### 2. kernel mean
        eff_mean_each <- v_alphas[Pij_eff$alpha_plastic]*v_thetas[Pij_eff$theta_plastic]
        # mean value (over the population) of the effective kernel mean at each timestep
        eff_mean_mean <- sum(eff_mean_each*previous_pop)/sum(previous_pop)
        output_list$df_eff$kernmean_mean[t] <- eff_mean_mean
        # variance (over the population) of the effective kernel mean at each timestep
        output_list$df_eff$kernmean_var[t] <- sum(previous_pop*(eff_mean_each-eff_mean_mean)^2)/sum(previous_pop)
        
        if("sp_struct" %in% keep){
          # some sort of measure of the spatial autocorrelation of fundamental kernel mean
          fund_mean_by_patch <- as.vector(rowSums(fund_mean_each*previous_pop))
          output_list$df_fund$kernmean_moran[t] <- Moran.I(fund_mean_by_patch,weight=moran_weights)$observed
          
          # some sort of measure of the spatial autocorrelation of effective kernel mean
          eff_mean_by_patch <- as.vector(rowSums(eff_mean_each*previous_pop))
          output_list$df_eff$kernmean_moran[t] <- Moran.I(eff_mean_by_patch,weight=moran_weights)$observed
          
          # slope of the line relating (among patches) distance to open water and kernel mode?
        }
      }
    } # if sum(previous_pop)>0
    
  } #t
  
  
  time_run <- proc.time()-starttime
  return(list(output_list=output_list,params=params,time_run=time_run))
}

######################## function to do competition step in parallel #######################
f_Competition <- function(i_patch,patch_abunds,patch_locations,temp_pop,ngroups){
  if(patch_abunds[i_patch]>0){ # if the patch isn't empty
    survivors=t(rmultinom(n=1,
                          size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose groups for min(abundance, K) survivors
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


############## process raw output of all-timesteps simulation function #############
f_ProcessLoopOutput <- function(params,Pij,group_index,time_run){
  list2env(x=params,envir=environment())
  # melt into a dataframe with columns patch, timestep, alpha, theta, p, popsize
  # add columns for param values at each time/patch/alpha/theta/p combo, scaled by the population size that has that combo
  group_index <- mutate(group_index,group=1:nrow(group_index))
  sim_melt <- reshape2::melt(Pij, varnames=c("patch","group",'t'),value.name="popsize") %>%
    left_join(group_index,by='group')
  sim_melt <- sim_melt[,c("patch", "t", "alpha", "theta", "p", "popsize")] %>%
    mutate(alpha_value=v_alphas[alpha],theta_value=v_thetas[theta],p_value=v_p[p],t=as.numeric(t))
  
  # mean param values at each timepoint
  # first take the sum across all cells of the param*popsize (numerator of the mean)
  by_t <- summarize(group_by(sim_melt,t),
                    alpha=sum(alpha_value*popsize),theta=sum(theta_value*popsize),p=sum(p_value*popsize),popsize=sum(popsize))
  # then divide by total popsize (denominator of the mean)
  by_t <- mutate(by_t, alpha=alpha/popsize, theta=theta/popsize, p=p/popsize)
  
  return(list(sim_melt=sim_melt,
              by_t=by_t,
              patch_locations=patch_locations,
              time_run=time_run,
              K=K,
              b=b))
}

############## process raw output of all-timesteps simulation function, faster #############
# use data.table instead of dplyr
f_ProcessLoopOutputDataTable <- function(params,Pij,group_index,time_run){
  list2env(x=params,envir=environment())
  npatch <- nrow(patch_locations)
  ngroups <- nrow(group_index)
  
  # melt into a datatable with columns patch, timestep, alpha, theta, p, popsize
  # add columns for param values at each time/patch/alpha/theta/p combo, scaled by the population size that has that combo
  sim_melt <- as.table(Pij)
  dimnames(sim_melt) <- list(1:npatch,1:ngroups,1:nsteps)
  sim_melt <- as.data.table(sim_melt)
  setnames(sim_melt, c("patch","group","t","popsize"))
  
  sim_melt[,`:=`(patch = as.numeric(patch),
                 group = as.numeric(group),
                 t = as.numeric(t))]
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
                    by=t]
  
  # then divide by total popsize (denominator of the mean)
  by_t[,`:=`(alpha = alpha/popsize,
             theta = theta/popsize,
             p = p/popsize)]
  
  return(list(sim_melt=sim_melt,
              by_t=by_t,
              patch_locations=patch_locations,
              time_run=time_run,
              K=K,
              b=b))
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
