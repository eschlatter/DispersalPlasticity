########################################################################
# first version of this function:
# generates parameter- (and K-) specific connectivity matrices as they're needed, and then stores them
f_RunMatrixLoop <- function(params){
  starttime <- proc.time()
  numCores <- detectCores()
  list2env(x=params,envir=environment())
  
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
      if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas, numCores) # calculate this matrix, if it hasn't been used before
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
  return(list(params=params,Pij=Pij,group_index=group_index,time_run=time_run))
}

########################################################################
# version that only stores specified info
f_RunMatrixLoopLite <- function(params, keep=list("p")){
  starttime <- proc.time()
  numCores <- detectCores()
  # make each run replicable
  seed <- Sys.time()
  set.seed(seed)
  
  list2env(x=params,envir=environment())
  
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
  
  ## 2. Population objects
  ## dimensions 1 = npatch, 2 = ngroups, 3 = nsteps
  previous_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  new_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  # initialize
  start_grp <- which(group_index$alpha==alpha_start & group_index$theta==theta_start & group_index$p==p_start)
  previous_pop[,start_grp] <- patch_locations$K_i
  
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
  by_t <- data.frame(t=1:nsteps,alpha=NA,theta=NA,p=NA,popsize=NA)
  
  ## 5. Initialize output data structures
  output_list=list()
  if("p" %in% keep){
    output_list <- c(output_list,list(v_pmeans=vector(mode="numeric",length=nsteps)))
  }
  
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
    pop_by_group <- colSums(previous_pop)
    
    for(g in which(pop_by_group>0)){
      # get parameter values for that parameter group
      v <- group_index[g,]
      
      # get population of each patch for that parameter group
      patch_pops <- previous_pop[,g]
      
      # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
      # (and accounting for the patch population x per capita output b_i from each patch)
      if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas, numCores) # calculate this matrix, if it hasn't been used before
      conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
      to_patch <- patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
      
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
      new_pop[i_patch,] <- survivors
    }
    
    previous_pop <- new_pop
    
    ################## Output ##################
    if(t %% max(1,round(nsteps/10)) == 0) print(t)
    
    ## mean value of p at each timestep
    if("p" %in% keep){
      p_means_by_group <- v_p[group_index$p]
      p_mean <- sum(p_means_by_group*colSums(previous_pop))/sum(previous_pop)
      output_list$v_pmeans[t] <- p_mean
    }
    
    
    
    # # melt into a datatable with columns patch, timestep, alpha, theta, p, popsize
    # sim_melt <- as.table(new_pop[,])
    # dimnames(sim_melt) <- list(1:npatch,1:ngroups)
    # sim_melt <- as.data.table(sim_melt)
    # setnames(sim_melt, c("patch","group","popsize"))
    # sim_melt <- sim_melt[sim_melt$popsize!=0,]
    # sim_melt[,`:=`(patch = as.numeric(patch),
    #                group = as.numeric(group),
    #                t = as.numeric(t))]
    # sim_melt[,`:=`(alpha = group_index$alpha[group],
    #                theta = group_index$theta[group],
    #                p = group_index$p[group])]
    # sim_melt[,`:=`(alpha_value = v_alphas[alpha],
    #                theta_value = v_thetas[theta],
    #                p_value = v_p[p])]
    # 
    # # mean param values at each timepoint
    # # first take the sum across all cells of the param*popsize (numerator of the mean)
    # by_t_i <- sim_melt[,.(alpha=sum(alpha_value*popsize),theta=sum(theta_value*popsize),p=sum(p_value*popsize),popsize=sum(popsize)),by=t]
    # 
    # # then divide by total popsize (denominator of the mean)
    # by_t_i[,`:=`(alpha = alpha/popsize,
    #            theta = theta/popsize,
    #            p = p/popsize)]
    # 
    # by_t[t,] <- by_t_i
  }
  
  time_run <- proc.time()-starttime
  return(c(output_list,list(params=params,seed=seed)))
}


########################################################################
# parallelize it a bit
f_RunMatrixLoopParallel <- function(params){
  starttime <- proc.time()
  
  list2env(x=params,envir=environment())
  numCores <- detectCores()
  
  # Data structures to describe space and dispersal
  # -------------------------------------------------------------------
  
  # see utility functions/f_MakeHabitat for details on what's in each object
  hab <- f_MakeHabitat(nx,ny,v_alphas,v_thetas,patch_locations)
  patch_locations <- hab$patch_locations
  patch_dists <- hab$patch_dists
  patch_angles <- hab$patch_angles
  patch_map <- hab$patch_map
  npatch <- hab$npatch
  rm(hab)
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
      
      # get population of each patch for that parameter group
      patch_pops <- Pij[,g,t-1]
      
      # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
      # (and accounting for the patch population x per capita output b_i from each patch)
      if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas) # calculate this matrix, if it hasn't been used before
      conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
      to_patch <- patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
      
      # Divide up to_patch among parameter groups that are the result of mutation
      temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
      
      for(mut_group in mutation_destinations[g,-1]){ # for each of the possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
        temp_pop[,mut_group] <- (mu/4)*to_patch+temp_pop[,mut_group]
      }
    }
    
    # results <- mclapply(which(pop_by_group>0), 
    #                     function(i) f_ReprodDispMut(i,group_index,Pij[,,t-1],mutation_destinations,all_conn_mats,patch_locations,patch_dists,patch_angles,v_p,v_alphas,v_thetas,mu),
    #                     mc.cores = numCores)
    # temp_pop <- Reduce('+',results) # adds all the individual temp_pop_g's together
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    patch_abunds <- rowSums(temp_pop)
    comp_results <- mclapply(1:npatch,function(i) f_Competition(i,patch_abunds,patch_locations,temp_pop),mc.cores = numCores)
    Pij[,,t] <- do.call(rbind,comp_results)
    
    if(t %% max(1,round(nsteps/10)) == 0) print(t)
  }
  
  time_run <- proc.time()-starttime
  return(list(params=params,Pij=Pij,group_index=group_index,time_run=time_run))
}


########################################################################
# function to do competition step in parallel
f_Competition <- function(i_patch,patch_abunds,patch_locations,temp_pop){
  if(patch_abunds[i_patch]>0){
    survivors=t(rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
                          size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                          prob = temp_pop[i_patch,])) # probability of each cell being chosen depends on its current abundance)
  }
  else survivors <- vector(mode='integer',length=ngroups)
  return(survivors)
}

########################################################################
# function to do reproduction, dispersal, and mutation steps in parallel
f_ReprodDispMut <- function(g,group_index,previous_step,mutation_destinations,all_conn_mats,patch_locations,patch_dists,patch_angles,v_p,v_alphas,v_thetas,mu){
  temp_pop_g <- matrix(0,nrow=nrow(patch_locations),ncol=nrow(group_index)) # hold intermediate population values from this group
  
  # get parameter values for that parameter group
  v <- group_index[g,] 
  
  # get population of each patch for that parameter group
  patch_pops <- previous_step[,g]
  if(sum(patch_pops)>0){
    
    # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
    # (and accounting for the patch population x per capita output b_i from each patch)
    if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas) # calculate this matrix, if it hasn't been used before
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




########################################################################
# take the raw output of simulation function and process it into something plottable
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

########################################################################
# take the raw output of simulation function and process it into something plottable
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
  by_t <- sim_melt[,.(alpha=sum(alpha_value*popsize),theta=sum(theta_value*popsize),p=sum(p_value*popsize),popsize=sum(popsize)),by=t]
  
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
