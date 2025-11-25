########################################################################
# first version of this function:
# generates parameter- (and K-) specific connectivity matrices as they're needed, and then stores them
f_RunMatrixLoop <- function(params){
  starttime <- proc.time()
  
  list2env(x=params,envir=environment())
  
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
    
    for(g in 1:ngroups){
      # get parameter values for that parameter group
      v <- group_index[g,] 
      
      # get population of each patch for that parameter group
      patch_pops <- Pij[,g,t-1]
      if(sum(patch_pops)>0){
        
        # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
        # (and accounting for the patch population x per capita output b_i from each patch)
        if(is.null(all_conn_mats[[g]])) all_conn_mats[[g]] <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas) # calculate this matrix, if it hasn't been used before
        conn_mat <- all_conn_mats[[g]] # otherwise, grab it from the list
        to_patch <- patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
        
        # Divide up to_patch among parameter groups that are the result of mutation
        temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
        
        for(mut_group in mutation_destinations[g,-1]){ # for each of the 4 possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
          temp_pop[,mut_group] <- (mu/4)*to_patch+temp_pop[,mut_group]
        }
      }
    }
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    patch_abunds <- rowSums(temp_pop)
    for(i_patch in 1:npatch){
      if(patch_abunds[i_patch]>0){
        survivors=rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
                            size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                            prob = temp_pop[i_patch,]) # probability of each cell being chosen depends on its current abundance)
        Pij[i_patch,,t] <- survivors
      }
      else Pij[i_patch,,t] <- 0L
    }
    
    if(t %% max(1,round(nsteps/10)) == 0) print(t)
  }
  
  time_run <- proc.time()-starttime
  return(list(params=params,Pij=Pij,group_index=group_index,time_run=time_run))
}


########################################################################
# second version of this function:
# calculate each connectivity matrix every time it's used, and don't store it
f_RunMatrixLoop2 <- function(params){
  starttime <- proc.time()
  
  list2env(x=params,envir=environment())
  
  if(!is.null(seed)) set.seed(seed)
  
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

  # Simulate
  # -------------------------------------------------------------------
  for(t in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    ################## Disturbance ##################
    
    if(t %% 10 == 0){ # disturbances last for 10 timesteps
      # first, reset all K's from any disturbance that occurred 10 timesteps previously
      if(prod(patch_locations$K_i == as.vector(K))==0){
        patch_locations$K_i = as.vector(K)
      }
      
      # then make a new disturbance (maybe)
      if(rbinom(n=1,size=1,p=disturb_prob)==1){
        disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
        disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
        patch_locations$K_i[disturb_patches] <- 0
      }
    }
    
    ################## Reproduction and Dispersal and Mutation ##################
    
    for(g in 1:ngroups){
      # get parameter values for that parameter group
      v <- group_index[g,] 
      
      # get population of each patch for that parameter group
      patch_pops <- Pij[,g,t-1]
      if(sum(patch_pops)>0){
        
        
        # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
        # (and accounting for the patch population x per capita output b_i from each patch)
        conn_mat <- f_GetPlasticConnMat(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas) # calculate the connectivity matrix for this group
        to_patch <- patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
        
        # Divide up to_patch among parameter groups that are the result of mutation
        temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
        
        for(mut_group in mutation_destinations[g,-1]){ # for each of the 4 possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
          temp_pop[,mut_group] <- (mu/4)*to_patch+temp_pop[,mut_group]
        }
      }
    }
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    patch_abunds <- rowSums(temp_pop)
    for(i_patch in 1:npatch){
      if(patch_abunds[i_patch]>0){
        survivors=rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
                            size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                            prob = temp_pop[i_patch,]) # probability of each cell being chosen depends on its current abundance)
        Pij[i_patch,,t] <- survivors
      }
      else Pij[i_patch,,t] <- 0L
    }
    
    if(t %% max(1,round(nsteps/10)) == 0) print(t)
  }
  
  time_run <- proc.time()-starttime
  return(list(params=params,Pij=Pij,group_index=group_index,time_run=time_run))
}

########################################################################
# third version of this function:
# don't store a single connectivity matrix
f_RunMatrixLoop3 <- function(params){
  starttime <- proc.time()
  
  list2env(x=params,envir=environment())

  # Data structures to describe space and dispersal
  # -------------------------------------------------------------------
  
  # list of patch locations and IDs
  # (dimensions: npatch x 3)
  # "location" is the center of the patch
  if(is.null(patch_locations)){ # if not specified by an input map, then make one
    patch_locations <- expand.grid(y=1:ny,x=1:nx) %>% # do y first so that patches are ordered columnwise, like the way R fills a matrix
      rowid_to_column(var='id')
  }
  npatch <- nrow(patch_locations)
  if(!"K_i" %in% colnames(patch_locations)) patch_locations$K_i <- as.vector(K)
  patch_locations$b_i <- as.vector(b)
  params$patch_locations <- patch_locations
  
  # a "map" of the patch numbers, spatially arranged
  # (dimensions: nx x ny)
  patch_map <- matrix(nrow=ny,ncol=nx)
  for(i in 1:npatch){
    patch_map[patch_locations$y[i],patch_locations$x[i]] <- patch_locations$id[i]
  }
  
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
  
  # Simulate
  # -------------------------------------------------------------------
  for(t in 2:nsteps){
    temp_pop[ ] <- 0 # reset temp_pop
    
    ################## Disturbance ##################

    if(t %% 10 == 0){ # disturbances last for 10 timesteps
      # first, reset all K's from any disturbance that occurred 10 timesteps previously
      if(prod(patch_locations$K_i == as.vector(K))==0){
        patch_locations$K_i = as.vector(K)
      }

      # then make a new disturbance (maybe)
      if(rbinom(n=1,size=1,p=disturb_prob)==1){
        disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
        disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
        patch_locations$K_i[disturb_patches] <- 0
      }
    }
    
    ################## Reproduction and Dispersal and Mutation ##################
    
    for(g in 1:ngroups){
      #print(paste0("t=",t,"g=",g))
      # get parameter values for that parameter group
      v <- group_index[g,] 
      
      # get population of each patch for that parameter group
      patch_pops <- Pij[,g,t-1]
      
      # get effective parameter values, given plasticity
      eff_params <- f_plasticityK(patch_locations$K_i,v_p[v$p],v$alpha,v$theta,n_alpha = length(v_alphas),n_theta = length(v_thetas))
      
      if(sum(patch_pops)>0){
        for(origin_patch in which(patch_pops>0)){
          # pick up alpha and theta values
          this_alpha <- v_alphas[eff_params$alpha_plastic[origin_patch]]
          this_theta <- v_thetas[eff_params$theta_plastic[origin_patch]]
          
          # calculate the contribution of the population of origin_patch with parameter group g to all other patches
          patch_dists <- sqrt((patch_locations$x[origin_patch]-patch_locations$x)^2+(patch_locations$y[origin_patch]-patch_locations$y)^2)
          patch_angles <- suppressWarnings(2*asin(1/(2*patch_dists))/(2*pi))
          patch_angles[is.nan(patch_angles)] <- 1
          patch_conn_mat <- (pgamma(patch_dists+0.5,shape=this_alpha,scale=this_theta)-
                               pgamma(patch_dists-0.5,shape=this_alpha,scale=this_theta))*patch_angles
          # account for the patch population x per capita output b_i from origin_patch         
          to_patch <- patch_locations$b_i[origin_patch]*patch_pops[origin_patch]*patch_conn_mat
          
          # Divide up to_patch among parameter groups that are the result of mutation
          temp_pop[,g] <- (1-mu)*to_patch+temp_pop[,g]
          
          for(mut_group in mutation_destinations[g,-1]){ # for each of the 4 possible mutations. This doesn't need to be a for loop, but let's do some error checking first.
            temp_pop[,mut_group] <- (mu/4)*to_patch+temp_pop[,mut_group]
          } # mut_group
          
        } # origin_patch
      } # if sum(patch_pops)>0
    } # g
    
    ################## Competition ##################
    
    # sample K (or current abundance, if <K) individuals per patch and distribute them among groups of parameter values
    # (with probability according to the current abundance of each group of param values in that patch)
    patch_abunds <- rowSums(temp_pop)
    for(i_patch in 1:npatch){
      if(patch_abunds[i_patch]>0){
        survivors=rmultinom(n=1, # there are this many cells (i.e., combos of parameter values) for the patch
                            size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                            prob = temp_pop[i_patch,]) # probability of each cell being chosen depends on its current abundance)
        Pij[i_patch,,t] <- survivors
      }
      else Pij[i_patch,,t] <- 0L
    } # i_patch
    
    if(t %% max(1,round(nsteps/10)) == 0) print(t)
  } # t
  
  time_run <- proc.time()-starttime
  return(list(params=params,Pij=Pij,group_index=group_index,time_run=time_run))
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
