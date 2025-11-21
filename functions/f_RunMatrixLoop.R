f_RunMatrixLoop <- function(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL,seed=NULL){
  starttime <- proc.time()
  
  if(!is.null(seed)) set.seed(seed)
  
  # -------------------------------------------------------------------
  # Data structures to describe space and dispersal
  # -------------------------------------------------------------------
  # see utility functions/f_MakeHabitat for details on what's in each object
  hab <- f_MakeHabitat(nx,ny,v_alphas,v_thetas,patch_locations)
  patch_locations <- hab$patch_locations
  patch_dists <- hab$patch_dists
  patch_angles <- hab$patch_angles
  npatch <- hab$npatch
  rm(hab)
  if(!"K_i" %in% colnames(patch_locations)) patch_locations$K_i <- as.vector(K)
  patch_locations$b_i <- as.vector(b)
  
  # -------------------------------------------------------------------
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
  Pij[,start_grp,1] <- K
  
  ## 3. mutation_destinations: for each parameter groups, what parameter groups can a single mutation reach?
  ## dimensions 1 = ngroups, 2 = number of types of mutation events (including no mutation)
  
  # first, list the possible mutations to each parameter
  # (each index represents a mutation event; only one parameter can change per mutation event)
  alpha_adds=c(0,1,-1,0,0)
  theta_adds=c(0,0,0,1,-1)
  p_adds=c(0,0,0,0,0)
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
  
  ## Initialize temporary data structures
  temp_pop <- matrix(0,nrow=npatch,ncol=ngroups) # hold intermediate population values for this timestep, before competition
  to_patch <- numeric(npatch) # hold numbers of immigrants to each patch during dispersal
  
  ## Precompute connectivity matrices
  ## Incorporates plasticity and K_i -- so needs to change if K_i changes
  all_conn_mats <- vector("list", ngroups)
  for (g in 1:ngroups) {
    v <- group_index[g,]
    # compute effective parameters for each patch with plasticity (once per group)
    eff_params <- f_plasticityK(patch_locations$K_i, 
                                v_p[v$p], 
                                v$alpha, 
                                v$theta,
                                n_alpha = length(v_alphas),
                                n_theta = length(v_thetas))
    # build matrix
    conn_mat <- f_GetConnectivityMatrix_vectorized(v_alphas[eff_params$alpha_plastic],
                                                   v_thetas[eff_params$theta_plastic],patch_dists,patch_angles)
    all_conn_mats[[g]] <- conn_mat
  }
  
  
  # -------------------------------------------------------------------
  # Simulate
  # -------------------------------------------------------------------
  for(t in 2:nsteps){
    temp_pop[,] <- 0 # reset temp_pop
    
    ################## Reproduction and Dispersal and Mutation ##################
    
    for(g in 1:ngroups){
      # get parameter values for that parameter group
      v <- group_index[g,] 
      
      # get population of each patch for that parameter group
      patch_pops <- Pij[,g,t-1]
      if(sum(patch_pops)>0){
        # get effective kernel parameters for each origin patch, given plasticity
        eff_params <- f_plasticityK(patch_locations$K_i,v_p[v$p],v$alpha,v$theta,n_alpha=length(v_alphas),n_theta = length(v_thetas))
        
        # calculate the connectivity matrix among patches, given the group parameter values and patch-level K's
        # (and accounting for the patch population x per capita output b_i from each patch)
        conn_mat <- all_conn_mats[[g]]
        to_patch <- patch_locations$b_i*(conn_mat %*% patch_pops) # vector of contribution of the population of this group to each patch
        to_patch_reverse <- patch_pops %*% conn_mat
        
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
        survivors=sample(x=ngroups, # there are this many cells (i.e., combos of parameter values) for the patch
                         size=min(patch_abunds[i_patch],patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                         prob = temp_pop[i_patch,], # probability of each cell being chosen depends on its current abundance
                         replace=TRUE) |> # same cell can be chosen by multiple individuals
          tabulate() # make into a vector with the number of individuals that chose each cell
        survivors <- c(survivors,rep(0, ngroups - length(survivors))) # add zeros to the end for cells past the last one chosen
        Pij[i_patch,,t] <- survivors
      }
      else Pij[i_patch,,t] <- 0L
    }
    
    if(t %% max(1,round(nsteps/10)) == 0) print(t)
  }
  
  # -------------------------------------------------------------------
  # Process data for output
  # -------------------------------------------------------------------
  
  # melt into a dataframe with columns patch, timestep, alpha, theta, p, popsize
  # add columns for param values at each time/patch/alpha/theta/p combo, scaled by the population size that has that combo
  group_index <- mutate(group_index,group=1:nrow(group_index))
  sim_melt <- melt(Pij, varnames=c("patch","group",'t'),value.name="popsize") %>%
    left_join(group_index)
  sim_melt <- sim_melt[,c("patch", "t", "alpha", "theta", "p", "popsize")] %>%
    mutate(alpha_value=v_alphas[alpha],theta_value=v_thetas[theta],p_value=v_p[p],t=as.numeric(t))
  
  # mean param values at each timepoint
  # first take the sum across all cells of the param*popsize (numerator of the mean)
  by_t <- summarize(group_by(sim_melt,t),
                    alpha=sum(alpha_value*popsize),theta=sum(theta_value*popsize),popsize=sum(popsize))
  # then divide by total popsize (denominator of the mean)
  by_t <- mutate(by_t, alpha=alpha/popsize, theta=theta/popsize)
  
  time_run <- proc.time()-starttime
  
  # -------------------------------------------------------------------
  # Output
  # -------------------------------------------------------------------
  return(list(sim_melt=sim_melt,
              by_t=by_t,
              patch_locations=patch_locations,
              time_run=time_run,
              K=K,
              b=b))
}
