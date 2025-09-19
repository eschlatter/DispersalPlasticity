f_RunMatrixSim <- function(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b){
  
  ########## Data structures to describe space and dispersal ##########
  hab <- f_MakeHabitat(nx,ny,v_alphas,v_thetas)
  patch_locations <- hab$patch_locations
  patch_map <- hab$patch_map
  patch_dists <- hab$patch_dists
  patch_angles <- hab$patch_angles
  conn_matrices <- hab$conn_matrices
  npatch <- hab$npatch
  rm(hab)
  
  ########## Data structure to describe population ##########
  
  # dimensions: npatch x length(v_alphas) x length(v_thetas) x length(v_p) x nsteps
  # 1: location (patch id)
  # 2: alpha (kernel shape parameter)
  # 3: theta (kernel scale parameter)
  # 4: p (plasticity parameter)
  # 5: timestep
  sim_array <- array(0,dim=c(npatch,length(v_alphas),length(v_thetas),length(v_p),nsteps))
  
  ########## Simulation ##########
  
  # initialize starting population
  # currently 1 individual (whatever that means) per patch
  # everybody starts with the same kernel and p
  sim_array[,alpha_start,theta_start,p_start,1] <- 1
  
  for(t in 2:nsteps){
    # reproduction and dispersal and mutation (each patch contributes to other patches)
    for(i_alpha in 1:length(v_alphas)){
      for(i_theta in 1:length(v_thetas)){
        cm <- conn_matrices[i_alpha,i_theta,,]
        for(i_patch in 1:npatch){
          cell_popsize <- sim_array[i_patch,i_alpha,i_theta,p_start,t-1]
          
          # no mutation
          sim_array[,i_alpha,i_theta,p_start,t] <- b*cell_popsize*(1-mu)*cm[i_patch,] + sim_array[,i_alpha,i_theta,p_start,t]
          
          # mutation
          # "absorbing boundaries" at the edge of allowable kernel parameters
          # surely there's a more efficient way to do this, but we'll go with this for now
          if(i_alpha!=length(v_alphas)){
            sim_array[,i_alpha+1,i_theta,p_start,t] <- b*cell_popsize*(mu/4)*cm[i_patch] + sim_array[,i_alpha+1,i_theta,p_start,t]
          }
          if(i_alpha!=1){
            sim_array[,i_alpha-1,i_theta,p_start,t] <- b*cell_popsize*(mu/4)*cm[i_patch] + sim_array[,i_alpha-1,i_theta,p_start,t]
          }
          # mutation to theta
          if(i_theta!=length(v_thetas)){
            sim_array[,i_alpha,i_theta+1,p_start,t] <- b*cell_popsize*(mu/4)*cm[i_patch] + sim_array[,i_alpha,i_theta+1,p_start,t]
          }
          if(i_theta!=1){
            sim_array[,i_alpha,i_theta-1,p_start,t] <- b*cell_popsize*(mu/4)*cm[i_patch] + sim_array[,i_alpha,i_theta-1,p_start,t]
          }
        } # i_patch
      } # i_theta
    } # i_alpha
    
    # competition
    # want to cap the population of each patch, probably at 1 for now.
    # so, if the patch has population greater than 1, scale the value in each box by 1/(sum of all boxes for that patch)
    pop_by_patch <- apply(sim_array[,,,,t],1,sum)
    pop_by_patch <- pmax(1,pop_by_patch) # scale by 1 (leave it alone) if population of a patch is less than 1
    sim_array[,,,,t] <- sweep(sim_array[,,,,t],MARGIN=1,FUN='/',STATS=pop_by_patch) # scale by patch population
  } # t
  
  ########## Make some plots ##########
  
  # process data for plotting
  
  # melt into a dataframe with columns patch, timestep, alpha, theta, p, popsize
  dimnames(sim_array) <- list(patch_locations$id,
                              v_alphas,
                              v_thetas,
                              v_p,
                              1:nsteps)
  sim_melt <- array2DF(sim_array) %>% 
    mutate(across(where(is.character), as.numeric))
  rm(sim_array)
  colnames(sim_melt) <- c('patch','alpha','theta','p','t','popsize')
  
  #param values at each time/patch/alpha/theta/p combo, scaled by the population size
  sim_melt <- mutate(sim_melt,
                     alpha_scale=alpha*popsize,
                     theta_scale=theta*popsize,
                     p_scale=p*popsize)
  
  # mean param values within each patch at each timepoint
  sim_melt_params <- summarize(group_by(sim_melt,t,patch),
                               alpha_mean=sum(alpha_scale),
                               theta_mean=sum(theta_scale),
                               popsize=sum(popsize))
  
  # mean param values, averaged over patches
  by_t <- summarize(group_by(sim_melt_params,t),
                    alpha=mean(alpha_mean),
                    theta=mean(theta_mean),
                    popsize=sum(popsize))
  return(by_t)
  
}