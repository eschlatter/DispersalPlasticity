f_RunMatrixSim <- function(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,heatmap_plot_int=NA,sleep_int=0, competition_method='sample'){
  starttime <- proc.time()
  if(!(competition_method %in% c('sample','rnorm'))) stop("competition method incorrectly specified")
  
  ########## Data structures to describe space and dispersal ##########
  hab <- f_MakeHabitat(nx,ny,v_alphas,v_thetas)
  patch_locations <- hab$patch_locations
  patch_map <- hab$patch_map
  patch_dists <- hab$patch_dists
  patch_angles <- hab$patch_angles
  conn_matrices <- hab$conn_matrices
  npatch <- hab$npatch
  rm(hab)
  patch_locations$K_i <- as.vector(K)
  patch_locations$b_i <- as.vector(b)
  
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
  sim_array[,alpha_start,theta_start,p_start,1] <- patch_locations$K_i
  
  for(t in 2:nsteps){
    # reproduction and dispersal and mutation (each patch contributes to other patches)
    for(i_alpha in 1:length(v_alphas)){
      for(i_theta in 1:length(v_thetas)){
        cm <- conn_matrices[i_alpha,i_theta,,]
        for(i_patch in 1:npatch){
          cell_popsize <- sim_array[i_patch,i_alpha,i_theta,p_start,t-1]
          b_i <- patch_locations$b_i[i_patch] # reproductive rate in the patch
          # no mutation
          sim_array[,i_alpha,i_theta,p_start,t] <- b_i*cell_popsize*(1-mu)*(cm[i_patch,]) + sim_array[,i_alpha,i_theta,p_start,t]
          
          # mutation
          # boundaries at the edge of allowable kernel params: collecting
          # (e.g., at minimum value of alpha, if mutation would lead to lower alpha, it just keeps the minimum)
          # surely there's a more efficient way to do this, but we'll go with this for now
          sim_array[,min(i_alpha+1,length(v_alphas)),i_theta,p_start,t] <- b_i*cell_popsize*(mu/4)*(cm[i_patch,]) + sim_array[,min(i_alpha+1,length(v_alphas)),i_theta,p_start,t]
          sim_array[,max(1,i_alpha-1),i_theta,p_start,t] <- b_i*cell_popsize*(mu/4)*(cm[i_patch,]) + sim_array[,max(1,i_alpha-1),i_theta,p_start,t]
          sim_array[,i_alpha,min(i_theta+1,length(v_thetas)),p_start,t] <- b_i*cell_popsize*(mu/4)*(cm[i_patch,]) + sim_array[,i_alpha,min(i_theta+1,length(v_thetas)),p_start,t]
          sim_array[,i_alpha,max(1,i_theta-1),p_start,t] <- b_i*cell_popsize*(mu/4)*(cm[i_patch,]) + sim_array[,i_alpha,max(1,i_theta-1),p_start,t]
          
        } # i_patch
      } # i_theta
    } # i_alpha
    
    ## Competition
    # Method 1: if a patch has population greater than K, sample K individuals and distribute them among cells in that patch
    # (with probability according to the current abundance of each cell)
    if(competition_method=='sample'){
      pop_by_patch <- apply(sim_array[,,,,t],1,sum)
      for(i_patch in 1:npatch){
          survivors <- sample(x=length(v_alphas)*length(v_thetas)*length(v_p),
                              size = min(pop_by_patch[i_patch],patch_locations$K_i[i_patch]),
                              prob = sim_array[i_patch,,,,t],
                              replace=TRUE)
          survivors <- as.data.frame(table(survivors)) %>% 
            mutate(cell=as.numeric(as.character(survivors)))
          
          new <- array(0, dim=dim(sim_array[i_patch,,,,t,drop=F]))
          new[survivors$cell] <- survivors$Freq
          sim_array[i_patch,,,,t] <- new
      } # i_patch
    }
    
    # Method 2: if a patch has population greater than K, scale the value in each cell by rnorm(mean = K/(sum of all boxes for that patch))
    if(competition_method=='rnorm'){
      pop_by_patch <- apply(sim_array[,,,,t],1,sum)
      scale_by_patch <- ifelse(pop_by_patch>K,pop_by_patch/K,1) # what to divide the patch population by, to reduce it to carrying capacity
      # store the values to scale each cell by
      scale_by_cell <- array(dim=c(npatch,length(v_alphas),length(v_thetas),length(v_p),1))
      for(i_patch in 1:npatch){
        scale_by_cell[i_patch,,,,1] <- pmax(rnorm(n=length(v_alphas)*length(v_thetas)*length(v_p),mean=scale_by_patch[i_patch],sd=mean(b)/10),0)
      }
      # do the scaling
      sim_array[,,,,t] <- sim_array[,,,,t,drop=F]/scale_by_cell
    }
    
    
    ## output
    if(t%%heatmap_plot_int==0){
      f_PlotHeatmaps(sim_array[,,,,t],patch_locations,t)
      Sys.sleep(sleep_int)
    }
    if(t %% round(nsteps/10) == 0) print(t)
    
  } # t
  
  
  ########## Process data for plotting ##########
  
  # melt into a dataframe with columns patch, timestep, alpha, theta, p, popsize
  dimnames(sim_array) <- list(patch_locations$id,
                              v_alphas,
                              v_thetas,
                              v_p,
                              1:nsteps)
  sim_melt <- array2DF(sim_array)
  rm(sim_array) # remove this to free up memory
  sim_melt <- mutate_all(sim_melt, as.numeric)
  colnames(sim_melt) <- c('patch','alpha','theta','p','t','popsize')
  
  # param values at each time/patch/alpha/theta/p combo, scaled by the population size that has that combo
  sim_melt <- mutate(sim_melt,
                     alpha_scale=alpha*popsize,
                     theta_scale=theta*popsize,
                     p_scale=p*popsize)
  
  # mean param values at each timepoint
  # first take the sum across all cells of the param*popsize (numerator of the mean)
  by_t <- summarize(group_by(sim_melt,t),alpha=sum(alpha_scale),theta=sum(theta_scale),popsize=sum(popsize))
  # then divide by total popsize (denominator of the mean)
  by_t <- mutate(by_t, alpha=alpha/popsize, theta=theta/popsize)
  
  time_run <- proc.time()-starttime
  ########## Output ##########
  return(list(sim_melt,by_t,patch_locations, time_run,K))
}
