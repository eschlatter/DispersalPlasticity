f_RunIBM <- function(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,heatmap_plot_int){
  starttime <- proc.time()
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
  start_sites=rep(patch_locations$id,K)
  # start_sites=rep(patch_locations$id[round(2*npatch/5),round(3*npatch/5)],K)
  # start_sites=rep(patch_locations$id[1:5],K)
  pop <- data.frame(t=1,
                    origin_site=NA, 
                    alpha=alpha_start, # note that these are the INDICES of the parameters in v_alphas and v_thetas, not the actual parameter values
                    theta=theta_start, 
                    dest_site=start_sites) # start at carrying capacity everywhere
  
  ########## Simulation ##########
  for(t_step in 1:nsteps){
    adults <- filter(pop,t==t_step)
    
    # reproduction
    larvae <- do.call("rbind", replicate(b, adults, simplify = FALSE)) %>%
      mutate(origin_site=dest_site) %>%
      mutate(dest_site=NA) %>%
      mutate(t=t_step+1)
    
    # dispersal
    for(i in 1:nrow(larvae)){
      dests <- conn_matrices[larvae[i,]$alpha, larvae[i,]$theta, larvae[i,]$origin_site, ]
      larvae[i,]$dest_site <- sample(1:npatch,size=1,prob=dests)
    } # i
    
    # mutation
    alpha_adds <- sample(c(0,1,-1),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
    theta_adds <- sample(c(0,1,-1),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
    larvae$alpha <- oob_squish(larvae$alpha + alpha_adds, c(1,length(v_alphas))) # bound each alpha index so it's >=1 and <= length(v_alphas)
    larvae$theta <- oob_squish(larvae$theta + theta_adds, c(1,length(v_thetas)))
    
    # competition
    # keep only K larvae (at random) at each destination site
    # (if there are fewer than K larvae at a site, they'll all be kept)
    larvae <- group_by(larvae,dest_site) %>%
      slice_sample(n=K)

    # cleanup for next timestep
    pop <- rbind(pop,larvae)
    
    if(t_step %% round(nsteps/10) == 0) print(t_step)
    
    # if(t_step%%heatmap_plot_int==0){
    #   f_PlotHeatmapsIBM(larvae,patch_locations,t_step)
    #   Sys.sleep(2)
    # }
    
  } # t_step
  
  # add values (not just indices) of parameters to dataframe
  pop <- mutate(pop,alpha_value=v_alphas[alpha],theta_value=v_thetas[theta])
  
  ########## Process data for plotting ##########
  by_t <- summarize(group_by(pop,t),alpha=mean(alpha_value),theta=mean(theta_value),popsize=n())
  
  ########## Output ##########
  time_run <- proc.time()-starttime
  return(list(pop,by_t, patch_locations, time_run))
}
