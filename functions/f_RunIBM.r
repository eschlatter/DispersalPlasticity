f_RunIBM <- function(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad=1,b_neutral=3,b_good=6,K,heatmap_plot_int){
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
  patch_locations$K_i <- as.vector(K)
  patch_locations$b_i <- as.vector(b)
  
  ########## Data structure to describe population ##########
  start_sites=rep(patch_locations$id,patch_locations$K_i)
  
  pop <- vector("list",length=nsteps) # list. List element t is the adults at timestep t.
  
  # initialize first generation of adults
  adults <- data.frame(t=1,
                       origin_site=NA, 
                       alpha=alpha_start, # note that these are the INDICES of the parameters in v_alphas and v_thetas, not the actual parameter values
                       theta=theta_start,
                       p=p_start,
                       dest_site=start_sites) # start at carrying capacity everywhere
  pop[[1]] <- adults
  
  ########## Simulation ##########
  for(t_step in 1:nsteps){
    # reproduction
    # create a dataframe with the same columns as adults and sum(K*b) number of rows, all empty
    larvae <- adults[FALSE,]
    larvae[1:sum(K*b),] <- NA
    
    on_row_larvae <- 1 #first row of larvae available for filling
    for(i_adult in 1:nrow(adults)){
      b_i <- patch_locations$b_i[adults$dest_site[i_adult]] # get the birth rate in that adult's patch
      larvae[on_row_larvae:(on_row_larvae+b_i-1),] <- adults[i_adult,]
      on_row_larvae <- on_row_larvae+b_i
    }
    larvae <- mutate(larvae, origin_site=dest_site) %>%
      mutate(dest_site=NA) %>%
      mutate(t=t_step+1) ### if we were going to add b_i to the dataframe, this would be the place to do it
    
    # dispersal
    for(i in 1:nrow(larvae)){
      b_i <- patch_locations$b_i[larvae$origin_site[i]] # habitat quality (b) in larva's origin patch. Maybe at this point it makes sense to just add this as a column.
      eff_params <- f_plasticity(b_i,larvae$p[i],larvae$alpha[i],larvae$theta[i],b_bad,b_neutral,b_good,length(v_alphas),length(v_thetas))
      dests <- conn_matrices[eff_params[[1]], eff_params[[2]], larvae$origin_site[i], ]
      larvae$dest_site[i] <- sample(1:npatch,size=1,prob=dests)
    } # i
    
    # mutation
    alpha_adds <- sample(c(0,1,-1),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
    theta_adds <- sample(c(0,1,-1),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
    larvae$alpha <- oob_squish(larvae$alpha + alpha_adds, c(1,length(v_alphas))) # bound each alpha index so it's >=1 and <= length(v_alphas)
    larvae$theta <- oob_squish(larvae$theta + theta_adds, c(1,length(v_thetas)))
    
    # competition
    # keep only K larvae (at random) at each destination site
    # (if there are fewer than K larvae at a site, they'll all be kept)
    # larvae <- group_by(larvae,dest_site) %>%
    #   slice_sample(n=K)
    adults <- larvae[FALSE,]
    adults[1:sum(K),] <- NA
    on_row_survivors <- 1
    for(site in unique(larvae$dest_site)){
      site_survivors <- filter(larvae,dest_site==site) %>%
        slice_sample(n=filter(patch_locations,id==site)$K_i)
      adults[on_row_survivors:(on_row_survivors+nrow(site_survivors)-1),] <- site_survivors
      on_row_survivors <- on_row_survivors+nrow(site_survivors)
    }
    
    # cleanup for next timestep
    pop[[t_step]] <- adults
    
    if(t_step %% round(nsteps/10) == 0) print(t_step)
    
    # if(t_step%%heatmap_plot_int==0){
    #   f_PlotHeatmapsIBM(larvae,patch_locations,t_step)
    #   Sys.sleep(2)
    # }
    
  } # t_step
  
  # combine timesteps into one dataframe
  pop_df <- do.call('rbind',pop)
  
  # add values (not just indices) of parameters to dataframe
  pop_df <- mutate(pop_df,alpha_value=v_alphas[alpha],theta_value=v_thetas[theta])
  
  ########## Process data for plotting ##########
  by_t <- summarize(group_by(pop_df,t),alpha=mean(alpha_value),theta=mean(theta_value),popsize=n())
  
  ########## Output ##########
  time_run <- proc.time()-starttime
  return(list(pop_df,by_t, patch_locations, time_run))
}