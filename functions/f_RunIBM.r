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

  pop <- data.frame(t=1,
                    origin_site=NA, 
                    alpha=alpha_start, # note that these are the INDICES of the parameters in v_alphas and v_thetas, not the actual parameter values
                    theta=theta_start,
                    p=p_start,
                    dest_site=start_sites) # start at carrying capacity everywhere
  
  ########## Simulation ##########
  for(t_step in 1:nsteps){
    adults <- filter(pop,t==t_step)
    
    # reproduction
    larvae <- adults[FALSE,]
    for(i_adult in 1:nrow(adults)){
      b_i <- patch_locations$b_i[adults$dest_site[i_adult]] # get the birth rate in that adult's patch
      larvae <- rbind(larvae,bind_rows(replicate(b_i, adults[i_adult,], simplify = FALSE))) 
    }
    larvae <- mutate(larvae, origin_site=dest_site) %>%
      mutate(dest_site=NA) %>%
      mutate(t=t_step+1)
    
    # dispersal
    for(i in 1:nrow(larvae)){
      b_i <- patch_locations$b_i[larvae$origin_site[i]] # habitat quality (b) in larva's origin patch. Maybe at this point it makes sense to just add this as a column.
      i_alpha <- larvae$alpha[i] # current larva's alpha index
      alpha_plastic <- case_when( # current larva's effective alpha index (given plasticity)
        b_i==b_neutral ~ i_alpha,
        b_i==b_bad ~ oob_squish(i_alpha+round(larvae$p[i]),c(1,length(v_alphas))), # we'll want something more sophisticated than round(v_p[i_p]) eventually
        b_i==b_good ~ oob_squish(i_alpha-round(larvae$p[i]),c(1,length(v_alphas)))
      )
      dests <- conn_matrices[alpha_plastic, larvae[i,]$theta, larvae[i,]$origin_site, ]
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
    # larvae <- group_by(larvae,dest_site) %>%
    #   slice_sample(n=K)
    for(site in unique(larvae$dest_site)){
      survivors <- filter(larvae,dest_site==site) %>%
        slice_sample(n=filter(patch_locations,id==site)$K_i)
      pop <- rbind(pop,survivors)
    }
    
    # cleanup for next timestep
    # pop <- rbind(pop,larvae)
    
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
