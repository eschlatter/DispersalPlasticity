############## after-the-fact heatmap function for IBM ###############
# inputs:
#   pop: population dataframe for the current timestep
#   patch_locations for mapmaking
#   plot_int to specify time intervals to plot
f_PlotAllHeatmapsIBM <- function(pop,patch_locations,plot_int=NA){
  if(is.na(plot_int)) plot_int <- round(max(pop$t)/10) # set the plotting interval, unless specified
  plot_ints <- seq(from=plot_int,to=max(pop$t),by=plot_int)
  
  by_patch <- pop %>%
    filter(t %in% plot_ints) %>%
    group_by(dest_site,t) %>%
    # variance of alpha and theta doesn't work at the moment, b/c just one individual per patch
    summarize(alpha_m=mean(alpha_value),alpha_v=var(alpha_value),theta_m=mean(theta_value),theta_v=var(theta_value),popsize=n()) %>%
    left_join(patch_locations,by=c("dest_site"="id"))
  
  for(t_i in plot_ints){
    by_patch_i <- filter(by_patch,t==t_i)
    
    plot_alpha <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=alpha_m))+
      geom_tile()+
      scale_x_continuous(breaks=0:10)+
      scale_y_continuous(breaks=0:10)+
      labs(x='x',y='y')+
      coord_fixed()
    
    plot_alpha_v <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=alpha_v))+
      geom_tile()+
      scale_x_continuous(breaks=0:10)+
      scale_y_continuous(breaks=0:10)+
      labs(x='x',y='y')+
      coord_fixed()
    
    plot_theta <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=theta_m))+
      geom_tile()+
      scale_x_continuous(breaks=0:10)+
      scale_y_continuous(breaks=0:10)+
      labs(x='x',y='y')+
      coord_fixed()
    
    plot_theta_v <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=theta_v))+
      geom_tile()+
      scale_x_continuous(breaks=0:10)+
      scale_y_continuous(breaks=0:10)+
      labs(x='x',y='y')+
      coord_fixed()
    
    plot_abund <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=popsize))+
      geom_tile()+
      scale_x_continuous(breaks=0:10)+
      scale_y_continuous(breaks=0:10)+
      labs(x='x',y='y')+
      coord_fixed()
    
    grid.arrange(plot_alpha, plot_alpha_v, plot_theta, plot_theta_v, plot_abund,ncol=2,top=paste0('t = ',t_i))
  }  
}

############## dynamic heatmap function for IBM ###############
# inputs:
#   pop_t: population dataframe for the current timestep
#   patch_locations for mapmaking
#   t current timestep
f_PlotHeatmapsIBM <- function(pop_t,patch_locations,t){
  by_patch <- pop_t %>%
    group_by(dest_site) %>%
    summarize(alpha=mean(alpha_value),theta=mean(theta_value),popsize=n()) %>%
    left_join(patch_locations,by=c("dest_site"="id"))
  
  plot_alpha <- ggplot(by_patch,aes(x=x-0.5,y=y-0.5,fill=alpha))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()
  
  plot_theta <- ggplot(by_patch,aes(x=x-0.5,y=y-0.5,fill=theta))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()
  
  plot_abund <- ggplot(by_patch,aes(x=x-0.5,y=y-0.5,fill=popsize))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()
  
  grid.arrange(plot_alpha, plot_theta, plot_abund,ncol=1,top=paste0('t = ',t))
}

############## after-the-fact bubble plot function for matrix model ##################
# inputs:
#   sim_melt: matrix output melted to df
#   patch_locations for mapmaking
#   plot_int: optionally, specify timestep intervals for plots
f_PlotBubbleMatrix <- function(sim_melt,patch_locations,plot_int=NA){
  if(is.na(plot_int)) plot_int <- round(max(sim_melt$t)/10) # set the plotting interval, unless specified
  plot_ints <- seq(from=plot_int,to=max(sim_melt$t),by=plot_int)
  
  ## alpha
  # create the coordinates for each bubble
  n_alphas <- length(unique(sim_melt$alpha))
  box_nrow <- round(sqrt(n_alphas))
  box_ncol <- ceiling(n_alphas/box_nrow)
  box_x <- seq(from=1/(2*box_ncol),to=1,by=1/box_ncol)
  box_y <- seq(from=1/(2*box_nrow),to=1,by=1/box_nrow)
  box_spots <- expand_grid(x_add=box_x,y_add=box_y)[1:n_alphas,] %>%
    mutate(alpha=round(v_alphas,digits=10))
  # process data
  alpha_df_plot <- group_by(sim_melt,patch,alpha,t) %>%
    mutate(alpha=round(alpha,digits=10)) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of theta
    left_join(box_spots,by='alpha') %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  ## theta
  # create the coordinates for each bubble
  n_thetas <- length(unique(sim_melt$theta))
  box_nrow <- round(sqrt(n_thetas))
  box_ncol <- ceiling(n_thetas/box_nrow)
  box_x <- seq(from=1/(2*box_ncol),to=1,by=1/box_ncol)
  box_y <- seq(from=1/(2*box_nrow),to=1,by=1/box_nrow)
  box_spots <- expand_grid(x_add=box_x,y_add=box_y)[1:n_thetas,] %>%
    mutate(theta=round(v_thetas,digits=10)) # need to round so the numbers match up when we join with theta_df_plot
  # process data
  theta_df_plot <- group_by(sim_melt,patch,theta,t) %>%
    mutate(theta=round(theta,digits=10)) %>% # need to round so the numbers match up when we join with box_spots
    summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of alpha
    left_join(box_spots,by='theta') %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  for(t_i in plot_ints){
    plot_alpha <- ggplot(filter(alpha_df_plot,t==t_i),aes(x=x+x_add-1,y=y+y_add-1)) +
      geom_point(aes(color=alpha,size=popsize),pch=19) + 
      geom_hline(yintercept=seq(from=0,to=ny,by=max(round(ny/10),1)))+
      geom_vline(xintercept=seq(from=0,to=nx,by=max(round(nx/10),1)))+
      labs(x='x',y='y')+
      coord_fixed()
    
    ggplot(filter(alpha_df_plot,t==t_i,alpha==v_alphas[1]),aes(x=x+x_add-1,y=y+y_add-1)) +
      geom_tile(aes(fill=popsize)) + 
      labs(x='x',y='y')+
      coord_fixed()
    
    plot_theta <- ggplot(filter(theta_df_plot,t==t_i),aes(x=x+x_add-1,y=y+y_add-1)) +
      geom_point(aes(color=theta,size=popsize),pch=19) + 
      geom_hline(yintercept=seq(from=0,to=ny,by=max(round(ny/10),1)))+
      geom_vline(xintercept=seq(from=0,to=nx,by=max(round(nx/10),1)))+
      labs(x='x',y='y')+
      coord_fixed()
    grid.arrange(plot_alpha, plot_theta,ncol=1,top=paste0('t = ',t_i))
  }
}

############## after-the-fact heatmap function for matrix model ##################
# inputs:
#   sim_array_t: just the portion of sim_array from timestep t
#   patch_locations for mapmaking
f_PlotAllHeatmaps <- function(sim_melt,patch_locations,plot_int=NA){
  if(is.na(plot_int)) plot_int <- round(max(sim_melt$t)/10) # set the plotting interval, unless specified
  plot_ints <- seq(from=plot_int,to=max(sim_melt$t),by=plot_int)
  
  # process data
  alpha_by_patch <- group_by(sim_melt,patch,alpha,t) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>%  # add up what's in the boxes with all values of theta
    group_by(patch,t) %>%
    summarize(alpha_m=sum(alpha*popsize)/sum(popsize), # at each patch, find the mean and variance of alpha values
              alpha_v=sum(popsize*(alpha-alpha_m)^2)/sum(popsize),
              .groups='drop') %>%   
    left_join(patch_locations,by=c("patch" = "id"))
  
  theta_by_patch <- group_by(sim_melt,patch,theta,t) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of alpha
    group_by(patch,t) %>%
    summarize(theta_m=sum(theta*popsize)/sum(popsize),
              theta_v=sum(popsize*(theta-theta_m)^2)/sum(popsize),,
              .groups='drop') %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  pop_by_patch <- group_by(sim_melt,patch,t) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  ## make plots
  for(t_i in plot_ints){
    plot_alpha <- filter(alpha_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=alpha_m))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed() 
    
    plot_alpha_v <- filter(alpha_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=alpha_v))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed() 
    
    plot_theta <- filter(theta_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=theta_m))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()
    
    plot_theta_v <- filter(theta_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=theta_v))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()
    
    plot_abund <- filter(pop_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=popsize))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()
    
    grid.arrange(plot_alpha, plot_alpha_v, plot_theta, plot_theta_v, plot_abund,ncol=2,top=paste0('t = ',t_i))
  }
}

############## dynamic heatmap function for matrix model ##################
# inputs:
#   sim_array_t: just the portion of sim_array from timestep t
#   patch_locations for mapmaking
f_PlotHeatmaps <- function(sim_array_t,patch_locations,t){
  ### data processing
  # melt into a dataframe with columns patch, timestep, alpha, theta, p, popsize
  dimnames(sim_array_t) <- list(patch_locations$id,
                                v_alphas,
                                v_thetas,
                                v_p)
  sim_melt <- array2DF(sim_array_t)
  sim_melt <- mutate_all(sim_melt, as.numeric)
  colnames(sim_melt) <- c('patch','alpha','theta','p','popsize')
  
  ### make maps
  ## map of alpha values
  alpha_by_patch <- group_by(sim_melt,patch,alpha) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>%  # add up what's in the boxes with all values of theta
    group_by(patch) %>%
    summarize(alpha=sum(alpha*popsize)/sum(popsize)) %>%   # at each patch, find the mean value of alpha
    left_join(patch_locations,by=c("patch" = "id"))
  
  plot_alpha <- ggplot(alpha_by_patch,aes(x=x-0.5,y=y-0.5,fill=alpha))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()
  
  ## map of theta values
  theta_by_patch <- group_by(sim_melt,patch,theta) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of alpha
    group_by(patch) %>%
    summarize(theta=sum(theta*popsize)/sum(popsize)) %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  plot_theta <- ggplot(theta_by_patch,aes(x=x-0.5,y=y-0.5,fill=theta))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()
  
  ## map of population size
  pop_by_patch <- group_by(sim_melt,patch) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  plot_abund <- ggplot(pop_by_patch,aes(x=x-0.5,y=y-0.5,fill=popsize))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()
  
  grid.arrange(plot_alpha, plot_theta, plot_abund,ncol=1,top=paste0('t = ',t))
}

############## quick diagnostic plotting function ##################

f_PlotOutput <- function(by_t,kern_timesteps,kern_xlim=25){
  p0 <- ggplot(by_t,aes(x=t))+
    geom_line(aes(y=alpha,color='alpha'))+
    geom_line(aes(y=theta,color='theta'))+
    labs(title='kernel parameters',y='value')+
    theme_minimal()+
    theme(legend.position = 'top')
  
  p1 <- ggplot(by_t,aes(x=alpha,y=theta))+
    geom_path(alpha=0.75,lwd=0.25)+
    theme_minimal()+
    geom_point(data=last(by_t),aes(x=alpha,y=theta),color='red')+
    labs(title='kernel parameters')
  
  p2 <- ggplot()+
    xlim(0,kern_xlim)+
    geom_function(fun=dgamma, args=list(shape=first(by_t$alpha),scale=first(by_t$theta)),aes(lty='first'))+
    lapply(kern_timesteps, function(i){geom_function(fun=dgamma,args=list(shape=by_t$alpha[i],scale=by_t$theta[i]),alpha=0.5,color='darkgray')})+
    geom_function(fun=dgamma,args=list(shape=median(by_t$alpha[kern_timesteps]),scale=median(by_t$theta[kern_timesteps])),color='black',lwd=0.75,aes(lty='last'))+
    scale_linetype_manual(values=c('first' = 'dashed',
                                   'last' = 'solid'))+
    theme_minimal()+
    theme(legend.position='top')+
    labs(x='distance',y='density',title='Kernel')
  
  p3 <- ggplot(by_t,aes(x=t,y=popsize))+
    geom_line()+
    theme_minimal()+
    labs(title='population size')
  
  grid.arrange(p0,p1,p2,p3,nrow=1)
}