# takes a dataframe sim_df (could be sim_melt or a subset of it, plus K) with columns popsize, alpha, theta, p_value, K, x, y
# plots the effective kernels (considering plasticity), with the thickness of the kernel line indicating abundance of that kernel in the dataset
f_PlotEffectiveKernels <- function(sim_df,v_alphas,v_thetas,patch_locations,plot_title=NULL){
  eff_pars <- f_plasticityK_new(K=sim_df$K,p=sim_df$p,alpha=sim_df$alpha,theta=sim_df$theta,n_alpha=length(v_alphas),n_theta=length(v_thetas),Kmin=min(patch_locations$K_i),Kmax=max(patch_locations$K_i))
  sim_df <- cbind(sim_df,eff_pars) 
  
  sim_df <- sim_df %>%
    group_by(alpha_plastic,theta_plastic) %>%
    summarize(abund=sum(popsize)) %>%
    ungroup()
  
  sim_df <- mutate(sim_df,abund_scale=(abund/max(abund))+.2)
  
  ggplot()+
    xlim(1,20)+
    lapply(1:nrow(sim_df), 
           function(i){geom_function(fun=dgamma,
                                     args=list(shape=sim_df$alpha_plastic[i],scale=sim_df$theta_plastic[i]),
                                     lwd=sim_df$abund_scale[i])} )+
    theme_minimal()+
    labs(x='distance',y='density',title=plot_title)
}

############## after-the-fact heatmap function for IBM ###############
# inputs:
#   pop: population dataframe for the current timestep
#   patch_locations for mapmaking
#   plot_int to specify time intervals to plot
f_PlotAllHeatmapsIBM <- function(pop,patch_locations,plot_int=NA){
  if(prod(is.na(plot_int))==1){
    plot_int <- round(max(pop$t)/10) # set the plotting interval, unless specified
    plot_ints <- seq(from=plot_int,to=max(pop$t),by=plot_int)
  } 
  else plot_ints <- plot_int
  
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
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_alpha_v <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=alpha_v))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_theta <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=theta_m))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_theta_v <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=theta_v))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_abund <- ggplot(by_patch_i,aes(x=x-0.5,y=y-0.5,fill=popsize))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    #grid.arrange(plot_alpha, plot_alpha_v, plot_theta, plot_theta_v, plot_abund,ncol=2,top=paste0('t = ',t_i))
    grid.arrange(arrangeGrob(plot_alpha,top="alpha"), arrangeGrob(plot_theta,top="theta"), ncol=2,top=paste0('t = ',t_i))
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
    coord_fixed()+
    scale_y_reverse()
  
  
  plot_theta <- ggplot(by_patch,aes(x=x-0.5,y=y-0.5,fill=theta))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()+
    scale_y_reverse()
  
  plot_abund <- ggplot(by_patch,aes(x=x-0.5,y=y-0.5,fill=popsize))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(x='x',y='y')+
    coord_fixed()+
    scale_y_reverse()
  
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
      coord_fixed()+
      scale_y_reverse()
    
    ggplot(filter(alpha_df_plot,t==t_i,alpha==v_alphas[1]),aes(x=x+x_add-1,y=y+y_add-1)) +
      geom_tile(aes(fill=popsize)) + 
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_theta <- ggplot(filter(theta_df_plot,t==t_i),aes(x=x+x_add-1,y=y+y_add-1)) +
      geom_point(aes(color=theta,size=popsize),pch=19) + 
      geom_hline(yintercept=seq(from=0,to=ny,by=max(round(ny/10),1)))+
      geom_vline(xintercept=seq(from=0,to=nx,by=max(round(nx/10),1)))+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    grid.arrange(plot_alpha, plot_theta,ncol=1,top=paste0('t = ',t_i))
  }
}

############## after-the-fact heatmap function for matrix model ##################
# inputs:
#   replot_data: if it's been done before for this sim, can reuse the processed data from last time
#   sim_array_t: just the portion of sim_array from timestep t
#   patch_locations for mapmaking
#   plot_int: specify timesteps to plot
f_PlotAllHeatmaps <- function(sim_melt,patch_locations,plot_int=NA,replot_data=NULL){
  if(prod(is.na(plot_int))==1){
    plot_int <- round(max(sim_melt$t)/10) # set the plotting interval, unless specified
    plot_ints <- seq(from=plot_int,to=max(sim_melt$t),by=plot_int)
  } 
  else plot_ints <- plot_int
  
  if(is.null(replot_data)){
    # process data
    alpha_by_patch <- group_by(sim_melt,patch,alpha_value,t) %>%
      summarize(popsize=sum(popsize),.groups='drop') %>%  # add up what's in the boxes with all values of theta and p
      group_by(patch,t) %>%
      summarize(alpha_m=sum(alpha_value*popsize)/sum(popsize), # at each patch, find the mean and variance of alpha values
                alpha_v=sum(popsize*(alpha_value-alpha_m)^2)/sum(popsize),
                .groups='drop') %>%   
      left_join(patch_locations,by=c("patch" = "id"))
    
    theta_by_patch <- group_by(sim_melt,patch,theta_value,t) %>%
      summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of alpha and p
      group_by(patch,t) %>%
      summarize(theta_m=sum(theta_value*popsize)/sum(popsize),
                theta_v=sum(popsize*(theta_value-theta_m)^2)/sum(popsize),
                .groups='drop') %>%
      left_join(patch_locations,by=c("patch" = "id"))
    
    p_by_patch <- group_by(sim_melt,patch,p_value,t) %>%
      summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of alpha and theta
      group_by(patch,t) %>%
      summarize(p_m=sum(p_value*popsize)/sum(popsize),
                p_v=sum(popsize*(p_value-p_m)^2)/sum(popsize),
                .groups='drop') %>%
      left_join(patch_locations,by=c("patch" = "id"))
    
    pop_by_patch <- group_by(sim_melt,patch,t) %>%
      summarize(popsize=sum(popsize),.groups='drop') %>%
      left_join(patch_locations,by=c("patch" = "id"))
  }
  else list2env(replot_data,envir=environment())
  
  ## make plots
  for(t_i in plot_ints){
    plot_alpha <- filter(alpha_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=alpha_m))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_alpha_v <- filter(alpha_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=alpha_v))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_theta <- filter(theta_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=theta_m))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_theta_v <- filter(theta_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=theta_v))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_p <- filter(p_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=p_m))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_p_v <- filter(p_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=p_v))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    plot_abund <- filter(pop_by_patch,t==t_i) %>%
      ggplot(aes(x=x-0.5,y=y-0.5,fill=popsize))+
      geom_tile()+
      labs(x='x',y='y')+
      coord_fixed()+
      scale_y_reverse()
    
    #grid.arrange(plot_alpha, plot_alpha_v, plot_theta, plot_theta_v, plot_abund,ncol=2,top=paste0('t = ',t_i))
    grid.arrange(plot_alpha, plot_theta, plot_abund, plot_p, plot_p_v, nrow=2,top=paste0('t = ',t_i))
  }
  replot_data <- list(alpha_by_patch,theta_by_patch,pop_by_patch)
  return(replot_data)
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
    labs(x='x',y='y')+
    coord_fixed()+
    scale_y_reverse()
  
  ## map of theta values
  theta_by_patch <- group_by(sim_melt,patch,theta) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of alpha
    group_by(patch) %>%
    summarize(theta=sum(theta*popsize)/sum(popsize)) %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  plot_theta <- ggplot(theta_by_patch,aes(x=x-0.5,y=y-0.5,fill=theta))+
    geom_tile()+
    labs(x='x',y='y')+
    coord_fixed()+
    scale_y_reverse()
  
  ## map of population size
  pop_by_patch <- group_by(sim_melt,patch) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  plot_abund <- ggplot(pop_by_patch,aes(x=x-0.5,y=y-0.5,fill=popsize))+
    geom_tile()+
    labs(x='x',y='y')+
    coord_fixed()+
    scale_y_reverse()
  
  grid.arrange(plot_alpha, plot_theta, plot_abund,ncol=1,top=paste0('t = ',t))
}

############## quick diagnostic plotting function ##################

f_PlotOutput <- function(by_t,kern_timesteps,kern_xlim=25,patch_locations=NULL,nx=NULL,ny=NULL){
  p0 <- ggplot(by_t,aes(x=t))+
    geom_line(aes(y=alpha,color='alpha'))+
    geom_line(aes(y=theta,color='theta'))+
    labs(title='kernel parameters',y='value')+
    theme_minimal()+
    theme(legend.position = 'top')
  
  fill_colors <- c("First" = "white", "Last" = "red")
  outline_colors <- c("First" = "red", "Last" = "red")
  p1 <- ggplot(by_t,aes(x=alpha,y=theta))+
    geom_path(alpha=0.75,lwd=0.25)+
    theme_minimal()+
    geom_point(data=first(by_t),aes(x=alpha,y=theta,color="First",fill="First"),pch=21)+
    geom_point(data=last(by_t),aes(x=alpha,y=theta,color="Last",fill="Last"),pch=21)+
    scale_color_manual(name = "Point", values=outline_colors)+
    scale_fill_manual(name = "Point", values=fill_colors)+
    labs(title='kernel parameters')+
    theme(legend.position = 'top')
  
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
  
  p4 <- ggplot(by_t,aes(x=t,y=p))+
    geom_line()+
    theme_minimal()+
    labs(title='plasticity parameter')
  
  if(!is.null(patch_locations)){
    p5 <- f_Plot_Landscape(patch_locations,nx,ny,do_now=FALSE)
    grid.arrange(p0,p1,p2,p4,p3,p5, nrow=2)    
  }  
  else grid.arrange(p0,p1,p2,p4,p3,nrow=2)
}

f_plot_gamma <- function(alpha,theta,kern_xlim=10,...){
  g <- ggplot()+
    xlim(0,kern_xlim)+
    geom_function(fun=dgamma, args=list(shape=alpha,scale=theta))+
    labs(title=paste('alpha =',round(alpha,3),', theta =',round(theta,3)))
  
  print(g)
}

f_Plot_Landscape <- function(patch_locations,nx,ny,do_now=TRUE){
  g <- ggplot(patch_locations,aes(x=x,y=y))+
    geom_tile(aes(fill=K_i))+
    geom_vline(data=data.frame(x=seq(from=0.5,to=nx+0.5,by=1)),aes(xintercept=x),alpha=0.8,color='darkgray')+
    geom_hline(data=data.frame(y=seq(from=0.5,to=ny+0.5,by=1)),aes(yintercept=y),alpha=0.8,color='darkgray')+
    scale_y_reverse()+
    theme_minimal()+
    labs(title="Carrying capacity (K) by patch")
  if(do_now==TRUE) print(g)
  else return(g)
}
