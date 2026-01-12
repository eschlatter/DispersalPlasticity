# takes a dataframe sim_df (should be sim_melt, joined with patch_locations)
# plots the effective kernels (considering plasticity), with the thickness of the kernel line indicating abundance of that kernel in the dataset
f_PlotEffectiveKernels <- function(sim_df,v_alphas,v_thetas,patch_locations,plot_title=NULL){
  eff_pars <- f_plasticityK_new(K=sim_df$K,p=sim_df$p,alpha=sim_df$alpha,theta=sim_df$theta,n_alpha=length(v_alphas),n_theta=length(v_thetas),Kmin=min(patch_locations$K_i),Kmax=max(patch_locations$K_i))
  sim_df <- cbind(sim_df,eff_pars) 
  
  sim_df <- sim_df[,.(abund=sum(popsize)),by=c("alpha_plastic,theta_plastic")]
  sim_df <- sim_df[,abund_scale := (abund/max(abund))+.05]
  
  ggplot()+
    xlim(0,20)+
    lapply(1:nrow(sim_df), 
           function(i){geom_function(fun=dgamma,
                                     args=list(shape=sim_df$alpha_plastic[i],scale=sim_df$theta_plastic[i]),
                                     lwd=sim_df$abund_scale[i],
                                     color='blue')} )+
    theme_minimal()+
    labs(x='distance',y='density',title=plot_title)
}

# takes a dataframe sim_df (should be sim_melt, joined with patch_locations)
# plots the effective kernels (considering plasticity), with the thickness of the kernel line indicating abundance of that kernel in the dataset
f_PlotKernels <- function(sim_df,v_alphas,v_thetas,patch_locations,plot_title=NULL,effective=TRUE){
  p_input <- ifelse(effective==TRUE,sim_df$p,rep(0,times=nrow(sim_df))) # if effective==FALSE, then we set p=0, so the "effective" parameters are just alpha and theta
  eff_pars <- f_plasticityK_new(K=sim_df$K,p=p_input,alpha=sim_df$alpha,theta=sim_df$theta,n_alpha=length(v_alphas),n_theta=length(v_thetas),Kmin=min(patch_locations$K_i),Kmax=max(patch_locations$K_i))
  sim_df <- cbind(sim_df,eff_pars) 
  
  sim_df <- sim_df[,.(abund=sum(popsize)),by=c("alpha_plastic,theta_plastic")]
  sim_df <- sim_df[,abund_scale := (abund/max(abund))+.05]
  
  ggplot()+
    xlim(0,20)+
    ylim(0,1)+
    lapply(1:nrow(sim_df), 
           function(i){geom_function(fun=dgamma,
                                     args=list(shape=sim_df$alpha_plastic[i],scale=sim_df$theta_plastic[i]),
                                     lwd=sim_df$abund_scale[i],
                                     color=ifelse(effective==TRUE,'blue','black')
           )
           } )+
    theme_minimal()+
    labs(x='distance',y='density',title=paste0(ifelse(effective==TRUE,"Effective ",""),"Kernels",plot_title))
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

############## after-the-fact heatmap function for matrix model: data.table version ##################
# inputs:
#   replot_data: if it's been done before for this sim, can reuse the processed data from last time
#   sim_melt: dataframe with columns patch, timestep, alpha, theta, p, popsize
#   patch_locations for mapmaking
#   plot_int: specify timesteps to plot
f_PlotAllHeatmapsDataTable <- function(sim_melt,patch_locations,plot_int=NA,replot_data=NULL){
  if(prod(is.na(plot_int))==1){
    plot_int <- round(max(sim_melt$t)/10) # set the plotting interval, unless specified
    plot_ints <- seq(from=plot_int,to=max(sim_melt$t),by=plot_int)
  } 
  else plot_ints <- plot_int
  
  if(is.null(replot_data)){
    # process data
    # alpha
    alpha_by_patch <- sim_melt[,.(popsize=sum(popsize)),by=c("patch","alpha_value","t")]
    alpha_by_patch <- alpha_by_patch[,.(alpha_m=sum(alpha_value*popsize)/sum(popsize),
                                        alpha_v=sum(popsize*(alpha_value-(sum(alpha_value*popsize)/sum(popsize)))^2)/sum(popsize)),by=c("patch","t")]
    alpha_by_patch <- alpha_by_patch[patch_locations,on=c(patch="id")]
    
    # theta
    theta_by_patch <- sim_melt[,.(popsize=sum(popsize)),by=c("patch","theta_value","t")]
    theta_by_patch <- theta_by_patch[,.(theta_m=sum(theta_value*popsize)/sum(popsize),
                                        theta_v=sum(popsize*(theta_value-(sum(theta_value*popsize)/sum(popsize)))^2)/sum(popsize)),by=c("patch","t")]
    theta_by_patch <- theta_by_patch[patch_locations,on=c(patch="id")]
    
    # p
    p_by_patch <- sim_melt[,.(popsize=sum(popsize)),by=c("patch","p_value","t")]
    p_by_patch <- p_by_patch[,.(p_m=sum(p_value*popsize)/sum(popsize),
                                p_v=sum(popsize*(p_value-(sum(p_value*popsize)/sum(popsize)))^2)/sum(popsize)),by=c("patch","t")]
    p_by_patch <- p_by_patch[patch_locations,on=c(patch="id")]
    
    # abundance
    pop_by_patch <- sim_melt[,.(popsize=sum(popsize)),by=c("patch","t")]
    pop_by_patch <- pop_by_patch[patch_locations,on=c(patch="id")]
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
  replot_data <- list(alpha_by_patch=alpha_by_patch,theta_by_patch=theta_by_patch,p_by_patch=p_by_patch,pop_by_patch=pop_by_patch)
  return(replot_data)
}

############## after-the-fact heatmap function for matrix model ##################
# inputs:
#   replot_data: if it's been done before for this sim, can reuse the processed data from last time
#   sim_melt: dataframe with columns patch, timestep, alpha, theta, p, popsize
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
  replot_data <- list(alpha_by_patch=alpha_by_patch,theta_by_patch=theta_by_patch,p_by_patch=p_by_patch,pop_by_patch=pop_by_patch)
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

f_PlotOutput <- function(by_t,kern_timesteps,kern_xlim=20,patch_locations=NULL,nx=NULL,ny=NULL,sim_melt=NULL){
  p0 <- ggplot(by_t,aes(x=t))+
    geom_line(aes(y=alpha,color='alpha'))+
    geom_line(aes(y=theta,color='theta'))+
    labs(title='kernel parameters',y='value')+
    theme_minimal()+
    theme(legend.position = 'inside')
  
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
    theme(legend.position = 'inside')
  
  p2 <- ggplot()+
    xlim(0,kern_xlim)+
    geom_function(fun=dgamma, args=list(shape=first(by_t$alpha),scale=first(by_t$theta)),aes(lty='first'))+
    lapply(kern_timesteps, function(i){geom_function(fun=dgamma,args=list(shape=by_t$alpha[i],scale=by_t$theta[i]),alpha=0.5,color='darkgray')})+
    geom_function(fun=dgamma,args=list(shape=median(by_t$alpha[kern_timesteps]),scale=median(by_t$theta[kern_timesteps])),color='black',lwd=0.75,aes(lty='last'))+
    scale_linetype_manual(values=c('first' = 'dashed',
                                   'last' = 'solid'))+
    theme_minimal()+
    theme(legend.position='inside')+
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
  }
  else p5 <- NULL
  
  if(!is.null(sim_melt)){
    # make a dataframe to give to f_PlotEffectiveKernels
    sim_df <- sim_melt[patch_locations,on=c(patch="id")] # join sim_melt and patch_locations
    sim_df <- sim_df[sim_df$popsize!=0,] # remove zero rows
    sim_df <- sim_df[(sim_df$t %in% kern_timesteps),] # choose only needed timesteps
    
    p6 <- f_PlotEffectiveKernels(sim_df,v_alphas,v_thetas,patch_locations,plot_title="Effective kernels, all patches")
    
    sim_df_island <- sim_df[(sim_df$x<15),]
    p7 <- f_PlotEffectiveKernels(sim_df_island,v_alphas,v_thetas,patch_locations,plot_title="Effective kernels, island")
    
    sim_df_mainland <- sim_df[(sim_df$x>15),]
    p8 <- f_PlotEffectiveKernels(sim_df_mainland,v_alphas,v_thetas,patch_locations,plot_title="Effective kernels, mainland")
    
    sim_df_last <- sim_df[(sim_df$t==nsteps),]
    sim_df_last <- sim_df_last[,.(abund=sum(popsize),x=first(x),y=first(y)),by=patch]
    p9 <- ggplot(sim_df_last)+
      geom_tile(aes(x=x,y=y,fill=abund))+
      scale_y_reverse()+
      labs(title=paste0("End population size, t=",nsteps))
  }
  else{
    p6 <- NULL
    p7 <- NULL
    p8 <- NULL
    p9 <- NULL
  } 
  
  grid.arrange(p0,p1,p2,p4,p3,p5,p6,p7,p8,p9,ncol=3)
}

############## diagnostic plotting function, reworked ##################

f_PlotOutputNew <- function(by_t,kern_timesteps=NULL,kern_xlim=20,patch_locations=NULL,nx=NULL,ny=NULL,sim_melt=NULL,
                            island_lims,mainland_lims){
  # if not specified, pick a reasonable sample of the ending timesteps to calculate patchwise kernels from
  if(is.null(kern_timesteps)) kern_timesteps=seq(from=round(0.75*nrow(by_t)), to=nrow(by_t),length.out=25)
  
  # make a dataframe to give to f_PlotEffectiveKernels
  sim_df <- sim_melt[patch_locations,on=c(patch="id")] # join sim_melt and patch_locations
  sim_df <- sim_df[sim_df$popsize!=0,] # remove zero rows
  
  # Page 1
  #-----------------------------
  # K across the seascape
  p1 <- f_Plot_Landscape(patch_locations,nx,ny,do_now=FALSE)
  dev.new(width=2,height=2,units="in")
  print(p1)
  
  # Pop size over time
  p2 <- ggplot(by_t,aes(x=t,y=popsize))+
    geom_line()+
    theme_minimal()+
    labs(title='population size')
  
  # By-patch pop size at the end of the sim
  sim_df_last <- sim_df[(sim_df$t==nsteps),]
  sim_df_last <- sim_df_last[,.(abund=sum(popsize),x=first(x),y=first(y)),by=patch]
  p3 <- ggplot(sim_df_last)+
    geom_tile(aes(x=x,y=y,fill=abund))+
    scale_y_reverse()+
    labs(title=paste0("End population size, t=",nsteps))
  
  # traces of kernel parameters
  p4 <- ggplot(by_t,aes(x=t))+
    geom_line(aes(y=alpha,color='alpha'))+
    geom_line(aes(y=theta,color='theta'))+
    labs(title='kernel parameters',y='value')+
    theme_minimal()+
    theme(legend.position = 'inside')
  
  # kernel parameters in alpha-theta space
  fill_colors <- c("First" = "white", "Last" = "red")
  outline_colors <- c("First" = "red", "Last" = "red")
  p5 <- ggplot(by_t,aes(x=alpha,y=theta))+
    geom_path(alpha=0.75,lwd=0.25)+
    theme_minimal()+
    geom_point(data=first(by_t),aes(x=alpha,y=theta,color="First",fill="First"),pch=21)+
    geom_point(data=last(by_t),aes(x=alpha,y=theta,color="Last",fill="Last"),pch=21)+
    scale_color_manual(name = "Point", values=outline_colors)+
    scale_fill_manual(name = "Point", values=fill_colors)+
    labs(title='kernel parameters')+
    theme(legend.position = 'inside')
  
  # trace of plasticity parameter
  p6 <- ggplot(by_t,aes(x=t,y=p))+
    geom_line()+
    theme_minimal()+
    labs(title='plasticity parameter')
  
  a <- dev.list()
  dev.set(which=as.numeric(a['RStudioGD']))
  grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
  
  # Page 2
  #-----------------------------
  # prep sim_df for by-patch stuff
  sim_df_start <- sim_df[(sim_df$t==1),]
  sim_df <- sim_df[(sim_df$t %in% kern_timesteps),] # choose only needed timesteps
  sim_df_island <- sim_df[(sim_df$x %in% island_lims$xmin:island_lims$xmax) & (sim_df$y %in% island_lims$ymin:island_lims$ymax),]
  sim_df_mainland <- sim_df[(sim_df$x %in% mainland_lims$xmin:mainland_lims$xmax) & (sim_df$y %in% mainland_lims$ymin:mainland_lims$ymax),]
  
  
  ## Non-effective kernels (actual parameter values)
  # starting
  p7 <- f_PlotKernels(sim_df=sim_df_start,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                      plot_title="\nAll patches \nt=1",effective=FALSE)
  
  # all
  p8 <- f_PlotKernels(sim_df=sim_df,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                      plot_title=paste0("\nAll patches \nt=",min(kern_timesteps),"-",max(kern_timesteps)),effective=FALSE)
  
  # island
  p9 <- f_PlotKernels(sim_df=sim_df_island,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                      plot_title=paste0("\nIsland \nt=",min(kern_timesteps),"-",max(kern_timesteps)),effective=FALSE)
  # mainland
  p10 <- f_PlotKernels(sim_df=sim_df_mainland,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                       plot_title=paste0("\nMainland \nt=",min(kern_timesteps),"-",max(kern_timesteps)),effective=FALSE)
  
  ## Effective kernels
  # starting
  p11 <- f_PlotKernels(sim_df=sim_df_start,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                       plot_title="\nAll patches \nt=1",effective=TRUE)
  # all
  p12 <- f_PlotKernels(sim_df=sim_df,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                       plot_title=paste0("\nAll patches \nt=",min(kern_timesteps),"-",max(kern_timesteps)),effective=TRUE)
  
  # island
  p13 <- f_PlotKernels(sim_df=sim_df_island,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                       plot_title=paste0("\nIsland \nt=",min(kern_timesteps),"-",max(kern_timesteps)),effective=TRUE)
  # mainland
  p14 <- f_PlotKernels(sim_df=sim_df_mainland,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                       plot_title=paste0("\nMainland \nt=",min(kern_timesteps),"-",max(kern_timesteps)),effective=TRUE)
  
  patch_locations$whichland <- ifelse((patch_locations$x %in% island_lims$xmin:island_lims$xmax) & (patch_locations$y %in% island_lims$ymin:island_lims$ymax),
                                      "island",ifelse((patch_locations$x %in% mainland_lims$xmin:mainland_lims$xmax) & (patch_locations$y %in% mainland_lims$ymin:mainland_lims$ymax),
                                                      "mainland","neither"))
  p15 <- ggplot(patch_locations)+
    geom_tile(aes(x=x,y=y,fill=whichland))+
    scale_y_reverse()+
    theme_minimal()
  
  grid.arrange(p7,p8,p9,p10,p11,p12,p13,p14,ncol=4)
  
  print(p15)
  
}

## make some dynamics plots based on output of the lite version of the simulation function (which just tracks summary statistics, not the whole population)
## input: output_lite, the list returned by f_RunMatrixLoopLite
## returns: none, but makes some graphs
f_PlotOutput_Lite <- function(output_lite){
  list2env(output_lite$params,environment())
  
  dat <- output_lite$output_df[2:nsteps,] # remove the first row, because we didn't calculate stats for it
  
  ## plot map with K
  g_map <- f_Plot_Landscape(patch_locations,nx,ny,do_now=FALSE)
  
  ## plot abundance
  g_abund <- ggplot(dat,aes(x=t,y=v_abund))+
    geom_line()+
    theme_minimal()+
    labs(title="Population size", x="time",y="population size")
  #  print(g_abund)
  
  ## plot p (mean +/- sd) at each timestep
  dat$p_lower <- dat$v_pmeans-sqrt(dat$v_pvars)
  dat$p_upper <- dat$v_pmeans+sqrt(dat$v_pvars)
  g_p <- ggplot(dat,aes(x=t,y=v_pmeans))+
    geom_line()+
    geom_ribbon(aes(ymin=p_lower,ymax=p_upper),alpha=0.15)+
    theme_minimal()+
    labs(title="Plasticity parameter",x='time',y='p (mean +/- sd)')
  #  print(g_p)
  
  ### Plot fundamental and effective kernel means (mean +/- sd) at each timestep
  dat$fund_mean_lower <- dat$fund_mean_mean-sqrt(dat$fund_mean_var)
  dat$fund_mean_upper <- dat$fund_mean_mean+sqrt(dat$fund_mean_var)
  dat$eff_mean_lower <- dat$eff_mean_mean-sqrt(dat$eff_mean_var)
  dat$eff_mean_upper <- dat$eff_mean_mean+sqrt(dat$eff_mean_var)
  g_kernmeans <- ggplot(dat,aes(x=t))+
    geom_line(aes(y=fund_mean_mean,color="Fundamental"))+
    geom_ribbon(aes(ymin=fund_mean_lower,ymax=fund_mean_upper,fill="Fundamental"),alpha=0.15)+
    geom_line(aes(y=eff_mean_mean,color="Effective"))+
    geom_ribbon(aes(ymin=eff_mean_lower,ymax=eff_mean_upper,fill="Effective"),alpha=0.15)+
    theme_minimal()+
    labs(title="Fundamental and effective kernel means",x='time',y='kernel mean (mean +/- sd at each timestep)')
  #  print(g_kernmeans)
  
  # dat$fund_mode_lower <- dat$fund_mode_mean-sqrt(dat$fund_mode_var)
  # dat$fund_mode_upper <- dat$fund_mode_mean+sqrt(dat$fund_mode_var)
  # dat$eff_mode_lower <- dat$eff_mode_mean-sqrt(dat$eff_mode_var)
  # dat$eff_mode_upper <- dat$eff_mode_mean+sqrt(dat$eff_mode_var)
  # ## plot properties of the fundamental kernel at each timestep
  # p_fundkern <- ggplot(dat,aes(x=t))+
  #   geom_line(aes(y=fund_mode_mean,color="mode"))+
  #   geom_ribbon(aes(ymin=fund_mode_lower,ymax=fund_mode_upper,fill="mode"),alpha=0.15)+
  #   geom_line(aes(y=fund_mean_mean,color="mean"))+
  #   geom_ribbon(aes(ymin=fund_mean_lower,ymax=fund_mean_upper,fill="mean"),alpha=0.15)+
  #   theme_minimal()+
  #   labs(title="Fundamental kernel",x='time',y='kernel property (mean +/- sd)')
  # ## plot properties of the effective kernel at each timestep
  # p_effkern <- ggplot(dat,aes(x=t))+
  #   geom_line(aes(y=eff_mode_mean,color="mode"))+
  #   geom_ribbon(aes(ymin=eff_mode_lower,ymax=eff_mode_upper,fill="mode"),alpha=0.15)+
  #   geom_line(aes(y=eff_mean_mean,color="mean"))+
  #   geom_ribbon(aes(ymin=eff_mean_lower,ymax=eff_mean_upper,fill="mean"),alpha=0.15)+
  #   theme_minimal()+
  #   labs(title="Effective kernel",x='time',y='kernel property (mean +/- sd)')
  # grid.arrange(p_fundkern,p_effkern,ncol=1)
  
  ### Plot spatial autocorrelation (global Moran's I) of both fundamental and effective kernel means
  g_autocorr <- ggplot(dat,aes(x=t))+
    geom_line(aes(y=fund_mean_moran,color='Fundamental'))+
    geom_line(aes(y=eff_mean_moran,color="Effective"))+
    theme_minimal()+
    labs(title="Spatial autocorrelation of kernel means",x = "time", y="Global Moran's I")
  
  grid.arrange(g_map,g_abund,g_p,g_kernmeans,g_autocorr,nrow=2)
}


## make some dynamics plots based on output of the lite version of the simulation function (which just tracks summary statistics, not the whole population)
## version for when the simulation tracks individual anemones, rather than habitat patches
## input: output_lite, the list returned by f_RunMatrixLoopLite
## returns: none, but makes some graphs
f_PlotOutput_Lite_Points <- function(output_lite){
  list2env(output_lite$params,environment())
  
  ## plot map with navigation radius circles
  circs=st_buffer(patch_sf,dist=nav_rad*1000) # dist is the radius of the circle, in meters
  if(!exists("map_sf")) map_sf=NULL
  g_map <- ggplot()+
    geom_sf(data=map_sf)+
    geom_sf(data=circs,fill=NA,color=alpha("gray",0.2))+
    geom_sf(data=patch_sf,size=0.2)+
    annotation_scale()+
    theme_minimal()+
    labs(title="Sites, with navigation radius")
  
  ## plot abundance
  abund_dat <- output_lite$output_list$df_abund[2:nsteps,] # remove the first row, because we didn't calculate stats for it
  g_abund <- ggplot(abund_dat,aes(x=t,y=abund))+
    geom_line()+
    theme_minimal()+
    labs(title="Population size", x="time",y="population size")
  
  ## plot p (mean +/- sd) at each timestep
  p_dat <- output_lite$output_list$df_p[2:nsteps,]
  p_dat$p_lower <- p_dat$mean-sqrt(p_dat$var)
  p_dat$p_upper <- p_dat$mean+sqrt(p_dat$var)
  g_p <- ggplot(p_dat,aes(x=t,y=mean))+
    geom_line()+
    geom_ribbon(aes(ymin=p_lower,ymax=p_upper),alpha=0.15)+
    theme_minimal()+
    labs(title="Plasticity parameter",x='time',y='p (mean +/- sd)')

  ### Plot fundamental and effective kernel means (mean +/- sd) at each timestep
  fund_dat <- output_lite$output_list$df_fund[2:nsteps,]
  fund_dat$fund_mean_lower <- fund_dat$kernmean_mean-sqrt(fund_dat$kernmean_var)
  fund_dat$fund_mean_upper <- fund_dat$kernmean_mean+sqrt(fund_dat$kernmean_var)
  
  eff_dat <- output_lite$output_list$df_eff[2:nsteps,]
  eff_dat$eff_mean_lower <- eff_dat$kernmean_mean-sqrt(eff_dat$kernmean_var)
  eff_dat$eff_mean_upper <- eff_dat$kernmean_mean+sqrt(eff_dat$kernmean_var)
  
  g_kernmeans <- ggplot(fund_dat,aes(x=t))+
    geom_line(aes(y=kernmean_mean,color="Fundamental"))+
    geom_ribbon(aes(ymin=fund_mean_lower,ymax=fund_mean_upper,fill="Fundamental"),alpha=0.15)+
    geom_line(data=eff_dat,aes(y=kernmean_mean,color="Effective"))+
    geom_ribbon(data=eff_dat,aes(ymin=eff_mean_lower,ymax=eff_mean_upper,fill="Effective"),alpha=0.15)+
    theme_minimal()+
    labs(title="Fundamental and effective kernel means",x='time',y='kernel mean (mean +/- sd at each timestep)')
  
  ### Arrange all plots
  grid.arrange(g_map,g_abund,g_p,g_kernmeans,nrow=2)
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
#    geom_vline(data=data.frame(x=seq(from=0.5,to=nx+0.5,by=1)),aes(xintercept=x),alpha=0.8,color='darkgray')+
#    geom_hline(data=data.frame(y=seq(from=0.5,to=ny+0.5,by=1)),aes(yintercept=y),alpha=0.8,color='darkgray')+
    scale_y_reverse()+
    theme_minimal()+
    labs(title="Carrying capacity (K) by patch")
  if(do_now==TRUE) print(g)
  else return(g)
}
