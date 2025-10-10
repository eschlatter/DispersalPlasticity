library(gridExtra)

############## heatmap function for IBM ###############
# inputs:
#   pop_t: population dataframe for the current timestep
#   patch_locations for mapmaking
#   t current timestep
f_PlotHeatmapsIBM <- function(pop_t,patch_locations,t){
  by_patch <- pop_t %>%
    group_by(dest_site) %>%
    summarize(alpha=mean(alpha),theta=mean(theta),popsize=n()) %>%
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

############## heatmap function ##################
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

f_PlotOutput <- function(by_t,kern_timesteps,kern_xlim=25){
  p0 <- ggplot(by_t,aes(x=t))+
    geom_line(aes(y=alpha,color='alpha'))+
    geom_line(aes(y=theta,color='theta'))+
    labs(title='kernel parameters',y='value')+
    theme_minimal()+
    theme(legend.position = 'top')
  
  p1 <- ggplot(by_t,aes(x=alpha,y=theta))+
    geom_line(alpha=0.75,lwd=0.25)+
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

f_MakeHabitat <- function(nx,ny,v_alphas,v_thetas){
  # list of patch locations and IDs
  # (dimensions npatch x 3)
  # "location" is the center of the patch
  patch_locations <- expand.grid(x=1:nx,y=1:ny) %>%
    rowid_to_column(var='id')
  # this is a simple version of patch_locations (all the squares of a grid); eventually we'll want to import a map.
  # I think we'll be able to just keep track of reef patches, not open ocean
  npatch <- nrow(patch_locations)
  
  # a "map" of the patch numbers, spatially arranged
  # (dimensions nx x ny)
  patch_map <- matrix(nrow=ny,ncol=nx)
  for(i in 1:npatch){
    patch_map[patch_locations$y[i],patch_locations$x[i]] <- patch_locations$id[i]
  }
  
  # a matrix of distances between the centers of each patch (i.e., r in polar coords)
  # (dimensions npatch x npatch)
  patch_dists <- matrix(nrow=npatch,ncol=npatch)
  colnames(patch_dists) <- patch_locations$id
  for(i in 1:npatch){
    patch_dists[i,] = sqrt((patch_locations$x[i]-patch_locations$x)^2+(patch_locations$y[i]-patch_locations$y)^2)
  }
  
  # a matrix of the size of the pie wedge between each patch (i.e., theta in polar coords -- not theta of the dispersal kernel)
  # assuming the width of the cell at the given distance is 1: not quite correct most of the time, but probably close enough
  # (dimensions npatch x npatch)
  patch_angles <- suppressWarnings(2*asin(1/(2*patch_dists))/(2*pi))
  patch_angles[is.nan(patch_angles)] <- 1
  
  #check: proportion of individuals from each patch that land in a patch (same configuration as patch_map).
  # Shouldn't be greater than 1 anywhere, and should be smaller where fewer patches are reachable.
  cm <- f_GetConnectivityMatrix(1,2,patch_dists,patch_angles)
  t(matrix(rowSums(cm),nrow=nrow(patch_map))) 
  
  # connectivity matrices
  # (dimensions nalpha x ntheta x npatch x npatch)
  # Should we calculate them for each possible kernel up front?
  # There are probably some kernels that won't get used, so this might be a bit wasteful. But, for now, let's do it. We can be more efficient later.
  # Should we allow kernels to evolve past the predefined ones? Maybe. (Probably?) Let's implement this later on.
  conn_matrices <- array(NA,dim=c(length(v_alphas),length(v_thetas),npatch,npatch))
  for(i_alpha in 1:length(v_alphas)){
    for(i_theta in 1:length(v_thetas)){
      conn_matrices[i_alpha,i_theta,,] <- f_GetConnectivityMatrix(v_alphas[i_alpha],v_thetas[i_theta],patch_dists,patch_angles)
    } # i_theta
  } # i_alpha
  
  return(list(patch_locations=patch_locations,
              patch_map=patch_map,
              patch_dists=patch_dists,
              patch_angles=patch_angles,
              conn_matrices=conn_matrices,
              npatch=npatch))
}


# inputs: kernel, seascape
# output: rates of dispersal from each patch to each other patch
  # matrix Connectivity: dimensions npatch x npatch
  # Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i
f_GetConnectivityMatrix <- function(alpha, theta, patch_dists, patch_angles){
  connectivity_matrix <- (pgamma(patch_dists+0.5,shape=alpha,scale=theta)-pgamma(patch_dists-0.5,shape=alpha,scale=theta))*patch_angles
}

# sample simple kernel function (though I think we'll just use the built-in gamma function eventually)
f_kernel <- function(distance){
  if(distance<=3) probability = 1/3
  else probability = 0
  probability
}