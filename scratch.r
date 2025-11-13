library(tidyverse)
library(scales)
library(gridExtra)

source('functions/f_RunIBM.R')
source('functions/f_RunMatrixSim.R')
#source('functions/f_RunMatrixSim2.R')
source('functions/utility_functions.R')

# if(plot_kernel_dynamic==TRUE & t %% round(nsteps/100) == 0){
#   if(sum(Pij[,t],na.rm=TRUE)>0){
#     # histogram of kernel means (alpha*theta)
#     thisstep <- cbind(Pij[,t],matrix_index) %>%
#       rename(popsize=paste0('Pij[, t]')) %>%
#       filter(popsize>0) %>%
#       mutate(kern_mean=v_alphas[alpha]*v_thetas[theta],
#              kern_mode=ifelse(v_alphas[alpha]<1,0,(v_alphas[alpha]-1)*v_thetas[theta]))
#     break_vec <- seq(from=0,to=1,length.out=50)
#     #break_vec <- seq(from=min(thisstep$kern_mode),to=max(thisstep$kern_mode),length.out=20)
#     break_inds <- data.frame(kern_mode_bin=1:length(break_vec),l_end=break_vec,r_end=lead(break_vec))[1:19,]
#     thisstep <- mutate(thisstep, kern_mode_bin = cut(kern_mode,breaks=break_vec, include.lowest=TRUE,labels=FALSE)) %>%
#       group_by(kern_mode_bin) %>%
#       summarize(popsize=sum(popsize))%>%
#       right_join(break_inds,by='kern_mode_bin')
#     thisstep$popsize[is.na(thisstep$popsize)] <- 0
#     
#     mode_ticks <- expand.grid(alpha=v_alphas,theta=v_thetas) %>%
#       mutate(mode=ifelse(alpha<1,0,(alpha-1)*theta))
#     
#     g <- ggplot(thisstep,aes(x=l_end,y=popsize))+
#       geom_bar(stat='identity')+
#       labs(x='kernel mode',title=paste("t =", t))+
#       geom_point(data=filter(mode_ticks,mode<=max(break_vec)),aes(x=mode),y=0)+
#       xlim(-.04,max(break_vec))+
#       ylim(0,sum(K))
#     print(g)
#   } else print(paste("t=",t,': pop=0'))
# }


sim_array_t <- sim_array[,,,,2]
############## heatmap function ##################
# inputs:
#   sim_array_t: just the portion of sim_array from timestep t
#   patch_locations for mapmaking
f_PlotHeatmaps <- function(sim_array_t,patch_locations){
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
    labs(title='alpha',x='x',y='y')
  
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
    labs(title='theta',x='x',y='y')
  
  ## map of population size
  pop_by_patch <- group_by(sim_melt,patch) %>%
    summarize(popsize=sum(popsize),.groups='drop') %>%
    left_join(patch_locations,by=c("patch" = "id"))
  
  plot_abund <- ggplot(pop_by_patch,aes(x=x-0.5,y=y-0.5,fill=popsize))+
    geom_tile()+
    scale_x_continuous(breaks=0:10)+
    scale_y_continuous(breaks=0:10)+
    labs(title='abundance',x='x',y='y')
  
  grid.arrange(plot_alpha, plot_theta, plot_abund,ncol=1)
}


##############  maps of alpha, theta, popsize values (from matrix data) ###############
nx <- 10 # size of space in the x dimension
ny <- 10 # size of space in the y dimension

nsteps <- 1000 # timesteps

# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- seq(from=0.1,to=5,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.1,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 1 # index (in v_alphas) of shape parameter initial value
theta_start <- 1 # index in (v_thetas) of scale parameter initial value

# plasticity parameter currently doesn't affect anything or evolve
v_p <- seq(from=0,to=0.9,length.out=10) # values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value

mu <- 0.01 # mutation frequency
# mutation increment specified by v_alphas, v_thetas

matrix_out <- f_RunMatrixSim(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu=0.01,b,K=6)
sim_melt <- matrix_out[[1]]
by_t <- matrix_out[[2]]

## map of alpha values
alpha_by_patch <- group_by(sim_melt,patch,alpha,t) %>%
  summarize(popsize=sum(popsize),.groups='drop') %>%  # add up what's in the boxes with all values of theta
  group_by(patch,t) %>%
  summarize(alpha=sum(alpha*popsize)/sum(popsize)) %>%   # at each timestep and patch, find the mean value of alpha
  left_join(patch_locations,by=c("patch" = "id"))

ggplot(filter(alpha_by_patch,t==100),aes(x=x-0.5,y=y-0.5,fill=alpha))+
  geom_tile()+
  scale_x_continuous(breaks=0:10)+
  scale_y_continuous(breaks=0:10)

## map of theta values
theta_by_patch <- group_by(sim_melt,patch,theta,t) %>%
  summarize(popsize=sum(popsize),.groups='drop') %>% # add up what's in the boxes with all values of alpha
  group_by(patch,t) %>%
  summarize(theta=sum(theta*popsize)/sum(popsize)) %>%
  left_join(patch_locations,by=c("patch" = "id"))

ggplot(filter(theta_by_patch,t==100),aes(x=x-0.5,y=y-0.5,fill=theta))+
  geom_tile()+
  scale_x_continuous(breaks=0:10)+
  scale_y_continuous(breaks=0:10)

## map of population size
pop_by_patch <- group_by(sim_melt,patch,t) %>%
  summarize(popsize=sum(popsize),.groups='drop') %>%
  left_join(patch_locations,by=c("patch" = "id"))

ggplot(filter(pop_by_patch,t==100),aes(x=x-0.5,y=y-0.5,fill=popsize))+
  geom_tile()+
  scale_x_continuous(breaks=0:10)+
  scale_y_continuous(breaks=0:10)

############# plot a couple kernels from inner and outer patches, from matrix sim data ################

alpha_in <- as.numeric(filter(alpha_by_patch,t==100,patch==46)$alpha)
theta_in <- as.numeric(filter(theta_by_patch,t==100,patch==46)$theta)
alpha_edge <- as.numeric(filter(alpha_by_patch,t==100,patch==1)$alpha)
theta_edge <- as.numeric(filter(theta_by_patch,t==100,patch==1)$theta)

ggplot()+
  xlim(0,25)+
  geom_function(fun=dgamma, args=list(shape=alpha_in,scale=theta_in),aes(lty='inner'))+
  geom_function(fun=dgamma, args=list(shape=alpha_edge,scale=theta_edge),aes(lty='edge'))+
  theme_minimal()+
  theme(legend.position='top')


############### sample ###################
library(dplyr)
cellpop <- 60
dests <- 1:5
probs <- c(0.8,0.1,0.05,0.04,0.01)
sum(probs)

a <- data.frame(patch=sample(dests,size=cellpop,probs,replace=TRUE)) %>%
  group_by(patch) %>%
  summarize(immigrants=n())


# get x numbers that add up to 1 from a uniform distribution of probability 1/x
