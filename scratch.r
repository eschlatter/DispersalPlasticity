library(tidyverse)
library(scales)
library(gridExtra)

source('functions/f_RunIBM.R')
source('functions/f_RunMatrixSim.R')
#source('functions/f_RunMatrixSim2.R')
source('functions/utility_functions.R')

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
