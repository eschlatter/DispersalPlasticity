library(tidyverse)
library(scales)

source('f_GetConnectivityMatrix.R')

########## Parameters ##########
nx <- 10 # size of space in the x dimension
ny <- 10 # size of space in the y dimension

nsteps <- 50 # timesteps

v_alphas <- seq(from=0.1,to=5,length.out=20) # values the shape parameter can take
v_thetas <- seq(from=0.1,to=5,length.out=20) # values the scale parameter can take
alpha_start <- 20 # index (in v_alphas) of shape parameter initial value
theta_start <- 20 # index in (v_thetas) of scale parameter initial value

v_p <- seq(from=0,to=0.9,length.out=10)
p_start <- 1 # index (in v_p) of plasticity parameter initial value

mu <- 0.01 # mutation frequency
delta <- 0.001 # mutation magnitude

b <- 2 # number of offspring per parent

########## Data structures to describe space and dispersal ##########

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
patch_angles <- 2*asin(1/(2*patch_dists))/(2*pi)
patch_angles[is.nan(patch_angles)] <- 1

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

########## Simulation ##########

pop <- data.frame(t=1,
                  origin_site=NA, 
                  alpha=alpha_start, # note that these are the INDICES of the parameters in v_alphas and v_thetas, not the actual parameter values
                  theta=theta_start, 
                  dest_site=patch_locations$id)

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
  # alpha_adds <- sample(c(0,delta,-delta),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
  # theta_adds <- sample(c(0,delta,-delta),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
  # larvae$alpha <- larvae$alpha + alpha_adds
  # larvae$theta <- larvae$theta + theta_adds
  
  alpha_adds <- sample(c(0,1,-1),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
  theta_adds <- sample(c(0,1,-1),size=nrow(larvae),replace=TRUE,prob=c(1-mu/2,mu/4,mu/4))
  larvae$alpha <- oob_squish(larvae$alpha + alpha_adds, c(1,length(v_alphas)))
  larvae$theta <- oob_squish(larvae$theta + theta_adds, c(1,length(v_thetas)))
  
  # competition
  larvae <- larvae[sample(nrow(larvae),nrow(larvae),replace=FALSE),]
  larvae <- distinct(larvae,dest_site,.keep_all=TRUE)
  
  # cleanup for next timestep
  pop <- rbind(pop,larvae)
  
} # t

pop <- mutate(pop,alpha_value=v_alphas[alpha],theta_value=v_thetas[theta])

########## Make some plots ##########

by_t <- summarize(group_by(pop,t),alpha=mean(alpha_value),theta=mean(theta_value),popsize=n())

ggplot(by_t,aes(x=t))+
  geom_line(aes(y=alpha,color='alpha'))+
  geom_line(aes(y=theta,color='theta'))+
  labs(title='kernel parameters')

ggplot()+
  xlim(0,25)+
  geom_function(fun=dgamma, args=list(shape=first(by_t$alpha),scale=first(by_t$theta)),aes(color='first'))+
  geom_function(fun=dgamma, args=list(shape=last(by_t$alpha),scale=last(by_t$theta)),aes(color='last'))+
  labs(title='kernels')

ggplot(by_t,aes(x=t,y=popsize))+
  geom_line()+
  labs(title='population size')




