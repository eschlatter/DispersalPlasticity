library(tidyverse)

source('f_GetConnectivityMatrix.R')

########## Parameters ##########
nx <- 5 # size of space in the x dimension
ny <- 5 # size of space in the y dimension

nsteps <- 5000 # timesteps

alphas <- seq(from=0.1,to=5,length.out=10) # values the shape parameter can take
thetas <- seq(from=0.1,to=5,length.out=10) # values the scale parameter can take

alpha_start <- 1 # index (in alphas) of shape parameter initial value
theta_start <- 1 # index in (thetas) of scale parameter initial value

ps <- seq(from=0,to=0.9,length.out=10)
p_start <- 1 # index (in ps) of plasticity parameter initial value

mu <- 0.01 # mutation frequency
delta <- 0.001 # mutation magnitude

b <- 1.1 # number of offspring per parent

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

#check: proportion of individuals from each patch that land in a patch (same configuration as patch_map).
# Shouldn't be greater than 1 anywhere, and should be smaller where fewer patches are reachable.
cm <- f_GetConnectivityMatrix(1,2,patch_dists,patch_angles)
t(matrix(rowSums(cm),nrow=5)) 

# connectivity matrices
# (dimensions nalpha x ntheta x npatch x npatch)
# Should we calculate them for each possible kernel up front?
# There are probably some kernels that won't get used, so this might be a bit wasteful. But, for now, let's do it. We can be more efficient later.
# Should we allow kernels to evolve past the predefined ones? Maybe. (Probably?) Let's implement this later on.
conn_matrices <- array(NA,dim=c(length(alphas),length(thetas),npatch,npatch))
for(i_alpha in 1:length(alphas)){
  for(i_theta in 1:length(thetas)){
    conn_matrices[i_alpha,i_theta,,] <- f_GetConnectivityMatrix(alphas[i_alpha],thetas[i_theta],patch_dists,patch_angles)
  } # i_theta
} # i_alpha

########## Data structure to describe population ##########

# population tracking array
# dimensions:
# 1: location
# 2: alpha (kernel shape parameter)
# 3: theta (kernel scale parameter)
# 4: p (plasticity parameter)
# 5: timestep
sim_array <- array(0,dim=c(npatch,length(alphas),length(thetas),length(ps),nsteps))

########## Simulation ##########

# initialize starting population
# currently 1 individual (whatever that means) per patch
# everybody starts with the same kernel and p
sim_array[,alpha_start,theta_start,p_start,1] <- 1

for(t in 2:nsteps){
  # reproduction and dispersal and mutation (each patch contributes to other patches)
  for(i_alpha in 1:length(alphas)){
    for(i_theta in 1:length(thetas)){
      cm <- conn_matrices[i_alpha,i_theta,,]
      for(i_patch in 1:npatch){
        patch_pop <- sim_array[i_patch,i_alpha,i_theta,p_start,t-1]
        # no mutation
        sim_array[,i_alpha,i_theta,p_start,t] <- b*patch_pop*(1-mu)*cm[i_patch,] + sim_array[,i_alpha,i_theta,p_start,t]
        
        # mutation
        # "absorbing boundaries" at the edge of allowable kernel parameters
        # surely there's a more efficient way to do this, but we'll go with this for now
        if(i_alpha!=length(alphas)){
          sim_array[,i_alpha+1,i_theta,p_start,t] <- b*patch_pop*(mu/4)*cm[i_patch] + sim_array[,i_alpha+1,i_theta,p_start,t]
        }
        if(i_alpha!=1){
          sim_array[,i_alpha-1,i_theta,p_start,t] <- b*patch_pop*(mu/4)*cm[i_patch] + sim_array[,i_alpha-1,i_theta,p_start,t]
        }
        # mutation to theta
        if(i_theta!=length(thetas)){
          sim_array[,i_alpha,i_theta+1,p_start,t] <- b*patch_pop*(mu/4)*cm[i_patch] + sim_array[,i_alpha,i_theta+1,p_start,t]
        }
        if(i_theta!=1){
          sim_array[,i_alpha,i_theta-1,p_start,t] <- b*patch_pop*(mu/4)*cm[i_patch] + sim_array[,i_alpha,i_theta-1,p_start,t]
        }
      } # i_patch
    } # i_theta
  } # i_alpha
  
  # competition
  # want to cap the population of each patch, probably at 1 for now.
  # so, if the patch has population greater than 1, scale the value in each box by 1/(sum of all boxes for that patch)
  pop_by_patch <- apply(sim_array[,,,,t],1,sum)
  for(i in 1:length(pop_by_patch)){
    pop_by_patch[i] <- max(1,pop_by_patch[i])
  }
  sim_array[,,,,t] <- sim_array[,,,,t]/pop_by_patch
  
} # t

########## A little test output ##########

# plot mean alpha and theta over time
alphas_over_time <- apply(sim_array,c(2,5),sum)
thetas_over_time <- apply(sim_array,c(3,5),sum)

test_alpha <- as.vector(alphas %*% alphas_over_time)/colSums(alphas_over_time)
test_theta <- as.vector(thetas %*% thetas_over_time)/colSums(thetas_over_time)

plot(1:nsteps,test_alpha,type='l')
lines(1:nsteps,test_theta,type='l',col='blue')

# plot kernel at the start and end
plot(seq(from=0,to=10,by=0.01),dgamma(seq(from=0,to=10,by=0.01),shape=first(test_alpha),scale=first(test_theta)),type='l',lty='dashed',main='kernels',xlab='distance',ylab='probability')
lines(seq(from=0,to=10,by=0.01),dgamma(seq(from=0,to=10,by=0.01),shape=last(test_alpha),scale=last(test_theta)),type='l')

# plot population size over time
pop_over_time <- apply(sim_array,5,sum)
plot(1:nsteps,pop_over_time,type='l')
