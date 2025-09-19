library(tidyverse)
library(scales)

source('sim_functions.R')

########## Parameters ##########
nx <- 10 # size of space in the x dimension
ny <- 10 # size of space in the y dimension

nsteps <- 1000 # timesteps

v_alphas <- seq(from=0.1,to=5,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.1,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 1 # index (in v_alphas) of shape parameter initial value
theta_start <- 1 # index in (v_thetas) of scale parameter initial value

v_p <- seq(from=0,to=0.9,length.out=10) # values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value

mu <- 0.01 # mutation frequency
delta <- 0.001 # mutation magnitude

b <- 2 # number of offspring per parent

########## Data structures to describe space and dispersal ##########
hab <- f_MakeHabitat(nx,ny,v_alphas,v_thetas)
patch_locations <- hab$patch_locations
patch_map <- hab$patch_map
patch_dists <- hab$patch_dists
patch_angles <- hab$patch_angles
conn_matrices <- hab$conn_matrices
npatch <- hab$npatch
rm(hab)

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

# process data for plotting (just summarize each timestep)
by_t <- summarize(group_by(pop,t),alpha=mean(alpha_value),theta=mean(theta_value),popsize=n())


# plots

ggplot(by_t,aes(x=t))+
  geom_line(aes(y=alpha,color='alpha'))+
  geom_line(aes(y=theta,color='theta'))+
  labs(title='kernel parameters')

ggplot()+
  xlim(0,25)+
  geom_function(fun=dgamma, args=list(shape=first(by_t$alpha),scale=first(by_t$theta)),aes(lty='first'))+
  geom_function(fun=dgamma, args=list(shape=last(by_t$alpha),scale=last(by_t$theta)),aes(lty='last'))+
  labs(title='kernels')

ggplot(by_t,aes(x=t,y=popsize))+
  geom_line()+
  labs(title='population size')

