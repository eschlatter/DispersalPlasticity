source('0_Setup.R')

########## Param set 1: 
nx <- 20 # size of space in the x dimension
ny <- 20 # size of space in the y dimension
nsteps <- 1000 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- seq(from=0.01,to=5,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.01,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 3 # index (in v_alphas) of shape parameter initial value
theta_start <- 3 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=0,to=0.9,length.out=5)
#v_p <- seq(from=0,to=0.9,length.out=2) # values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
b <- 8 # reproductive output
K <- 10
#K <- round(runif(n=nx*ny,min=3,max=15)) # carrying capacity (also what plasticity responds to)

disturb_prob=0
patch_locations=NULL

params <- list(nx=nx,ny=ny,nsteps=nsteps,
               v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,b=b,K=K,
               disturb_prob=disturb_prob,patch_locations=NULL,seed=NULL)

#save(params,file='params/ParSet1.RData')

########## Param set 2:
source('0_Setup.R')
# simulate 33x33 seascape (use fractal base map)
nx <- 33 # size of space in the x dimension
ny <- 33 # size of space in the y dimension
nsteps <- 1000 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- seq(from=0.01,to=5,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.01,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 3 # index (in v_alphas) of shape parameter initial value
theta_start <- 3 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=0,to=0.9,length.out=5)
#v_p <- seq(from=0,to=0.9,length.out=2) # values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
b <- 8 # reproductive output
K <- 10
#K <- round(runif(n=nx*ny,min=3,max=15)) # carrying capacity (also what plasticity responds to)

disturb_prob=0

## new K overlay with existing base map
patch_map <- matrix(0,nrow=33,ncol=33)
for(i in 1:(33*33)){
  patch_map[params$patch_locations$y[i],params$patch_locations$x[i]] <- params$patch_locations$id[i]
}
patch_map[patch_map>0] <- 1
frac_out <- f_GenerateMapWithK(base_map=patch_map,K_range=c(15,15), h=1.8, k=5, p=0.5, h_base=0.8, plot_flag=TRUE)

frac_out <- f_GenerateMapWithK(K_range=c(15,15), h=1.8, k=5, p=0.5, h_base=0.8, plot_flag=TRUE)
list2env(frac_out,envir=environment())

params <- list(nx=nx,ny=ny,nsteps=nsteps,
               v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,b=b,K=K,
               disturb_prob=0,patch_locations=patch_locations,seed=NULL)

save(params,file='params/ParSet3.RData')

# -----------------------------------------------------
# # load parameters
# load('params/ParSet1.RData')
# list2env(pars_save,envir=environment())

# #simulate seascape (use real base map)
# load('seascapes/bze_map_sub.RData')
# bze_out <- f_GenerateMapWithK(base_map, K_range=c(20,35), h = 0.6, plot_flag=TRUE)
# list2env(bze_out,envir=environment())
# 