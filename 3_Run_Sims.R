source('0_Setup.R')

#### Load a parameter set (or skip this, and just set them directly below)
# load('params/ParSet4.RData')
# list2env(x=params,envir=environment())

#### Set parameters (or adjust them as desired, if they were loaded above)
nsteps <- 1000 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- seq(from=0.01,to=5,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.01,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 1 # index (in v_alphas) of shape parameter initial value
theta_start <- 1 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=0,to=0.9,length.out=5)# values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
disturb_prob=0 # probability of disturbance (per 10 timesteps)
nav_rad=0.5 # navigation radius (in km)
hab_type="points" # "points" or "grid"
keep=list("abund","p","kern") # types of summary stats to collect at each timestep

params <- list(nsteps=nsteps,v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,disturb_prob=disturb_prob,nav_rad=nav_rad,
               keep=keep,seed=NULL,hab_type=hab_type)


#### Run sim
list2env(x=params,envir=environment())

### set up anemone locations
# use a map of Kimbe Bay:
  pts_output <- f_SimPtsOnMap(map_extent="large",n_anems=100,show_map=TRUE)
  list2env(x=pts_output,envir=environment())
# use a grid:
  # patch_locations <- expand.grid(x=seq(from=1,by=0.1,length.out=10),y=seq(from=1,by=0.1,length.out=10))
  # patch_locations$id <- 1:nrow(patch_locations)
  # dists_mat=NULL

### set b for each anemone
patch_locations$b_i <- sample(1:1000,size=nrow(patch_locations),replace = TRUE) # number of offspring per parent

### set up habitat data structures
hab_output <- f_MakeHabitat(v_alphas=v_alphas,v_thetas=v_thetas, patch_locations=patch_locations,hab_type=hab_type,nav_rad=nav_rad,dists_mat=dists_mat,overlap_method=2)


### run sim
# Note: I've been working with f_RunMatrixLoopLite (in functions/f_RunMatrixLoop.R), which doesn't store all the population info at each timestep; it only stores summary stats.
# We may want to go back to keeping everything at each timestep (that's f_RunMatrixLoop), at least temporarily. But I haven't updated that version of the function recently.
sim_loop1 <- f_RunMatrixLoopLite(params,show_plot = FALSE,makehab_output=hab_output)
save(sim_loop1,file="output/20250120_sim4.RData")

### quick diagnostic plots
load(file="output/20250120_sim3.RData")
f_PlotOutput_Lite_Points(sim_loop1)






#### Code for checking speed and memory (haven't updated this recently)
# ## time profile
# Rprof(filename='profile_loop.out')
# sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
# Rprof(NULL)
# summaryRprof('profile_loop.out')
# 
# ## memory profile
# Rprofmem(filename='profile_loop.out',threshold=10000)
# sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
# Rprofmem(NULL)
# noquote(readLines("profile_loop.out"))
