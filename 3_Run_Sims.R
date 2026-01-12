source('0_Setup.R')

#### Load a parameter set
load('params/ParSet4.RData')

#### Make adjustments to it as desired
# params$v_p <- -2:2 # possible plasticity values
# params$p_start <- 0 # starting plasticity index (0 = random distribution of starting values)
# params$alpha_start <- 0 # starting alpha index (0 = random distribution of starting values)
# params$theta_start <- 0 # starting theta index (0 = random distribution of starting values)
# params$mu=0.01 # mutation rate
# params$hab_type="grid" # grid-style versus point-style spatial sim
params$nsteps=10000 # number of simulation timesteps
params$nav_rad=1 # navigation radius (in km)
params$patch_locations$b_i <- sample(1:1000,size=nrow(params$patch_locations),replace = TRUE) # number of offspring per parent

#### Run sim
list2env(x=params,envir=environment())
# set up habitat data structures
makehab_output <- f_MakeHabitat(nx=nx,ny=ny,v_alphas=v_alphas,v_thetas=v_thetas,patch_locations=patch_locations,
                                hab_type=hab_type,nav_rad=nav_rad,numCores=parallelly::availableCores())
# run sim
# I've been working with f_RunMatrixLoopLite (in functions/f_RunMatrixLoop.R), which doesn't store all the population info at each timestep; it only stores summary stats.
# We'll probably want to go back to keeping everything at each timestep (that's f_RunMatrixLoop), at least temporarily. But I haven't updated that version of the function yet.
sim_loop1 <- f_RunMatrixLoopLite(params,show_plot = TRUE,makehab_output=makehab_output)
# quick diagnostic plots
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
