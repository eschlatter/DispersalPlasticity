source('0_Setup.R')
load("experiments/Array1_SimIndex.RData") #load the dataframe of habitats to use

run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("Run",run_i))

# load the parameters
load('params/ParSetArray1.RData')
list2env(params,envir=environment())
p_start <- 1
alpha_start <- 3
theta_start <- 3
v_p <- c(0,0)

# load the habitat
#load(paste0("seascapes/kmaps/kmap_constant",run_i,".RData"))
load(paste0("seascapes/kmaps/kmap_",sim_index$Kmap[run_i],".RData"))
list2env(kmap_i,envir=environment())
rm(kmap_i)

# put everything in a list to pass to the simulation function
params <- list(nx=nx,ny=ny,nsteps=nsteps,
               v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,b=b,K=K,
               disturb_prob=disturb_prob,patch_locations=patch_locations)

# run the simulation
sim_out <- f_RunMatrixLoopLite(params,keep=list("p","abund","kern","sp_struct"))

# save the output
save(sim_out,file=paste0('output/array4/',run_i,'_array4.RData'))
