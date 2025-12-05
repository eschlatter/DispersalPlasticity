## run the model
source('0_Setup.R')
load('params/VariableKPars.RData')
params$v_p <- -3:3
params$p_start <- 1
params$nsteps=2000
params$mu=0.05
list2env(x=params,envir=environment())
sim_loop1 <- f_RunMatrixLoop(params)
# sim_loop2 <- f_RunMatrixLoop(params)
# all.equal(sim_loop1$Pij,sim_loop2$Pij)

# save(sim_loop1,file='UniformK.RData')
# load(file='UniformK.RData')
sim_loop_out1 <- f_ProcessLoopOutputDataTable(sim_loop1$params,sim_loop1$Pij,sim_loop1$group_index, sim_loop1$time_run)
# sim_loop_out2 <- f_ProcessLoopOutputDataTable(sim_loop2$params,sim_loop2$Pij,sim_loop2$group_index, sim_loop2$time_run)
# all.equal(sim_loop1$by_t[2:nsteps,],sim_loop_out2$by_t[2:nsteps,])

# ## other versions
# sim_loop_new <- f_RunMatrixLoopNew(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=40,file_out="model_runs/20251121_1.h5")
# sim_ibm <- f_RunIBM(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad,b_neutral,b_good,K,heatmap_plot_int=5000)
# sim_loop2 <- f_RunMatrixLoop2(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=40)
# sim_sparse <- f_RunMatrixSparse(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad,b_neutral,b_good,K,patch_locations = patch_locations, seed=NULL)
# # check that different model versions give the same output
# all.equal(sim_loop_out$by_t,sim_loop_out2$by_t)

## plot output
f_PlotOutput(sim_loop_out1$by_t, kern_timesteps=round(seq(from=0.75*nsteps,to=nsteps,length.out=25)),patch_locations=sim_loop_out1$patch_locations,nx=params$nx,ny=params$ny,sim_melt=sim_loop_out1$sim_melt)

# f_PlotOutput(sim_loop1$by_t, kern_timesteps=seq(from=0.75*nsteps,to=nsteps,length.out=25),patch_locations=patch_locations,nx=nx,ny=ny)
# f_PlotOutput(sim_loop_out2$by_t, kern_timesteps=seq(from=0.75*nsteps,to=nsteps,length.out=25),patch_locations=sim_loop_out2$patch_locations,nx=params$nx,ny=params$ny)
f_PlotAllHeatmaps(sim_loop_out1$sim_melt,sim_loop_out1$patch_locations,plot_int=nsteps)

# ## time profile
# Rprof(filename='profile_loop.out')
# sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
# Rprof(NULL)
# summaryRprof('profile_loop.out')
# 
# 
# ## memory profile
# Rprofmem(filename='profile_loop.out',threshold=10000)
# sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
# Rprofmem(NULL)
# noquote(readLines("profile_loop.out"))
