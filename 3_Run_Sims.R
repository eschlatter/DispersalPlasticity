source('0_Setup.R')

# # load parameters
# load('params/ParSet1.RData')
# list2env(pars_save,envir=environment())

# #simulate seascape (use real base map)
# load('seascapes/bze_map_sub.RData')
# bze_out <- f_GenerateMapWithK(base_map, K_range=c(20,35), h = 0.6, plot_flag=TRUE)
# list2env(bze_out,envir=environment())
# 
# or simulate seascapes (use fractal base map)
frac_out <- f_GenerateMapWithK(K_range=c(1,15), h=1.8, k=5, p=0.5, h_base=0.8, plot_flag=TRUE)
list2env(frac_out,envir=environment())

## run the model
source('1_Generate_Param_Sets.R')
patch_locations=NULL
nsteps=1000
sim_sparse <- f_RunMatrixSparse(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad,b_neutral,b_good,K,patch_locations = patch_locations, seed=NULL)
sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0.1,patch_locations=NULL, seed=NULL)
sim_loop_new <- f_RunMatrixLoopNew(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=40,file_out="model_runs/20251121_1.h5")
#sim_ibm <- f_RunIBM(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad,b_neutral,b_good,K,heatmap_plot_int=5000)

sim_loop1 <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=40)
sim_loop2 <- f_RunMatrixLoop2(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=40)

# check any changes didn't affect output
all.equal(sim_loop1$by_t,sim_loop2$by_t)

# timing
sim_sparse$time_run
sim_loop$time_run
sim_loop2$time_run


## plot output
f_PlotOutput(sim_sparse$by_t, kern_timesteps=seq(from=0.75*nsteps,to=nsteps,length.out=25))
f_PlotAllHeatmaps(sim_sparse$sim_melt,sim_sparse$patch_locations)

f_PlotOutput(sim_loop$by_t, kern_timesteps=seq(from=0.75*nsteps,to=nsteps,length.out=25))
f_PlotAllHeatmaps(sim_loop$sim_melt,sim_loop$patch_locations,plot_int=seq(from=489,to=500,by=1))

f_PlotOutput(sim_loop_new$by_t, kern_timesteps=seq(from=0.75*nsteps,to=nsteps,length.out=25))
f_PlotAllHeatmaps(sim_loop_new$sim_melt,sim_loop_new$patch_locations)

#f_PlotOutput(ibm_sim[[2]], kern_timesteps=seq(from=0.75*nsteps,to=nsteps,length.out=25))

## time profile
Rprof(filename='profile_loop.out')
sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
Rprof(NULL)
summaryRprof('profile_loop.out')


## memory profile
Rprofmem(filename='profile_loop.out',threshold=10000)
sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
Rprofmem(NULL)
noquote(readLines("profile_loop.out"))
