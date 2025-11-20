source('0_Setup.R')

# load parameters
load('params/ParSet1.RData')
list2env(pars_save,envir=environment())

# #simulate seascape (use real base map)
# load('seascapes/bze_map_sub.RData')
# bze_out <- f_GenerateMapWithK(base_map, K_range=c(20,35), h = 0.6, plot_flag=TRUE)
# list2env(bze_out,envir=environment())
# 
# # or simulate seascapes (use fractal base map)
# frac_out <- f_GenerateMapWithK(K_range=c(1,15), h=1.8, k=1, p=0.3, h_base=0.8, plot_flag=TRUE)
# list2env(frac_out,envir=environment())

## run the model
nsteps=300
patch_locations=NULL
sim_out <- f_RunMatrix(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad,b_neutral,b_good,K,patch_locations = patch_locations, disturb_prob=0.1)

## plot output
f_PlotOutput(sim_out$by_t, kern_timesteps=seq(from=0.75*nsteps,to=nsteps,length.out=25))
f_PlotAllHeatmaps(sim_out$sim_melt,sim_out$patch_locations)


a <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL)
