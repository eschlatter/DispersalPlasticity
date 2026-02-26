source('0_Setup.R')

#### Load a habitat (or skip this, and run the relevant sections of 2_Generate_Seascapes.R to make your own)
hab_file="seascapes/2026_02_26/hab_4_2_pt"
load(file=paste0(hab_file,".RData"))
hab_rast <- rast(paste0(hab_file,".tif")) # load hab_rast
hab_params <- c(hab_params,hab_rast=hab_rast)

#### Load simulation parameters (or set them below)

#### Set simulation parameters (or load them above)
nsteps <- 1000 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- c(1) # values the shape parameter can take
v_thetas <- seq(from=0.01,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 1 # index (in v_alphas) of shape parameter initial value
theta_start <- 3 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=-1,to=1,length.out=3)# values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
disturb_prob=0 # probability of disturbance (per 10 timesteps)

params <- list(nsteps=nsteps,v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,disturb_prob=disturb_prob)
#save(params,file="seascapes/2026_02_26/params_1.RData")

#### Run sim
sim_out <- f_RunSim(params,hab_params=hab_params,output_flag="all",show_plot = FALSE,output_thin=10)
# plotting for output_flag="all"
sim_processed <- f_ProcessPij(sim_out$Pij,patch_locations = hab_params$patch_locations,group_index = sim_out$group_index,output_thin=10)
sim_melt <- sim_processed$sim_melt
by_t <- sim_processed$by_t
save(sim_melt,by_t,hab_file,params,file="seascapes/2026_02_26/sim_out_4.RData")

load("seascapes/2026_02_26/sim_out_4.RData")

a <- f_PlotAllHeatmaps(sim_melt = sim_melt,patch_locations = hab_params$patch_locations,reef_sf = hab_params$reef_sf,hab_type = hab_params$hab_type)
ggplot(by_t,aes(x=t_i,y=p))+geom_line()

# plotting for output_flag="lite"
ggplot(filter(sim_out$df_out,t_i %in% seq(from=1,to=max(sim_out$df_out$t_i),by=1)),aes(x=t_i,y=mean))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sqrt(var),ymax=mean+sqrt(var)),alpha=0.15)+
  facet_wrap(vars(metric),scales="free")
