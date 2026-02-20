source('0_Setup.R')

# test commit

#### Load a habitat (or skip this, and run the relevant sections of 2_Generate_Seascapes.R to make your own)
load(file="seascapes/hab_params_2.RData")


#### Set parameters (or adjust them as desired, if they were loaded above)
nsteps <- 100 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- seq(from=0.01,to=5,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.01,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 1 # index (in v_alphas) of shape parameter initial value
theta_start <- 1 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=0,to=0.9,length.out=5)# values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
disturb_prob=0 # probability of disturbance (per 10 timesteps)

params <- list(nsteps=nsteps,v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,disturb_prob=disturb_prob)

#### Run sim
sim_out <- f_RunSim(params,hab_params=make_hab_out,output_flag="all",show_plot = TRUE)
sim_melt <- f_ProcessPij(sim_out$Pij,patch_locations = make_hab_out$patch_locations,group_index = sim_out$group_index)$sim_melt
by_t <- f_ProcessPij(sim_out$Pij,patch_locations = make_hab_out$patch_locations,group_index = sim_out$group_index)$by_t
f_PlotAllHeatmapsDataTable(sim_melt,make_hab_out$patch_locations)

ggplot(by_t,aes(x=t_i,y=p))+
  geom_line()+
  labs(title="population mean plasticity parameter")
