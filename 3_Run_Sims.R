source('0_Setup.R')

#### Load a habitat (or skip this, and run the relevant sections of 2_Generate_Seascapes.R to make your own)
basemap_file="experiments/Exp1_20260413"
hab_file="hab_test_full"
sim_file="test_full"

# q_autocorr=0.37
# rep_i=60
# hab_file <- paste0("set2/hab_",q_autocorr,"_sim",rep_i)
load(file=paste0(basemap_file,"/",hab_file,".RData"))

# hab_rast <- rast(paste0(basemap_file,"/",hab_file,".tif")) # load hab_rast
# hab_params <- c(hab_params,hab_rast=hab_rast)
# plot(hab_rast)

list2env(x=hab_params,envir=environment())

# circs=st_buffer(sfc_patches,dist=nav_rad)
# ggplot(reef_sf)+
#   geom_sf()+
#   geom_sf(data=sfc_patches)+
#   geom_sf(data=circs,alpha=0)+
#   annotation_scale()

#### Set (or load) simulation parameters
nsteps <- 10000 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- c(1) # values the shape parameter can take
v_thetas <- seq(from=0.25,to=4.75,by=0.5) # values the scale parameter can take
alpha_start <- 1 # index (in v_alphas) of shape parameter initial value
theta_start <- 0 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=-5,to=5,by=1)# values the plasticity parameter can take
p_start <- 0 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
disturb_prob=0 # probability of disturbance (per 10 timesteps)

params <- list(nsteps=nsteps,v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,disturb_prob=disturb_prob)
#save(params,file="seascapes/2026_02_26/params_1.RData")
#load("seascapes/2026_02_26/params_1.RData")

#### Run sim
f_RunSimNew(params,hab_params=hab_params,output_flag="all",show_plot = FALSE,output_thin=25,
            output_file=paste0(basemap_file,"/",sim_file))

#### Post-processing
popmat <- fread(paste0(basemap_file,"/",sim_file,"_raw.csv")) # population output file
load(paste0(basemap_file,"/",sim_file,"_metadata.RData")) # metadata file

sim_melt <- f_ProcessPopmat(popmat,metadata_list$patch_locations,metadata_list$group_index,
                                params=metadata_list$params)
save(sim_melt,file=paste0(basemap_file,"/",sim_file,"_processed.RData"))


# #### Make some quick check plots
# source('0_Setup.R')
# load(paste0(basemap_file,"/",sim_file,"_processed.RData")) # processed data
# load(paste0(basemap_file,"/",sim_file,"_metadata.RData")) # metadata file
# load(paste0(basemap_file,"/",metadata_list$hab_file,".RData")) # habitat file
# by_t <- sim_melt$by_t
# sim_melt <- sim_melt$sim_melt
# 
# a <- f_PlotAllHeatmaps(sim_melt = sim_melt,patch_locations = hab_params$patch_locations,reef_sf = hab_params$reef_sf,hab_type = hab_params$hab_type)
# ggplot(by_t,aes(x=t_i,y=p))+geom_line(aes(color="p"))+geom_line(aes(y=alpha*theta,color="kernmean"))
# ggplot(by_t,aes(x=t_i))+geom_line(aes(y=theta,color="theta"))+geom_line(aes(y=alpha,color="alpha"))+geom_line(aes(y=p,color="p"))
# 
# 
# 





#load("seascapes/2026_03_30/1x25km_res=10m/sim_out_hab_3_t1000.RData")
# load("seascapes/2026_03_04/sim_out_8.RData")
# load(file=paste0(hab_file,".RData"))
# 
# # plotting for output_flag="lite"
# ggplot(filter(sim_out$df_out,t_i %in% seq(from=1,to=max(sim_out$df_out$t_i),by=1)),aes(x=t_i,y=mean))+
#   geom_line()+
#   geom_ribbon(aes(ymin=mean-sqrt(var),ymax=mean+sqrt(var)),alpha=0.15)+
#   facet_wrap(vars(metric),scales="free")
