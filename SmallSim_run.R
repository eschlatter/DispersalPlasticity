source('0_Setup.R')
experiment_folder <- "experiments/SmallSim_Testing"
basemap_id=1
popmap_id=2
qmap_id=5
hab_id=5
temp_path <- paste0("/scratch.global/schla103/SmallSim_Testing/map",hab_id)
source("functions/sim_fns_parallel.R")
source("functions/f_RunSimComboStorageLite.R")

########## Load data ###########
# params
load(file=paste0(experiment_folder,"/params_1.RData")) 
list2env(x=params,envir=environment())
disturb_prob <- 0
params$disturb_prob <- 0

# hab_params
load(file=paste0(experiment_folder,"/habparams_",hab_id,".RData")) 

# patch_dists
patchdist_file <- paste0(experiment_folder,"/b",basemap_id,"/patchdists_b",basemap_id,
                         "_p",popmap_id,".RData")
nsteps=500
params$nsteps=500
disturb_prob <- 1
params$disturb_prob <- 1
output_flag="lite"
show_plot=FALSE
output_thin=1
output_file=paste0(experiment_folder,"/output/map",hab_id)
run_i=1
connmat_folder=temp_path
connmat_format="rds"
connmat_size_GB=0.33
jobmem_GB=10
start_location="everywhere"
normalize_connmats=TRUE
disturb_freq=nsteps
disturb_radius=0

# #### Run sim
# profvis({
  f_RunSimComboStorageLite(params,hab_params=hab_params,output_flag=output_flag,show_plot = show_plot,output_thin=output_thin,
                           output_file=output_file,run_i=run_i,
                           connmat_folder=connmat_folder,connmat_format=connmat_format,patchdist_file=patchdist_file,
                           connmat_size_GB = connmat_size_GB, jobmem_GB = jobmem_GB,normalize_connmats=normalize_connmats,
                           disturb_freq=disturb_freq,disturb_radius=disturb_radius)
# })

# # f_RunSimComboStorage(params,hab_params=hab_params,output_flag=output_flag,show_plot = show_plot,output_thin=output_thin,
# #                      output_file=output_file,run_i=run_i,
# #                      connmat_folder=connmat_folder,connmat_format=connmat_format,
# #                      patchdist_file=patchdist_file)
# # 
# prof_sim <- profvis({out_sim <- f_RunSimComboStorage_Parallel(params=params, hab_params=hab_params, keep=list("abund","p","kern","sp_struct"),
#                               output_flag=output_flag,show_plot=FALSE,output_thin=output_thin,output_file=output_file,run_i=run_i,
#                               connmat_folder=connmat_folder,connmat_format=connmat_format,patchdist_file=patchdist_file,
#                               connmat_size_GB=connmat_size_GB,jobmem_GB=jobmem_GB)
# })
# htmlwidgets::saveWidget(prof_sim,"profile.html")
# browseURL("profile.html")


#### Quick output
dat_out <- read.csv(paste0(experiment_folder,"/output/map",hab_id,"_raw.csv"))

abund_df <- filter(dat_out,metric=="larval_abund")
dat_out <- filter(dat_out,metric %in% c("p","theta","efftheta")) %>%
  mutate(metric=factor(metric,levels=c("p","theta","efftheta")))

g_mean <- ggplot(dat_out,aes(x=t_i))+
  geom_line(aes(y=median,color="median"))+
  geom_line(aes(y=mean,color="mean"))+
  geom_ribbon(aes(ymin=mean-sqrt(var),ymax=mean+sqrt(var)),alpha=0.2)+
  # geom_ribbon(aes(ymin=q05,ymax=q95),alpha=0.2)+
  # geom_ribbon(aes(ymin=q25,ymax=q75),alpha=0.2)+
  facet_wrap(vars(metric),scales="free",nrow=4)
g_moran <- ggplot(filter(dat_out,metric!="abund"),aes(x=t_i))+
  geom_line(aes(y=MoranI))+
  facet_wrap(vars(metric),nrow=4)
g_corr <- ggplot(filter(dat_out,metric!="abund"),aes(x=t_i))+
  geom_line(aes(y=corr_q))+
  geom_hline(yintercept = 0,lty="dashed")+
  facet_wrap(vars(metric),nrow=4)

grid.arrange(g_mean,g_corr,nrow=1)
ggplot(abund_df,aes(x=t_i,y=median))+geom_line()

dat_out |> filter(t_i>100) |> summarise(median=median(median), .by=metric)
