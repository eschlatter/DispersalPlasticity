##### Square qmaps

source('0_Setup.R')
experiment_index <- read.csv("experiments/Exp3_20260528/experiment_index_disturbset.csv")

for(i in 1:4){
  qmap_id <- experiment_index$q_ID[i]
  autocorr_range <- round(experiment_index$range[i])
  qrast <- rast(paste0("experiments/Exp3_20260528/b1/set1/qmap_b1_q",qmap_id,".tif"))
  qrast_plot <- ggplot()+
    ggspatial::layer_spatial(qrast$q)+
    scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "#d2f2f7")+
    #geom_sf(data=hab_params$sfc_patches,size=0.05)+
    #annotation_scale()+
    labs(title=paste0(autocorr_range," m"))+
    theme(legend.position = "none",axis.text=element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust=0.5))
  print(qrast_plot)
}

##### Kimbe Bay qmap

source('0_Setup.R')
experiment_name="Exp4_20260611"
experiment_folder <- "experiments/Exp4_20260611"
basemap_id=3
popmap_id=3
qmap_id=3
hab_id=3
temp_path <- paste0("/scratch.global/schla103/",experiment_name,"/map_",hab_id)
source("functions/sim_fns_parallel.R")
source("functions/f_RunSimComboStorageLite.R")

# params
load(file=paste0(experiment_folder,"/params_2.RData")) 
list2env(x=params,envir=environment())
disturb_prob <- 1
params$disturb_prob <- 1

# hab_params
load(file=paste0(experiment_folder,"/habfiles/habparams_",hab_id,".RData")) 
hab_params$patch_locations$b <- hab_params$patch_locations$b/10

qrast <- rast(paste0("experiments/Exp4_20260611/b3/set1/qmap_b3_q",qmap_id,".tif"))
qrast_plot <- ggplot()+
  ggspatial::layer_spatial(qrast$q)+
  scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "#d2f2f7")+
  #geom_sf(data=hab_params$sfc_patches,size=0.05)+
  #annotation_scale()+
  labs(title=paste0(2893," m"))+
  theme(legend.position = "none",axis.text=element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust=0.5))
print(qrast_plot)


####### Kernels
source('0_Setup.R')
library(grid)
# output_folder <- "experiments/Exp3_20260528/output/smallerkern/rep2"
experiment_name <- "Exp3_20260528"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index_disturbset.csv"))
plotdynamics <- FALSE

load(file=paste0(experiment_folder,"/params_3.RData"))
v_thetas <- params$v_thetas
#v_thetas <- c(0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28) # for smallerkern
#v_thetas <- 10^(-4:1) # for smallkern
v_p <- -5:5
v_p_val <- (v_thetas[2]/v_thetas[1])^-v_p
df_p <- data.frame(p=v_p,pval=v_p_val)

df_plast <- data.frame(b=0:10)
df_plast$p1_2 <- (v_thetas[2]/v_thetas[1])^(f_plasticityb(b=0:10,p=1,alpha=1,theta=3,n_alpha=1,n_theta=8,bmin=1,bmax=9)$theta_plastic-3)
df_plast$p1_4 <- (v_thetas[2]/v_thetas[1])^(f_plasticityb(b=0:10,p=2,alpha=1,theta=3,n_alpha=1,n_theta=8,bmin=1,bmax=9)$theta_plastic-3)
df_plast$p2 <- (v_thetas[2]/v_thetas[1])^(f_plasticityb(b=0:10,p=-1,alpha=1,theta=3,n_alpha=1,n_theta=8,bmin=1,bmax=9)$theta_plastic-3)
df_plast$p4 <- (v_thetas[2]/v_thetas[1])^(f_plasticityb(b=0:10,p=-2,alpha=1,theta=3,n_alpha=1,n_theta=8,bmin=1,bmax=9)$theta_plastic-3)


ggplot(df_plast,aes(x=b))+
  geom_hline(aes(yintercept=1),linetype='dashed')+
  geom_line(aes(y=p1_4,color="1/4"))+
  geom_line(aes(y=p1_2,color="1/2"))+
  geom_line(aes(y=p2,color="2"))+
  geom_line(aes(y=p4,color="4"))+
  labs(x="habitat quality",y="adjustment to fundamental kernel mean\n(multiplication factor)",color="p")+
  scale_y_log10()


ggplot()+
  xlim(0,10)+
  geom_function(fun=dgamma, args=list(shape=1,scale=1),aes(color="1"))+
  geom_function(fun=dgamma, args=list(shape=1,scale=2),aes(color="2"))+
  labs(x="distance (km)",y="density",color="kernel mean (km)")+
  theme_minimal()

ggplot()+
  xlim(0,5)+
  lapply(1:length(v_thetas), 
         function(i){geom_function(fun=dgamma,
                                   args=list(shape=1,scale=v_thetas[i]),
                                   aes(color=factor(v_thetas[i])),
                                   n=5001)})+
  theme_minimal()+
  labs(x='distance (km)',y='density',color="theta",title="Dispersal Kernels")
