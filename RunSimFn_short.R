# library and working directory paths
if(version$major=="4" & version$minor=="4.0") .libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old") else .libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")

# load packages and functions
library(data.table)
library(parallel)
library(terra)
library(sf)
library(calculus)
library(dplyr)
source("functions/Function_simulate.R")
source("functions/Functions_auxiliary.R")

repID <- ifelse(is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))),1,as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
output_flag="all" # "lite" or "all"; isn't actually doing anything at this point -- we're outputting summary info at interval output_thin, and full pop info at an interval hard-coded in (I know, I know, but I don't want to add another parameter right now)
output_thin=1

# Experiment info
experiment_folder <- ""  #contains info on batch of simulations
experimentID <- 1
notes <- "test"   #can include a character string here with notes on the sim

# bio parameters
mutation_type <- "frequent"
mu <- 0.01 # for plasticity
mu_theta <- 1 # in km
nav_rad <- 0.05 # navigation radius, in km
adult_survival_prob <- 0
base_fecund <- 1000
habID <- 9 # this is a small habitat (population 200), so sims run fast

# sim parameters
theta_start_min <- 0.005
theta_start_max <- 140
p_start_min <- -10
p_start_max <- 10
nsteps <- 300

# sim options
normalize_offspring=FALSE
plasticity_on="multiplicative"
larval_output_by_theta=FALSE # save generic data about the relationship between theta and larval output

# disturbance parameters
Dp=0 # probability of disturbance (per timestep)
De=0.25 # extent of disturbance
Dl=5 # duration of disturbance
disturb_method="fractal" # "fractal" or "circle"
Dm=0.1 # magnitude of disturbance (disturbed anemones' fecundity multiplied by this amt; smaller is more impactful)

simID <- SimFn(repID=repID,
               experimentID=experimentID,
               output_flag=output_flag,
               output_thin=output_thin,
               experiment_folder=experiment_folder,
               notes=notes,
               
               # bio parameters
               mutation_type=mutation_type,
               mu=mu,
               mu_theta=mu_theta,
               nav_rad=nav_rad,
               adult_survival_prob=adult_survival_prob,
               base_fecund=base_fecund,
               habID=habID,
               
               # sim parameters
               theta_start_min=theta_start_min,
               theta_start_max=theta_start_max,
               p_start_min=p_start_min,
               p_start_max=p_start_max,
               nsteps=nsteps,
               
               # sim options
               normalize_offspring=normalize_offspring,
               plasticity_on=plasticity_on,
               larval_output_by_theta=larval_output_by_theta,
               
               # disturbance parameters
               Dp=Dp,
               De=De,
               Dl=Dl,
               disturb_method=disturb_method,
               Dm=Dm)

#### Some output plots ####
library(ggplot2)
library(colorspace)
library(gridExtra)
df_summary <- read.csv(paste0("output/raw/",simID,"_summary.csv"))
df_all <- if(file.exists(paste0("output/raw/",simID,"_all.csv"))) read.csv(paste0("output/raw/",simID,"_all.csv")) else NULL

#### 1. Dynamics of key variables over time
ggplot(df_summary,aes(x=t_i,y=median))+
  geom_ribbon(aes(ymin=q25,ymax=q75),alpha=0.2)+
  geom_line()+
  facet_wrap(vars(metric),scales='free')

#### 2. Map of the population at the end
experiment_i <- read.csv("maps/index_habs.csv") |> dplyr::filter(hab_id==habID)
qrast <- rast(paste0("maps/b",experiment_i$basemap_id,"/qmap_b",
                     experiment_i$basemap_id,"_q",experiment_i$qmap_id,".tif"))
load(paste0("maps/habfiles/hab_",experiment_i$hab_id,".RData"))
qmap <- ggplot()+
  ggspatial::layer_spatial(qrast$q)+ #aggregate(qrast$q,2): doesn't work great for kimbe, but okay for others
  scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
  #geom_sf(data=hab_params$sfc_patches,size=0.5,color='red')+ # include anemones
  ggspatial::annotation_scale()+
  #    labs(title=paste0("Sim ", simID_i, "\n",gsub("Export full, ","",experiment_i$notes)))+
  theme(legend.position = "bottom",axis.text=element_blank(),axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5))
load(paste0("output/raw/",simID,"_popsnapshot.RData"))
pop_df <- left_join(pop_df,patch_locations,by=join_by("patch"=="id"))
pop_df$eff_theta <- pop_df$theta*exp(2*pop_df$p*(pop_df$q-0.5))
pop_df$eff_theta <- pmax(pmin(pop_df$eff_theta,140),0.005)
pmap <- ggplot(pop_df,aes(x=x,y=y,color=p))+
  geom_point(alpha=0.5)+
  labs(title=paste0("t = ",t_i),x=NULL,y=NULL)+
  coord_fixed()+
  theme(legend.position="top",axis.text=element_blank(),axis.ticks = element_blank())+
  scale_color_continuous_diverging(palette="Purple-Brown",mid=0)
thetamap <- ggplot(pop_df,aes(x=x,y=y,color=theta))+
  geom_point(alpha=0.5)+
  labs(title=paste0("t = ",t_i),x=NULL,y=NULL)+
  coord_fixed()+
  theme(legend.position="top",axis.text=element_blank(),axis.ticks = element_blank())+
  scale_color_continuous(palette="Blues")
effthetamap <- ggplot(pop_df,aes(x=x,y=y,color=eff_theta))+
  geom_point(alpha=0.5)+
  labs(title=paste0("t = ",t_i),x=NULL,y=NULL)+
  coord_fixed()+
  theme(legend.position="top",axis.text=element_blank(),axis.ticks = element_blank())+
  scale_color_continuous(palette="Oranges",trans="reverse")+
  guides(color = guide_colorbar(reverse = TRUE))
q_v_theta <- ggplot(pop_df,aes(x=q,y=theta))+
  geom_point(alpha=0.5)+
  geom_smooth()+
  labs(title=paste0("Corr(q, theta) = ",round(cor(pop_df$q,pop_df$theta),2)))
q_v_efftheta <- ggplot(pop_df,aes(x=q,y=eff_theta))+
  geom_point(alpha=0.5)+
  geom_smooth()+
  labs(title=paste0("Corr(q, efftheta) = ",round(cor(pop_df$q,pop_df$eff_theta),2)))
q_v_p <- ggplot(pop_df,aes(x=q,y=p))+
  geom_point(alpha=0.5)+
  geom_smooth()+
  labs(title=paste0("Corr(q, p) = ",round(cor(pop_df$q,pop_df$p),2)))
grid.arrange(qmap,pmap,thetamap,effthetamap,q_v_p,q_v_theta,q_v_efftheta,
             layout_matrix=rbind(c(1,1,2,3,4),
                                 c(1,1,2,3,4),
                                 c(1,1,5,6,7)),
             top=paste("Last timestep, ",simID))
