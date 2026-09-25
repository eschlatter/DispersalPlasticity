if(version$major=="4" & version$minor=="4.0") .libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old") else .libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
source("New/functions.R")
library(terra)
library(sf)
library(data.table)
library(dplyr)
experiment_folder <- "New/Maps"

#### Basemap 1: 5x5 full ####
basemapID=1
popmapID=2

# Variables
#### Spatial variation
v_qmaps <- seq(from=219, to=417,by=2)
df_qmaps <- data.frame(qmap_id=v_qmaps,hab_id=NA)
# index_habs <- read.csv("New/Maps/index_habs.csv")
for(q_i in 1:length(v_qmaps)){
  qmapID <- v_qmaps[q_i]
  # hab_i <- filter(index_habs,basemap_id==basemapID,popmap_id==popmapID,qmap_id==qmapID)
  habID <- f_MakeHabitat(qmapID=qmapID,popmapID = popmapID,experiment_folder = experiment_folder)
  df_qmaps$hab_id[q_i] <- habID
}

#### Temporal variation (disturbance)
#     disturbance extent (fraction of habitat)
v_De <- c(0.2,0.4,0.6,0.8) 
mid_De <- 0.4
#     disturbance duration
v_Dl <- c(1,5,10,15)  
mid_Dl <- 5
#     disturbance magnitude (factor to reduce fecundity)
v_Dm <- c(0.05,0.1,0.2,0.4)  
mid_Dm <- 0.1
df_disturbances <- data.frame(De=c(0,v_De,rep(mid_De,8)),
                              Dl=c(0,rep(mid_Dl,4),v_Dl,rep(mid_Dl,4)),
                              Dm=c(1,rep(mid_Dm,8),v_Dm),
                              Dp=c(0,rep(0.1,12)),
                              De_seq=c(1,rep(1,4),rep(0,8)),
                              Dl_seq=c(1,rep(0,4),rep(1,4),rep(0,4)),
                              Dm_seq=c(1,rep(0,8),rep(1,4))) |>
  group_by(De,Dl,Dm,Dp) |>
  # columns to indicate whether each row is part of the univariate sequence for each variable
  summarize(De_seq=sum(De_seq),Dl_seq=sum(Dl_seq),Dm_seq=sum(Dm_seq)) |>
  # start with all the De's, so I have something to look at if things take a long time
  arrange(-De_seq) |>
  # give a disturb_id
  tibble::rownames_to_column(var="disturb_id")
v_disturb <- df_disturbances$disturb_id

df_experiments <- expand.grid(qmap_id=v_qmaps,disturb_id=v_disturb) |>
  mutate(basemap_id=basemapID,popmap_id=popmapID) |>
  left_join(df_disturbances) |>
  left_join(df_qmaps)
df_experiments$experiment_id <- 1:nrow(df_experiments)

save(df_experiments,df_qmaps,file=paste0("New/df_experiments_b",basemapID,".RData"))


## choose a few qmaps to use
qmap_index <- read.csv("New/Maps/index_qmaps.csv")
df_qmaps <- left_join(df_qmaps,qmap_index,by=join_by(qmap_id))
df_qmaps_sub <- slice_sample(df_qmaps,by=h,n=5)
habs_to_use <- df_qmaps_sub$hab_id

save(habs_to_use,file="New/habs_to_use_090926.RData")


##### Some mutation rate experiments ########
basemapID=1
popmapID=2
load("New/df_experiments_b1.RData")
load("New/habs_to_use_090926.RData")
qmap_index <- read.csv("New/Maps/index_qmaps.csv")

# Start with the qmaps we're using for the larger experiment (9/9/26), then reduce further
# Two of each h-value x 4 h-values = 8 maps
df_qmaps <- filter(df_qmaps,hab_id %in% habs_to_use) |>
  left_join(qmap_index,by=join_by(qmap_id)) |>
  slice_sample(by=h,n=2)
v_qmaps <- df_qmaps$qmap_id

# And we're just going to do two disturbance regimes: no disturbance and the intermediate value for each param
mid_De <- 0.4 #     disturbance extent (fraction of habitat)
mid_Dl <- 5 #     disturbance duration
mid_Dm <- 0.1 #     disturbance magnitude (factor to reduce fecundity)
df_disturbances <- data.frame(disturb_id=1:2,
                              De=c(0,mid_De),
                              Dl=c(0,mid_Dl),
                              Dm=c(1,mid_Dm),
                              Dp=c(0,0.1),
                              De_seq=c(1,1),
                              Dl_seq=c(1,1),
                              Dm_seq=c(1,1))
v_disturb <- df_disturbances$disturb_id

# Now let's get some mutation magnitudes
v_mu_theta <- exp(seq(from=log(0.001),to=log(1),length.out=5))
v_mu_p <- exp(seq(from=log(0.001),to=log(1),length.out=5))
df_experiments_freq <- expand.grid(qmap_id=v_qmaps,disturb_id=v_disturb,mu_theta=v_mu_theta,mu_p=v_mu_p) |>
  mutate(basemap_id=basemapID,popmap_id=popmapID,mutation_type="frequent") |>
  left_join(df_disturbances) |>
  left_join(df_qmaps)

# And we'll do a few more mutation experiments, using the "rare" method of mutation. mu_p and mu_theta refer to the frequency of mutations. Magnitudes are ~N(0,1).
v_mu_theta <- c(0.001,0.01,0.1)
v_mu_p <- c(0.001,0.01,0.1)
df_experiments_rare <- expand.grid(qmap_id=v_qmaps,disturb_id=v_disturb,mu_theta=v_mu_theta,mu_p=v_mu_p) |>
  mutate(basemap_id=basemapID,popmap_id=popmapID,mutation_type="rare") |>
  left_join(df_disturbances) |>
  left_join(df_qmaps)

df_experiments <- rbind(df_experiments_freq,df_experiments_rare)
df_experiments$experiment_id <- 1:nrow(df_experiments)

save(df_experiments,df_qmaps,file=paste0("New/df_experiments_b1_mutation.RData"))

par(mfrow=c(4,2))
for(qmap_id in df_qmaps$qmap_id){
  qrast <- rast(paste0("New/Maps/b1/qmap_b1_q",qmap_id,".tif"))
  plot(qrast$q,main=qmap_id)
}


par(mfrow=c(4,5))
for(hab_i in habs_to_use){
  row_i <- df_qmaps[df_qmaps$hab_id==hab_i,]
  qrast <- rast(paste0("New/Maps/b1/qmap_b1_q",row_i$qmap_id,".tif"))
  plot(qrast$q,main=row_i$qmap_id)
}

##### Old habmaps ######

# Want to run some sims on the medium Kimbe map at low density. That's b4, p6
# We'll use same habitat quality everywhere. That's q13. (Could also use q12, or make some more, but they're slow.)
f_MakeHabitat(qmapID=13,popmapID = 6, experiment_folder = experiment_folder)
# Then run MakeRefMat for popmapID=6.

# Also want to run some sims on the full-density square map. That's b1, p2.
# Various q options: 1-7
for(q_i in 1:7){
  f_MakeHabitat(qmapID=q_i,popmapID = 2,experiment_folder = experiment_folder)
}
f_MakeHabitat(qmapID=14,popmapID = 2,experiment_folder = experiment_folder)
f_MakeHabitat(qmapID=14,popmapID = 1,experiment_folder = experiment_folder)
# Then run MakeRefMat for popmapID=2.

# Okay, that full-density map will be good later, but it's taking way too long on an interactive node.
# Let's get b1, p1. All the same q's.
for(q_i in 1:7){
  f_MakeHabitat(qmapID=q_i,popmapID = 1,experiment_folder = experiment_folder)
}
# Then run MakeRefMat for popmapID=1.