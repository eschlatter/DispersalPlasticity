.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
library(sf)
library(terra)
#library(raster)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggspatial)
library(tibble)
library(parallel)
source("New/functions.R")


# From a basemap:
# Decide what autocorr ranges you'd like to have.
# Pick a pretty closely-sampled range of h values across the entire meaningful range.
# Simulate some number (100?) landscapes for each, and get the median autocorr range.
# Pick the h with the median closest to each of your target autocorr ranges.
# Use the middle 10 (or however many you want) seascapes from those h values.

#source('0_Setup.R')
experiment_folder <- "New/Maps"
basemap_id <- 4
popmap_id <- 7 # popmap to use for measuring variogram
set_id <- 1
seascape_set <- paste0("set",set_id) # in case we do this for a basemap and don't like the outcome, everything's in its own folder
seascapeset_folder <- paste0(experiment_folder,"/Kimbe/",seascape_set)
numCores <- parallelly::availableCores()

# load data collection points for fitting variogram
load(file=paste0(experiment_folder,"/b",basemap_id,"/pop_b",basemap_id,"_p",popmap_id,".RData"))
spdf1 <- as_Spatial(sfc_patches)

v_h <- seq(from=-1.98,to=1.98,by=0.2) # set of h values to try
# variogram fitting parameters
variog_width <- 10
variog_cutoff <- 20000

#### Do this part once per seascape set: ######

# ## index of sims
# index_step1 <- expand.grid(h=v_h,rep_i=1:100) %>%
#   rownames_to_column(var="simID") %>%
#   mutate(simID=as.integer(simID)) %>%
#   mutate(target_dist="A",basemap_id=basemap_id,popmap_id=popmap_id) # target dist (A) is uniform
# write.csv(x = index_step1,file = paste0(seascapeset_folder,"/index_step1.csv"),row.names = FALSE)
# ## file to hold info on each sim (just create header row, and append rows as we go)
# df_step1 <- data.frame(simID=numeric(),
#                        range=numeric(),sill=numeric(),SSErr=numeric(),
#                        model=factor())
# write.table(df_step1, file = paste0(seascapeset_folder,"/df_step1.csv"), sep = ",", append = FALSE,
#             quote = FALSE, col.names = TRUE, row.names = FALSE)

######### Do this part with a separate slurm task for each value of h ##############

## pick current value of h (based on slurm_i), and get index of sims to do
slurm_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(slurm_i)) slurm_i=1

# # create empty file to write results to
# df_step1 <- data.frame(simID=numeric(),range=numeric(),sill=numeric(),SSErr=numeric(),model=factor())
# write.table(df_step1, file = paste0(seascapeset_folder,"/df_step1_",slurm_i,".csv"), 
#             sep = ",", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)

h_i <- v_h[slurm_i] # value of h for this task
print(paste("h_i =",h_i))

index_runs_h <- read.csv(paste0(seascapeset_folder,"/index_step1.csv")) %>%
  filter(abs(h-h_i)<0.001) # the equality statement gives bugs because of rounding

## for each sim, generate the qmap, save it, fit the variogram, and add the variogram stats to df_dists
## maybe we'll make this an mclapply, so all 100 reps can run simultaneously
df_list <- mclapply(1:nrow(index_runs_h),function(current_sim){
  simID_i <- index_runs_h$simID[current_sim] # ID for this sim in the index of all sims for this basemap
  target_dist <- index_runs_h$target_dist[current_sim]


  q_out <- f_GenerateHabQual(q_autocorr=h_i,target_dist=target_dist,plot_flag=FALSE,
                             experiment_folder=experiment_folder,basemap_id=basemap_id,qmap_id=NULL,
                             variog_fit=TRUE, variog_width=variog_width,variog_cutoff=variog_cutoff)
  q_rast <- q_out$q_rast
  vgmf <- q_out$vgm_fit

  # save the map
  writeRaster(q_rast,filename=paste0(seascapeset_folder,"/qmap_b",basemap_id,"_q",simID_i,".tif"),overwrite=TRUE)

  # get variogram info
  spdf1$q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  # empirical variogram
  vgm1 <- gstat::variogram(q~1,data=spdf1,cressie=TRUE,width=variog_width,cutoff=variog_cutoff)
  # run gaussian and spherical models, and pick the better one
  vgmf <- gstat::fit.variogram(vgm1,gstat::vgm(c("Gau","Sph","Exp")))


  # store output
  df_single <- data.frame(simID=simID_i,
                          range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
                          model=vgmf$model[2])
},
mc.cores = numCores) #mclapply

df_all <- do.call(rbind,df_list)
write.table(df_all, file = paste0(seascapeset_folder,"/df_step1_",slurm_i,".csv"),
            sep = ",", append = FALSE, quote = FALSE, row.names = FALSE)

####### Didn't give the slurm task quite enough time, so now I have to go back and calculate the variogram stats again. Oops. ####
all_files <- list.files(seascapeset_folder)
df_files <- all_files[grep("df_step",all_files)]
qmap_files <- all_files[grep("qmap_b4",all_files)]

# First, make a data frame of all the ones whose variogram stats we have
df_step1 <- read.csv(paste0(seascapeset_folder,"/",df_files[1]))
for(i in 2:length(df_files)){
  df_file_i <- read.csv(paste0(seascapeset_folder,"/",df_files[i]))
  df_step1 <- rbind(df_step1,df_file_i)
}
df_step1 <- filter(df_step1,model!="")
write.table(df_step1, file = paste0(seascapeset_folder,"/df_step1.csv"),
            sep = ",", append = FALSE, quote = FALSE, row.names = FALSE)

# Then, list all the qmaps in the folder
df_qmaps <- data.frame(file=qmap_files) %>%
  tidyr::separate_wider_delim(file,"q",names=c("throw","throw2","keep")) %>%
  mutate(qmap_id=as.numeric(gsub("[^0-9]", "", keep))) %>%
  dplyr::select(qmap_id) %>%
  mutate(done=qmap_id %in% as.numeric(df_step1$simID))

qmaps_not_done <- df_qmaps$qmap_id[df_qmaps$done==FALSE]
for(qmap_i in qmaps_not_done){
  q_rast <- rast(paste0(seascapeset_folder,"/qmap_b4_q",qmap_i,".tif"))
  
  # get variogram info
  spdf1$q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  # empirical variogram
  vgm1 <- gstat::variogram(q~1,data=spdf1,cressie=TRUE,width=variog_width,cutoff=variog_cutoff)
  # run gaussian and spherical models, and pick the better one
  vgmf <- gstat::fit.variogram(vgm1,gstat::vgm(c("Gau","Sph","Exp")))

  df_single <- data.frame(simID=qmap_i,
                          range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
                          model=vgmf$model[2])
  write.table(df_single, file = paste0(seascapeset_folder,"/df_step1.csv"),
              sep = ",", append = TRUE, quote = FALSE, col.names=FALSE, row.names = FALSE)
  
}

############ Pick which qmaps to use ############
df_step1 <- read.csv(paste0(seascapeset_folder,"/df_step1.csv"))
index_step1 <- read.csv(paste0(seascapeset_folder,"/index_step1.csv"))
df_dists <- df_step1 %>%
  mutate(simID=as.integer(simID)) %>%
  left_join(index_step1,by="simID") %>%
  filter(model!="Sph",range<50000)

ggplot(df_dists,aes(x=h,y=range/1000,color=model))+
  geom_point()+
  labs(y="range(km)")

median_by_h <- df_dists %>% group_by(h) %>% summarize(med_range_km=median(range)/1000)
ggplot(median_by_h,aes(x=h,y=med_range_km))+
  geom_point()

# Pick which h values to use, and make a dataframe
df_dists_choose <- df_dists[0,]
n_ranges <- 4
range_targets <- seq(from=min(median_by_h$med_range_km),to=max(median_by_h$med_range_km),length.out=4)
h_to_use <- vector(length=n_ranges)
for(i in 1:n_ranges){
  h_to_use[i] <- median_by_h$h[which.min(abs(median_by_h$med_range_km-range_targets[i]))]
}
n_maps_per_range <- 25
for(h_i in h_to_use){
  df_dists_i <- filter(df_dists,abs(h-h_i)<0.001)
  med_h_i <- median_by_h$med_range_km[abs(median_by_h$h-h_i)<0.001]
  df_dists_choose <- rbind(df_dists_choose,slice_min(df_dists_i,abs(range/1000-med_h_i),n=n_maps_per_range))
}

# Look at the choices:
ggplot(df_dists_choose,aes(x=h,y=range))+
  geom_point()+
  labs(title="Spatial Autocorrelation Scales, Kimbe maps",y="Scale of autocorrelation (m)",x="Fractal landscape algorithm aggregation parameter")
ggsave("New/Maps/Kimbe/AutocorrScalePlot.png",width=1500,height = 1500,units="px")
save(df_dists_choose,file="New/Maps/Kimbe/index.RData") # save the index

## plot all the maps.
for(h_now in h_to_use){
  par(mfrow=c(5,5))
  choose_h <- filter(df_dists_choose,h==h_now)
  for(map_i in 1:n_maps_per_range){
    q_rast <- rast(paste0(seascapeset_folder,"/qmap_b4_q",choose_h$simID[map_i],".tif")) # load q_rast
    plot(q_rast$q,main=paste0(choose_h$simID[map_i],",A=",round(choose_h$range[map_i]),"m"),legend=FALSE,cex.main=0.75)
  }
}

index_start <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
df_dists_choose$qmap_id <- index_start:(index_start+nrow(df_dists_choose)-1)
save(df_dists_choose,file="New/Maps/Kimbe/index.RData") # The index now contains the new and old map numbers
index_qmaps_add <- df_dists_choose %>%
  mutate(description=paste0("Kimbe, ",
                            recode_values(h,from=sort(unique(df_dists_choose$h)),
                                          to=c("very low","low","moderate","high"))," autocorrelation")) %>%
  dplyr::select(qmap_id,basemap_id,description,h,target_dist,range,sill,SSErr,model) %>%
  mutate(original_range=range)
index_qmaps_add$target_dist <- "uniform"
fwrite(index_qmaps_add,file="New/Maps/index_qmaps.csv",append=TRUE)
# Save the .tif for each map
# It's the same map in the Kimbe folder and the q4 folder.
# This is redundant, but I want the file structure to match the other basemaps, and they have an extra copy because they needed scale alteration.
for(i in 1:nrow(df_dists_choose)){
  file.copy(from=paste0(seascapeset_folder,"/qmap_b4_q",df_dists_choose$simID[i],".tif"),to=paste0("New/Maps/Kimbe/q",df_dists_choose$simID[i],".tif"))
  file.copy(from=paste0(seascapeset_folder,"/qmap_b4_q",df_dists_choose$simID[i],".tif"),to=paste0("New/Maps/b4/qmap_b4_q",df_dists_choose$qmap_id[i],".tif"))
}

