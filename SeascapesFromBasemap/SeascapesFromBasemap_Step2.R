source('0_Setup.R')
experiment_folder <- "experiments/Exp2_20260423"
basemap_id <- 2
popmap_id <- 2 # popmap to use for measuring variogram
set_id <- 1
seascape_set <- paste0("set",set_id) # in case we do this for a basemap and don't like the outcome, everything's in its own folder
seascapeset_folder <- paste0(experiment_folder,"/b",basemap_id,"/",seascape_set)


# # ############ Pick which qmaps to use ############
# # Pick the h with the median closest to each of your target autocorr ranges.
# v_ranges <- c(18,64,223,780,2700) # vector of the target scales of spatial autocorrelation (range of semivariogram)
# # Use the middle 10 (or however many you want) seascapes from those h values.
# # For each chosen qmap, generate a hab_params, and save it in the basemap folder

# First, combine all the csv files
files <- list.files(path = seascapeset_folder,
                    pattern = "df_step1_", full.names = TRUE)
list_of_dfs <- lapply(files, read.csv)
df_step1 <- do.call(rbind, list_of_dfs)

index_runs <- read.csv(paste0(seascapeset_folder,"/index_step1.csv"))
df_dists <- df_step1 %>%
  mutate(simID=as.integer(simID)) %>%
  left_join(index_runs,by="simID") %>%
  filter(model!="Sph",range<100000)

ggplot(df_dists,aes(x=h,y=range/1000,color=model))+
  geom_point()+
  lims(x=c(-1,0.5),y=c(0,0.1))+
#  lims(y=c(0,50))+
  labs(y="range(km)")
# median by h
median_by_h <- df_dists %>% group_by(h) %>% summarize(med_range_km=median(range)/1000)
ggplot(median_by_h,aes(x=h,y=med_range_km))+
  geom_point()+
  lims(x=c(-1,0.2),y=c(0,0.1))

# Pick which h values to use, and make a dataframe
df_dists_choose <- df_dists[0,]
# # for basemap1 (5x5km):
# h_to_use <- c(-1,-0.5,-0.1,0.25,0.65)
# for basemap2 (1x25km):
h_to_use <- c(-0.75,-0.25,0.5,1.25,1.95)
#h_to_use <- seq(from=-0.4, to=0.8, by=0.2)
for(h_i in h_to_use){
  df_dists_i <- filter(df_dists,abs(h-h_i)<0.001)
  med_h_i <- median_by_h$med_range_km[abs(median_by_h$h-h_i)<0.001]
  df_dists_choose <- rbind(df_dists_choose,slice_min(df_dists_i,range-med_h_i,n=10))
}

# Look at the choices:
ggplot(df_dists_choose,aes(x=h,y=range))+
  geom_point()
par(mfrow=c(2,5))
for(h_now in h_to_use){
  choose_h <- filter(df_dists_choose,h==h_now)
  for(map_i in 1:10){
    q_rast <- rast(paste0(seascapeset_folder,"/qmap_b2_q",choose_h$simID[map_i],".tif")) # load q_rast 
    plot(q_rast$q)
  }
}


# Save them
save(df_dists_choose,file=paste0(seascapeset_folder,"/df_step2.RData"))


############ Put them all together into an index of sims ############
# This index should have a row for every unique combo of basemap, popmap, qmap, target_dist, and params.
# Might actually do multiple runs per row (with different or the same initial conditions)
source('0_Setup.R')
experiment_folder <- "experiments/Exp2_20260423"
basemaps <- c(1,2) # set of basemaps to include
sets <- c(1,1) # goes along with the basemaps vector: ith element of sets should be the set used by the ith element of basemaps
seascape_set <- paste0("set",sets)
target_dists <- c("identity","C","D","E")

# put qmaps from all the basemaps together
experiment_index <- data.frame(simID=integer(),range=numeric(),sill=numeric(),SSErr=numeric(),model=character(),
                               h=numeric(),rep_i=integer(),target_dist=character(),
                               basemap_id=integer(),popmap_id=integer(),param_id=integer())
for(i in 1:length(basemaps)){
  basemap_id <- basemaps[i]
  set_name <- seascape_set[i]
  load(paste0(experiment_folder,"/b",basemap_id,"/",set_name,"/df_step2.RData"))
  experiment_index <- rbind(experiment_index,df_dists_choose)
}

# multiply out by target_dists
n_maps <- nrow(experiment_index)
experiment_index <- slice(experiment_index,rep(1:n_maps,each=length(target_dists))) %>%
  mutate(target_dist = rep(target_dists,times=n_maps)) %>%
  rename(q_ID=simID) %>%
  mutate(mapID=row_number(),range_ident=range) %>%
  dplyr::select(-rep_i) # don't need this anymore, and it'll be confusing with reps later
experiment_index[experiment_index$target_dist!="identity",c("range","sill","SSErr","model")] <- NA

# (we could multiply out by param set here, if we needed)

write.csv(experiment_index,file=paste0(experiment_folder,"/experiment_index.csv"),row.names=FALSE)

######### Also make a basemap index
basemap_index <- data.frame(basemap_id=1:2,set_id=1,variog_width=c(10,10),variog_cutoff=c(5000,10000))
save(basemap_index,file=paste0(experiment_folder,"/basemap_index.RData"))
