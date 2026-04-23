source('0_Setup.R')
experiment_folder <- "experiments/Exp1_20260413"
basemap_folder <- "basemap_5x5"
seascapeset_folder <- "set1" # in case we do this for a basemap and don't like the outcome, everything's in its own folder

load(file=paste0(experiment_folder,"/",basemap_folder,"/pop_density800.RData")) 

v_width <- 10
v_cutoff <- 5000

# From a basemap:
# Decide what autocorr ranges you'd like to have.
# Pick a pretty closely-sampled range of h values across the entire meaningful range.
# Simulate some number (100?) landscapes for each, and get the median autocorr range.
# Pick the h with the median closest to each of your target autocorr ranges.
# Use the middle 10 (or however many you want) seascapes from those h values.

v_ranges <- c(18,64,223,780,2700) # vector of the target scales of spatial autocorrelation (range of semivariogram)
v_h <- seq(from=-1,to=1.98,by=0.05)
# done_ones <- filter(df_dists,rep_i==1)$simID
# still_need <- (1:60)[(1:60) %notin% done_ones]
# still_need <- c(8, 12, 13, 15, 17, 18, 19, 20, 22, 23, 24, 25, 27, 28, 29, 30, 33, 34,
#                 35, 38, 39, 40, 42, 45, 47, 49, 50, 52, 54, 55, 57, 59, 60)

# # ## index of sims
# index_runs <- expand.grid(h=v_h,rep_i=1:100) %>%
#   rownames_to_column(var="simID") %>%
#   mutate(simID=as.integer(simID)) %>%
#   mutate(target_dist="identity")
# write.csv(x = index_runs,file = paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/index_runs.csv"),row.names = FALSE)

# ## file to hold info on each sim (just create header row, and append rows as we go)
# df_dists <- data.frame(simID=numeric(),
#                        range=numeric(),sill=numeric(),SSErr=numeric(),
#                        model=factor())
# write.table(df_dists, file = paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/df_dists.csv"), sep = ",", append = FALSE,
#             quote = FALSE, col.names = TRUE, row.names = FALSE)

######### Do this part with a separate slurm task for each value of h ##############

## pick current value of h (based on slurm_i), and get index of sims to do
slurm_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# create empty file to write results to
df_dists <- data.frame(simID=numeric(),range=numeric(),sill=numeric(),SSErr=numeric(),model=factor())
write.table(df_dists, file = paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/df_dists_",slurm_i,".csv"), 
            sep = ",", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)

h_i <- v_h[slurm_i] # value of h for this task
print(paste("h_i =",h_i))

index_runs_h <- read.csv(paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/index_runs.csv")) %>%
  filter(abs(h-h_i)<0.001) # the equality statement gives bugs because of rounding

## for each sim, generate the qmap, save it, fit the variogram, and add the variogram stats to df_dists
## maybe we'll make this an mclapply, so all 100 reps can run simultaneously
for(current_sim in 1:nrow(index_runs_h)){
  simID_i <- index_runs_h$simID[current_sim] # ID for this sim in the index of all sims for this basemap
  target_dist <- index_runs_h$target_dist[current_sim]
  
  basemap_file <- paste0(experiment_folder,"/",basemap_folder)
  qmap_file <- paste0(seascapeset_folder,"/sim_",simID_i)
  
  q_rast <- f_GenerateHabQual(base_rast=basemap_file,q_autocorr=h_i,target_dist=target_dist,
                              plot_flag=FALSE,qmap_file = qmap_file)$q_rast
  
  # get variogram info
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  # empirical variogram
  vgm1 <- variogram(q~1,data=spdf1,cressie=TRUE,width=v_width,cutoff=v_cutoff)
  # run gaussian and spherical models, and pick the better one
  vgmf <- fit.variogram(vgm1,vgm(c("Gau","Sph","Exp")))
  
  # store output
  df_all <- data.frame(simID=simID_i,
                       range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
                       model=vgmf$model[2])
  write.table(df_all, file = paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/df_dists_",slurm_i,".csv"),
              sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
} #rep_i


# # ############ Pick which qmaps to use ############
# # Pick the h with the median closest to each of your target autocorr ranges.
# # Use the middle 10 (or however many you want) seascapes from those h values.
# # For each chosen qmap, generate a hab_params, and save it in the basemap folder

# # First, combine all the csv files
# files <- list.files(path = paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder),
#                     pattern = "df_dists*", full.names = TRUE)
# list_of_dfs <- lapply(files, read.csv)
# df_dists <- do.call(rbind, list_of_dfs)
# 
# v_ranges <- c(18,64,223,780,2700) # vector of the target scales of spatial autocorrelation (range of semivariogram)
# index_runs <- read.csv(paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/index_runs.csv"))
# df_dists <- read.csv(paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/df_dists.csv")) %>%
#   mutate(simID=as.integer(simID)) %>%
#   left_join(index_runs,by="simID") 
# %>%
#   filter(model!="Sph")
# 
# ggplot(df_dists,aes(x=h,y=range,color=model))+
#   geom_point()+
#   lims(y=c(0,5000))
# # median by h
# median_by_h <- df_dists %>% group_by(h) %>% summarize(median(range))
# plot(median_by_h$h,median_by_h$`median(range)`)
