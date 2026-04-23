source('0_Setup.R')
experiment_folder <- "experiments/Exp1_20260413"
basemap_folder <- "basemap_5x5"
seascapeset_folder <- "set1" # in case we do this for a basemap and don't like the outcome, everything's in its own folder


# # ############ Pick which qmaps to use ############
# # Pick the h with the median closest to each of your target autocorr ranges.
# # Use the middle 10 (or however many you want) seascapes from those h values.
# # For each chosen qmap, generate a hab_params, and save it in the basemap folder

# First, combine all the csv files
files <- list.files(path = paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder),
                    pattern = "df_dists*", full.names = TRUE)
list_of_dfs <- lapply(files, read.csv)
df_dists <- do.call(rbind, list_of_dfs)

v_ranges <- c(18,64,223,780,2700) # vector of the target scales of spatial autocorrelation (range of semivariogram)
index_runs <- read.csv(paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/index_runs.csv"))
df_dists <- df_dists %>%
  mutate(simID=as.integer(simID)) %>%
  left_join(index_runs,by="simID") %>%
  filter(model!="Gau")

ggplot(df_dists,aes(x=h,y=range/1000,color=model))+
  geom_point()+
  #lims(x=c(-1,0.5),y=c(0,0.1))+
  lims(y=c(0,50))+
  labs(y="range(km)")
# median by h
median_by_h <- df_dists %>% group_by(h) %>% summarize(med_range_km=median(range)/1000)
plot(median_by_h$h,median_by_h$med_range_km,ylim=c(0,5))

## for 1x25, and for 5x5
df_dists_choose <- df_dists[0,]
h_to_use <- seq(from=-0.4, to=0.8, by=0.2)
for(h_i in h_to_use){
  df_dists_i <- filter(df_dists,abs(h-h_i)<0.001)
  med_h_i <- median_by_h$med_range_km[abs(median_by_h$h-h_i)<0.001]
  df_dists_choose <- rbind(df_dists_choose,slice_min(df_dists_i,range-med_h_i,n=10))
}
save(df_dists_choose,file=paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/df_choose.RData"))


