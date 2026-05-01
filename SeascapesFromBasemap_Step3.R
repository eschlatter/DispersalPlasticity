######## Generate the maps and variogram stats for all target_dists ###########
source('0_Setup.R')
experiment_name <- "Exp2_20260423"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
load(paste0(experiment_folder,"/basemap_index.RData")) # basemap index

# Create new folders in the project directory and MSI scratch directory
temp_dir <- paste0("connmats_",experiment_name)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
if(!dir.exists(temp_path)) dir.create(temp_path)
keep_path <- paste0(experiment_folder,"/habfiles")
if(!dir.exists(keep_path)) dir.create(keep_path)

## Do one slurm task per q_ID
slurm_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#slurm_i=2
qID_i <- unique(experiment_index$q_ID)[slurm_i] # get current value of q_ID
index_i <- filter(experiment_index,q_ID==qID_i)
numCores <- parallelly::availableCores()

variog_data <- data.frame(mapID=integer(),range=numeric(),sill=numeric(),SSErr=numeric(),model=character())

# for each q_ID, process each map:
for(i_map in 1:nrow(index_i)){
  # get the info related to the map
  mapID_i <- index_i$mapID[i_map]
  qID_i <- index_i$q_ID[i_map]
  basemap_i <- index_i$basemap_id[i_map]
  popmap_i <- index_i$popmap_id[i_map]
  param_i <- index_i$param_id[i_map]
  dist_i <- index_i$target_dist[i_map]
  set_i <- basemap_index[basemap_index$basemap_id==basemap_i,"set_id"]
  variog_width_i <- basemap_index[basemap_index$basemap_id==basemap_i,"variog_width"]
  variog_cutoff_i <- basemap_index[basemap_index$basemap_id==basemap_i,"variog_cutoff"]
  seascapeset_folder <- paste0(experiment_folder,"/b",basemap_i,"/set",set_i)
  
  # load the qmap
  q_rast <- rast(paste0(seascapeset_folder,"/qmap_b",basemap_i,"_q",qID_i,".tif"))
  # load popmap
  load(paste0(experiment_folder,"/b",basemap_i,"/pop_b",basemap_i,"_p",popmap_i,".RData"))
  
  ##### transform qmap to the target distribution
  q_rast_new <- q_rast
  q_rast_new$q <- f_TransformDist(matrix(values(q_rast$q),nrow=nrow(q_rast),byrow=TRUE),dist_i)
  q_rast <- q_rast_new # use q_rast from here on
  
  ##### get variogram stats
  # get variogram info
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  # empirical variogram
  vgm1 <- variogram(q~1,data=spdf1,cressie=TRUE,width=variog_width_i,cutoff=variog_cutoff_i)
  # run gaussian and spherical models, and pick the better one
  vgmf <- fit.variogram(vgm1,vgm(c("Gau","Sph","Exp")))
  
  # store output
  variog_data <- rbind(variog_data,data.frame(mapID=mapID_i,range=vgmf$range[2],
                                              sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
                                              model=vgmf$model[2]))
  
  ##### get hab_params, and write to keep_path under the mapID
  load(paste0(experiment_folder,"/params_",param_i,".RData"))
  nav_rad <- params$nav_rad
  units(nav_rad) <- 'km'
  npatch <- length(sfc_patches)
  
  ## put q_rast and K_rast together
  K_rast <- rast(paste0(experiment_folder,"/b",basemap_i,"/pop_b",basemap_i,"_p",popmap_i,".tif"))
  hab_rast <- c(q_rast,K_rast)
  
  ## create df_patches (important: ID should be in the same order as in sfc_patches, or distance matrix will be wrong)
  q_vect <- terra::extract(hab_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)
  K_vect <- terra::extract(hab_rast$K,vect(sfc_patches),xy=TRUE,search_radius=500)
  patch_coords <- st_coordinates(sfc_patches)
  df_patches <- cbind(q_vect[,c("ID","q")],patch_coords)
  df_patches$K <- K_vect$K[df_patches$ID]
  df_patches$id <- df_patches$ID
  df_patches$x <- df_patches$X
  df_patches$y <- df_patches$Y
  df_patches <- df_patches[,c("id","x","y","q","K")]
  df_patches$b <- f_q_to_b(df_patches$q) # calculate reproductive rate (b) from habitat quality (q)
  
  
  ## create patch_angles (save it separately)
  # first load patch_dists
  load(file=paste0(experiment_folder,"/b",basemap_i,"/patchdists_b",basemap_i,"_p",popmap_i,".RData"))
  # create patch_angles
  patch_angles <- (2*nav_rad)/(2*pi*patch_dists)
  patch_angles[is.nan(patch_angles)] <- 1
  patch_angles[drop_units(patch_angles)>1] <- 1
  # save patch_angles
  if(!dir.exists(paste0(temp_path,"/map_",mapID_i))) dir.create(paste0(temp_path,"/map_",mapID_i))
  # write_fst(as.data.frame(patch_angles),paste0(temp_path,"/map_",mapID_i,"/patch_angles"),compress=0)
  saveRDS(patch_angles,file=paste0(temp_path,"/map_",mapID_i,"/patch_angles",".rds"),compress=FALSE) #seems noticeably faster, but let's play around with it
  
  # calculate overlap_discount
  n_neighbors <- rowSums(patch_dists<nav_rad) # number of points within distance nav_rad of focal point (including focal point)
  overlap_discount <- 1/n_neighbors
  
  hab_params <- list(npatch=npatch,
                     patch_locations=df_patches,
                     overlap_discount=overlap_discount,
                     reef_sf=reef_sf,
                     sfc_patches=sfc_patches,
                     hab_type=hab_type,
                     nav_rad=nav_rad,
                     mapID=mapID_i)
  save(hab_params,file=paste0(keep_path,"/habparams_",mapID_i,".RData"))
}

write.csv(variog_data,file=paste0(experiment_folder,"/variog_fit_",slurm_i,".csv"),row.names=FALSE)

### Then, afterward, put all the variog_fit_i.csv's together
# files <- list.files(path = experiment_folder,
#                     pattern = "variog_fit_", full.names = TRUE)
# list_of_dfs <- lapply(files, read.csv)
# variog_fits <- do.call(rbind, list_of_dfs)
# write.csv(variog_fits,"range_vals.csv")
# file.remove(files)