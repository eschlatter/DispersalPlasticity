.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
experiment_folder <- "New/Maps"
library(foreach)
library(doParallel)
library(sf)
library(terra)
library(calculus)
library(parallel)

#popmapID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
popmapID <- 2
popmap_row <- read.csv(paste0(experiment_folder,"/index_popmaps.csv")) |> 
  dplyr::filter(popmap_id==popmapID)
npatch <- popmap_row$npatch

#### SET UP OBJECTS ####
numCores <- parallelly::availableCores()
print(numCores)

# get patch and map objects
basemapID <- popmap_row$basemap_id
patchdist_file <- paste0(experiment_folder,"/b",basemapID,"/patchdists_b",basemapID,"_p",popmapID,".RData")
output_file <- paste0(experiment_folder,"/b",basemapID,"/PopmapObjects_b",basemapID,"_p",popmapID,"_new.RData")
load(patchdist_file)
patch_dists <- collapse::qM(patch_dists) # remove units, but it's in km
patch_dists <- patch_dists[1:npatch,1:npatch]

### make a reference matrix of distance bins, theta values, where each site is the kernel evaluated at the distance
# (this is fast and small)
v_dist_bins <- c(seq(from=0,to=1,by=0.001),seq(from=1,to=max(patch_dists)+0.01,by=0.01))
v_theta_bins <- exp(seq(from=log(0.005),to=log(20*max(patch_dists)),length.out=50)) # let theta get much bigger
#v_theta_bins <- c(seq(from=0.005,to=0.095,by=0.005),seq(from=0.1,to=max(patch_dists)/2,by=0.1))

### figure out which row of ref_mat each distance corresponds to
# (returns the lower bound distance)
# this takes some time, and is the same size as patch_dists.
patch_dists <- asplit(patch_dists,2) # split patch_dists into a list of columns
gc()
cluster <- makeForkCluster(nnodes=numCores)
registerDoParallel(cluster)
patch_dists_ref <- foreach(col_i = patch_dists,.combine="cbind") %dopar% {
  vapply(col_i,function(cell_i) match(TRUE,cell_i<v_dist_bins)-1,FUN.VALUE=double(1),USE.NAMES = FALSE)
  # gc()
}
stopCluster(cl = cluster)
storage.mode(patch_dists_ref) <- "integer"
colnames(patch_dists_ref) <- NULL
rownames(patch_dists_ref) <- NULL

rm(patch_dists)

##### offmap corrections, the precise way
# offmap corrections
offmap_corrections <- matrix(nrow=npatch,ncol=length(v_theta_bins))

base_rast <- rast(paste0("New/Maps/b",basemapID,"/base_b",basemapID,".tif"))
map_ext <- terra::ext(base_rast)
map_bounds <- list(x=c(as.numeric(map_ext[c(1,2)])/1000),
                   y=c(as.numeric(map_ext[c(3,4)])/1000))
load(paste0("New/Maps/b",basemapID,"/pop_b",basemapID,"_p",popmap_row$popmap_id,".RData"))
origin_pts <- st_as_sf(data.frame(st_coordinates(sfc_patches)),coords=c("X","Y"),crs=crs(base_rast))

# offmap_corrections (row for each origin anemone, col for each theta val)
offmap_corrections <- mclapply(1:npatch,function(i){
  # get x and y coordinates of anemone
  central_pt <- st_coordinates(origin_pts[i,])
  centr_x <- central_pt[1,1]/1000
  centr_y <- central_pt[1,2]/1000
  offmap_corrections_j <- vector(length=length(v_theta_bins))
  for(j in seq_along(v_theta_bins)){
    th_j <- v_theta_bins[j]
    if(th_j>0.01){ # if theta<0.01, the integral gets weird, and we're going to assume nobody falls off the map. If there's a little edge effect at those tiny values, I don't mind anyway.
      vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                            bounds=map_bounds,
                            params=list(theta=th_j),
                            method="divonne",
                            relTol=0.0001)$value
      offmap_corrections_j[j] <- (2*pi*th_j)/vol_inmap
    } else offmap_corrections_j[j] <- 1 
  } # for j in seq_along(v_theta_bins)
  return(offmap_corrections_j)
},mc.cores=numCores)
offmap_corrections <- do.call(rbind,offmap_corrections)

save(patch_dists_ref,v_dist_bins,v_theta_bins,offmap_corrections,file=paste0(output_file))
