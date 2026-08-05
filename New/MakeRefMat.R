.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
experiment_folder <- "New/Maps"
library(foreach)
library(doParallel)
library(sf)
library(calculus)

popmapID <- 6
popmap_row <- read.csv(paste0(experiment_folder,"/index_popmaps.csv")) |> 
  dplyr::filter(popmap_id==popmapID)
npatch <- popmap_row$npatch
calc_offmap_corr <- TRUE

#### SET UP OBJECTS ####
numCores <- parallelly::availableCores()
print(numCores)

# get patch and map objects
basemapID <- popmap_row$basemap_id
patchdist_file <- paste0(experiment_folder,"/b",basemapID,"/patchdists_b",basemapID,"_p",popmapID,".RData")
reefsf_file <- paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".RData")
output_file <- paste0(experiment_folder,"/b",basemapID,"/RefMat_b",basemapID,"_p",popmapID,".RData")
load(reefsf_file)
load(patchdist_file)
patch_dists <- collapse::qM(patch_dists) # remove units, but it's in km
patch_dists <- patch_dists[1:npatch,1:npatch]

# make a reference matrix of distance bins, theta values, where each site is the kernel evaluated at the distance
# (this is fast and small)
v_dist_bins <- c(seq(from=0,to=1,by=0.001),seq(from=1,to=max(patch_dists)+0.01,by=0.01))
v_theta_bins <- c(seq(from=0.005,to=0.095,by=0.005),seq(from=0.1,to=max(patch_dists)/2,by=0.1))
ref_mat <- matrix(nrow=length(v_dist_bins),ncol=length(v_theta_bins))
for(th_i in seq_along(v_theta_bins)){
  ref_mat[,th_i] <- (1/(2*pi*v_theta_bins[th_i]))*dgamma(v_dist_bins,shape=1,scale=v_theta_bins[th_i])
}

# offmap corrections
if(calc_offmap_corr==TRUE){
  map_ext <- terra::ext(reef_sf)
  map_bounds <- list(x=c(as.numeric(map_ext[c(1,2)])/1000),
                     y=c(as.numeric(map_ext[c(3,4)])/1000))
  central_pt <- st_coordinates(st_centroid(reef_sf))
  centr_x <- mean(map_bounds$x)
  centr_y <- mean(map_bounds$y)
  offmap_correction <- vector(length=length(v_theta_bins)) # object to hold the offmap correction for each theta
  
  # could parallelize this, but is it worth it?
  for(i in seq_along(v_theta_bins)){
    th_i <- v_theta_bins[i]
    if(th_i>0.01){
      vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                            bounds=map_bounds,
                            params=list(theta=th_i),
                            method="divonne",
                            relTol=0.0001)$value
      offmap_correction[i] <- (2*pi*th_i)/vol_inmap
    } else offmap_correction[i] <- 1  
  } # for i in seq_along(v_theta_bins)
} else offmap_correction <- NULL # if calc_offmap_corr==TRUE

# figure out which row of ref_mat each distance corresponds to
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
save(patch_dists_ref,ref_mat,v_dist_bins,v_theta_bins,offmap_correction,file=paste0(output_file))
