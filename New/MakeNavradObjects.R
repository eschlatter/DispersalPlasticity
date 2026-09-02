# For a population and nav_rad:
# Take the existing RefMat and convert the kernel values to kernel integrals,
#   by some combination of multiplying by pi*nav_rad^2 and actually computing integrals

.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
experiment_folder <- "New/Maps"
library(foreach)
library(doParallel)
library(sf)
library(terra)
library(calculus)
library(parallel)

nav_rad <- 0.05
popmapID <- 2
popmap_row <- read.csv(paste0(experiment_folder,"/index_popmaps.csv")) |> 
  dplyr::filter(popmap_id==popmapID)
npatch <- popmap_row$npatch

#### SET UP OBJECTS ####
numCores <- parallelly::availableCores()
print(numCores)

# get patch and map objects
basemapID <- popmap_row$basemap_id
load(paste0("New/Maps/b",basemapID,"/PopmapObjects_b",basemapID,"_p",popmapID,"_new.RData"))

# what if we just do them all?
ref_mat_new <- mclapply(seq_along(v_theta_bins),function(th_i){
  ref_mat_th <- vector(length=length(v_dist_bins))
  theta_i <- v_theta_bins[th_i]
  dists_to_do <- seq_along(v_dist_bins)
  for(d_i in dists_to_do){
    dist_i <- v_dist_bins[d_i]
    square_bounds <- list(x=c(dist_i-nav_rad,dist_i+nav_rad),y=c(-nav_rad,nav_rad))
    # calculate the volume under the integral in a square centered at dist with side length 2*nav_rad
    vol_square <- integral(function(x,y,theta){(1/(2*pi*theta))*(1/theta)*exp(-sqrt(x^2+y^2)/theta)},
                           bounds=square_bounds,
                           params=list(theta=theta_i),
                           method="divonne",
                           relTol=0.0001)$value
    vol_circ <- (pi/4)*vol_square
    ref_mat_th[d_i] <- vol_circ
  }
  return(ref_mat_th)
},mc.cores=numCores)
ref_mat <- do.call(cbind,ref_mat_new)
save(ref_mat,file=paste0("New/Maps/b",basemapID,"/NavradObjects_b",basemapID,"_p",popmapID,"_n=",nav_rad,".RData"))

# plot(v_dist_bins,ref_mat_new[,1],type='l')
# plot(v_dist_bins,ref_mat[,1],type='l')


