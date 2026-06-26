## script to run that just generates the connectivity matrices and saves them to global scratch
source('0_Setup.R')
library(calculus)
library(collapse)
experiment_name <- "Exp3_20260528"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index_smallkern.csv"))
load(paste0(experiment_folder,"/basemap_index.RData")) # basemap index

# Create new folders in the project directory and MSI scratch directory
temp_dir <- paste0("connmats_smallkern_",experiment_name)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
if(!dir.exists(temp_path)) dir.create(temp_path)

# identify which map we're on (as indexed by experiment_index)
#run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run_i=5
map_info <- experiment_index[run_i,]
mapID_i <- experiment_index$mapID[run_i]
if(!dir.exists(paste0(temp_path,"/map_",mapID_i))) dir.create(paste0(temp_path,"/map_",mapID_i))

# params
load(paste0(experiment_folder,"/params_",map_info$param_id,".RData"))
list2env(x=params,envir=environment())

# hab_params
load(file=paste0(experiment_folder,"/habfiles/habparams_",mapID_i,".RData")) 
list2env(x=hab_params,envir=environment()) 

# load patch_dists
load(file=paste0(experiment_folder,"/b",map_info$basemap_id,"/patchdists_b",map_info$basemap_id,"_p",map_info$popmap_id,".RData"))
units(patch_dists) <- NULL
units(nav_rad) <- NULL

group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)

# map info
map_ext <- ext(reef_sf)
map_bounds <- list(x=c(as.numeric(map_ext[1])/1000,as.numeric(map_ext[2])/1000),
                   y=c(as.numeric(map_ext[3])/1000,as.numeric(map_ext[4])/1000))
# central point, to use for off-map correction and self-recruitment calculations
central_pt <- st_coordinates(st_centroid(reef_sf))
centr_x <- central_pt[1]/1000
centr_y <- central_pt[2]/1000
square_bounds <- list(x=c(centr_x-nav_rad,centr_x+nav_rad),y=c(centr_y-nav_rad,centr_y+nav_rad))

# overlap discount (multiply destination patch by this factor)
patch_locations$overlap_discount <- 1/colSums(patch_dists<nav_rad)

# index of alpha/theta combinations
eff_kern_index <- matrix(1:(length(v_alphas)*length(v_thetas)),ncol=length(v_alphas),nrow=length(v_thetas))
eff_kern_df <- data.frame(index=as.vector(eff_kern_index),
                          alpha_ind=rep(1:length(v_alphas),each=length(v_thetas)),
                          theta_ind=rep(1:length(v_thetas),times=length(v_alphas)))
eff_kern_df$alpha_val <- v_alphas[eff_kern_df$alpha_ind]
eff_kern_df$theta_val <- v_thetas[eff_kern_df$theta_ind]

# make the off-the-map corrections and self-recruitment calculations
for(i in 1:nrow(eff_kern_df)){
  al_i <- eff_kern_df$alpha_val[i]
  th_i <- eff_kern_df$theta_val[i]
  
  # off-the-map correction
  vol_inmap <- max(integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                            bounds=map_bounds,
                            params=list(theta=th_i),
                            relTol=0.00001)$value,
                   2*pi*th_i)
  eff_kern_df$offmap_corr[i] <- 1+(2*pi*th_i-vol_inmap)/vol_inmap
  
  # self-recruitment values
  vol_selfrecruit <- integral(function(r,theta){(1/th_i)*exp(-r/th_i)},
                              bounds=list(r=c(0,nav_rad),theta=c(0,2*pi)),
                              coordinates="polar",
                              relTol=0.00001,
                              method="divonne"
  )$value
  eff_kern_df$selfrecruit[i] <- vol_selfrecruit
}

######### Make and save connectivity matrices ##########
for(g in 1:nrow(group_index)){
  if(!file.exists(paste0(temp_path,"/map_",mapID_i,"/grp_",g))){
    # get the connectivity matrix among patches, given the group parameter values and patch-level q's
    # (and accounting for the patch population x per capita output b_i from each patch)
    
    # calculate this matrix
    ## This is what uses multiple cores: f_GetPlasticConnMat uses f_GetConnectivityMatrix_parallel
    conn_mat <- f_GetConnMat(g,
                             group_index,patch_locations,v_p,v_alphas,v_thetas,eff_kern_index,eff_kern_df,patch_dists,nav_rad,
                             numCores=parallelly::availableCores())
    # then store it
    write_fst(qDF(conn_mat),paste0(temp_path,"/map_",mapID_i,"/grp_",g),compress=0)
    #saveRDS(object=conn_mat,file=paste0(saveto_folder,"/grp_",g,".rds"),compress=FALSE)
  }
} # g

