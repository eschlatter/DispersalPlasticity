## script to run that just generates the connectivity matrices and saves them to global scratch
source('0_Setup.R')
library(calculus)
library(collapse)
library(profvis)
experiment_name <- "Exp4_20260611"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
#load(paste0(experiment_folder,"/basemap_index.RData")) # basemap index

calc_offmap_corr=FALSE

# Create new folders in the project directory and MSI scratch directory
temp_dir <- paste0(experiment_name)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
if(!dir.exists(temp_path)) dir.create(temp_path)

# identify which map we're on (as indexed by experiment_index)
#run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run_i=4
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
patch_dists <- collapse::qM(patch_dists)
#units(patch_dists) <- NULL
units(nav_rad) <- NULL

# make distance bins, and find the nearest bin for each dist in patch_dists
# this will return NA's for distance 0
# smalldists <- exp(seq(from=log(min(0.9*patch_dists[patch_dists!=0])),to=log(-10*v_thetas[1]*log(v_thetas[1])),length.out=10))
# dist_list <- data.frame(lowbound=unique(c(smalldists,seq(from=nav_rad,to=max(patch_dists)+nav_rad,by=nav_rad))))
dist_list <- data.frame(lowbound=seq(from=0,to=max(patch_dists)+nav_rad,by=nav_rad))
dist_list$integral_eval <- dist_list$lowbound+nav_rad/2
patch_dists_bin <- cut(patch_dists,breaks=dist_list$lowbound,labels=FALSE)

# group index
group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)

# map info
map_ext <- ext(reef_sf)
map_bounds <- list(x=c(as.numeric(map_ext[1])/1000,as.numeric(map_ext[2])/1000),
                   y=c(as.numeric(map_ext[3])/1000,as.numeric(map_ext[4])/1000))
# central point, to use for off-map correction and self-recruitment calculations
central_pt <- st_coordinates(st_centroid(reef_sf))
centr_x <- mean(map_bounds$x)
centr_y <- mean(map_bounds$y)
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

# hold the integrals at each distance for each param combo
patch_dists_integrals <- matrix(nrow=nrow(dist_list),ncol=nrow(eff_kern_df))

eff_kern_df$offmap_corr <- NA
eff_kern_df$selfrecruit <- NA
# make the off-the-map corrections and self-recruitment calculations
for(i in 1:nrow(eff_kern_df)){
  al_i <- eff_kern_df$alpha_val[i]
  th_i <- eff_kern_df$theta_val[i]
  
  # off-the-map correction
  if(calc_offmap_corr==TRUE & th_i>0.01){
    vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                          bounds=map_bounds,
                          params=list(theta=th_i),
                          method="divonne",
                          relTol=0.0001)$value
    eff_kern_df$offmap_corr[i] <- 1+(2*pi*th_i-vol_inmap)/vol_inmap
  } else eff_kern_df$offmap_corr[i] <- 1
  
  # self-recruitment values
  vol_selfrecruit <- integral(function(r,theta){(1/th_i)*exp(-r/th_i)},
                              bounds=list(r=c(0,nav_rad),theta=c(0,2*pi)),
                              coordinates="polar",
                              relTol=0.00001,
                              method="divonne"
  )$value
  eff_kern_df$selfrecruit[i] <- vol_selfrecruit
  
  for(dd in 1:nrow(dist_list)){
    # calculate an integral
    dist_integral <- (pi)/(4)*integral(function(x,y){(1/th_i)*exp(-sqrt(x^2+y^2)/th_i)},
                              bounds=list(x=c(dist_list$integral_eval[dd]-nav_rad,dist_list$integral_eval[dd]+nav_rad),
                                          y=c(-nav_rad,nav_rad)))$value
    patch_dists_integrals[dd,i] <- dist_integral
  }
}

######### Make and save connectivity matrices ##########
for(g in 1:nrow(group_index)){
  if(!file.exists(paste0(temp_path,"/map_",mapID_i,"/grp_",g))){
    # get the connectivity matrix among patches, given the group parameter values and patch-level q's
    # (and accounting for the patch population x per capita output b_i from each patch)
    
    # calculate this matrix
    v <- group_index[g,]
    # compute effective parameters for each patch with plasticity (once per group)
    eff_params <- f_plasticityb(patch_locations$b, 
                                v_p[v$p], 
                                v$alpha, 
                                v$theta,
                                n_alpha = length(v_alphas),
                                n_theta = length(v_thetas))
    alpha <- eff_params$alpha_plastic
    theta <- eff_params$theta_plastic
    npatch=length(alpha)
    # row of eff_kern_df that corresponds to each origin patch
    kernparam_rows <- sapply(1:npatch,function(i) eff_kern_index[eff_params$theta_plastic[i],eff_params$alpha_plastic[i]])
    
    connectivity_matrix <- vapply(1:npatch,function(patch_i){ # do once for each origin patch
      # Evaluate volume under kernel in a circle of radius nav_rad around destination site
      cm_i <- patch_dists_integrals[patch_dists_bin[,patch_i],kernparam_rows[patch_i]]
      cm_i[is.na(cm_i)] <- eff_kern_df$selfrecruit[kernparam_rows[patch_i]] # any site with distance 0
      # Do all three steps below in one
      cm_i <- (cm_i*patch_locations$overlap_discount)*eff_kern_df$offmap_corr[kernparam_rows[patch_i]]/(2*pi*v_thetas[theta[patch_i]])
      # # Divide by 2pitheta so all integrals sum to one
      # cm_i <- cm_i/(2*pi*theta[i])
      # # Apply a correction factor to account for the fact that broader kernels lose more offspring to edge effects
      # cm_i <- cm_i*eff_kern_df$offmap_corr[kernparam_rows[i]]
      # # Apply a correction for close neighbors
      # cm_i <- cm_i*patch_locations$overlap_discount
    },numeric(npatch))
    
    
    # ## This is what uses multiple cores: f_GetPlasticConnMat uses f_GetConnectivityMatrix_parallel
    # conn_mat <- f_GetConnMat(g,
    #                          group_index,patch_locations,v_p,v_alphas,v_thetas,eff_kern_index,eff_kern_df,patch_dists,nav_rad,
    #                          numCores=parallelly::availableCores())
    # # then store it
    write_fst(qDF(connectivity_matrix),paste0(temp_path,"/map_",mapID_i,"/grp_",g),compress=0)
    #saveRDS(object=conn_mat,file=paste0(saveto_folder,"/grp_",g,".rds"),compress=FALSE)
  }
} # g
