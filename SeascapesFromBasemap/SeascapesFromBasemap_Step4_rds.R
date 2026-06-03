## script to run that just generates the connectivity matrices and saves them to global scratch
## use a bunch of cores, so this goes fast
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

# identify which map we're on (as indexed by experiment_index)
run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
map_info <- experiment_index[run_i,]
mapID_i <- experiment_index$mapID[run_i]

# params
load(paste0(experiment_folder,"/params_",map_info$param_id,".RData"))
list2env(x=params,envir=environment())

# hab_params
load(file=paste0(experiment_folder,"/habfiles/habparams_",mapID_i,".RData")) 
list2env(x=hab_params,envir=environment()) 

map_info <- list(basemap_id=1,popmap_id=2)
mapID_i=5

# patch_dists and patch_angles
load(file=paste0(experiment_folder,"/b",map_info$basemap_id,"/patchdists_b",map_info$basemap_id,"_p",map_info$popmap_id,".RData"))
patch_angles <- readRDS(paste0(temp_path,"/map_",mapID_i,"/patch_angles",".rds"))

if(hab_type=="points") patch_locations$K <- 1 # sometimes these end up as 0 from map resolution issues
units(patch_dists) <- NULL
units(nav_rad) <- NULL
units(patch_angles) <- NULL

## group_index: all unique combinations of parameters alpha, theta, and p
group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)

######### Make and save connectivity matrices ##########

for(g in 1:ngroups){
  # get the connectivity matrix among patches, given the group parameter values and patch-level q's
  # (and accounting for the patch population x per capita output b_i from each patch)
  
  # calculate this matrix
  ## This is what uses multiple cores: f_GetPlasticConnMat uses f_GetConnectivityMatrix_parallel
  conn_mat <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
                                  patch_dists=patch_dists, patch_angles=patch_angles,
                                  overlap_discount=overlap_discount,
                                  v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas,
                                  nav_rad=nav_rad, numCores=parallelly::availableCores())
  # then store it
  #write_fst(as.data.frame(conn_mat),paste0(temp_path,"/map_",mapID_i,"/grp_",g),compress=0)
  saveRDS(object=conn_mat,file=paste0(temp_path,"/grp_",g,".rds"),compress=FALSE)
  
} # g
