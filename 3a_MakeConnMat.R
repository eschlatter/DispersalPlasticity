## script to run that just generates the connectivity matrices and saves them to global scratch
## use a bunch of cores, so this goes fast
source('0_Setup.R')
experiment_folder <- "experiments/Exp1_20260413"
numCores <- parallelly::availableCores()

# identify which run we're on (as indexed by sim_design)
run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(run_i)) run_i <- 1

# Create a new folder in the MSI scratch directory
temp_dir <- paste0("sim_",run_i)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
dir.create(temp_path)

########## Load data ###########
# sim_design df (to show which params and hab_params to load)
load(paste0(experiment_folder,"/sim_design.RData")) 

# params
load(file=paste0(experiment_folder,"/params_",sim_design$params[run_i],".RData")) 
list2env(x=params,envir=environment())

# hab_params
load(file=paste0(experiment_folder,"/",sim_design$base[run_i],"/hab_",sim_design$hab[run_i],".RData")) 
list2env(x=hab_params,envir=environment()) 
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
  conn_mat <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
                                  patch_dists=patch_dists, patch_angles=patch_angles,
                                  overlap_discount=overlap_discount,
                                  v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas,
                                  nav_rad=nav_rad, numCores=numCores)
  # then store it
  write_fst(as.data.frame(conn_mat),paste0(temp_path,"/grp_",g),compress=0)
  #saveRDS(conn_mat,file=paste0(temp_path,"/grp_",g,".rds"),compress=FALSE)
  
} # g