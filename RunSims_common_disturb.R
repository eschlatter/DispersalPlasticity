## After connectivity matrices have been made, run the sims
source('0_Setup.R')
library(profvis)
experiment_name <- "Exp3_20260528"
experiment_subname <- "common_disturb_b01"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index_disturbset.csv"))

# experiment_index <- rbind(experiment_index,experiment_index,experiment_index)
# experiment_index$repID <- rep(1:3,each=4)
experiment_index$repID <- 1

## only did this once
# experiment_index <- experiment_index[c(17,1,5,13),]
# experiment_index$param_id <- 3
# write.csv(experiment_index,file=paste0(experiment_folder,"/experiment_index_disturbset.csv"),row.names=FALSE)

load(paste0(experiment_folder,"/basemap_index.RData")) # basemap index

# identify which run we're on (as indexed by sim_design)
run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run_info <- experiment_index[run_i,]

# Identify folders in the project directory and MSI scratch directory
temp_dir <- paste0("connmats_smallerkern2_",experiment_name)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
#temp_path <- paste0("/dev/shm")
keep_path <- paste0(experiment_folder,"/habfiles")

########## Load data ###########
# params
load(file=paste0(experiment_folder,"/params_",run_info$param_id,".RData"))
params$disturb_prob <- 1
list2env(x=params,envir=environment())

# hab_params
load(file=paste0(keep_path,"/habparams_",run_info$mapID,".RData")) 
hab_params$patch_locations$b <- hab_params$patch_locations$b/10

# patch_dists
patchdist_file <- paste0(experiment_folder,"/b",run_info$basemap_id,"/patchdists_b",run_info$basemap_id,
                         "_p",run_info$popmap_id,".RData")

jobmem_GB <- as.numeric(Sys.getenv("SLURM_MEM_PER_NODE"))/1000

#### Run sim
f_RunSimComboStorageLite(params,hab_params=hab_params,output_flag="lite",show_plot = FALSE,output_thin=20,
                         output_file=paste0(experiment_folder,"/output/",experiment_subname,"/rep",run_info$repID,"/map_",run_info$mapID),run_i=run_i,
                         connmat_folder=paste0(temp_path,"/map_",run_info$mapID),connmat_format="fst",
                         patchdist_file=patchdist_file,
                         connmat_size_GB = 3.2, jobmem_GB = jobmem_GB,normalize_connmats=TRUE,
                         disturb_freq=1,disturb_radius=sqrt(0.03*25/pi))

# prof_sim <- profvis({f_RunSimComboStorageLite(params,hab_params=hab_params,output_flag="lite",show_plot = FALSE,output_thin=25,
#                          output_file=paste0(experiment_folder,"/output/reptest/map_",run_info$mapID),run_i=run_i,
#                          connmat_folder=paste0(temp_path,"/map_",run_info$mapID),connmat_format="fst",
#                          patchdist_file=patchdist_file,
#                          connmat_size_GB = 3.2, jobmem_GB = 300)
# })
# htmlwidgets::saveWidget(prof_sim,"profile.html")
#browseURL("profile.html")

