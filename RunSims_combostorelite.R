## After connectivity matrices have been made, run the sims
source('0_Setup.R')
library(profvis)
experiment_name <- "Exp3_20260528"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
load(paste0(experiment_folder,"/basemap_index.RData")) # basemap index

# identify which run we're on (as indexed by sim_design)
run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run_info <- experiment_index[run_i,]

# Identify folders in the project directory and MSI scratch directory
temp_dir <- paste0("connmats_",experiment_name)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
#temp_path <- paste0("/dev/shm")
keep_path <- paste0(experiment_folder,"/habfiles2")

source("functions/f_RunSimComboStorageLite.R")

########## Load data ###########
# params
load(file=paste0(experiment_folder,"/params_",run_info$param_id,".RData")) 
list2env(x=params,envir=environment())

# hab_params
load(file=paste0(keep_path,"/habparams_",run_info$mapID,".RData")) 

# patch_dists
patchdist_file <- paste0(experiment_folder,"/b",run_info$basemap_id,"/patchdists_b",run_info$basemap_id,
                         "_p",run_info$popmap_id,".RData")

jobmem_GB <- as.numeric(Sys.getenv("SLURM_MEM_PER_NODE"))/1000

#### Run sim
f_RunSimComboStorageLite(params,hab_params=hab_params,output_flag="lite",show_plot = FALSE,output_thin=20,
                                              output_file=paste0(experiment_folder,"/output/rep2/map_",run_info$mapID),run_i=run_i,
                                              connmat_folder=paste0(temp_path,"/map_",run_info$mapID),connmat_format="fst",
                                              patchdist_file=patchdist_file,
                                              connmat_size_GB = 3.2, jobmem_GB = jobmem_GB)


# prof_sim <- profvis({f_RunSimComboStorageLite(params,hab_params=hab_params,output_flag="lite",show_plot = FALSE,output_thin=25,
#                          output_file=paste0(experiment_folder,"/output/reptest/map_",run_info$mapID),run_i=run_i,
#                          connmat_folder=paste0(temp_path,"/map_",run_info$mapID),connmat_format="fst",
#                          patchdist_file=patchdist_file,
#                          connmat_size_GB = 3.2, jobmem_GB = 300)
# })
# htmlwidgets::saveWidget(prof_sim,"profile.html")
#browseURL("profile.html")

# f_RunSimNewScratch(params,hab_params=hab_params,output_flag="all",show_plot = FALSE,output_thin=25,
#                    output_file=paste0(experiment_folder,"/output/map_",run_info$mapID),run_i=run_i,
#                    connmat_folder=paste0(temp_path,"/map_",run_info$mapID),connmat_format="fst")

# #### Post-processing
# popmat <- fread(paste0(experiment_folder,"/output/map_",run_info$mapID,"_raw.csv")) # population output file
# load(paste0(experiment_folder,"/output/map_",run_info$mapID,"_metadata.RData")) # metadata file
# 
# sim_melt <- f_ProcessPopmat(popmat,metadata_list$patch_locations,metadata_list$group_index,
#                             params=metadata_list$params)
# save(sim_melt,file=paste0(experiment_folder,"/output/map_",run_info$mapID,"_processed.RData"))

# #### Make some quick check plots
# source('0_Setup.R')
# load(paste0(experiment_folder,"/output/map_",run_info$mapID,"_processed.RData")) # processed data
# load(paste0(experiment_folder,"/output/map_",run_info$mapID,"_metadata.RData")) # metadata file
# load(paste0(keep_path,"/habparams_",run_info$mapID,".RData")) # habitat file
# by_t <- sim_melt$by_t
# sim_melt <- sim_melt$sim_melt
# 
# # a <- f_PlotAllHeatmaps(sim_melt = sim_melt,patch_locations = hab_params$patch_locations,reef_sf = hab_params$reef_sf,hab_type = hab_params$hab_type)
# ggplot(by_t,aes(x=t_i,y=p))+geom_line(aes(color="p"))+geom_line(aes(y=alpha*theta,color="kernmean"))
# ggplot(by_t,aes(x=t_i))+geom_line(aes(y=theta,color="theta"))+geom_line(aes(y=alpha,color="alpha"))+geom_line(aes(y=p,color="p"))
# ggplot(by_t,aes(x=t_i,y=popsize))+geom_line()
