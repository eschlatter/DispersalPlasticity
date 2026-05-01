## After connectivity matrices have been made, run the sims
source('0_Setup.R')
experiment_name <- "Exp2_20260423"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
#load(paste0(experiment_folder,"/basemap_index.RData")) # basemap index

# Identify folders in the project directory and MSI scratch directory
temp_dir <- paste0("connmats_",experiment_name)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
#temp_path <- paste0("/dev/shm")
keep_path <- paste0(experiment_folder,"/habfiles")


#### Post-processing: for loop

# identify which run we're on (as indexed by sim_design)
#run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run_i=2
run_info <- experiment_index[run_i,]

popmat <- fread(paste0(experiment_folder,"/output/map_",run_info$mapID,"_raw.csv")) # population output file
load(paste0(experiment_folder,"/output/map_",run_info$mapID,"_metadata.RData")) # metadata file

sim_melt <- f_ProcessPopmat(popmat,metadata_list$patch_locations,metadata_list$group_index,
                            params=metadata_list$params)
by_t <- sim_melt$by_t
sim_melt <- sim_melt$sim_melt


#save(sim_melt,file=paste0(experiment_folder,"/output/map_",run_info$mapID,"_processed.RData"))

#### Make some quick check plots
# source('0_Setup.R')
# load(paste0(experiment_folder,"/output/map_",run_info$mapID,"_processed.RData")) # processed data
# load(paste0(experiment_folder,"/output/map_",run_info$mapID,"_metadata.RData")) # metadata file
# load(paste0(keep_path,"/habparams_",run_info$mapID,".RData")) # habitat file


# a <- f_PlotAllHeatmaps(sim_melt = sim_melt,patch_locations = hab_params$patch_locations,reef_sf = hab_params$reef_sf,hab_type = hab_params$hab_type)
ggplot(by_t,aes(x=t_i,y=p))+geom_line(aes(color="p"))+geom_line(aes(y=alpha*theta,color="kernmean"))
#ggplot(by_t,aes(x=t_i))+geom_line(aes(y=theta,color="theta"))+geom_line(aes(y=alpha,color="alpha"))+geom_line(aes(y=p,color="p"))
ggplot(by_t,aes(x=t_i,y=popsize))+geom_line()
