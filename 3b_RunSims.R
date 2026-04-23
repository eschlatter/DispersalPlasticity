## After connectivity matrices have been made, run the sims
source('0_Setup.R')
experiment_folder <- "experiments/Exp1_20260413"

# identify which run we're on (as indexed by sim_design)
run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(run_i)) run_i <- 1

# Create a new folder in the MSI scratch directory
temp_dir <- paste0("sim_",run_i)
temp_path <- file.path("/scratch.global","schla103",temp_dir)

########## Load data ###########
# sim_design df (to show which params and hab_params to load)
load(paste0(experiment_folder,"/sim_design.RData")) 

# params
load(file=paste0(experiment_folder,"/params_",sim_design$params[run_i],".RData")) 
list2env(x=params,envir=environment())

# hab_params
load(file=paste0(experiment_folder,"/",sim_design$base[run_i],"/hab_",sim_design$hab[run_i],".RData")) 

#### Run sim
f_RunSimNew(params,hab_params=hab_params,output_flag="all",show_plot = FALSE,output_thin=25,
            output_file=paste0(experiment_folder,"/output/sim_",sim_file))

#### Post-processing
popmat <- fread(paste0(experiment_folder,"/output/sim_",sim_file,"_raw.csv")) # population output file
load(paste0(experiment_folder,"/output/sim_",sim_file,"_metadata.RData")) # metadata file

sim_melt <- f_ProcessPopmat(popmat,metadata_list$patch_locations,metadata_list$group_index,
                            params=metadata_list$params)
save(sim_melt,file=paste0(experiment_folder,"/output/sim_",sim_file,"_processed.RData"))