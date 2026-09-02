.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
source("New/SimFn.R")
library(data.table)
library(parallel)

repID <- ifelse(is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))),1,as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
output_flag="lite"
output_thin=1
experiment_folder <- "New"  #name of directory where everything is stored
notes <- "dont normalize, q=0.5"   #can include a character string here with notes on the sim

# bio parameters
mu <- 0.001
nav_rad <- 0.05
adult_survival_prob <- 0
base_fecund <- 20
habID <- 16

# sim parameters
theta_start_min <- 0.005
theta_start_max <- 0.005
p_start_min <- 0
p_start_max <- 0
nsteps <- 5000

# sim options
normalize_offspring=FALSE
plasticity_on=FALSE
larval_output_by_theta=FALSE

NewSimFn2(repID=repID,
         output_flag=output_flag,
         output_thin=output_thin,
         experiment_folder=experiment_folder,
         notes=notes,
         
         # bio parameters
         mu=mu,
         nav_rad=nav_rad,
         adult_survival_prob=adult_survival_prob,
         base_fecund=base_fecund,
         habID=habID,
         
         # sim parameters
         theta_start_min=theta_start_min,
         theta_start_max=theta_start_max,
         p_start_min=p_start_min,
         p_start_max=p_start_max,
         nsteps=nsteps,
         
         # sim options
         normalize_offspring=normalize_offspring,
         plasticity_on=plasticity_on,
         larval_output_by_theta=larval_output_by_theta)
