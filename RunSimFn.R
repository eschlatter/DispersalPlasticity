if(version$major=="4" & version$minor=="4.0") .libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old") else .libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
library(data.table)
library(parallel)
library(terra)
library(sf)
library(calculus)
source("functions/Function_simulate.R")
source("functions/Functions_auxiliary.R")

repID <- ifelse(is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))),1,as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
v_De <- c(0.9,0.9,0.8,0.8,0.7,0.7,0.6,0.6)
print(paste0("Experiment3: limits of disturbance while maintaining population persistence,  De=",v_De[repID],", theta=1km init condit"))

output_flag="lite" # "lite" or "all"
output_thin=1
experiment_folder <- "New"  #name of directory where everything is stored
notes <- paste0("Finding top De value (try ",v_De[repID],")")   #can include a character string here with notes on the sim

# bio parameters
mu <- 0.001 # for plasticity
mu_theta <- 0.01 # in km
nav_rad <- 0.05
adult_survival_prob <- 0
base_fecund <- 10
habID <- 8

# sim parameters
# theta_start_min <- 0.005
# theta_start_max <- 140
theta_start_min <- 1
theta_start_max <- 1
p_start_min <- -10
p_start_max <- 10
nsteps <- 400

# sim options
normalize_offspring=FALSE
plasticity_on="multiplicative"
larval_output_by_theta=FALSE # save generic data about the relationship between theta and larval output

# disturbance parameters
Dp=0.1 # probability of disturbance (per timestep)
De=v_De[repID] # extent of disturbance
Dl=16 # duration of disturbance
disturb_method="fractal" # "fractal" or "circle"
Dm=0 # magnitude of disturbance (disturbed anemones' fecundity multiplied by this amt; smaller is more impactful)

SimFn(repID=repID,
      output_flag=output_flag,
      output_thin=output_thin,
      experiment_folder=experiment_folder,
      notes=notes,
      
      # bio parameters
      mu=mu,
      mu_theta=mu_theta,
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
      larval_output_by_theta=larval_output_by_theta,
      
      # disturbance parameters
      Dp=Dp,
      De=De,
      Dl=Dl,
      disturb_method=disturb_method,
      Dm=Dm)