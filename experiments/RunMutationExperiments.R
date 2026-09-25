if(version$major=="4" & version$minor=="4.0") .libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old") else .libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
library(data.table)
library(parallel)
library(terra)
library(sf)
library(calculus)
library(dplyr)
source("New/SimFn.R")
source("New/functions.R")

# Running one experiment per basemap: 1 (5x5), 2 (1x25), 7 (5x5 patchy), 6 (1x25 patchy), 4 (Kimbe)
basemapID <- 1

repID <- ifelse(is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))),1,as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
output_flag="lite" # "lite" or "all"; isn't actually doing anything at this point -- we're outputting summary info at interval output_thin, and full pop info at an interval hard-coded in (I know, I know, but I don't want to add another parameter right now)
output_thin=1

# Experiment ID
experiment_folder <- "New"  #name of directory where everything is stored
notes <- paste0("Basemap ",basemapID," experiment, mutation magnitude and rate")   #can include a character string here with notes on the sim
load(paste0("New/df_experiments_b1_mutation.RData"))
exp_i <- df_experiments[repID,]
rm(df_experiments,df_qmaps)
experimentID <- exp_i$experiment_id

# bio parameters
mutation_type <- exp_i$mutation_type
mu <- exp_i$mu_p # for plasticity
mu_theta <- exp_i$mu_theta # in km
nav_rad <- 0.05
adult_survival_prob <- 0
base_fecund <- 10
habID <- exp_i$hab_id

# sim parameters
theta_start_min <- 0.005
theta_start_max <- 140
p_start_min <- -10
p_start_max <- 10
nsteps <- 3000

# sim options
normalize_offspring=FALSE
plasticity_on="multiplicative"
larval_output_by_theta=FALSE # save generic data about the relationship between theta and larval output

# disturbance parameters
Dp=exp_i$Dp # probability of disturbance (per timestep)
De=exp_i$De # extent of disturbance
Dl=exp_i$Dl # duration of disturbance
disturb_method="fractal" # "fractal" or "circle"
Dm=exp_i$Dm # magnitude of disturbance (disturbed anemones' fecundity multiplied by this amt; smaller is more impactful)

NewSimFn3(repID=repID,
          experimentID=experimentID,
          output_flag=output_flag,
          output_thin=output_thin,
          experiment_folder=experiment_folder,
          notes=notes,
          
          # bio parameters
          mutation_type=mutation_type,
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