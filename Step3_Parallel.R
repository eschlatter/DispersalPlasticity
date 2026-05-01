######## Generate the hab_params needed for sims ###########
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

## One slurm task, lots of cores
numCores <- parallelly::availableCores()

## First: Parallelize over q_ID; make the target_dist maps and calculate the variogram stats
i <- 1
qID_i <- unique(experiment_index$q_ID)[i] # get current value of h


## Second: Parallelize over mapID; make and store the patch_angles, conn_mats, hab_params



index_i <- filter(experiment_index,q_ID==qID_i)


for(i_map in 1:nrow(index_i)){
  
}

## Get params for the current maps
basemap_i <- first(index_i$basemap_id)
popmap_i <- first(index_i$popmap_id)
param_i <- first(index_i$param_id)

# variogram-fitting and navigation radius parameters
nav_rad <- as_units(0.05,'km')
variog_width <- 10
variog_cutoff <- 5000

# load index of maps chosen from previous step
load(paste0(seascapeset_folder,"/df_step2.RData"))

# pick the maps for the current value of h
h_i <- unique(df_dists_choose$h)[slurm_i] # get current value of h
df_step2_i <- filter(df_dists_choose,abs(h-h_i)<0.001) # make an index of each map we're starting with
df_step3_i <- df_step2_i # create new index to hold the existing map files, plus their target-distribution variants

for(sim_i in 1:nrow(df_step2_i)){
  q_id <- df_step2_i$simID[sim_i] # ID for this qmap
  
  # first, make and save the patch_angles and conn_mats for that map
  # make directory to put them in
  sub_fold <- paste0("b",basemap_id,"_p",popmap_id,"_q",q_id,"_dI")
  dir.create(file.path(temp_path,sub_fold))
  
  # generate conn_mats (save in temp folder)
  # generate patch_angles (save in temp folder)
  # generate hab_params (save locally)
  
  
  qmap_file <- paste0(seascapeset_folder,"/sim_",df_i$simID[sim_i])
  hab_file <- paste0("habfiles/hab_",df_i$simID[sim_i],"_identity")
  make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,qmap_file=qmap_file,popmap_file = popmap_file,
                                basemap_file=file.path(experiment_folder,basemap_folder),
                                overlap_method="simple",hab_file = hab_file)
  
  q_rast <- rast(paste0(experiment_folder,"/",basemap_folder,"/",qmap_file,".tif")) # load q_rast
  load(file=paste0(experiment_folder,"/",basemap_folder,"/pop_density800.RData")) # load popmap
  
  # then, for target distributions, make the map, save the variogram stats, and make and save the conn_mats
  for(target_dist in c("C","D","E")){
    qmap_file <- paste0(seascapeset_folder,"/sim_",df_i$simID[sim_i],"_",target_dist)
    
    q_rast_new <- q_rast
    q_rast_new$q <- f_TransformDist(matrix(values(q_rast$q),nrow=nrow(q_rast),byrow=TRUE),target_dist)
    writeRaster(q_rast_new,filename=paste0(experiment_folder,"/",basemap_folder,"/",qmap_file,".tif"),overwrite=TRUE)
    
    # get variogram info
    spdf1 <- as_Spatial(sfc_patches)
    spdf1$q <- terra::extract(q_rast_new$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
    # empirical variogram
    vgm1 <- variogram(q~1,data=spdf1,cressie=TRUE,width=v_width,cutoff=v_cutoff)
    # run gaussian and spherical models, and pick the better one
    vgmf <- fit.variogram(vgm1,vgm(c("Gau","Sph","Exp")))
    
    # store output
    df_i_all <- rbind(df_i_all,data.frame(simID=df_i$simID[sim_i],
                                          range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
                                          model=vgmf$model[2],h=h_i,rep_i=NA,target_dist=target_dist))
    
    
    
    # make directory for the conn_mats and patch_angles
    sub_fold <- paste0("b",basemap_id,"_p",popmap_id,"_q",q_id,"_d",target_dist)
    dir.create(temp_path,sub_fold)
    
    hab_file <- paste0("habfiles/hab_",df_i$simID[sim_i],"_",target_dist)
    
    make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,qmap_file=qmap_file,popmap_file = popmap_file,
                                  basemap_file=file.path(experiment_folder,basemap_folder),
                                  overlap_method="simple",hab_file = hab_file)
  } 
}

write.table(df_i_all,file=file.path(experiment_folder,basemap_folder,"habfiles",paste0("df_",slurm_i,".csv")),
            row.names=FALSE)
