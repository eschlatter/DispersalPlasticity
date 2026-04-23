######## Generate the hab_params needed for sims ###########
source('0_Setup.R')
experiment_folder <- "experiments/Exp1_20260413"
basemap_folder <- "basemap_1x25"
seascapeset_folder <- "set1"
popmap_file <- "pop_density800"
nav_rad <- as_units(0.05,'km')
load(paste0(experiment_folder,"/",basemap_folder,"/",seascapeset_folder,"/df_choose.RData"))
v_width=10
v_cutoff=5000

slurm_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
h_i <- unique(df_dists_choose$h)[slurm_i]
df_i <- filter(df_dists_choose,abs(h-h_i)<0.001)
df_i_all <- df_i # hold all the hab files

for(sim_i in 1:nrow(df_i)){
  qmap_file <- paste0(seascapeset_folder,"/sim_",df_i$simID[sim_i])
  hab_file <- paste0("habfiles/hab_",df_i$simID[sim_i],"_identity")
  make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,qmap_file=qmap_file,popmap_file = popmap_file,
                                basemap_file=file.path(experiment_folder,basemap_folder),
                                overlap_method="simple",hab_file = hab_file)
  
  q_rast <- rast(paste0(experiment_folder,"/",basemap_folder,"/",qmap_file,".tif")) # load q_rast
  load(file=paste0(experiment_folder,"/",basemap_folder,"/pop_density800.RData")) # load popmap
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
    
    hab_file <- paste0("habfiles/hab_",df_i$simID[sim_i],"_",target_dist)
    
    make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,qmap_file=qmap_file,popmap_file = popmap_file,
                                  basemap_file=file.path(experiment_folder,basemap_folder),
                                  overlap_method="simple",hab_file = hab_file)
  } 
}

write.table(df_i_all,file=file.path(experiment_folder,basemap_folder,"habfiles",paste0("df_",slurm_i,".csv")),
            row.names=FALSE)
