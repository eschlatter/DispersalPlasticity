f_MakeHabitat <- function(nav_rad,hab_type="points",
                          qmap_file=NULL,popmap_file=NULL,hab_file=NULL,patchdist_file=NULL,
                          qmap_id=NULL,popmap_id=NULL,habmap_id=NULL,experiment_folder=NULL){
  # load in data
  if(!is.null(qmap_id)) qmap_file <- paste0(experiment_folder,"/b",basemap_id,"/set",set_id,"/qmap_b",basemap_id,"_q",qmap_id)
  if(!is.null(popmap_id)) popmap_file <- paste0(experiment_folder,"/b",basemap_id,"/pop_b",basemap_id,"_p",popmap_id)
  if(!is.null(hab_id)) hab_file <- paste0(experiment_folder,"/habfiles/habparams_",hab_id)
  if(is.null(patchdist_file)) patchdist_file <- paste0(experiment_folder,"/b",basemap_id,"/patchdists_b",basemap_id,"_p",popmap_id,".RData")
  
  load(paste0(popmap_file,".RData")) # loads reef_sf,patch_dists,sfc_patches,hab_type
  q_rast <- rast(paste0(qmap_file,".tif")) # load q_rast
  K_rast <- rast(paste0(popmap_file,".tif")) # load K_rast
  load(patchdist_file) # load patch_dists
  
  units(nav_rad) <- 'km'
  npatch <- length(sfc_patches)
  
  ## put q_rast and K_rast together
  hab_rast <- c(q_rast,K_rast)
  
  ## create df_patches (important: ID should be in the same order as in sfc_patches, or distance matrix will be wrong)
  q_vect <- terra::extract(hab_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)
  K_vect <- terra::extract(hab_rast$K,vect(sfc_patches),xy=TRUE,search_radius=500)
  patch_coords <- st_coordinates(sfc_patches)
  df_patches <- cbind(q_vect[,c("ID","q")],patch_coords)
  df_patches$K <- K_vect$K[df_patches$ID]
  df_patches$id <- df_patches$ID
  df_patches$x <- df_patches$X
  df_patches$y <- df_patches$Y
  df_patches <- df_patches[,c("id","x","y","q","K")]
  df_patches$b <- f_q_to_b(df_patches$q) # calculate reproductive rate (b) from habitat quality (q)
  
  ## make overlap_discount
  overlap_discount <- 1/rowSums(patch_dists<nav_rad)
  
  hab_params <- list(npatch=npatch,
                     patch_locations=df_patches,
                     overlap_discount=overlap_discount,
                     reef_sf=reef_sf,
                     sfc_patches=sfc_patches,
                     hab_type=hab_type,
                     nav_rad=nav_rad,
                     mapID=hab_id)
  
  save(hab_params,file=paste0(hab_file,".RData"))
  writeRaster(hab_rast,filename=paste0(hab_file,".tif"),overwrite=TRUE)
}