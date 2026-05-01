####### For any row of the sim index (i.e., basemap/popmap/qmap/target_dist/params combo),
####### generate the remaining files it needs to run a sim, and store them appropriately

f_MakeSeascapeFiles <- function(basemap_id,popmap_id,qmap_id,target_dist,param_id,target_dist,
                                experiment_folder,seascapeset_folder,temp_path,sub_fold,keep_path,
                                ){

# params
load(file=paste0(experiment_folder,"/params_",param_id,".RData")) 
list2env(x=params,envir=environment())

######### Load patch_dists, make patch_angles, save patch_angles ##########
# load patch_dists
load(file=paste0(experiment_folder,"/b",basemap_id,"/patchdists_b",basemap_id,"_p",popmap_id,".RData")) 
# make patch_angles
patch_angles <- (2*nav_rad)/(2*pi*patch_dists)
patch_angles[is.nan(patch_angles)] <- 1
patch_angles[drop_units(patch_angles)>1] <- 1
# save patch_angles
write_fst(as.data.frame(patch_angles),paste0(temp_path,"/",sub_fold,"/patch_angles"),compress=0)

######### Generate hab_params (without patch_angles or patch_dists) and save locally ##########
#### (basically, this is f_MakeHabitat)
load(paste0(experiment_folder,"/b",basemap_id,"/pop_b",basemap_id,"_p",popmap_id,".RData")) # loads reef_sf,patch_dists,sfc_patches,hab_type
q_rast <- rast(paste0(seascapeset_folder,"/qmap_b",basemap_id,"_q",qmap_id,".tif")) # load q_rast
K_rast <- rast(paste0(experiment_folder,"/b",basemap_id,"/pop_b",basemap_id,"_p",popmap_id,".tif")) # load K_rast

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

# make overlap_discount
n_neighbors <- rowSums(patch_dists<nav_rad) # number of points within distance nav_rad of focal point (including focal point)
overlap_discount <- 1/n_neighbors

hab_file <- paste0("/habfile_b",basemap_id,"_p",popmap_id,"_q",qmap_id,"_d",target_dist)
hab_params <- list(npatch=npatch,
                   patch_locations=df_patches,
                   overlap_discount=overlap_discount,
                   reef_sf=reef_sf,
                   sfc_patches=sfc_patches,
                   hab_type=hab_type,
                   nav_rad=nav_rad,
                   hab_file=hab_file)

save(hab_params,file=paste0(keep_path,hab_file,".RData"))
writeRaster(hab_rast,filename=paste0(keep_path,hab_file,".tif"),overwrite=TRUE)






######### Make and save connectivity matrices ##########
## group_index: all unique combinations of parameters alpha, theta, and p
group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)

for(g in 1:ngroups){
  # get the connectivity matrix among patches, given the group parameter values and patch-level q's
  # (and accounting for the patch population x per capita output b_i from each patch)
  
  # calculate this matrix
  conn_mat <- f_GetPlasticConnMat(g=g, group_index=group_index, patch_locations=patch_locations,
                                  patch_dists=patch_dists, patch_angles=patch_angles,
                                  overlap_discount=overlap_discount,
                                  v_p=v_p, v_alphas=v_alphas, v_thetas=v_thetas,
                                  nav_rad=nav_rad, numCores=numCores)
  # then store it
  write_fst(as.data.frame(conn_mat),paste0(temp_path,"/grp_",g),compress=0)
  #saveRDS(conn_mat,file=paste0(temp_path,"/grp_",g,".rds"),compress=FALSE)
  
} # g




# # hab_params
# load(file=paste0(experiment_folder,"/",sim_design$base[run_i],"/hab_",sim_design$hab[run_i],".RData")) 
# list2env(x=hab_params,envir=environment()) 
# if(hab_type=="points") patch_locations$K <- 1 # sometimes these end up as 0 from map resolution issues
# units(patch_dists) <- NULL
# units(nav_rad) <- NULL
# units(patch_angles) <- NULL


}