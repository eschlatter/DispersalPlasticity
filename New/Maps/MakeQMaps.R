.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
library(sf)
library(terra)
library(raster)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggspatial)
source("New/functions.R")

experiment_folder <- "New/Maps"

# First: use SeascapesFromBasemap/SeascapesFromBasemap_Step2.R to
#  choose qmaps from among the very many I've already generated and 
#  save them in Maps/5x5 or Maps/1x25

# Next: for each map:
#  convert it to the correct scale (0-1 instead of 1-9)
#  give it a new number, save it in the full-pop basemap and add it to the index
#  edit it to match the relevant patchy basemap
#  give that a new number, save it in the patchy-pop basemap and add it to the index

#### 1x25 -- only do this once! ####
load("New/Maps/1x25/index.RData")
df_dists_choose$qID_full <- NA
df_dists_choose$qID_patchy <- NA
qmap_list <- list.files("New/Maps/1x25") 
qmap_list <- qmap_list[grep(".tif",qmap_list)]

variog_width=10
variog_cutoff=25000

for(qmap in qmap_list){
  basemapID <- 2 # full-occupancy map
  qmap_num <- as.numeric(gsub("[^0-9]", "", qmap))
  df_dists_i <- df_dists_choose[which(df_dists_choose$simID==qmap_num),]
  # transform to uniform distribution and convert to correct scale
  qrast <- rast(paste0("New/Maps/1x25/",qmap))
  values(qrast$q) <- f_TransformDist(values(qrast$q),"A") # transform to uniform distribution
  minq <- min(values(qrast$q))
  maxq <- max(values(qrast$q))
  values(qrast$q) <- (values(qrast$q)-minq)/(maxq-minq)
  
  # fit variogram for this map
  load(paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".RData")) # get reef_sf
  sfc_patches <- sf::st_sample(reef_sf,units::drop_units(st_area(reef_sf)/1000))
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(qrast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  vgm1 <- gstat::variogram(q~1,data=spdf1,cressie=TRUE,width=variog_width,cutoff=variog_cutoff) # empirical variogram
  vgmf <- gstat::fit.variogram(vgm1,gstat::vgm(c("Gau","Sph","Exp"))) # run several models, and pick the best one
  #plot(vgm1,vgmf)
  vgm_fit=list(range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
               model=vgmf$model[2])
  
  # give it a new number
  qID_full <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
  # save in the full-pop basemap
  writeRaster(qrast,filename=paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qID_full,".tif"),overwrite=TRUE)
  # add it to the index
  index_qmaps <- data.frame(qmap_id=qID_full,basemap_id=basemapID,description=paste0("autocorr=",round(vgm_fit$range),", 1x25km full"),
                            h=df_dists_i$h,target_dist="identity",
                            range=vgm_fit$range,sill=vgm_fit$sill,
                            SSErr=vgm_fit$SSErr,model=vgm_fit$model,original_range=df_dists_i$range)
  fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)
  
  # edit it to the relevant patchy basemap, and convert that back to the appropriate scale
  basemapID <- 6
  basemap_file <- paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID)
  qrast_patchy <- qrast
  base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
  # block out non-reef areas
  qmat <- values(qrast_patchy$q)
  qmat[is.na(values(base_rast))] <- NA # if they're NA's
  qmat[values(base_rast)==0] <- NA # if they're zeros
  values(qrast_patchy$q) <- qmat
  values(qrast_patchy$reef) <- qmat
  values(qrast_patchy$reef)[values(qrast_patchy$reef)>0] <- 1
  # transform to uniform distribution
  values(qrast_patchy$q) <- f_TransformDist(values(qrast_patchy$q),"A") 
  # convert to 0-1 scale
  minq <- min(values(qrast_patchy$q),na.rm=TRUE)
  maxq <- max(values(qrast_patchy$q),na.rm=TRUE)
  values(qrast_patchy$q) <- (values(qrast_patchy$q)-minq)/(maxq-minq)
  #  give that a new number, save it in the patchy-pop basemap and add it to the index
  qID_patchy <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
  
  
  # fit variogram for this map
  load(paste0(basemap_file,".RData")) # get reef_sf
  sfc_patches <- sf::st_sample(reef_sf,round(units::drop_units(st_area(reef_sf)/1000)))
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(qrast_patchy$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  vgm1 <- gstat::variogram(q~1,data=spdf1,width=variog_width,cutoff=variog_cutoff) # empirical variogram
  vgmf <- gstat::fit.variogram(vgm1,gstat::vgm(c("Gau","Sph","Exp"))) # run several models, and pick the best one
  #plot(vgm1,vgmf)
  vgm_fit=list(range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
               model=vgmf$model[2])
  
  
  # save in the patchy-pop basemap
  writeRaster(qrast_patchy,filename=paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qID_patchy,".tif"),overwrite=TRUE)
  # add it to the index
  index_qmaps <- data.frame(qmap_id=qID_patchy,basemap_id=basemapID,description=paste0("autocorr=",round(vgm_fit$range),", 1x25km patchy"),
                            h=df_dists_i$h,target_dist="identity",
                            range=vgm_fit$range,sill=vgm_fit$sill,
                            SSErr=vgm_fit$SSErr,model=vgm_fit$model,original_range=df_dists_i$range)
  fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)
  # add the new qID numbers to df_dists
  df_dists_i$qID_full <- qID_full
  df_dists_i$qID_patchy <- qID_patchy
  df_dists_choose[which(df_dists_choose$simID==qmap_num),] <- df_dists_i
}

save(df_dists_choose,file="New/Maps/1x25/index.RData")


#### 5x5 -- only do this once! ####
load("New/Maps/5x5/index.RData")
df_dists_choose$qID_full <- NA
df_dists_choose$qID_patchy <- NA
qmap_list <- list.files("New/Maps/5x5") 
qmap_list <- qmap_list[grep(".tif",qmap_list)]

variog_width=10
variog_cutoff=5000

for(qmap in qmap_list){
  basemapID <- 1 # full-occupancy map
  qmap_num <- as.numeric(gsub("[^0-9]", "", qmap))
  df_dists_i <- df_dists_choose[which(df_dists_choose$simID==qmap_num),]
  # transform to uniform distribution and convert to correct scale
  qrast <- rast(paste0("New/Maps/5x5/",qmap))
  values(qrast$q) <- f_TransformDist(values(qrast$q),"A") # transform to uniform distribution
  minq <- min(values(qrast$q))
  maxq <- max(values(qrast$q))
  values(qrast$q) <- (values(qrast$q)-minq)/(maxq-minq)
  
  # fit variogram for this map
  load(paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".RData")) # get reef_sf
  sfc_patches <- sf::st_sample(reef_sf,units::drop_units(st_area(reef_sf)/1000))
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(qrast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  vgm1 <- gstat::variogram(q~1,data=spdf1,cressie=TRUE,width=variog_width,cutoff=variog_cutoff) # empirical variogram
  vgmf <- gstat::fit.variogram(vgm1,gstat::vgm(c("Gau","Sph","Exp"))) # run several models, and pick the best one
  #plot(vgm1,vgmf)
  vgm_fit=list(range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
               model=vgmf$model[2])
  
  # give it a new number
  qID_full <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
  # save in the full-pop basemap
  writeRaster(qrast,filename=paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qID_full,".tif"),overwrite=TRUE)
  # add it to the index
  index_qmaps <- data.frame(qmap_id=qID_full,basemap_id=basemapID,description=paste0("autocorr=",round(vgm_fit$range),", 5x5km full"),
                            h=df_dists_i$h,target_dist="identity",
                            range=vgm_fit$range,sill=vgm_fit$sill,
                            SSErr=vgm_fit$SSErr,model=vgm_fit$model,original_range=df_dists_i$range)
  fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)
  
  # edit it to the relevant patchy basemap, and convert that back to the appropriate scale
  basemapID <- 5
  basemap_file <- paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID)
  qrast_patchy <- qrast
  base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
  # block out non-reef areas
  qmat <- values(qrast_patchy$q)
  qmat[is.na(values(base_rast))] <- NA # if they're NA's
  qmat[values(base_rast)==0] <- NA # if they're zeros
  values(qrast_patchy$q) <- qmat
  values(qrast_patchy$reef) <- qmat
  values(qrast_patchy$reef)[values(qrast_patchy$reef)>0] <- 1
  # transform to uniform distribution
  values(qrast_patchy$q) <- f_TransformDist(values(qrast_patchy$q),"A") 
  # convert to 0-1 scale
  minq <- min(values(qrast_patchy$q),na.rm=TRUE)
  maxq <- max(values(qrast_patchy$q),na.rm=TRUE)
  values(qrast_patchy$q) <- (values(qrast_patchy$q)-minq)/(maxq-minq)
  #  give that a new number, save it in the patchy-pop basemap and add it to the index
  qID_patchy <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
  
  
  # fit variogram for this map
  load(paste0(basemap_file,".RData")) # get reef_sf
  sfc_patches <- sf::st_sample(reef_sf,round(units::drop_units(st_area(reef_sf)/500)))
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(qrast_patchy$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  vgm1 <- gstat::variogram(q~1,data=spdf1,width=variog_width,cutoff=variog_cutoff) # empirical variogram
  vgmf <- gstat::fit.variogram(vgm1,gstat::vgm(c("Gau","Sph","Exp"))) # run several models, and pick the best one
  #plot(vgm1,vgmf)
  vgm_fit=list(range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
               model=vgmf$model[2])
  
  
  # save in the patchy-pop basemap
  writeRaster(qrast_patchy,filename=paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qID_patchy,".tif"),overwrite=TRUE)
  # add it to the index
  index_qmaps <- data.frame(qmap_id=qID_patchy,basemap_id=basemapID,description=paste0("autocorr=",round(vgm_fit$range),", 5x5km patchy"),
                            h=df_dists_i$h,target_dist="identity",
                            range=vgm_fit$range,sill=vgm_fit$sill,
                            SSErr=vgm_fit$SSErr,model=vgm_fit$model,original_range=df_dists_i$range)
  fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)
  # add the new qID numbers to df_dists
  df_dists_i$qID_full <- qID_full
  df_dists_i$qID_patchy <- qID_patchy
  df_dists_choose[which(df_dists_choose$simID==qmap_num),] <- df_dists_i
}

save(df_dists_choose,file="New/Maps/5x5/index.RData")





