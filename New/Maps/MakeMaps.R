.libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
library(sf)
#library(raster)
library(data.table)
library(terra)
library(ggplot2)
library(ggspatial)
library(dplyr)
source("New/functions.R")

experiment_folder <- "New/Maps"

# # Make indices for basemaps, popmaps, and qmaps. Each time a new one is generated, put it in the index.
# index_basemaps <- data.frame(basemapID=integer(),description=character(),area_km2=double())
# fwrite(index_basemaps,file="New/Maps/index_basemaps.csv")
# 
# index_popmaps <- data.frame(popmapID=integer(),basemapID=integer(),description=character(),density=double(),npatch=integer())
# fwrite(index_popmaps,file="New/Maps/index_popmaps.csv")
# 
# index_qmaps <- data.frame(qmapID=integer(),basemapID=integer(),description=character(),
#                           h=double(),target_dist=character(),
#                           range=double(),sill=double(),SSErr=double(),model=character())
# fwrite(index_qmaps,file="New/Maps/index_qmaps.csv")
# 
# index_habs <- data.frame(hab_id=integer(),basemap_id=integer(),qmap_id=integer(),
#                          popmap_id=integer(),npatch=integer(),q_autocorr_scale=double())
# fwrite(index_habs,file=paste0(experiment_folder,"/index_habs.csv"),append = TRUE)

#### b1: 5x5km, uniform ####
basemapID=1
hab_sim <- f_GenerateBasemap(x_dist=5000,y_dist=5000,resol=c(10,10),method="uniform",h=NA,prop_hab=NA,
                             make_dist_mat = FALSE,plot_flag=TRUE,experiment_folder=experiment_folder,basemap_id = basemapID)
index_basemaps <- data.frame(basemapID=basemapID,description="5x5km uniform",area_km2=25)
fwrite(index_basemaps,file="New/Maps/index_basemaps.csv",append=TRUE)

##### b1 popmaps #####
# 1: small-density (8 anems/km2)
popmapID=1
pts_out <- f_SimPtsOnMap(n_anems=25*8,inwater_dist=FALSE,samp_type = "random",plot_flag=FALSE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmapID=popmapID,basemapID=basemapID,
                            description="low density",density=8,npatch=200)
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

# 2: full-density (800 anems/km2)
popmapID=2
pts_out <- f_SimPtsOnMap(n_anems=25*800,inwater_dist=FALSE,samp_type = "random",plot_flag=FALSE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmapID=popmapID,basemapID=basemapID,
                            description="full density",density=800,npatch=20000)
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

##### b1 Qmaps #####
# 1: no autocorrelation
#qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
qmapID=1
h=-2
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="white noise",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 2: very weak autocorrelation
qmapID=2
h=-0.5
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="very weak autocorrelation",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 3: weak autocorrelation
qmapID=3
h=0
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="weak autocorrelation",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 4: strong autocorrelation
qmapID=4
h=0.5
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="strong autocorrelation",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 5: strong autocorrelation
qmapID=5
h=1
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="strong autocorrelation",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 6: very strong autocorrelation
qmapID=6
h=1.3
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="strong autocorrelation",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 7: very strong autocorrelation
qmapID=7
h=1.9
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="strong autocorrelation",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# q=0.5 everywhere
# constant q value
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
h=NA
target_dist=NA
basemap_file <- paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID)
base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
temp_rast <- 0.5*base_rast
q_rast <- c(base_rast,temp_rast)
names(q_rast) <- c("reef","q")
# save
writeRaster(q_rast,filename=paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qmapID,".tif"),overwrite=TRUE)
index_qmaps <- data.frame(qmap_id=qmapID,basemap_id=basemapID,description="same habitat quality (0.5) everywhere, 5x5km",
                          h=NA,target_dist=NA,
                          range=NA,sill=NA,
                          SSErr=NA,model=NA)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

#### b2: 1x25km, uniform ####
basemapID=2
hab_sim <- f_GenerateBasemap(x_dist=1000,y_dist=25000,resol=c(10,10),method="uniform",h=NA,prop_hab=NA,
                             make_dist_mat = FALSE,plot_flag=TRUE,experiment_folder=experiment_folder,basemap_id = basemapID)
index_basemaps <- data.frame(basemapID=basemapID,description="1x25km uniform",area_km2=25)
fwrite(index_basemaps,file="New/Maps/index_basemaps.csv",append=TRUE)

##### b2 popmaps #####
# 1: small-density (8 anems/km2)
popmapID=max(read.csv("New/Maps/index_popmaps.csv")$popmapID,0)+1
pts_out <- f_SimPtsOnMap(n_anems=25*8,inwater_dist=FALSE,samp_type = "random",plot_flag=TRUE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmapID=popmapID,basemapID=basemapID,
                            description="low density",density=8,npatch=200)
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

# 2: full-density (800 anems/km2)
popmapID=max(read.csv("New/Maps/index_popmaps.csv")$popmapID,0)+1
pts_out <- f_SimPtsOnMap(n_anems=25*800,inwater_dist=FALSE,samp_type = "random",plot_flag=FALSE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmapID=popmapID,basemapID=basemapID,
                            description="full density",density=800,npatch=20000)
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

##### b2 Qmaps ##### Need to do these later (or poach ones I've done); they take a long time
# # 1: no autocorrelation
# qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
# h=-2
# target_dist="identity"
# q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
#                            experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)
# 
# index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="white noise",
#                           h=h,target_dist=target_dist,
#                           range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
#                           SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
# fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

#### b3: Kimbe, small ####
basemapID=3
# Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
kimbe_reef_WGS <- read_sf('seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson') |>
  filter(class=="Coral/Algae")
# Project to UTM
# calculate the UTM zone
bounds <- st_bbox(kimbe_reef_WGS)
meanlon <- (bounds$xmax+bounds$xmin)/2
meanlat <- (bounds$ymax+bounds$ymin)/2
UTM_zone <- floor((meanlon+180)/6)+1 # formula for finding UTM zone based on lat/lon
UTM_hem <- ifelse(meanlat>0,32600,32700) # crs starts with 326 for northern hemisphere, 327 for southern
UTM_crs <- as.numeric(UTM_hem+UTM_zone) # build the crs to pass to st_transform
# project from WGS84 to UTM
reef_sf <- st_transform(kimbe_reef_WGS,crs=UTM_crs)
ggplot(reef_sf)+geom_sf()
st_bbox(reef_sf) # now the coordinates are in meters
kimbe_reef_area <- sum(st_area(reef_sf))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# Rasterize the shapefile to get base_rast
resol=c(10,10)
units(resol) <- "m"
template <- rast(extent=raster::extent(reef_sf),resolution=resol,crs=crs(reef_sf))
base_rast <- terra::rasterize(x=reef_sf,y=template)
# plot(base_rast,xlim=c(177000,177500),ylim=c(9400600,9400700))
# plot(base_rast)
names(base_rast) <- "lyr.1"
#plot(base_rast,col="black",xlim=c(177000,177500),ylim=c(9400600,9400700))

# Get bathymetry from marmap
kimbe_bathy <- rast("seascapes/Field data/Kimbe2-20251219174118/Bathymetry---composite-depth/bathymetry_0.tif") # RasterLayer
kimbe_bathy <- project(kimbe_bathy,base_rast) # project to UTM, to match the other files
kimbe_bathy <- raster(kimbe_bathy) # convert to older RasterLayer format to use with marmap
bathy_rast <- aggregate(kimbe_bathy,fact=2) # decrease resolution of bathymetry file
bathy_rast <- marmap::as.bathy(bathy_rast)
marmap_transmat <- marmap::trans.mat(-bathy_rast) # depths are positive, so need to take the negative of the bathymetry object (or change the range of values with arguments to trans.mat)

# save reef_sf and bathy_rast in an .RData file (base_b1.RData), and base_rast as a .tif (base_b1.tif). All in a folder called b1.
# added marmap_transmat, because we'll need it to calculate in-water distances
dir.create(path=paste0(experiment_folder,"/b",basemapID),recursive=TRUE)
save(reef_sf,bathy_rast,marmap_transmat,file=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".RData"))
writeRaster(base_rast,filename=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".tif"),overwrite=TRUE)
index_basemaps <- data.frame(basemapID=basemapID,description="Kimbe small",area_km2=as.numeric(kimbe_reef_area))
fwrite(index_basemaps,file="New/Maps/index_basemaps.csv",append=TRUE)

##### b3 popmaps #####
# full-density (800 anems/km2)
popmapID=max(read.csv("New/Maps/index_popmaps.csv")$popmapID,0)+1
pts_out <- f_SimPtsOnMap(n_anems=as.numeric(kimbe_reef_area)*800,inwater_dist=FALSE,samp_type = "random",plot_flag=TRUE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmapID=popmapID,basemapID=basemapID,
                            description="full density, small Kimbe",density=800,npatch=round(as.numeric(kimbe_reef_area)*800))
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

##### b3 qmaps #####
# low autocorrelation
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
h=-2
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="low autocorrelation, small Kimbe",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# moderate autocorrelation
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
h=0
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="moderate autocorrelation, small Kimbe",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# high autocorrelation
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
h=0.5
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="moderate autocorrelation, small Kimbe",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# very high autocorrelation
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
h=1.5
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="moderate autocorrelation, small Kimbe",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

#### b4: Kimbe, medium ####
basemapID=4
# Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
kimbe_reef_WGS <- read_sf('seascapes/Field data/Kimbe_large-20260114171835/Benthic-Map/benthic.geojson') %>%
  filter(class=="Coral/Algae")
kimbe_reef_WGS <- st_crop(kimbe_reef_WGS,xmin=150.015,xmax=150.25,ymin=-5.5,ymax=-5.05)
plot(kimbe_reef_WGS)
# kimbe_reef_WGS <- st_crop(kimbe_reef_WGS,xmin=150.015,xmax=150.25,ymin=-5.559,ymax=-5.05)
# Project to UTM
# calculate the UTM zone
bounds <- st_bbox(kimbe_reef_WGS)
meanlon <- (bounds$xmax+bounds$xmin)/2
meanlat <- (bounds$ymax+bounds$ymin)/2
UTM_zone <- floor((meanlon+180)/6)+1 # formula for finding UTM zone based on lat/lon
UTM_hem <- ifelse(meanlat>0,32600,32700) # crs starts with 326 for northern hemisphere, 327 for southern
UTM_crs <- as.numeric(UTM_hem+UTM_zone) # build the crs to pass to st_transform
# project from WGS84 to UTM
reef_sf <- st_transform(kimbe_reef_WGS,crs=UTM_crs) |>
  filter(st_geometry_type(geometry)=="POLYGON")
ggplot(reef_sf)+geom_sf()
st_bbox(reef_sf) # now the coordinates are in meters
kimbe_reef_area <- sum(st_area(reef_sf))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# Rasterize the shapefile to get base_rast
resol=c(10,10)
units(resol) <- "m"
template <- rast(extent=raster::extent(reef_sf),resolution=resol,crs=crs(reef_sf))
base_rast <- terra::rasterize(x=reef_sf,y=template)
names(base_rast) <- "lyr.1"

# Get bathymetry from marmap
kimbe_bathy <- rast("seascapes/Field data/Kimbe_large-20260114171835/Bathymetry---composite-depth/bathymetry_0.tif") # RasterLayer
kimbe_bathy <- crop(kimbe_bathy,ext(kimbe_reef_WGS))
kimbe_bathy <- project(kimbe_bathy,base_rast) # project to UTM, to match the other files
kimbe_bathy <- raster(kimbe_bathy) # convert to older RasterLayer format to use with marmap
bathy_rast <- aggregate(kimbe_bathy,fact=10) # decrease resolution of bathymetry file
bathy_rast <- marmap::as.bathy(bathy_rast)
marmap_transmat <- marmap::trans.mat(-bathy_rast) # depths are positive, so need to take the negative of the bathymetry object (or change the range of values with arguments to trans.mat)

# save reef_sf and bathy_rast in an .RData file (base_b1.RData), and base_rast as a .tif (base_b1.tif). All in a folder called b1.
# added marmap_transmat, because we'll need it to calculate in-water distances
dir.create(path=paste0(experiment_folder,"/b",basemapID),recursive=TRUE)
save(reef_sf,bathy_rast,marmap_transmat,file=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".RData"))
writeRaster(base_rast,filename=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".tif"),overwrite=TRUE)
index_basemaps <- data.frame(basemapID=basemapID,description="Kimbe medium",area_km2=as.numeric(kimbe_reef_area))
fwrite(index_basemaps,file="New/Maps/index_basemaps.csv",append=TRUE)

##### b4 popmaps #####
# low-density (16 anems/km2)
popmapID=max(read.csv("New/Maps/index_popmaps.csv")$popmapID,0)+1
pts_out <- f_SimPtsOnMap(n_anems=as.numeric(kimbe_reef_area)*16,inwater_dist=FALSE,samp_type = "random",plot_flag=TRUE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmapID=popmapID,basemapID=basemapID,
                            description="low density, medium Kimbe",density=16,npatch=round(as.numeric(kimbe_reef_area)*16))
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

# full-density (800 anems/km2)
popmapID=max(read.csv("New/Maps/index_popmaps.csv")$popmapID,0)+1
pts_out <- f_SimPtsOnMap(n_anems=as.numeric(kimbe_reef_area)*800,inwater_dist=FALSE,samp_type = "random",plot_flag=FALSE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmapID=popmapID,basemapID=basemapID,
                            description="full density, medium Kimbe",density=800,npatch=round(as.numeric(kimbe_reef_area)*800))
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

##### b4 qmaps #####
# low autocorrelation
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
h=-2
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="(supposedly) low autocorrelation, medium Kimbe",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)


# constant q value
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmapID,0)+1
h=NA
target_dist=NA

# if a filepath was specified, load the saved base map. Otherwise it's ready to go.
basemap_file <- paste0(experiment_folder,"/b",basemap_id,"/base_b",basemap_id)
base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast

nx=ncol(base_rast)
ny=nrow(base_rast)
frac_map <- matrix(data=0.5,nrow=ny,ncol=nx)
# mask out non-habitat locations
frac_map[is.na(matrix(values(base_rast),nrow=ny,byrow=TRUE))] <- NA # if they're NA's
frac_map[matrix(values(base_rast),nrow=ny,byrow=TRUE)==0] <- NA # if they're zeros

temp_rast <- rast(ext(base_rast), resolution=res(base_rast), crs = crs(base_rast))
values(temp_rast) <- frac_map
# put together with base_rast in new object
q_rast <- c(base_rast,temp_rast)
names(q_rast) <- c("reef","q")

# save
writeRaster(q_rast,filename=paste0(experiment_folder,"/b",basemap_id,"/qmap_b",basemap_id,"_q",qmap_id,".tif"),overwrite=TRUE)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="same habitat quality (0.5) everywhere, medium Kimbe",
                          h=NA,target_dist=NA,
                          range=NA,sill=NA,
                          SSErr=NA,model=NA)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)


# ##### Add min_dist and max_dist to the first 7 popmap files -- didn't do it when they were made.
# #I've updated the function now so it'll do it when they're generated.
# index_popmaps <- read.csv(file="New/Maps/index_popmaps.csv")
# for(i in 1:nrow(index_popmaps)){
#   row_i <- index_popmaps[i,]
#   # load patchdists
#   load(paste0("New/Maps/b",row_i$basemap_id,"/patchdists_b",row_i$basemap_id,"_p",row_i$popmap_id,".RData"))
#   max_dist <- max(patch_dists)
#   min_dist <- min(patch_dists)
#   # load pop.RData
#   load(paste0("New/Maps/b",row_i$basemap_id,"/pop_b",row_i$basemap_id,"_p",row_i$popmap_id,".RData"))
#   save(reef_sf,sfc_patches,hab_type,basemap_file,min_dist,max_dist,
#        file=paste0("New/Maps/b",row_i$basemap_id,"/pop_b",row_i$basemap_id,"_p",row_i$popmap_id,".RData"))
# }

#### b5: 5x5km, squares in the corners ####
basemapID=5
x_dist=5000
y_dist=5000
resol=c(10,10)
corner_length <- 1000 # side length of the corner squares, in m

# base_rast
base_rast <- rast(xmin=250000,xmax=250000+x_dist,ymin=0,ymax=y_dist,crs="EPSG:32631",resol=resol)
nx=ncol(base_rast)
ny=nrow(base_rast)
base_mat <- matrix(0,nrow=ny,ncol=nx)
corner_length_rast_x <- corner_length/resol[1]
corner_length_rast_y <- corner_length/resol[2]
base_mat[1:corner_length_rast_y,1:corner_length_rast_x] <- 1
base_mat[((ny-corner_length_rast_y):ny),((nx-corner_length_rast_x):nx)] <- 1
values(base_rast) <- base_mat

# reef_sf
bb <- as.numeric(ext(base_rast)[c(1,3,2,4)])
corner1 <- st_sfc(st_polygon(list(rbind(c(bb[1],bb[4]), # xmin, ymax
                                        c(bb[1]+corner_length,bb[4]),
                                        c(bb[1]+corner_length,bb[4]-corner_length),
                                        c(bb[1],bb[4]-corner_length),
                                        c(bb[1],bb[4])))),crs="EPSG:32631")
corner2 <- st_sfc(st_polygon(list(rbind(c(bb[3],bb[2]), # xmax,ymin
                                        c(bb[3]-corner_length,bb[2]),
                                        c(bb[3]-corner_length,bb[2]+corner_length),
                                        c(bb[3],bb[2]+corner_length),
                                        c(bb[3],bb[2])))),crs="EPSG:32631")
reef_sf <- st_sf(st_union(corner1,corner2),crs="EPSG:32631")
patch_dists <- NULL
sfc_patches <- NULL
bathy_rast <- NULL

ggplot()+
#  ggspatial::layer_spatial(base_rast)+
  geom_sf(data=reef_sf)

dir.create(path=paste0(experiment_folder,"/b",basemapID),recursive=TRUE)
save(reef_sf,patch_dists,sfc_patches,bathy_rast,file=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".RData"))
writeRaster(base_rast,filename=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".tif"),overwrite=TRUE)

index_basemaps <- data.frame(basemap_id=basemapID,description="5x5km corner patches",area_km2=2)
fwrite(index_basemaps,file="New/Maps/index_basemaps.csv",append=TRUE)

##### b5 popmaps #####
# 1: full-density (800 anems/km2)
popmapID=max(read.csv("New/Maps/index_popmaps.csv")$popmap_id,0)+1
pts_out <- f_SimPtsOnMap(n_anems=2*800,inwater_dist=FALSE,samp_type = "random",plot_flag=FALSE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmap_id=popmapID,basemap_id=basemapID,
                            description="full density",density=800,npatch=1600)
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

##### b5 Qmaps #####
# q=0.5 everywhere
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
h=NA
target_dist=NA
basemap_file <- paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID)
base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
temp_rast <- 0.5*base_rast
q_rast <- c(base_rast,temp_rast)
names(q_rast) <- c("reef","q")
# save
writeRaster(q_rast,filename=paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qmapID,".tif"),overwrite=TRUE)
index_qmaps <- data.frame(qmap_id=qmapID,basemap_id=basemapID,description="q=0.5 everywhere, 5x5km corners",
                          h=NA,target_dist=NA,
                          range=NA,sill=NA,
                          SSErr=NA,model=NA)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 1: no autocorrelation
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
h=-2
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmap_id=qmapID,basemap_id=basemapID,description="white noise",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

#### b6: 1x25km, 5 regular patches of length 1.5km ####
basemapID=6
x_dist=1000
y_dist=25000
resol=c(10,10)
patch_length <- 1500 # length of the habitable patches, in m
patch_num <- 5
gap_length <- (y_dist-(patch_length*patch_num))/(patch_num-1)

# base_rast
base_rast <- rast(xmin=250000,xmax=250000+x_dist,ymin=0,ymax=y_dist,crs="EPSG:32631",resol=resol)
nx=ncol(base_rast)
ny=nrow(base_rast)
base_mat <- matrix(0,nrow=ny,ncol=nx)
patch_length_rast <- patch_length/resol[2]
gap_length_rast <- gap_length/resol[2]

for(i in 1:patch_num){
  base_mat[(1+(i-1)*(patch_length_rast+gap_length_rast)):(((i-1)*(patch_length_rast+gap_length_rast))+patch_length_rast),1:nx] <- 1
  }
values(base_rast) <- base_mat

# reef_sf
basemap_stars <- stars::st_as_stars(base_rast[[1]])
basemap_contour <- stars::st_contour(basemap_stars,breaks=c(0.5))
reef_sf <- basemap_contour[basemap_contour$Min>0,] # pick out just the reef part for the shapefile
areareef <- st_area(reef_sf)
units(areareef) <- "km^2"

ggplot()+
  ggspatial::layer_spatial(base_rast)+
  geom_sf(data=reef_sf,alpha=0.1,color='white')

patch_dists <- NULL
sfc_patches <- NULL
bathy_rast <- NULL

dir.create(path=paste0(experiment_folder,"/b",basemapID),recursive=TRUE)
save(reef_sf,patch_dists,sfc_patches,bathy_rast,file=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".RData"))
writeRaster(base_rast,filename=paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID,".tif"),overwrite=TRUE)

index_basemaps <- data.frame(basemap_id=basemapID,description="1x25km regular patches (5, 1.5km each)",area_km2=as.numeric(areareef))
fwrite(index_basemaps,file="New/Maps/index_basemaps.csv",append=TRUE)

##### b6 popmaps #####
# 1: full-density (800 anems/km2)
popmapID=max(read.csv("New/Maps/index_popmaps.csv")$popmap_id,0)+1
pts_out <- f_SimPtsOnMap(n_anems=round(as.numeric(areareef)*800),inwater_dist=FALSE,samp_type = "random",plot_flag=FALSE,
                         experiment_folder=experiment_folder,basemap_id=basemapID,popmap_id=popmapID)
index_popmaps <- data.frame(popmap_id=popmapID,basemap_id=basemapID,
                            description="full density",density=800,npatch=round(as.numeric(areareef)*800))
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv",append=TRUE)

##### b6 Qmaps #####
# q=0.5 everywhere
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
h=NA
target_dist=NA
basemap_file <- paste0(experiment_folder,"/b",basemapID,"/base_b",basemapID)
base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
temp_rast <- 0.5*base_rast
q_rast <- c(base_rast,temp_rast)
names(q_rast) <- c("reef","q")
# save
writeRaster(q_rast,filename=paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qmapID,".tif"),overwrite=TRUE)
index_qmaps <- data.frame(qmap_id=qmapID,basemap_id=basemapID,description="q=0.5 everywhere, 1x25km patches",
                          h=NA,target_dist=NA,
                          range=NA,sill=NA,
                          SSErr=NA,model=NA)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)

# 1: no autocorrelation
qmapID <- max(read.csv("New/Maps/index_qmaps.csv")$qmap_id,0)+1
h=-2
target_dist="identity"
q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)
q_out$vgm_fit
index_qmaps <- data.frame(qmap_id=qmapID,basemap_id=basemapID,description="white noise, 1x25km patches",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)
