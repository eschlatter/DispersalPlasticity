.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
source("New/functions.R")
library(data.table)
library(terra)
library(raster)
library(sf)
library(ggplot2)
library(ggspatial)
library(dplyr)

experiment_folder <- "New/Maps"

# Make indices for basemaps, popmaps, and qmaps. Each time a new one is generated, put it in the index.
index_basemaps <- data.frame(basemapID=integer(),description=character(),area_km2=double())
fwrite(index_basemaps,file="New/Maps/index_basemaps.csv")

index_popmaps <- data.frame(popmapID=integer(),basemapID=integer(),description=character(),density=double(),npatch=integer())
fwrite(index_popmaps,file="New/Maps/index_popmaps.csv")

index_qmaps <- data.frame(qmapID=integer(),basemapID=integer(),description=character(),
                          h=double(),target_dist=character(),
                          range=double(),sill=double(),SSErr=double(),model=character())
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv")

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



q_out <- f_GenerateHabQual(q_autocorr=h,target_dist=target_dist,plot_flag=TRUE,
                           experiment_folder=experiment_folder,basemap_id=basemapID,qmap_id=qmapID)

index_qmaps <- data.frame(qmapID=qmapID,basemapID=basemapID,description="(supposedly) low autocorrelation, medium Kimbe",
                          h=h,target_dist=target_dist,
                          range=q_out$vgm_fit$range,sill=q_out$vgm_fit$sill,
                          SSErr=q_out$vgm_fit$SSErr,model=q_out$vgm_fit$model)
fwrite(index_qmaps,file="New/Maps/index_qmaps.csv",append=TRUE)
