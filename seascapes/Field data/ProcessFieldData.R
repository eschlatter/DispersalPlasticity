.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(sfheaders)
library(ggspatial)
library(terra)
library(raster)
library(gdistance)
library(marmap)

# ############# Use actual anemone locations ################
# 
# # Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
# kimbe_reef <- read_sf('seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson') %>%
#   filter(class=="Coral/Algae")
# kimbe_reef_area <- sum(st_area(kimbe_reef))
# units(kimbe_reef_area) <- "km^2" # and convert to km^2
# 
# # Convert reef shapefile to raster for simulation
# template <- rast(extent=raster::extent(kimbe_reef),resolution=0.0002)
# kimbe_raster <- terra::rasterize(x=kimbe_reef,y=template)
# kimbe_raster <- raster::raster(kimbe_raster)
# n_patch <- sum(kimbe_raster@data@values,na.rm=TRUE) # how many patches
# 
# plot(kimbe_raster,col='black')
# points(x=anemones$Longitude,y=anemones$Latitude,col='red',pch=20,cex=0.2)
# text(x=150.07,y=-5.45,labels=paste0(n_patch," patches"))
# 
# # Anemone data from Zoe
# anemones <- read.csv('seascapes/Field data/DispersalPlasticity_anemones_metadata.csv')
# gps_pts <- group_by(anemones,Tag) %>%
#   summarize(Latitude=first(Latitude),Longitude=first(Longitude),Reef=first(Reef))
# gps_pts <- st_as_sf(gps_pts,coords=c("Longitude","Latitude"))
# st_crs(gps_pts) <- 4326 # this looks right on the map, so I'm going with it
# reefs <- group_by(anemones,Reef) %>%
#   summarize(Latitude=mean(Latitude),Longitude=mean(Longitude))
# 
# # Plot anemone locations
# ggplot(kimbe_reef)+
#   geom_sf()+
#   geom_sf(data=gps_pts,aes(color=Reef),size=0.75)+
#   geom_text_repel(data=reefs,aes(x=Longitude,y=Latitude,label=Reef,color=Reef),box.padding=0.75,size=3)+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   annotate("text",x=150.085,y=-5.458,label=paste0("Reef area = ",round(kimbe_reef_area,2),units(kimbe_reef_area)),hjust="left",size=3)+
#   annotation_scale() # ggspatial
# 
# # Empirical dispersal distances for A. percula (from Almany et al 2017)
# mean_disp_dist_09 <- 10 # 15km in 2009
# mean_disp_dist_11 <- 10 # 10km in 2011
# lambda_09 <- 18.9 # exponential distribution with mean=18.9km
# lambda_11 <- 13.3 # exponential distribution with mean = 13.3km
# # that study used all of Kimbe Bay, diameter ~150km
# # our study has a long dimension of ~7km

############# Real anemone locations on restricted reef map ################
# Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
kimbe_reef <- read_sf('seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson') %>%
  filter(class=="Coral/Algae")
kimbe_reef_area <- sum(st_area(kimbe_reef))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# Get bathymetry, also from Allen Coral Atlas (depth in CENTIMETERS)
kimbe_bathy <- raster("seascapes/Field data/Kimbe2-20251219174118/Bathymetry---composite-depth/bathymetry_0.tif") # RasterLayer
kimbe_bathy <- marmap::as.bathy(kimbe_bathy)

# Anemone data from Zoe
anemones <- read.csv('seascapes/Field data/DispersalPlasticity_anemones_metadata.csv')
gps_pts <- group_by(anemones,Tag) %>%
  summarize(Latitude=first(Latitude),Longitude=first(Longitude),Reef=first(Reef))
reefs <- group_by(anemones,Reef) %>%
  summarize(Latitude=mean(Latitude),Longitude=mean(Longitude))
anemone_spots_df <- gps_pts %>%
  mutate(depth=marmap::get.depth(kimbe_bathy,gps_pts$Longitude,gps_pts$Latitude,locator=FALSE)$depth) 
gps_pts <- st_as_sf(gps_pts,coords=c("Longitude","Latitude"))
st_crs(gps_pts) <- 4326 # this looks right on the map, so I'm going with it

ggplot(gps_pts)+
  geom_contour_filled(data=marmap::as.xyz(kimbe_bathy),aes(x=V1,y=V2,z=V3),alpha=0.5)+
  geom_sf(data=kimbe_reef,fill='black',color='black',alpha=0.2)+  
  geom_sf(color='red',size=1)+
  annotation_scale()+
  theme_minimal()

############# Simulate anemone locations on restricted reef map ################
# Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
kimbe_reef <- read_sf('seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson') %>%
  filter(class=="Coral/Algae")
kimbe_reef_area <- sum(st_area(kimbe_reef))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# Get bathymetry, also from Allen Coral Atlas (depth in CENTIMETERS)
kimbe_bathy <- raster("seascapes/Field data/Kimbe2-20251219174118/Bathymetry---composite-depth/bathymetry_0.tif") # RasterLayer
kimbe_bathy <- marmap::as.bathy(kimbe_bathy)

# Place anemones randomly on reef
anemone_spots <- st_sample(kimbe_reef,5)
anemone_spots_df <- sfc_to_df(anemone_spots)%>%
  mutate(depth=marmap::get.depth(kimbe_bathy,x,y,locator=FALSE)$depth)

ggplot(anemone_spots)+
  geom_contour_filled(data=marmap::as.xyz(kimbe_bathy),aes(x=V1,y=V2,z=V3),alpha=0.5)+
  geom_sf(data=kimbe_reef,fill='black',color='black',alpha=0.2)+  
  geom_sf(color='red',size=1)+
  annotation_scale()+
  theme_minimal()

# distances among anemones
trans1 <- trans.mat(-kimbe_bathy) # depths are positive, so need to take the negative of the bathymetry object (or change the range of values with arguments to trans.mat)
anemone_dists <- lc.dist(trans1,anemone_spots_df[,c("x","y")],res='dist') #distances are in km
anemone_dists_mat <- matrix(0,nrow(anemone_spots_df),nrow(anemone_spots_df)) # convert dist into matrix form
anemone_dists_mat[lower.tri(anemone_dists_mat,diag=FALSE)] <- anemone_dists
# visualize an in-water path
anemone_paths <- lc.dist(trans1,anemone_spots_df[,c("x","y")],res='path')
# Get one trajectory for plotting
a <- as.data.frame(anemone_paths[[3]])
b <- sf_linestring(a,x="x",y="y")
st_crs(b) <- st_crs(anemone_spots)
ggplot(anemone_spots)+
  geom_contour_filled(data=marmap::as.xyz(kimbe_bathy),aes(x=V1,y=V2,z=V3),alpha=0.5)+
  geom_sf(data=kimbe_reef,fill='black',color='black',alpha=0.2)+  
  geom_sf(color='red',size=1)+
  annotation_scale()+
  theme_minimal()+
  geom_sf(data=b,color="yellow")


############# Simulate anemone locations on larger reef map ################
# Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
kimbe_reef <- read_sf('seascapes/Field data/Kimbe_large-20260114171835/Benthic-Map/benthic.geojson') %>%
  filter(class=="Coral/Algae")
kimbe_reef_area <- sum(st_area(kimbe_reef))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# Get bathymetry from marmap
kimbe_bathy <- raster("seascapes/Field data/Kimbe_large-20260114171835/Bathymetry---composite-depth/bathymetry_0.tif") # RasterLayer
kimbe_bathy <- marmap::as.bathy(kimbe_bathy)

# Place anemones randomly on reef
anemone_spots <- st_sample(kimbe_reef,5)
anemone_spots_df <- sfc_to_df(anemone_spots)%>%
  mutate(depth=marmap::get.depth(kimbe_bathy,x,y,locator=FALSE)$depth) %>%
  filter(depth>0)

g <- ggplot(anemone_spots)+
  #geom_contour_filled(data=marmap::as.xyz(kimbe_bathy),aes(x=V1,y=V2,z=V3),alpha=0.5)+
  geom_sf(data=kimbe_reef,fill='black',color='black',alpha=0.2)+  
  geom_sf(color='red',size=1)+
  annotation_scale()+
  theme_minimal()

# distances among anemones
trans1 <- trans.mat(-kimbe_bathy) # depths are positive, so need to take the negative of the bathymetry object (or change the range of values with arguments to trans.mat)
save(trans1,file="kimbe_trans_mat.RData")
sponge_paths <- lc.dist(trans1,anemone_spots_df[,c("x","y")],res='path')

# Get one trajectory for plotting
a <- as.data.frame(sponge_paths[[3]])
b <- sf_linestring(a,x="x",y="y")
st_crs(b) <- st_crs(anemone_spots)

ggplot(anemone_spots)+
  geom_contour_filled(data=marmap::as.xyz(kimbe_bathy),aes(x=V1,y=V2,z=V3),alpha=0.5)+
  geom_sf(data=kimbe_reef,fill='black',color='black',alpha=0.2)+  
  geom_sf(color='red',size=1)+
  annotation_scale()+
  theme_minimal()+
  geom_sf(data=b,color="yellow")

