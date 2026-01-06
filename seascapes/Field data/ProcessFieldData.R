library(dplyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(ggspatial)
library(terra)

# Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
kimbe_reef <- read_sf('seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson') %>%
  filter(class=="Coral/Algae")
kimbe_reef_area <- sum(st_area(kimbe_reef))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# Convert reef shapefile to raster for simulation
template <- rast(extent=raster::extent(kimbe_reef),resolution=0.0002)
kimbe_raster <- terra::rasterize(x=kimbe_reef,y=template)
kimbe_raster <- raster(kimbe_raster)
n_patch <- sum(kimbe_raster@data@values,na.rm=TRUE) # how many patches

plot(kimbe_raster,col='black')
points(x=anemones$Longitude,y=anemones$Latitude,col='red',pch=20,cex=0.2)
text(x=150.07,y=-5.45,labels=paste0(n_patch," patches"))

# Anemone data from Zoe
anemones <- read.csv('seascapes/Field data/DispersalPlasticity_anemones_metadata.csv')
gps_pts <- group_by(anemones,Tag) %>%
  summarize(Latitude=first(Latitude),Longitude=first(Longitude),Reef=first(Reef))
reefs <- group_by(anemones,Reef) %>%
  summarize(Latitude=mean(Latitude),Longitude=mean(Longitude))

# Plot anemone locations
ggplot(kimbe_reef)+
  geom_sf()+
  geom_point(data=gps_pts,aes(x=Longitude,y=Latitude,color=Reef),size=0.75)+
  geom_text_repel(data=reefs,aes(x=Longitude,y=Latitude,label=Reef,color=Reef),box.padding=0.75,size=3)+
  theme_minimal()+
  theme(legend.position = "none")+
  annotate("text",x=150.085,y=-5.458,label=paste0("Reef area = ",round(kimbe_reef_area,2),units(kimbe_reef_area)),hjust="left",size=3)+
  annotation_scale() # ggspatial

# Empirical dispersal distances for A. percula (from Almany et al 2017)
mean_disp_dist_09 <- 10 # 15km in 2009
mean_disp_dist_11 <- 10 # 10km in 2011
lambda_09 <- 18.9 # exponential distribution with mean=18.9km
lambda_11 <- 13.3 # exponential distribution with mean = 13.3km
# that study used all of Kimbe Bay, diameter ~150km
# our study has a long dimension of ~7km
