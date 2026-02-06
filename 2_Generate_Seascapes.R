source('0_Setup.R')
library(png)

## landscape configuration from real map, with simulated (fractal) K's
# ## import Belize map and take a little subset
# bzemap <- readPNG("map_SLiM.png")
# bzemap_sub <- bzemap[,,1]
# bzemap_sub[bzemap_sub==1] <- 1000
# bzemap_sub[bzemap_sub<1] <- TRUE
# bzemap_sub[bzemap_sub==1000] <- FALSE
# base_map <- bzemap_sub[1250:1280,555:585]
# image.plot(base_map)
# sum(base_map)
# save(base_map,file='seascapes/bze_map_sub.RData')
load('seascapes/bze_map_sub.RData')
bze_out <- f_GenerateMapWithK(base_map, K_range=c(1,15), h = 0.8, plot_flag=TRUE)


## fractal landscapes, 33x33 (~325 patches)
frac_out <- f_GenerateMapWithK(K_range=c(1,15), h=1.8, k=5, p=0.3, h_base=0.8, plot_flag=TRUE)
frac_out <- f_GenerateMapWithK(K_range=c(1,15), h=0.2, k=5, p=0.3, h_base=0.8, plot_flag=TRUE)


## Create a basemap. We'll keep this constant for awhile, and vary habitat on it.
a <- f_GenerateMapWithK(base_map=NULL,K_range=c(3,10),h=0.2,k=5,p=0.3,h_base=0.8,plot_flag=TRUE)
basemap_1 <- matrix(0,nrow=a$ny,ncol=a$nx)
for(i in 1:nrow(a$patch_locations)) basemap_1[a$patch_locations$y[i],a$patch_locations$x[i]] <- 1
save(basemap_1,file='seascapes/basemap1.RData')

#################### Create base map (to pass to function f_SimPtsOnMap) ##########################

# Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
#'seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson'
kimbe_reef <- read_sf('seascapes/Field data/Kimbe_large-20260114171835/Benthic-Map/benthic.geojson') %>%
  filter(class=="Coral/Algae")
kimbe_reef_area <- sum(st_area(kimbe_reef))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# Get bathymetry from marmap
# 'seascapes/Field data/Kimbe2-20251219174118/Bathymetry---composite-depth/bathymetry_0.tif'
kimbe_bathy1 <- raster("seascapes/Field data/Kimbe_large-20260114171835/Bathymetry---composite-depth/bathymetry_0.tif") # RasterLayer
kimbe_bathy1 <- aggregate(kimbe_bathy1,fact=20) # decrease resolution of bathymetry file
kimbe_bathy <- marmap::as.bathy(kimbe_bathy1)

# create base_matrix
template_rast <- raster(ext=kimbe_bathy1@extent,res=c(0.0005,0.0005))
base_rast <- st_rasterize(kimbe_reef,template=st_as_stars(template_rast),align=TRUE)
base_matrix <- base_rast$ID
base_matrix_coords <- st_coordinates(base_rast) %>% rename(lon=x,lat=y)
coord_inds <- expand.grid(y=1:nrow(base_matrix),x=1:ncol(base_matrix))
base_matrix_coords <- cbind(base_matrix_coords,coord_inds)

transmat_20x <- trans.mat(-kimbe_bathy) # depths are positive, so need to take the negative of the bathymetry object (or change the range of values with arguments to trans.mat)

reef_sf <- kimbe_reef
bathy_raster <- kimbe_bathy
marmap_transmat <- transmat_20x

hab_sim=list(base_matrix=base_matrix,base_matrix_coords=base_matrix_coords,reef_sf=reef_sf,bathy_rast=bathy_raster,crs=crs(kimbe_bathy1))

save(base_matrix,reef_sf,bathy_raster,marmap_transmat,file="seascapes/kimbe_large_20x.RData")