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

# reef_sf: Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
sf_path <- 'seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson'
#sf_path <- 'seascapes/Field data/Kimbe_large-20260114171835/Benthic-Map/benthic.geojson'
kimbe_reef <- read_sf(sf_path) %>%
  filter(class=="Coral/Algae")
kimbe_reef_area <- sum(st_area(kimbe_reef))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# base_rast: make SpatRaster from sf
er <- rast(ext(kimbe_reef), resolution=c(0.0005,0.0005), crs = crs(kimbe_reef))
base_rast <- rasterize(vect(kimbe_reef),er)

# bathy_raster: Get bathymetry for marmap
bathy_path <- 'seascapes/Field data/Kimbe2-20251219174118/Bathymetry---composite-depth/bathymetry_0.tif'
#bathy_path <- "seascapes/Field data/Kimbe_large-20260114171835/Bathymetry---composite-depth/bathymetry_0.tif"
kimbe_bathy1 <- raster(bathy_path) # RasterLayer
kimbe_bathy1 <- aggregate(kimbe_bathy1,fact=20) # decrease resolution of bathymetry file
kimbe_bathy <- marmap::as.bathy(kimbe_bathy1)

#transmat_20x <- trans.mat(-kimbe_bathy) # depths are positive, so need to take the negative of the bathymetry object (or change the range of values with arguments to trans.mat)

reef_sf <- kimbe_reef
base_rast <- base_rast
bathy_raster <- kimbe_bathy
#marmap_transmat <- transmat_20x

hab_sim=list(base_rast=base_rast,reef_sf=reef_sf,bathy_rast=bathy_raster)

save(base_matrix,reef_sf,bathy_raster,marmap_transmat,file="seascapes/kimbe_large_20x.RData")



################ how reef-like can we make the fractal landscapes? #########################
source('0_Setup.R')
fn_slice_map <- function(frac_map,min_val,max_val){
  peaks <- which(frac_map>max_val)
  valleys <- which(frac_map<(min_val))
  frac_map_cut <- frac_map
  frac_map_cut[peaks] <- NA
  frac_map_cut[valleys] <- NA
  frac_map_cut[!is.na(frac_map_cut)] <- 1
  image.plot(frac_map_cut)
  print(paste("Number of habitat patches: ",sum(frac_map_cut,na.rm=TRUE)))
  return(frac_map_cut)
}

k=10 # dimensions of landscape
h=1.2 # how clumped is it? 0<h<2
frac_map <- fracland(k=k,h=h,binary=FALSE,plotflag=TRUE)
a <- fn_slice_map(frac_map,quantile(frac_map,0.3),quantile(frac_map,0.4))
fn_slice_map(frac_map,quantile(frac_map,0.6),quantile(frac_map,0.8))
fn_slice_map(frac_map,quantile(frac_map,0.2),quantile(frac_map,0.4))

h=0.7 # how clumped is it? 0<h<2
frac_map <- fracland(k=k,h=h,binary=FALSE,plotflag=FALSE)
fn_slice_map(frac_map,quantile(frac_map,0.3),quantile(frac_map,0.7))

################# all parts of seascape generating process ############################

# f_GenerateBasemap
x_dist=34000
y_dist=34000
units(x_dist) <- 'm'
units(y_dist) <- 'm'
resol=c(0.005,0.005)
h=0.9
prop_hab=0.2
hab_sim <- f_GenerateBasemap(x_dist=x_dist,y_dist=y_dist,resol=resol,method="fractal",h=0.7,prop_hab=prop_hab,make_dist_mat = TRUE)
plot(hab_sim$base_rast)
base_rast <- hab_sim$base_rast
bathy_rast <- hab_sim$bathy_rast
reef_sf <- hab_sim$reef_sf
patch_dists <- hab_sim$patch_dists
sfc_patches <- hab_sim$sfc_patches
plot(base_rast)
plot.bathy(bathy_rast)
ggplot(reef_sf)+geom_sf()

# f_GenerateHabQual
q_range=c(2,11)
q_autocorr=0.7
qual_out <- f_GenerateHabQual(base_rast,q_range,q_autocorr,plot_flag=TRUE)
qual_out$q_rast -> q_rast
plot(q_rast)

# f_GenerateK
# use this for hab_type=grid
K_range <- c(3,10)
K_autocorr <- 0
K_out <- f_GenerateK(base_rast,K_range=K_range,K_autocorr=K_autocorr,plot_flag = TRUE)
K_rast <- K_out$K_rast

# f_SimPtsOnMap
# use this for hab_type=points
pts_out <- f_SimPtsOnMap(reef_sf,base_rast,n_anems=50,inwater_dist=FALSE,show_map=TRUE)
K_rast=pts_out$K_rast
sfc_patches=pts_out$sfc_patches
patch_dists=pts_out$patch_dists

# f_MakeHabitat
nav_rad <- 0.5
units(nav_rad) <- 'km'
circs=st_buffer(sfc_patches,dist=nav_rad)
make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,q_rast,K_rast,
                              patch_dists,sfc_patches,reef_sf,overlap_method="simple")
list2env(x=make_hab_out,envir=environment())

ggplot(reef_sf)+geom_sf()+geom_sf(data=sfc_patches)+geom_sf(data=circs,alpha=0)