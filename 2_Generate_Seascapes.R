source('0_Setup.R')
# library(png) # check if we need this

#################### Create base map (to pass to the rest of the seascape-generating pipeline) ##########################

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

################# all parts of seascape-generating pipeline ############################
source('0_Setup.R')
# f_GenerateBasemap
# simulated basemap (or we can use a real one; see above)
# x_dist=34000 --> raster dimensions 61x61
# x_dist=17000 --> raster dimensions 31x31

x_dist=17000
y_dist=17000
units(x_dist) <- 'm'
units(y_dist) <- 'm'
resol=c(0.005,0.005)
h=0.9
prop_hab=0.2
hab_sim <- f_GenerateBasemap(x_dist=x_dist,y_dist=y_dist,resol=resol,method="fractal",h=h,prop_hab=prop_hab,make_dist_mat = TRUE,plot_flag=TRUE)
base_rast <- hab_sim$base_rast
bathy_rast <- hab_sim$bathy_rast
reef_sf <- hab_sim$reef_sf
patch_dists <- hab_sim$patch_dists
sfc_patches <- hab_sim$sfc_patches

# f_GenerateHabQual
q_range=c(5,50)
q_autocorr=0.9
qual_out <- f_GenerateHabQual(base_rast,q_range,q_autocorr,plot_flag=TRUE)
q_rast <- qual_out$q_rast

# f_GenerateK
# use this for hab_type=grid
K_range <- c(3,10)
K_autocorr <- 1.4
K_out <- f_GenerateK(base_rast,K_range=K_range,K_autocorr=K_autocorr,plot_flag = TRUE)
K_rast <- K_out$K_rast
hab_type=K_out$hab_type

# f_SimPtsOnMap
# use this for hab_type=points
pts_out <- f_SimPtsOnMap(reef_sf,base_rast,n_anems=30,inwater_dist=FALSE,plot_flag=TRUE)
K_rast=pts_out$K_rast
sfc_patches=pts_out$sfc_patches
patch_dists=pts_out$patch_dists
hab_type=pts_out$hab_type

# f_MakeHabitat
nav_rad <- 0.5
units(nav_rad) <- 'km'
circs=st_buffer(sfc_patches,dist=nav_rad)
make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,q_rast,K_rast,
                              patch_dists,sfc_patches,reef_sf,hab_type,overlap_method="simple")
list2env(x=make_hab_out,envir=environment())

# plot(hab_rast)
ggplot(reef_sf)+geom_sf()+geom_sf(data=sfc_patches)+geom_sf(data=circs,alpha=0)+annotation_scale()

alpha <- 0.9
theta <- 1.5
f_plot_gamma(alpha,theta,kern_xlim = 17)
f_GetConnectivityMatrix_parallel(rep(alpha,3),rep(theta,3),drop_units(patch_dists),drop_units(patch_angles),overlap_discount,drop_units(nav_rad),numCores=1)
