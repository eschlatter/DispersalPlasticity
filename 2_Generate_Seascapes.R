source('0_Setup.R')

################# create a hab_params object to pass to simulation, do each of the following: ############################
# 1. Generate a base map
x_dist=as_units(2000,'m')
y_dist=as_units(2000,'m')
resol=c(0.00008,0.00008)
base_h=1.1
prop_hab=.1
basemap_file="seascapes/2026_02_26/basemap_3"
hab_sim <- f_GenerateBasemap(x_dist=x_dist,y_dist=y_dist,resol=resol,method="fractal",h=base_h,prop_hab=prop_hab,
                             make_dist_mat = FALSE,plot_flag=TRUE,basemap_file=basemap_file)
reef_area <- st_area(hab_sim$reef_sf)
units(reef_area)='km^2'

# 2. Create a habitat quality layer
qmap_file="seascapes/2026_02_26/qmap_3_1"
q_range=c(2,30)
q_autocorr=0.8
qual_out <- f_GenerateHabQual(base_rast=basemap_file,q_range,q_autocorr,plot_flag=TRUE,qmap_file = qmap_file)

# 3. Do EITHER 3a (hab_type=grid) OR 3b (hab_type=points)

# 3a. Create a layer that assigns each reef grid square a carrying capacity (# of individuals it can hold)
popmap_file="seascapes/2026_02_26/popmap_2_1_grid"
K_range <- c(2,10)
K_autocorr <- 1.8
K_out <- f_GenerateK(base_rast=basemap_file,K_range=K_range,K_autocorr=K_autocorr,
                     plot_flag = TRUE,popmap_file = popmap_file)

# 3b. Place points at random on the reef, each of which represents the habitat of a single individual
popmap_file="seascapes/2026_02_26/popmap_3_1_pt"
n_anems=230
inwater_dist=FALSE
pts_out <- f_SimPtsOnMap(basemap_file = basemap_file,n_anems=n_anems,inwater_dist=inwater_dist,popmap_file=popmap_file,plot_flag=TRUE)

# 4. Put everything together
hab_file="seascapes/2026_02_26/hab_3_1_pt"
nav_rad <- as_units(0.5,'km')
make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,qmap_file=qmap_file,popmap_file = popmap_file,overlap_method="simple",hab_file = hab_file)


## look at the output
load(file=paste0(hab_file,".RData"))
list2env(x=hab_params,envir=environment())
hab_rast <- rast(paste0(hab_file,".tif")) # load hab_rast
plot(hab_rast)

circs=st_buffer(sfc_patches,dist=nav_rad)
ggplot(reef_sf)+
  geom_sf()+
  geom_sf(data=sfc_patches)+
  geom_sf(data=circs,alpha=0)+
  annotation_scale()

# for points-type habitat:
ggplot()+
  ggspatial::layer_spatial(hab_rast$q)+
  scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
  annotation_scale()+
  geom_sf(data=sfc_patches,pch=21,fill='white')+
  labs(title="Habitat quality")



# in progress:
#################### Load in base map (instead of simulating one in f_GenerateBasemap) ##########################
# need to update this so it outputs all the necessary objects

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