source('0_Setup.R')

################# create a hab_params object to pass to simulation, do each of the following: ############################
# 1. Generate a base map
x_dist=500000
y_dist=500000
resol=c(100000,100000)
#resol=c(0.00008,0.00008)
base_h=0.8
prop_hab=0.2
basemap_file="seascapes/2026_03_24/50x50km_res=1km"
basemap_file=NULL
hab_sim <- f_GenerateBasemap(x_dist=x_dist,y_dist=y_dist,resol=resol,method="uniform",h=base_h,prop_hab=prop_hab,
                             make_dist_mat = TRUE,plot_flag=TRUE,basemap_file=basemap_file)
reef_area <- st_area(hab_sim$reef_sf)
units(reef_area)='km^2'

# 2. Create a habitat quality layer
qmap_file="seascapes/2026_03_24/50x50_test"
q_range=c(2,30)
q_autocorr=0.9
qual_out <- f_GenerateHabQual(base_rast=basemap_file,q_range,q_autocorr,binary=TRUE,plot_flag=TRUE,qmap_file = qmap_file)

# 3. Do EITHER 3a (hab_type=grid) OR 3b (hab_type=points)

# 3a. Create a layer that assigns each reef grid square a carrying capacity (# of individuals it can hold)
popmap_file="seascapes/2026_03_04/popmap_1_grid"
K_range <- c(1,1)
K_autocorr <- 1.8
K_out <- f_GenerateK(base_rast=basemap_file,K_range=K_range,K_autocorr=K_autocorr,
                     plot_flag = TRUE,popmap_file = popmap_file)

# 3b. Place points at random on the reef, each of which represents the habitat of a single individual
popmap_file="seascapes/2026_03_24/popmap_50x50km"
#n_anems=round(drop_units(reef_area*638))
n_anems=1000
inwater_dist=FALSE
pts_out <- f_SimPtsOnMap(basemap_file = basemap_file,n_anems=n_anems,inwater_dist=inwater_dist,
                         samp_type = "random",popmap_file=popmap_file,plot_flag=TRUE)

# 4. Put everything together
hab_file="seascapes/2026_03_04/hab_7_pt"
nav_rad <- as_units(0.05,'km')
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
  scale_fill_continuous(palette = 'BluGrn',na.value = "grey")+
  annotation_scale()+
  #geom_sf(data=sfc_patches,pch=21,fill='white',size=1)+
  labs(fill="Habitat\nquality")



# in progress:
#################### Load in base map (instead of simulating one in f_GenerateBasemap) ##########################
# need to update this so it outputs all the necessary objects

# reef_sf: Get reef shapefile(s) of field data collection area (downloaded from Allen Coral Atlas)
sf_path <- 'seascapes/Field data/Kimbe2-20251219174118/Benthic-Map/benthic.geojson'
sf_path <- 'seascapes/Field data/Kimbe_large-20260114171835/Benthic-Map/benthic.geojson'
kimbe_reef <- read_sf(sf_path) %>%
  filter(class=="Coral/Algae")
kimbe_reef_area <- sum(st_area(kimbe_reef))
units(kimbe_reef_area) <- "km^2" # and convert to km^2

# base_rast: make SpatRaster from sf
er <- rast(ext(kimbe_reef), resolution=c(0.0008,0.0008), crs = crs(kimbe_reef))

### how to rasterize?

## 1. default
base_rast_default <- rasterize(vect(kimbe_reef),er)

ggplot()+
  ggspatial::layer_spatial(base_rast_default$layer)+
  scale_fill_continuous(name="reef",na.value = "lightblue")+
  annotation_scale()+
  lims(x=c(150.06,150.1),y=c(-5.45,-5.38))+
  labs(title=paste0("default, cells=",sum(!is.na(as.matrix(base_rast_default$layer)))))

## 2. a cell is reef if the polygon touches it at all
base_rast_touch <- rasterize(vect(kimbe_reef),er,touches=TRUE)

ggplot()+
  ggspatial::layer_spatial(base_rast_touch$layer)+
  scale_fill_continuous(name="reef",na.value = "lightblue")+
  annotation_scale()+
  lims(x=c(150.06,150.1),y=c(-5.45,-5.38))+
  labs(title=paste0("touch, cells=",sum(!is.na(as.matrix(base_rast_touch$layer)))))

## 3. a cell is reef if x% of the cell is covered by reef polygon
## with frac_inhabit=0.15, this is pretty good.
frac_inhabit=0.15
base_rast_frac <- rasterize(vect(kimbe_reef),er,cover=TRUE)
base_rast_frac$layer[base_rast_frac$layer<frac_inhabit] <- NA
base_rast_frac$layer[base_rast_frac$layer>=frac_inhabit] <- 1

ggplot()+
  ggspatial::layer_spatial(base_rast_frac$layer)+
  scale_fill_continuous(name="reef",na.value = "lightblue")+
  annotation_scale()+
  lims(x=c(150.06,150.1),y=c(-5.45,-5.38))+
  labs(title=paste0("frac=",frac_inhabit,", cells=",sum(!is.na(as.matrix(base_rast_frac$layer)))))

### don't forget to do this at the end, for input into f_GenerateHabQual
base_rast$layer[is.na(base_rast$layer)] <- 0

ggplot()+
  ggspatial::layer_spatial(base_rast$layer)+
  scale_fill_continuous(name="reef",na.value = "lightblue")+
  annotation_scale()+
  lims(x=c(150.06,150.1),y=c(-5.45,-5.38))

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