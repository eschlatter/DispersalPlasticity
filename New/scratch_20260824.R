.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
library(data.table)
library(parallel)
library(terra)
library(sf)
library(calculus)

habID=2
load(paste0("New/Maps/habfiles/hab_",habID,".RData"))
base_rast <- rast(paste0("New/Maps/habfiles/hab_",habID,".tif"))
raster_box <- st_as_sfc(st_bbox(base_rast))

circs <- st_buffer(hab_params$sfc_patches,dist=50)
st_area(circs[[1]]) # this is approximately pi*50^2, so the circle around the points is correct
circs_union <- st_union(circs) # put them all together
circs_crop <- st_intersection(circs_union,raster_box) # and crop at the edge of the map

# calculate reachable area and plot the map
area_prop <- st_area(circs_crop)/as.numeric(expanse(base_rast$reef)$area)
ggplot(circs_crop)+geom_sf()+theme_minimal()+labs(title=paste0("habID=",habID,": ",round(area_prop,3)," of reef area"))


