overlap_area <- seq(from=0,to=5,by=0.1)

low_est <- 1/(1+overlap_area)
high_est <- 1-(overlap_area/2)

plot(overlap_area,low_est,type='l',col='red')
lines(overlap_area,high_est,col='blue')

library(sf)
anemones <- read.csv('seascapes/Field data/DispersalPlasticity_anemones_metadata.csv')
patch_locations <- group_by(anemones,Tag) %>%
  summarize(Latitude=first(Latitude),Longitude=first(Longitude)) 
patch_locations <- patch_locations %>%
  mutate(id=1:nrow(patch_locations)) %>%
  dplyr::select(id,Latitude,Longitude)
patch_sf <- st_as_sf(patch_locations,coords=c("Longitude","Latitude"))
st_crs(patch_sf) <- 4326

patch_dists <- st_distance(patch_sf)
str(patch_dists)
units(patch_dists) <- "km"
image.plot(matrix(as.numeric(patch_dists),nrow=230))

## get the intersection areas with sf
circs=st_buffer(patch_sf,dist=100) # dist is the diameter of the circle, in meters
onecirc_area=st_area(circs[1,])

all_overlaps <- lapply(1:3,f_FindOverlapAreas)
all_overlaps <- unlist(all_overlaps)

ggplot(patch_sf)+
  geom_sf(size=0.2)+
  geom_sf(data=circs,fill=NA)+
  annotation_scale()

intersect_areas <- st_intersection(circs[1,],circs[-1,])

ggplot(intersect_areas)+
  geom_sf()

ggplot(intersect_areas[1,])+
  geom_sf(data=circs[1,])+
  geom_sf(data=intersect_areas[1,],fill='black')








patch_locations$K_i <- 1

nx=max(patch_locations$x)-min(patch_locations$x)
ny=max(patch_locations$y)-min(patch_locations$y)
v_alphas <- 1:5
v_thetas <- 1:5
conn_out=FALSE