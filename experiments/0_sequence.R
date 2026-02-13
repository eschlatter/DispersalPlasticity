source('0_Setup.R')

# f_GenerateBasemap
x_dist=5000
y_dist=5000
resol=c(0.0005,0.0005)
h=0.9
prop_hab=0.2
hab_sim <- f_GenerateBasemap(x_dist=x_dist,y_dist=y_dist,resol=resol,method="fractal",h=0.7,prop_hab=prop_hab)
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
K_range <- c(3,10)
K_autocorr <- 0
K_out <- f_GenerateK(base_rast,K_range=K_range,K_autocorr=K_autocorr,plot_flag = TRUE)
K_rast <- K_out$K_rast

# f_SimPtsOnMap
pts_out <- f_SimPtsOnMap(reef_sf,base_rast,n_anems=50,inwater_dist=FALSE,show_map=TRUE)
K_rast=pts_out$K_rast
sfc_patches=pts_out$sfc_patches
patch_dists=pts_out$patch_dists

# f_MakeHabitat
make_hab_out <- f_MakeHabitat(nav_rad=quantile(patch_dists,0.1),q_rast,K_rast,
                               patch_dists,sfc_patches,reef_sf,overlap_method="simple")
list2env(x=make_hab_out,envir=environment())
