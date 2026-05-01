source('0_Setup.R')
basemap_file="seascapes/2026_03_30/1x25km_res=10m"
qmap_file="q0.1"
popmap_file="pop_density800"
hab_file="hab_2"
nav_rad <- as_units(0.05,'km')
make_hab_out <- f_MakeHabitat(nav_rad=nav_rad,qmap_file=qmap_file,popmap_file = popmap_file,basemap_file=basemap_file,
                              overlap_method="simple",hab_file = hab_file)
