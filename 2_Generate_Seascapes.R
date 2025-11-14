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