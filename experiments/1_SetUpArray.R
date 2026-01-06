source('0_Setup.R')
## set up runs

## h values
h_vals <- round(seq(from=-1.8,to=1.8,by=0.6),2)
n_basemaps <- 3 # number of basemaps to use
n_k_per_base <- 3 # number of Kmaps to produce, for each basemap + h value combination
n_rep_per_k <- 3 # number of identical (up to stochasticity) reps to do of each map


## BaseMap:
## reef vs open water maps. For now, all generated with fractal landscape, h=0.9. Could also try varying h and/or using real maps or another method.
## Let's make 3 of them.
BaseMaps <- list()
for(i in 1:n_basemaps){
  hab <- f_GenerateMapWithK(base_map=NULL,K_range=c(1,1),h=0,k=6,p=0.3,h_base=0.9,plot_flag=TRUE)
  basemap_i <- matrix(0,nrow=hab$ny,ncol=hab$nx)
  for(i in 1:nrow(hab$patch_locations)) basemap_i[hab$patch_locations$y[i],hab$patch_locations$x[i]] <- 1
  BaseMaps <- c(BaseMaps,list(basemap_i))
}

## Kmap:
## habitat maps 
KMaps <- list()
BaseKMapIndex <- data.frame(Kmap=1:(length(BaseMaps)*length(h_vals)*n_k_per_base+length(BaseMaps)),BaseMap=NA,h=NA) #keep track of which base map belongs with each kmap
ind <- 0
for(basemap in 1:length(BaseMaps)){  ## For each base map
  for(h_i in h_vals){  ## Loop over all values of h we want to try
    for(i in 1:n_k_per_base){     ## Let's do 3 Kmaps of each h value per base map
      ind <- ind+1
      kmap_i <- f_GenerateMapWithK(base_map=BaseMaps[[basemap]],K_range=c(3,18),h=h_i,k=7,p=0.3,h_base=0.8,plot_flag=FALSE)
      save(kmap_i,file=paste0("seascapes/kmaps/kmap_",ind,".RData"))
      BaseKMapIndex$BaseMap[ind] <- basemap
      BaseKMapIndex$h[ind] <- h_i
    }
  }
}
# also add one map with constant K for each basemap
for(basemap in 1:length(BaseMaps)){
  ind <- ind+1
  kmap_i <- f_GenerateMapWithK(base_map=BaseMaps[[basemap]],K_range=c(10,10),h=0,k=7,p=0.3,h_base=0.8,plot_flag = TRUE)
  BaseKMapIndex$BaseMap[ind] <- basemap
  BaseKMapIndex$h[ind] <- "constant"
}

## All runs
## do 3 reps per Kmap
sim_index <- do.call(rbind,replicate(n_rep_per_k,BaseKMapIndex,simplify=FALSE))
sim_index$Run <- 1:nrow(sim_index)

save(sim_index,file="experiments/Array1_SimIndex.RData")

# test different size maps for memory and speed
# ----------------------------------------------------

# source('0_Setup.R')
# for(k_i in 5:7){
#   kmap_i <- f_GenerateMapWithK(base_map=NULL,K_range=c(3,18),h=0.9,k=k_i,p=0.3,h_base=0.9,plot_flag=TRUE)
#   save(kmap_i,file=paste0("seascapes/sizetest/kmap_",k_i,".RData"))  
# }
# 
# load(paste0("seascapes/sizetest/kmap_",9,".RData")) #78000 patches
# load(paste0("seascapes/sizetest/kmap_",8,".RData")) #20000 patches
# load(paste0("seascapes/sizetest/kmap_",7,".RData")) #5000 patches
# load(paste0("seascapes/sizetest/kmap_",6,".RData")) #1200 patches
# load(paste0("seascapes/sizetest/kmap_",5,".RData")) #325 patches
