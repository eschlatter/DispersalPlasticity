source('0_Setup.R')
library(gstat)
library(spdep)

basemap_file="seascapes/2026_03_24/50x50km_res=1km"
popmap_file="seascapes/2026_03_24/popmap_50x50km"

# method="gaussian"
# v_qs <- 1:4
# noise=0.01

method="fractal"
v_qs <- seq(from=0,to=1.99,length.out=6)

load(paste0(popmap_file,".RData"))
dists_mat <- drop_units(patch_dists)

df_all <- data.frame(q=numeric(),range=numeric(),sill=numeric())
for(q_i in v_qs){
  print(q_i)
  df_q <- data.frame(st_coordinates(sfc_patches))
  npatch <- nrow(df_q)
  
  # create habitat quality layers
  if(method=="fractal"){
    qmap_file=NULL
    q_range=c(1,10)
    q_autocorr=q_i
    for(i in 1:100){
      q_rast <- f_GenerateHabQual(base_rast=basemap_file,q_range,q_autocorr,binary=FALSE,plot_flag=FALSE,qmap_file = qmap_file)$q_rast
      df_q <- cbind(df_q,terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q)
    }
  } else if(method=="gaussian"){
    for(i in 1:100){
      cov_Sigma <- exp(-q_i*dists_mat) # let covariance decay exponentially with distance
      q_vec <- mvrnorm(n=1, mu=rep(5,times=npatch), Sigma=cov_Sigma+noise*diag(npatch))
      df_q <- cbind(df_q,q_vec)
      
      ggplot(df_q,aes(x=X,y=Y))+geom_point(aes(color=q_vec))
    }
  }
  qnames <- paste0("q",1:100)
  names(df_q) <- c("X","Y",qnames)
  
  
  # fit sample variograms
  for(i in 1:(length(df_q)-2)){
    print(i)
    spdf1 <- as_Spatial(sfc_patches)
    spdf1$q <- df_q[,i+2]
    vgm1 <- variogram(q~1,data=spdf1,cutoff=50000)
    vgmf <- fit.variogram(vgm1,vgm("Gau"))
    # g <- plot(vgm1,vgmf,main=paste(i))
    # ggplot(df_q,aes(x=X,y=Y))+geom_point(aes(color=q2))
    # print(g)
    df_all <- rbind(df_all,data.frame(q=q_i,range=vgmf$range[2],sill=vgmf$psill[2]))
  }
}


#save(df_all,file="experiments/autocorr/Experiment5.RData")

ggplot(df_all,aes(x=q,y=range/1000,group=q))+
  geom_boxplot()+
  labs(y="semivariogram range (km)",x="fracland parameter h")
+
  lims(y=c(0,4))
ggplot(df_all,aes(x=log(q),y=sill,group=q))+
  geom_boxplot()+
  labs(y="semivariogram sill",x="fracland parameter h")+
  lims(y=c(0,4))








######## Incremental Moran's I plot
# calculate incremental spatial autocorrelation (like Robin's paper)
f_incr_Moran <- function(dists_mat,v_q,max_bin){
  dists_mat_NA <- dists_mat
  diag(dists_mat_NA) <- NA
  nearest_dists <- do.call(pmin,c(as.data.frame(dists_mat_NA),na.rm=TRUE))
  #  min_band <- min(nearest_dists) # minimum distance so everybody has at least one neighbor
  min_band <- median(nearest_dists)
  band_incr <- mean(nearest_dists) # average distance to each feature's nearest neighboring feature
  farthest_dists <- do.call(pmax,c(as.data.frame(dists_mat_NA),na.rm=TRUE)) 
  #max_band <- min(farthest_dists) # max dist so everybody still has at least one neighbor
  max_band <- max_bin
  
  incr_moran <- data.frame(band=seq(from=min_band,by=band_incr,to=max_band),z=NA, MI=NA, p=NA)
  
  for(i in 1:nrow(incr_moran)){
    dist_band <- incr_moran$band[i]
    
    # set to 1 for all points within the range (dist_band,dist_band+band_incr), and 0 for all points outside
    s.dist <- dists_mat
    s.dist[dists_mat<dist_band] <- 1
    s.dist[dists_mat>dist_band] <- 0
    # s.dist[s.dist>(dist_band+band_incr)] <- 0
    # s.dist[s.dist>0] <- 1
    
    # spatial weights matrix (normalized version of s.dist)
    w.dist <- s.dist/colSums(s.dist)
    nonzero_neigh <- which(colSums(w.dist,na.rm=TRUE)!=0)
    w.dist <- w.dist[nonzero_neigh,][,nonzero_neigh] # get rid of sites without neighbors in this band (is this ok?)
    
    # store it
    if(length(nonzero_neigh)>0){
      MC <- Moran.I(x=v_q[nonzero_neigh],weight=w.dist)
      z <- (MC$observed-MC$expected)/MC$sd # calculate the z-score
      incr_moran$z[i] <- z
      incr_moran$MI[i] <- MC$observed
      incr_moran$p[i] <- 2*pnorm(abs(z),lower.tail = FALSE)
    } else incr_moran$MC[i] <- NA
  }  
  
  g <- ggplot(incr_moran,aes(x=band,y=MI))+
    geom_line(alpha=0.3)+
    geom_point(aes(color=(p<0.05)),size=0.5)+
    scale_color_manual(values=c("black","grey"),breaks=c(TRUE, FALSE))+
    #labs(x="distance band (km)",y="Moran's I",title=paste0("Incremental Spatial Autocorrelation, phi=",phi))+
    geom_hline(yintercept=0,alpha=0.3,lty="dashed")+
    theme_minimal()
  print(g)
  
  return(incr_moran)  
}

sp_aut_dist <- incr_moran$band[which.max(incr_moran$z)]

# # 4. Make objects for use below
# q_rast <- rast(paste0(qmap_file,".tif")) # load q_rast
# v_q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
max_bin=3

ggplot()+
  ggspatial::layer_spatial(q_rast$q)+
  scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
  geom_sf(data=sfc_patches)+
  annotation_scale()

moran_out <- f_incr_Moran(dists_mat,v_q,max_bin)

####### scratch
d <- seq(from=0,to=5,by=0.2)
H=0.5
C=3
sv <- (C/2)*d^(2*H)
lines(d,sv)
