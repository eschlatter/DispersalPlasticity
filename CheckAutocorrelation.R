source('0_Setup.R')

######## Look at output from array job ##########
basemap_file="seascapes/2026_03_30/1x25km_res=10m"
df_all <- read.csv(paste0(basemap_file,"/set1/df_all.csv")) %>%
  filter(!is.na(range))

### plots
group_by(df_all,q_autocorr,model) %>%
  ggplot(aes(x=factor(q_autocorr)))+
  geom_boxplot(aes(y=range/1000,color=model))+
  lims(y=c(0,10))

group_by(df_all,q_autocorr) %>%
  ggplot(aes(x=factor(q_autocorr)))+
  geom_boxplot(aes(y=range/1000))+
  lims(y=c(0,5))

table(df_all$model,df_all$q_autocorr)

centers <- group_by(df_all,q_autocorr) |>
  summarize(l40=quantile(range,0.48),u60=quantile(range,0.52),median=median(range))

df_all <- left_join(df_all,centers,by=join_by(q_autocorr))
df_all <- mutate(df_all,in_middle=(range>l40 & range<u60))

group_by(df_all,q_autocorr) %>%
  filter(in_middle==TRUE) %>%
  ggplot(aes(x=factor(q_autocorr)))+
  geom_boxplot(aes(y=range/1000))+
  lims(y=c(0,5))

df_all_mid <- filter(df_all,in_middle==TRUE) %>%
  arrange(range)

par(mfrow=c(1,2))

keeps_5x5 <- c(48,44,36,34,26,25,20,13,6)
keeps_1x25 <- c()
for(i in 1:nrow(df_all_mid)){
  rep_i <- df_all_mid$rep_i[i]
  q_autocorr <- df_all_mid$q_autocorr[i]
  qmap_file <- paste0("set1/q",q_autocorr,"_sim",rep_i)
  q_rast <- rast(paste0(basemap_file,"/",qmap_file,".tif")) # load q_rast
  # g <- ggplot()+ggspatial::layer_spatial(q_rast$q)+
  #   scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
  #   annotation_scale(pad_x=unit(0,"cm"))+
  #   theme_minimal()+
  #   theme(legend.position = "none",axis.text=element_blank())+
  #   labs(title=paste0("h=",q_autocorr,", sim",rep_i,", range=",round(df_all_mid$range[i])))
  # print(g)
  plot(q_rast$q,main=paste0("range=",round(df_all_mid$range[i]/1000,1),"km"),
       axes=FALSE,legend=FALSE,col=colorspace::sequential_hcl(256,palette="BluGrn",rev=TRUE))
#  plot(q_rast$q,main=paste0("h=",q_autocorr,"_sim",rep_i,", range=",round(df_all_mid$range[i])),
#       legend=FALSE,axes=FALSE,col=colorspace::sequential_hcl(256,palette="Purples2",rev=TRUE))
  sbar(d=1000,labels=c(0,1),below="km",type="bar",divs=2)
}

rep_i=87
q_autocorr=0.35
qmap_file <- paste0("set1/q",q_autocorr,"_sim",rep_i)
q_rast <- rast(paste0(basemap_file,"/",qmap_file,".tif")) # load q_rast
plot(q_rast$q)

###### look at variation by target distribution

df_dists <- data.frame(rep_i=numeric(),q_autocorr=numeric(),target_dist=factor(),
                       range=numeric(),sill=numeric(),SSErr=numeric(),
                       model=factor())
write.table(df_dists, file = paste0(basemap_file,"/set2/df_dists.csv"), sep = ",", append = FALSE,
                         quote = FALSE, col.names = TRUE, row.names = FALSE)

basemap_file="seascapes/2026_03_30/1x25km_res=10m"
popmap_file="pop_density800"
load(paste0(basemap_file,"/",popmap_file,".RData"))
dists_mat <- drop_units(patch_dists)
spdf1 <- as_Spatial(sfc_patches)
rep_i=1
q_autocorr=0.35
qmap_file <- paste0("set1/q",q_autocorr,"_sim",rep_i)
q_rast <- rast(paste0(basemap_file,"/",qmap_file,".tif")) # load q_rast

for(target_dist in c("identity","A","B","C","D","E")){
  q_rast$newdist <- f_TransformDist(matrix(values(q_rast$q),nrow=500,byrow=TRUE),target_dist)
  plot(q_rast$newdist,main=target_dist,axes=FALSE,col=colorspace::sequential_hcl(256,palette="BluGrn",rev=TRUE))
  
  # spdf1$q <- terra::extract(q_rast$newdist,vect(sfc_patches),xy=TRUE,search_radius=500)$newdist
  # # empirical variogram
  # vgm1 <- variogram(q~1,data=spdf1,cressie=TRUE,width=100,cutoff=10000) # choose a relevant bin width and cutoff here
  # # run exponential, gaussian and spherical models, and pick the best one (lowest std error)
  # vgmf <- fit.variogram(vgm1,vgm(c("Gau","Sph","Exp")))
  # 
  # g <- plot(vgm1,vgmf,main=paste0(target_dist))
  # print(g)
  # 
  # df_all <- data.frame(rep_i=rep_i, q_autocorr=q_autocorr,target_dist=target_dist,
  #                      range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
  #                      model=vgmf$model[2])
  # write.table(df_all, file = paste0(basemap_file,"/set2/df_dists.csv"), sep = ",", append = TRUE, 
  #             quote = FALSE, col.names = FALSE, row.names = FALSE)
}

df_dists <- read.csv(paste0(basemap_file,"/set2/df_dists.csv"))
########

## load data
basemap_file="seascapes/2026_03_27/10x10km_res=2m"
popmap_file="pop_1000"
load(paste0(basemap_file,"/",popmap_file,".RData"))
dists_mat <- drop_units(patch_dists)
qmap_file <- "set1/q=0_sim3"
q_rast <- rast(paste0(basemap_file,"/",qmap_file,".tif")) # load q_rast

## params
list_dists <- c('identity','A','B','C','D','E')
q_vec <- seq(from=0,to=0.45,length.out=10)
loop_over <- q_vec
target_dist="identity"

## data structures
df_all <- data.frame(rep_i=numeric(),loop_i=numeric(),range_gau=numeric(),sill_gau=numeric(),SSErr_gau=numeric(),
                     range_sph=numeric(),sill_sph=numeric(),SSErr_sph=numeric())

### quick little function wrapper to make this work in parallel
f_mc_qmaps <- function(i){
  # generate a qmap
  qmap_file <- paste0("set1/q=",q_autocorr,"_sim",i) # save the SpatRaster for later
  q_rast <- f_GenerateHabQual(base_rast=basemap_file,q_autocorr=q_autocorr,target_dist=target_dist,
                              plot_flag=TRUE,qmap_file = qmap_file)$q_rast
  
  # get variogram info
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  # empirical variogram
  vgm1 <- variogram(q~1,data=spdf1,cressie=TRUE)
  # run gaussian and spherical models, and pick the better one
  vgm_gau <- fit.variogram(vgm1,vgm(c("Gau")))
  vgm_sph <- fit.variogram(vgm1,vgm(c("Sph")))
  
  # store output
  df_all <- data.frame(rep_i=i, loop_i=loop_i,range_gau=vgm_gau$range[2],sill_gau=vgm_gau$psill[2],SSErr_gau=attr(vgm_gau,"SSErr"),
                       range_sph=vgm_sph$range[2],sill_sph=vgm_sph$psill[2],SSErr_sph=attr(vgm_sph,"SSErr"))
}

index_runs <- expand.grid(loop_i=loop_over,rep_i=1:100)
run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
loop_i <- index_runs$loop_i[run_i]
rep_i <- index_runs$rep_i[run_i]

print(paste0("h=",loop_i,", rep=",rep_i))


### generate all the qmaps and calculate variogram range
for(loop_i in loop_over){
  print(loop_i)
  q_autocorr=loop_i

  a <- lapply(1:20, f_mc_qmaps)
  df_loop <- do.call(rbind,a)
  df_all <- rbind(df_all,df_loop)
} # loop_i

save(df_all,file=paste0(basemap_file,"/set1/df_all.RData"))

df_all <- read.csv(paste0(basemap_file,"/set2/df_all.csv"))

### plots
group_by(df_all,q_autocorr,model) %>%
  ggplot(aes(x=factor(q_autocorr)))+
  geom_boxplot(aes(y=range/1000,color=model))+
  lims(y=c(0,5))

group_by(df_all,q_autocorr) %>%
  ggplot(aes(x=factor(q_autocorr)))+
  geom_boxplot(aes(y=range/1000))+
  lims(y=c(0,5))

table(df_all$model,df_all$q_autocorr)

qmap_file <- "set1/q0.45_sim15"
q_rast <- rast(paste0(basemap_file,"/",qmap_file,".tif")) # load q_rast

ggplot(df_all,aes(x=loop_i,y=range/1000,group=loop_i))+
  geom_boxplot()+
  labs(y="semivariogram range (km)",x="h")

ggplot(df_all,aes(x=loop_i,y=sill,group=loop_i))+
  geom_boxplot()+
  labs(y="semivariogram sill",x="fracland parameter h")



df_all <- df_all |>
  rownames_to_column(var="map_id") |>
  mutate(range_bin=cut(range,breaks=seq(from=0,to=30000,by=5000))) |>
  mutate(range_bin_lower=as.numeric( sub("\\((.+),.*", "\\1", range_bin) ),
         range_bin_upper=as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", range_bin) )) |>
  mutate(range_bin=cut(range,breaks=seq(from=0,to=30000,by=5000),labels=FALSE))
table(df_all$loop_i,df_all$range_bin)

rplo <- function(i){
  plot(rast_list[[i]]$q,main=paste(i,", h=",df_all$loop_i[i]))
}





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
