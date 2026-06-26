source('0_Setup.R')
library(calculus)
experiment_folder <- "experiments/SmallSim_Testing"
hab_id=5
param_id=1

load(file.path(experiment_folder,paste0("params_",param_id,".RData"))) # params
list2env(x=params,envir=environment())
load(paste0(experiment_folder,"/habparams_",hab_id,".RData")) # hab params
list2env(x=hab_params,envir=environment())

group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)
map_ext <- ext(reef_sf)
map_bounds <- list(x=c(as.numeric(map_ext[1])/1000,as.numeric(map_ext[2])/1000),
                   y=c(as.numeric(map_ext[3])/1000,as.numeric(map_ext[4])/1000))

# Divide the map up into squares of side length 2*nav_rad
# (numbering starts in lower left corner and goes rowwise)
grid_sf <- st_make_grid(reef_sf, cellsize = c(2*as.numeric(nav_rad)*1000, 2*as.numeric(nav_rad)*1000)) %>%
  st_sf() %>%
  mutate(grid_id = row_number()) # Unique ID for each cell

# Find which square each anemone is in
points_with_id <- st_join(st_sf(sfc_patches), grid_sf)

# Make overlap_discount: a vector with an entry for each square, 1/how many anems are in it
overlap_discount <- tabulate(points_with_id$grid_id,nbins=nrow(grid_sf))
occupied_squares <- which(overlap_discount>0)

# check:
ip <- 51
ggplot()+
  geom_sf(data=grid_sf)+
  geom_sf(data=grid_sf[occupied_squares,],color='blue',lwd=1.2)+
  geom_sf(data=points_with_id)+
  # geom_sf(data=grid_sf[which.max(overlap_discount),],color='red',lwd=1.2)+ # the square with the most anemones
  # geom_sf(data=filter(points_with_id,grid_id==which.max(overlap_discount)),color='red')+ # and the anemones it contains
  geom_sf(data=points_with_id[ip,],color='red') # and the anemones it contains


## For a single parameter set (no plasticity, small kernel),
# get the exact values of the integral under the kernel between each pair of anemones
#g=16 # narrow kernel
#g=20 # or wide kernel
g=35

kern_vols <- matrix(nrow=npatch,ncol=npatch) # hold volume under kernel
conn_mat_35 <- matrix(nrow=npatch,ncol=npatch)
vols_inmap <- vector()

grp_data <- group_index[g,]
# plastic values of alpha and theta
plast_params <- f_plasticityb(patch_locations$b,v_p[grp_data$p],grp_data$alpha,grp_data$theta,n_alpha=length(v_alphas),n_theta=length(v_thetas))
alpha_val_p <- v_alphas[plast_params$alpha_plastic]
theta_val_p <- v_thetas[plast_params$theta_plastic]

# For each origin anemone (and its habitat quality):
for(anem_i in 1:npatch){
  patch_data <- patch_locations[anem_i,]
  xpt <- patch_data$x/1000 # coordinates should be in km, to match the kernel parameters
  ypt <- patch_data$y/1000
  al_i <- alpha_val_p[anem_i]
  th_i <- theta_val_p[anem_i]
  
  sq_probs <- rep(NA,length=nrow(grid_sf)) # making it the length of the total number of squares, rather than occupied squares, so indexing is easier later
  # Calculate the centroid-to-area dispersal probability from the anemone to each occupied square, and multiply by overlap_discount for that square
  for(dest_square in occupied_squares){
    sq_i <- grid_sf[dest_square,]
    sq_bounds <- list(x=as.numeric(st_bbox(sq_i))[c(1,3)]/1000,y=as.numeric(st_bbox(sq_i))[c(2,4)]/1000) # coords should be in km, to match the kernel params
    kern_int <- integral(function(x,y){(1/th_i)*exp(-sqrt((xpt-x)^2+(ypt-y)^2)/th_i)},
                         bounds=sq_bounds, #,absTol=1e-12 # might want to mess with this
                         relTol=0.00001
    )$value
    sq_probs[dest_square] <- kern_int
    # Try to vectorize this -- the documentation says that can lead to significant performance improvements
    # (Not vectorizing now. It's not clear how/whether vectorization works over multiple sets of bounds, rather than multiple variable/param inputs.)
  } # dest_square
  #sq_probs <- sq_probs*(1/overlap_discount)
  # Distribute these values to the appropriate destination anemones to get a vector of dispersal probabilities to each anemone
  kern_vols[,anem_i] <- sq_probs[points_with_id$grid_id]
  
  # Get the off-the-map correction for this anemone/param combo:
  # Calculate the volume of the kernel within the map (vol_inmap).
  vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((xpt-x)^2+(ypt-y)^2)/theta)},
                        bounds=map_bounds,
                        params=list(theta=th_i),
                        relTol=0.00001
  )$value
  vols_inmap[anem_i] <- vol_inmap
  # Correction factor is 1+(2*pi*theta-vol_inmap)/vol_inmap.
  inmap_correction <- 1+(2*pi*th_i - vol_inmap)/vol_inmap
  # Multiply the dispersal probability vector by the off-the-map correction factor.
  dest_probs <- kern_vols[,anem_i]*inmap_correction
  # Divide by 2*pi*theta, because the total integral of the dispersal kernel should be 2pitheta, and we want to correct for that too
  # dest_probs <- dest_probs/(2*pi*th_i)
  # Fill this in as the row (column) of destination probabilities for that origin anemone
  conn_mat_35[,anem_i] <- dest_probs
} # anem_i (origin anemone)
save(conn_mat_35,kern_vols,vols_inmap,file="scratch_conn_mat_35.RData")
#load("scratch_conn_mat_35.RData")

#### the estimations
# patch distances, in km
patch_dists <- st_distance(sfc_patches)
units(patch_dists) <- 'km'
units(patch_dists) <- NULL
units(nav_rad) <- NULL

f_kernvolest <- function(i,method="a"){
  if(method=="a"){
    kernvol <- pgamma(patch_dists[i,]+nav_rad,shape=alpha_val_p[i],scale=theta_val_p[i])-
      pgamma(pmax(patch_dists[i,]-nav_rad,0),shape=alpha_val_p[i],scale=theta_val_p[i])+
      ifelse(nav_rad>patch_dists[i,],
             pgamma(nav_rad-patch_dists[i,],shape=alpha_val_p[i],scale=theta_val_p[i]),
             0) # when patch_dists<nav_rad, correct for it    
  }
  if(method=="b"){ # just use the kernel value at the destination point, rather than integrating over a range of +/-nav_rad
    kernvol <- dgamma(patch_dists[i,],shape=alpha_val_p[i],scale=theta_val_p[i])
  }
  if(method=="c"){ # method a, but without correcting for patch_dists<nav_rad.
    kernvol <- pgamma(patch_dists[i,]+nav_rad,shape=alpha_val_p[i],scale=theta_val_p[i])-
      pgamma(pmax(patch_dists[i,]-nav_rad,0),shape=alpha_val_p[i],scale=theta_val_p[i])
  }
  return(kernvol)
}
a_conn_mat <- 2*nav_rad*sapply(1:npatch,f_kernvolest,method="a")
b_conn_mat <- 4*nav_rad^2*sapply(1:npatch,f_kernvolest,method="b")
#c_conn_mat <- 2*nav_rad*sapply(1:npatch,f_kernvolest,method="c")

integral_ests <- data.frame(kernvol=as.vector(kern_vols),
                            aconnmat=as.vector(a_conn_mat),
                            bconnmat=as.vector(b_conn_mat),
                            patchdist=as.vector(patch_dists),
                            originpatch=rep(1:npatch,each=npatch),
                            destpatch=rep(1:npatch,times=npatch))
integral_ests$q_origin <- patch_locations$q[integral_ests$originpatch]
integral_ests$theta_origin <- theta_val_p[integral_ests$originpatch]
integral_ests$theta_ind_origin <- plast_params$theta_plastic[integral_ests$originpatch]

integral_ests <- mutate(integral_ests,
                        step1=bconnmat/(2*pi*theta_origin))

ggplot(integral_ests,aes(x=patchdist,y=step1*100,color=factor(theta_origin)))+
  geom_point(alpha=0.4)+
#  scale_color_viridis_b(na.value = "red")+
#  geom_abline(slope=1,intercept=0,color="red")+
  theme_minimal()#+
  # geom_function(fun=dgamma, args=list(shape=1,scale=2.25),color='black')+
  # geom_function(fun=dgamma, args=list(shape=1,scale=4.25),color='black')


group_by(integral_ests,theta_origin) |>
  summarize(avg_bconnmat=mean(bconnmat),avg_step1=mean(step1))


central_pt <- st_coordinates(st_centroid(reef_sf))
centr_x <- central_pt[1]/1000
centr_y <- central_pt[2]/1000
out_over_in <- vector()
for(i in 1:length(v_thetas)){
  th_i <- v_thetas[i]
  # calculate the proportion of the kernel volume that falls within the map for a kernel originating at central_pt
  vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                        bounds=map_bounds,
                        params=list(theta=th_i),
                        relTol=0.00001
  )$value
  out_over_in[i] <- (2*pi*th_i-vol_inmap)/vol_inmap
}


integral_ests$step2 <- integral_ests$step1*(1+out_over_in[integral_ests$theta_ind_origin])

a <- group_by(integral_ests,theta_origin) |>
  summarize(avg_bconnmat=mean(bconnmat),avg_step1=mean(step1),avg_step2=mean(step2))
ggplot(pivot_longer(a,cols=-theta_origin),aes(x=theta_origin,y=value,color=name))+geom_line()


ggplot(filter(integral_ests,originpatch!=51),aes(x=patchdist,y=step2,color=factor(theta_origin)))+
  geom_point(alpha=0.4)+
  #  scale_color_viridis_b(na.value = "red")+
  #  geom_abline(slope=1,intercept=0,color="red")+
  theme_minimal()#+

par(mfrow=c(1,2))
plot(kern_vols,a_conn_mat,main="method a")
lines(c(0,1),c(0,1),col='red')

plot(kern_vols,b_conn_mat,main="method b")
lines(c(0,1),c(0,1),col='red')

# plot(conn_mat_20,c_conn_mat,main="method c")
# lines(c(0,1),c(0,1),col='red')


