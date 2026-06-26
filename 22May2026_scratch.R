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

# patch distances and nav_rad, in km
patch_dists <- st_distance(sfc_patches)
units(patch_dists) <- 'km'
units(patch_dists) <- NULL
units(nav_rad) <- NULL

# overlap discount (multiply destination patch by this factor)
overlap_discount <- 1/colSums(patch_dists<nav_rad)

# off-the-map correction
central_pt <- st_coordinates(st_centroid(reef_sf))
centr_x <- central_pt[1]/1000
centr_y <- central_pt[2]/1000
offmap_corr <- matrix(ncol=length(v_alphas),nrow=length(v_thetas))
for(al_i in 1:ncol(offmap_corr)){
  for(th_i in 1:nrow(offmap_corr)){
    # calculate the proportion of the kernel volume that falls within the map for a kernel originating at central_pt
    vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                          bounds=map_bounds,
                          params=list(theta=v_thetas[th_i]), ###### BUG! I was using a theta index here before.
                          relTol=0.00001
    )$value
    offmap_corr[th_i,al_i] <- 1+(2*pi*v_thetas[th_i]-vol_inmap)/vol_inmap ###### BUG! Was also using theta index here.
  }
}

conn_mats <- list()
step1 <- list()
step2 <- list()
step3 <- list()

integral_ests <- data.frame(g=integer(),
                              connmat=numeric(),
                              step1=numeric(),
                              step2=numeric(),
                              step3=numeric(),
                              alpha_val=numeric(),
                              theta_val=numeric(),
                            patchdist=numeric(),
                            originpatch=integer(),
                            destpatch=integer())
for(g in 1:ngroups){
  # plastic values of alpha and theta
  grp_data <- group_index[g,]
  plast_params <- f_plasticityb(patch_locations$b,v_p[grp_data$p],grp_data$alpha,grp_data$theta,n_alpha=length(v_alphas),n_theta=length(v_thetas))
  alpha_val_p <- v_alphas[plast_params$alpha_plastic]
  theta_val_p <- v_thetas[plast_params$theta_plastic]
  mat_alphas <- matrix(rep(alpha_val_p,each=npatch),nrow=npatch)
  mat_thetas <- matrix(rep(theta_val_p,each=npatch),nrow=npatch)
  
  # approximate the area under the integral in a square of side length 2*nav_rad centered on each destination patch
  # conn_mat[i,j] is the density of propagules from origin site j landing in destination site i
  conn_mats[[g]] <- 4*(nav_rad^2)*dgamma(patch_dists,shape=mat_alphas,scale=mat_thetas)
  
  # apply corrections:
  # First, divide by 2pitheta so all integrals sum to one
  step1[[g]] <- conn_mats[[g]]/(2*pi*mat_thetas)
  # Second, apply a correction factor to account for the fact that broader kernels lose more offspring to edge effects
  step2mat <- matrix(sapply(1:npatch,function(i) offmap_corr[plast_params$theta_plastic[i],plast_params$alpha_plastic[i]]),
                     nrow=npatch,ncol=npatch,byrow=TRUE)
  step2[[g]] <- step1[[g]]*step2mat
  # Third, apply a correction for close neighbors
  step3[[g]] <- step2[[g]]*overlap_discount
  
  integral_ests_i <- data.frame(g=g,
                                   connmat=as.vector(conn_mats[[g]]),
                                   step1=as.vector(step1[[g]]),
                                   step2=as.vector(step2[[g]]),
                                   step3=as.vector(step3[[g]]),
                                   alpha_val=as.vector(mat_alphas),
                                   theta_val=as.vector(mat_thetas),
                                patchdist=as.vector(patch_dists),
                                originpatch=rep(1:npatch,each=npatch),
                                destpatch=rep(1:npatch,times=npatch))
  integral_ests <- rbind(integral_ests,integral_ests_i)
}

ggplot(integral_ests,aes(x=theta_val,y=step3,group=interaction(g,originpatch),color=factor(g)))+geom_point(stat="summary",fun="sum")


ggplot(pivot_longer(integral_ests,cols=c(connmat,step1,step2,step3)),
       aes(x=theta_val,y=value,group=interaction(g,originpatch),color=factor(g)))+
  geom_point(stat="summary",fun="sum")+
  facet_wrap(~name,scales = "free")



ggplot(filter(integral_ests,originpatch!=destpatch),aes(x=theta_val,y=step3,group=interaction(g,originpatch),color=factor(g)))+geom_point(stat="summary",fun="sum")
