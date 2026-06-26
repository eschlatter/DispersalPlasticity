source('0_Setup.R')
library(calculus)
experiment_folder <- "experiments/SmallSim_Testing"
hab_id=5
param_id=1

load(file.path(experiment_folder,paste0("params_",param_id,".RData"))) # params
list2env(x=params,envir=environment())
load(paste0(experiment_folder,"/habparams_",hab_id,".RData")) # hab params
list2env(x=hab_params,envir=environment())

# patch distances and nav_rad, in km
patch_dists <- st_distance(sfc_patches)
#patch_angles <- suppressWarnings(2*asin(nav_rad/patch_dists)/(2*pi))
patch_angles <- (2*nav_rad)/(2*pi*patch_dists)
patch_angles[is.nan(patch_angles)] <- 1
units(patch_angles) <- NULL
patch_angles[patch_angles==Inf] <- 1
units(patch_dists) <- 'km'
units(patch_dists) <- NULL
units(nav_rad) <- NULL

group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)
map_ext <- ext(reef_sf)
map_bounds <- list(x=c(as.numeric(map_ext[1])/1000,as.numeric(map_ext[2])/1000),
                   y=c(as.numeric(map_ext[3])/1000,as.numeric(map_ext[4])/1000))

# info on a central point, to use for off-map correction and self-recruitment calculations
central_pt <- st_coordinates(st_centroid(reef_sf))
centr_x <- central_pt[1]/1000
centr_y <- central_pt[2]/1000
square_bounds <- list(x=c(centr_x-nav_rad,centr_x+nav_rad),y=c(centr_y-nav_rad,centr_y+nav_rad))

# overlap discount (multiply destination patch by this factor)
overlap_discount <- 1/colSums(patch_dists<nav_rad)

# index of alpha/theta combinations
kern_param_index <- matrix(1:(length(v_alphas)*length(v_thetas)),ncol=length(v_alphas),nrow=length(v_thetas))
kern_params <- data.frame(index=as.vector(kern_param_index),
                          alpha_ind=rep(1:length(v_alphas),each=length(v_thetas)),
                          theta_ind=rep(1:length(v_thetas),times=length(v_alphas)))
kern_params$alpha_val <- v_alphas[kern_params$alpha_ind]
kern_params$theta_val <- v_thetas[kern_params$theta_ind]

# make the off-the-map corrections and self-recruitment calculations
for(i in 1:nrow(kern_params)){
  al_i <- kern_params$alpha_val[i]
  th_i <- kern_params$theta_val[i]
  
  # off-the-map correction
  vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                        bounds=map_bounds,
                        params=list(theta=th_i),
                        relTol=0.00001
  )$value
  kern_params$offmap_corr[i] <- 1+(2*pi*th_i-vol_inmap)/vol_inmap
  
  # self-recruitment values
  vol_selfrecruit <- integral(function(r,theta){(1/th_i)*exp(-r/th_i)},
                              bounds=list(r=c(0,nav_rad),theta=c(0,2*pi)),
                              coordinates="polar",
                              relTol=0.00001
  )$value
  # vol_selfrecruit <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
  #                       bounds=square_bounds,
  #                       params=list(theta=th_i),
  #                       relTol=0.00001
  # )$value
  kern_params$selfrecruit[i] <- vol_selfrecruit
}


## Generate the conn_mats -------------
conn_mats <- list()
oldway <- list()
step0 <- list()
step1 <- list()
step2 <- list()
step3 <- list()

integral_ests <- data.frame(g=integer(),
                            oldway=numeric(),
                            connmat=numeric(),
                            step0=numeric(),
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
  
  # row of kernparam that corresponds to each origin patch
  kernparam_rows <- sapply(1:length(alpha_val_p),function(i) kern_param_index[plast_params$theta_plastic[i],plast_params$alpha_plastic[i]])
  
  
  # oldway[[g]] <- overlap_discount*patch_angles*
  #   (pgamma(patch_dists+nav_rad,shape=mat_alphas,scale=mat_thetas)-
  #      pgamma(pmax(patch_dists-nav_rad,0),shape=mat_alphas,scale=mat_thetas)+
  #      ifelse(nav_rad>patch_dists,pgamma(nav_rad-patch_dists,shape=mat_alphas,scale=mat_thetas),0))
  
  # approximate the area under the integral in a circle of radius nav_rad centered on each destination patch
  # conn_mat[i,j] is the density of propagules from origin site j landing in destination site i
  inprogress <- pi*(nav_rad^2)*dgamma(patch_dists,shape=mat_alphas,scale=mat_thetas)
  conn_mats[[g]] <- inprogress
  
  # apply corrections:
  # First first, fix the self-recruitment values
  selfrecruits <- kern_params$selfrecruit[kernparam_rows]
  diag(inprogress) <- selfrecruits
  #diag(inprogress) <- 0
  step0[[g]] <- inprogress
  # First, divide by 2pitheta so all integrals sum to one
  inprogress <- inprogress/(2*pi*mat_thetas)
  step1[[g]] <- inprogress
  # Second, apply a correction factor to account for the fact that broader kernels lose more offspring to edge effects
  step2mat <- matrix(kern_params$offmap_corr[kernparam_rows],nrow=npatch,ncol=npatch,byrow=TRUE)
  inprogress <- inprogress*step2mat
  step2[[g]] <- inprogress
  # Third, apply a correction for close neighbors
  inprogress <- inprogress*overlap_discount
  step3[[g]] <- inprogress
  conn_mats[[g]] <- inprogress
  integral_ests_i <- data.frame(g=g,
                                oldway=as.vector(oldway[[g]]),
                                connmat=as.vector(conn_mats[[g]]),
                                step0=as.vector(step0[[g]]),
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

df_all <- data.frame(g=integer(),originpatch=integer(),repr_out=numeric())
for(i in 1:nrow(group_index)){
  group_index$surv_total[i] <- sum(step3[[i]])
  group_index$surv_total_old[i] <- sum(oldway[[i]])
  temp_df <- data.frame(g=i,originpatch=1:npatch,repr_out=colSums(conn_mats[[i]]),repr_out_old=colSums(oldway[[i]]))
  df_all <- rbind(df_all,temp_df)
}
ggplot(group_index,aes(x=theta,color=p))+geom_point(aes(y=surv_total,color="New Way"))+geom_point(aes(y=surv_total_old,color="Old Way"))


df_all <- left_join(df_all,mutate(group_index,g=row_number()),by=join_by(g))
df_all$alpha_val <- v_alphas[df_all$alpha]
df_all$theta_val <- v_thetas[df_all$theta]
df_all <- left_join(df_all,patch_locations,by=join_by(originpatch==id))
df_all$plastic_theta <- v_thetas[f_plasticityb(df_all$b,v_p[df_all$p],df_all$alpha,df_all$theta,
                                               n_alpha=length(v_alphas),n_theta=length(v_thetas))$theta_plastic]

#load(file="2kpop_noself.RData")

# effective theta value of origin patch vs reproductive output
g1 <- df_all |> group_by(g,originpatch) |> 
  summarize(reprod_out=sum(repr_out),theta_val=first(plastic_theta),.groups="drop") |> 
  ggplot(aes(x=theta_val,y=reprod_out))+geom_point()+geom_smooth(method="loess",se=FALSE)+
  labs(title="New way",x="Effective theta",y="Reproductive output")+
  lims(y=c(0,0.55))

g2 <- df_all |> group_by(g,originpatch) |> 
  summarize(reprod_out=sum(repr_out),reprod_out_old=sum(repr_out_old),theta_val=first(plastic_theta),.groups="drop") |> 
  ggplot(aes(x=theta_val,y=reprod_out_old))+geom_point()+geom_smooth(method="loess",se=FALSE)+
  labs(title="Old way",x="Effective theta",y="Reproductive output")+
  lims(y=c(0,0.55))

grid.arrange(g2,g1,nrow=1)

# fundamental theta value of origin patch vs reproductive output
df_all |> group_by(g,originpatch) |>
  summarize(reprod_out=sum(repr_out),theta_val=first(theta_val),.groups="drop") |> 
  ggplot(aes(x=theta_val,y=reprod_out))+geom_point()+geom_smooth(method="loess",se=FALSE)+
  labs(title="Fundamental theta of origin patch vs reproductive output")


#save(df_all,file="2kpop_withself.RData")







p2 <- ggplot(integral_ests,aes(x=theta_val,y=step3,group=interaction(g,originpatch),color=factor(g)))+geom_point(stat="summary",fun="sum")

grid.arrange(p1,p2,nrow=2)

integral_ests |> 
  group_by(g,originpatch) |> 
  summarize(reprod_out=sum(step3),theta_val=first(theta_val),.groups="drop") |> 
  mutate(th_val=v_thetas[group_index$theta[g]]) |>
  ggplot(aes(x=th_val,y=reprod_out))+geom_point()+geom_smooth(method="loess")

integral_ests |> group_by(g,originpatch) |> summarize(reprod_out=sum(step3),theta_val=first(theta_val),.groups="drop") |>
  ggplot(aes(x=theta_val,y=reprod_out))+geom_point()+geom_smooth(method="loess")


ggplot(pivot_longer(integral_ests,cols=c(connmat,step0,step1,step2,step3)),
       aes(x=theta_val,y=value,group=interaction(g,originpatch),color=factor(g)))+
  geom_point(stat="summary",fun="sum")+
  facet_wrap(~name,scales = "free")


g=20
patch_locations$a20_export <- colSums(step3[[g]])
patch_locations$a20_import <- rowSums(step3[[g]])
ggplot(patch_locations,aes(x=x,y=y,color=a20_import))+geom_point()





# Same thing, but a bigger dataset ---------
source('0_Setup.R')
library(calculus)
experiment_folder <- "experiments/SmallSim_Testing"
hab_id=1
param_id=1

load(file.path(experiment_folder,paste0("params_",param_id,".RData"))) # params
list2env(x=params,envir=environment())
load(paste0(experiment_folder,"/habparams_",hab_id,".RData")) # hab params
list2env(x=hab_params,envir=environment())

q_rast <- rast(paste0(experiment_folder,"/b1/set1/qmap_b1_q",hab_id,".tif"))
ggplot()+
  ggspatial::layer_spatial(q_rast$q)+
  scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
  geom_sf(data=sfc_patches)+
  annotation_scale()

# patch distances and nav_rad, in km
patch_dists <- st_distance(sfc_patches)
units(patch_dists) <- 'km'
units(patch_dists) <- NULL
units(nav_rad) <- NULL

group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)
map_ext <- ext(reef_sf)
map_bounds <- list(x=c(as.numeric(map_ext[1])/1000,as.numeric(map_ext[2])/1000),
                   y=c(as.numeric(map_ext[3])/1000,as.numeric(map_ext[4])/1000))

# info on a central point, to use for off-map correction and self-recruitment calculations
central_pt <- st_coordinates(st_centroid(reef_sf))
centr_x <- central_pt[1]/1000
centr_y <- central_pt[2]/1000
square_bounds <- list(x=c(centr_x-nav_rad,centr_x+nav_rad),y=c(centr_y-nav_rad,centr_y+nav_rad))

# overlap discount (multiply destination patch by this factor)
overlap_discount <- 1/colSums(patch_dists<nav_rad)

# index of alpha/theta combinations
kern_param_index <- matrix(1:(length(v_alphas)*length(v_thetas)),ncol=length(v_alphas),nrow=length(v_thetas))
kern_params <- data.frame(index=as.vector(kern_param_index),
                          alpha_ind=rep(1:length(v_alphas),each=length(v_thetas)),
                          theta_ind=rep(1:length(v_thetas),times=length(v_alphas)))
kern_params$alpha_val <- v_alphas[kern_params$alpha_ind]
kern_params$theta_val <- v_thetas[kern_params$theta_ind]

# make the off-the-map corrections and self-recruitment calculations
for(i in 1:nrow(kern_params)){
  al_i <- kern_params$alpha_val[i]
  th_i <- kern_params$theta_val[i]
  
  # off-the-map correction
  vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                        bounds=map_bounds,
                        params=list(theta=th_i),
                        relTol=0.00001
  )$value
  kern_params$offmap_corr[i] <- 1+(2*pi*th_i-vol_inmap)/vol_inmap
  
  # self-recruitment values
  vol_selfrecruit <- integral(function(r,theta){(1/th_i)*exp(-r/th_i)},
                              bounds=list(r=c(0,nav_rad),theta=c(0,2*pi)),
                              coordinates="polar",
                              relTol=0.00001
  )$value
  # vol_selfrecruit <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
  #                       bounds=square_bounds,
  #                       params=list(theta=th_i),
  #                       relTol=0.00001
  # )$value
  kern_params$selfrecruit[i] <- vol_selfrecruit
}


## Generate the conn_mats
conn_mats <- list()
for(g in 1:ngroups){
  # plastic values of alpha and theta
  grp_data <- group_index[g,]
  plast_params <- f_plasticityb(patch_locations$b,v_p[grp_data$p],grp_data$alpha,grp_data$theta,n_alpha=length(v_alphas),n_theta=length(v_thetas))
  alpha_val_p <- v_alphas[plast_params$alpha_plastic]
  theta_val_p <- v_thetas[plast_params$theta_plastic]
  mat_alphas <- matrix(rep(alpha_val_p,each=npatch),nrow=npatch)
  mat_thetas <- matrix(rep(theta_val_p,each=npatch),nrow=npatch)
  
  # row of kernparam that corresponds to each origin patch
  kernparam_rows <- sapply(1:length(alpha_val_p),function(i) kern_param_index[plast_params$theta_plastic[i],plast_params$alpha_plastic[i]])
  
  # approximate the area under the integral in a circle of radius nav_rad centered on each destination patch
  # conn_mat[i,j] is the density of propagules from origin site j landing in destination site i
  inprogress <- pi*(nav_rad^2)*dgamma(patch_dists,shape=mat_alphas,scale=mat_thetas)
  
  # apply corrections:
  # First first, fix the self-recruitment values
  selfrecruits <- kern_params$selfrecruit[kernparam_rows]
  diag(inprogress) <- selfrecruits
# diag(inprogress) <- 0
  # First, divide by 2pitheta so all integrals sum to one
  inprogress <- inprogress/(2*pi*mat_thetas)
  # Second, apply a correction factor to account for the fact that broader kernels lose more offspring to edge effects
  step2mat <- matrix(kern_params$offmap_corr[kernparam_rows],nrow=npatch,ncol=npatch,byrow=TRUE)
  inprogress <- inprogress*step2mat
  # Third, apply a correction for close neighbors
  inprogress <- inprogress*overlap_discount
  conn_mats[[g]] <- inprogress
}

df_all <- data.frame(g=integer(),originpatch=integer(),repr_out=numeric())
for(i in 1:nrow(group_index)){
  group_index$surv_total[i] <- sum(conn_mats[[i]])
  temp_df <- data.frame(g=i,originpatch=1:npatch,repr_out=colSums(conn_mats[[i]]))
  df_all <- rbind(df_all,temp_df)
}

ggplot(group_index,aes(x=theta,y=surv_total,color=p))+geom_point()


df_all <- left_join(df_all,mutate(group_index,g=row_number()),by=join_by(g))
df_all$alpha_val <- v_alphas[df_all$alpha]
df_all$theta_val <- v_thetas[df_all$theta]
df_all <- left_join(df_all,patch_locations,by=join_by(originpatch==id))
df_all$plastic_theta <- v_thetas[f_plasticityb(df_all$b,v_p[df_all$p],df_all$alpha,df_all$theta,
                                               n_alpha=length(v_alphas),n_theta=length(v_thetas))$theta_plastic]

#save(df_all,file="hab2_noself.RData")
load(file="2kpop_withself.RData")

# effective theta value of origin patch vs reproductive output
df_all |> group_by(g,originpatch) |> 
  summarize(reprod_out=sum(repr_out),theta_val=first(plastic_theta),p=v_p[first(p)],.groups="drop") |> 
  ggplot(aes(x=theta_val,y=reprod_out))+geom_point()+geom_smooth(method="loess",se=FALSE)+
  labs(title="Effective theta of origin patch vs reproductive output")

# fundamental theta value of origin patch vs reproductive output
df_all |> group_by(g,originpatch) |>
  summarize(reprod_out=sum(repr_out),theta_val=first(theta_val),.groups="drop") |> 
  ggplot(aes(x=theta_val,y=reprod_out))+geom_point()+geom_smooth(method="loess",se=FALSE)+
  labs(title="Fundamental theta of origin patch vs reproductive output")




# scratch ------------
# offmap_corr <- matrix(ncol=length(v_alphas),nrow=length(v_thetas))
# for(al_i in 1:ncol(offmap_corr)){
#   for(th_i in 1:nrow(offmap_corr)){
#     # calculate the proportion of the kernel volume that falls within the map for a kernel originating at central_pt
# vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
#                       bounds=map_bounds,
#                       params=list(theta=th_i),
#                       relTol=0.00001
# )$value
#     offmap_corr[th_i,al_i] <- 1+(2*pi*th_i-vol_inmap)/vol_inmap
#   }
# }
# for(al_i in 1:ncol(offmap_corr)){
#   for(th_i in 1:nrow(offmap_corr)){
#     v_inmapvols <- vector()
#     for(origin_pt in 1:npatch){
#       # calculate the proportion of the kernel volume that falls within the map for a kernel originating at central_pt
#       coord_x <- patch_locations$x[origin_pt]/1000
#       coord_y <- patch_locations$y[origin_pt]/1000
#       v_inmapvols[origin_pt] <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((coord_x-x)^2+(coord_y-y)^2)/theta)},
#                             bounds=map_bounds,
#                             params=list(theta=th_i),
#                             relTol=0.00001
#       )$value
#     }
#     v_inmapvols <- median(v_inmapvols)
#     v_offmapcorr <- 1+(2*pi*th_i-v_inmapvols)/v_inmapvols
#     offmap_corr[th_i,al_i] <- v_offmapcorr
#   }
# }


# # self-recruitment values
# square_bounds <- list(x=c(centr_x-nav_rad,centr_x+nav_rad),y=c(centr_y-nav_rad,centr_y+nav_rad))
# for(i in 1:nrow(kern_params)){
#   th_i <- kern_params$theta_val[i]
#   al_i <- kern_params$alpha_val[i]
#   vol_selfrecruit <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
#                               bounds=square_bounds,
#                               params=list(theta=th_i),
#                               relTol=0.00001
#   )$value
#   kern_params$selfrecruit[i] <- vol_selfrecruit
# }
# 

