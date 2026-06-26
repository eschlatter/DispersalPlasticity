# Same thing, but a bigger dataset ---------
source('0_Setup.R')
library(calculus)
experiment_folder <- "experiments/SmallSim_Testing"
hab_id=5
param_id=1
saveto_folder=paste0("/scratch.global/schla103/",experiment_folder,"/hab",hab_id)

load(file.path(experiment_folder,paste0("params_",param_id,".RData"))) # params
list2env(x=params,envir=environment())
#load(paste0(experiment_folder,"/habfiles/habparams_",hab_id,".RData")) # hab params
load(paste0(experiment_folder,"/habparams_",hab_id,".RData")) # hab params
list2env(x=hab_params,envir=environment())

# q_rast <- rast(paste0(experiment_folder,"/b1/set1/qmap_b1_q",hab_id,".tif"))
# ggplot()+
#   ggspatial::layer_spatial(q_rast$q)+
#   scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
#   geom_sf(data=sfc_patches)+
#   annotation_scale()

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
patch_locations$overlap_discount <- 1/colSums(patch_dists<nav_rad)

# index of alpha/theta combinations
eff_kern_index <- matrix(1:(length(v_alphas)*length(v_thetas)),ncol=length(v_alphas),nrow=length(v_thetas))
eff_kern_df <- data.frame(index=as.vector(eff_kern_index),
                          alpha_ind=rep(1:length(v_alphas),each=length(v_thetas)),
                          theta_ind=rep(1:length(v_thetas),times=length(v_alphas)))
eff_kern_df$alpha_val <- v_alphas[eff_kern_df$alpha_ind]
eff_kern_df$theta_val <- v_thetas[eff_kern_df$theta_ind]

# make the off-the-map corrections and self-recruitment calculations
for(i in 1:nrow(eff_kern_df)){
  al_i <- eff_kern_df$alpha_val[i]
  th_i <- eff_kern_df$theta_val[i]
  
  # off-the-map correction
  vol_inmap <- integral(function(x,y,theta){(1/theta)*exp(-sqrt((centr_x-x)^2+(centr_y-y)^2)/theta)},
                        bounds=map_bounds,
                        params=list(theta=th_i),
                        relTol=0.00001
  )$value
  eff_kern_df$offmap_corr[i] <- 1+(2*pi*th_i-vol_inmap)/vol_inmap
  
  # self-recruitment values
  vol_selfrecruit <- integral(function(r,theta){(1/th_i)*exp(-r/th_i)},
                              bounds=list(r=c(0,nav_rad),theta=c(0,2*pi)),
                              coordinates="polar",
                              relTol=0.00001
  )$value
  eff_kern_df$selfrecruit[i] <- vol_selfrecruit
}

dir.create(saveto_folder,recursive=TRUE)
save(group_index,patch_locations,v_p,v_alphas,v_thetas,eff_kern_index,eff_kern_df,patch_dists,nav_rad,
     file=paste0(saveto_folder,"/patch_dists_etc.RData"))

load(paste0(saveto_folder,"/patch_dists_etc.RData"))

# #### function for connmats
# # output: rates of dispersal from each patch to each other patch
# # Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i
# f_GetConnMat <- function(g,group_index,patch_locations,v_p,v_alphas,v_thetas,eff_kern_index,eff_kern_df,patch_dists,nav_rad,numCores){
#   v <- group_index[g,]
#   # compute effective parameters for each patch with plasticity (once per group)
#   eff_params <- f_plasticityb(patch_locations$b, 
#                               v_p[v$p], 
#                               v$alpha, 
#                               v$theta,
#                               n_alpha = length(v_alphas),
#                               n_theta = length(v_thetas))
#   alpha <- eff_params$alpha_plastic
#   theta <- eff_params$theta_plastic
#   npatch=length(alpha)
#   # row of kernparam that corresponds to each origin patch
#   kernparam_rows <- sapply(1:npatch,function(i) eff_kern_index[eff_params$theta_plastic[i],eff_params$alpha_plastic[i]])
#   
#   connectivity_matrix <- mclapply(1:npatch,function(i){ # do once for each origin patch
#     # Evaluate volume under kernel in a circle of radius nav_rad around destination site
#     cm_i <- pi*nav_rad^2*dgamma(patch_dists[,i],shape=alpha[i],scale=theta[i])
#     # Improve volume estimate for self-recruitment value
#     cm_i[i] <- eff_kern_df$selfrecruit[kernparam_rows[i]]
#     # Divide by 2pitheta so all integrals sum to one
#     cm_i <- cm_i/(2*pi*theta[i])
#     # Apply a correction factor to account for the fact that broader kernels lose more offspring to edge effects
#     cm_i <- cm_i*eff_kern_df$offmap_corr[kernparam_rows[i]]
#     # Apply a correction for close neighbors
#     cm_i <- cm_i*patch_locations$overlap_discount
#   }, mc.cores=numCores)
#   connectivity_matrix <- do.call(cbind,connectivity_matrix)
#   return(connectivity_matrix)
# }
#   

g=1
numCores=1

a <- f_GetConnMat(g,
                  group_index,patch_locations,v_p,v_alphas,v_thetas,eff_kern_index,eff_kern_df,patch_dists,nav_rad,
                  numCores=parallelly::availableCores())


######### Make and save connectivity matrices ##########

for(g in 1:nrow(group_index)){
  # get the connectivity matrix among patches, given the group parameter values and patch-level q's
  # (and accounting for the patch population x per capita output b_i from each patch)
  
  # calculate this matrix
  ## This is what uses multiple cores: f_GetPlasticConnMat uses f_GetConnectivityMatrix_parallel
  conn_mat <- f_GetConnMat(g,
                           group_index,patch_locations,v_p,v_alphas,v_thetas,eff_kern_index,eff_kern_df,patch_dists,nav_rad,
                           numCores=parallelly::availableCores())
  # then store it
  #write_fst(as.data.frame(conn_mat),paste0(temp_path,"/map_",mapID_i,"/grp_",g),compress=0)
  saveRDS(object=conn_mat,file=paste0(saveto_folder,"/grp_",g,".rds"),compress=FALSE)
  
} # g





# ## Looking at it afterwards --------
# df_all <- data.frame(g=integer(),originpatch=integer(),repr_out=numeric())
# for(i in 1:nrow(group_index)){
#   group_index$surv_total[i] <- sum(conn_mats[[i]])
#   temp_df <- data.frame(g=i,originpatch=1:npatch,repr_out=colSums(conn_mats[[i]]))
#   df_all <- rbind(df_all,temp_df)
# }
# 
# ggplot(group_index,aes(x=theta,y=surv_total,color=p))+geom_point()
# 
# 
# df_all <- left_join(df_all,mutate(group_index,g=row_number()),by=join_by(g))
# df_all$alpha_val <- v_alphas[df_all$alpha]
# df_all$theta_val <- v_thetas[df_all$theta]
# df_all <- left_join(df_all,patch_locations,by=join_by(originpatch==id))
# df_all$plastic_theta <- v_thetas[f_plasticityb(df_all$b,v_p[df_all$p],df_all$alpha,df_all$theta,
#                                                n_alpha=length(v_alphas),n_theta=length(v_thetas))$theta_plastic]
# 
# #save(df_all,file="hab2_noself.RData")
# load(file="2kpop_withself.RData")
# 
# # effective theta value of origin patch vs reproductive output
# df_all |> group_by(g,originpatch) |> 
#   summarize(reprod_out=sum(repr_out),theta_val=first(plastic_theta),p=v_p[first(p)],.groups="drop") |> 
#   ggplot(aes(x=theta_val,y=reprod_out))+geom_point()+geom_smooth(method="loess",se=FALSE)+
#   labs(title="Effective theta of origin patch vs reproductive output")
# 
# # fundamental theta value of origin patch vs reproductive output
# df_all |> group_by(g,originpatch) |>
#   summarize(reprod_out=sum(repr_out),theta_val=first(theta_val),.groups="drop") |> 
#   ggplot(aes(x=theta_val,y=reprod_out))+geom_point()+geom_smooth(method="loess",se=FALSE)+
#   labs(title="Fundamental theta of origin patch vs reproductive output")
