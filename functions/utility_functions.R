library(gridExtra)

# generate a simulated basemap
# inputs:
#   x_dist,y_dist: size of map in the x and y directions (in meters)
#   resol: vector of resolution (in degrees) in the x and y directions
#   method: "fractal" (fractal landscape); "uniform" (all sites habitable)
#   h: habitat aggregation value (between -2 and 2; higher = more autocorrelated)
#   prop_hab: proportion of the map that should be habitat
# output:
#   reef_sf: sfc_multipolygon of the reef area
#   bathy_rast: marmap::bathy object, for getting in-water distances
#   base_rast: SpatRaster with 0 for open water, 1 for reef
f_GenerateBasemap <- function(x_dist=500,y_dist=500,resol=c(0.00005,0.00005),method="fractal",h=NA,prop_hab=NA,make_dist_mat=TRUE){
  # create empty raster
  endpt_lat <- as.numeric(geosphere::destPoint(c(0,0),b=0,d=y_dist)[,'lat'])
  endpt_lon <- as.numeric(geosphere::destPoint(c(0,0),b=90,d=x_dist)[,'lon'])
  base_rast <- rast(xmin=0,xmax=endpt_lon,ymin=0,ymax=endpt_lat,resolution=c(0.00005,0.00005))
  crs(base_rast) <- "epsg:4326"
  nx=ncol(base_rast)
  ny=nrow(base_rast)
  
  # generate habitat configuration and add it to base_rast
  if(method=="fractal"){
    #fractal landscape
    dimens <- (2^(1:15)+1)
    k <- first(which(dimens>=max(nx,ny)))
    base_map=fracland(k=k,h=h,p=1-prop_hab,binary=TRUE,plotflag=FALSE)
    base_map <- base_map[1:ny,1:nx]
  } else if(method=="uniform"){
    # full grid is habitable
    base_map <- matrix(1,nrow=ny,ncol=nx)
  } else stop("method incorrectly specified")
  values(base_rast) <- base_map   # add habitat configuration to raster
  
  # create bathy_raster (for getting in-water distances)
  bathy_rast <- marmap::as.bathy(raster(base_rast)) # this rotates it! Why???
  bathy_rast[bathy_rast==1] <- -20 # reef
  bathy_rast[bathy_rast==0] <- -100 # open water
  # bathy_rast[bathy_rast==-1] <- 100 # land 
  ## autoplot.bathy(bathy_rast,geom="raster")
  
  # create reef_sf (for simulating point locations and plotting)
  if(method=="uniform"){
    # generate reef_sf manually
    reef_sf <- st_multipolygon(x = list(list(rbind(c(0,0),c(endpt_lon,0),c(endpt_lon,endpt_lat),c(0,endpt_lat),c(0,0)))))
  } else{
    basemap_stars <- st_as_stars(base_rast[[1]])
    basemap_contour <- st_contour(basemap_stars,breaks=c(0.5))
    reef_sf <- basemap_contour[basemap_contour$Min>0,] # pick out just the reef part for the shapefile
  }
  
  # create distance matrix
  if(make_dist_mat==TRUE){
    # take base_rast and make an sf object with a point in the center of each reef cell
    yn_reef <- values(base_rast,dataframe=TRUE) # reef/water values
    coords <- crds(base_rast,df=TRUE) # coordinates
    patches <- cbind(coords,yn_reef) # put them together 
    patches <- subset(patches,lyr.1!=0) # remove open water patches
    st_patches <- st_multipoint(x=as.matrix(patches[,c('x','y')]))
    sfc_patches <- st_sfc(st_patches,crs=crs(base_rast))
    sfc_patches <- st_cast(sfc_patches,'POINT')
    
    # # also make a dataframe, so we can connect distances/angles to patches
    # df_patches <- data.frame(st_coordinates(sfc_patches))
    # npatch <- nrow(df_patches)
    # df_patches$id <- 1:npatch
    
    # make patch_dists
    patch_dists <- st_distance(sfc_patches)
    units(patch_dists) <- 'km'
  } else {
    patch_dists=NA
    sfc_patches=NA
  }
  
  return(list(base_rast=base_rast,bathy_rast=bathy_rast,reef_sf=reef_sf,
              patch_dists=patch_dists,sfc_patches=sfc_patches))
}

# Generate habitat quality map
# Inputs:
#   base_rast: SpatRaster with 0 for open water, 1 for reef
#   q_range: c(qmin,qmax), where q=habitat quality
#   q_autocorr: some measure of the spatial autocorrelation in q.
#               For the fractal landscape method, it's h in the fracland function (higher values = more autocorrelated, range=(-2,2)(technically (-infinity,2) but don't worry about it))
# Outputs:
#   q_rast: SpatRaster object with two layers: reef (0 for open water and 1 for reef) and q (habitat quality)
f_GenerateHabQual <- function(base_rast,q_range,q_autocorr,plot_flag=FALSE){
  # add habitat values on top of base map
  nx=ncol(base_rast)
  ny=nrow(base_rast)
  dimens <- (2^(1:15)+1)   # find the k value to use in fracland function, given the dimensions of the base map
  k <- first(which(dimens>=max(nx,ny)))
  # generate a fractal layer
  frac_map <- fracland(k=k,h=q_autocorr,binary=FALSE,plotflag=FALSE)
  # convert to desired range of q values
  frac_map <- (frac_map-min(frac_map))/(max(frac_map)-min(frac_map)) # first to 0-1
  frac_map <- frac_map*(q_range[2]-q_range[1])+q_range[1] # then to q_range
  frac_map <- frac_map[1:ny,1:nx] 
  empty_rast <- rast(ext(base_rast), resolution=res(base_rast), crs = crs(base_rast))
  values(empty_rast) <- frac_map
  q_rast <- c(base_rast,empty_rast*base_rast)
  names(q_rast) <- c("reef","q")
  
  # # put things in the right format for simulation inputs
  # a <- subset(values(hab_rast,dataframe=TRUE),!is.na(reef))
  # b <- crds(hab_rast,df=TRUE)
  # reef_locations <- cbind(a,b)
  # 
  # if(plot_flag==TRUE){
  #   g <- ggplot(reef_locations,aes(x=x,y=y,fill=q))+geom_tile()
  #   print(g)
  # }
  
  return(list(q_rast=q_rast))
}

# Generate carrying capacity map
# Inputs:
#   df_patches: data frame with a row for every reef cell, and its x and y values
#   base_rast: SpatRaster object with one layer: reef (0 for open water and 1 for reef)
#   K_range: c(Kmin,Kmax)
#   K_autocorr: some measure of the spatial autocorrelation in K.
#               For the fractal landscape method, it's h in the fracland function (higher values = more autocorrelated, range=(-2,2))
# Outputs:
#   df_patches: df with columns for x, y, and K
#   hab_rast: SpatRaster object (matching base_rast in extent and resolution) with 1 layer: K (carrying capacity)
f_GenerateK <- function(base_rast,K_range,K_autocorr,plot_flag=FALSE){
  # add habitat values on top of base map
  nx=ncol(base_rast)
  ny=nrow(base_rast)
  dimens <- (2^(1:15)+1)   # find the k value to use in fracland function, given the dimensions of the base map
  k <- first(which(dimens>=max(nx,ny)))
  # generate a fractal layer
  frac_map <- fracland(k=k,h=K_autocorr,binary=FALSE,plotflag=FALSE)
  # convert to desired range of q values
  frac_map <- (frac_map-min(frac_map))/(max(frac_map)-min(frac_map)) # first to 0-1
  frac_map <- frac_map*(K_range[2]-K_range[1])+K_range[1] # then to q_range
  frac_map <- frac_map[1:ny,1:nx] 
  empty_rast <- rast(ext(base_rast), resolution=res(base_rast), crs = crs(base_rast))
  values(empty_rast) <- frac_map
  K_rast <- empty_rast*base_rast
  names(K_rast) <- c("K")
  
  if(plot_flag==TRUE){
    plot(K_rast)
  }
  
  return(list(K_rast=K_rast))
}

# # put habitable sites into a dataframe
# a <- values(K_rast,dataframe=TRUE) # K values
# b <- crds(K_rast,df=TRUE) # coordinates
# df_patches <- cbind(b,a) # put them together
# df_patches <- subset(df_patches,K!=0) # remove open water patches


# Generates the specified number of anemones at random locations on the map
# Calculates distance matrix between points
# Inputs:
#   reef_sf: sfc_multipolygon of the reef area
#   hab_rast: SpatRaster object with two layers: reef (0 for open water and 1 for reef) and q (habitat quality)
#   n_anems = number of anemone locations to simulate
# Outputs:
#   df_patches: dataframe with columns for id, x, y, and K (K=1 always for point version of sim)
#   sfc_patches
#   patch_dists
#   K_rast
f_SimPtsOnMap <- function(reef_sf,base_rast,n_anems=50,inwater_dist=FALSE,show_map=TRUE){
  # sample the anemones
  sfc_patches <- st_sample(reef_sf,n_anems)
  # df_patches <- data.frame(st_coordinates(sfc_patches))
  # npatch <- nrow(df_patches)
  # df_patches$id <- 1:npatch
  # df_patches$K <- 1
  
  # make patch_dists
  patch_dists <- st_distance(sfc_patches)
  units(patch_dists) <- 'km'
  
  # K_rast is 1 everywhere
  K_rast <- base_rast
  names(K_rast) <- "K"
  
  if(show_map==TRUE){
    g <- ggplot(reef_sf)+geom_sf()+geom_sf(data=sfc_patches)
    print(g)
  }
  
  return(list(K_rast=K_rast,sfc_patches=sfc_patches,patch_dists=patch_dists))
  
}

# nav_rad = navigation radius (in km). When hab_type="grid", should be set to 1.
f_MakeHabitat2 <- function(nav_rad,q_rast,K_rast,patch_dists,sfc_patches){
  nx=dim(q_rast)[2]
  ny=dim(q_rast)[1]
  npatch <- nrow(df_patches)

  ## put q_rast and K_rast together
  hab_rast <- c(q_rast,K_rast)
    
  ## create df_patches (important: ID should be in the same order as in sfc_patches, or distance matrix will be wrong)
  q_vect <- terra::extract(hab_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)
  K_vect <- terra::extract(hab_rast$K,vect(sfc_patches),xy=TRUE,search_radius=500)
  df_patches <- q_vect
  df_patches$K <- K_vect$K[df_patches$ID]
  df_patches$id <- df_patches$ID
  df_patches$ID <- NULL
  
  ## make patch_angles
  patch_angles <- suppressWarnings(2*asin(nav_rad/patch_dists)/(2*pi))
  patch_angles[is.nan(patch_angles)] <- 1
  
  ## make overlap_discount
  
  
  return(list(npatch=npatch,nx=nx,ny=ny,hab_rast=hab_rast,patch_locations=df_patches,
              patch_dists=patch_dists,patch_angles=patch_angles))
}


# takes existing map (patch_locations)
# or generates uniform grid (if patch_locations==NULL)
# returns connectivity matrices and related objects for use in simulation
# hab_type = "grid" (grid of habitat patches with x-y coordinates), "points" (anemone locations with gps points)
# nav_rad = navigation radius (in km). Used when hab_type="points". When hab_type="grid", it's set to 1 so the patch_angle calculation still works.
# patch_locations can be:
#   1) NULL: generates an nx-x-ny grid
#   2) if hab_type=="grid", a dataframe with columns id, x, y
#   3) if hab_type=="points", a shapefile
# if not given a distance matrix, it'll calculate as the crow flies. If you want in-water distance, do it beforehand and pass the distance matrix.
f_MakeHabitat <- function(nx,ny,v_alphas,v_thetas,patch_locations=NULL,conn_out=FALSE,hab_type="grid",nav_rad=1,dists_mat=NULL,overlap_method=1){
  numCores <- parallelly::availableCores()
  
  # list of patch locations and IDs
  # (dimensions: npatch x 3)
  # if it's being generated here, "location" is the center of the grid square
  if(is.null(patch_locations) & hab_type=="grid"){ # if not specified by an input map, then make one
    patch_locations <- expand.grid(y=1:ny,x=1:nx) %>% # do y first so that patches are ordered columnwise, like the way R fills a matrix
      rowid_to_column(var='id')
  }
  npatch <- nrow(patch_locations)
  
  if(hab_type=="grid"){
    # a "map" of the patch numbers, spatially arranged
    # (dimensions: nx x ny)
    patch_map <- matrix(nrow=ny,ncol=nx)
    for(i in 1:npatch){
      patch_map[patch_locations$y[i],patch_locations$x[i]] <- patch_locations$id[i]
    }
    
    # a matrix of distances between the centers of each patch (i.e., r in polar coords)
    # (dimensions: npatch x npatch)
    if(is.null(dists_mat)){
      patch_dists <- matrix(nrow=npatch,ncol=npatch)
      colnames(patch_dists) <- patch_locations$id
      for(i in 1:npatch){
        patch_dists[i,] = sqrt((patch_locations$x[i]-patch_locations$x)^2+(patch_locations$y[i]-patch_locations$y)^2)
      }
    } else {
      patch_dists <- dists_mat
    }
    
    
    # set objects that are mainly used in points mode
    overlap_discount=rep(1,times=npatch)
    patch_sf=NULL
  }
  
  if(hab_type=="points"){
    # convert formats so patch_locations is a dataframe, and patch_sf is a shapefile
    if(inherits(patch_locations,"sf")==FALSE){     # if patch_locations is not a shapefile
      # create a shapefile of points from patch_locations
      patch_sf <- st_as_sf(patch_locations,coords=c("x","y"))
      st_crs(patch_sf) <- 4326
    } else { # if patch_locations is a shapefile
      patch_sf <- patch_locations
      patch_locations <- cbind(st_drop_geometry(patch_sf),st_coordinates(patch_sf))
    }
    
    # calculate distances between all points, in km
    if(is.null(dists_mat)){
      patch_dists <- st_distance(patch_sf)
      units(patch_dists) <- "km"
      patch_dists <- matrix(as.numeric(patch_dists),nrow=nrow(patch_dists))      
    } else {
      patch_dists <- dists_mat
    }
    
    # ## OVERLAP METHOD 1:
    if(overlap_method==1){
      # if patches aren't on a grid,
      # find the overlap of each patch's basin of attraction with other basins
      # first define the basins
      circs=st_buffer(patch_sf,dist=set_units(nav_rad,km)) # nav_rad is specified in km
      onecirc_area=st_area(circs[1,])
      # then calculate the overlaps (this is slow; should use mclapply)
      all_overlaps <- mclapply(1:npatch,function(i) f_FindOverlapAreas(i,circs,onecirc_area),mc.cores = numCores)
      all_overlaps <- unlist(all_overlaps)
      overlap_discount <- 1/(1+all_overlaps)
    }
    
    ## OVERLAP METHOD 2:
    if(overlap_method==2){
      n_neighbors <- rowSums(patch_dists<nav_rad) # number of points within distance nav_rad of focal point (including focal point)
      overlap_discount <- 1/n_neighbors
    }
    
    
    # set objects that are mainly used in grid mode
    patch_map=NULL
  }
  
  # a matrix of the angle of the pie wedge between each patch
  # (dimensions: npatch x npatch)
  patch_angles <- suppressWarnings(2*asin(nav_rad/patch_dists)/(2*pi))
  patch_angles[is.nan(patch_angles)] <- 1
  
  return(list(patch_locations=patch_locations,
              patch_map=patch_map,
              patch_dists=patch_dists,
              patch_angles=patch_angles,
              npatch=npatch,
              overlap_discount=overlap_discount,
              patch_sf=patch_sf))
}


###### other stuff from f_SimPtsOnMap
  # # get their habitat quality values; store in patch_locations
  # df_patches <- terra::extract(hab_rast[[2]],vect(anemone_spots),xy=TRUE,search_radius=500)
  # # add a column for K (which will be 1 everywhere)
  # df_patches$K <- 1
  # 
  # 
  # 
  # plot(hab_rast[[2]])
  # points(patch_locations$x,patch_locations$y,col='red',pch=19,cex=0.2*patch_locations$q)
  # 
  # if(inwater_dist==TRUE){ # calculate in-water distance, if desired
  #   dists <- lc.dist(marmap_transmat,patch_locations[,c("x","y")],res='dist',meters=TRUE) #distances are in meters
  #   dists_mat <- matrix(0,nrow(patch_locations),nrow(patch_locations)) # convert into matrix form
  #   dists_mat[lower.tri(dists_mat,diag=FALSE)] <- dists
  #   dists_mat <- dists_mat/1000 # in km. Convert after the fact, because if lc.dist works in km it rounds to the nearest km.
  #   dists_mat[upper.tri(dists_mat,diag=FALSE)] <- t(dists_mat)[upper.tri(t(dists_mat),diag=FALSE)] # convert from lower tri to full
  # } else{ # otherwise, calculate euclidean distance
  #   
  # }
  # 
  # # check for infinite distances
  # inf_dists <- which(dists_mat[,1]==Inf)
  # if(length(inf_dists)>0){
  #   print(paste0("Removed ",length(inf_dists)," points at distance infinity"))
  #   patch_locations <- patch_locations[-inf_dists,]
  #   dists_mat <- dists_mat[-inf_dists,][,-inf_dists]
  # }
  # 
  # # check for points too close together
  # v_remove=c()
  # for(pt_i in 1:nrow(patch_locations)){
  #   too_close <- which(dists_mat[-pt_i,pt_i]<0.001) # points less than 1m away
  #   if(length(too_close)>0) {v_remove <- c(v_remove,too_close)}
  # }
  # if(length(v_remove)>0){
  #   patch_locations <- patch_locations[-v_remove,]
  #   dists_mat <- dists_mat[-v_remove,][,-v_remove]
  #   print(paste0("Points too close: ",length(v_remove)," removed"))
  # } 
  # 
  # # optionally, plot
  # if(show_map==TRUE){
  #   g <- ggplot(anemone_spots)+
  #     geom_contour_filled(data=marmap::as.xyz(bathy_raster),aes(x=V1,y=V2,z=V3),alpha=0.5)+
  #     geom_sf(data=reef_sf,fill='black',color='black',alpha=0.2)+  
  #     geom_sf(color='red',size=1)+
  #     annotation_scale()+
  #     theme_minimal()
  #   print(g)
  # }
###### other stuff from f_SimPtsOnMap



# ## plasticity function: plasticity in dispersal kernel in response to K (carrying capacity) of patch
# # b, p, alpha, theta: current set of parameter values
# # b_bad, b_good, b_neutral: values of b that trigger plastic responses
# # n_alpha, n_theta: length of v_alpha and v_theta; used to avoid exceeding the maximum parameter values with plastic response
# f_plasticity <- function(b_i, p_i, alpha_i, theta_i, b_bad=1, b_neutral=5, b_good=9, n_alpha=5, n_theta=5){
#   if(b_good!=b_neutral){ # check for the possibility that there isn't variation in b. Assuming there is:
#     alpha_add <- round((b_neutral-b_i)/(b_good-b_neutral))} # calculate what to add to the alpha index, based on plasticity. It's -1, 0, or +1.
#   else alpha_add <- 0
#   alpha_plastic <- oob_squish(alpha_i+round(p_i)*alpha_add, c(1,n_alpha))
#   theta_plastic <- theta_i
#   return(list(alpha_plastic=alpha_plastic, theta_plastic=theta_plastic))
# }
# f_plasticity2 <- function(b_i, p_i, alpha_i, theta_i, b_bad=1, b_neutral=5, b_good=9, n_alpha=5, n_theta=5){
#   if(b_good!=b_neutral){ # check for the possibility that there isn't variation in b. Assuming there is:
#     alpha_add <- round((b_neutral-b_i)/(b_good-b_neutral))} # calculate what to add to the alpha index, based on plasticity. It's -1, 0, or +1.
#   else alpha_add <- 0
#   alpha_plastic <- oob_squish(alpha_i+round(p_i)*alpha_add, c(1,n_alpha))
#   theta_plastic <- theta_i
#   return(data.frame(alpha_plastic=alpha_plastic,theta_plastic=theta_plastic))
# }

## plasticity function: plasticity in dispersal kernel in response to K (carrying capacity) of patch
## inputs: vectors of values for K, p, alpha, and theta. 
## p should be the actual value of p; alpha and theta should be indices in v_alpha and v_theta
f_plasticityK <- function(K, p, alpha, theta, n_alpha=5, n_theta=5, Kmin=NULL, Kmax=NULL){
  # if alpha and theta are scalars, recycle them to vectors of same length as K
  if(length(alpha)==1) alpha=rep_len(alpha,length(K))
  if(length(theta)==1) theta=rep_len(theta,length(K))
  # define plasticity thresholds of K, if not given
  if(is.null(Kmin)) Kmin=min(K)
  if(is.null(Kmax)) Kmax=max(K)
  # if no variation in K, no plasticity
  if(Kmin==Kmax) alpha_plastic <- alpha
  else {
    alpha_add <- round(ifelse(K<Kmin, p, ifelse(K>Kmax, -p, p-2*p*(K-Kmin)/(Kmax-Kmin))))
    alpha_plastic <- oob_squish(alpha+alpha_add, c(1,n_alpha))
  }
  theta_plastic <- theta
  return(list(alpha_plastic=alpha_plastic,theta_plastic=theta_plastic))
}

## plasticity function: plasticity in dispersal kernel in response to b (reproductive output) of patch
## inputs: vectors of values for b, p, alpha, and theta. 
## p should be the actual value of p; alpha and theta should be indices in v_alpha and v_theta
f_plasticityb <- function(b, p, alpha, theta, n_alpha=5, n_theta=5, bmin=NULL, bmax=NULL){
  # if alpha and theta are scalars, recycle them to vectors of same length as b
  if(length(alpha)==1) alpha=rep_len(alpha,length(b))
  if(length(theta)==1) theta=rep_len(theta,length(b))
  # define plasticity thresholds of b, if not given
  if(is.null(bmin)) bmin=min(b)
  if(is.null(bmax)) bmax=max(b)
  # if no variation in K, no plasticity
  if(bmin==bmax) alpha_plastic <- alpha
  else {
    alpha_add <- round(ifelse(b<bmin, p, ifelse(b>bmax, -p, p-2*p*(b-bmin)/(bmax-bmin))))
    alpha_plastic <- oob_squish(alpha+alpha_add, c(1,n_alpha))
  }
  theta_plastic <- theta
  return(list(alpha_plastic=alpha_plastic,theta_plastic=theta_plastic))
}

# Inputs:
#   base_map: a matrix representing habitat configuration (0=open ocean, 1=reef)
#   if no base map is provided, specify:
#     k (to get dimension of habitat)
#     p (proportion of map covered)
#     h_base (aggregation of habitat area)
#   K_range (vector with 2 elements): range of carrying capacity values the final map should have 
#   h: level of aggregation of K values (K is generated via fractal landscape method)
# Returns:
#   patch_locations: a dataframe of habitable patches, their locations, and their K values
#   K in vector form
#   nx, ny: dimensions of the map
f_GenerateMapWithK <- function(base_map=NULL,K_range,h=0.8,k=5,p=0.3,h_base=0.7,plot_flag=FALSE){
  # if no base map is provided, simulate one
  if(is.null(base_map)){
    base_map=fracland(k=k,h=h_base,p=1-p,binary=TRUE,plotflag=FALSE)
  }
  
  # add habitat values on top of base map
  nx=ncol(base_map)
  ny=nrow(base_map)
  # find the k value to use in fracland function, given the dimensions of the base map
  dimens <- (2^(1:15)+1)
  k <- first(which(dimens>=max(nx,ny)))
  # generate a fractal layer
  frac_map <- fracland(k=k,h=h,binary=FALSE,plotflag=FALSE)
  # convert to desired range of K values
  frac_map <- (frac_map-min(frac_map))/(max(frac_map)-min(frac_map)) # first to 0-1
  frac_map <- frac_map*(K_range[2]-K_range[1])+K_range[1]
  # create full map
  full_map <- base_map
  full_map[which(base_map==TRUE)] <- frac_map[which(base_map==TRUE)]
  
  # put things in the right format for simulation inputs
  patch_locations <- as.data.frame(which(full_map!=0,arr.ind=TRUE)) %>%
    rename(y=row,x=col) %>%
    rowid_to_column(var='id')
  patch_locations$K_i <- round(full_map[full_map!=0])
  K <- patch_locations$K_i
  
  if(plot_flag==TRUE){
    #    dev.new()
    f_Plot_Landscape(patch_locations,nx,ny)
    #    a <- dev.list()
    #    dev.set(which=as.numeric(a['RStudioGD']))
  }
  
  return(list(patch_locations=patch_locations,K=K,nx=nx,ny=ny))
}

f_PlotDecayFn <- function(phi){
  plot(1:100,exp(-phi^2*1:100),type='l',main=paste0('phi=',phi))
}

library(spdep)
# generates spatially-autocorrelated b_i values
# takes a set of points (patch_locations) and the associated distance matrix (dists_mat)
# returns patch_locations, with a column added for b_i
f_SimBValues <- function(patch_locations,dists_mat,phi,show_plot=FALSE,noise=0){
  sp_aut_dist=NA
  # simulate b values
  cov_Sigma <- exp(-phi*dists_mat) # let covariance decay exponentially with distance
  # check in case covariance matrix isn't positive definite (maybe because phi is too small)
  if(min(eigen(cov_Sigma)$values)<(-1e-14)){
    print("Warning: covariance matrix is not positive definite. No b_i values simulated.")
    return(list(patch_locations=patch_locations,sp_aut_dist=sp_aut_dist))
  }
  patch_locations$b_i <- mvrnorm(n=1, mu=rep(5,times=nrow(patch_locations)), Sigma=cov_Sigma+noise*diag(nrow(patch_locations)))
  
  # optionally, plot
  if(show_plot==TRUE){
    g <- ggplot(patch_locations, aes(x=x,y=y,color=b_i))+
      geom_point()+
      labs(title=paste0('phi=',phi))+
      theme_minimal()
    print(g)
  }
  
  # calculate incremental spatial autocorrelation (like Robin's paper)
  dists_mat_NA <- dists_mat
  diag(dists_mat_NA) <- NA
  nearest_dists <- do.call(pmin,c(as.data.frame(dists_mat_NA),na.rm=TRUE))
  #  min_band <- min(nearest_dists) # minimum distance so everybody has at least one neighbor
  min_band <- median(nearest_dists)
  band_incr <- mean(nearest_dists) # average distance to each feature's nearest neighboring feature
  farthest_dists <- do.call(pmax,c(as.data.frame(dists_mat_NA),na.rm=TRUE)) 
  max_band <- min(farthest_dists) # max dist so everybody still has at least one neighbor
  
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
      MC <- Moran.I(x=patch_locations$b_i[nonzero_neigh],weight=w.dist)
      z <- (MC$observed-MC$expected)/MC$sd # calculate the z-score
      incr_moran$z[i] <- z
      incr_moran$MI[i] <- MC$observed
      incr_moran$p[i] <- 2*pnorm(abs(z),lower.tail = FALSE)
    } else incr_moran$MC[i] <- NA
  }  
  
  h <- ggplot(incr_moran,aes(x=band,y=MI))+
    geom_line(alpha=0.3)+
    geom_point(aes(color=(p<0.05)),size=0.5)+
    scale_color_manual(values=c("black","grey"),breaks=c(TRUE, FALSE))+
    labs(x="distance band (km)",y="Moran's I",title=paste0("Incremental Spatial Autocorrelation, phi=",phi))+
    geom_hline(yintercept=0,alpha=0.3,lty="dashed")+
    theme_minimal()
  print(h)
  sp_aut_dist <- incr_moran$band[which.max(incr_moran$z)]
  
  return(list(patch_locations=patch_locations,sp_aut_dist=sp_aut_dist,incr_moran=incr_moran))
}

# inputs: kernel, seascape
# output: rates of dispersal from each patch to each other patch
# matrix Connectivity: dimensions npatch x npatch
# Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i
f_GetConnectivityMatrix <- function(alpha, theta, patch_dists, patch_angles){
  connectivity_matrix <- (pgamma(patch_dists+0.5,shape=alpha,scale=theta)-pgamma(patch_dists-0.5,shape=alpha,scale=theta))*patch_angles
}

# output: rates of dispersal from each patch to each other patch
# Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i
f_GetConnectivityMatrix_parallel <- function(alpha,theta,patch_dists,patch_angles,overlap_discount,nav_rad,numCores){
  npatch=length(alpha)
  connectivity_matrix <- mclapply(1:npatch,function(i) cm_i <- overlap_discount[i]*patch_angles[i,]*
                                    (pgamma(patch_dists[i,]+nav_rad,shape=alpha[i],scale=theta[i])-
                                       pgamma(pmax(patch_dists[i,]-nav_rad,0),shape=alpha[i],scale=theta[i])),
                                  mc.cores=numCores)
  connectivity_matrix <- do.call(rbind,connectivity_matrix)
  return(connectivity_matrix)
}

# inputs:
#   matrices patch_dists and patch_angles
#   vectors of patch-specific values of alpha and theta (if scalar, recycled for use across all patches)
# output:
#   matrix of proportion of individuals that disperse from each patch to each other patch, for the given patch-specific alphas and thetas
f_GetConnectivityMatrix_vectorized <- function(alpha, theta, patch_dists, patch_angles,numCores) {
  np=nrow(patch_dists)
  if(length(alpha)==1) alpha <- rep(alpha,np)
  if(length(theta)==1) theta <- rep(theta,np)
  # create matrices of alpha and theta of the FROM patch
  alpha_mat <- matrix(alpha, nrow=np, ncol=np, byrow=FALSE)
  theta_mat <- matrix(theta, nrow=np, ncol=np, byrow=FALSE)
  # Compute connectivity (row=FROM patch, col=TO patch)
  connectivity_matrix <- patch_angles*(pgamma(patch_dists+0.5, shape=alpha_mat, scale=theta_mat)-
                                         pgamma(patch_dists-0.5, shape=alpha_mat, scale=theta_mat))
  return(connectivity_matrix)
}

## function to run within f_RunMatrixLoop that gets the plastic connectivity matrix for a given parameter group, g, defined by its index
f_GetPlasticConnMat <- function(g, group_index, patch_locations, patch_dists, patch_angles, overlap_discount, v_p, v_alphas, v_thetas,nav_rad,numCores){
  v <- group_index[g,]
  # compute effective parameters for each patch with plasticity (once per group)
  eff_params <- f_plasticityb(patch_locations$b_i, 
                              v_p[v$p], 
                              v$alpha, 
                              v$theta,
                              n_alpha = length(v_alphas),
                              n_theta = length(v_thetas))
  # build matrix
  # conn_mat <- f_GetConnectivityMatrix_vectorized(v_alphas[eff_params$alpha_plastic],
  #                                              v_thetas[eff_params$theta_plastic],patch_dists,patch_angles,numCores)
  conn_mat <- f_GetConnectivityMatrix_parallel(alpha=v_alphas[eff_params$alpha_plastic],
                                               theta=v_thetas[eff_params$theta_plastic],
                                               patch_dists=patch_dists,patch_angles=patch_angles,overlap_discount=overlap_discount,nav_rad=nav_rad,numCores=numCores)
}

# # Take a vector of a time series, and identify the point when it has reached equilibrium
# # this is an ad hoc method, but it seems like it kind of works:
# #   it's at least good at cutting out enough of the early portion of the time series.
# #   less good at detecting smaller changes in equilibrium level after the first move away from initial conditions.
# # inputs: v_t, vector of time series data
# # returns: equil_start, timepoint where equilibrium has been reached
# f_FindEquil <- function(v_t,showplot=FALSE){
#   nt=length(v_t)
#   equil_pt <- NULL
#   
#   for(i in 1:19){ # try up to 19 possible timepoints
#     # split data into training and testing
#     end_train=i*nt/20
#     train <- v_t[1:end_train]
#     test <- v_t[(end_train+1):nt]
#     
#     # make an ARIMA model with the training data
#     m1 <- auto.arima(train,seasonal=FALSE)
#     # use it to forecast for the length of the test data
#     m1_fore <- forecast(m1,nt-end_train,level=c(95))
#     
#     # check if the test data lies in the model's bounds of prediction
#     # if it has reached equilibrium, then we're just predicting the mean with increasingly large confidence intervals -- that's ok.
#     # if it hasn't reached equilibrium yet, then either 1) there's a trend in the training data that doesn't continue in the test data,
#     # or 2) there's a high-order autoregression signal in the training data that means the confidence bounds are really huge
#     # we test for both possibilities below
#     in_bounds <- (test>m1_fore$lower)&(test<m1_fore$upper)
#     pct_in_bounds <- sum(in_bounds)/length(in_bounds)
#     max_CI_width <- as.numeric(tail(m1_fore$upper,1))-as.numeric(tail(m1_fore$lower,1))
#     #test_autocorr_strength <- m1$arma[1]<3
#     test_autocorr_strength <- max_CI_width<(100*diff(range(v_t)))
#     
#     # if so, and if the bounds of prediction aren't enormous, then choose this timepoint and stop the loop
#     # (i.e., the time series has reached equilibrium by the end of the training segment. We can use the test segment as equilibrium data.)
#     if(pct_in_bounds>0.95 & (test_autocorr_strength==TRUE)){
#       equil_pt <- end_train
#       break
#     } 
#   }
#   
#   if(showplot==TRUE & !is.null(equil_pt)){
#     par(mfrow=c(2,1))
#     plot(1:nt,v_t,type='l')
#     lines(1:length(train),train,col='blue')
#     plot(m1_fore)
#     lines(1:nt,v_t)
#     par(mfrow=c(1,1))
#   }
#   
#   if(is.null(equil_pt)) print("Warning: no equilibrium found")
#   return(equil_pt) # if no equilibrium point is found, then the returned value is nt.
# }

# used in lapply in f_MakeHabitat:
# find area of overlap of one circle (j) with all other circles
# assumes a sfc_POLYGON object called circs, with a row for each circle
f_FindOverlapAreas <- function(j,circs,onecirc_area){
  all_intersects <- st_intersection(circs[j,],circs[-j,])
  all_intersects <- st_collection_extract(all_intersects,type="POLYGON")
  all_intersects <- st_make_valid(all_intersects)
  all_intersects <- st_area(all_intersects)
  
  overlap_area <- sum(all_intersects)/onecirc_area
  return(overlap_area)
}
