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
  base_rast <- rast(xmin=0,xmax=endpt_lon,ymin=0,ymax=endpt_lat,resolution=resol)
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
f_MakeHabitat <- function(nav_rad,q_rast,K_rast,patch_dists,sfc_patches,reef_sf,overlap_method="simple"){
  units(nav_rad) <- 'km'
  npatch <- length(sfc_patches)

  ## put q_rast and K_rast together
  hab_rast <- c(q_rast,K_rast)
    
  ## create df_patches (important: ID should be in the same order as in sfc_patches, or distance matrix will be wrong)
  q_vect <- terra::extract(hab_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)
  K_vect <- terra::extract(hab_rast$K,vect(sfc_patches),xy=TRUE,search_radius=500)
  df_patches <- q_vect
  df_patches$K <- K_vect$K[df_patches$ID]
  df_patches$id <- df_patches$ID
  df_patches$ID <- NULL
  df_patches$b <- f_q_to_b(df_patches$q) # calculate reproductive rate (b) from habitat quality (q)
  
  ## make patch_angles
  patch_angles <- suppressWarnings(2*asin(nav_rad/patch_dists)/(2*pi))
  patch_angles[is.nan(patch_angles)] <- 1
  
  ## make overlap_discount
  if(overlap_method=="complicated"){
    # if patches aren't on a grid,
    # find the overlap of each patch's basin of attraction with other basins
    # first define the basins
    circs=st_buffer(sfc_patches,dist=nav_rad) # nav_rad is specified in km
    onecirc_area=st_area(circs[1,])
    # then calculate the overlaps (this is slow; should use mclapply)
    all_overlaps <- mclapply(1:npatch,function(i) f_FindOverlapAreas(i,circs,onecirc_area),mc.cores = parallelly::availableCores())
    all_overlaps <- unlist(all_overlaps)
    overlap_discount <- 1/(1+all_overlaps)
  } else{
    n_neighbors <- rowSums(patch_dists<nav_rad) # number of points within distance nav_rad of focal point (including focal point)
    overlap_discount <- 1/n_neighbors
  }
  
  return(list(npatch=npatch,
              hab_rast=hab_rast,
              patch_locations=df_patches,
              patch_dists=patch_dists,
              patch_angles=patch_angles,
              overlap_discount=overlap_discount,
              reef_sf=reef_sf))
}

# function to calculate reproductive rate (b) from habitat quality (q)
# right now this is boring, but maybe we'll want it to do something more interesting at some point
# input: vector q
# output: vector b
f_q_to_b <- function(q){
  b=as.integer(q)
  return(b)
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

## function to run within f_RunMatrixLoop that gets the plastic connectivity matrix for a given parameter group, g, defined by its index
f_GetPlasticConnMat <- function(g, group_index, patch_locations, patch_dists, patch_angles, overlap_discount, v_p, v_alphas, v_thetas,nav_rad,numCores){
  v <- group_index[g,]
  # compute effective parameters for each patch with plasticity (once per group)
  eff_params <- f_plasticityb(patch_locations$b, 
                              v_p[v$p], 
                              v$alpha, 
                              v$theta,
                              n_alpha = length(v_alphas),
                              n_theta = length(v_thetas))
  # build matrix
  conn_mat <- f_GetConnectivityMatrix_parallel(alpha=v_alphas[eff_params$alpha_plastic],
                                               theta=v_thetas[eff_params$theta_plastic],
                                               patch_dists=patch_dists,patch_angles=patch_angles,overlap_discount=overlap_discount,nav_rad=nav_rad,numCores=numCores)
}


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
