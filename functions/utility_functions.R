library(gridExtra)

# generate a simulated basemap
# inputs:
#   x_dist,y_dist: size of map in the x and y directions (in meters)
#   resol: vector of resolution (in meters) in the x and y directions
#   method: "fractal" (fractal landscape); "uniform" (all sites habitable)
#   h: habitat aggregation value (between -2 and 2; higher = more autocorrelated)
#   prop_hab: proportion of the map that should be habitat
#   make_dist_mat: flag to generate distance matrix and sfc_patches here.
#     Not necessary (and potentially slow) if the basemap will only be used for "points"-type simulations,
#     but better to do it here if it'll be used for "patch"-type.
# output:
#   reef_sf: sfc_multipolygon of the reef area
#   bathy_rast: marmap::bathy object, for getting in-water distances
#   base_rast: SpatRaster with 0 for open water, 1 for reef
f_GenerateBasemap <- function(x_dist=500,y_dist=500,resol=c(100,100),
                              method="fractal",h=NA,prop_hab=NA,make_dist_mat=TRUE,
                              plot_flag=FALSE,basemap_file=NULL){
  
  # create empty raster of the appropriate size
  base_rast <- rast(xmin=250000,xmax=250000+x_dist,ymin=0,ymax=y_dist,crs="EPSG:32631",resol=resol)
  nx=ncol(base_rast)
  ny=nrow(base_rast)
  
  # generate habitat configuration (reef vs open water) and add it to base_rast
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
  
  values(base_rast) <- base_map
  
  # create bathy_rast (for getting in-water distances)
  bathy_rast <- marmap::as.bathy(raster(base_rast)) # this rotates it! Why???
  bathy_rast[bathy_rast==1] <- -20 # reef
  bathy_rast[bathy_rast==0] <- -100 # open water
  
  # create reef_sf (for simulating point locations and plotting)
  if(method=="uniform"){
    # generate reef_sf manually
    reef_sf <- st_multipolygon(x = list(list(rbind(c(250000,0),c(250000+x_dist,0),c(250000+x_dist,y_dist),c(250000,y_dist),c(250000,0)))))
    reef_sf <- st_sf(geom = st_sfc(reef_sf),crs="EPSG:32631")
  } else{
    basemap_stars <- st_as_stars(base_rast[[1]])
    basemap_contour <- st_contour(basemap_stars,breaks=c(0.5))
    reef_sf <- basemap_contour[basemap_contour$Min>0,] # pick out just the reef part for the shapefile
  }
  
  if(plot_flag==TRUE){
    g_rast <- ggplot()+ggspatial::layer_spatial(as.factor(base_rast$lyr.1))+annotation_scale()+labs(title="Base Map")+
      scale_fill_manual(values=c("#a6cee3","#d95f02"),name=NULL,labels=c("Water","Reef"))+theme(legend.position = "bottom")
    print(g_rast)
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
    
    # make patch_dists
    patch_dists <- st_distance(sfc_patches,which="Euclidean") # distances in meters
    patch_dists <- drop_units(patch_dists/1000) # convert to km
    units(patch_dists) <- 'km' # re-add units, because there's a bug somewhere when I convert directly
  } else {
    patch_dists=NA
    sfc_patches=NA
  }
  
  if(!is.null(basemap_file)){
    save(reef_sf,patch_dists,sfc_patches,bathy_rast,file=paste0(basemap_file,".RData"))
    writeRaster(base_rast,filename=paste0(basemap_file,".tif"),overwrite=TRUE)
  }
  
  return(list(base_rast=base_rast,bathy_rast=bathy_rast,reef_sf=reef_sf,
              patch_dists=patch_dists,sfc_patches=sfc_patches))
}

# Generate habitat quality map
# Inputs:
#   base_rast: SpatRaster with 0 for open water, 1 for reef. Or the filepath of an existing SpatRaster.
#   q_range: c(qmin,qmax), where q=habitat quality
#   q_autocorr: some measure of the spatial autocorrelation in q.
#               For the fractal landscape method, it's h in the fracland function (higher values = more autocorrelated, range=(-2,2)(technically (-infinity,2) but don't worry about it))
# Outputs:
#   q_rast: SpatRaster object with two layers: reef (0 for open water and 1 for reef) and q (habitat quality)
f_GenerateHabQual <- function(base_rast,q_range,q_autocorr,binary=FALSE,target_dist='identity',
                              plot_flag=FALSE,qmap_file=NULL){
  # if a filepath was specified, load the saved base map. Otherwise it's ready to go.
  if(typeof(base_rast)=="character"){   
    basemap_file <- base_rast
    base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
    print(paste0("using saved base map: ",basemap_file))
  } 
  nx=ncol(base_rast)
  ny=nrow(base_rast)
  # find the k value to use in fracland function, given the dimensions of the base map
  dimens <- (2^(1:15)+1)
  k <- first(which(dimens>=max(nx,ny)))
  # generate a fractal layer
  frac_map <- fracland(k=k,h=q_autocorr,binary=FALSE,plotflag=FALSE)
  frac_map <- frac_map[1:ny,1:nx] 
  temp_rast <- rast(ext(base_rast), resolution=res(base_rast), crs = crs(base_rast))
  # convert to desired distribution of q values
  values(temp_rast) <- matrix(f_TransformDist(frac_map,target_dist),nrow=nrow(frac_map))
  # mask out non-habitat locations
  temp_rast[base_rast==0] <- NA

  # put together with base_rast in new object
  q_rast <- c(base_rast,temp_rast)
  names(q_rast) <- c("reef","q")
  
  if(plot_flag==TRUE){
    q_rast_plot <- q_rast
    q_rast_plot$q[q_rast$q==0] <- NA
    g <- ggplot()+ggspatial::layer_spatial(q_rast_plot$q)+
      scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
      annotation_scale()+labs(title="Habitat quality")
    print(g)
  }
  
  # only save data if 1) qmap_file is given, and 2) the basemap was previously saved
  if(!is.null(qmap_file) & exists("basemap_file")){
    writeRaster(q_rast,filename=paste0(qmap_file,".tif"),overwrite=TRUE)
    save(basemap_file,file=paste0(qmap_file,".RData")) # save the path to the basemap data
  }
  
  return(list(q_rast=q_rast))
}

# Generate carrying capacity map (for "grid" type simulation)
# Inputs:
#   df_patches: data frame with a row for every reef cell, and its x and y values
#   base_rast: SpatRaster object with one layer: reef (0 for open water and 1 for reef)
#   K_range: c(Kmin,Kmax)
#   K_autocorr: some measure of the spatial autocorrelation in K.
#               For the fractal landscape method, it's h in the fracland function (higher values = more autocorrelated, range=(-2,2))
# Outputs:
#   df_patches: df with columns for x, y, and K
#   hab_rast: SpatRaster object (matching base_rast in extent and resolution) with 1 layer: K (carrying capacity)
f_GenerateK <- function(base_rast,K_range,K_autocorr,plot_flag=FALSE,popmap_file=NULL){
  # if a filepath was specified, load the saved base map. Otherwise it's ready to go.
  if(typeof(base_rast)=="character"){   
    basemap_file <- base_rast
    base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
    print(paste0("using saved base map: ",basemap_file))
    
    load(paste0(basemap_file,".RData")) # load reef_sf, patch_dists, sfc_patches
  } else {
    reef_sf <- NULL
    patch_dists <- NULL
    sfc_patches <- NULL
  }
  # add habitat values on top of base map
  nx=ncol(base_rast)
  ny=nrow(base_rast)
  dimens <- (2^(1:15)+1)   # find the k value to use in fracland function, given the dimensions of the base map
  k <- first(which(dimens>=max(nx,ny)))
  # generate a fractal layer
  frac_map <- fracland(k=k,h=K_autocorr,binary=FALSE,plotflag=FALSE)
  # convert to desired range of q values
  # frac_map <- (frac_map-min(frac_map))/(max(frac_map)-min(frac_map)) # first to 0-1
  # frac_map <- frac_map*(K_range[2]-K_range[1])+K_range[1] # then to q_range
  frac_map <- frac_map[1:ny,1:nx] 
  K_rast <- rast(ext(base_rast), resolution=res(base_rast), crs = crs(base_rast))
  values(K_rast) <- frac_map
  K_rast[base_rast==0] <- NA
  # convert to desired range of q values
  K_rast <- (K_rast-minmax(K_rast)[1])/(minmax(K_rast)[2]-minmax(K_rast)[1]) # first to 0-1
  K_rast <- K_rast*(K_range[2]-K_range[1])+K_range[1] # then to K_range
  names(K_rast) <- c("K")
  
  if(plot_flag==TRUE){
    K_rast_plot <- K_rast
    K_rast_plot$K[K_rast$K==0] <- NA
    g <- ggplot()+ggspatial::layer_spatial(K_rast_plot$K)+scale_fill_continuous(palette = 'RedOr',name="K",na.value = "grey")+annotation_scale()+labs(title="Carrying capacity by patch")
    print(g)
  }
  
  hab_type="grid"
  if(!is.null(popmap_file)){
    writeRaster(K_rast,filename=paste0(popmap_file,".tif"),overwrite=TRUE)
    save(reef_sf,patch_dists,sfc_patches,hab_type,basemap_file,file = paste0(popmap_file,".RData"))
  }
  
  return(list(K_rast=K_rast,
              reef_sf=reef_sf,patch_dists=patch_dists,sfc_patches=sfc_patches,
              hab_type=hab_type))
}


# Generates the specified number of anemones at random locations on the map (for "points" type simulation)
# Calculates distance matrix between points
# Inputs:
#   basemap_file: load reef_sf and base_rast from file (this will overwrite objects specified directly)
#   reef_sf: sfc_multipolygon of the reef area
#   base_rast: SpatRaster object with one layer: reef (0 for open water and 1 for reef)
#   n_anems = number of anemone locations to simulate
#   samp_type = "random","regular"
#   inwater_dist = whether to use in-water method to calculate distance, instead of Euclidean
# Outputs:
#   sfc_patches
#   patch_dists
#   K_rast
f_SimPtsOnMap <- function(basemap_file=NULL,reef_sf=NULL,base_rast=NULL,
                          n_anems=50,samp_type="random",inwater_dist=FALSE,
                          plot_flag=FALSE,popmap_file=NULL){
  if(!is.null(basemap_file)){
    load(paste0(basemap_file,".RData")) # load reef_sf, bathy_rast (also sfc_patches and patch_dists, but these will be overwritten)
    base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
    print(paste0("using saved base map: ",basemap_file))
  } else{
    reef_sf <- NULL
  }
  # sample the anemones
  if(samp_type=="regular"){
    sfc_patches <- st_sample(reef_sf,size=round(n_anems),type="regular")
    #st_crs(sfc_patches) <- 4326
  }
  if(samp_type=="random"){
    sfc_patches <- st_sample(reef_sf,size=round(n_anems*1.2),type="random")
    #st_crs(sfc_patches) <- 4326
  }
  # make sure they're all on the reef
  on_reef <- extract(base_rast,sfc_to_df(sfc_patches)[,c("x","y")])
  sfc_patches <- sfc_patches[on_reef$lyr.1==1 & !is.na(on_reef$lyr.1)]
  sfc_patches <- sfc_patches[1:min(length(sfc_patches),n_anems)]
  
  # make patch_dists
  patch_dists <- st_distance(sfc_patches)
  patch_dists <- drop_units(patch_dists/1000) # again, for some reason this is much faster than converting directly
  units(patch_dists) <- 'km'
  
  # K_rast is 1 everywhere
  K_rast <- base_rast
  names(K_rast) <- "K"
  
  if(plot_flag==TRUE){
    g <- ggplot(reef_sf)+geom_sf()+geom_sf(data=sfc_patches)+labs(title="Anemone locations")+theme(panel.background=element_rect(fill="lightblue"),panel.grid = element_blank())
    print(g)
  }
  
  ## output
  hab_type="points"
  if(!is.null(popmap_file)){
    save(reef_sf,sfc_patches,patch_dists,hab_type,basemap_file,file=paste0(popmap_file,".RData"))
    writeRaster(K_rast,filename=paste0(popmap_file,".tif"),overwrite=TRUE)
  }
  
  return(list(K_rast=K_rast,reef_sf=reef_sf,sfc_patches=sfc_patches,patch_dists=patch_dists,hab_type=hab_type))
  
}

# Inputs:
#   nav_rad = navigation radius (in km).
#   overlap_method: how to calculate discount for sites within nav_rad of each other.
#     "simple" (divide by number of sites within nav_rad)
#     or "complicated" (draw all the circles and calculate area of overlap. This is slow.)
#   qmap_file, popmap_file: if both are given, load in from saved data.
#   otherwise, q_rast, K_rast, patch_dists, sfc_patches, reef_sf, and hab_type must be given directly
#   hab_file: output filepath
# Outputs:
#   
f_MakeHabitat <- function(nav_rad,overlap_method="simple",qmap_file=NULL,popmap_file=NULL,
                          q_rast=NULL,K_rast=NULL,patch_dists=NULL,sfc_patches=NULL,reef_sf=NULL,hab_type=NULL,
                          hab_file=NULL){
  # load in data, if necessary
  if(!is.null(qmap_file) & !is.null(popmap_file)){
    load(paste0(qmap_file,".RData"))
    q_basemap <- basemap_file
    load(paste0(popmap_file,".RData")) # loads reef_sf,patch_dists,sfc_patches,hab_type
    if(basemap_file==q_basemap){
      q_rast <- rast(paste0(qmap_file,".tif")) # load q_rast
      K_rast <- rast(paste0(popmap_file,".tif")) # load q_rast
    } else stop("error: habitat quality map and population maps are not compatible")
  }
  
  units(nav_rad) <- 'km'
  npatch <- length(sfc_patches)
  
  ## put q_rast and K_rast together
  hab_rast <- c(q_rast,K_rast)
  
  ## create df_patches (important: ID should be in the same order as in sfc_patches, or distance matrix will be wrong)
  q_vect <- terra::extract(hab_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)
  K_vect <- terra::extract(hab_rast$K,vect(sfc_patches),xy=TRUE,search_radius=500)
  patch_coords <- st_coordinates(sfc_patches)
  df_patches <- cbind(q_vect[,c("ID","q")],patch_coords)
  df_patches$K <- K_vect$K[df_patches$ID]
  df_patches$id <- df_patches$ID
  df_patches$x <- df_patches$X
  df_patches$y <- df_patches$Y
  df_patches <- df_patches[,c("id","x","y","q","K")]
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
  
  hab_params <- list(npatch=npatch,
                     patch_locations=df_patches,
                     patch_dists=patch_dists,
                     patch_angles=patch_angles,
                     overlap_discount=overlap_discount,
                     reef_sf=reef_sf,
                     sfc_patches=sfc_patches,
                     hab_type=hab_type,
                     nav_rad=nav_rad,
                     hab_file=hab_file)
  
  if(!is.null(hab_file)){
    save(hab_params,file=paste0(hab_file,".RData"))
    writeRaster(hab_rast,filename=paste0(hab_file,".tif"),overwrite=TRUE)
  }
  return(hab_params) # note that this doesn't include hab_rast, so if you want this later, you'll need to load it from hab_params$hab_file
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
  if(bmin==bmax){
    alpha_plastic <- alpha
    theta_plastic <- theta
  } 
  else {
    increment_add <- round(ifelse(b<bmin, p, ifelse(b>bmax, -p, p-2*p*(b-bmin)/(bmax-bmin))))
    #alpha_plastic <- oob_squish(alpha+increment_add, c(1,n_alpha))
    theta_plastic <- oob_squish(theta+increment_add, c(1,n_theta))
    alpha_plastic <- alpha
    #theta_plastic <- theta
  }
  return(list(alpha_plastic=alpha_plastic,theta_plastic=theta_plastic))
}

# output: rates of dispersal from each patch to each other patch
# Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i
f_GetConnectivityMatrix_parallel <- function(alpha,theta,patch_dists,patch_angles,overlap_discount,nav_rad,numCores){
  npatch=length(alpha)
  connectivity_matrix <- mclapply(1:npatch,function(i) cm_i <- overlap_discount[i]*patch_angles[i,]*
                                    (pgamma(patch_dists[i,]+nav_rad,shape=alpha[i],scale=theta[i])-
                                       pgamma(pmax(patch_dists[i,]-nav_rad,0),shape=alpha[i],scale=theta[i])+
                                       ifelse(nav_rad>patch_dists[i,],pgamma(nav_rad-patch_dists[i,],shape=alpha[i],scale=theta[i]),0)# when patch_dists<nav_rad, correct for it
                                    ),
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

# Function to transform distributions
# target_dist can be A, B, C, D, E, or "identity"
# A = uniform, B = unimodal intermediate, C = bimodal, D = unimodal small, E = unimodal large
# all distributions have a range of about 1-9
f_TransformDist <- function(starting_dist,target_dist){
  if(target_dist=="identity"){
    return(starting_dist)
  } else{
    load(paste0('data/target_dists/dist_',target_dist,'.RData'))
    # cdf values of each element of the starting distribution
    start_cdf <- as.numeric(as.factor(starting_dist))/length(starting_dist)
    # new value is the quantile of the new distribution at the cdf value from the original
    new_dist <- as.numeric(quantile(target_dist,start_cdf))
    return(new_dist)
  }
}