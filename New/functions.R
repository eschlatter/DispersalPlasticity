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
                              plot_flag=FALSE,experiment_folder=NULL,basemap_id=NULL){
  
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
  
  if(!is.null(basemap_id)){
    dir.create(path=paste0(experiment_folder,"/b",basemap_id),recursive=TRUE)
    save(reef_sf,patch_dists,sfc_patches,bathy_rast,file=paste0(experiment_folder,"/b",basemap_id,"/base_b",basemap_id,".RData"))
    writeRaster(base_rast,filename=paste0(experiment_folder,"/b",basemap_id,"/base_b",basemap_id,".tif"),overwrite=TRUE)
  }
  
  return(list(base_rast=base_rast,bathy_rast=bathy_rast,reef_sf=reef_sf,
              patch_dists=patch_dists,sfc_patches=sfc_patches))
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
f_SimPtsOnMap <- function(reef_sf=NULL,base_rast=NULL,n_anems=50,samp_type="random",inwater_dist=FALSE,
                          plot_flag=FALSE,experiment_folder=NULL,basemap_id=NULL,popmap_id=NULL,marmap_transmat=NULL){
  basemap_file <- paste0(experiment_folder,"/b",basemap_id,"/base_b",basemap_id)
  load(paste0(basemap_file,".RData")) # load reef_sf, bathy_rast (also sfc_patches and patch_dists, but these will be overwritten)
  base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
  print(paste0("using saved base map: ",basemap_file))
  
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
  on_reef <- extract(base_rast,sfheaders::sfc_to_df(sfc_patches)[,c("x","y")])
  sfc_patches <- sfc_patches[on_reef$lyr.1==1 & !is.na(on_reef$lyr.1)]
  sfc_patches <- sfc_patches[1:min(length(sfc_patches),n_anems)]
  
  # make patch_dists
  if(inwater_dist==TRUE){
    # first, remove any points that don't have depths in the bathymetry map (more rounding errors)
    # this means the number of anemones may be slightly lower than the input n_anem
    depths <- get.depth(marmap_transmat@history[[3]],x=st_coordinates(sfc_patches),locator=FALSE)
    sfc_patches <- sfc_patches[!is.na(depths$depth)]
    # calculate distances
    anemone_dists <- lc.dist(marmap_transmat,st_coordinates(sfc_patches),res='dist',meters=TRUE) #distances are in km
    anemone_dists_mat <- matrix(0,length(sfc_patches),length(sfc_patches)) # convert dist into matrix form
    anemone_dists_mat[lower.tri(anemone_dists_mat,diag=FALSE)] <- anemone_dists
    anemone_dists_mat[upper.tri(anemone_dists_mat)] <- t(anemone_dists_mat)[upper.tri(anemone_dists_mat)]
    patch_dists <- anemone_dists_mat/1000
    units(patch_dists) <- "km"
  } else{
    patch_dists <- st_distance(sfc_patches)
    patch_dists <- collapse::qM(patch_dists)/1000 # again, for some reason this is much faster than converting directly
    units(patch_dists) <- 'km'
  }
  
  # K_rast is 1 everywhere
  K_rast <- base_rast
  names(K_rast) <- "K"
  
  if(plot_flag==TRUE){
    g <- ggplot(reef_sf)+
      geom_sf()+
      geom_sf(data=sfc_patches)+
      labs(title="Anemone locations")+
      theme(panel.background=element_rect(fill="lightblue"),panel.grid = element_blank())
    print(g)
  }
  
  ## output
  hab_type="points"
  if(!is.null(popmap_id)){
    save(reef_sf,sfc_patches,hab_type,basemap_file,
         file=paste0(experiment_folder,"/b",basemap_id,"/pop_b",basemap_id,"_p",popmap_id,".RData"))
    save(patch_dists,file=paste0(experiment_folder,"/b",basemap_id,"/patchdists_b",basemap_id,"_p",popmap_id,".RData"))
    writeRaster(K_rast,
                filename=paste0(experiment_folder,"/b",basemap_id,"/pop_b",basemap_id,"_p",popmap_id,".tif"),
                overwrite=TRUE)
  }
  
  return(list(K_rast=K_rast,reef_sf=reef_sf,sfc_patches=sfc_patches,patch_dists=patch_dists,hab_type=hab_type))
}

# Generate habitat quality map
# Inputs:
#   base_rast: SpatRaster with 0 for open water, 1 for reef. Or the filepath of an existing SpatRaster.
#   q_autocorr: some measure of the spatial autocorrelation in q.
#               For the fractal landscape method, it's h in the fracland function (higher values = more autocorrelated, range=(-2,2)(technically (-infinity,2) but don't worry about it))
#   target_dist: set a statistical distribution for the q values in the map to follow.
#               "identity" = leave as it is; 'A' through 'E' are distribution choices, all with range q=1-9, as follows:
#               A = uniform, B = intermediate unimodal, C = bimodal, D = low unimodal, E = high unimodal.
#   qmap_file: where to save the output
# Outputs:
#   q_rast: SpatRaster object with two layers: reef (0 for open water and 1 for reef) and q (habitat quality)
f_GenerateHabQual <- function(q_autocorr,target_dist='identity',plot_flag=FALSE,experiment_folder=NULL,
                              basemap_id=NULL,qmap_id=NULL){
  # if a filepath was specified, load the saved base map. Otherwise it's ready to go.
  basemap_file <- paste0(experiment_folder,"/b",basemap_id,"/base_b",basemap_id)
  base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
  print(paste0("using saved base map: ",basemap_file))
  
  nx=ncol(base_rast)
  ny=nrow(base_rast)
  # find the k value to use in fracland function, given the dimensions of the base map
  dimens <- (2^(1:15)+1)
  k <- first(which(dimens>=max(nx,ny)))
  # generate a fractal layer
  frac_map <- fracland(k=k,h=q_autocorr,binary=FALSE,plotflag=FALSE)
  frac_map <- frac_map[1:ny,1:nx] 
  # mask out non-habitat locations
  frac_map[is.na(matrix(values(base_rast),nrow=ny,byrow=TRUE))] <- NA # if they're NA's
  frac_map[matrix(values(base_rast),nrow=ny,byrow=TRUE)==0] <- NA # if they're zeros
  temp_rast <- rast(ext(base_rast), resolution=res(base_rast), crs = crs(base_rast))
  # convert to desired distribution of q values
  values(temp_rast) <- f_TransformDist(frac_map,target_dist)
  
  # put together with base_rast in new object
  q_rast <- c(base_rast,temp_rast)
  names(q_rast) <- c("reef","q")
  
  if(plot_flag==TRUE){
    q_rast_plot <- q_rast
    q_rast_plot$q[q_rast$reef==0] <- NA
    g <- ggplot()+ggspatial::layer_spatial(q_rast_plot$q)+
      scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
      annotation_scale()+labs(title="Habitat quality")
    print(g)
  }
  
  # get variogram info
  load(paste0(basemap_file,".RData")) # get reef_sf
  sfc_patches <- sf::st_sample(reef_sf,units::drop_units(st_area(reef_sf)/1000))
  spdf1 <- as_Spatial(sfc_patches)
  spdf1$q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
  vgm1 <- gstat::variogram(q~1,data=spdf1,cressie=TRUE) # empirical variogram
  vgmf <- gstat::fit.variogram(vgm1,gstat::vgm(c("Gau","Sph","Exp"))) # run several models, and pick the best one
  #plot(vgm1,vgmf)
  vgm_fit=list(range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
               model=vgmf$model[2])
  
  # only save data if 1) qmap_file is given, and 2) the basemap was previously saved
  if(!is.null(qmap_id)){
    writeRaster(q_rast,filename=paste0(experiment_folder,"/b",basemap_id,"/qmap_b",basemap_id,"_q",qmap_id,".tif"),overwrite=TRUE)
  }
  
  return(list(q_rast=q_rast,vgm_fit=vgm_fit))
}

# Function to transform distributions
# target_dist can be A, B, C, D, E, or "identity"
# A = uniform, B = unimodal intermediate, C = bimodal, D = unimodal small, E = unimodal large
# all distributions have a range of about 0-1
# returns an object with the same dimensions as starting_dist
f_TransformDist <- function(starting_dist,target_dist){
  # get just the non-NA cells to transform
  good_inds <- which(!is.na(starting_dist))
  starting_dist_use <- starting_dist[good_inds]
  
  # do the transformation
  if(target_dist=="identity"){ 
    # just converts it to range 0-1; doesn't change the distribution otherwise
    new_dist_use <- (starting_dist_use-min(starting_dist_use))/(max(starting_dist_use)-min(starting_dist_use))
  } else{
    # convert to new distribution
    load(paste0('data/target_dists/dist_',target_dist,'.RData'))
    # cdf values of each element of the starting distribution
    start_cdf <- as.numeric(as.factor(starting_dist_use))/length(starting_dist_use)
    # new value is the quantile of the new distribution at the cdf value from the original
    new_dist_use <- as.numeric(quantile(target_dist,start_cdf))
    new_dist_use <- (new_dist_use-min(new_dist_use))/(max(new_dist_use)-min(new_dist_use))
  }
  
  # store back in object with the same dimensions and NA's as starting_dist
  new_dist <- rep(NA,length(starting_dist))
  new_dist[good_inds] <- new_dist_use
  dim(new_dist) <- dim(starting_dist)
  return(new_dist)
}

f_MakeHabitat <- function(qmapID=NULL,popmapID=NULL,experiment_folder=NULL){
  # load in data
  qmap_row <- read.csv(paste0(experiment_folder,"/index_qmaps.csv")) |> 
    dplyr::filter(qmap_id==qmapID)
  popmap_row <- read.csv(paste0(experiment_folder,"/index_popmaps.csv")) |> 
    dplyr::filter(popmap_id==popmapID)
  if(qmap_row$basemap_id==popmap_row$basemap_id){basemapID <- qmap_row$basemap_id
  } else stop(paste0("Incompatible base maps: Popmap ",popmapID," basemap ",
                     popmap_row$basemap_id," vs. qmap ", qmapID, " basemap ",
                     qmap_row$basemap_id))
  
  # check if this qmap/popmap combo has been done already
  hab_row <- read.csv(paste0(experiment_folder,"/index_habs.csv")) |>
    dplyr::filter(qmap_id==qmapID & popmap_id==popmapID)
  if(nrow(hab_row)>0){
    stop(paste0("Already done: Hab ",hab_row$hab_id,", Qmap ",qmapID,", Popmap ",popmapID))
  } else{
    habID <- max(read.csv("New/Maps/index_habs.csv")$hab_id,0)+1
    print(paste0("Generating: Hab ",habID,", Qmap ",qmapID,", Popmap ",popmapID))
  }

  # load files
  qmap_file <- paste0(experiment_folder,"/b",basemapID,"/qmap_b",basemapID,"_q",qmapID)
  popmap_file <- paste0(experiment_folder,"/b",basemapID,"/pop_b",basemapID,"_p",popmapID)
  hab_file <- paste0(experiment_folder,"/habfiles/hab_",habID)
  #patchdist_file <- paste0(experiment_folder,"/b",basemapID,"/patchdists_b",basemapID,"_p",popmapID,".RData")
  
  load(paste0(popmap_file,".RData")) # loads reef_sf,sfc_patches,hab_type
  q_rast <- rast(paste0(qmap_file,".tif")) # load q_rast
  K_rast <- rast(paste0(popmap_file,".tif")) # load K_rast
  #load(patchdist_file) # load patch_dists
  
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
  
  hab_params <- list(npatch=npatch,
                     patch_locations=df_patches,
                     reef_sf=reef_sf,
                     sfc_patches=sfc_patches,
                     habID=habID,
                     q_autocorr_scale=qmap_row$range)
  
  save(hab_params,file=paste0(hab_file,".RData"))
  writeRaster(hab_rast,filename=paste0(hab_file,".tif"),overwrite=TRUE)
  
  index_habs <- data.frame(hab_id=habID,basemap_id=basemapID,qmap_id=qmapID,
                              popmap_id=popmapID,npatch=npatch,q_autocorr_scale=qmap_row$range)
  fwrite(index_habs,file=paste0(experiment_folder,"/index_habs.csv"))
}


#' Create neutral landscape maps
#' 
#' Use standard methods to generate fractal maps. Binary and continuous surfaces may be produced.
#' 
#' @param k integer. The extent of the map (2^k+1)^2 pixels
#' @param h numeric. Level of aggregation in the map.
#' @param p numeric (0,1). The proportion of map in habitat 1
#' @param binary logical. If TRUE, a 0/1 categorical landscape is produced.
#' @param plotflag logical. If TRUE, the map will be plotted.
#' @param rasterflag logical. If TRUE, a spatially-referenced raster will be returned.
#' @param minx numeric. The minimum x coordinate of the raster.
#' @param miny numeric. The minimum y coordinate of the raster.
#' @param cellsize integer. The number of spatial units per raster pixel. 
#' @param projstr string. The proj4string describing the spatial reference of the raster.
#' @author Shannon Pittman, James Forester
#' @export
#' @example examples/neutral.landscape_example.R
fracland <- function(k, h, p, binary = TRUE, plotflag = FALSE, rasterflag = FALSE, minx = 562569, 
                     miny = 5230469, cellsize = 30, projstr = "+proj=utm +zone=15 +datum=NAD83",  range= NULL, ...) {
  ## Function for creating neutral landscapes Shannon Pittman University of Minnesota May, 2013 k = the extent of the map (2^k+1)^2 pixels h =
  ## how clumped the map should be (ranging from ?? to ??) -- weird behavior at higher values p = proportion of map in habitat 1 binary =
  ## plotflag == if TRUE will plot a filled contour version of the matrix
  
  ## function call: testmap=land(6,1,.5,FALSE,TRUE)
  library(sp)
  
  # k <- 6 # Scalar-determines size of landscape
  A <- 2^k + 1  # Scalar-determines length of landscape matrix
  
  # if(algorithm=="Gaussian"){
  #   if(is.null(range)) range = A*h
  #   B=ideal.map(A,A,p=p, range=range, binmap=binary, ...)
  #   
  # }else{
  
  B <- matrix(0, A, A)  # Creates landscape matrix
  #HabValue <- 0.5  # Determines value assigned to 'Habitat' points
  # h <- 0.9 # how clumped the landscape is
  #PixelsPerMeter <- 1  #set scale of landscape
  # p <- 0.5 #percentage of landscape that is 'Habitat'
  B[1, 1] <- 0
  B[1, A] <- 0
  B[A, 1] <- 0
  B[A, A] <- 0
  
  
  iter <- 1
  for (iter in 1:k) {
    scalef <- (0.5 + (1 - h)/2)^(iter)
    
    d <- 2^(k - iter)
    
    # ALL SQUARE STEPS#
    for (i in seq(d + 1, A - d, 2 * d)) {
      for (j in seq(d + 1, A - d, 2 * d)) {
        B[i, j] <- mean(c(B[i - d, j - d], B[i - d, j + d], B[i + d, j - d], B[i + d, j + d])) + scalef * rnorm(n = 1)
      }
    }
    
    # OUTSIDE DIAMOND STEP#
    for (j in seq(d + 1, A - d, 2 * d)) {
      B[1, j] <- mean(c(B[1, j - d], B[1, j + d], B[1 + d, j])) + scalef * rnorm(n = 1)
      B[A, j] <- mean(c(B[A, j - d], B[A, j + d], B[A - d, j])) + scalef * rnorm(n = 1)
    }
    
    for (i in seq(d + 1, A - d, 2 * d)) {
      B[i, 1] <- mean(c(B[i - d, 1], B[i + d, 1], B[i, 1 + d])) + scalef * rnorm(n = 1)
      B[i, A] <- mean(c(B[i - d, A], B[i + d, A], B[i, A - d])) + scalef * rnorm(n = 1)
    }
    
    # INSIDE DIAMOND STEP#
    if (2 * d + 1 <= A - 2 * d) {
      for (i in seq(d + 1, A - d, 2 * d)) {
        for (j in seq(2 * d + 1, A - 2 * d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
        }
      }
      
      for (i in seq(2 * d + 1, A - 2 * d, 2 * d)) {
        for (j in seq(d + 1, A - d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
        }
      }
    }
    
    iter <- iter + 1
  }
  
  if (binary == T) {
    R <- sort(B)
    PosR <- (1 - p) * length(R)  #larger values become habitat, designated as 0
    pval <- R[PosR]
    T1 <- which(B > pval)
    T2 <- which(B <= pval)
    B[T1] <- 0  #habitat is 0
    B[T2] <- 1
    if (plotflag) 
      filled.contour(B, levels = c(0, 0.5, 1), col = c("black", "white"))
  } else {
    if (plotflag) 
      filled.contour(B)
  }
  #  }
  
  if (rasterflag) {
    library(raster)
    mapdat <- list()
    mapdat$x <- seq(minx+0.5*cellsize, by = cellsize, len = ncol(B))
    mapdat$y <- seq(miny+0.5*cellsize, by = cellsize, len = nrow(B))
    mapdat$z <- B
    B <- raster(mapdat$z, xmn = range(mapdat$x)[1]-0.5*cellsize, xmx = range(mapdat$x)[2]+0.5*cellsize, ymn = range(mapdat$y)[1]-0.5*cellsize, ymx = range(mapdat$y)[2]+0.5*cellsize, crs = CRS(projstr))
    
  }
  return(B)
  
}
