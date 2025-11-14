library(gridExtra)

# for current p, b, alpha, and theta, returns effective alpha and theta values, given plasticity
# b, p, alpha, theta: current set of parameter values
# b_bad, b_good, b_neutral: values of b that trigger plastic responses
# n_alpha, n_theta: length of v_alpha and v_theta; used to avoid exceeding the maximum parameter values with plastic response
f_plasticity <- function(b_i, p_i, alpha_i, theta_i, b_bad=1, b_neutral=5, b_good=9, n_alpha=5, n_theta=5){
  if(b_good!=b_neutral){ # check for the possibility that there isn't variation in b. Assuming there is:
    alpha_add <- round((b_neutral-b_i)/(b_good-b_neutral))} # calculate what to add to the alpha index, based on plasticity. It's -1, 0, or +1.
  else alpha_add <- 0
  alpha_plastic <- oob_squish(alpha_i+round(p_i)*alpha_add, c(1,n_alpha))
  theta_plastic <- theta_i
  return(list(alpha_plastic=alpha_plastic, theta_plastic=theta_plastic))
}

f_plasticity2 <- function(b_i, p_i, alpha_i, theta_i, b_bad=1, b_neutral=5, b_good=9, n_alpha=5, n_theta=5){
  if(b_good!=b_neutral){ # check for the possibility that there isn't variation in b. Assuming there is:
    alpha_add <- round((b_neutral-b_i)/(b_good-b_neutral))} # calculate what to add to the alpha index, based on plasticity. It's -1, 0, or +1.
  else alpha_add <- 0
  alpha_plastic <- oob_squish(alpha_i+round(p_i)*alpha_add, c(1,n_alpha))
  theta_plastic <- theta_i
  return(data.frame(alpha_plastic=alpha_plastic,theta_plastic=theta_plastic))
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
    dev.new()
    f_Plot_Landscape(patch_locations,nx,ny)
    a <- dev.list()
    dev.set(which=as.numeric(a['RStudioGD']))
  }
  
  return(list(patch_locations=patch_locations,K=K,nx=nx,ny=ny))
}

# takes existing map (patch_locations)
# or generates uniform grid (if patch_locations==NULL)
# returns connectivity matrices and related objects for use in simulation
f_MakeHabitat <- function(nx,ny,v_alphas,v_thetas,patch_locations=NULL){
  # list of patch locations and IDs
  # (dimensions: npatch x 3)
  # "location" is the center of the 
  if(is.null(patch_locations)){ # if not specified by an input map, then make one
    patch_locations <- expand.grid(y=1:ny,x=1:nx) %>% # do y first so that patches are ordered columnwise, like the way R fills a matrix
      rowid_to_column(var='id')
  }
  npatch <- nrow(patch_locations)
  
  # a "map" of the patch numbers, spatially arranged
  # (dimensions: nx x ny)
  patch_map <- matrix(nrow=ny,ncol=nx)
  for(i in 1:npatch){
    patch_map[patch_locations$y[i],patch_locations$x[i]] <- patch_locations$id[i]
  }
  
  # a matrix of distances between the centers of each patch (i.e., r in polar coords)
  # (dimensions: npatch x npatch)
  patch_dists <- matrix(nrow=npatch,ncol=npatch)
  colnames(patch_dists) <- patch_locations$id
  for(i in 1:npatch){
    patch_dists[i,] = sqrt((patch_locations$x[i]-patch_locations$x)^2+(patch_locations$y[i]-patch_locations$y)^2)
  }
  
  # a matrix of the size of the pie wedge between each patch (i.e., theta in polar coords -- not theta of the dispersal kernel)
  # assuming the width of the cell at the given distance is 1: not quite correct most of the time, but probably close enough
  # (dimensions: npatch x npatch)
  patch_angles <- suppressWarnings(2*asin(1/(2*patch_dists))/(2*pi))
  patch_angles[is.nan(patch_angles)] <- 1
  
  # connectivity matrices
  # (dimensions: nalpha x ntheta x npatch x npatch)
  # Should we calculate them for each possible kernel up front?
  # There are probably some kernels that won't get used, so this might be a bit wasteful. But, for now, let's do it. We can be more efficient later.
  conn_matrices <- array(NA,dim=c(length(v_alphas),length(v_thetas),npatch,npatch))
  for(i_alpha in 1:length(v_alphas)){
    for(i_theta in 1:length(v_thetas)){
      conn_matrices[i_alpha,i_theta,,] <- f_GetConnectivityMatrix(v_alphas[i_alpha],v_thetas[i_theta],patch_dists,patch_angles)
    } # i_theta
  } # i_alpha
  
  return(list(patch_locations=patch_locations,
              patch_map=patch_map,
              patch_dists=patch_dists,
              patch_angles=patch_angles,
              conn_matrices=conn_matrices,
              npatch=npatch))
}


# inputs: kernel, seascape
# output: rates of dispersal from each patch to each other patch
# matrix Connectivity: dimensions npatch x npatch
# Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i
f_GetConnectivityMatrix <- function(alpha, theta, patch_dists, patch_angles){
  connectivity_matrix <- (pgamma(patch_dists+0.5,shape=alpha,scale=theta)-pgamma(patch_dists-0.5,shape=alpha,scale=theta))*patch_angles
}
