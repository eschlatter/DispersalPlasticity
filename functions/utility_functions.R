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

f_MakeHabitat <- function(nx,ny,v_alphas,v_thetas){
  # list of patch locations and IDs
  # (dimensions npatch x 3)
  # "location" is the center of the patch
  patch_locations <- expand.grid(y=1:ny,x=1:nx) %>% # do y first so that patches are ordered columnwise, like the way R fills a matrix
    rowid_to_column(var='id')
  # this is a simple version of patch_locations (all the squares of a grid); eventually we'll want to import a map.
  # I think we'll be able to just keep track of reef patches, not open ocean
  npatch <- nrow(patch_locations)
  
  # a "map" of the patch numbers, spatially arranged
  # (dimensions nx x ny)
  patch_map <- matrix(nrow=ny,ncol=nx)
  for(i in 1:npatch){
    patch_map[patch_locations$y[i],patch_locations$x[i]] <- patch_locations$id[i]
  }
  
  # a matrix of distances between the centers of each patch (i.e., r in polar coords)
  # (dimensions npatch x npatch)
  patch_dists <- matrix(nrow=npatch,ncol=npatch)
  colnames(patch_dists) <- patch_locations$id
  for(i in 1:npatch){
    patch_dists[i,] = sqrt((patch_locations$x[i]-patch_locations$x)^2+(patch_locations$y[i]-patch_locations$y)^2)
  }
  
  # a matrix of the size of the pie wedge between each patch (i.e., theta in polar coords -- not theta of the dispersal kernel)
  # assuming the width of the cell at the given distance is 1: not quite correct most of the time, but probably close enough
  # (dimensions npatch x npatch)
  patch_angles <- suppressWarnings(2*asin(1/(2*patch_dists))/(2*pi))
  patch_angles[is.nan(patch_angles)] <- 1
  
  #check: proportion of individuals from each patch that land in a patch (same configuration as patch_map).
  # Shouldn't be greater than 1 anywhere, and should be smaller where fewer patches are reachable.
  cm <- f_GetConnectivityMatrix(1,2,patch_dists,patch_angles)
  t(matrix(rowSums(cm),nrow=nrow(patch_map))) 
  
  # connectivity matrices
  # (dimensions nalpha x ntheta x npatch x npatch)
  # Should we calculate them for each possible kernel up front?
  # There are probably some kernels that won't get used, so this might be a bit wasteful. But, for now, let's do it. We can be more efficient later.
  # Should we allow kernels to evolve past the predefined ones? Maybe. (Probably?) Let's implement this later on.
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