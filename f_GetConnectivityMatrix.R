# inputs: kernel, seascape
# output: rates of dispersal from each patch to each other patch
  # matrix Connectivity: dimensions npatch x npatch
  # Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i

f_GetConnectivityMatrix <- function(alpha, theta, patch_dists, patch_angles){
  connectivity_matrix <- (pgamma(patch_dists+0.5,shape=alpha,scale=theta)-pgamma(patch_dists-0.5,shape=alpha,scale=theta))*patch_angles
}

# sample simple kernel function (though I think we'll just use the built-in gamma function eventually)
f_kernel <- function(distance){
  if(distance<=3) probability = 1/3
  else probability = 0
  probability
}