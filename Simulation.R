library(tidyverse)

source('f_GetConnectivityMatrix.R')

########## Parameters ##########
nx <- 5 # size of space in the x dimension
ny <- 5 # size of space in the y dimension

########## Data structures to describe space and dispersal ##########

# list of patch locations and IDs
# (dimensions npatch x 3)
# "location" is the center of the patch
patch_locations <- expand.grid(x=1:nx,y=1:ny) %>%
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

# a matrix of distances between the centers of each patch (i.e., r)
# (dimensions npatch x npatch)
patch_dists <- matrix(nrow=npatch,ncol=npatch)
colnames(patch_dists) <- patch_locations$id
for(i in 1:npatch){
  patch_dists[i,] = sqrt((patch_locations$x[i]-patch_locations$x)^2+(patch_locations$y[i]-patch_locations$y)^2)
}

# a matrix of the size of the pie wedge between each patch (i.e., theta)
# assuming the width of the cell at the given distance is 1: not quite correct most of the time, but probably close enough
# (dimensions npatch x npatch)
patch_thetas <- 2*asin(1/(2*patch_dists))/(2*pi)
patch_thetas[is.nan(patch_thetas)] <- 1

cm <- f_GetConnectivityMatrix(1,patch_dists,patch_thetas)
t(matrix(rowSums(cm),nrow=5)) #check: proportion of individuals from each patch that land in a patch (same configuration as patch_map).
# Shouldn't be greater than 1 anywhere, and should be smaller where fewer patches are reachable.

