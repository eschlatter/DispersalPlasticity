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

f_plasticityK <- function(K_i, p, alpha, theta, n_alpha=5, n_theta=5){
  if(length(unique(K_i))>1){ # check for the possibility that there isn't variation in K. Assuming there is:
    bad_thresh <- min(K_i)+(1/3)*(max(K_i)-min(K_i))
    good_thresh <- min(K_i)+(2/3)*(max(K_i)-min(K_i))
    alpha_add <- case_when(K_i<bad_thresh ~ -1, #  decide what to add to the alpha index, based on plasticity. It's -1, 0, or +1.
                           K_i<=good_thresh ~ 0,
                           K_i>good_thresh ~ 1)
  }
  else alpha_add <- rep(0,length(K_i))
  alpha_plastic <- oob_squish(alpha+round(p)*alpha_add, c(1,n_alpha))
  theta_plastic <- theta
  return(data.frame(alpha_plastic=alpha_plastic,theta_plastic=theta_plastic))
}

## improve the plasticity function
## it still takes vectors of values for K, p, alpha, and theta. 
## p should be the actual value of p; alpha and theta should be indices in v_alpha and v_theta
f_plasticityK_new <- function(K, p, alpha, theta, n_alpha=5, n_theta=5, Kmin=NULL, Kmax=NULL){
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

# takes existing map (patch_locations)
# or generates uniform grid (if patch_locations==NULL)
# returns connectivity matrices and related objects for use in simulation
f_MakeHabitat <- function(nx,ny,v_alphas,v_thetas,patch_locations=NULL,conn_out=FALSE){
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
  if(conn_out==TRUE){
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
  
  return(list(patch_locations=patch_locations,
              patch_map=patch_map,
              patch_dists=patch_dists,
              patch_angles=patch_angles,
              npatch=npatch))
}


# inputs: kernel, seascape
# output: rates of dispersal from each patch to each other patch
# matrix Connectivity: dimensions npatch x npatch
# Connectivity[i,j] = the proportion of dispersers from patch j that land in patch i
f_GetConnectivityMatrix <- function(alpha, theta, patch_dists, patch_angles){
  connectivity_matrix <- (pgamma(patch_dists+0.5,shape=alpha,scale=theta)-pgamma(patch_dists-0.5,shape=alpha,scale=theta))*patch_angles
}


f_GetConnectivityMatrix_parallel <- function(alpha,theta,patch_dists,patch_angles,numCores){
  npatch=length(alpha)
  list_cms <- list(length=npatch)
  connectivity_matrix <- mclapply(1:npatch,function(i) cm_i <- patch_angles[i,]*(pgamma(patch_dists[i,]+0.5,shape=alpha[i],scale=theta[i])-
                                                                                   pgamma(patch_dists[i,]-0.5,shape=alpha[i],scale=theta[i])),
                                  mc.cores=numCores)
  connectivity_matrix <- do.call(rbind,connectivity_matrix)
  colnames(connectivity_matrix) <- 1:npatch
  return(connectivity_matrix)
}

# for(i in 1:length(alpha)){ # for each from patch
#   # vector of proportion of individuals from patch i that land in each patch
#   cm_i <- patch_angles[i,]*(pgamma(patch_dists[i,]+0.5,shape=alpha[i],scale=theta[i])-
#                               pgamma(patch_dists[i,]-0.5,shape=alpha[i],scale=theta[i]))
#   list_cms[[i]] <- cm_i
# }
# connectivity_matrix <- matrix(unlist(list_cms),nrow=length(alpha))
#comp_results <- mclapply(1:npatch,function(i) f_Competition(i,patch_abunds,patch_locations,temp_pop),mc.cores = numCores)

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
f_GetPlasticConnMat <- function(g, group_index, patch_locations, patch_dists, patch_angles, v_p, v_alphas, v_thetas,numCores){
  v <- group_index[g,]
  # compute effective parameters for each patch with plasticity (once per group)
  eff_params <- f_plasticityK_new(patch_locations$K_i, 
                                  v_p[v$p], 
                                  v$alpha, 
                                  v$theta,
                                  n_alpha = length(v_alphas),
                                  n_theta = length(v_thetas))
  # build matrix
  conn_mat <- f_GetConnectivityMatrix_parallel(v_alphas[eff_params$alpha_plastic],
                                               v_thetas[eff_params$theta_plastic],patch_dists,patch_angles,numCores)
}

# Take a vector of a time series, and identify the point when it has reached equilibrium
# this is an ad hoc method, but it seems like it kind of works:
#   it's at least good at cutting out enough of the early portion of the time series.
#   less good at detecting smaller changes in equilibrium level after the first move away from initial conditions.
# inputs: v_t, vector of time series data
# returns: equil_start, timepoint where equilibrium has been reached
f_FindEquil <- function(v_t,showplot=FALSE){
  nt=length(v_t)
  equil_pt <- NULL
  
  for(i in 1:19){ # try up to 19 possible timepoints
    # split data into training and testing
    end_train=i*nt/20
    train <- v_t[1:end_train]
    test <- v_t[(end_train+1):nt]
    
    # make an ARIMA model with the training data
    m1 <- auto.arima(train,seasonal=FALSE)
    # use it to forecast for the length of the test data
    m1_fore <- forecast(m1,nt-end_train,level=c(95))
    
    # check if the test data lies in the model's bounds of prediction
    # if it has reached equilibrium, then we're just predicting the mean with increasingly large confidence intervals -- that's ok.
    # if it hasn't reached equilibrium yet, then either 1) there's a trend in the training data that doesn't continue in the test data,
    # or 2) there's a high-order autoregression signal in the training data that means the confidence bounds are really huge
    # we test for both possibilities below
    in_bounds <- (test>m1_fore$lower)&(test<m1_fore$upper)
    pct_in_bounds <- sum(in_bounds)/length(in_bounds)
    max_CI_width <- as.numeric(tail(m1_fore$upper,1))-as.numeric(tail(m1_fore$lower,1))
    #test_autocorr_strength <- m1$arma[1]<3
    test_autocorr_strength <- max_CI_width<(100*diff(range(v_t)))
    
    # if so, and if the bounds of prediction aren't enormous, then choose this timepoint and stop the loop
    # (i.e., the time series has reached equilibrium by the end of the training segment. We can use the test segment as equilibrium data.)
    if(pct_in_bounds>0.95 & (test_autocorr_strength==TRUE)){
      equil_pt <- end_train
      break
    } 
  }
  
  if(showplot==TRUE & !is.null(equil_pt)){
    par(mfrow=c(2,1))
    plot(1:nt,v_t,type='l')
    lines(1:length(train),train,col='blue')
    plot(m1_fore)
    lines(1:nt,v_t)
    par(mfrow=c(1,1))
  }
  
  if(is.null(equil_pt)) print("Warning: no equilibrium found")
  return(equil_pt) # if no equilibrium point is found, then the returned value is nt.
}
