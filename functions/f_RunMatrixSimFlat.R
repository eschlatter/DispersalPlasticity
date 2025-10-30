f_RunMatrixSimFlat <- function(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad=1,b_neutral=3,b_good=6,K,sleep_int=0, competition_method='sample'){
  starttime <- proc.time()
  if(!(competition_method %in% c('sample','rnorm'))) stop("competition method incorrectly specified")
  
  ########## Data structures to describe space and dispersal ##########
  # see utility functions/f_MakeHabitat for details on what's in each object
  hab <- f_MakeHabitat(nx,ny,v_alphas,v_thetas)
  patch_locations <- hab$patch_locations
  patch_map <- hab$patch_map
  patch_dists <- hab$patch_dists
  patch_angles <- hab$patch_angles
  conn_matrices <- hab$conn_matrices
  npatch <- hab$npatch
  rm(hab)
  patch_locations$K_i <- as.vector(K)
  patch_locations$b_i <- as.vector(b)
  
  ########## Data structures to describe population ##########
  # 1. matrix_index: "key" to which parameter values and patch are represented by rows of the following two objects
  # 2. Tij: transition matrix that performs reproduction, dispersal, and mutation steps of simulation
  # 3. Pij: holds vectors of abundance for each patch/parameter combo (rows), at each timestep (columns)
  
  ###### 1. matrix_index
  ## list of all combinations of patch, alpha, theta, p -- in the order they're represented in Tij and Pij
    # patch (spatial location)
    # alpha index (alpha is the kernel shape parameter; this is the index of the value found in v_alphas)
    # theta (theta is the kernel scale parameter; this is the index of the value found in v_thetas)
    # p (p is the plasticity parameter; this is the index of the value found in v_p)
  matrix_index <- expand.grid(patch=1:npatch,alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  
  ###### 2. Tij
  # only need to build this once at the beginning (or as often as patch-specific birth rate changes).
  
  # initialize it
  Tij <- matrix(0,nrow=nrow(matrix_index),ncol=nrow(matrix_index))
    
  # first, list out all the "groups" -- i.e., unique combinations of parameters.
  # each of these represents a contiguous (npatch x npatch) section of the transition matrix Tij
  # e.g., the first npatch x npatch cells in Tij represent dispersal from patch i to j,
  #   starting with the first set of parameters in group_index (first "group") and experiencing no mutation.
  # and the cells in the first npatch rows and columns (npatch+1):(2*npatch) represent dispersal from patch i to (npatch+j),
  #   starting with the parameters in the first row of group_index (first "group") and ending with the parameters in the second row of group_index (second "group").
  group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  
  # do the following for every combination of parameters (i.e., "group")
  for(group_row in 1:nrow(group_index)){
    v <- group_index[group_row,] # get parameter values for that row
    
    # make the patch-wise connectivity matrix
    # build it one row at a time because each row represents a different origin patch
    # origin patches may have different patch qualities, and therefore different dispersal kernels
    patchwise_cmat <- matrix(0,nrow=npatch,ncol=npatch)
    for(patch_i in 1:npatch){
      # get the patch quality for the current origin patch
      b_i <- b[patch_i]
      # get the effective kernel parameters, given any plastic response to origin patch quality
      eff_params <- f_plasticity(b_i,v_p[v$p],v$alpha,v$theta,b_bad,b_neutral,b_good,n_alpha=length(v_alphas),n_theta = length(v_thetas))
      # dispersal frequency matrix from the origin patch to everywhere
      cmat <- conn_matrices[eff_params$alpha_plastic,eff_params$theta_plastic,,patch_i]
      # multiply by b_i to get number of larvae (per adult in origin patch) dispersing from origin patch to everywhere
      patchwise_cmat[patch_i,] <- b_i*cmat
    }
    
    # #### remove this later! Testing if rounding makes a difference.
    # patchwise_cmat <- signif(patchwise_cmat,digits=6)
    
    # distribute the larvae among their new parameter values (with mutation):
    # place patchwise_cmat multiplied by (1-mu) into Tij for the current parameter group,
    # and (multiplied by mu/4) for the parameter group of each mutation possibility
    Tij_group_inds <- which(matrix_index$alpha==v$alpha & matrix_index$theta==v$theta & matrix_index$p==v$p) # indices in Tij for that parameter group
    
    # no mutation
    Tij[Tij_group_inds,Tij_group_inds] <- (1-mu)*patchwise_cmat
    
    # mutation
    alpha_adds=c(1,-1,0,0)
    theta_adds=c(0,0,1,-1)
    # find the parameter groups that are the endpoints of mutation
    mutated_param_groups <- mapply(function(x1,x2) which(group_index$alpha==v$alpha+x1 & group_index$theta==v$theta+x2 & group_index$p==v$p)
                                   ,alpha_adds,theta_adds) # (can't vectorize that which statement, so use mapply)
    # identify mutated parameter groups that are out of bounds (i.e., a mutation that decreases alpha when alpha is already at its minimum)
    out_of_bounds <- is.na(mutated_param_groups>0)
    # for each parameter group that's out of bounds, we'll add those who would have undergone that mutation back to the original parameter group
    Tij[Tij_group_inds,Tij_group_inds] <- Tij[Tij_group_inds,Tij_group_inds]+sum(out_of_bounds)*(mu/4)*patchwise_cmat
    # for the mutated parameter groups that are in-bounds:
    # list of the j-indices of each of those parameter groups in Tij
    mut_Tij_inds <- lapply(mutated_param_groups[!out_of_bounds],FUN=function(x) seq(from=((x-1)*npatch+1),to=x*npatch,by=1))
    # store everything where it goes (this shouldn't be a for loop. But it is for now because I'm less likely to make a mistake.)
    for(mut_group in 1:length(mut_Tij_inds)){
      Tij[Tij_group_inds,mut_Tij_inds[[mut_group]]] <- (mu/4)*patchwise_cmat
    }
  }
  
  # Store Tij as a sparse matrix, since it's got tons of zeroes
  # (Can I just build it as sparse to begin with, and avoid using all the memory to create Tij? Probably. Look into this.)
  Tij_sparse <- as(Tij, "sparseMatrix")
  rm(Tij)
  
  ###### 3. Pij
  # One row for each row of matrix_index and Tij (i.e. each combo of parameter values and patch); one column for each timestep. 
  # (it's faster for R to grab columns than rows)
  Pij <- matrix(data=0,nrow=nrow(matrix_index),ncol=nsteps)
  
  ########## Simulation ##########
  
  # initialize starting population
  # everybody starts with same parameter values and each site at its carrying capacity
  inds_to_fill <- which(matrix_index$alpha==alpha_start & matrix_index$theta==theta_start & matrix_index$p==p_start)
  Pij[inds_to_fill,1] <- patch_locations$K_i[matrix_index$patch[inds_to_fill]]
  
  for(t in 2:nsteps){
    # reproduction and dispersal and mutation
    temp_newpop <- Pij[,t-1] %*% Tij_sparse # technically the Pij vector should be a row vector rather than a column vector, but R doesn't care
    
    ## Competition
    # Method 1: if a patch has population greater than K, sample K individuals and distribute them among cells in that patch
    # (with probability according to the current abundance of each cell)
    if(competition_method=='sample'){
      for(i_patch in 1:npatch){
        patch_inds_Pij <- which(matrix_index$patch==i_patch) # get the row numbers in Pij corresponding to the focal patch
        patch_abund <- sum(temp_newpop[patch_inds_Pij]) # calculate abundance in the patch
        if(patch_abund>0){ # if abundance is (less than or) equal to 0, the relevant Pij entries for the next timestep are already 0, so we don't need to do anything else
          # choose min(abundance,K) survivors, assigned to parameter values according to the current abundances in those cells
          survivors <- sample(x=length(patch_inds_Pij), # there are this many cells (i.e., combos of parameter values) for the patch
                              size=min(patch_abund,patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                              prob = temp_newpop[patch_inds_Pij], # probability of each cell being chosen depends on its current abundance
                              replace=TRUE) |> # same cell can be chosen by multiple individuals
            tabulate() # make into a vector with the number of individuals that chose each cell
          survivors <- c(survivors,rep(0, length(patch_inds_Pij) - length(survivors))) # add zeros to the end for cells past the last one chosen
          Pij[patch_inds_Pij,t] <- survivors # put those surviving larvae in the right place in Pij at the new timestep
        }
      }
    }
    
    # # Method 2: if a patch has population greater than K, scale the value in each cell by rnorm(mean = K/(sum of all boxes for that patch))
    # if(competition_method=='rnorm'){
    #   pop_by_patch <- apply(sim_array[,,,,t],1,sum,na.rm=TRUE)
    #   scale_by_patch <- ifelse(pop_by_patch>K,pop_by_patch/K,1) # what to divide the patch population by, to reduce it to carrying capacity
    #   # store the values to scale each cell by
    #   scale_by_cell <- array(dim=c(npatch,length(v_alphas),length(v_thetas),length(v_p),1))
    #   for(i_patch in 1:npatch){
    #     scale_by_cell[i_patch,,,,1] <- pmax(rnorm(n=length(v_alphas)*length(v_thetas)*length(v_p),mean=scale_by_patch[i_patch],sd=mean(b)/10),0)
    #   }
    #   # do the scaling
    #   sim_array[,,,,t] <- sim_array[,,,,t,drop=F]/scale_by_cell
    # }

    if(t %% round(nsteps/10) == 0) print(t)
    
  } # t
  
  ########## Process data for plotting ##########
  
  # melt into a dataframe with columns patch, timestep, alpha, theta, p, popsize
  # add columns for param values at each time/patch/alpha/theta/p combo, scaled by the population size that has that combo
  sim_melt <- cbind(Pij,matrix_index) %>%
    pivot_longer(cols=1:nsteps,names_to='t',values_to = 'popsize') %>%
    mutate(alpha_value=v_alphas[alpha],theta_value=v_thetas[theta],p_value=v_p[p],t=as.numeric(t))
  
  # mean param values at each timepoint
  # first take the sum across all cells of the param*popsize (numerator of the mean)
  by_t <- summarize(group_by(sim_melt,t),
                    alpha=sum(alpha_value*popsize),theta=sum(theta_value*popsize),popsize=sum(popsize))
  # then divide by total popsize (denominator of the mean)
  by_t <- mutate(by_t, alpha=alpha/popsize, theta=theta/popsize)

  time_run <- proc.time()-starttime
  # ######### Output ##########
  return(list(sim_melt=sim_melt,
              by_t=by_t,
              patch_locations=patch_locations,
              time_run=time_run,
              K=K,
              b=b))
}
