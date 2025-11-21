f_RunMatrixSparse <- function(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,b_bad=1,b_neutral=3,b_good=6,K,disturb_prob=0,patch_locations=NULL){
  starttime <- proc.time()
  
  ########## Data structures to describe space and dispersal ##########
  # see utility functions/f_MakeHabitat for details on what's in each object
  hab <- f_MakeHabitat(nx,ny,v_alphas,v_thetas,patch_locations)
  patch_locations <- hab$patch_locations
  patch_map <- hab$patch_map
  patch_dists <- hab$patch_dists
  patch_angles <- hab$patch_angles
  conn_matrices <- hab$conn_matrices
  npatch <- hab$npatch
  rm(hab)
  if(!"K_i" %in% colnames(patch_locations)) patch_locations$K_i <- as.vector(K)
  patch_locations$b_i <- b_neutral
  
  ########## Data structures to describe population ##########
  # 1. matrix_index: "key" to which parameter values and patch are represented by rows of the following two objects
  # 2. Tij: transition matrix that performs reproduction, dispersal, and mutation steps of simulation
  # 3. Pij: holds vectors of abundance for each patch/parameter combo (rows), at each timestep (columns)
  
  ###### 1. matrix_index, group_index
  ## list of all combinations of patch, alpha, theta, p -- in the order they're represented in Tij and Pij
  # patch (spatial location)
  # alpha index (alpha is the kernel shape parameter; this is the index of the value found in v_alphas)
  # theta (theta is the kernel scale parameter; this is the index of the value found in v_thetas)
  # p (p is the plasticity parameter; this is the index of the value found in v_p)
  matrix_index <- expand.grid(patch=1:npatch,alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  
  ## list out all the "groups" -- i.e., unique combinations of parameters.
  # each of these represents a contiguous (npatch x npatch) section of the transition matrix Tij
  # e.g., the first npatch x npatch cells in Tij represent dispersal from patch i to j,
  #   starting with the first set of parameters in group_index (first "group") and experiencing no mutation.
  # and the cells in the first npatch rows and columns (npatch+1):(2*npatch) represent dispersal from patch i to (npatch+j),
  #   starting with the parameters in the first row of group_index (first "group") and ending with the parameters in the second row of group_index (second "group").
  group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
  
  ###### 2. Tij (only need to build this once at the beginning -- or as often as patch-specific birth rate changes)
  
  # mutation options: each element represents a separate possibility
  # first index represents no mutation
  # index i represents alpha=alpha+alpha_adds[i],theta=theta+theta_adds[i],p=p+p_adds[i]
  alpha_adds=c(0,1,-1,0,0)
  theta_adds=c(0,0,0,1,-1)
  p_adds=c(0,0,0,0,0)
  
  # mutation_destinations: a matrix with dimensions nrow(group_index) x length(alpha_adds)
  # for each row (parameter group), the indices of the parameter groups that could be the result of (+ or -) mutation in alpha or theta
  # (we'll add mutation in p later)
  f_mut_finder <- function(grp,mut_num){
    a <- which(group_index$alpha==group_index$alpha[grp]+alpha_adds[mut_num] & 
                 group_index$theta==group_index$theta[grp]+theta_adds[mut_num] & 
                 group_index$p==group_index$p[grp]+p_adds[mut_num])
    return(ifelse(length(a)!=0,a,grp))
  }
  mutation_destinations <- matrix(nrow=nrow(group_index),ncol=length(alpha_adds))
  for(mut_num in 1:length(alpha_adds)){
    mutation_destinations[,mut_num] <- sapply(1:nrow(group_index),f_mut_finder,mut_num=mut_num)
  }
  
  # timing
  phase1_total <- 0
  phase2_total <- 0

  # fill up Tij
  Tij_list <- list() # we'll add elements of Tij to a list, and then put everything in the sparse matrix at once
  # Tij_df <- data.frame(i=integer(length=nrow(matrix_index)*npatch*8),j=integer(length=nrow(matrix_index)*npatch*8),x=0)
  # Tij_df_at <- 1
  # do the following for every combination of parameters (i.e., "group")
  for(group_row in 1:nrow(group_index)){
    v <- group_index[group_row,] # get parameter values for that row

    start_phase1 <- proc.time()
    
    # get effective kernel parameters, given plasticity
    eff_params <- f_plasticity2(patch_locations$b_i,v_p[v$p],v$alpha,v$theta,b_bad,b_neutral,b_good,n_alpha=length(v_alphas),n_theta = length(v_thetas))
#    eff_params <- f_plasticityK(as.vector(K),v_p[v$p],v$alpha,v$theta,n_alpha=length(v_alphas),n_theta = length(v_thetas))
    
    # patchwise_cmat: a matrix of per-parent larval dispersal rates among patches, specific to the given parameter group
    # build the matrix one row at a time because each row represents a different origin patch
    # origin patches may have different patch quality, and therefore different dispersal kernels
    cmat <- list()
    for(patch_i in 1:npatch){
      # dispersal frequency matrix from the origin patch to everywhere
      # multiplied by b_i to get number of larvae (per adult in origin patch) dispersing from origin patch to everywhere
      cmat[[patch_i]] <- patch_locations$b_i[patch_i]*conn_matrices[eff_params$alpha_plastic[patch_i],eff_params$theta_plastic[patch_i],,patch_i]
    }
    patchwise_cmat <- do.call(rbind,cmat)
    rm(cmat)
    
    phase1_total <- (proc.time()-start_phase1) + phase1_total
    start_phase2 <- proc.time()
    
    # distribute the larvae among their new parameter values (with mutation):
    # place patchwise_cmat multiplied by (1-mu) into Tij for the current parameter group,
    # and (multiplied by mu/4) for the parameter group of each mutation possibility
    Tij_group_inds <- (1:npatch)+npatch*(group_row-1)
    temp_df <- expand.grid(i=Tij_group_inds,j=Tij_group_inds)
    temp_df$x <- as.vector(patchwise_cmat)
    temp_df$i <- as.integer(temp_df$i)
    temp_df$j <- as.integer(temp_df$j)
    
    
    # no mutation
    Tij_list <- append(Tij_list,list(mutate(temp_df,x=(1-mu)*x)))
    # Tij_df[Tij_df_at:(Tij_df_at+nrow(temp_df)-1),] <- mutate(temp_df,x=(1-mu)*x)
    # Tij_df_at <- Tij_df_at+nrow(temp_df)

    # with mutation
    for(mut_group in mutation_destinations[group_row,-1]){ # for each of the 4 possible mutations
      Tij_list <- append(Tij_list,list(mutate(temp_df,x=(mu/4)*x,j=j+npatch*(mut_group-group_row))))
      # Tij_df[Tij_df_at:(Tij_df_at+nrow(temp_df)-1),] <- mutate(temp_df,x=(mu/4)*x,j=j+npatch*(mut_group-group_row))
      # Tij_df_at <- Tij_df_at+nrow(temp_df)
    }
    
    phase2_total <- (proc.time()-start_phase2) + phase2_total
  }
  
  # store Tij in sparse matrix form, and remove everything else
  # for(Tij_part in Tij_list){
  #   write_csv(Tij_part,file="temp_Tij.csv",append=TRUE,col_names=TRUE)
  # }
  # rm(Tij_list)
  # Tij <- sparseMatrix(i=read_csv('temp_Tij.csv',col_select=i),
  #                     j=read_csv('temp_Tij.csv',col_select=j),
  #                     x=read_csv('temp_Tij.csv',col_select=x),
  #                     dims=c(nrow(matrix_index),nrow(matrix_index)))  
  Tij_df <- bind_rows(Tij_list)
  Tij_df <- Tij_df[which(Tij_df$x!=0),]
  rm(Tij_list)
  Tij <- sparseMatrix(i=Tij_df$i,j=Tij_df$j,x=Tij_df$x,dims=c(nrow(matrix_index),nrow(matrix_index)))
  rm(Tij_df)
  
  # print("Tij finished")
  # print(proc.time()-starttime)
  # print("Phase 1:")
  # print(phase1_total)
  # print("Phase 2:")
  # print(phase2_total)
  
  ###### 3. Pij
  # One row for each row of matrix_index and Tij (i.e. each combo of parameter values and patch); one column for each timestep. 
  # (it's faster for R to grab columns than rows)
  # Pij <- matrix(data=0,nrow=nrow(matrix_index),ncol=nsteps)
  Pij <- rep(list(vector(mode="numeric",length=nrow(matrix_index))),nsteps)
  
  ########## Simulation ##########
  
  # get these here, to avoid all the which() statements in the competition loop
  patch_inds_Pij_list <- lapply(1:npatch,function(i_patch)which(matrix_index$patch==i_patch))
  
  # initialize starting population
  # everybody starts with same parameter values and each site at its carrying capacity
  #inds_to_fill <- which(matrix_index$alpha==alpha_start & matrix_index$theta==theta_start & matrix_index$p==p_start)
  # Pij[inds_to_fill,1] <- patch_locations$K_i[matrix_index$patch[inds_to_fill]]
  inds_to_fill <- which(matrix_index$alpha==alpha_start & matrix_index$theta==theta_start & matrix_index$p==p_start)
  Pij[[1]][inds_to_fill] <- patch_locations$K_i[matrix_index$patch[inds_to_fill]]
  
  for(t in 2:nsteps){
    
    # ## Disturbance
    # if(t %% 10 == 0){
    #   # first, reset all K's from any disturbance that occurred 10 timesteps previously
    #   patch_locations$K_i <- as.vector(K)
    #   
    #   # then make a new disturbance (maybe)
    #   if(rbinom(n=1,size=1,p=disturb_prob)==1){
    #     disturb <- ideal.map(ny, nx, p = 0.2, nshape = 1, type = "circle", maxval = 1, minval = 0, binmap = TRUE, rasterflag = FALSE, plotflag=FALSE)
    #     disturb_patches <- as.numeric(na.omit(patch_map[disturb!=0]))
    #     patch_locations$K_i[disturb_patches] <- 0
    #   }
    # }
    
    ## Reproduction and Dispersal and Mutation
    temp_newpop <- Pij[[t-1]] %*% Tij # technically the Pij vector should be a row vector rather than a column vector, but R doesn't care
    
    ## Competition
    # if a patch has population greater than K, sample K individuals and distribute them among cells in that patch
    # (with probability according to the current abundance of each cell)
    for(i_patch in 1:npatch){
      patch_inds_Pij <- patch_inds_Pij_list[[i_patch]] # get the row numbers in Pij corresponding to the focal patch
      patch_abund <- sum(temp_newpop[patch_inds_Pij]) # calculate abundance in the patch
      if(patch_abund>0){ # if abundance is (less than or) equal to 0, the relevant Pij entries for the next timestep are already 0, so we don't need to do anything else
        # choose min(abundance,K) survivors, assigned to parameter values according to the current abundances in those cells
        survivors <- sample(x=length(patch_inds_Pij), # there are this many cells (i.e., combos of parameter values) for the patch
                            size=min(patch_abund,patch_locations$K_i[i_patch]), # choose cells for min(abundance, K) survivors
                            prob = temp_newpop[patch_inds_Pij], # probability of each cell being chosen depends on its current abundance
                            replace=TRUE) |> # same cell can be chosen by multiple individuals
          tabulate() # make into a vector with the number of individuals that chose each cell
        survivors <- c(survivors,rep(0, length(patch_inds_Pij) - length(survivors))) # add zeros to the end for cells past the last one chosen
        Pij[[t]][patch_inds_Pij] <- survivors # put those surviving larvae in the right place in Pij at the new timestep
      }
    }
    
    if(t %% round(nsteps/10) == 0) print(t)
    
  } # t
  
  ########## Process data for plotting ##########
  
  Pij <- do.call(cbind,Pij)
  
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
