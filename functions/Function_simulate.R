### Core simulation function.

SimFn <- function(repID, # for multi-simulation experiments, to track repetitions of the same parameter set
                  experimentID, # for multi-simulation experiments, to track parameter sets
                  output_flag="lite", # "all": keep data on all individuals at each timestep (can result in very large data files), or "lite": keep summary statistics
                  output_thin=1, # for output_flag="all", interval between recorded timesteps (e.g., 1=every timestep is recorded; 10=every 10th timestep is recorded)
                  experiment_folder,  # for multi-simulation experiments, name of relevant directory
                  notes=NA,   # can include a character string here with notes on the sim
                  
                  # bio parameters
                  mutation_type, # "frequent" or "rare"
                  mu, # mutation magnitude (for plasticity)
                  mu_theta, # mutation magnitude (for dispersal kernel parameter)
                  nav_rad, # navigation radius (in km)
                  adult_survival_prob, # adult survival probability
                  base_fecund, # base fecundity; altered for each individual based on habitat quality
                  habID, # habitat ID: specifies the base map and habitat quality layer to use for this simulation
                  
                  # sim parameters
                  theta_start_min, # minimum initial value (of uniform distribution) for dispersal kernel parameter
                  theta_start_max, # maximum initial value (of uniform distribution) for dispersal kernel parameter
                  p_start_min, # minimum initial value (of uniform distribution) for plasticity parameter
                  p_start_max, # maximum initial value (of uniform distribution) for plasticity parameter
                  nsteps, # number of timesteps to simulate
                  
                  # sim options
                  normalize_offspring=FALSE, # option to force all individuals to produce the same number of offspring (useful for testing)
                  plasticity_on=TRUE, # option to turn off plasticity (useful for testing)
                  larval_output_by_theta=FALSE, # save generic data about the relationship between kernel parameter and larval output
                  
                  # disturbance parameters
                  Dp=0, # probability of disturbance (per timestep)
                  De=0, # extent of disturbance
                  Dl=1, # duration of disturbance
                  disturb_method="fractal", # disturbance type: fractal landscape shape (new shape generated for each disturbance) or circular (new center point chosen at random for each disturbance)
                  Dm=0.1 # magnitude of disturbance (fraction each habitat quality is multiplied by)
){
  
  #### ------------------------Simulation metadata -----------------------------
  # All params
  df_params <- data.frame(mu=mu,nav_rad=nav_rad,adult_survival_prob=adult_survival_prob,
                          base_fecund=base_fecund,
                          theta_start_min=theta_start_min,theta_start_max=theta_start_max,
                          p_start_min=p_start_min,p_start_max=p_start_max,
                          nsteps=nsteps,output_thin=output_thin,
                          normalize_offspring=normalize_offspring,plasticity_on=plasticity_on)
  # generate simID
  while(exists("simID")==FALSE){
    simID <- paste0(sprintf("%08d", sample(10000000:99999999,size=1)))
    if(simID %in% read.csv("_index_sims.csv")$sim_id) rm(simID) # check that it hasn't been used already
  }
  print(paste0("simID: ",simID))
  simDate <- format(Sys.time(),"%Y_%m_%d %H:%M:%S")
  slurmJob <- as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID")) # as.numeric(Sys.getenv("SLURM_JOB_ID"))
  slurmTask <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if(is.na(slurmTask)) slurmTask <- 1
  df_index <- data.frame(sim_id=paste(simID),simDate=simDate,slurmJob=paste0(slurmJob,"_",slurmTask),
                         rep_id=repID,notes=notes)
  # Habitat info
  hab_row <- read.csv("maps/index_habs.csv") |> 
    dplyr::filter(hab_id==habID)
  # Disturbance parameters
  disturb_params <- data.frame(Dp=Dp,De=De,Dl=Dl,disturb_method=disturb_method, Dm=Dm,
                               mu_theta=mu_theta,experiment_id=experimentID,mutation_type=mutation_type)
  # Write all sim info (params, metadata, and habitat) to index
  index_sims <- cbind(df_index,hab_row,df_params,disturb_params)
  fwrite(index_sims,file="_index_sims.csv",append = TRUE)
  print(t(index_sims))
  
  # output file name
  output_file <- paste0("output/raw/",simID)
  
  # set seed with simID
  set.seed(as.numeric(simID),kind="L'Ecuyer-CMRG")
  
  #### -------------------Set up and precalculate objects ----------------------
  starttime <- proc.time()
  numCores=max(1,parallelly::availableCores())
  
  # kernel lookup table
  load(paste0("maps/b",hab_row$basemap_id,"/PopmapObjects_b",hab_row$basemap_id,"_p",hab_row$popmap_id,"_new.RData")) # includes offmap correction 
  load(paste0("maps/b",hab_row$basemap_id,"/NavradObjects_b",hab_row$basemap_id,"_p",hab_row$popmap_id,"_n=",nav_rad,".RData"))
  theta_max <- max(v_theta_bins)
  theta_min <- min(v_theta_bins)
  
  # patch_locations
  load(paste0("maps/habfiles/hab_",habID,".RData"))
  patch_locations <- hab_params$patch_locations
  patch_locations$b <- round(base_fecund*patch_locations$q)
  patch_locations$q_original <- patch_locations$q # store original habitat quality, to return to after disturbance
  
  # disturbance objects
  disturb_length <- 0 # initialize. how long has the disturbance been going?
  disturbed_anems <- c()
  if(disturb_method=="fractal"){
    sfc_patches <- hab_params$sfc_patches
    basemap_file <- paste0("maps/b",hab_row$basemap_id,"/base_b",hab_row$basemap_id)
    base_rast <- rast(paste0(basemap_file,".tif")) # load base_rast
    while(prod(dim(base_rast))>1e4){    # make it lower-resolution, so calculating the disturbances goes faster
      base_rast <- aggregate(base_rast,2)
    }
    nx=ncol(base_rast)
    ny=nrow(base_rast)
    dimens <- (2^(1:15)+1)
    k <- first(which(dimens>=max(nx,ny)))
    
  } else{
    xmin=min(patch_locations$x)
    xmax=max(patch_locations$x)
    ymin=min(patch_locations$y)
    ymax=max(patch_locations$y)
    Dr <- sqrt(((xmax-xmin)*(ymax-ymin)*De)/pi) # disturbance radius (in meters)
  }
  rm(hab_params)
  
  # overlap discount
  nav_rad_ref <- which.max(v_dist_bins>nav_rad) # first one larger than nav_rad
  overlap_discount <- 1/rowSums(patch_dists_ref<nav_rad_ref)
  # hist(overlap_discount,breaks=0.5*(1:(2*max(overlap_discount)+1)))
  
  #### ------------------Initialize population ---------------------------------
  n_init <- round(0.75*nrow(patch_locations))
  pop_df <- data.table(patch=sample(patch_locations$id,n_init,replace = FALSE),
                       theta=runif(n = n_init,min = theta_start_min,max = theta_start_max),
                       p=runif(n = n_init,min = p_start_min,max = p_start_max))
  pop_df$ancestor <- pop_df$patch # column to track which patch the individual's ancestor came from
  # calculate effective theta for all adults
  if(plasticity_on=="multiplicative"){
    pop_df$eff_theta <- pop_df$theta*exp(2*pop_df$p*(patch_locations$q[pop_df$patch]-0.5)) # p has a multiplicative effect
  } else if(plasticity_on=="additive"| plasticity_on==TRUE){
    pop_df$eff_theta <- pop_df$theta+pop_df$p*(patch_locations$q[pop_df$patch]-0.5) # p has an additive effect
  } else pop_df$eff_theta <- pop_df$theta
  
  
  #### ------------Save generic larval output data, if desired -----------------
  if(larval_output_by_theta==TRUE){
    df_thetas <- data.frame(theta=v_theta_bins,output_mean=NA,output_sd=NA)
    larvae_out <- matrix(nrow=nrow(patch_locations),ncol=nrow(pop_df))
    
    for(th_i in seq_along(v_theta_bins)){
      eff_th <- v_theta_bins[th_i]
      for(adult_i in 1:nrow(pop_df)){
        # get its origin patch and effective theta value
        origin_patch <- pop_df$patch[adult_i]
        # evaluate the kernel from parent to each destination patch
        # rows are destination patches (indexed in patch_locations); cols are parents (indexed in pop_df)
        larvae_out[,adult_i] <- overlap_discount*offmap_corrections[adult_i,th_i]*ref_mat[patch_dists_ref[,origin_patch],th_i]
      } # adult_i
      outs <- colSums(larvae_out)
      df_thetas$output_mean[th_i] <- mean(outs)
      df_thetas$output_sd[th_i] <- sd(outs)
    }
    save(df_thetas,file=paste0(output_file,"_larval_output_by_theta.RData"))
  }
  
  #### ----------------------------Simulation ----------------------------------
  print("start simulation")
  print(proc.time()-starttime)
  interval_starttime <- proc.time()
  
  for(t_i in 1:nsteps){
    # check the population isn't zero
    if(nrow(pop_df)==0){
      print(paste0("Timestep ",t_i,": Population = 0"))
      break
    }
    
    # create matrix to hold larvae
    larvae_out <- matrix(nrow=nrow(patch_locations),ncol=nrow(pop_df))
    
    ##### ---------------------Disturbance -------------------------------------
    # In anemones affected by disturbance, habitat quality is reduced (affects fecundity and plastic dispersal response)
    
    # check for an active disturbance. If there is one:
    #   if its time hasn't run out, keep it going
    #   if its time has run out, end it
    if(disturb_length>0){
      disturb_length <- disturb_length+1
      if(disturb_length>Dl){  # end disturbance and reset
        disturb_length <- 0
        disturbed_anems <- c()
        patch_locations$q <- patch_locations$q_original
        patch_locations$b <- round(base_fecund*patch_locations$q)
      }
    }
    
    # if there's not an active disturbance already, one occurs with probability Dp
    # disturbed_anems is a vector, indexed as rows of patch_locations, with 1 = UNDISTURBED, 0 = DISTURBE
    # switching this around. Now, disturbed_anems is the indices of disturbed anemones in patch_locations.
    if(disturb_length==0 & rbinom(1,1,Dp)==TRUE){  
      print(paste0("Disturbance: ",t_i))
      fwrite(data.frame(t_i=t_i),file=paste0(output_file,"_disturb.csv"),append=TRUE) # record start time of new disturbance
      disturb_length <- 1
      # identify anemones affected by disturbance
      if(disturb_method=="fractal"){ ## fracland method (quite slow for maps of large extent)
        disturb_map <- fracland(k=k,h=1.5,p=De,binary=FALSE,plotflag=FALSE) 
        disturb_map <- disturb_map[1:ny,1:nx] 
        disturb_rast <- rast(ext(base_rast), resolution=res(base_rast), crs = crs(base_rast))
        values(disturb_rast) <- f_TransformDist(disturb_map,"A")<De
        disturbed_anems <- which(terra::extract(disturb_rast,vect(sfc_patches),xy=TRUE,search_radius=500)$lyr.1==1)
      } else{ ## circle method
        disturbance_center <- c(runif(1,xmin,xmax),runif(1,ymin,ymax))
        disturb_dists <- sqrt((patch_locations$x-disturbance_center[1])^2+(patch_locations$y-disturbance_center[2])^2)
        disturbed_anems <- which(disturb_dists<Dr)
      }
      # reduce disturbed anemones' habitat quality and fecundity
      patch_locations$q[disturbed_anems] <- Dm*patch_locations$q[disturbed_anems]
      patch_locations$b[disturbed_anems] <- Dm*patch_locations$b[disturbed_anems]
    }
    
    #### ------------------------ Plasticity -----------------------------------
    if(plasticity_on=="multiplicative"){
      pop_df$eff_theta <- pop_df$theta*exp(2*pop_df$p*(patch_locations$q[pop_df$patch]-0.5)) # p has a multiplicative effect
    } else if(plasticity_on=="additive"| plasticity_on==TRUE){
      pop_df$eff_theta <- pop_df$theta+pop_df$p*(patch_locations$q[pop_df$patch]-0.5) # p has an additive effect
    } else pop_df$eff_theta <- pop_df$theta
    
    #### ------------------ Output adult trait values --------------------------
    if(t_i %% output_thin == 0){
      # p
      p_quants <- quantile(pop_df$p, probs=c(0.05,0.25,0.5,0.75,0.95))
      df_p <- data.frame(t_i=t_i,metric="p",mean=mean(pop_df$p),var=var(pop_df$p),
                         q05=p_quants[1],q25=p_quants[2],median=p_quants[3],q75=p_quants[4],q95=p_quants[5],
                         corr_q=cor(pop_df$p,patch_locations$q[pop_df$patch]))
      fwrite(df_p,file=paste0(output_file,"_summary.csv"),append=TRUE)
      
      # theta (fundamental)
      theta_quants <- quantile(pop_df$theta, probs=c(0.05,0.25,0.5,0.75,0.95))
      df_theta <- data.frame(t_i=t_i,metric="theta",mean=mean(pop_df$theta),var=var(pop_df$theta),
                             q05=theta_quants[1],q25=theta_quants[2],median=theta_quants[3],q75=theta_quants[4],q95=theta_quants[5],
                             corr_q=cor(pop_df$theta,patch_locations$q[pop_df$patch]))
      fwrite(df_theta,file=paste0(output_file,"_summary.csv"),append=TRUE)
      
      efftheta_quants <- quantile(pop_df$eff_theta, probs=c(0.05,0.25,0.5,0.75,0.95))
      df_efftheta <- data.frame(t_i=t_i,metric="efftheta",mean=mean(pop_df$eff_theta),var=var(pop_df$eff_theta),
                                q05=efftheta_quants[1],q25=efftheta_quants[2],median=efftheta_quants[3],q75=efftheta_quants[4],q95=efftheta_quants[5],
                                corr_q=cor(pop_df$eff_theta,patch_locations$q[pop_df$patch]))
      fwrite(df_efftheta,file=paste0(output_file,"_summary.csv"),append=TRUE)
    }
    
    ### --------------------------- Dispersal ----------------------------------
    for(adult_i in 1:nrow(pop_df)){
      # get its origin patch and effective theta value
      origin_patch <- pop_df$patch[adult_i]
      #eff_th <- min(pop_df$eff_theta[adult_i],max(v_theta_bins))
      eff_th <- max(min(pop_df$eff_theta[adult_i],theta_max),theta_min)
      #comp_results$eff_theta <- pmax(pmin(comp_results$eff_theta,theta_max),theta_min)
      
      # find the closest values in v_theta_bins, and their corresponding weights
      u_th_ind <- max(match(TRUE,eff_th<v_theta_bins),2) # set to smallest bin if it's smaller than that
      if(is.na(u_th_ind)) u_th_ind <- length(v_theta_bins) # set to largest bin if it's larger than that
      l_th_ind <- u_th_ind-1
      l_th <- v_theta_bins[l_th_ind]
      u_th <- v_theta_bins[u_th_ind]
      wt_upper <- max((eff_th - l_th)/(u_th - l_th),0)
      wt_lower <- max(1-wt_upper,0)
      
      # evaluate the kernel from parent to each destination patch
      # rows are destination patches (indexed in patch_locations); cols are parents (indexed in pop_df)
      larvae_out[,adult_i] <- overlap_discount*(wt_lower*offmap_corrections[adult_i,l_th_ind]+wt_upper*offmap_corrections[adult_i,u_th_ind])*
        (wt_lower*ref_mat[patch_dists_ref[,origin_patch],l_th_ind]+wt_upper*ref_mat[patch_dists_ref[,origin_patch],u_th_ind])
    } # adult_i
    
    # normalize larvae_out so each parent produces an equal number of offspring (before accounting for differences in fecundity).
    # this removes any differences based on location or kernel
    if(normalize_offspring==TRUE){
      larvae_out <- larvae_out/rep(colSums(larvae_out),each=nrow(larvae_out))
    }
    
    #### --------------------------- Reproduction ------------------------------
    # multiply number of larvae from each origin patch by that patch's fecundity (disturbed patches have zero fecundity)
    larvae_out <- larvae_out*matrix(rep(patch_locations$b[pop_df$patch],each=nrow(larvae_out)),nrow=nrow(larvae_out))
    
    #### -------------------------- Adult mortality ----------------------------
    # decide which adults will survive to the next generation: 1 for survival, 0 for death
    # INDEXED AS ROWS OF pop_df!
    adults_survive <- rbinom(nrow(pop_df),1,adult_survival_prob)
    # and decide which anemones will have vacancies: 1 for vacancy, 0 for no vacancy
    # INDEXED AS ROWS OF patch_locations!
    will_have_vacancy <- rep(1,length=nrow(patch_locations))
    will_have_vacancy[pop_df$patch] <- !adults_survive
    will_have_vacancy <- which(will_have_vacancy==1)
    
    #### --------------------------- Competition -------------------------------
    patch_abunds <- rowSums(larvae_out)
    nadults=ncol(larvae_out)
    
    comp_results <- mclapply(will_have_vacancy,function(i_patch){
      i_abund=patch_abunds[i_patch]
      K_i=patch_locations$K[i_patch]
      temp_pop_i=larvae_out[i_patch,]
      # decide the maximum number of larvae that can settle in each patch, based on how many arrived (patch_abunds[i_patch])
      if(i_abund>0 & i_abund<1){
        n_setts <- rbinom(1,1,prob=i_abund)
      } else if(i_abund>=1){
        n_setts <- round(i_abund)
      } else n_setts=0
      
      # do the competition
      if(n_setts>0){ # if the patch isn't empty
        survivors=rmultinom(n=1,
                            size=min(n_setts,K_i), # choose groups for min(n_setts, K) settlers
                            prob = temp_pop_i) # probability of each group being chosen depends on its current abundance
        survivors <- which(survivors>0) #### NEED TO CHANGE THIS IF WE EVER HAVE K>1!
        # survivors <- sample.int(length(temp_pop_i),size = min(n_setts,K_i),prob=temp_pop_i) # this performs worse, even though it's simpler
        new_pop_i <- cbind(as.vector(survivors),rep(i_patch,length(survivors)))
      } else new_pop_i <- NULL
      
      return(new_pop_i)
    },mc.cores=numCores)
    
    comp_results <- do.call(rbind,comp_results)
    comp_results <- as.data.table(comp_results)
    if(nrow(comp_results)==0) comp_results <- data.table(parent=integer(),patch=integer())
    colnames(comp_results) <- c("parent","patch")
    
    #### ------------------ Mutation (and trait inheritance) -------------------
    if(mutation_type=="rare"){
      # # in this version, mu_theta and mu_p are mutation probabilities. The magnitudes of mutations are normally-distributed with mean=0 and sd=1. 
      theta_muts <- rbinom(nrow(comp_results),1,prob=mu_theta)
      theta_muts[theta_muts==1] <- rnorm(n=sum(theta_muts),mean=0,sd=1)
      p_muts <- rbinom(nrow(comp_results),1,prob=mu)
      p_muts[p_muts==1] <- rnorm(n=sum(p_muts),mean=0,sd=1)
      comp_results$theta <- pmax(pmin(pop_df$theta[comp_results$parent] + theta_muts,
                                      theta_max),theta_min) # cap theta
      comp_results$p <- pop_df$p[comp_results$parent] + p_muts
      
    } else{
      # # in this version, mu_theta and mu_p are mutation magnitudes. The probability of mutation for each individual is 1. 
      comp_results$theta <- pmax(pmin(pop_df$theta[comp_results$parent]+rnorm(nrow(comp_results),mean=0,sd=sqrt(mu_theta)),
                                      theta_max),theta_min) # cap theta
      comp_results$p <- pop_df$p[comp_results$parent]+ rnorm(nrow(comp_results),mean=0,sd=sqrt(mu))
    }
    comp_results$ancestor <- pop_df$ancestor[comp_results$parent]
    comp_results$eff_theta <- NA
    pop_df <- rbind(pop_df[adults_survive,c("patch","theta","p","eff_theta","ancestor")],comp_results[,c("patch","theta","p","eff_theta","ancestor")])
    
    #### --------------------------------- Output ------------------------------
    if(t_i %% 100 == 0){
      print(t_i)
      print(proc.time()-interval_starttime)
      interval_starttime <- proc.time()  
    }
    
    if(output_flag=="all" & (t_i<=25 | (t_i %% output_thin) == 0)){
      # if(t_i %in% c(seq(from=300, to=2000, by=300),seq(from=2000,to=3000,by=100))){
      export_mat <- pop_df
      export_mat$t_i <- t_i
      export_mat$q <- patch_locations$q[export_mat$patch]
      export_mat$disturb <- 0
      export_mat$disturb[export_mat$patch %in% disturbed_anems] <- 1
      
      fwrite(export_mat,file=paste0(output_file,"_all.csv"),append=TRUE)
    }
    
    if(t_i %% 100 == 0){
      save(pop_df,t_i,patch_locations,disturb_length,file=paste0(output_file,"_popsnapshot.RData"))
      # save(comp_results,t_i,file=paste0(output_file,"_comp_results.RData"))
    }
    
    # adult and larval abundance
    if(t_i %% output_thin == 0){
      if(nrow(pop_df)>0){
        # adult abundance
        df_adult <- data.frame(t_i=t_i,metric="abund",mean=NA,var=NA,q05=NA,q25=NA,median=nrow(pop_df),q75=NA,q95=NA,corr_q=NA)
        fwrite(df_adult,file=paste0(output_file,"_summary.csv"),append=TRUE)
        
        # number of larvae produced this timestep
        df_larv <- data.frame(t_i=t_i,metric="larval_abund",mean=NA,var=NA,
                              q05=NA,q25=NA,median=sum(patch_abunds),q75=NA,q95=NA,
                              corr_q=NA)
        fwrite(df_larv,file=paste0(output_file,"_summary.csv"),append=TRUE)
      }
    }
  } # t_i
  
  return(simID)
}