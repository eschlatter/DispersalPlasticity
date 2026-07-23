#source('0_Setup.R')
.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
library(profvis)
library(foreach)
library(doParallel)
library(data.table)
#setDTthreads(1, restore_after_fork=FALSE)
#### PARAMETERS ####
mu=0.01
output_flag="lite"
output_thin=1
output_file="output/TryNew_1/TryNew_MultiRows"
nsteps=100
numCores=parallelly::availableCores()
print(numCores)
npatch=1000

prof_sim <- profvis({
  
  #### SET UP OBJECTS ####
  starttime <- proc.time()
  
  # patch_locations
  load("experiments/Exp4_20260611/habfiles/habparams_4.RData")
  patch_locations <- hab_params$patch_locations[1:npatch,]
  
  # get patch dists
  patchdist_file <- "experiments/Exp4_20260611/b3/patchdists_b3_p3.RData"
  load(patchdist_file)
  patch_dists <- collapse::qM(patch_dists) # remove units, but it's in km
  patch_dists <- patch_dists[1:npatch,1:npatch]
  
  # make a reference matrix of distance bins, theta values, where each site is the kernel evaluated at the distance
  # (this is fast and small)
  v_dist_bins <- c(seq(from=0,to=1,by=0.001),seq(from=1,to=max(patch_dists)+0.01,by=0.01))
  v_theta_bins <- c(seq(from=0.005,to=0.1,by=0.005),seq(from=0.1,to=max(patch_dists)/2,by=0.1))
  ref_mat <- matrix(nrow=length(v_dist_bins),ncol=length(v_theta_bins))
  for(th_i in seq_along(v_theta_bins)){
    ref_mat[,th_i] <- dgamma(v_dist_bins,shape=1,scale=v_theta_bins[th_i])
  }
  
  # figure out which row of ref_mat each distance corresponds to
  # (returns the lower bound distance)
  # this takes some time, and is the same size as patch_dists.
  patch_dists <- asplit(patch_dists,2) # split patch_dists into a list of columns
  gc()
  cluster <- makeForkCluster(nnodes=numCores)
  registerDoParallel(cluster)
  patch_dists_ref <- foreach(col_i = patch_dists,.combine="cbind") %dopar% {
    vapply(col_i,function(cell_i) match(TRUE,cell_i<v_dist_bins)-1,FUN.VALUE=double(1),USE.NAMES = FALSE)
    # gc()
  }
  stopCluster(cl = cluster)
  storage.mode(patch_dists_ref) <- "integer"
  colnames(patch_dists_ref) <- NULL
  rownames(patch_dists_ref) <- NULL
  save(patch_dists_ref,ref_mat,file=paste0(output_file,"_patch_dists_ref.RData"))
  rm(patch_dists)
  
  # initialize population
  pop_df <- data.frame(patch=sample(patch_locations$id,0.75*nrow(patch_locations),replace = FALSE),
                       theta=1,
                       p=0)
  # calculate effective theta for all adults
  pop_df$eff_theta <- pop_df$theta+pop_df$p*patch_locations$q[pop_df$patch]
  
  #### SIMULATION ####
  print("start simulation")
  print(proc.time()-starttime)
  interval_starttime <- proc.time()
  for(t_i in 1:nsteps){
    # create matrix to hold larvae
    larvae_out <- matrix(nrow=nrow(patch_locations),ncol=nrow(pop_df))
    
    ##### DISTURBANCE ####
    # decide which adults will survive to the next generation
    # INDEXED AS ROWS OF pop_df!
    adults_survive <- sample(1:nrow(pop_df),size = 0.5*nrow(pop_df),replace = FALSE) # rows of pop_df to keep for next time
    
    # and decide which anemones will have vacancies
    # (currently this depends precisely on who survives above, but it could be more complicated)
    # INDEXED AS ROWS OF patch_locations!
    will_have_vacancy <- which(patch_locations$id %notin% pop_df$patch[adults_survive])
    
    #### REPRODUCTION ####
    # for each adult (in parallel):
    for(adult_i in 1:nrow(pop_df)){
      pop_df_i <- pop_df[adult_i,]
      # get its origin patch
      origin_patch <- pop_df_i$patch
      
      # grab its theta value (eff_th)
      eff_th <- pop_df_i$eff_theta
      # find the closest values in v_theta_bins, and their corresponding weights
      u_th_ind <- max(match(TRUE,eff_th<v_theta_bins),2) # set to smallest bin if it's smaller than that
      if(is.na(u_th_ind)) u_th_ind <- length(v_theta_bins) # set to largest bin if it's larger than that
      l_th_ind <- u_th_ind-1
      l_th <- v_theta_bins[l_th_ind]
      u_th <- v_theta_bins[u_th_ind]
      wt_upper <- (eff_th - l_th)/(u_th - l_th)
      wt_lower <- 1-wt_upper
      
      # evaluate the kernel from parent to each destination patch
      # rows are destination patches (indexed in patch_locations); cols are parents (indexed in pop_df)
      larvae_out[,adult_i] <- pmax(patch_locations$b[pop_df_i$patch]*(wt_lower*ref_mat[patch_dists_ref[,origin_patch],l_th_ind]+
                                                                   wt_upper*ref_mat[patch_dists_ref[,origin_patch],u_th_ind]),0)
    } # adult_i
    
    #### COMPETITION ####
    patch_abunds <- rowSums(larvae_out)
    nadults=ncol(larvae_out)
    
    # for anemones with a vacancy (as indexed in will_have_vacancy), decide which larva survives
    cluster <- makeForkCluster(nnodes=numCores)
    registerDoParallel(cluster)
    
    comp_results <- vector("list",length(will_have_vacancy))
    comp_results <- foreach(i_patch = will_have_vacancy) %dopar% {
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
        survivors=t(rmultinom(n=1,
                              size=min(n_setts,K_i), # choose groups for min(n_setts, K) settlers
                              prob = temp_pop_i)) # probability of each group being chosen depends on its current abundance
        survivors <- which(survivors>0) #### NEED TO CHANGE THIS IF WE EVER HAVE K>1!
        new_pop_i <- data.frame(parent=survivors, patch=i_patch)
      } else new_pop_i <- NULL

      # if(n_setts>0){ # if the patch isn't empty
      #   survivors=sample(1:nadults,size=min(n_setts,K_i),replace = TRUE, prob = temp_pop_i)
      #   new_pop_i <- data.frame(parent=survivors, patch=i_patch)
      # } else new_pop_i <- NULL
      
      comp_results[[i_patch]] <- new_pop_i
    }
    stopCluster(cl = cluster)
    comp_results <- rbindlist(comp_results)
    
    
    # comp_results <- mclapply(will_have_vacancy,
    #                                  function(ii) f_CompetitionNew(i_patch=ii,i_abund=patch_abunds[ii],K_i=patch_locations$K[ii],
    #                                                               temp_pop_i=larvae_out[ii,],nadults=ncol(larvae_out)),
    #                          mc.cores=numCores)
    # comp_results <- rbindlist(comp_results)
    
    #### MUTATION ####
    comp_results$theta <- pop_df$theta[comp_results$parent] + rnorm(nrow(comp_results),mean=0,sd=mu)
    comp_results$p <- pop_df$p[comp_results$parent]+ rnorm(nrow(comp_results),mean=0,sd=mu)
    
    #### PLASTICITY #### 
    comp_results$eff_theta <- comp_results$theta+comp_results$p*patch_locations$q[comp_results$patch]
    
    pop_df <- rbind(pop_df[adults_survive,],comp_results[,c("patch","theta","p","eff_theta")])
    
    
    #### OUTPUT ####
    print(t_i)
    print(proc.time()-interval_starttime)
    interval_starttime <- proc.time()
    
    if((output_flag=="all") & (t_i %% output_thin == 0)){
      export_mat <- pop_df
      export_mat$t_i <- t_i
      fwrite(export_mat,file=paste0(output_file,"_raw.csv"),append=TRUE)
    }

    if(t_i %% 100 == 0){
      save(pop_df,t_i,file=paste0(output_file,"_popsnapshot.RData"))
    }

    if(output_flag=="lite"){
      if(t_i %% output_thin == 0){
        if(nrow(pop_df)>0){

          df_i <- data.frame(t_i=rep(t_i,5),metric=NA,mean=NA,var=NA,q05=NA,q25=NA,median=NA,q75=NA,q95=NA,corr_q=NA)

          # adult abundance
          df_i[1,] <- data.frame(t_i=t_i,metric="abund",mean=NA,var=NA,q05=NA,q25=NA,median=nrow(pop_df),q75=NA,q95=NA,corr_q=NA)
          # p
          p_quants <- quantile(pop_df$p, probs=c(0.05,0.25,0.5,0.75,0.95))
          df_i[2,] <- data.frame(t_i=t_i,metric="p",mean=mean(pop_df$p),var=var(pop_df$p),
                                 q05=p_quants[1],q25=p_quants[2],median=p_quants[3],q75=p_quants[4],q95=p_quants[5],
                                 corr_q=cor(pop_df$p,patch_locations$q[pop_df$patch]))
          # theta (fundamental)
          theta_quants <- quantile(pop_df$theta, probs=c(0.05,0.25,0.5,0.75,0.95))
          df_i[3,] <- data.frame(t_i=t_i,metric="theta",mean=mean(pop_df$theta),var=var(pop_df$theta),
                                 q05=theta_quants[1],q25=theta_quants[2],median=theta_quants[3],q75=theta_quants[4],q95=theta_quants[5],
                                 corr_q=cor(pop_df$theta,patch_locations$q[pop_df$patch]))
          # theta (effective)
          efftheta_quants <- quantile(pop_df$eff_theta, probs=c(0.05,0.25,0.5,0.75,0.95))
          df_i[4,] <- data.frame(t_i=t_i,metric="efftheta",mean=mean(pop_df$eff_theta),var=var(pop_df$eff_theta),
                                 q05=efftheta_quants[1],q25=efftheta_quants[2],median=efftheta_quants[3],q75=efftheta_quants[4],q95=efftheta_quants[5],
                                 corr_q=cor(pop_df$eff_theta,patch_locations$q[pop_df$patch]))
          # number of larvae produced this timestep
          df_i[5,] <- data.frame(t_i=t_i,metric="larval_abund",mean=NA,var=NA,
                                 q05=NA,q25=NA,median=sum(patch_abunds),q75=NA,q95=NA,
                                 corr_q=NA)

          fwrite(df_i,file=paste0(output_file,"_raw.csv"),append=TRUE)
        }
      }
    } # if output_flag=="lite"

    
  } # t_i
  
  
}) # profvis
#
 save(prof_sim,file=paste0(output_file,"_profile.RData"))
 # 
# # load(paste0(output_file,"_profile.RData"))
# htmlwidgets::saveWidget(prof_sim,"profile_mono.html")
# browseURL("profile_mono.html")

# #### PROCESS OUTPUT ####
# dat_out <- read.csv(paste0(output_file,"_raw.csv")) %>%
#   filter(metric %in% c("p","theta","efftheta"))
# 
# abund_out <- read.csv(paste0(output_file,"_raw.csv")) %>%
#   filter(metric %in% c("abund","larval_abund"))
# 
# 
# ggplot(dat_out,aes(x=t_i,y=median))+
#   geom_line()+
#   facet_wrap(vars(metric),nrow=3,scales="free", strip.position = "left")
# 
