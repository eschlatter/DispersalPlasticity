.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
library(data.table)
library(parallel)

#### User inputs ####
habID <- 17
repID <- 1
output_flag="lite"
output_thin=1
experiment_folder <- "New"  #name of directory where everything is stored
seed_i <- NULL   #for re-running sims
notes <- NA   #can include a character string here with notes on the sim

# bio parameters
mu <- 0.001
nav_rad <- 0.05
adult_survival_prob <- 0
base_fecund <- 20

# sim parameters
theta_start_min <- 0.005
theta_start_max <- 0.005
p_start_min <- 0
p_start_max <- 0
nsteps <- 6000

# sim options
normalize_offspring=TRUE
plasticity_on=FALSE

#### Initial outputs ####
# All params
df_params <- data.frame(mu=mu,nav_rad=nav_rad,adult_survival_prob=adult_survival_prob,
                        base_fecund=base_fecund,
                        theta_start_min=theta_start_min,theta_start_max=theta_start_max,
                        p_start_min=p_start_min,p_start_max=p_start_max,
                        nsteps=nsteps,output_thin=output_thin,
                        normalize_offspring=normalize_offspring,plasticity_on=plasticity_on)
# Sim metadata
# generate simID
while(exists("simID")==FALSE){
  simID <- paste0(sprintf("%08d", sample(10000000:99999999,size=1)))
  if(simID %in% read.csv(paste0(experiment_folder,"/_index_sims.csv"))$sim_id) rm(simID) # check that it hasn't been used already
}
simDate <- format(Sys.time(),"%Y_%m_%d %H:%M:%S")
slurmJob <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
slurmTask <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(slurmTask)) slurmTask <- 1
df_index <- data.frame(sim_id=paste(simID),simDate=simDate,slurmJob=paste0(slurmJob,"_",slurmTask),
                       rep_id=repID,notes=notes)
# Habitat info
hab_row <- read.csv(paste0(experiment_folder,"/Maps/index_habs.csv")) |> 
  dplyr::filter(hab_id==habID)
# Write all sim info (params, metadata, and habitat) to index
index_sims <- cbind(df_index,hab_row,df_params)
fwrite(index_sims,file=paste0(experiment_folder,"/_index_sims.csv"),append = TRUE)

# output file name
output_file <- paste0(experiment_folder,"/output/",simID)

# set and save seed
RNGkind("L'Ecuyer-CMRG")
if(!exists("seed_i")) seed_i <- .Random.seed
set.seed(seed_i)
save(seed_i,file=paste0(output_file,"_seed.RData"))

#### Set up objects ####
starttime <- proc.time()
numCores=max(1,parallelly::availableCores())

# patch_locations
load(paste0(experiment_folder,"/Maps/habfiles/hab_",habID,".RData"))
patch_locations <- hab_params$patch_locations
patch_locations$b <- round(base_fecund*patch_locations$q)
rm(hab_params)

# kernel lookup table
load(paste0(experiment_folder,"/Maps/b",index_sims$basemap_id,"/RefMat_b",index_sims$basemap_id,
            "_p",index_sims$popmap_id,".RData"))
theta_max <- max(v_theta_bins)
theta_min <- min(v_theta_bins)

# multiply the kernel by pi*nav_rad^2 once here, because it's the same throughout the sim
# (and replace the self-recruitment values with real integrals, if we're doing that)
if(file.exists(paste0(experiment_folder,"/Maps/b",index_sims$basemap_id,"/RefMat_selfrecruit_b",index_sims$basemap_id,
                      "_p",index_sims$popmap_id,"_navrad",nav_rad,".RData"))){
  load(paste0(experiment_folder,"/Maps/b",index_sims$basemap_id,"/RefMat_selfrecruit_b",index_sims$basemap_id,
              "_p",index_sims$popmap_id,"_navrad",nav_rad,".RData"))
} else{
  ref_mat <- pi*nav_rad^2*ref_mat
  for(i in seq_along(v_theta_bins)){
    th_i <- v_theta_bins[i]
    ref_mat[1,i] <- calculus::integral(function(r,theta){(1/th_i)*exp(-r/th_i)},
                                       bounds=list(r=c(0,nav_rad),theta=c(0,2*pi)),
                                       coordinates="polar",
                                       relTol=0.00001,
                                       method="divonne"
    )$value
  }
  save(ref_mat,file=paste0(experiment_folder,"/Maps/b",index_sims$basemap_id,"/RefMat_selfrecruit_b",index_sims$basemap_id,
                           "_p",index_sims$popmap_id,"_navrad",nav_rad,".RData"))
}


# overlap discount
nav_rad_ref <- which.max(v_dist_bins>nav_rad) # first one larger than nav_rad
overlap_discount <- 1/rowSums(patch_dists_ref<nav_rad_ref)
# hist(overlap_discount,breaks=0.5*(1:(2*max(overlap_discount)+1)))

# initialize population
n_init <- round(0.75*nrow(patch_locations))
pop_df <- data.table(patch=sample(patch_locations$id,n_init,replace = FALSE),
                     theta=runif(n = n_init,min = theta_start_min,max = theta_start_max),
                     p=runif(n = n_init,min = p_start_min,max = p_start_max))
# calculate effective theta for all adults
if(plasticity_on==TRUE){pop_df$eff_theta <- pop_df$theta+pop_df$p*patch_locations$q[pop_df$patch]
} else pop_df$eff_theta <- pop_df$theta

#### Simulation ####
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
  
  ##### Disturbance ####
  # decide which adults will survive to the next generation
  # INDEXED AS ROWS OF pop_df!
  adults_survive <- sample(1:nrow(pop_df),size = adult_survival_prob*nrow(pop_df),replace = FALSE) # rows of pop_df to keep for next time
  
  # and decide which anemones will have vacancies
  # (currently this depends precisely on who survives above, but it could be more complicated)
  # INDEXED AS ROWS OF patch_locations!
  will_have_vacancy <- which(patch_locations$id %notin% pop_df$patch[adults_survive])
  
  # # For checking if the number of offspring produced varies based on theta
  # df_thetas <- data.frame(theta=v_theta_bins,output_mean=NA,output_sd=NA)
  # 
  # for(th_i in seq_along(v_theta_bins)){
  #   eff_th <- v_theta_bins[th_i]
  #   # find the closest values in v_theta_bins, and their corresponding weights
  #   u_th_ind <- max(match(TRUE,eff_th<v_theta_bins),2) # set to smallest bin if it's smaller than that
  #   if(is.na(u_th_ind)) u_th_ind <- length(v_theta_bins) # set to largest bin if it's larger than that
  #   l_th_ind <- u_th_ind-1
  #   l_th <- v_theta_bins[l_th_ind]
  #   u_th <- v_theta_bins[u_th_ind]
  #   wt_upper <- max((eff_th - l_th)/(u_th - l_th),0)
  #   wt_lower <- max(1-wt_upper,0)
  # 
  #   for(adult_i in 1:nrow(pop_df)){
  #     # get its origin patch and effective theta value
  #     origin_patch <- pop_df$patch[adult_i]
  # 
  #     # evaluate the kernel from parent to each destination patch
  #     # rows are destination patches (indexed in patch_locations); cols are parents (indexed in pop_df)
  #     larvae_out[,adult_i] <- overlap_discount*(pi*nav_rad^2)*(wt_lower*offmap_correction[l_th_ind]+wt_upper*offmap_correction[u_th_ind])*
  #       (wt_lower*ref_mat[patch_dists_ref[,origin_patch],l_th_ind]+wt_upper*ref_mat[patch_dists_ref[,origin_patch],u_th_ind])
  #   } # adult_i
  # 
  #   selfrecruit_inds <- which(patch_dists_ref[,pop_df$patch]==1)
  #   larvae_out[selfrecruit_inds] <- larvae_out[selfrecruit_inds]/(pi*nav_rad^2)
  # 
  #   outs <- colSums(larvae_out)
  #   df_thetas$output_mean[th_i] <- mean(outs)
  #   df_thetas$output_sd[th_i] <- sd(outs)
  # }
  # 
  # plot(df_thetas$theta,df_thetas$output_mean,type='o-')
  
  ### Dispersal ####
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
    larvae_out[,adult_i] <- overlap_discount*(wt_lower*offmap_correction[l_th_ind]+wt_upper*offmap_correction[u_th_ind])*
      (wt_lower*ref_mat[patch_dists_ref[,origin_patch],l_th_ind]+wt_upper*ref_mat[patch_dists_ref[,origin_patch],u_th_ind])
  } # adult_i
  
  # normalize larvae_out so each parent produces an equal number of offspring (before accounting for differences in fecundity).
  # this removes any differences based on location or kernel
  if(normalize_offspring==TRUE){
    larvae_out <- larvae_out/rep(colSums(larvae_out),each=nrow(larvae_out))
  }
  
  #### Reproduction ####
  # multiply number of larvae from each origin patch by that patch's fecundity
  larvae_out <- larvae_out*matrix(rep(patch_locations$b[pop_df$patch],each=nrow(larvae_out)),nrow=nrow(larvae_out))
  
  #### Competition ####
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
  #comp_results <- comp_results[comp_results[,1]>0,]
  comp_results <- as.data.table(comp_results)
  if(nrow(comp_results)==0) comp_results <- data.table(parent=integer(),patch=integer())
  colnames(comp_results) <- c("parent","patch")
  
  #### Mutation ####
  comp_results$theta <- pmax(pmin(pop_df$theta[comp_results$parent]+rnorm(nrow(comp_results),mean=0,sd=sqrt(mu)),
                                theta_max),theta_min) # cap theta
#  comp_results$theta <- pop_df$theta[comp_results$parent] + rnorm(nrow(comp_results),mean=0,sd=sqrt(mu))
  comp_results$p <- pop_df$p[comp_results$parent]+ rnorm(nrow(comp_results),mean=0,sd=sqrt(mu))
  
  #### Plasticity ####
  if(plasticity_on==TRUE){comp_results$eff_theta <- comp_results$theta+comp_results$p*(patch_locations$q[comp_results$patch]-0.5)
  } else comp_results$eff_theta <- comp_results$theta
  
  pop_df <- rbind(pop_df[adults_survive,c("patch","theta","p","eff_theta")],comp_results[,c("patch","theta","p","eff_theta")])
  
  #### Output ####
  print(t_i)
  print(proc.time()-interval_starttime)
  interval_starttime <- proc.time()
  
  if(t_i %% output_thin == 0){
    export_mat <- pop_df
    export_mat$t_i <- t_i
    fwrite(export_mat,file=paste0(output_file,"_raw.csv"),append=TRUE)
  }
  
  if(t_i %% 100 == 0){
    save(pop_df,t_i,file=paste0(output_file,"_popsnapshot.RData"))
  }
  
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
      
      fwrite(df_i,file=paste0(output_file,"_summary.csv"),append=TRUE)
    }
  }
  
  
} # t_i
# 
#### Process output ####
library(dplyr)
library(ggplot2)
dat_out <- read.csv(paste0(output_file,"_summary.csv")) %>%
  filter(metric %in% c("p","theta","efftheta"))

abund_out <- read.csv(paste0(output_file,"_summary.csv")) %>%
  filter(metric %in% c("abund","larval_abund"))

if(exists("v_theta_bins")){
  linemin=min(v_theta_bins)
  linemax=max(v_theta_bins)
} else{
  linemin=0
  linemax=3.5
}

g1 <- ggplot(filter(dat_out,metric=="theta"),aes(x=t_i,y=median))+
  geom_line()+
  geom_ribbon(aes(ymin=q05,ymax=q95),alpha=0.2)+
  labs(y="Theta (km)\n(median,5-95%)",x=NULL)+
  geom_hline(yintercept=c(linemin,linemax),linetype='dashed')

g2 <- ggplot(filter(dat_out,metric=="p"),aes(x=t_i,y=median))+
  geom_line()+
  geom_ribbon(aes(ymin=q05,ymax=q95),alpha=0.2)+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y="Plasticity\n(median,5-95%)",x=NULL)

g3 <- ggplot(filter(abund_out,metric=="abund"),aes(x=t_i,y=median))+
  geom_line()+
  labs(y="Adult\nabundance")

gridExtra::grid.arrange(g1,g2,g3,ncol=1)

# ggplot(abund_out,aes(x=t_i,y=median))+geom_line()+
#   facet_wrap(vars(metric),nrow=2,scales="free", strip.position = "left")

all_out <- read.csv(paste0(output_file,"_raw.csv"))
all_out_last <- filter(all_out,t_i==max(all_out$t_i))
