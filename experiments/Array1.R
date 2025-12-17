source('0_Setup.R')
load("experiments/Test_Array1_SimIndex.RData") #load the dataframe of habitats to use

### set up data structures
output_all_list=list()

### Run sims
run_i <- c(5,6,7)[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
#for(run_i in sim_index$Run[1:3]){
  print(paste("Run",run_i))
  
  load('params/ParSetArray1.RData')
  list2env(params,envir=environment())
  
  #load(paste0("seascapes/kmaps/kmap_",run_i,".RData"))
  load(paste0("seascapes/sizetest/kmap_",run_i,".RData"))
  list2env(kmap_i,envir=environment())
  rm(kmap_i)
 
  # g <- f_Plot_Landscape(patch_locations,nx,ny,do_now=FALSE)+
  #   labs(title=paste("Run",run_i))
  # print(g)
  
  nsteps <- 12
  params <- list(nx=nx,ny=ny,nsteps=nsteps,
                 v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
                 alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
                 mu=mu,b=b,K=K,
                 disturb_prob=disturb_prob,patch_locations=patch_locations)
  
#  sim_out=list(test=paste("this is run",run_i))
  sim_out <- f_RunMatrixLoopLite(params,keep=list("p","abund","kern", "sp_struct"))
#  output_all_list[[run_i]] <- sim_out #store everything for checking later
# } # run_i

save(sim_out,file=paste0('output/sizetest/',run_i,'_sizetest.RData'))
  #save(output_all_list,file='output/Test_Array1Out.RData')