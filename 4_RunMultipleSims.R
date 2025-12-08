## Check if spatial variation in the (effective) kernel mean is the same across stochastic realizations on the same seascape

source('0_Setup.R')
load('params/VariableKPars.RData')
params$v_p <- -2:2
params$p_start <- 3
params$nsteps=1000
params$mu=0.01
list2env(x=params,envir=environment())

nrep=3
kernmeans_df <- data.frame(all_kernmean=numeric(),island_kernmean=numeric(),mainland_kernmean=numeric())
for(i in 1:nrep){
  sim_loop1 <- f_RunMatrixLoop(params)
  sim_loop_out1 <- f_ProcessLoopOutputDataTable(sim_loop1$params,sim_loop1$Pij,sim_loop1$group_index, sim_loop1$time_run)
  by_t <- sim_loop_out1$by_t
  sim_melt <- sim_loop_out1$sim_melt
  
  kern_timesteps=seq(from=round(0.75*nrow(by_t)), to=nrow(by_t),length.out=25)
  sim_df <- sim_melt[patch_locations,on=c(patch="id")] # join sim_melt and patch_locations
  sim_df <- sim_df[sim_df$popsize!=0,] # remove zero rows
  sim_df <- sim_df[(sim_df$t %in% kern_timesteps),] # choose only needed timesteps
  eff_pars <- f_plasticityK_new(K=sim_df$K,p=sim_df$p,alpha=sim_df$alpha,theta=sim_df$theta,n_alpha=length(v_alphas),n_theta=length(v_thetas),Kmin=min(patch_locations$K_i),Kmax=max(patch_locations$K_i))
  sim_df <- cbind(sim_df,eff_pars)
  sim_df$alpha_plastic_val <- v_alphas[sim_df$alpha_plastic]
  sim_df$theta_plastic_val <- v_thetas[sim_df$theta_plastic]
  
  sim_df_island <- sim_df[(sim_df$x %in% island_lims$xmin:island_lims$xmax) & (sim_df$y %in% island_lims$ymin:island_lims$ymax),]
  sim_df_mainland <- sim_df[(sim_df$x %in% mainland_lims$xmin:mainland_lims$xmax) & (sim_df$y %in% mainland_lims$ymin:mainland_lims$ymax),]
  
  # average effective parameter everywhere
  all_alpha <- sum(sim_df$popsize*sim_df$alpha_plastic_val)/sum(sim_df$popsize)
  all_theta <- sum(sim_df$popsize*sim_df$theta_plastic_val)/sum(sim_df$popsize)
  all_kernmean <- all_alpha*all_theta
  
  # average effective parameter island
  island_alpha <- sum(sim_df_island$popsize*sim_df_island$alpha_plastic_val)/sum(sim_df_island$popsize)
  island_theta <- sum(sim_df_island$popsize*sim_df_island$theta_plastic_val)/sum(sim_df_island$popsize)
  island_kernmean <- island_alpha*island_theta
  
  # average effective parameter mainland
  mainland_alpha <- sum(sim_df_mainland$popsize*sim_df_mainland$alpha_plastic_val)/sum(sim_df_mainland$popsize)
  mainland_theta <- sum(sim_df_mainland$popsize*sim_df_mainland$theta_plastic_val)/sum(sim_df_mainland$popsize)
  mainland_kernmean <- mainland_alpha*mainland_theta
  
  kernmeans_df <- rbind(kernmeans_df,data.frame(all_kernmean=all_kernmean,island_kernmean=island_kernmean,mainland_kernmean=mainland_kernmean))
}

kernmeans_df <- pivot_longer(kernmeans_df,everything())
ggplot(kernmeans_df)+
  geom_boxplot(aes(name,value))+
  theme_minimal()+
  labs(x='type of reef patch',y='kernel mean')