source('0_Setup.R')
## set up runs
# h_i=1.0 # value of h (aggregation level of K's) to use in this set of sims
# h can maybe meaningfully go from 0 to 1.8ish? I'm not sure. 

## 1. Seascape(s)
load('seascapes/basemap1.RData')
load('params/ParSet1.RData')
list2env(params,envir=environment())

## 2. Starting condition(s)
param_starts <- expand.grid(alpha_start=c(1,5),theta_start=c(1,5),p_start=c(1,5))
param_starts <- rbind(param_starts,c(0,0,0)) # distribute starting individuals among all parameter values

## 3. Habitat quality map(s)
v_h <- c(0.1*1:12)
nhabs=length(v_h)
habitats=list() # hold all the K-scapes generated
all_params <- param_starts[NULL,]
for(i in 1:nhabs){
  habitats[[i]] <- f_GenerateMapWithK(base_map=basemap_1,K_range=c(3,18),h=v_h[i],k=5,p=0.3,h_base=0.8,plot_flag=FALSE)
  new_params <- param_starts
  new_params$hab_num <- i
  all_params <- rbind(all_params,new_params)
}

# Run a bunch of sims
# --------------------------------------
### set up data structures
output_all_list=list()
output_all_df <- data.frame(p_mean=numeric(),p_var_among_t=numeric(),p_var_within_t=numeric(),
                            fund_kernmean_mean=numeric(),fund_kernmean_var_among_t=numeric(),fund_kernmean_var_within_t=numeric(),
                            eff_kernmean_mean=numeric(),eff_kernmean_var_among_t=numeric(),eff_kernmean_var_within_t=numeric(),
                            fund_kernmoran_mean=numeric(),fund_kernmoran_var_among_t=numeric(),
                            eff_kernmoran_mean=numeric(),eff_kernmoran_var_among_t=numeric())

for(hab_i in 1:nhabs){
  for(pars_i in 1:nrow(param_starts)){
    sim_num=(hab_i-1)*nrow(param_starts)+pars_i
    print(paste("Sim num:",sim_num,"hab:",hab_i))
    print(param_starts[pars_i,])
    
    load('seascapes/basemap1.RData')
    load('params/ParSet1.RData')
    list2env(params,envir=environment())
    list2env(habitats[[hab_i]],envir=environment())
    alpha_start <- param_starts$alpha_start[pars_i]
    theta_start <- param_starts$theta_start[pars_i]
    p_start <- param_starts$p_start[pars_i]
    nsteps <- 10000
    params <- list(nx=nx,ny=ny,nsteps=nsteps,
                   v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
                   alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
                   mu=mu,b=b,K=K,
                   disturb_prob=0,patch_locations=patch_locations)
    
    sim_out <- f_RunMatrixLoopLite(params,keep=list("p","abund","kern", "sp_struct"))
    output_all_list[[sim_num]] <- sim_out #store everything for checking later
    output_all_df <- rbind(output_all_df,f_GetOutputStats(sim_out$output_df,nsteps))
  } # pars_i
} # hab_i

output_all_df <- cbind(all_params,output_all_df)
output_all_df$h <- 0.1*output_all_df$hab_num
output_all_df$start_grp <- rep(1:nrow(param_starts),nhabs)

save(output_all_list,output_all_df,habitats,file="model_runs/variable_h.RData")

# Look at the output
# ----------------------------------------------

# check output from a single sim
f_PlotOutput_Lite(output_all_list[[75]])

# demonstrate gradient in aggregation
h1 <- f_Plot_Landscape(habitats[[1]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.1")
h2 <- f_Plot_Landscape(habitats[[2]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.2")
h3 <- f_Plot_Landscape(habitats[[3]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.3")
h4 <- f_Plot_Landscape(habitats[[4]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.4")
h5 <- f_Plot_Landscape(habitats[[5]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.5")
h6 <- f_Plot_Landscape(habitats[[6]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.6")
h7 <- f_Plot_Landscape(habitats[[7]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.7")
h8 <- f_Plot_Landscape(habitats[[8]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.8")
h9 <- f_Plot_Landscape(habitats[[9]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 0.9")
h10 <- f_Plot_Landscape(habitats[[10]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 1.0")
h11 <- f_Plot_Landscape(habitats[[11]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 1.1")
h12 <- f_Plot_Landscape(habitats[[12]]$patch_locations,nx,ny,do_now = FALSE)+
  labs(title="h = 1.2")
grid.arrange(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,nrow=3)

# h vs p
ggplot(output_all_df,aes(x=h,group=start_grp))+
  geom_point(aes(y=p_mean),position=position_dodge(width=0.05))+
  geom_errorbar(aes(ymin=p_mean-sqrt(p_var_among_t),ymax=p_mean+sqrt(p_var_among_t)),position=position_dodge(width=0.05))+
  labs(x='Habitat quality aggregation (h)',y="plasticity parameter (mean +/- sd over time)",title="Plasticity vs Habitat aggregation")+
  theme_minimal()

# h vs fundamental and effective kernels
##  mean
ggplot(output_all_df,aes(x=h,group=start_grp))+
  geom_point(aes(y=fund_kernmean_mean,color='Fundamental Kernel'),position=position_dodge(width=0.05))+
  geom_errorbar(aes(ymin=fund_kernmean_mean-sqrt(fund_kernmean_var_among_t),ymax=fund_kernmean_mean+sqrt(fund_kernmean_var_among_t),color='Fundamental Kernel'),
                position=position_dodge(width=0.05))+
  geom_point(aes(y=eff_kernmean_mean,color='Effective Kernel'),position=position_dodge(width=0.05))+
  geom_errorbar(aes(ymin=eff_kernmean_mean-sqrt(eff_kernmean_var_among_t),ymax=eff_kernmean_mean+sqrt(eff_kernmean_var_among_t),color='Effective Kernel'),
                position=position_dodge(width=0.05))+
  labs(x='Habitat quality aggregation (h)',y="dispersal kernel mean\n(mean +/- sd over time)",title="Dispersal kernel mean vs Habitat aggregation")+
  theme_minimal()
##  within-timestep variance
ggplot(output_all_df,aes(x=h,group=start_grp))+
  geom_point(aes(y=fund_kernmean_var_within_t,color='Fundamental Kernel'),position=position_dodge(width=0.05))+
  geom_point(aes(y=eff_kernmean_var_within_t,color='Effective Kernel'),position=position_dodge(width=0.05))+
  labs(x='Habitat quality aggregation (h)',y="Within-timestep variance in kernel mean",title="Dispersal kernel variation vs Habitat aggregation")+
  theme_minimal()

# h vs spatial structure in kernel
ggplot(output_all_df,aes(x=h,group=start_grp))+
  geom_point(aes(y=fund_kernmoran_mean,color='Fundamental Kernel'),position=position_dodge(width=0.05))+
  geom_errorbar(aes(ymin=fund_kernmoran_mean-sqrt(fund_kernmoran_var_among_t),ymax=fund_kernmoran_mean+sqrt(fund_kernmoran_var_among_t),color='Fundamental Kernel'),
                position=position_dodge(width=0.05))+
  geom_point(aes(y=eff_kernmoran_mean,color='Effective Kernel'),position=position_dodge(width=0.05))+
  geom_errorbar(aes(ymin=eff_kernmoran_mean-sqrt(eff_kernmoran_var_among_t),ymax=eff_kernmoran_mean+sqrt(eff_kernmoran_var_among_t),color='Effective Kernel'),
                position=position_dodge(width=0.05))+
  labs(x='Habitat quality aggregation (h)',y="Moran's I of dispersal kernel mean\n(mean +/- sd over time)",title="Dispersal kernel spatial structure vs Habitat aggregation")+
  theme_minimal()

# look at dynamics of a) effective kernel and b) plasticity from the various starting conditions in the same habitat
sim_nums <- which(output_all_df$h==0.9)
df_by_t_selected <- output_all_list[[1]]$output_df[NULL,]
df_by_t_selected$sim_num <- factor()
for(sim in sim_nums){
  sim_list <- output_all_list[[sim]]$output_df
  sim_list$sim_num <- sim
  df_by_t_selected <- rbind(df_by_t_selected,sim_list)
}
df_by_t_selected$sim_num <- as.factor(df_by_t_selected$sim_num)
df_by_t_selected <- df_by_t_selected[!is.na(df_by_t_selected$v_abund),]

ggplot(df_by_t_selected,aes(x=t))+
  geom_ribbon(aes(ymin=eff_mean_mean-sqrt(eff_mean_var),ymax=eff_mean_mean+sqrt(eff_mean_var),fill=sim_num),alpha=0.1)+
  geom_line(aes(y=eff_mean_mean,color=sim_num))+
  labs(x="time",y="Effective kernel mean\n(mean +/- sd within timesteps)",title="Simulation dynamics, h=0.9: Effective kernel mean")+
  theme_minimal()

ggplot(df_by_t_selected,aes(x=t))+
  geom_ribbon(aes(ymin=fund_mean_mean-sqrt(fund_mean_var),ymax=fund_mean_mean+sqrt(fund_mean_var),fill=sim_num),alpha=0.1)+
  geom_line(aes(y=fund_mean_mean,color=sim_num))+
  labs(x="time",y="Fundamental kernel mean\n(mean +/- sd within timesteps)",title="Simulation dynamics, h=0.9: Fundamental kernel mean")+
  theme_minimal()

ggplot(df_by_t_selected,aes(x=t))+
  geom_ribbon(aes(ymin=v_pmeans-sqrt(v_pvars),ymax=v_pmeans+sqrt(v_pvars),fill=sim_num),alpha=0.1)+
  geom_line(aes(y=v_pmeans,color=sim_num))+
  labs(x="time",y="plasticity\n(mean +/- sd within timesteps)",title="Simulation dynamics, h=0.9: Plasticity parameter")+
  theme_minimal()

ggplot(filter(df_by_t_selected,t<100),aes(x=t))+
  geom_line(aes(y=v_abund,color=sim_num))+
  labs(x="time",y="abundance",title="Simulation dynamics, h=0.9: Abundance")+
  theme_minimal()

print(all_params[sim_nums,1:3])
