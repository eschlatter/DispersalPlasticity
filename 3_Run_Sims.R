source('0_Setup.R')

#### Load a parameter set (or skip this, and just set them directly below)
# load('params/ParSet4.RData')
# list2env(x=params,envir=environment())

#### Set parameters (or adjust them as desired, if they were loaded above)
nsteps <- 100 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- seq(from=0.01,to=1,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.01,to=0.9,length.out=5) # values the scale parameter can take
alpha_start <- 1 # index (in v_alphas) of shape parameter initial value
theta_start <- 1 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=0,to=0.9,length.out=5)# values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
disturb_prob=0 # probability of disturbance (per 10 timesteps)
keep=list("abund","p","kern") # types of summary stats to collect at each timestep

params <- list(nsteps=nsteps,v_alphas=v_alphas,v_thetas=v_thetas,v_p=v_p,
               alpha_start=alpha_start,theta_start=theta_start,p_start=p_start,
               mu=mu,disturb_prob=disturb_prob,
               keep=keep,seed=NULL)

#### Run the seascape-generating pipeline (see the last section of 2_GenerateSeascapes.R) to obtain make_hab_out

#### Run sim
# Note: I've been working with f_RunMatrixLoopLite (in functions/f_RunMatrixLoop.R), which doesn't store all the population info at each timestep; it only stores summary stats.
# We may want to go back to keeping everything at each timestep (that's f_RunMatrixLoop), at least temporarily. But I haven't updated that version of the function recently.
sim_loop1 <- f_RunSimLite(params,hab_params=make_hab_out,show_plot = TRUE)
f_PlotOutput_Lite_Points(sim_loop1)

patch_locations <- make_hab_out$patch_locations
last_ts <- sim_loop1$last_timestep
group_index <- sim_loop1$group_index

npatch <- nrow(patch_locations)
ngroups <- nrow(group_index)

# melt into a datatable with columns patch, timestep, alpha, theta, p, popsize
# add columns for param values at each time/patch/alpha/theta/p combo, scaled by the population size that has that combo
sim_melt <- as.table(last_ts)
dimnames(sim_melt) <- list(1:npatch,1:ngroups)
sim_melt <- as.data.table(sim_melt)
setnames(sim_melt, c("patch","group","popsize"))

sim_melt[,`:=`(patch = as.numeric(patch),
               group = as.numeric(group))]
sim_melt[,`:=`(alpha = group_index$alpha[group],
               theta = group_index$theta[group],
               p = group_index$p[group])]
sim_melt[,`:=`(alpha_value = v_alphas[alpha],
               theta_value = v_thetas[theta],
               p_value = v_p[p])]

# mean param values at each patch
# first take the sum across all cells of the param*popsize (numerator of the mean)
by_patch <- sim_melt[ ,.(alpha=sum(alpha_value*popsize),
                     theta=sum(theta_value*popsize),
                     p=sum(p_value*popsize),
                     popsize=sum(popsize)),
                  by=patch]

# then divide by total popsize (denominator of the mean)
by_patch[,`:=`(alpha = alpha/popsize,
           theta = theta/popsize,
           p = p/popsize)]

patch_locations$patch <- patch_locations$id

by_patch_xy <- merge(by_patch,patch_locations,by='patch')
by_patch_xy <- mutate(by_patch_xy,kernmean=alpha*theta)

alph <- ggplot(by_patch_xy,aes(x=x,y=y,fill=kernmean))+
  geom_tile()+
  coord_fixed()+
  labs(title="dispersal kernel mean\n(at end of simulation)")

p_plot <- ggplot(by_patch_xy,aes(x=x,y=y,fill=p))+
  geom_tile()+
  coord_fixed()+
  labs(title="plasticity\n(at end of simulation)")

popsize_plot <- ggplot(by_patch_xy,aes(x=x,y=y,fill=popsize))+
  geom_tile()+
  coord_fixed()+
  labs(title="population size\n(at end of simulation)")

q_plot <- ggplot(by_patch_xy,aes(x=x,y=y,fill=q))+
  geom_tile()+
  coord_fixed()+
  labs(title="habitat quality")

K_plot <- ggplot(by_patch_xy,aes(x=x,y=y,fill=round(K)))+
  geom_tile()+
  coord_fixed()+
  labs(title="carrying capacity")

grid.arrange(alph,popsize_plot,p_plot,q_plot,K_plot,ncol=2)



df_params <- summarize(by_patch_xy,popsize=sum(popsize),.by=c(alpha,theta))
df_params$popsize <- scale(df_params$popsize)+0.5

ggplot()+
  xlim(0,20)+
  lapply(1:nrow(df_all), 
         function(i){geom_function(fun=dgamma,
                                   args=list(shape=df_all$Var1[i],scale=df_all$Var2[i]),
                                   color='black',alpha=0.3)} )+
  lapply(1:nrow(df_params), 
         function(i){geom_function(fun=dgamma,
                                   args=list(shape=df_params$alpha[i],scale=df_params$theta[i]),
                                   lwd=df_params$popsize[i],
                                   color='blue',alpha=0.5)} )+  
  theme_minimal()+
  labs(x='distance',y='density')

# save(sim_loop1,file="output/20250120_sim4.RData")
# 
# ### quick diagnostic plots
# load(file="output/20250120_sim3.RData")
# f_PlotOutput_Lite_Points(sim_loop1)






#### Code for checking speed and memory (haven't updated this recently)
# ## time profile
# Rprof(filename='profile_loop.out')
# sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
# Rprof(NULL)
# summaryRprof('profile_loop.out')
# 
# ## memory profile
# Rprofmem(filename='profile_loop.out',threshold=10000)
# sim_loop <- f_RunMatrixLoop(nx,ny,nsteps,v_alphas,v_thetas,v_p,alpha_start,theta_start,p_start,mu,b,K,disturb_prob=0,patch_locations=NULL, seed=NULL)
# Rprofmem(NULL)
# noquote(readLines("profile_loop.out"))
