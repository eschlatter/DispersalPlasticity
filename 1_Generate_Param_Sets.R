source('0_Setup.R')

########## Param set 1: 
nx <- 20 # size of space in the x dimension
ny <- 20 # size of space in the y dimension
nsteps <- 10 # timesteps
# dispersal kernel is a gamma distribution, shape=alpha, scale=theta
v_alphas <- seq(from=0.01,to=5,length.out=5) # values the shape parameter can take
v_thetas <- seq(from=0.01,to=5,length.out=5) # values the scale parameter can take
alpha_start <- 3 # index (in v_alphas) of shape parameter initial value
theta_start <- 3 # index in (v_thetas) of scale parameter initial value
v_p <- seq(from=0,to=0.9,length.out=2) # values the plasticity parameter can take
p_start <- 1 # index (in v_p) of plasticity parameter initial value
mu <- 0.01 # mutation frequency
b <- 10 # reproductive output
K <- round(runif(n=nx*ny,min=3,max=15)) # carrying capacity (also what plasticity responds to)
#pars_save <- list(nx=nx,ny=ny,nsteps=nsteps,v_alphas=v_alphas,alpha_start=alpha_start,v_thetas=v_thetas,theta_start=theta_start,v_p=v_p,p_start=p_start,mu=mu,b=b,b_bad=b_bad,b_neutral=b_neutral,b_good=b_good,K=K)
#save(pars_save,file='params/ParSet1.RData')

patch_locations=NULL
b_neutral=b
b_bad=4
b_good=16
