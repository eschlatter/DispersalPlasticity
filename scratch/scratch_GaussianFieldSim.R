if(version$major=="4" & version$minor=="4.0") .libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old") else .libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
library(gstat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
#library(sp)

# unconditional simulation on a 100 x 100 grid
var.range=500
xy <- expand.grid(x=seq(from=1,to=100,by=1), y=seq(from=1,to=100,by=1))
g.dummy <- gstat(formula = z~1,
                 locations = ~x+y,
                 dummy = TRUE,
                 beta = 0,
                 model = vgm(psill=1, model="Exp", range=var.range),
                 nmax = 5)
yy <- predict(g.dummy, newdata = xy, nsim = 1) |>
  pivot_longer(cols=starts_with("sim"),names_to="sim",values_to="z")

p1 <- ggplot(yy,aes(x=x,y=y))+
  geom_tile(aes(fill=z))+
  facet_wrap(vars(sim))

p2 <- ggplot(yy,aes(x=z,fill=sim,group=sim))+
  geom_histogram(alpha=0.5)+
  facet_wrap(vars(sim))+
  geom_vline(aes(xintercept = 0),linetype='dashed')

grid.arrange(p1,p2,top=paste0("Range=",var.range),nrow=1)

# when the range gets too large relative to the map diameter,
#   we start getting distributions that don't look normal and zero-centered anymore.
# if we center and normalize them, will they still have the same range?

# first, let's fit variograms to the existing ones:

v_ranges <- seq(from=10,to=100,by=10)
df_fits <- data.frame(actual_range=vector(mode="numeric",length=length(v_ranges)*10),
                      fit_range=vector(mode="numeric",length=length(v_ranges)*10))
df_row <- 0
for(i in 1:length(v_ranges)){
  var.range=v_ranges[i]
  print(var.range)
  xy <- expand.grid(x=seq(from=1,to=100,by=1), y=seq(from=1,to=100,by=1))
  g.dummy <- gstat(formula = z~1,
                   locations = ~x+y,
                   dummy = TRUE,
                   beta = 0,
                   model = vgm(psill=1, model="Exp", range=var.range),
                   nmax = 5)
  for(rep_i in 1:10){
    df_row <- df_row+1
    yy <- predict(g.dummy, newdata = xy, nsim = 1)
    coordinates(yy) <- ~x+y
    yy_gstat <- gstat(id="sim1",formula=sim1~1, data=yy)
    yy_variog <- variogram(yy_gstat,width=5)
    yy_variog_fit <- fit.variogram(yy_variog,vgm(c("Gau","Sph","Exp")))
    df_fits[df_row,] <- data.frame(actual_range=var.range,
                                   fit_range=yy_variog_fit$range[2])  
  }
}

# even these aren't that good!
# need to pick the right variogram width and cutoff, probably

# first one: width=10, cutoff=500
df_fits_1 <- df_fits
ggplot(df_fits_1,aes(x=actual_range,y=fit_range))+geom_point()+geom_abline(aes(slope=1,intercept=0))+
  lims(y=c(0,200))+
  labs(title="Width=10,cutoff=500")

# second one: width=1, cutoff=500
df_fits_2 <- df_fits
ggplot(df_fits_2,aes(x=actual_range,y=fit_range))+geom_point()+geom_abline(aes(slope=1,intercept=0))+
  lims(y=c(0,200))+
  labs(title="Width=1,cutoff=500")

# third one: width=1, cutoff=100
df_fits_3 <- df_fits
ggplot(df_fits_3,aes(x=actual_range,y=fit_range))+geom_point()+geom_abline(aes(slope=1,intercept=0))+
  lims(y=c(0,500))+
  labs(title="Width=1,cutoff=100")

# fourth one: width=1, cutoff=sqrt(200000)/3 (default)
df_fits_4 <- df_fits
ggplot(df_fits_4,aes(x=actual_range,y=fit_range))+geom_point()+geom_abline(aes(slope=1,intercept=0))+
  lims(y=c(0,200))+
  labs(title="Width=1,cutoff=149 (default)")

## Okay. Based on this (don't know if it generalizes to other sizes/resolutions of map), I've learned:
# Variogram fits aren't that accurate.
# They're only reliable at all up to about 1/2 the diameter of the map (up to 50, in this case).

v_ranges <- seq(from=10,to=50,by=10)
df_fits <- data.frame(actual_range=vector(mode="numeric",length=length(v_ranges)*10),
                      fit_range=vector(mode="numeric",length=length(v_ranges)*10),
                      normalized_fit_range=vector(mode="numeric",length=length(v_ranges)*10))
df_row <- 0
for(i in 1:length(v_ranges)){
  var.range=v_ranges[i]
  print(var.range)
  xy <- expand.grid(x=seq(from=1,to=100,by=1), y=seq(from=1,to=100,by=1))
  g.dummy <- gstat(formula = z~1,
                   locations = ~x+y,
                   dummy = TRUE,
                   beta = 0,
                   model = vgm(psill=1, model="Exp", range=var.range),
                   nmax = 5)
  for(rep_i in 1:10){
    df_row <- df_row+1
    yy <- predict(g.dummy, newdata = xy, nsim = 1)
    coordinates(yy) <- ~x+y
    yy_gstat <- gstat(id="sim1",formula=sim1~1, data=yy)
    yy_variog <- variogram(yy_gstat,width=5)
    yy_variog_fit <- fit.variogram(yy_variog,vgm(c("Gau","Sph","Exp")))
    
    yy$sim1_norm <- (yy$sim1-mean(yy$sim1))/var(yy$sim1)
    nyy_gstat <- gstat(id="sim1",formula=sim1_norm~1, data=yy)
    nyy_variog <- variogram(nyy_gstat,width=5)
    nyy_variog_fit <- fit.variogram(nyy_variog,vgm(c("Gau","Sph","Exp")))
    df_fits[df_row,] <- data.frame(actual_range=var.range,
                                   fit_range=yy_variog_fit$range[2],
                                   normalized_fit_range=nyy_variog_fit$range[2])  
  }
}

# Let's try to fit the normalized data, too
ggplot(df_fits,aes(x=fit_range,y=normalized_fit_range))+
  geom_point(aes(color=actual_range))+
  #lims(x=c(0,100),y=c(0,100))+
  geom_abline(aes(slope=1,intercept=0))+
  labs(title="Width=5,cutoff=149 (default)")

# Okay! So, normalizing isn't a problem. Doesn't change the fit at all. Not an iota.
# But, the range we put in still doesn't match the range we get out particularly well.

