library(dplyr)
library(urca)
library(zoo)

load('ts_sets.RData')
v_t <- ts_sets[[6]]
v_t <- sim_loop1$v_pmeans

nt=length(v_t)
dat=data.frame(t=1:nt,p=v_t)
dat <- mutate(dat, lag1 = lag(p), diff=p-lag1)

a <- ur.za(dat$diff) # Zivot-Andrews test on diffs
# (Zivot-Andrews is running a bunch of ADFs, I think)
# reject the null ==> stationary
# using the difference, these reject the null by a LOT
summary(a)

par(mfrow=c(2,1))
plot(v_t,type='l')
abline(v=a@bpoint,col='red')
plot(dat$diff,type='l')
abline(v=a@bpoint,col='red')
par(mfrow=c(1,1))

f_FindEquil3 <- function(v_t,showplot=FALSE){
  nt=length(v_t)
  
  # use Zivot-Andrews test for stationarity with a breakpoint
  # (our time series isn't stationary, but a 1-lagged difference is)
  v_t_diff <- diff(v_t)
  za_mod <- ur.za(v_t_diff) 
  breakpt <- za_mod@bpoint
  breakpt_pad <- breakpt+min(breakpt/2,(nt-breakpt)/2) # add a little cushion

  if(za_mod@teststat>za_mod@cval[2]) return("not stationary") # check that it's stationary
  
  if(breakpt_pad>(nt/2)) breakpt_pad=NA # return NA if it's more than halfway through the time series
  
  if(showplot==TRUE){
    par(mfrow=c(2,1))
    plot(v_t,type='l')
    abline(v=breakpt,col='blue')
    abline(v=breakpt_pad,col='red')
    plot(v_t_diff,type='l')
    abline(v=breakpt,col='blue')
    abline(v=breakpt_pad,col='red')
    par(mfrow=c(1,1))
  }
  
  return(breakpt_pad)
}


# can we do something simpler with the diffs? Look at the moving average?
ma <- rollmean(dat$diff,20)
plot(ma,type='l')



## Look at when the mean stabilizes
nt=length(v_t)
dat=data.frame(t=1:nt,p=v_t)
dat <- mutate(dat, lag1 = lag(p), diff=p-lag1)

breakpts <- seq(from=0,to=nt,length.out=min(50,nt))


# Take a vector of a time series, and identify the point when it has reached equilibrium
# this is an ad hoc method, but it seems like it kind of works:
#   it's at least good at cutting out enough of the early portion of the time series.
#   less good at detecting smaller changes in equilibrium level after the first move away from initial conditions.
# inputs: v_t, vector of time series data
# returns: equil_start, timepoint where equilibrium has been reached
f_FindEquil2 <- function(v_t,showplot=FALSE){
  nt=length(v_t)
  equil_pt <- NULL
  
  for(i in 1:19){ # try up to 19 possible start points
    # split data into training and testing
    start_pt=i*nt/20
    test <- v_t[start_pt:nt]
    
    # make an ARIMA model with the test data
    m1 <- auto.arima(test,seasonal=FALSE)
    resids <- checkresiduals(m1,plot=FALSE)
    
    # if so, and if the bounds of prediction aren't enormous, then choose this timepoint and stop the loop
    # (i.e., the time series has reached equilibrium by the end of the training segment. We can use the test segment as equilibrium data.)
    if(resids$p.value>0.05){
      equil_pt <- start_pt
      break
    } 
  }
  
  if(showplot==TRUE & !is.null(equil_pt)){
    plot(1:nt,v_t,type='l')
    abline(v=equil_pt,col='red')
    abline(h=mean(test),col='blue')
  }
  
  if(is.null(equil_pt)) print("Warning: no equilibrium found")
  return(equil_pt) # if no equilibrium point is found, then the returned value is nt.
}

v_t <- ts_sets[[6]]
f_FindEquil2(v_t,showplot = TRUE)


# Zivot-Andrews test
library(urca)
a <- ur.za(v_t)
summary(a)


# stationarity tests generally don't find these time series are stationary ever

# Augmented Dickey-Fuller test for stationarity
# null hypothesis: there is a unit root
# rejection of the null ==> stationary
library(tseries)
a <- adf.test(v_t[5000:10000])

# KPSS for stationarity
# null hypothesis is stationarity
a <- ur.kpss(v_t[5000:10000])

# keep this one!
# --------------------------------------
### another idea is something to do with forecasting. Use some early part of the dataset to forecast the rest.
# If it's totally wrong (i.e., it's predicting based on a trend in the early data that goes away), then use more, until it's in the ballpark.
# Then throw away whatever you had to use, and keep the rest.
load('ts_sets.RData')
library(forecast)
this_set=ts_sets[[8]]
plot(this_set,type='l')

nt=length(this_set)

for(i in 1:20){
  end_train=i*nt/20
  train <- this_set[1:end_train]
  test <- this_set[(end_train+1):nt]
  
  m1 <- auto.arima(train,seasonal=FALSE)
  m1_fore <- forecast(m1,nt-end_train,level=c(95))
  
  in_bounds <- (test>m1_fore$lower)&(test<m1_fore$upper)
  pct_in_bounds <- sum(in_bounds)/length(in_bounds)
  if(pct_in_bounds>0.95 & ((as.numeric(tail(m1_fore$upper,1))-as.numeric(tail(m1_fore$lower,1)))<(100*diff(range(this_set))))){
    equil_start <- end_train
    break
  }
}

equil_start

par(mfrow=c(2,1))
end_train=equil_start
train <- this_set[1:end_train]
test <- this_set[(end_train+1):nt]
plot(1:nt,this_set,type='l')
lines(1:length(train),train,col='blue')

m1 <- auto.arima(train,seasonal=FALSE)
m1_fore <- forecast(m1,nt-end_train,level=c(95))
plot(m1_fore)
lines(1:nt,this_set)




### playing around with time series
ts_sets=list()
ts_sets <- c(ts_sets,list(sim_loop1$v_pmeans))
save(ts_sets,file='ts_sets.RData')

library(forecast)
dat <- data.frame(t=1:nsteps,p=sim_loop1$v_pmeans)
dat <- as.ts(dat)

tsdisplay(sim_loop1$v_pmeans,lag.max=400)

model1 <- Arima(sim_loop1$v_pmeans,c(1,0,0))
summary(model1)
checkresiduals(model1)
plot(1:nsteps,sim_loop1$v_pmeans,type='l')
lines(fitted(model1),col='blue')

autoarima <- auto.arima(sim_loop1$v_pmeans)
summary(autoarima)
checkresiduals(autoarima)
plot(1:nsteps,sim_loop1$v_pmeans,type='l')
lines(fitted(autoarima),col='blue')

dat_subsamp <- sim_loop1$v_pmeans[seq(from=1000,to=10000,by=100)]
tsdisplay(dat_subsamp)
autoarima_subsamp <- auto.arima(dat_subsamp)
summary(autoarima_subsamp)
checkresiduals(autoarima_subsamp)
plot(1:length(dat_subsamp),dat_subsamp,type='l')
lines(fitted(autoarima_subsamp),col='blue')

dat_end <- sim_loop1$v_pmeans[7000:10000]
autoarima_end <- auto.arima(dat_end)
summary(autoarima_end)
checkresiduals(autoarima_end)


## what if we just try a moving average
load('ts_sets.RData')
library(TTR)
lagval=5
subsamp=50
this_set=ts_sets[[5]]
nt=length(this_set)
dat <- data.frame(p=this_set,t=1:nt) %>%
  filter(t %in% seq(from=1,to=nt,by=subsamp)) %>% # first, subsample
  mutate(lag20=lag(p,lagval),diff20=p-lag20) # then take a lagged difference

crit <- vector()
for(i in (1+5):(nrow(dat))){
  # model the differenced time series
  marima <- auto.arima(dat$diff20[i:200])
  # check whether there is a nonzero mean in the fitted model
  crit[i] <- "intercept" %in% names(coef(marima)) 
}
# timestep to start is the first one without a nonzero mean (i.e., where mean=0)
crit_val <- dat$t[min(which(crit==FALSE))]

p1 <- ggplot(data.frame(t=1:nt,p=this_set),aes(x=t,y=p))+
  geom_line()+
  theme_minimal()+
  geom_vline(xintercept=crit_val,color='red')+
  labs(title='original time series')

p2 <- ggplot(dat,aes(x=t))+
  geom_line(aes(y=diff20))+
  geom_vline(xintercept=crit_val,color='red')+
  theme_minimal()+
  labs(title=paste0('lagged difference, subsamp=',subsamp,", lag=",lagval))

grid.arrange(p1,p2,ncol=1)







# first, subsample
dat_sub <- data.frame(t=1:length(dat),p=dat)[seq(from=1,to=10000,by=50),]
# then take a moving average
dat_sub$mavg <- SMA(dat_sub$p,n=10)

plot(dat_sub$mavg,type='l')

nt <- length(dat_sub)
isflat <- vector(length=nt)
for(i in 15:nt){
  m1 <- lm(mavg~t,data=dat_sub[i:nt,])
  isflat[i] <- m1$coefficients[2] # p-value of the linear relationship between t(time) and p(plasticity). When it's not distinguishable from 0, the time series is flat enough.
}

plot(dat_sub$mavg,type='l')
trendmodel <- lm(mavg~t,data=dat_sub[20:200,])
summary(trendmodel)

plot(dat_sub$mavg,type='l')
lines(20:200,fitted(trendmodel))

m1 <- gls(p~t,
          data=dat,
          correlation=corAR1(form= ~ t))
summary(m1)

acf(dat$p,lag.max=9000)
c <- pacf(dat$p)

dat <- mutate(dat, lag1 = lag(p), diff=p-lag1)
plot(dat$t,dat$diff,type='l')
abline(h=0,col='red')
acf(dat$diff[2:10000])

library(tseries)
adf <- adf.test(sim_loop1$v_pmeans[100:10000])$p.value

adf <- vector()
test_steps <- seq(from=100,to=nsteps,by=100)


for(i in 1:length(test_steps)){
  adf[i] <- adf.test(sim_loop1$v_pmeans[test_steps[i]:10000])$p.value
}

acf(sim_loop1$v_pmeans[seq(from=100,to=10000,by=100)])

adf.test(dat$diff[5000:10000])
