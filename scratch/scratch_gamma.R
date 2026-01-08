# Plot a bunch of gamma functions with the same mean but different alpha and theta values
library(ggplot2)

mean=1
alpha=c(0.01,1,1.01,1.5,2,5,10,100)
theta=mean/alpha

v_alpha=seq(from=0.01,to=5,length.out=5)
v_theta=seq(from=0.01,to=5,length.out=5)
pars <- expand.grid(v_alpha,v_theta)
colnames(pars)=c('alpha','theta')

pars <- subset(pars,theta==v_theta[4])

ggplot()+
  xlim(0,10*mean)+
  lapply(1:nrow(pars),function(i) geom_function(fun=dgamma,args=list(shape=pars$alpha[i],scale=pars$theta[i]),aes(color=paste0('alpha=',pars$alpha[i],'theta=',pars$theta[i]))))+
  geom_vline(aes(xintercept=0.5))
