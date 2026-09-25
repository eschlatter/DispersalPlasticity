library(ggplot2)

q <- seq(from=0,to=1,length.out=100)
plasticity_fn <- function(theta_f,p,q){
  theta_e <- theta_f*(exp(2*p*(q-0.5)))
}

v_p <- seq(from=-1,by=0.2,to=1)

ggplot()+
  xlim(0,1)+
  lapply(1:length(v_p),
         function(i){geom_function(fun=plasticity_fn,
                                   args=list(theta_f=1,p=v_p[i]),
                                   aes(color=factor(v_p[i])),
                                   n=100)})+
  theme_minimal()+
  labs(x="q",y="multiplication factor",color="p")