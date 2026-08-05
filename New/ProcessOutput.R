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
  labs(y="Theta\n(median,5-95%)",x=NULL)+
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

all_out <- read.csv(paste0(output_file,"_raw.csv"))
all_out_last <- filter(all_out,t_i==max(all_out$t_i))
