source('0_Setup.R')
library(grid)
output_folder <- "experiments/Exp3_20260528/output/rep2"
experiment_name <- "Exp3_20260528"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
plotdynamics <- FALSE

theta_increment <- 0.5

raw_files <- grep(".csv",list.files(output_folder),value=TRUE)
v_maps <- numeric()
for(i in raw_files){
  v_maps <- c(v_maps,as.numeric(strsplit(i,"_")[[1]][2]))
}

## Put results of all sims together and plot summaries ----

df_out <- data.frame(mapID=v_maps,repID=NA,mean_abund=NA,mean_p=NA,mean_theta=NA,
                     mean_efftheta=NA,moran_p=NA,moran_theta=NA,moran_efftheta=NA,
                     corr_p=NA,corr_theta=NA,corr_efftheta=NA)
max_t <- vector()
for(i in 1:length(v_maps)){
  mapID <- v_maps[i]
  dat_out <- read.csv(paste0(output_folder,"/map_",mapID,"_raw.csv"))
  abund_df <- filter(dat_out,metric=="abund")
  dat_out <- filter(dat_out,metric!="abund") %>%
    mutate(metric=factor(metric,levels=c("p","theta","efftheta")))
  
  max_t[i] <- max(dat_out$t_i)
  ## results to return for summary
  min_timestep <- min(500,max(dat_out$t_i))
  summ_df <- dat_out %>%
    filter(t_i>=min_timestep) %>%
    group_by(metric) %>%
    summarize(mean=mean(mean),moran=mean(MoranI),corr=mean(corr_q)) %>%
    pivot_wider(names_from=c("metric"),values_from=c("mean","moran","corr"))
  if(nrow(summ_df)!=0){
    df_out_i <- cbind(data.frame(mapID=mapID,repID=1,mean_abund=mean(filter(abund_df,t_i>=min_timestep)$mean)),summ_df)
    df_out[i,] <- df_out_i
  }
}
df_results <- left_join(df_out,experiment_index,by=join_by(mapID))

g_plast <- ggplot(df_results,aes(x=range_ident))+
  geom_hline(aes(yintercept=0),linetype="dashed",alpha=0.5)+
  geom_point(aes(y=mean_p*theta_increment,color=target_dist),position = position_dodge(width=.01))+ # multiply p*theta interval to get plasticity in km
  labs(x=NULL,y="Plasticity value (km)",color="Habitat Quality\nDistribution")+
  scale_color_discrete(breaks=c("identity","C","D","E"),labels = c("Normal", "Bimodal", "Mostly bad","Mostly good"))+
  theme_minimal()+
  scale_x_log10()+
  theme(legend.position = "top",panel.border = element_rect(colour = "gray", fill = NA, linewidth = 0.75))

g_fundkern <- ggplot(df_results,aes(x=range_ident))+
  geom_point(aes(y=mean_theta,color=target_dist),position = position_dodge(width=.01))+
  labs(x=NULL,y="Fundamental Kernel mean (km)",color="Habitat Quality\nDistribution")+
  scale_color_discrete(breaks=c("identity","C","D","E"),labels = c("Normal", "Bimodal", "Mostly bad","Mostly good"))+
  theme_minimal()+
  scale_x_log10()+
  theme(legend.position = "none",panel.border = element_rect(colour = "gray", fill = NA, linewidth = 0.75))

g_effkern <- ggplot(df_results,aes(x=range_ident))+
  geom_point(aes(y=mean_efftheta,color=target_dist),position = position_dodge(width=.01))+
  labs(x="Spatial scale of habitat quality autocorrelation (m)",y="Effective Kernel mean (km)",color="Habitat Quality\nDistribution")+
  scale_color_discrete(breaks=c("identity","C","D","E"),labels = c("Normal", "Bimodal", "Mostly bad","Mostly good"))+
  theme_minimal()+
  scale_x_log10()+
  theme(legend.position = "none",panel.border = element_rect(colour = "gray", fill = NA, linewidth = 0.75))

allplot <- grid.arrange(g_plast,g_fundkern,g_effkern)
ggsave(filename=paste0(output_folder,"/plots/allplot.png"),
       plot=allplot,
       width=1600,
       height=2400,
       units="px")

experiment_index$target_dist_name <- case_match(experiment_index$target_dist,"identity"~"Normal","C"~"Bimodal","D"~"Mostly bad","E"~"Mostly good")
## Plot individual sim dynamics ----
#sims_to_plot <- sample(v_maps,5)
sims_to_plot <- v_maps
sims_to_plot <- 17:20
for(i in sims_to_plot){
  mapID <- i
  
  dat_out <- read.csv(paste0(output_folder,"/map_",mapID,"_raw.csv"))
  abund_df <- filter(dat_out,metric=="abund")
  dat_out <- filter(dat_out,metric!="abund") %>%
    mutate(metric=factor(metric,levels=c("p","theta","efftheta")))
  prows <- dat_out$metric=="p"
  dat_out$mean[prows] <- theta_increment*dat_out$mean[prows]
  dat_out$var[prows] <- theta_increment^2*dat_out$var[prows]
  
  ## qmap
  experiment_i <- filter(experiment_index,mapID==i)
  qrast <- rast(paste0(experiment_folder,"/b",experiment_i$basemap_id,"/set1/qmap_b",
                       experiment_i$basemap_id,"_q",experiment_i$q_ID,".tif"))
  temp_rast <- rast(ext(qrast), resolution=res(qrast), crs = crs(qrast))
  values(temp_rast) <- f_TransformDist(matrix(values(qrast$q),nrow=dim(qrast)[1],byrow=TRUE),experiment_i$target_dist)
  load(file=paste0(experiment_folder,"/habfiles2/habparams_",i,".RData"))
  reef_sf <- hab_params$reef_sf
  temp_rast_cropped <- crop(temp_rast,vect(reef_sf),mask=TRUE)
  qrast_plot <- ggplot()+
    ggspatial::layer_spatial(temp_rast_cropped)+
    scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
    #geom_sf(data=hab_params$sfc_patches,size=0.05)+
    annotation_scale()+
    labs(title="Habitat quality")+
    theme(legend.position = "bottom")
  
  ## dynamics plots
  g_mean <- ggplot(dat_out,aes(x=t_i))+
    geom_hline(aes(yintercept=0),linetype='dashed',alpha=0.5)+
    geom_line(aes(y=mean))+
    geom_ribbon(aes(ymin=mean-sqrt(var),ymax=mean+sqrt(var)),alpha=0.2)+
    facet_wrap(vars(metric),scales="free",nrow=3, strip.position = "left",
               labeller = as_labeller(c("p"="Plasticity\n(km)","theta"="Kernel\nMean (km)","efftheta"="Effective Kernel\nMean (km)")))+
    labs(y=NULL)+
    theme(strip.placement = "outside",strip.text = element_text(size = 12))
  g_moran <- ggplot(filter(dat_out,metric!="abund"),aes(x=t_i))+
    geom_line(aes(y=MoranI))+
    facet_wrap(vars(metric),nrow=3)+
    labs(y=NULL)+
    theme(strip.background = element_blank(),strip.text.x = element_blank())
  g_corr <- ggplot(filter(dat_out,metric!="abund"),aes(x=t_i))+
    geom_hline(aes(yintercept=0),linetype='dashed',alpha=0.5)+
    geom_line(aes(y=corr_q))+
    facet_wrap(vars(metric),nrow=3)+
    labs(y=NULL)+
    theme(strip.background = element_blank(),strip.text.x = element_blank())
  
  #g_abund <- ggplot(abund_df,aes(x=t_i,y=mean))+geom_line()
  
  gplots <- grid.arrange(
    arrangeGrob(g_mean, top = "Value (mean +/- sd)"),
    #arrangeGrob(g_moran, top = "Spatial Autocorrelation\n(Moran's I)"),
    arrangeGrob(g_corr, top = "Correlation with Habitat Quality"),
    qrast_plot,
    ncol = 3,top=textGrob(paste0("Map ",mapID,", ",experiment_index$target_dist_name[i]), gp = gpar(fontsize = 14, fontface = "bold"))
  )
  
  ggsave(filename=paste0(output_folder,"/plots/map_",i,".png"),
         plot=gplots,
         width=2100,
         height=1500,
         units="px")
}
