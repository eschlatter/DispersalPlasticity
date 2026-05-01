source('0_Setup.R')
output_folder <- "experiments/Exp2_20260423/output/rep2"
experiment_name <- "Exp2_20260423"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
plotdynamics <- FALSE

raw_files <- grep(".csv",list.files(output_folder),value=TRUE)
v_maps <- numeric()
for(i in raw_files){
  v_maps <- c(v_maps,as.numeric(strsplit(i,"_")[[1]][2]))
}

## Put results of all sims together and plot summaries ----

df_out <- data.frame(mapID=v_maps,repID=NA,mean_abund=NA,mean_p=NA,mean_theta=NA,
                     mean_efftheta=NA,moran_p=NA,moran_theta=NA,moran_efftheta=NA,
                     corr_p=NA,corr_theta=NA,corr_efftheta=NA)

for(i in 1:length(v_maps)){
  mapID <- v_maps[i]
  dat_out <- read.csv(paste0(output_folder,"/map_",mapID,"_raw.csv"))
  abund_df <- filter(dat_out,metric=="abund")
  dat_out <- filter(dat_out,metric!="abund") %>%
    mutate(metric=factor(metric,levels=c("p","theta","efftheta")))
  
  ## results to return for summary
  min_timestep <- min(200,max(dat_out$t_i))
  summ_df <- dat_out %>%
    filter(t_i>=min_timestep) %>%
    group_by(metric) %>%
    summarize(mean=mean(mean),moran=mean(MoranI),corr=mean(corr_q)) %>%
    pivot_wider(names_from=c("metric"),values_from=c("mean","moran","corr"))
  if(nrow(summ_df)!=0){
    df_out_i <- cbind(data.frame(mapID=mapID,repID=2,mean_abund=mean(filter(abund_df,t_i>=min_timestep)$mean)),summ_df)
    df_out[i,] <- df_out_i
  }
}
df_results <- left_join(df_out,experiment_index,by=join_by(mapID))


ggplot(df_results,aes(x=range_ident))+
  geom_hline(aes(yintercept=0),linetype="dashed",alpha=0.5)+
  geom_point(aes(y=mean_p*0.5,color=target_dist))+ # multiply p*theta interval to get plasticity in km
  labs(x="Spatial scale of habitat quality autocorrelation (m)",y="Plasticity value (km)",color="Habitat Quality\nDistribution")+
  scale_color_discrete(breaks=c("identity","C","D","E"),labels = c("Normal", "Bimodal", "Mostly bad","Mostly good"))

ggplot(df_results,aes(x=range_ident))+
  geom_point(aes(y=mean_theta,color=target_dist))+
  labs(x="Spatial scale of habitat quality autocorrelation (m)",y="Fundamental Kernel mean (km)",color="Habitat Quality\nDistribution")+
  scale_color_discrete(breaks=c("identity","C","D","E"),labels = c("Normal", "Bimodal", "Mostly bad","Mostly good"))

ggplot(df_results,aes(x=range_ident))+
  geom_point(aes(y=mean_efftheta,color=target_dist))+
  labs(x="Spatial scale of habitat quality autocorrelation (m)",y="Effective Kernel mean (km)",color="Habitat Quality\nDistribution")+
  scale_color_discrete(breaks=c("identity","C","D","E"),labels = c("Normal", "Bimodal", "Mostly bad","Mostly good"))


## Plot individual sim dynamics ----
sims_to_plot <- sample(v_maps,5)
for(i in sims_to_plot){
  mapID <- i
  dat_out <- read.csv(paste0(output_folder,"/map_",mapID,"_raw.csv"))
  abund_df <- filter(dat_out,metric=="abund")
  dat_out <- filter(dat_out,metric!="abund") %>%
    mutate(metric=factor(metric,levels=c("p","theta","efftheta")))
  
  ## qmap
  experiment_i <- filter(experiment_index,mapID==i)
  qrast <- rast(paste0(experiment_folder,"/b",experiment_i$basemap_id,"/set1/qmap_b",
                       experiment_i$basemap_id,"_q",experiment_i$q_ID,".tif"))
  temp_rast <- rast(ext(qrast), resolution=res(qrast), crs = crs(qrast))
  values(temp_rast) <- f_TransformDist(matrix(values(qrast$q),nrow=dim(qrast)[1],byrow=TRUE),experiment_i$target_dist)
  qrast_plot <- ggplot()+ggspatial::layer_spatial(temp_rast)+
    scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
    annotation_scale()+labs(title="Habitat quality")+theme(legend.position = "bottom")
  
  ## dynamics plots
  g_mean <- ggplot(dat_out,aes(x=t_i))+
    geom_hline(aes(yintercept=0),linetype='dashed',alpha=0.5)+
    geom_line(aes(y=mean))+
    geom_ribbon(aes(ymin=mean-sqrt(var),ymax=mean+sqrt(var)),alpha=0.2)+
    facet_wrap(vars(metric),scales="free",nrow=3, strip.position = "left",
               labeller = as_labeller(c("p"="\nPlasticity","theta"="Kernel\nMean","efftheta"="Effective\nKernel Mean")))+
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
  
  gplots <- grid.arrange(
    arrangeGrob(g_mean, top = "Value\n(mean +/- sd)"),
    arrangeGrob(g_moran, top = "Spatial Autocorrelation\n(Moran's I)"),
    arrangeGrob(g_corr, top = "Correlation with\nHabitat Quality"),
    qrast_plot,
    ncol = 4,top=paste0("Map ",mapID)
  )
}
