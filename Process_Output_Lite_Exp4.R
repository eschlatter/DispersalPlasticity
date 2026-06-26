source('0_Setup.R')
library(grid)
experiment_name <- "Exp4_20260611"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
plotdynamics <- FALSE

load(file=paste0(experiment_folder,"/params_2.RData"))
v_thetas <- params$v_thetas
#v_thetas <- c(0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28) # for smallerkern
#v_thetas <- 10^(-4:1) # for smallkern
v_p <- -5:5
v_p_val <- (v_thetas[2]/v_thetas[1])^-v_p
df_p <- data.frame(p=v_p,pval=v_p_val)

## Put results of all sims together and plot summaries ----
experiment_index$target_dist_name <- case_match(experiment_index$target_dist,"identity"~"Normal","C"~"Bimodal","D"~"Mostly bad","E"~"Mostly good")

df_all <- list()
abund_df <- list()
elt_i <- 1
max_t <- vector()

subexperiment_names <- c("disturb_frequent_b01")
#subexperiment_names <- c("common_disturb","rare_disturb","no_disturb")

for(subexp_i in subexperiment_names){
  output_folder <- paste0(experiment_folder,"/output/",subexp_i)
  
  for(repID_i in 1){
    raw_files <- grep(".csv",list.files(paste0(output_folder,"/rep",repID_i)),value=TRUE)
    v_maps <- numeric()
    for(file_i in raw_files){
      v_maps <- c(v_maps,as.numeric(gsub("\\D", "", file_i)))
      # v_maps <- c(v_maps,substr(strsplit(file_i,"_")[[1]][1],start=4,stop=10))
      # v_maps <- c(v_maps,as.numeric(strsplit(file_i,"_")[[1]][2]))
    }
    
    for(i in 1:length(v_maps)){
      mapID_i <- v_maps[i]
      if(file.exists(paste0(output_folder,"/rep",repID_i,"/map",mapID_i,"_raw.csv"))){
        dat_out <- read.csv(paste0(output_folder,"/rep",repID_i,"/map",mapID_i,"_raw.csv"))
        dat_out$mapID <- mapID_i
        dat_out$repID <- repID_i
        dat_out$subex <- subexp_i
        abund_df[[elt_i]] <- filter(dat_out,metric=="larval_abund")
        dat_out <- filter(dat_out,metric %in% c("p","theta","efftheta")) %>%
          mutate(metric=factor(metric,levels=c("p","theta","efftheta")))
        prows <- dat_out$metric=="p"
        dat_out[prows,c("mean","var","q05","q25","median","q75","q95")] <- (v_thetas[2]/v_thetas[1])^-dat_out[prows,c("mean","var","q05","q25","median","q75","q95")]
        
        df_all[[elt_i]] <- dat_out
        elt_i <- elt_i+1
        
        max_t[i] <- max(dat_out$t_i)
      }
    }
  }
}

#}
df_all <- do.call(rbind,df_all) |> left_join(experiment_index,by=join_by(mapID))
abund_df <- do.call(rbind,abund_df) |> left_join(experiment_index,by=join_by(mapID))

plot_metrics <- filter(df_all,t_i>50) |>
  group_by(mapID,repID,subex,metric) |> # one value per sim (and metric)
  summarize(mean=mean(mean),q25=mean(q25),q75=mean(q75),range=2893,subex=first(subex)) |> # mean values across timesteps
  ggplot(aes(x=range,y=mean,color=subex))+
  geom_point(position = position_dodge(width=0.1))+
  geom_errorbar(aes(ymin=q25,ymax=q75),position = position_dodge(width=0.1),width=0.2)+
  geom_hline(data=filter(dat_out,metric=="p"),aes(yintercept=1),linetype="dashed")+
  scale_x_log10()+
  facet_wrap(vars(metric),nrow=3,scales="free", strip.position = "left",
             labeller = as_labeller(c("p"="Plasticity\n(kern response\nto best habitat)","theta"="Fundamental\nKernel Mean\n(km)","efftheta"="Effective\nKernel Mean\n(km)")))+
  ggh4x::scale_y_facet(PANEL %in% c(1),trans="log10",breaks = breaks_log(base=2),labels = label_log(base=2))+
  scale_color_discrete(labels=c("common_disturb"="Small and Frequent","rare_disturb"="Large and Infrequent","no_disturb"="None"))+
  theme(legend.position="bottom",legend.direction="vertical",strip.placement = "outside",strip.text = element_text(size = 11))+
  labs(x="scale of habitat quality autocorrelation (m)",
       y="mean (within a simulation over time) \nand 50% IQR (among patches within a sim and timestep)\n(IQR is averaged across timesteps)",
       title="Plasticity and kernel means",
       color="Disturbance")

plot_corrs <- filter(df_all,t_i>50) |>
  group_by(mapID,repID,subex,metric) |> # one value per sim (and metric)
  summarize(corr_q=mean(corr_q),sd_corr_q=sd(corr_q),range=2893,subex=first(subex)) |>
  ggplot(aes(x=range,y=corr_q,color=subex))+
  geom_point(position = position_dodge(width=0.1))+
  geom_errorbar(aes(ymin=corr_q-sd_corr_q,ymax=corr_q+sd_corr_q),position = position_dodge(width=0.1),width=0.2)+
  geom_hline(aes(yintercept=0),linetype='dashed',alpha=0.5)+
  scale_x_log10()+
  facet_wrap(vars(metric),nrow=3, strip.position = "left",
             labeller = as_labeller(c("p"="Plasticity\n(kern response\nto best habitat)","theta"="Fundamental\nKernel Mean\n(km)","efftheta"="Effective\nKernel Mean\n(km)")))+
  scale_color_discrete(labels=c("common_disturb"="Small and Frequent","rare_disturb"="Large and Infrequent","no_disturb"="None"))+
  theme(legend.position="bottom",legend.direction="vertical",strip.placement = "outside",strip.text = element_text(size = 11))+
  labs(x="scale of habitat quality autocorrelation (m)",
       y="mean correlation of variable with habitat quality (within a simulation over time)",
       title="Correlation with habitat quality",color="Disturbance")

allplot <- grid.arrange(plot_metrics,plot_corrs,ncol=2)
if(!dir.exists(paste0(experiment_folder,"/output/plots"))) dir.create(paste0(output_folder,"/plots"))
ggsave(filename=paste0(experiment_folder,"/output/plots/allplot_",format(Sys.time(),"%Y%b%d_%H%M"),".png"),
       plot=allplot,
       width=2400,
       height=2400,
       units="px")

## Plot individual sim dynamics ----
#subexp_i <- "common_disturb"
#subexp_i <- "rare_disturb"
#subexp_label <- "Small and Frequent Disturbance"
#subexp_label <- "Large and Infrequent Disturbance"

subexp_i <- "no_disturb"
subexp_label <- "No Disturbance"

subexp_i <- "disturb_frequent"
subexp_label <- "Frequent Disturbance"

subexp_i <- "disturb_frequent_b01"
subexp_label <- "Frequent Disturbance, Low Fecundity"

subexp_i <- "disturb_frequent_b1000"
subexp_label <- "Frequent Disturbance, High Fecundity"


# get the map numbers
output_folder <- paste0(experiment_folder,"/output/",subexp_i)
v_maps <- numeric()
for(repID_i in 1){
  raw_files <- grep(".csv",list.files(paste0(output_folder,"/rep",repID_i)),value=TRUE)
  for(file_i in raw_files){
    v_maps <- c(v_maps,as.numeric(gsub("\\D", "", file_i)))
  }
}
v_maps <- unique(v_maps)

sims_to_plot <- seq_along(v_maps)
for(i in sims_to_plot){
  map_i <- v_maps[i]
  
  dat_out <- filter(df_all,mapID==map_i & subex==subexp_i)
  abund_df_i <- filter(abund_df,mapID==map_i & subex==subexp_i)
  dat_out <- filter(dat_out,metric %in% c("p","theta","efftheta"))
  
  ## qmap
  experiment_i <- filter(experiment_index,mapID==map_i)
  qrast <- rast(paste0(experiment_folder,"/b",experiment_i$basemap_id,"/set1/qmap_b",
                       experiment_i$basemap_id,"_q",experiment_i$qmap_id,".tif"))
  # qrast <- rast(paste0(experiment_folder,"/b",experiment_i$basemap_id,"/set1/qmap_b",
  # experiment_i$basemap_id,"_q",experiment_i$q_ID,".tif"))
  temp_rast <- rast(ext(qrast), resolution=res(qrast), crs = crs(qrast))
  values(temp_rast) <- f_TransformDist(matrix(values(qrast$q),nrow=dim(qrast)[1],byrow=TRUE),experiment_i$target_dist)
  load(file=paste0(experiment_folder,"/habfiles/habparams_",map_i,".RData"))
  reef_sf <- hab_params$reef_sf
  temp_rast_cropped <- crop(temp_rast,vect(reef_sf),mask=TRUE)
  qrast_plot <- ggplot()+
    ggspatial::layer_spatial(temp_rast_cropped)+
    scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
    #geom_sf(data=hab_params$sfc_patches,size=0.05)+
    annotation_scale()+
    labs(title="Habitat quality")+
    theme(legend.position = "bottom",axis.text=element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust=0.5))
  
  ## dynamics plots
  g_median <- ggplot(dat_out,aes(x=t_i,group=repID,color=repID))+
    geom_hline(data=filter(dat_out,metric=="p"),aes(yintercept=1),linetype="dashed")+
    geom_line(aes(y=mean),alpha=0.75)+
    geom_ribbon(aes(ymin=q05,ymax=q95,color=NULL,fill=repID),alpha=0.2)+
    facet_wrap(vars(metric),scales="free",nrow=3, strip.position = "left",
               labeller = as_labeller(c("p"="Plasticity\n(kern response\nto best habitat)","theta"="Fundamental\nKernel Mean\n(km)","efftheta"="Effective\nKernel Mean\n(km)")))+
    labs(y=NULL)+
    # lims(x=c(0,max(max_t)))+
    ggh4x::scale_y_facet(PANEL %in% c(1),trans="log10",breaks = breaks_log(base=4),labels = label_log(base=4))+
    theme(strip.placement = "outside",strip.text = element_text(size = 11),legend.position="none")
  g_corr <- ggplot(filter(dat_out,metric!="abund"),aes(x=t_i,group=repID,color=repID))+
    geom_hline(aes(yintercept=0),linetype='dashed',alpha=0.5)+
    geom_line(aes(y=corr_q),alpha=0.75)+
    facet_wrap(vars(metric),nrow=3)+
    labs(y=NULL)+
    # lims(x=c(0,max(max_t)))+
    theme(strip.background = element_blank(),strip.text.x = element_blank(),legend.position="none")
  
  g_abund <- ggplot(abund_df_i,aes(x=t_i,y=median,group=repID,color=repID))+
    geom_line(alpha=0.75)+
    labs(title="Larval Abundance",y=NULL)+
    scale_y_continuous(labels = label_number(accuracy = 1))+
    theme(legend.position="none")
  
  gplots <- grid.arrange(
    arrangeGrob(g_median, top = "Mean and\n50% interquartile range"),
    arrangeGrob(g_corr, top = "Correlation with\nHabitat Quality"),
    arrangeGrob(qrast_plot,g_abund,layout_matrix = matrix(c(1,1,2),ncol=1)),
    ncol = 3,top=textGrob(paste0("Map ",map_i,", ",subexp_label), gp = gpar(fontsize = 14, fontface = "bold"))
  )
  
  ggsave(filename=paste0(experiment_folder,"/output/plots/",subexp_i,"_map",map_i,".png"),
         plot=gplots,
         width=2100,
         height=1500,
         units="px")
}


## Look at pop snapshots ---------------
load("experiments/Exp3_20260528/output/common_disturb/rep1/map_13_popsnapshot.R")
load("experiments/Exp3_20260528/habfiles/habparams_13.RData")
patch_locations <- hab_params$patch_locations

v_alphas <- 1
v_thetas <- c(0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28) # for smallerkern
#v_thetas <- 10^(-4:1) # for smallkern
v_p <- -5:5
v_p_val <- (v_thetas[2]/v_thetas[1])^-v_p
group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))

group_index$theta_val <- v_thetas[group_index$theta]
group_index$p_val <- v_p_val[group_index$p]

patch_locations$p_ind <- previous_pop %*% group_index$p
patch_locations$p_val <- previous_pop %*% group_index$p_val

ggplot(filter(patch_locations,p_val>10),aes(x=x,y=y,color=p_val))+
  geom_point()

head(patch_locations)
ggplot(patch_locations,aes(x=q,y=p_ind))+geom_point(alpha=0.2)+geom_smooth()

cor(patch_locations$q,patch_locations$p_val)
