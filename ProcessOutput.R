if(version$major=="4" & version$minor=="4.0") .libPaths("/projects/standard/mrunj/shared/Rlibs_schla103_old") else .libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
setwd("/projects/standard/mrunj/shared/Dispersal_plasticity")
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(terra)
experiment_folder <- "New"
output_folder <- "New/output"
sim_index <- read.csv(paste0(experiment_folder,"/_index_sims.csv"))

#### Get and update progress of all running sims ####
if(file.exists(paste0(experiment_folder,"/_index_progress.csv"))){
  progress_index <- read.csv(paste0(experiment_folder,"/_index_progress.csv"))
} else progress_index <- data.frame(sim_id=numeric(),progress=double())
sim_index <- full_join(progress_index,sim_index,by=join_by(sim_id))
sims_to_check <- which(sim_index$progress<1 | is.na(sim_index$progress)) # the ones that aren't done yet
for(sim_row in sims_to_check){
  sim_i <- sim_index$sim_id[sim_row]
  if(file.exists(paste0(experiment_folder,"/output/",sim_i,"_summary.csv"))){
    sim_summ <- read.csv(paste0(experiment_folder,"/output/",sim_i,"_summary.csv"))
    sim_index$progress[sim_row] <- round(max(sim_summ$t_i)/sim_index$nsteps[sim_row],2)
  }
}
progress_index <- sim_index[,c("sim_id","progress")]
write.csv(progress_index,file=paste0(experiment_folder,"/_index_progress.csv"),row.names = FALSE)

#### Compare multiple sims ####

## First, pick which ones

# Pick by slurm run ID
#slurmID <- "15329060|15347441|15347872|15367600|15368129|15368266" # set of 6 full-pop sims from Model Checks Aug6: looking at qautocorr and normalizing larval output
slurmID <- "17253917"

sims_to_pick <- grep(slurmID,sim_index$slurmJob)
sim_index_sub <- sim_index[sims_to_pick,]

##### Runs from 8/27 ####
library(data.table)
sims_to_pick <- which(sim_index$simDate>="2026_08_27")
sim_index_sub <- filter(sim_index,simDate>="2026_08_27")
for(run_i in 1:nrow(sim_index_sub)){
  simID_i <- sim_index_sub$sim_id[run_i]
  pop_df <- fread(paste0(output_folder,"/",simID_i,"_all.csv"))
  ts_keep <- c(1:25,25*(2:(ceiling(max(pop_df$t_i)/25))))
  pop_df_thin <- filter(pop_df,t_i %in% ts_keep)
  fwrite(pop_df_thin,file=paste0(output_folder,"/",simID_i,"_all_thin.csv"))
}

# # Or pick by filtering sim_index
# sim_index_sub <- filter(sim_index,simDate>="2026_08_10",npatch==200,plasticity_on==FALSE,nsteps>100,!grepl("NA",slurmJob))
# sims_to_pick <- which(sim_index$sim_id %in% sim_index_sub$sim_id)

## Then process the data
keep_dynamics <- TRUE # might want to turn this off if data files get too big
# objects to hold everything
df_out <- data.frame(sim_id=rep(NA,length(sims_to_pick)),median_larv_abund=NA,
                     mean_p=NA,mean_theta=NA,mean_efftheta=NA,
                     var_p=NA,var_theta=NA,var_efftheta=NA,
                     median_p=NA,median_theta=NA,median_efftheta=NA,
                     q05_p=NA,q05_theta=NA,q05_efftheta=NA,
                     q25_p=NA,q25_theta=NA,q25_efftheta=NA,
                     q75_p=NA,q75_theta=NA,q75_efftheta=NA,
                     q95_p=NA,q95_theta=NA,q95_efftheta=NA,
                     corr_p=NA,corr_theta=NA,corr_efftheta=NA)
l_dat_out <- vector(mode="list",length=length(sims_to_pick))
l_abund_out <- vector(mode="list",length=length(sims_to_pick))

for(i in seq_along(sims_to_pick)){
  sim_i <- sims_to_pick[i]
  simID_i <- sim_index$sim_id[sim_i]
  if(file.exists(paste0(output_folder,"/",simID_i,"_summary.csv"))){
    dat_out <- read.csv(paste0(output_folder,"/",simID_i,"_summary.csv"))
    abund_df <- filter(dat_out,metric %in% c("abund","larval_abund")) %>% mutate(sim_id=simID_i) %>% dplyr::select(sim_id,t_i,metric,median)
    dat_out <- filter(dat_out,metric %in% c("p","theta","efftheta")) %>%
      mutate(metric=factor(metric,levels=c("p","theta","efftheta")),sim_id=simID_i)
    
    if(keep_dynamics==TRUE){
      l_dat_out[[i]] <- dat_out
      l_abund_out[[i]] <- abund_df
    }
    
    ## results to return for summary
    min_timestep <- min(1500,max(dat_out$t_i))
    summ_df <- dat_out %>%
      filter(t_i>=min_timestep) %>%
      group_by(metric) %>%
      summarize(mean=median(mean),var=median(var),median=median(median),q05=median(q05),q25=median(q25),q75=median(q75),q95=median(q95),corr=median(corr_q)) %>%
      pivot_wider(names_from=c("metric"),values_from=c("mean","var","median","q05","q25","q75","q95","corr"))
    
    if(nrow(summ_df)!=0){
      df_out_i <- cbind(data.frame(sim_id=simID_i,median_larv_abund=median(filter(abund_df,t_i>=min_timestep)$median)),summ_df)
      df_out[i,] <- df_out_i
    } 
    
  } # if file.exists
} # i in seq_along(sims_to_pick)

dat_out_all <- do.call(rbind,l_dat_out)
abund_out_all <- do.call(rbind,l_abund_out)
df_out <- drop_na(df_out,sim_id)
df_out <- left_join(df_out,sim_index,by=join_by(sim_id))

##### Dynamics #####
if(keep_dynamics==TRUE){
  g_theta <- ggplot(filter(dat_out_all,metric=="theta"),aes(x=t_i,y=median,group=factor(sim_id)))+
    geom_line(aes(color=factor(sim_id)))+
    geom_ribbon(aes(ymin=q05,ymax=q95,fill=factor(sim_id)),alpha=0.1)+
    labs(y="Theta\n(median,5-95%)",x=NULL)+
    theme_minimal()+
    theme(legend.position = "none")
  
  g_p <- ggplot(filter(dat_out_all,metric=="p"),aes(x=t_i,y=median,group=factor(sim_id)))+
    geom_line(aes(color=factor(sim_id)))+
    geom_ribbon(aes(ymin=q05,ymax=q95,fill=factor(sim_id)),alpha=0.1)+
    geom_hline(yintercept = 0,linetype='dashed')+
    labs(y="Plasticity\n(median,5-95%)",x=NULL)+
    theme_minimal()+
    theme(legend.position = "none")
  
  g_th_q_corr <- ggplot(filter(dat_out_all,metric=="theta"),aes(x=t_i,y=corr_q,group=factor(sim_id)))+
    geom_line(aes(color=factor(sim_id)))+
    geom_hline(yintercept = 0,linetype='dashed')+
    labs(y="Correlation,\nTheta vs q",x=NULL)+
    theme_minimal()+
    theme(legend.position = "none")
  
  g_effth_q_corr <- ggplot(filter(dat_out_all,metric=="efftheta"),aes(x=t_i,y=corr_q,group=factor(sim_id)))+
    geom_line(aes(color=factor(sim_id)))+
    geom_hline(yintercept = 0,linetype='dashed')+
    labs(y="Correlation,\nEffective Theta vs q",x=NULL)+
    theme_minimal()+
    theme(legend.position = "none")
  
  g_abund <- ggplot(filter(abund_out_all,metric=="abund"),aes(x=t_i,y=median,group=factor(sim_id)))+
    geom_line(aes(color=factor(sim_id)))+
    labs(y="Adult\nabundance")+
    theme_minimal()
  
  gplots <- grid.arrange(g_theta,g_p,g_th_q_corr,g_effth_q_corr,g_abund,ncol=1,top=paste(slurmID,": ", gsub("Export full, ","",sim_index_sub$notes[1])))
}

##### p-theta plots in final timestep #####
plotlist <- vector(mode="list",length=nrow(sim_index_sub))
for(simrep in 1:nrow(sim_index_sub)){
  simID_i <- sim_index_sub$sim_id[simrep]
  habID_i <- sim_index_sub$hab_id[simrep]
  load(paste0(output_folder,"/",simID_i,"_popsnapshot.RData"))
  pop_df <- left_join(pop_df,patch_locations,by=join_by(patch==id))
  plotlist[[simrep]] <- ggplot(pop_df,aes(x=theta,y=p,color=factor(ancestor)))+
    geom_hline(aes(yintercept=0),linetype='dashed')+
    geom_point(size=0.5)+
    lims(x=c(0,sim_index_sub$theta_start_max[1]),y=c(-55,10))+
    labs(title=paste0(simID_i,": ",gsub("Export full, ","",sim_index_sub$notes[simrep]),": ",length(unique(pop_df$ancestor))," lineages"))+
    theme(legend.position = "none")
  
  # plotlist[[simrep]] <- ggplot(pop_df,aes(x=x,y=y,color=factor(ancestor)))+
  #   geom_point(size=0.2)+
  #   labs(title=paste0(simID_i,": ",gsub("Export full, ","",sim_index_sub$notes[simrep]),": ",length(unique(pop_df$ancestor))," lineages"))+
  #   theme(legend.position = "none")+
  #   coord_fixed()
}
grid.arrange(grobs=plotlist)

##### Multisim comparison #####
df_out$q_autocorr_scale[is.na(df_out$q_autocorr_scale)] <- 0
df_out$q_autocorr_scale <- factor(df_out$q_autocorr_scale,levels=sort(unique(df_out$q_autocorr_scale)),labels=c("short-scale","long-scale"))
df_out$disturbance <- df_out$Dp>0

# create a dispersal column based on the notes column
# sim_index_sub <- mutate(sim_index_sub,dispersal=as.character(1+(grepl("dont normalize",notes)+2*grepl("new dispersal",notes))))
# sim_index_sub$dispersal <- replace_values(sim_index_sub$dispersal,"1"~"normalized","2"~"old method","3"~"new method")
# df_out <- left_join(df_out,dplyr::select(sim_index_sub,sim_id,dispersal),by=join_by(sim_id))

p1 <- ggplot(df_out,aes(x=factor(round(q_autocorr_scale)/1000),y=median_theta,group=disturbance,color=disturbance))+
  geom_crossbar(aes(ymin=q25_theta,ymax=q75_theta,fill=factor(disturbance)),position = position_jitterdodge(),width=0.15,alpha=0.2)+ # median over time of population median and population 90% quantiles
  #geom_errorbar(aes(ymin=median_theta-var_theta,ymax=median_theta+var_theta))+ # median +/- variance over time of population median
  #geom_point(size=2.2,position = position_dodge(width = 0.4))+
  labs(x=NULL,y="dispersal kernel median (km)\n(median over time of population median and 50% quantiles)",color="Disturbance")+
  theme(legend.position="top")
p2 <- ggplot(df_out,aes(x=factor(round(q_autocorr_scale)/1000),y=median_p,group=factor(disturbance),color=factor(disturbance)))+
  geom_crossbar(aes(ymin=q25_p,ymax=q75_p),position = position_dodge(width = 0.4))+ # median over time of population median and population 90% quantiles
  #geom_errorbar(aes(ymin=median_theta-var_theta,ymax=median_theta+var_theta))+ # median +/- variance over time of population median
  geom_point(size=2.2,position = position_dodge(width = 0.4))+
  labs(x=NULL,
       y="plasticity: multiplicative change in kernel in best habitat\n(median over time of population median and 90% quantiles)")+
  theme(legend.position = "none")

p3 <- ggplot(df_out,aes(x=factor(round(q_autocorr_scale)/1000),y=corr_theta,group=factor(disturbance),color=factor(disturbance)))+
  #geom_errorbar(aes(ymin=q05_theta,ymax=q95_theta),position = position_dodge(width = 0.4))+ # median over time of population median and population 90% quantiles
  #geom_errorbar(aes(ymin=median_theta-var_theta,ymax=median_theta+var_theta))+ # median +/- variance over time of population median
  geom_point(size=2.2,position = position_dodge(width = 0.4))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  labs(x="spatial autocorrelation in habitat quality",y="correlation of habitat quality and kernel mean")+
  theme(legend.position = "bottom")

grid.arrange(p1,p2,p3,ncol=1)

ggplot(filter(df_out,normalize_offspring==FALSE),aes(x=factor(q_autocorr_scale),y=mean_theta,color=normalize_offspring))+
  geom_point()+
  labs(x="spatial autocorrelation in habitat quality",y="dispersal kernel mean")+
  theme(legend.position = "none")


##### Map #####
simID_i=11802483
experiment_i <- filter(sim_index,sim_id==simID_i)
qrast <- rast(paste0(experiment_folder,"/Maps/b",experiment_i$basemap_id,"/qmap_b",
                     experiment_i$basemap_id,"_q",experiment_i$qmap_id,".tif"))
load(paste0(experiment_folder,"/Maps/habfiles/hab_",experiment_i$hab_id,".RData"))
map_plot <- ggplot()+
  ggspatial::layer_spatial(qrast$q)+ #aggregate(qrast$q,2): doesn't work great for kimbe, but okay for others
  scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
  #geom_sf(data=hab_params$sfc_patches,size=0.5,color='red')+ # include anemones
  ggspatial::annotation_scale()+
  labs(title="Habitat quality")+
  theme(legend.position = "bottom",axis.text=element_blank(),axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5))
map_plot

##### Single-timestep detail #####
simID_i <- sim_index$sim_id[sims_to_pick[5]]
load(paste0(experiment_folder,"/output/",simID_i,"_popsnapshot.RData")) # load pop snapshot
load(paste0(experiment_folder,"/Maps/habfiles/hab_",sim_index$hab_id[sims_to_pick[1]],".RData")) # load hab params
patch_locations <- hab_params$patch_locations
pop_df <- left_join(pop_df,patch_locations,by=join_by(patch==id))
plot(pop_df$q,pop_df$theta,main=paste(cor(pop_df$q,pop_df$theta)))

load(paste0(experiment_folder,"/output/",simID_i,"_comp_results.RData")) # load comp results
comp_results$q_origin <- patch_locations$q[comp_results$parent]
comp_results$q_dest <- patch_locations$q[comp_results$patch]
par(mfrow=c(2,1))
plot(comp_results$q_origin,comp_results$theta,main=paste("After: corr = ",round(cor(comp_results$q_origin,comp_results$theta),3)))
plot(comp_results$q_dest,comp_results$theta,main=paste("Before: corr = ",round(cor(comp_results$q_dest,comp_results$theta),3)))

g1 <- ggplot()+
  ggspatial::layer_spatial(qrast$q)+ #aggregate(qrast$q,2): doesn't work great for kimbe, but okay for others
  scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
  geom_point(data=pop_df,aes(x=x,y=y,color=theta))+ # include anemones
  scale_color_continuous(palette = 'Blues',name="theta",na.value = "#d2f2f7")+
  ggspatial::annotation_scale()+
  labs(title=paste("Theta-q corr = ", round(cor(pop_df$q,pop_df$theta),2)))+
  theme(legend.position = "bottom",axis.text=element_blank(),axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5))
g2 <- ggplot(pop_df,aes(x=q,y=theta))+geom_point()
grid.arrange(g1,g2,nrow=1)

#### One sim ####
.libPaths("/projects/standard/mrunj/shared/Rlib_schla103")
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(terra)
experiment_folder <- "New"
output_folder <- "New/output"
sim_index <- read.csv(paste0(experiment_folder,"/_index_sims.csv"))

# choose the sim
simID <- sim_index$sim_id[1]
simID <- 32413396
output_file <- paste0(experiment_folder,"/output/",simID)
experiment_i <- filter(sim_index,sim_id==simID)
if(experiment_i$basemap_id %in% c(3,4)) plot_q=FALSE else plot_q=TRUE

dat_out <- read.csv(paste0(output_file,"_summary.csv")) %>%
  filter(metric %in% c("p","theta","efftheta"))

abund_out <- read.csv(paste0(output_file,"_summary.csv")) %>%
  filter(metric %in% c("abund","larval_abund"))

load(paste0(experiment_folder,"/Maps/b",experiment_i$basemap_id,"/pop_b",experiment_i$basemap_id,"_p",experiment_i$popmap_id,".RData"))


##### Dynamics #####
g_theta <- ggplot(filter(dat_out,metric=="theta"),aes(x=t_i,y=median))+
  geom_line()+
  geom_ribbon(aes(ymin=q05,ymax=q95),alpha=0.2)+
  labs(y="Theta\n(median,5-95%)",x=NULL)+
  geom_hline(yintercept=as.numeric(c(min_dist,max_dist/2)),linetype='dashed') # bounds of theta in the lookup table

g_p <- ggplot(filter(dat_out,metric=="p"),aes(x=t_i,y=median))+
  geom_line()+
  geom_ribbon(aes(ymin=q05,ymax=q95),alpha=0.2)+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y="Plasticity\n(median,5-95%)",x=NULL)

g_th_q_corr <- ggplot(filter(dat_out,metric=="theta"),aes(x=t_i,y=corr_q))+
  geom_line()+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y="Correlation,\nTheta vs q",x=NULL)

g_abund <- ggplot(filter(abund_out,metric=="abund"),aes(x=t_i,y=median))+
  geom_line()+
  labs(y="Adult\nabundance")

if(plot_q==TRUE){
  qrast <- rast(paste0(experiment_folder,"/Maps/b",experiment_i$basemap_id,"/qmap_b",
                       experiment_i$basemap_id,"_q",experiment_i$qmap_id,".tif"))
  map_plot <- ggplot()+
    ggspatial::layer_spatial(qrast$q)+ #aggregate(qrast$q,2): doesn't work great for kimbe, but okay for others
    scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
    #geom_sf(data=hab_params$sfc_patches,size=0.05)+ # include anemones
    ggspatial::annotation_scale()+
    labs(title="Habitat quality")+
    theme(legend.position = "bottom",axis.text=element_blank(),axis.ticks = element_blank(),
          plot.title = element_text(hjust=0.5))
} else{
  map_plot <- ggplot(reef_sf)+geom_sf()+ggspatial::annotation_scale()+labs(title="Habitat shape")+
    theme(axis.text=element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust=0.5))
}

gplots <- grid.arrange(
  # arrangeGrob(g_theta),
  # arrangeGrob(g_p),
  arrangeGrob(g_theta,g_p,g_th_q_corr,layout_matrix=matrix(c(1,2,3),ncol=1)),
  arrangeGrob(map_plot,g_abund,layout_matrix = matrix(c(1,1,2),ncol=1)),
  ncol = 2,top=grid::textGrob(simID, gp = gpar(fontsize = 14, fontface = "bold"))
)

group_by(dat_out,metric) |> summarize(median=median(median))

all_out <- read.csv(paste0(output_file,"_raw.csv"))
all_out_last <- filter(all_out,t_i==max(all_out$t_i))

##### Animation
library(gganimate)
simID_i <- 99214704
sim_index_i <- filter(sim_index,sim_id==simID_i)
pop_df <- read.csv(paste0(output_folder,"/",simID_i,"_all.csv"))
g1 <- ggplot(pop_df,aes(x=theta,y=p,color=as.factor(ancestor),group=t_i))+
  geom_point()+
  transition_time(t_i)+
  theme(legend.position = "none")+
  labs(title = 'Timestep: {frame_time}')
animate(g1,renderer = gifski_renderer(),duration=20,end_pause=15)
anim_save(paste0(output_folder,"/",simID_i,"_gganimate.gif"))

##### Single-timestep #####
simID_i <- 22655209
sim_index_i <- filter(sim_index,sim_id==simID_i)
qrast <- rast(paste0("New/Maps/b",sim_index_i$basemap_id,"/qmap_b",sim_index_i$basemap_id,"_q",sim_index_i$qmap_id,".tif"))
load(paste0(output_folder,"/",simID_i,"_popsnapshot.RData"))
ggplot(pop_df)+geom_histogram(aes(x=theta))
pop_df <- left_join(pop_df,patch_locations,by=join_by(patch==id))

thetaplot <- ggplot(pop_df,aes(x=x,y=y,color=theta))+geom_point(size=0.5)+coord_fixed()+labs(title="Dispersal kernel mean")
pplot <- ggplot(pop_df,aes(x=x,y=y,color=p))+geom_point(size=0.5)+coord_fixed()+labs(title="Plasticity")
qplot <- ggplot()+ggspatial::layer_spatial(qrast$q)+
  scale_fill_continuous(palette = 'BluGrn',name="q",na.value = "grey")+
  annotation_scale()+labs(title="Habitat quality")
grid.arrange(thetaplot,pplot,qplot,nrow=1,top=paste0("simID=",simID_i))

ggplot(pop_df,aes(x=theta,y=p,color=q))+geom_hline(aes(yintercept=0),linetype='dashed')+geom_point()+lims(x=c(0,160))

##### Larval output by theta #####
output_folder <- "New/output"
simID_i <- 63819350
load(paste0(output_folder,"/",simID_i,"_larval_output_by_theta.RData"))
ggplot(df_thetas,aes(x=theta,y=output_mean))+
  geom_ribbon(aes(ymin=output_mean-output_sd,ymax=output_mean+output_sd),alpha=0.15)+
  geom_point(size=0.5)+
  geom_line()+
  labs(#title=paste0("max at theta=",df_thetas$theta[which.max(df_thetas$output_mean)]," km"),
    y="larval output (mean +/- sd among sites)")+
  ylim(c(0,1.1))

##### Animation
library(gganimate)
simID <- 32413396
pop_df <- read.csv(paste0(output_folder,"/",simID,"_all.csv"))
g1 <- ggplot(pop_df,aes(x=theta,y=p))+
  geom_point()+
  transition_time(t_i)+
  ease_aes('linear')
animate(g1,renderer = gifski_renderer())
anim_save(paste0(output_folder,"/",simID,"_gganimate.gif"), anim)
