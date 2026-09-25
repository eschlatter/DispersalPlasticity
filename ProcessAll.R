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
sim_index <- read.csv(paste0(experiment_folder,"/_index_sims.csv")) |> filter(!is.na(hab_id))

basemapID=1
#load(paste0("New/df_experiments_b",basemapID,"_mutation.RData"))
#load(paste0("New/df_experiments_b",basemapID,"_090426.RData"))
load(paste0("New/df_experiments_b",basemapID,".RData"))
df_experiments$rep_id <- 1:nrow(df_experiments)
df_experiments <- df_experiments[,c("experiment_id","De_seq","Dl_seq","Dm_seq","disturb_id")]
sim_index <- left_join(sim_index,df_experiments,by=join_by(experiment_id))

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
slurmID <- "18113615|18114187" # all 5x5 full
slurmID <- "176034" # 5x5 full, larger mutation rates
#slurmID <- "184714" # mutation experiment
sims_to_pick <- grep(slurmID,sim_index$slurmJob)
sim_index_sub <- sim_index[sims_to_pick,]
# # Or pick by filtering sim_index
#sim_index_sub <- filter(sim_index,simDate>"2026_08_27" & simDate<"2026_08_28")
sim_index_sub <- filter(sim_index_sub,Dl==0)
sims_to_pick <- which(sim_index$sim_id %in% sim_index_sub$sim_id)


## Download dynamics to make animations
save(sim_index_sub,file=paste0("toDownload/20260909_largemut_nodisturb_simindex.RData"))
for(i in 1:nrow(sim_index_sub)){
  sim_i <- sim_index_sub$sim_id[i]
  #file.copy(from=paste0("New/output/",sim_i,"_all.csv"),to=paste0("toDownload/",sim_i,"_all.csv"))
  alldat <- read.csv(paste0("New/output/",sim_i,"_all.csv"))
  print(unique(alldat$t_i))
}

alldat <- read.csv(paste0("New/output/",sim_i,"_all.csv"))
unique(alldat$t_i)


##### Add to the index, for the running sims that somehow got missed (overwriting, maybe?)
repID=7
file_lines <- readr::read_tsv(paste0(output_folder,"/slurm/18114187_",repID,".output"),col_names=FALSE)
first_line_number <- grep("simID: ", file_lines$X1)
simID <- readr::parse_number(strsplit(file_lines$X1[first_line_number],":")[[1]][2])
sim_index_i <- file_lines[(first_line_number+2):(first_line_number+30),1] |>
  separate_wider_delim(cols = X1,delim=" ",names=c("metric","value"),too_many="merge") |>
  mutate(value=stringr::str_trim(value)) |>
  mutate(value=gsub("\\\\","",value)) |>
  t()
colnames(sim_index_i) <- sim_index_i[1,]
sim_index_i <- as.data.frame(sim_index_i) 
sim_index_i <- sim_index_i[2,]


target_types <- purrr::map_chr(sim_index, ~ class(.x)[1])
sim_index_i <- sim_index_i %>%
  readr::type_convert()

sim_index_i <- mutate(sim_index_i,across(c(sim_id,), as.numeric))

##### Runs from 8/27 ####
# library(data.table)
# sims_to_pick <- which(sim_index$simDate>="2026_08_27" & sim_index$simDate<"2026_08_28")
# sim_index_sub <- filter(sim_index,simDate>="2026_08_27",simDate<"2026_08_28")
# for(run_i in 2:nrow(sim_index_sub)){
#   simID_i <- sim_index_sub$sim_id[run_i]
#   pop_df <- fread(paste0(output_folder,"/",simID_i,"_all.csv"))
#   ts_keep <- c(1:25,25*(2:(ceiling(max(pop_df$t_i)/25))))
#   pop_df_thin <- filter(pop_df,t_i %in% ts_keep)
#   fwrite(pop_df_thin,file=paste0(output_folder,"/",simID_i,"_all_thin.csv"))
#   # file.remove(paste0(output_folder,"/",simID_i,"_all.csv"))
#   # file.rename(from=paste0(output_folder,"/",simID_i,"_all_thin.csv"),to=paste0(output_folder,"/",simID_i,"_all.csv"))
# }


## Then process the data
keep_dynamics <- FALSE # might want to turn this off if data files get too big
# objects to hold everything
df_out <- data.frame(sim_id=rep(NA,length(sims_to_pick)),median_larv_abund=NA,
                     mean_p=NA,mean_theta=NA,mean_efftheta=NA,
                     var_p=NA,var_theta=NA,var_efftheta=NA,
                     median_p=NA,median_theta=NA,median_efftheta=NA,
                     q05_p=NA,q05_theta=NA,q05_efftheta=NA,
                     q25_p=NA,q25_theta=NA,q25_efftheta=NA,
                     q75_p=NA,q75_theta=NA,q75_efftheta=NA,
                     q95_p=NA,q95_theta=NA,q95_efftheta=NA,
                     corr_p=NA,corr_theta=NA,corr_efftheta=NA,
                     var_corr_p=NA,var_corr_theta=NA,var_corr_efftheta=NA)
if(keep_dynamics==TRUE){
  l_dat_out <- vector(mode="list",length=length(sims_to_pick))
  l_abund_out <- vector(mode="list",length=length(sims_to_pick))
}

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
    min_timestep <- 2800 # min(1500,max(dat_out$t_i))
    summ_df <- dat_out %>%
      filter(t_i>=min_timestep) %>%
      group_by(metric) %>%
      summarize(mean=median(mean),var=median(var),median=median(median),q05=median(q05),q25=median(q25),q75=median(q75),q95=median(q95),corr=median(corr_q),var_corr=var(corr_q)) %>%
      pivot_wider(names_from=c("metric"),values_from=c("mean","var","median","q05","q25","q75","q95","corr","var_corr"))
    
    if(nrow(summ_df)!=0){
      df_out_i <- cbind(data.frame(sim_id=simID_i,median_larv_abund=median(filter(abund_df,t_i>=min_timestep)$median)),summ_df)
      df_out[i,] <- df_out_i
    } 
    
  } # if file.exists
} # i in seq_along(sims_to_pick)

if(keep_dynamics==TRUE){
  dat_out_all <- do.call(rbind,l_dat_out)
  abund_out_all <- do.call(rbind,l_abund_out)
}
df_out <- drop_na(df_out,sim_id)
df_out <- left_join(df_out,sim_index,by=join_by(sim_id))
save(df_out,sim_index_sub,sims_to_pick,file="New/output/df_out_b1_090426.RData")


abund_out_all <- left_join(abund_out_all,sim_index_sub)
g_abund <- ggplot(filter(abund_out_all,metric=="abund"),aes(x=t_i,y=median,group=factor(sim_id)))+
  geom_line(aes(color=factor(De)))+
  labs(y="Adult\nabundance")+
  theme_minimal()
g_abund

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
sim_index_sub_sub <- filter(sim_index_sub,mutation_type=="frequent")
sim_index_sub_sub <- filter(sim_index_sub,disturb_id==11,q_autocorr_scale>2000)
sim_index_sub_sub <- filter(sim_index_sub,sim_id %in% c(88195989, 20176438, 83533336, 47091439))
plotlist <- vector(mode="list",length=nrow(sim_index_sub_sub))
i=0
for(simrep in 1:nrow(sim_index_sub_sub)){
  i=i+1
  simID_i <- sim_index_sub_sub$sim_id[simrep]
  habID_i <- sim_index_sub_sub$hab_id[simrep]
  load(paste0(output_folder,"/",simID_i,"_popsnapshot.RData"))
  pop_df <- left_join(pop_df,patch_locations,by=join_by(patch==id))
  plotlist[[i]] <- ggplot(pop_df,aes(x=theta,y=p, color=q))+
    geom_hline(aes(yintercept=0),linetype='dashed')+
    geom_point(size=0.5)+
    lims(x=c(0,sim_index_sub_sub$theta_start_max[1]),y=c(-30,25))+
    #labs(title=paste0(simID_i,": ",gsub("Export full, ","",sim_index_sub_sub$notes[simrep]),": ",length(unique(pop_df$ancestor))," lineages"))+
    # labs(title=paste0(simID_i,", t = ",t_i,"\nmu_theta = ",round(sim_index_sub_sub$mu_theta[simrep],3),", mu_p = ",round(sim_index_sub_sub$mu[simrep],3)))+
    labs(title=paste0(simID_i,", t = ",t_i))+
    theme(legend.position = "right")
  
  # plotlist[[simrep]] <- ggplot(pop_df,aes(x=x,y=y,color=factor(ancestor)))+
  #   geom_point(size=0.2)+
  #   labs(title=paste0(simID_i,": ",gsub("Export full, ","",sim_index_sub$notes[simrep]),": ",length(unique(pop_df$ancestor))," lineages"))+
  #   theme(legend.position = "none")+
  #   coord_fixed()
}
grid.arrange(grobs=plotlist)

##### Playing around with clustering #####

dist_mat <- dist(pop_df[,c("theta","p")])
hc <- hclust(dist_mat)


##### Multisim comparison #####
load("New/output/df_out_b1_090926.RData")

sim_index_sub <- filter(sim_index_sub,Dm_seq==1)
sims_to_pick <- which(sim_index$sim_id %in% sim_index_sub$sim_id)
df_out <- filter(df_out,sim_id %in% sim_index_sub$sim_id)
# df_out$q_autocorr_scale[is.na(df_out$q_autocorr_scale)] <- 0
# df_out$q_autocorr_scale <- factor(df_out$q_autocorr_scale,levels=sort(unique(df_out$q_autocorr_scale)),labels=c("short-scale","long-scale"))
# df_out$disturbance <- df_out$Dp>0

# create a dispersal column based on the notes column
# sim_index_sub <- mutate(sim_index_sub,dispersal=as.character(1+(grepl("dont normalize",notes)+2*grepl("new dispersal",notes))))
# sim_index_sub$dispersal <- replace_values(sim_index_sub$dispersal,"1"~"normalized","2"~"old method","3"~"new method")
# df_out <- left_join(df_out,dplyr::select(sim_index_sub,sim_id,dispersal),by=join_by(sim_id))

p1 <- ggplot(df_out,aes(x=q_autocorr_scale/1000,y=median_theta,group=1-Dm,color=1-Dm))+
  geom_errorbar(aes(ymin=q25_theta,ymax=q75_theta),position=position_dodge(width=0.1))+ # median over time of population median and population 90% quantiles
  #geom_errorbar(aes(ymin=median_theta-var_theta,ymax=median_theta+var_theta))+ # median +/- variance over time of population median
  geom_point(position=position_dodge(width=0.1))+
  labs(x=NULL,y="dispersal kernel median (km)\n(median over time of population median and 50% quantiles)",color="Disturbance\nSeverity")+
  theme(legend.position="top")
p2 <- ggplot(df_out,aes(x=q_autocorr_scale/1000,y=median_p,group=1-Dm,color=1-Dm))+
  geom_errorbar(aes(ymin=q25_p,ymax=q75_p),position=position_dodge(width=0.1))+ # median over time of population median and population 90% quantiles
  #geom_errorbar(aes(ymin=median_theta-var_theta,ymax=median_theta+var_theta))+ # median +/- variance over time of population median
  geom_point(position=position_dodge(width=0.1))+
  labs(x=NULL,
       y="plasticity: kernel multiplied by e^p in best habitat\n(median over time of population median and 50% quantiles)",
       color="Disturbance")+
  theme(legend.position="none")
p3 <- ggplot(df_out,aes(x=q_autocorr_scale/1000,y=corr_theta,group=1-Dm,color=1-Dm))+
  geom_errorbar(aes(ymin=corr_theta-var_corr_theta,ymax=corr_theta+var_corr_theta),position = position_dodge(width = 0.1))+ # median over time of population median and population 90% quantiles
  geom_point(position = position_dodge(width = 0.1))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  labs(x="spatial autocorrelation in habitat quality",
       y="correlation of habitat quality and kernel mean\n(median over time of population median and population 90% quantiles)")+
  theme(legend.position = "none")

grid.arrange(p1,p2,p3,ncol=1)

ggplot(filter(df_out,normalize_offspring==FALSE),aes(x=factor(q_autocorr_scale),y=mean_theta,color=normalize_offspring))+
  geom_point()+
  labs(x="spatial autocorrelation in habitat quality",y="dispersal kernel mean")+
  theme(legend.position = "none")


##### GAMs #####
library(mgcv)
load("New/output/df_out_b1_090926.RData")
load("New/output/df_out_b1_090426.RData")
df_out$Ds <- 1-df_out$Dm

## Kernel prediction models
# Disturbance severity
kern_vs_Ds <- gam(median_theta~te(q_autocorr_scale, Ds),data=filter(df_out,Dm_seq==1))
#vis.gam(kern_vs_Ds,plot.type = "contour")
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),Ds=seq(0,1,by=0.025))
predictions <- predict(kern_vs_Ds,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g1 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=Ds,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=Ds,z=predictions),color='black')+
  geom_point(data=kern_vs_Ds$model,aes(x=q_autocorr_scale,y=Ds),color='black')+
  theme_minimal()+
  labs(y="Disturbance Severity",fill="Theta")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(0,max(df_out$median_theta)),palette = "Blues")
# Disturbance duration
kern_vs_Dl <- gam(median_theta~te(q_autocorr_scale, Dl),data=filter(df_out,Dl_seq==1))
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),Dl=seq(0,15,by=0.25))
predictions <- predict(kern_vs_Dl,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g2 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=Dl,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=Dl,z=predictions),color='black')+
  geom_point(data=kern_vs_Dl$model,aes(x=q_autocorr_scale,y=Dl),color='black')+
  theme_minimal()+
  labs(y="Disturbance Duration",fill="Dispersal\nKernel")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(0,max(df_out$median_theta)),palette = "Blues")
# Disturbance extent
kern_vs_De <- gam(median_theta~te(q_autocorr_scale, De),data=filter(df_out,De_seq==1))
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),De=seq(0,1,by=0.025))
predictions <- predict(kern_vs_De,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g3 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=De,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=De,z=predictions),color='black')+
  geom_point(data=kern_vs_De$model,aes(x=q_autocorr_scale,y=De),color='black')+
  theme_minimal()+
  labs(y="Disturbance Extent",fill="Dispersal\nKernel")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(0,max(df_out$median_theta)),palette = "Blues")

grid.arrange(g1,g2,g3)


## Plasticity prediction models
# Disturbance severity
p_vs_Ds <- gam(median_p~te(q_autocorr_scale, Ds),data=filter(df_out,Dm_seq==1))
#vis.gam(p_vs_Ds,plot.type = "contour")
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),Ds=seq(0,1,by=0.025))
predictions <- predict(p_vs_Ds,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g1 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=Ds,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=Ds,z=predictions),color='black')+
  geom_point(data=p_vs_Ds$model,aes(x=q_autocorr_scale,y=Ds),color='black')+
  theme_minimal()+
  labs(y="Disturbance Severity",fill="Plasticity")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(min(df_out$median_p),max(df_out$median_p)),palette = "Blues")
# Disturbance duration
p_vs_Dl <- gam(median_p~te(q_autocorr_scale, Dl),data=filter(df_out,Dl_seq==1))
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),Dl=seq(0,15,by=0.25))
predictions <- predict(p_vs_Dl,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g2 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=Dl,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=Dl,z=predictions),color='black')+
  geom_point(data=p_vs_Dl$model,aes(x=q_autocorr_scale,y=Dl),color='black')+
  theme_minimal()+
  labs(y="Disturbance Duration",fill="Plasticity")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(min(df_out$median_p),max(df_out$median_p)),palette = "Blues")
# Disturbance extent
p_vs_De <- gam(median_p~te(q_autocorr_scale, De),data=filter(df_out,De_seq==1))
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),De=seq(0,1,by=0.025))
predictions <- predict(p_vs_De,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g3 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=De,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=De,z=predictions),color='black')+
  geom_point(data=p_vs_De$model,aes(x=q_autocorr_scale,y=De),color='black')+
  theme_minimal()+
  labs(y="Disturbance Extent",fill="Plasticity")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(min(df_out$median_p),max(df_out$median_p)),palette = "Blues")

grid.arrange(g1,g2,g3)

## Local adaptation prediction models
# Disturbance severity
LA_vs_Ds <- gam(corr_theta~te(q_autocorr_scale, Ds),data=filter(df_out,Dm_seq==1))
#vis.gam(p_vs_Ds,plot.type = "contour")
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),Ds=seq(0,1,by=0.025))
predictions <- predict(LA_vs_Ds,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g1 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=Ds,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=Ds,z=predictions),color='black')+
  geom_point(data=LA_vs_Ds$model,aes(x=q_autocorr_scale,y=Ds),color='black')+
  theme_minimal()+
  labs(y="Disturbance Severity",fill="Local adaptation\n(correlation between hab quality and fundamental kernel)")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(min(df_out$corr_theta),max(df_out$corr_theta)),palette = "Blues")
# Disturbance duration
LA_vs_Dl <- gam(corr_theta~te(q_autocorr_scale, Dl),data=filter(df_out,Dl_seq==1))
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),Dl=seq(0,15,by=0.25))
predictions <- predict(LA_vs_Dl,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g2 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=Dl,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=Dl,z=predictions),color='black')+
  geom_point(data=LA_vs_Dl$model,aes(x=q_autocorr_scale,y=Dl),color='black')+
  theme_minimal()+
  labs(y="Disturbance Duration",fill="Local adaptation\n(correlation between hab quality and fundamental kernel)")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "none")+
  scale_fill_continuous(limits=c(min(df_out$corr_theta),max(df_out$corr_theta)),palette = "Blues")
# Disturbance extent
LA_vs_De <- gam(corr_theta~te(q_autocorr_scale, De),data=filter(df_out,De_seq==1))
testdata <- expand.grid(q_autocorr_scale=seq(0,3500,by=0.5),De=seq(0,1,by=0.025))
predictions <- predict(LA_vs_De,newdata=testdata,type='response',se=FALSE)
df_preds <- data.frame(testdata,predictions)
g3 <- ggplot(df_preds)+
  geom_tile(aes(x=q_autocorr_scale,y=De,fill=predictions))+
  geom_contour(aes(x=q_autocorr_scale,y=De,z=predictions),color='black')+
  geom_point(data=LA_vs_De$model,aes(x=q_autocorr_scale,y=De),color='black')+
  theme_minimal()+
  labs(y="Disturbance Extent",fill="Local\nAdaptation")+
  theme(axis.title.x = element_blank(),axis.text.x=element_blank(),legend.position = "right")+
  scale_fill_continuous(limits=c(min(df_out$corr_theta),max(df_out$corr_theta)),palette = "Blues")

grid.arrange(g1,g2,g3)


##### Some kernel plots ######
th_j=50
int_trunc <- calculus::integral(function(x,theta){dgamma(x,shape=1,scale=theta)},
                                   bounds=list(x=c(0,5)),
                                   params=list(theta=th_j))$value
ggplot()+
  lims(x=c(0,5),y=c(0,0.5))+
  geom_function(fun=function(x,theta){(1/int_trunc)*dgamma(x,shape=1,scale=theta)}, args=list(theta=th_j), linewidth=4)+
  labs(x="Distance (km)")+
  geom_text(aes(x=4,y=0.45,label=paste0("theta = ",th_j)),size=14)+
  theme(axis.text.y=element_blank(),axis.text.x=element_text(size=20),
        axis.title.y=element_blank(),axis.title.x=element_text(size=20),
        panel.grid=element_blank())

##### Sample disturbance maps #####

frac_map <- fracland(k=9,h=1.5,p=0.8,binary=TRUE, plotflag=TRUE)


##### Another type of plot: heatmaps? #####
qmap_index <- read.csv("New/Maps/index_qmaps.csv")
df_out2 <- left_join(df_out,qmap_index,by=join_by(qmap_id)) |>
  mutate(Dm_ord=factor(1-Dm, levels=sort(unique(1-Dm)),ordered=TRUE),
         De_ord=factor(De, levels=sort(unique(De)),ordered=TRUE),
         Dl_ord=factor(Dl, levels=sort(unique(Dl)),ordered=TRUE),
         h_ord=factor(h, levels=sort(unique(h)),ordered=TRUE))

ggplot(filter(df_out2,Dm_seq==1))+
  geom_tile(aes(x=Dm_ord,y=h_ord,fill=log(median_theta)))+
  coord_fixed()+
  labs(x="Disturbance Severity",y="Spatial Autocorrelation Scale")
ggplot(filter(df_out2,De_seq==1))+
  geom_tile(aes(x=De_ord,y=h_ord,fill=log(median_theta)))+
  coord_fixed()+
  labs(x="Disturbance Extent",y="Spatial Autocorrelation Scale")
ggplot(filter(df_out2,Dl_seq==1))+
  geom_tile(aes(x=Dl_ord,y=h_ord,fill=log(median_theta)))+
  coord_fixed()+
  labs(x="Disturbance Duration",y="Spatial Autocorrelation Scale")

ggplot(filter(df_out2,Dm_seq==1))+
  geom_tile(aes(x=Dm_ord,y=h_ord,fill=median_p))+
  coord_fixed()+
  labs(x="Disturbance Severity",y="Plasticity")
ggplot(filter(df_out2,De_seq==1))+
  geom_tile(aes(x=De_ord,y=h_ord,fill=median_p))+
  coord_fixed()+
  labs(x="Disturbance Extent",y="Plasticity")
ggplot(filter(df_out2,Dl_seq==1))+
  geom_tile(aes(x=Dl_ord,y=h_ord,fill=median_p))+
  coord_fixed()+
  labs(x="Disturbance Duration",y="Plasticity")

ggplot(filter(df_out2,Dm_seq==1))+
  geom_tile(aes(x=Dm_ord,y=h_ord,fill=corr_theta))+
  coord_fixed()+
  labs(x="Disturbance Severity",y="Local Adaptation")
ggplot(filter(df_out2,De_seq==1))+
  geom_tile(aes(x=De_ord,y=h_ord,fill=corr_theta))+
  coord_fixed()+
  labs(x="Disturbance Extent",y="Local Adaptation")
ggplot(filter(df_out2,Dl_seq==1))+
  geom_tile(aes(x=Dl_ord,y=h_ord,fill=corr_theta))+
  coord_fixed()+
  labs(x="Disturbance Duration",y="Local Adaptation")


##### Map #####
maplist=vector(mode="list",length=nrow(sim_index_sub))
for(i in 1:nrow(sim_index_sub)){
  simID_i <- sim_index_sub$sim_id[i]
  experiment_i <- filter(sim_index,sim_id==simID_i)
  qrast <- rast(paste0(experiment_folder,"/Maps/b",experiment_i$basemap_id,"/qmap_b",
                       experiment_i$basemap_id,"_q",experiment_i$qmap_id,".tif"))
  load(paste0(experiment_folder,"/Maps/habfiles/hab_",experiment_i$hab_id,".RData"))
  maplist[[i]] <- ggplot()+
    ggspatial::layer_spatial(qrast$q)+ #aggregate(qrast$q,2): doesn't work great for kimbe, but okay for others
    scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
    #geom_sf(data=hab_params$sfc_patches,size=0.5,color='red')+ # include anemones
    ggspatial::annotation_scale()+
#    labs(title=paste0("Sim ", simID_i, "\n",gsub("Export full, ","",experiment_i$notes)))+
    labs(title=paste0("Sim ", simID_i, "\nAutocorr=",round(experiment_i$q_autocorr_scale)))+
    theme(legend.position = "none",axis.text=element_blank(),axis.ticks = element_blank(),
          plot.title = element_text(hjust=0.5))
}

grid.arrange(grobs=maplist)


library(colorspace)
for(i in 1:nrow(sim_index_sub)){
  simID_i <- sim_index_sub$sim_id[i]
  experiment_i <- filter(sim_index,sim_id==simID_i)
  
  # qmap
  qrast <- rast(paste0(experiment_folder,"/Maps/b",experiment_i$basemap_id,"/qmap_b",
                       experiment_i$basemap_id,"_q",experiment_i$qmap_id,".tif"))
  load(paste0(experiment_folder,"/Maps/habfiles/hab_",experiment_i$hab_id,".RData"))
  qmap <- ggplot()+
    ggspatial::layer_spatial(qrast$q)+ #aggregate(qrast$q,2): doesn't work great for kimbe, but okay for others
    scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
    #geom_sf(data=hab_params$sfc_patches,size=0.5,color='red')+ # include anemones
    ggspatial::annotation_scale()+
    #    labs(title=paste0("Sim ", simID_i, "\n",gsub("Export full, ","",experiment_i$notes)))+
    theme(legend.position = "bottom",axis.text=element_blank(),axis.ticks = element_blank(),
          plot.title = element_text(hjust=0.5))
  
  load(paste0(experiment_folder,"/output/",experiment_i$sim_id,"_popsnapshot.RData"))
  pop_df <- left_join(pop_df,patch_locations,by=join_by("patch"=="id"))
  pop_df$eff_theta <- pop_df$theta*exp(2*pop_df$p*(pop_df$q-0.5))
  pop_df$eff_theta <- pmax(pmin(pop_df$eff_theta,140),0.005)
  
  pmap <- ggplot(pop_df,aes(x=x,y=y,color=p))+
    geom_point(size=0.2,alpha=0.5)+
    labs(title=paste0("t = ",t_i),x=NULL,y=NULL)+
    coord_fixed()+
    theme(legend.position="top",axis.text=element_blank(),axis.ticks = element_blank())+
    scale_color_continuous_diverging(palette="Purple-Brown",mid=0)
  
  thetamap <- ggplot(pop_df,aes(x=x,y=y,color=theta))+
    geom_point(size=0.2,alpha=0.5)+
    labs(title=paste0("t = ",t_i),x=NULL,y=NULL)+
    coord_fixed()+
    theme(legend.position="top",axis.text=element_blank(),axis.ticks = element_blank())+
    scale_color_continuous(palette="Blues")
  
  effthetamap <- ggplot(pop_df,aes(x=x,y=y,color=eff_theta))+
    geom_point(size=0.2,alpha=0.5)+
    labs(title=paste0("t = ",t_i),x=NULL,y=NULL)+
    coord_fixed()+
    theme(legend.position="top",axis.text=element_blank(),axis.ticks = element_blank())+
    scale_color_continuous(palette="Oranges",trans="reverse")+
    guides(color = guide_colorbar(reverse = TRUE))
  

  q_v_theta <- ggplot(pop_df,aes(x=q,y=theta))+
    geom_point(size=0.2,alpha=0.5)+
    geom_smooth()+
    labs(title=paste0("Corr(q, theta) = ",round(cor(pop_df$q,pop_df$theta),2)))
  q_v_efftheta <- ggplot(pop_df,aes(x=q,y=eff_theta))+
    geom_point(size=0.2,alpha=0.5)+
    geom_smooth()+
    labs(title=paste0("Corr(q, efftheta) = ",round(cor(pop_df$q,pop_df$eff_theta),2)))
  q_v_p <- ggplot(pop_df,aes(x=q,y=p))+
    geom_point(size=0.2,alpha=0.5)+
    geom_smooth()+
    labs(title=paste0("Corr(q, p) = ",round(cor(pop_df$q,pop_df$p),2)))
  
  plot_to_save <- grid.arrange(qmap,pmap,thetamap,effthetamap,q_v_p,q_v_theta,q_v_efftheta,
               layout_matrix=rbind(c(1,1,2,3,4),
                                   c(1,1,2,3,4),
                                   c(1,1,5,6,7)),
               top=paste0("Sim ",simID_i, ", Autocorr=",round(experiment_i$q_autocorr_scale),
                          " m, Ds=",1-experiment_i$Dm,
                          ", Dl=",experiment_i$Dl,
                          ", De=",experiment_i$De)
               )
  ggsave(filename=paste0(experiment_folder,"/plots/",simID_i,"_laststep.png"),plot = plot_to_save,width = 35.2,height=16,units="cm")
}

##### Single-timestep detail #####
simID_i <- sim_index$sim_id[sims_to_pick[5]]
simID_i <- 56944889
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
simID <- 56944889
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
  labs(y="Theta\n(median,5-95%)",x=NULL)

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
