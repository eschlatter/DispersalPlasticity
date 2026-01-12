source('0_Setup.R')
library(data.table)

########### Load files #################

## Info on habitats
load("experiments/Array1_SimIndex.RData")
sim_index <- as.data.table(sim_index)

sim_index_order <- arrange(sim_index,Kmap) %>%
  arrange(BaseMap) %>%
  arrange(h)
sim_index_order$adds <- rep(-0.2+0.005*1:27,7)
sim_index_order <- sim_index_order[,c('Run','adds')]
sim_index <- left_join(sim_index,sim_index_order,by='Run')
#sim_index$Run <- as.character(sim_index$Run)

## get the names of files containing the data
filenames <- data.frame(name=list.files(path="output/array4"))%>%
  filter(grepl(".RData",name))
files <- filenames$name[1:189]
files_constant <- filenames$name[190:192]

## set up storage structures
list_p <- list()
list_abund <- list()
list_fund <- list()
list_eff <- list()
run_times <- data.frame(user.self=numeric(),sys.self=numeric(),elapsed=numeric(),user.child=numeric(),sys.child=numeric(),run=integer())

## get each file, put its data into the storage structures
for(file_i in 1:length(files)){
  load(paste0('output/array4/',files[file_i]))
  runnum <- base::strsplit(files[file_i],"_")[[1]][1] %>% as.integer()
  list_p[[file_i]] <- cbind(sim_out$output_list$df_p,run=runnum)
  list_abund[[file_i]] <- cbind(sim_out$output_list$df_abund,run=runnum)
  list_fund[[file_i]] <- cbind(sim_out$output_list$df_fund,run=runnum)
  list_eff[[file_i]] <- cbind(sim_out$output_list$df_eff,run=runnum)
  run_times <- rbind(run_times,cbind(t(sim_out$time_run),run=runnum))
}

for(file_i in 1:length(files_constant)){
  load(paste0('output/array4/',files_constant[file_i]))
  runnum <- base::strsplit(files[file_i],"_")[[1]][1] %>% as.integer()
  list_p[[189+file_i]] <- cbind(sim_out$output_list$df_p,run=paste0('constant',file_i))
  list_abund[[189+file_i]] <- cbind(sim_out$output_list$df_abund,run=paste0('constant',file_i))
  list_fund[[189+file_i]] <- cbind(sim_out$output_list$df_fund,run=paste0('constant',file_i))
  list_eff[[189+file_i]] <- cbind(sim_out$output_list$df_eff,run=paste0('constant',file_i))
  run_times <- rbind(run_times,cbind(t(sim_out$time_run),run=paste0('constant',file_i)))
}

########### Process data #################

# store all timesteps
df_p_all <- as.data.table(do.call(rbind,list_p))
df_abund_all <- as.data.table(do.call(rbind,list_abund))
df_fund_all <- as.data.table(do.call(rbind,list_fund))
df_eff_all <- as.data.table(do.call(rbind,list_eff))

# post-equilibrium timesteps only
cutoff_t=2500
df_p_cut <- df_p_all[df_p_all$t>cutoff_t,]
df_abund_cut <- df_abund_all[df_abund_all$t>cutoff_t,]
df_fund_cut <- df_fund_all[df_fund_all$t>cutoff_t,]
df_eff_cut <- df_eff_all[df_eff_all$t>cutoff_t,]

## Get summary stats by sim
#   p
df_p_summary <- df_p_cut[ ,.(mean=mean(mean),
                             var_among_t=var(mean),
                             var_within_t=mean(var),
                             var_within_t_var=var(var)),
                          by=run]
df_p_summary <- sim_index[df_p_summary,on=c(Run="run")]
#   fundamental kernel
df_fund_summary <- df_fund_cut[ ,.(kernmean=mean(kernmean_mean),
                                   kernmean_var_among_t=var(kernmean_mean),
                                   kernmean_var_within_t=mean(kernmean_var),
                                   kernmean_var_within_t_var=var(kernmean_var)),
                                by=run]
df_fund_summary <- sim_index[df_fund_summary,on=c(Run="run")]
#   effective kernel
df_eff_summary <- df_eff_cut[ ,.(kernmean=mean(kernmean_mean),
                                 kernmean_var_among_t=var(kernmean_mean),
                                 kernmean_var_within_t=mean(kernmean_var),
                                 kernmean_var_within_t_var=var(kernmean_var)),
                              by=run]
df_eff_summary <- sim_index[df_eff_summary,on=c(Run="run")]

save(df_abund_all,df_abund_cut,df_eff_all,df_eff_cut,df_eff_summary,
     df_fund_all,df_fund_cut,df_fund_summary,df_p_all,df_p_cut,df_p_summary,sim_index,
     file="output/array4/array4_processed.RData")

########### Output plots #################
load('output/array4/array4_processed.RData')

# h vs p
ggplot(df_p_summary,aes(x=h,group=Run))+
  geom_point(aes(y=mean,color=as.factor(BaseMap)),position=position_dodge(width=0.1))+
  geom_errorbar(aes(ymin=mean-sqrt(var_among_t),ymax=mean+sqrt(var_among_t),color=as.factor(BaseMap)),position=position_dodge(width=0.1))+
  labs(x='Habitat quality aggregation (h)',y="plasticity parameter (mean +/- sd over time)",title="Plasticity vs Habitat aggregation")+
  theme_minimal()

p_constants <- filter(df_p_summary,is.na(Kmap))
p_constants$BaseMap <- c(1,2,3)

ggplot(filter(df_p_summary,!is.na(Kmap)),aes(x=h+adds,group=Run))+
  geom_point(aes(y=mean,color=as.factor(BaseMap)))+
  geom_errorbar(aes(ymin=mean-sqrt(var_among_t),ymax=mean+sqrt(var_among_t),color=as.factor(BaseMap)))+
  geom_rect(data=p_constants,aes(xmin=-2,xmax=2,ymin=mean-sqrt(var_among_t),ymax=mean+sqrt(var_among_t),fill=as.factor(BaseMap)),alpha=0.1)+
  geom_hline(data=p_constants,aes(yintercept=mean,color=as.factor(BaseMap)),lty='dashed')+
  labs(x='Habitat quality aggregation (h)',y="plasticity parameter (mean +/- sd over time)",title="Plasticity vs Habitat aggregation")+
  theme_minimal()

# h vs fundamental and effective kernels
##  mean
ggplot(df_fund_summary,aes(x=h+adds,group=Run))+
  geom_point(aes(y=kernmean,color='Fundamental Kernel'),position=position_dodge(width=0.05))+
  geom_errorbar(aes(ymin=kernmean-sqrt(kernmean_var_among_t),ymax=kernmean+sqrt(kernmean_var_among_t),color='Fundamental Kernel'),
                position=position_dodge(width=0.05))+
  # geom_point(data=df_eff_summary,aes(y=kernmean,color='Effective Kernel'),position=position_dodge(width=0.05))+
  # geom_errorbar(data=df_eff_summary,aes(ymin=kernmean-sqrt(kernmean_var_among_t),ymax=kernmean+sqrt(kernmean_var_among_t),color='Effective Kernel'),
  #               position=position_dodge(width=0.05))+
  labs(x='Habitat quality aggregation (h)',y="dispersal kernel mean\n(mean +/- sd over time)",title="Dispersal kernel mean vs Habitat aggregation")+
  theme_minimal()
##  within-timestep variance
ggplot(df_fund_summary,aes(x=h+adds,group=Run))+
  geom_point(aes(y=kernmean_var_within_t,color='Fundamental Kernel'))+
  geom_errorbar(aes(ymin=kernmean_var_within_t-sqrt(kernmean_var_within_t_var),ymax=kernmean_var_within_t+sqrt(kernmean_var_within_t_var),color='Fundamental Kernel'))+
  geom_point(data=df_eff_summary,aes(y=kernmean_var_within_t,color='Effective Kernel'))+
  geom_errorbar(data=df_eff_summary,aes(ymin=kernmean_var_within_t-sqrt(kernmean_var_within_t_var),ymax=kernmean_var_within_t+sqrt(kernmean_var_within_t_var),color='Effective Kernel'))+
  labs(x='Habitat quality aggregation (h)',y="Within-timestep variance in kernel mean",title="Dispersal kernel variation vs Habitat aggregation")+
  theme_minimal()

########### Dynamics plots #################
# pick what runs to plot
plotruns <- filter(sim_index,Kmap==54)$Run
plotruns <- c(180)

ggplot(filter(df_abund_all,run %in% plotruns),aes(x=t,y=abund,group=as.factor(run)))+
  geom_line(aes(color=as.factor(run)))

# get subsets of the data
df_eff_selected <- filter(df_eff_all,run %in% plotruns)
df_eff_selected <- sim_index[df_eff_selected,on=c(Run="run")]
df_eff_selected$Run <- as.factor(df_eff_selected$Run)
df_eff_selected$Kmap <- as.factor(df_eff_selected$Kmap)

df_fund_selected <- filter(df_fund_all,run %in% plotruns)
df_fund_selected <- sim_index[df_fund_selected,on=c(Run="run")]
df_fund_selected$Run <- as.factor(df_fund_selected$Run)
df_fund_selected$Kmap <- as.factor(df_fund_selected$Kmap)

df_p_selected <- filter(df_p_all,run %in% plotruns)
df_p_selected <- sim_index[df_p_selected,on=c(Run="run")]
df_p_selected$Run <- as.factor(df_p_selected$Run)
df_p_selected$Kmap <- as.factor(df_p_selected$Kmap)

# plot plasticity
ggplot(df_p_selected,aes(x=t,group=Run))+
  geom_ribbon(aes(ymin=mean-sqrt(var),ymax=mean+sqrt(var),fill=Run),alpha=0.1)+
  geom_line(aes(y=mean,color=Run))+
  labs(x="time",y="Plasticity parameter\n(mean +/- sd within timesteps)",title="Simulation dynamics: Plasticity parameter")+
  theme_minimal()

# plot effective kernel
ggplot(df_eff_selected,aes(x=t,group=Run))+
  geom_ribbon(aes(ymin=kernmean_mean-sqrt(kernmean_var),ymax=kernmean_mean+sqrt(kernmean_var),fill=Run),alpha=0.1)+
  geom_line(aes(y=kernmean_mean,color=Run))+
  labs(x="time",y="Effective kernel mean\n(mean +/- sd within timesteps)",title="Simulation dynamics: Effective kernel mean")+
  theme_minimal()

ggplot(df_eff_selected,aes(x=t,group=Run))+
  geom_line(aes(y=kernmean_moran,color=Run))+
  labs(x="time",y="Moran's I of effective kernel mean",title="Simulation dynamics:\nSpatial structure of effective kernel mean")+
  theme_minimal()

# plot fundamental kernel
ggplot(df_fund_selected,aes(x=t,group=Run))+
  geom_ribbon(aes(ymin=kernmean_mean-sqrt(kernmean_var),ymax=kernmean_mean+sqrt(kernmean_var),fill=Run),alpha=0.1)+
  geom_line(aes(y=kernmean_mean,color=Run))+
  labs(x="time",y="Fundamental kernel mean\n(mean +/- sd within timesteps)",title="Simulation dynamics: Fundamental kernel mean")+
  theme_minimal()

ggplot(df_fund_selected,aes(x=t,group=Run))+
  geom_line(aes(y=kernmean_moran,color=Run))+
  labs(x="time",y="Moran's I of fundamental kernel mean",title="Simulation dynamics:\nSpatial structure of fundamental kernel mean")+
  theme_minimal()


# find Moran's I of the habitat

load(paste0("seascapes/kmaps/kmap_",21,".RData"))
npatch=nrow(kmap_i$patch_locations)
patch_dists <- matrix(nrow=npatch,ncol=npatch)
colnames(patch_dists) <- kmap_i$patch_locations$id
for(i in 1:npatch){
  patch_dists[i,] = sqrt((kmap_i$patch_locations$x[i]-kmap_i$patch_locations$x)^2+
                           (kmap_i$patch_locations$y[i]-kmap_i$patch_locations$y)^2)
}
moran_weights <- 1/(patch_dists)
diag(moran_weights) <- 0
moran_hab <- Moran.I(kmap_i$patch_locations$K_i,weight=moran_weights)$observed


########### Habitat plots #################
basemap_list <- unique(sim_index$BaseMap)
plotlist=list()

for(base_i in 1:length(basemap_list)){
  # get a Kmap with the correct base map
  kmap_ind <- min(which(sim_index$BaseMap==basemap_list[base_i]))
  load(paste0("seascapes/kmaps/kmap_",sim_index$Kmap[kmap_ind],".RData"))
  
  # make a homogenous basemap
  basemap_1 <- matrix(0,nrow=kmap_i$ny,ncol=kmap_i$nx)
  for(i in 1:nrow(kmap_i$patch_locations)) basemap_1[kmap_i$patch_locations$y[i],kmap_i$patch_locations$x[i]] <- 1
  
  # put it in Kmap format and save it
  kmap_i <- f_GenerateMapWithK(base_map=basemap_1,K_range=c(10,10),h=0,k=7,p=0.3,h_base=0.8,plot_flag = FALSE)
  save(kmap_i,file=paste0("seascapes/kmaps/kmap_constant",base_i,".RData"))
  
  # store the plot
  g <- f_Plot_Landscape(kmap_i$patch_locations,kmap_i$nx,kmap_i$ny)+
    labs(title=paste("Base Map",base_i))+
    theme(legend.position = "none")+
    coord_fixed()+
    lims(x=c(-1,66),y=c(66,-1)) # 2^k+1
  print(g)
}

hneg2 <- filter(sim_index,h==1.8) %>%
  group_by(BaseMap) %>%
  summarize(Run=first(Run))

hneg2 <- data.frame(Run=c(180))

kmap_ind <- hneg2$Run[1]
load(paste0("seascapes/kmaps/kmap_",sim_index$Kmap[kmap_ind],".RData"))
f_Plot_Landscape(kmap_i$patch_locations,kmap_i$nx,kmap_i$ny)+
  labs(title='h=1.8')+
  coord_fixed()+
  lims(x=c(-1,66),y=c(66,-1)) # 2^k+1

kmap_ind <- hneg2$Run[2]
load(paste0("seascapes/kmaps/kmap_",sim_index$Kmap[kmap_ind],".RData"))
f_Plot_Landscape(kmap_i$patch_locations,kmap_i$nx,kmap_i$ny)+
  #labs(title='h=1.8')+
  coord_fixed()+
  lims(x=c(-1,66),y=c(66,-1)) # 2^k+1

kmap_ind <- hneg2$Run[3]
load(paste0("seascapes/kmaps/kmap_",sim_index$Kmap[kmap_ind],".RData"))
f_Plot_Landscape(kmap_i$patch_locations,kmap_i$nx,kmap_i$ny)+
  labs(title='h=1.8')+
  coord_fixed()+
  lims(x=c(-1,66),y=c(66,-1)) # 2^k+1

load(paste0("seascapes/kmaps/kmap_",57,".RData"))
f_Plot_Landscape(kmap_i$patch_locations,kmap_i$nx,kmap_i$ny)+
  coord_fixed()+
  labs(title='h=0.6')+
  lims(x=c(-1,66),y=c(66,-1)) # 2^k+1
