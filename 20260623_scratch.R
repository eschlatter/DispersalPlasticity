source('0_Setup.R')
experiment_index <- read.csv("experiments/Exp3_20260528/experiment_index_disturbset.csv")

i <- 2

mapID <- experiment_index$mapID[i]
load(paste0("experiments/Exp3_20260528/habfiles/habparams_",mapID,".RData"))
patch_abunds <- read.csv(file=paste0("experiments/Exp3_20260528/output/no_disturb_b01/rep1/map_",mapID,"_patch_abunds.csv"),header=TRUE)
model_df <- data.frame(v_b=hab_params$patch_locations$b,patch_abunds=NA)
model_out <- data.frame(t_i=integer(),intercept=numeric(),slope=numeric())

for(t_i in 1:nrow(patch_abunds)){
  
  model_df$patch_abunds <- t(patch_abunds[t_i,])
  lm_i <- lm(patch_abunds~v_b,data = model_df)
  model_out[t_i,] <- data.frame(t_i=t_i,intercept=lm_i$coefficients[1],slope=lm_i$coefficients[2])
}

g1 <- ggplot(model_out,aes(x=t_i,y=slope))+geom_line()+labs(title=paste0("Map ",mapID),y="slope of number of arriving larvae vs habitat quality")
g2 <- ggplot(model_df,aes(x=v_b,y=patch_abunds))+geom_point()+geom_smooth(method="lm")+
  labs(x="Habitat quality",y="Number of larvae arriving at anemone",title=paste0("Map ",mapID,", t=",t_i))

qID <- experiment_index$q_ID[i]
qrast <- rast(paste0("experiments/Exp3_20260528/b1/set1/qmap_b1_q",qID,".tif"))
qrast_plot <- ggplot()+
  ggspatial::layer_spatial(qrast$q)+
  scale_fill_continuous(palette = 'Greens',name="q",na.value = "#d2f2f7")+
  #geom_sf(data=hab_params$sfc_patches,size=0.05)+
  annotation_scale()+
  theme(legend.position = "bottom",axis.text=element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust=0.5))


grid.arrange(g2,g1,qrast_plot,nrow=1)

patch_abunds_adult <- read.csv(file=paste0("experiments/Exp3_20260528/output/no_disturb_b01/rep1/map_",mapID,"_patch_abunds_adult.csv"),header=TRUE)
avg_patch_abunds_adult <- colSums(patch_abunds_adult)/nrow(patch_abunds_adult)
df_plot <- data.frame(b=hab_params$patch_locations$b/10,adults=avg_patch_abunds_adult)
ggplot(df_plot,aes(x=b,y=adults))+geom_point()
