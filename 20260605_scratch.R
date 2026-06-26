source('0_Setup.R')
#library(profvis)
experiment_name <- "Exp4_20260611"
experiment_folder <- paste0("experiments/",experiment_name)
experiment_index <- read.csv(paste0(experiment_folder,"/experiment_index.csv"))
load(paste0(experiment_folder,"/basemap_index.RData")) # basemap index

# identify which run we're on (as indexed by sim_design)
run_i <- 3
run_info <- experiment_index[run_i,]

# Identify folders in the project directory and MSI scratch directory
temp_dir <- experiment_name
#temp_dir <- paste0("connmats_smallerkern2_",experiment_name)
temp_path <- file.path("/scratch.global","schla103",temp_dir)
connmat_folder=paste0(temp_path,"/map_",run_info$mapID)

# params
load(file=paste0(experiment_folder,"/params_",run_info$param_id,".RData")) 
list2env(x=params,envir=environment())

# hab_params
load(file=paste0(experiment_folder,"/habfiles/habparams_",run_info$mapID,".RData"))

# group index
group_index <- expand.grid(alpha=1:length(v_alphas),theta=1:length(v_thetas),p=1:length(v_p))
ngroups <- nrow(group_index)

v_repr_in <- matrix(nrow=hab_params$npatch,ncol=ngroups)

group_index$repr_mean <- NA
group_index$repr_var <- NA
for(g in 36:42){
  conn_mat <- collapse::qM(read_fst(paste0(connmat_folder,"/grp_",g)))
  v_repr_in[,g] <- rowSums(conn_mat)
  
  v_repr_out <- colSums(conn_mat)
  # conn_mat <- sweep(conn_mat,2,v_repr_out,FUN='/')
  # v_repr_out <- colSums(conn_mat)
  group_index$repr_mean[g] <- mean(v_repr_out,na.rm=TRUE)
  group_index$repr_var[g] <- var(v_repr_out,na.rm=TRUE)
  rm(conn_mat)
  gc()
}

df_repr_in <- data.frame(b=hab_params$patch_locations$b) |> cbind(v_repr_in)

for(g in 36:42){
  plotcorr <- plot(hab_params$patch_locations$b,v_repr_in[,g],main=paste(g))
  print(plotcorr)
}

plot1 <- ggplot(hab_params$patch_locations,aes(x=x,y=y,color=v_repr_in[,36]))+
  geom_point()
plot2 <- ggplot(hab_params$patch_locations,aes(x=x,y=y,color=b))+
  geom_point()
grid.arrange(plot1,plot2,nrow=1)


ggplot(group_index[36:42,],aes(x=1000*v_thetas[theta]))+
  geom_point(aes(y=repr_mean))+
  geom_errorbar(aes(ymin=repr_mean-sqrt(repr_var),ymax=repr_mean+sqrt(repr_var)))+
  scale_x_log10()+
  labs(x="Dispersal kernel mean (m)", y=" Dispersal output (mean +/- sd across patches)")

plot(group_index$theta,group_index$repr_mean)

g=31
conn_mat <- collapse::qM(read_fst(paste0(connmat_folder,"/grp_",g)))
conn_mat[1:10,1:10]
hist(colSums(conn_mat))
mean(colSums(conn_mat))
sd(colSums(conn_mat),na.rm=TRUE)
sum(conn_mat,na.rm=TRUE)

g=36
conn_mat <- collapse::qM(read_fst(paste0(connmat_folder,"/grp_",g)))
conn_mat[1:10,1:10]
hist(colSums(conn_mat))
mean(colSums(conn_mat),na.rm=TRUE)
sd(colSums(conn_mat),na.rm=TRUE)
sum(conn_mat,na.rm=TRUE)
