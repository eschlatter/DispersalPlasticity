source('0_Setup.R')

## load data
basemap_file="seascapes/2026_03_30/5x5km_res=10m"
popmap_file="pop_density800"
load(paste0(basemap_file,"/",popmap_file,".RData"))
dists_mat <- drop_units(patch_dists)

# ## first, create the file to hold output values. (Just do this once -- not each time through the array.)
# df_all <- data.frame(rep_i=numeric(),q_autocorr=numeric(),range=numeric(),sill=numeric(),SSErr=numeric(),
#                      model=factor())
# write.table(df_all, file = paste0(basemap_file,"/set2/df_all.csv"), sep = ",", append = FALSE,
#             quote = FALSE, col.names = TRUE, row.names = FALSE)

## get params
v_autocorrs <- c(-1,0,0.07,0.37,1.99)
index_runs <- expand.grid(q_autocorr=v_autocorrs,rep_i=1:100)
run_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
q_autocorr <- index_runs$q_autocorr[run_i]
rep_i <- index_runs$rep_i[run_i]

target_dist='identity'

# generate a qmap
qmap_file <- paste0("set2/q",q_autocorr,"_sim",rep_i) # save the SpatRaster for later
q_rast <- f_GenerateHabQual(base_rast=basemap_file,q_autocorr=q_autocorr,target_dist=target_dist,
                            plot_flag=FALSE,qmap_file = qmap_file)$q_rast

# get variogram info
spdf1 <- as_Spatial(sfc_patches)
spdf1$q <- terra::extract(q_rast$q,vect(sfc_patches),xy=TRUE,search_radius=500)$q
# empirical variogram
vgm1 <- variogram(q~1,data=spdf1,cressie=TRUE,width=10,cutoff=5000) # choose a relevant bin width and cutoff here
# run exponential, gaussian and spherical models, and pick the best one (lowest std error)
vgmf <- fit.variogram(vgm1,vgm(c("Gau","Sph","Exp")))

# store output
df_all <- data.frame(rep_i=rep_i, q_autocorr=q_autocorr,range=vgmf$range[2],sill=vgmf$psill[2],SSErr=attr(vgmf,"SSErr"),
                     model=vgmf$model[2])
write.table(df_all, file = paste0(basemap_file,"/set2/df_all.csv"), sep = ",", append = TRUE, 
            quote = FALSE, col.names = FALSE, row.names = FALSE)
