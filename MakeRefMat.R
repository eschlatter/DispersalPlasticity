.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
#library(profvis)
library(foreach)
library(doParallel)

#### PARAMETERS ####
patchdist_file <- "experiments/Exp4_20260611/b3/patchdists_b3_p3.RData"
output_file="output/TryNew_1/RefMat_10k"
numCores=parallelly::availableCores()
print(numCores)
npatch=10000

#### SET UP OBJECTS ####

# get patch dists
load(patchdist_file)
patch_dists <- collapse::qM(patch_dists) # remove units, but it's in km
patch_dists <- patch_dists[1:npatch,1:npatch]

# make a reference matrix of distance bins, theta values, where each site is the kernel evaluated at the distance
# (this is fast and small)
v_dist_bins <- c(seq(from=0,to=1,by=0.001),seq(from=1,to=max(patch_dists)+0.01,by=0.01))
v_theta_bins <- c(seq(from=0.005,to=0.1,by=0.005),seq(from=0.1,to=max(patch_dists)/2,by=0.1))
ref_mat <- matrix(nrow=length(v_dist_bins),ncol=length(v_theta_bins))
for(th_i in seq_along(v_theta_bins)){
  ref_mat[,th_i] <- dgamma(v_dist_bins,shape=1,scale=v_theta_bins[th_i])
}

# figure out which row of ref_mat each distance corresponds to
# (returns the lower bound distance)
# this takes some time, and is the same size as patch_dists.
patch_dists <- asplit(patch_dists,2) # split patch_dists into a list of columns
gc()
cluster <- makeForkCluster(nnodes=numCores)
registerDoParallel(cluster)
patch_dists_ref <- foreach(col_i = patch_dists,.combine="cbind") %dopar% {
  vapply(col_i,function(cell_i) match(TRUE,cell_i<v_dist_bins)-1,FUN.VALUE=double(1),USE.NAMES = FALSE)
  # gc()
}
stopCluster(cl = cluster)
storage.mode(patch_dists_ref) <- "integer"
colnames(patch_dists_ref) <- NULL
rownames(patch_dists_ref) <- NULL
save(patch_dists_ref,ref_mat,v_dist_bins,v_theta_bins,file=paste0(output_file,".RData"))