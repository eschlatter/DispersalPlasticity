.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
source("New/functions.R")
library(terra)
library(sf)
library(data.table)
experiment_folder <- "New/Maps"

# Want to run some sims on the medium Kimbe map at low density. That's b4, p6
# We'll use same habitat quality everywhere. That's q13. (Could also use q12, or make some more, but they're slow.)
f_MakeHabitat(qmapID=13,popmapID = 6, experiment_folder = experiment_folder)
# Then run MakeRefMat for popmapID=6.

# Also want to run some sims on the full-density square map. That's b1, p2.
# Various q options: 1-7
for(q_i in 1:7){
  f_MakeHabitat(qmapID=q_i,popmapID = 2,experiment_folder = experiment_folder)
}
f_MakeHabitat(qmapID=14,popmapID = 2,experiment_folder = experiment_folder)
f_MakeHabitat(qmapID=14,popmapID = 1,experiment_folder = experiment_folder)
# Then run MakeRefMat for popmapID=2.

# Okay, that full-density map will be good later, but it's taking way too long on an interactive node.
# Let's get b1, p1. All the same q's.
for(q_i in 1:7){
  f_MakeHabitat(qmapID=q_i,popmapID = 1,experiment_folder = experiment_folder)
}
# Then run MakeRefMat for popmapID=1.