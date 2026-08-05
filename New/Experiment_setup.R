.libPaths("/projects/standard/mrunj/shared/Rlibs_schla103")
source("New/functions.R")
experiment_folder <- "New/Maps"

# Want to run some sims on the medium Kimbe map at low density. That's b4, p6
# We'll use same habitat quality everywhere. That's q13. (Could also use q12, or make some more, but they're slow.)
f_MakeHabitat(qmapID=13,popmapID = 6, experiment_folder = experiment_folder)

# Then run MakeRefMat for popmapID=6. 