source('functions/neutral_landscape.R')
source('functions/plotting_functions.R')

## fractal landscape, 17x17 (~87 patches)
K_neutral=5
amt_covered <- 0.3
testmap=fracland(k=4,h=0.6,p=1-amt_covered,binary=FALSE,plotflag=FALSE)
nx=ncol(testmap)
ny=nrow(testmap)
testmap[testmap<quantile(testmap,1-amt_covered)] <- NA
patch_locations <- as.data.frame(which(!is.na(testmap),arr.ind=TRUE)) %>%
  rename(y=row,x=col) %>%
  rowid_to_column(var='id')
patch_locations <- mutate(patch_locations,K_i=testmap[y+(x-1)*nrow(testmap)])
patch_locations$K_i <- round((K_neutral/3)*scale(patch_locations$K_i)+K_neutral)
f_Plot_Landscape(patch_locations,nx,ny)

## fractal landscape, 33x33 (~325 patches)
K_neutral=5
amt_covered <- 0.3
testmap=fracland(k=5,h=0.5,p=1-amt_covered,binary=FALSE,plotflag=FALSE)
nx=ncol(testmap)
ny=nrow(testmap)
testmap[testmap<quantile(testmap,1-amt_covered)] <- NA
patch_locations <- as.data.frame(which(!is.na(testmap),arr.ind=TRUE)) %>%
  rename(y=row,x=col) %>%
  rowid_to_column(var='id')
patch_locations <- mutate(patch_locations,K_i=testmap[y+(x-1)*nrow(testmap)])
patch_locations$K_i <- round((K_neutral/3)*scale(patch_locations$K_i)+K_neutral)
K <- patch_locations$K_i

f_Plot_Landscape(patch_locations,nx,ny) # plot it
save(patch_locations,K,nx,ny,file='seascapes/Fractal_33x33_n325.RData')
