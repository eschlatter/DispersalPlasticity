# Create the target habitat quality distributions
# A: uniform
set_A <- runif(1000,min=1,max=9)
target_dist <- set_A
save(target_dist,file='data/target_dists/dist_A.RData')
# B: normal
set_B <- rnorm(1000,mean=5,sd=1)
target_dist <- set_B
save(target_dist,file='data/target_dists/dist_B.RData')
# C: bimodal
set_C <- c(rnorm(1000,mean=2,sd=0.5),rnorm(1000,mean=8,sd=0.5))
target_dist <- set_C
save(target_dist,file='data/target_dists/dist_C.RData')
# D: bad habitat common
set_D <- rgamma(1000,shape=2,scale=1)
set_D <- set_D[set_D<9 & set_D>=1]
target_dist <- set_D
save(target_dist,file='data/target_dists/dist_D.RData')
# E: good habitat common
# let's just flip D around
set_E <- 10-set_D
target_dist <- set_E
save(target_dist,file='data/target_dists/dist_E.RData')



# 
# 
# starting_dist <- rnorm(10000,15)
# hist(f_TransformDist(starting_dist,set_E))
# 
# 
# df_dists <- data.frame(original=rnorm(1000,mean=-0.2,sd=3))
# ggplot(df_dists,aes(x=original_cdf))+geom_density()
# 
# # get cdf
# theor <- pnorm(df_dists$original,mean=-0.2,sd=3)
# # or, if we don't have a theoretical cdf, the ordinal value/total number of values should do a pretty good job
# original_ord <- as.numeric(as.factor(df_dists$original)) # get ordinal values for each of the originals
# emp <- original_ord/length(original_ord)
# plot(emp,theor) # yep!
# df_dists$original_cdf <- emp
# 
# # if we have a theoretical target function
# df_dists$transformed_data <- qexp(emp,rate=5)
# 
# # if we have an empirical target function
# dist_dataset <- c(rnorm(1000,mean=1,sd=1),rnorm(1000,mean=10,sd=1))
# #target_fn <- ecdf(dist_dataset)
# 
# df_dists$transformed_data <- as.numeric(quantile(dist_dataset,df_dists$original_cdf))
# ggplot(df_dists,aes(x=transformed_data))+geom_density()
# hist(df_dists$transformed_data)
# hist(dist_dataset)
