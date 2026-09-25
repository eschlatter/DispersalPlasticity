v_phis=c(0.1,0.5,0.8,1,1.2,2,5,10)
v_phis <- seq(from=0.4,by=0.05,to=1)

for(phi_i in 1:length(v_phis)){
  a <- f_SimBValues(patch_locations,dists_mat,v_phis[phi_i],show_plot=TRUE)
  print(a$sp_aut_dist)
}
