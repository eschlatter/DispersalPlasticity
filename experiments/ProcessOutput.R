
## get the names of files containing the data
filenames <- data.frame(name=list.files(path="output"))%>%
  filter(grepl(".RData",name))
files <- filenames$name

## set up storage structures
list_p <- list()
list_abund <- list()
list_fund <- list()
list_eff <- list()
run_times <- data.frame(user.self=numeric(),sys.self=numeric(),elapsed=numeric(),user.child=numeric(),sys.child=numeric(),run=integer())

## grab each file, put its data into the storage structures
for(file_i in 1:length(files)){
  load(paste0('output/',files[file_i]))
  runnum <- base::strsplit(files[file_i],"_")[[1]][1] %>% as.integer()
  list_p[[file_i]] <- cbind(sim_out$output_list$df_p,run=runnum)
  list_abund[[file_i]] <- cbind(sim_out$output_list$df_abund,run=runnum)
  list_fund[[file_i]] <- cbind(sim_out$output_list$df_fund,run=runnum)
  list_eff[[file_i]] <- cbind(sim_out$output_list$df_eff,run=runnum)
  run_times <- rbind(run_times,cbind(t(sim_out$time_run),run=runnum)) # this isn't working
}

df_p_all <- do.call(rbind,list_p)
df_abund_all <- do.call(rbind,list_abund)
df_fund_all <- do.call(rbind,list_fund)
df_eff_all <- do.call(rbind,list_eff)

df_p_all <- df_p_all[!is.na(df_p_all$mean),]
