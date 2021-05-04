# clear workspace
rm(list = ls())

for (run in 1:10){
  tree <- list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined"),
                     pattern="\\.typed.network.trees$", full.names = TRUE)
  system(paste0("java -jar strNetworkSummarizer.jar -burnin 0 -strategy mean ",
                tree[1], paste0(" ../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined/h5n1_mcc.typed.network.trees")))
  
  tree_rootFilter <-list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max"),
                               pattern="\\combined.rootFilter.trees$", full.names = TRUE)
  
  system(paste0("java -jar strNetworkSummarizer.jar -burnin 0 -strategy mean ",
                tree_rootFilter[1], paste0(" ../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_mcc.typed.network.trees")))
  
  
  tree_segmentRootFilter <-list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max"),
                                      pattern="\\combined.rootFilter.trees$", full.names = TRUE)
  
  system(paste0("java -jar strNetworkSummarizer.jar -burnin 0 -strategy mean ",
                tree_segmentRootFilter[1], paste0(" ../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max/h5n1_mcc.typed.network.trees")))
  
  
  
}