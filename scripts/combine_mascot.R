######################################################
######################################################
# combine the gene tree runs and run the mcc trees
# some parts adapted from N.F. Muller scripts for CoalRe
# analysis
######################################################
######################################################
library(ggplot2)
library(coda)
library(stringr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
# this.dir <- dirname(parent.frame(2)$ofile)
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)
post_ess =data.frame(run=character(),group=integer(), ess=double())
post_ess2 =data.frame(run=character(),group=integer(), ess=double())
burnIn = 10

for (run in seq(1,10)){
  # system(paste0("rm -r ../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/combined"))
  system(paste0("mkdir ../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/combined"))
  
  logs <- list.files(path=paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/output/trimmed"), pattern="h5n1_rep[0-9]*\\.log$", full.names = TRUE)
  
  logs2 <- list.files(path=paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/output/root_filter_max"), pattern="h5n1_rep[0-9]*\\.reNum.rootFilter.log$", full.names = TRUE)
  
  in_command <- paste0(" -b ",burnIn)
  for (j in seq(1,3)){
    if (run==2 & j==3)
      break
    in_command = paste(in_command, " -log ", logs[j])
    t <- read.table(logs[[j]], header=TRUE, sep="\t")
    # take a burnin
    t <- t[-seq(1,ceiling(length(t$posterior)/burnIn)), ]
    # calculate ess values for posterior
    ess <- effectiveSize(mcmc(t$posterior))
    post_ess <- rbind(post_ess, data.frame(run=paste0("run_", run, j), group=run, ess=as.numeric(ess)))
    
    t2 <- read.table(logs2[[j]], header=TRUE, sep="\t")
    # take a burnin
    t2 <- t2[-seq(1,ceiling(length(t2$posterior)/burnIn)), ]
    # calculate ess values for posterior
    ess2 <- effectiveSize(mcmc(t2$posterior))
    post_ess2 <- rbind(post_ess2, data.frame(run=paste0("run_", run, j), group=run, ess=as.numeric(ess2)))
  }
  
  out_command = gsub("rep0", "combined", logs[1])
  out_command = gsub("output/trimmed", "combined", out_command)
  
  combined_command = gsub(".log",".log", out_command)
  combined_command = paste(" -o ", gsub("_rep0", "_combined",combined_command), sep="")

  if (!file.exists(paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.log"))){
    system(paste("/Applications/BEAST\\ 2.6.2/bin/logcombiner", in_command, combined_command, sep=" "))
  }
  
  t <- read.table(paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.log"), header=TRUE, sep="\t")
  ess <- effectiveSize(mcmc(t$posterior))
  post_ess <- rbind(post_ess, data.frame(run=paste0("run_", run, 1), group=run, ess=as.numeric(ess)))
  post_ess <- rbind(post_ess, data.frame(run=paste0("run_", run, 2), group=run, ess=as.numeric(ess)))
  if (run != 2){
    post_ess <- rbind(post_ess, data.frame(run=paste0("run_", run, 3), group=run, ess=as.numeric(ess)))
  }
  
  ess_hist <- ggplot(post_ess, aes(x=run, y=ess, fill=group)) +
    geom_bar(stat="identity", position ="identity", alpha=0.5, width=1) + 
    labs(fill="run", x = "runs and repeats") + labs(y = "ESS for posterior") +
    theme(legend.position = "none", axis.ticks = element_blank()) +
    scale_x_discrete(labels=  c("", "run_1", "", "   run_2", "", "", "run_3", "",
                                "", "run_4", "",
                                "", "run_5", "",
                                "", "run_6", "",
                                "", "run_7", "",
                                "", "run_8", "",
                                "", "run_9", "",
                                "", "run_10", ""))
  ess_hist <-ess_hist + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 13))
  
  
  ## network root filtering

  path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_', run, "/")
  # setwd(paste0(path,"output"))
  out <- paste0(path,"output","/root_filter_max")
  ha_files <- list.files(path=paste0(path, "output/trimmed"), pattern="*.ha.str.trees")
  s <- ""
  ss_ha <- ""
  ss_na <- ""
  for (file in ha_files){
    if (run==2 & file==ha_files[3])
      break
    s <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".ha.str.trees", ".rootFilter.log")))
    ss_ha <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".ha.str.trees", ".ha.str.rootFilter.trees")))
    ss_na <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".na.str.trees", ".na.str.rootFilter.trees")))
  }
  x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', s,' -o ',
              paste0(out, "/", str_replace(ha_files[3], "rep2.ha.str.trees", "combined.rootFilter.log")), ' -b 10')
  if (!file.exists(paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.trees"))){
    system(x)
  }

  x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', ss_ha,' -o ',
              paste0(out, "/", str_replace(ha_files[2], "rep1.ha.str.trees", "combined.rootFilter.trees")), ' -b 10')
  if (!file.exists(paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.log"))){
    system(x)
  }

  t2 <- read.table(paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.log"), header=TRUE, sep="\t")
  ess2 <- effectiveSize(mcmc(t2$posterior))
  post_ess2 <- rbind(post_ess2, data.frame(run=paste0("run_", run, 1), group=run, ess=as.numeric(ess2)))
  post_ess2 <- rbind(post_ess2, data.frame(run=paste0("run_", run, 2), group=run, ess=as.numeric(ess2)))
  if (run != 2){
    post_ess2 <- rbind(post_ess2, data.frame(run=paste0("run_", run, 3), group=run, ess=as.numeric(ess2)))
  }


  ess_hist2 <- ggplot(post_ess2, aes(x=run, y=ess, fill=group)) +
    geom_bar(stat="identity", position ="identity", alpha=0.5, width=1) +
    labs(fill="run", x = "runs and repeats") + labs(y = "ESS for posterior") +
    theme(legend.position = "none", axis.ticks = element_blank()) +
    scale_x_discrete(labels=  c("", "run_1", "", "   run_2", "", "", "run_3", "",
                                "", "run_4", "",
                                "", "run_5", "",
                                "", "run_6", "",
                                "", "run_7", "",
                                "", "run_8", "",
                                "", "run_9", "",
                                "", "run_10", ""))
  ess_hist2 <-ess_hist2 + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 13))
}

ggsave("mascot_posterior_ess_full.pdf",ess_hist, width = 15)
ggsave("mascot_posterior_ess_root_filtering.pdf",ess_hist2, width = 15)