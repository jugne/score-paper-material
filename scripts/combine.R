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
post_ess3 =data.frame(run=character(),group=integer(), ess=double())
burnIn = 10

for (run in seq(1,1)){
  # system(paste0("rm -r ../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined"))
  # system(paste0("mkdir ../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined"))

  trees <- list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output"), pattern="\\.typed.network.trees$", full.names = TRUE)
  logs <- list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output"), pattern="h5n1_rep[0-9]*\\.log$", full.names = TRUE)

  logs2 <- list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max"), pattern="h5n1_rep[0-9]*\\.reNum.rootFilter.log$", full.names = TRUE)
  logs3 <- list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max"), pattern="h5n1_rep[0-9]*\\.reNum.rootFilter.log$", full.names = TRUE)


    in_command <- paste0(" -b ",burnIn," -log")
    for (j in seq(1,3)){
      in_command = paste(in_command, trees[j])
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

      t3 <- read.table(logs3[[j]], header=TRUE, sep="\t")
      # take a burnin
      t3 <- t3[-seq(1,ceiling(length(t3$posterior)/burnIn)), ]
      # calculate ess values for posterior
      ess3 <- effectiveSize(mcmc(t3$posterior))
      post_ess3 <- rbind(post_ess3, data.frame(run=paste0("run_", run, j), group=run, ess=as.numeric(ess3)))
      # post_ess = c(post_ess, as.numeric(ess["posterior"]))
    }

    out_command = gsub("rep0_", "", trees[1])
    out_command = gsub("output", "combined", out_command)

    combined_command = gsub(".trees",".trees", out_command)
    combined_command = paste(" -o ", gsub("_rep0", "_combined",combined_command), sep="")

    # system(paste("/Applications/BEAST\\ 2.6.2/bin/logcombiner", in_command, combined_command, "", sep=" "))
    # system(paste("/Applications/BEAST\\ 2.6.2/bin/logcombiner", gsub(".typed.network.trees",".log", in_command), gsub(".typed.network.trees",".log", combined_command), sep=" "))





    # t <- read.table(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.log"), header=TRUE, sep="\t")
    # ess <- effectiveSize(mcmc(t$posterior))
    # post_ess <- rbind(post_ess, data.frame(run=paste0("run_", run, 1), group=run, ess=as.numeric(ess)))
    # post_ess <- rbind(post_ess, data.frame(run=paste0("run_", run, 2), group=run, ess=as.numeric(ess)))
    # post_ess <- rbind(post_ess, data.frame(run=paste0("run_", run, 3), group=run, ess=as.numeric(ess)))
    # 
    # ess_hist <- ggplot(post_ess, aes(x=run, y=ess, fill=group)) +
    #   geom_bar(stat="identity", position ="identity", alpha=0.5, width=1) +
    #   labs(fill="run", x = "runs and repeats") + labs(y = "ESS for posterior") +
    #   theme(legend.position = "none", axis.ticks = element_blank()) +
    #   scale_x_discrete(labels=  c("", "run_1", "", "", "run_2", "", "", "run_3", "",
    #                               "", "run_4", "",
    #                               "", "run_5", "",
    #                               "", "run_6", "",
    #                               "", "run_7", "",
    #                               "", "run_8", "",
    #                               "", "run_9", "",
    #                               "", "run_10", ""))
    # ess_hist <-ess_hist + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 13)) + scale_y_continuous(limits = c(0,4000))


    ## network root filtering

    path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_', run, "/")
    # setwd(paste0(path,"output"))
    out <- paste0(path,"output","/root_filter_max")
    system(paste0("rm -r ",out,"/combined"))
    system(paste0("mkdir ",out,"/combined"))

    files <- list.files(path=paste0(path, "output"), pattern="*.typed.network.trees")
    s <- ""
    ss <- ""
    for (file in files){
      s <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".typed.network.trees", ".rootFilter.log")))
      ss <- paste0(ss, " -log ", paste0(out, "/", str_replace(file, ".typed.network.trees", ".typed.network.rootFilter.trees")))
    }
    x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', s,' -o ',
                paste0(out, "/", str_replace(files[3], "rep2.typed.network.trees", "combined.rootFilter.log")), ' -b 10')
    # system(x)

    x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', ss,' -o ',
                paste0(out, "/", str_replace(files[1], "rep0.typed.network.trees", "combined.rootFilter.trees")), ' -b 10')
    # system(x)

    # t2 <- read.table(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.log"), header=TRUE, sep="\t")
    # ess2 <- effectiveSize(mcmc(t2$posterior))
    # post_ess2 <- rbind(post_ess2, data.frame(run=paste0("run_", run, 1), group=run, ess=as.numeric(ess2)))
    # post_ess2 <- rbind(post_ess2, data.frame(run=paste0("run_", run, 2), group=run, ess=as.numeric(ess2)))
    # post_ess2 <- rbind(post_ess2, data.frame(run=paste0("run_", run, 3), group=run, ess=as.numeric(ess2)))
    # 
    # ess_hist2 <- ggplot(post_ess2, aes(x=run, y=ess, fill=group)) +
    #   geom_bar(stat="identity", position ="identity", alpha=0.5, width=1) +
    #   labs(fill="run", x = "runs and repeats") + labs(y = "ESS for posterior") +
    #   theme(legend.position = "none", axis.ticks = element_blank()) +
    #   scale_x_discrete(labels=  c("", "run_1", "", "", "run_2", "", "", "run_3", "",
    #                               "", "run_4", "",
    #                               "", "run_5", "",
    #                               "", "run_6", "",
    #                               "", "run_7", "",
    #                               "", "run_8", "",
    #                               "", "run_9", "",
    #                               "", "run_10", ""))
    # ess_hist2 <-ess_hist2 + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 13)) + scale_y_continuous(limits = c(0,4000))



    ## filtering by nodes corresponding to segment roots

    path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_', run, "/")
    # setwd(paste0(path,"output"))
    out <- paste0(path,"output","/segment_root_filter_max")
    system(paste0("rm -r ",out,"/combined"))
    system(paste0("mkdir ",out,"/combined"))
    files <- list.files(path=paste0(path, "output"), pattern="*.typed.network.trees")
    s <- ""
    ss <- ""
    for (file in files){
      s <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".typed.network.trees", ".rootFilter.log")))
      ss <- paste0(ss, " -log ", paste0(out, "/", str_replace(file, ".typed.network.trees", ".typed.network.rootFilter.trees")))
    }
    x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', s,' -o ',
                paste0(out, "/", str_replace(files[3], "rep2.typed.network.trees", "combined.rootFilter.log")), ' -b 10')
    # system(x)

    x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', ss,' -o ',
                paste0(out, "/", str_replace(files[1], "rep0.typed.network.trees", "combined.rootFilter.trees")), ' -b 10')
    system(x)


    # t3 <- read.table(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max/h5n1_combined.rootFilter.log"), header=TRUE, sep="\t")
    # ess3 <- effectiveSize(mcmc(t3$posterior))
    # post_ess3 <- rbind(post_ess3, data.frame(run=paste0("run_", run, 1), group=run, ess=as.numeric(ess3)))
    # post_ess3 <- rbind(post_ess3, data.frame(run=paste0("run_", run, 2), group=run, ess=as.numeric(ess3)))
    # post_ess3 <- rbind(post_ess3, data.frame(run=paste0("run_", run, 3), group=run, ess=as.numeric(ess3)))
    # 
    # ess_hist3 <- ggplot(post_ess3, aes(x=run, y=ess, fill=group)) +
    #   geom_bar(stat="identity", position ="identity", alpha=0.5, width=1) +
    #   labs(fill="run", x = "runs and repeats") + labs(y = "ESS for posterior") +
    #   theme(legend.position = "none", axis.ticks = element_blank()) +
    #   scale_x_discrete(labels=  c("", "run_1", "", "", "run_2", "", "", "run_3", "",
    #                               "", "run_4", "",
    #                               "", "run_5", "",
    #                               "", "run_6", "",
    #                               "", "run_7", "",
    #                               "", "run_8", "",
    #                               "", "run_9", "",
    #                               "", "run_10", ""))
    # ess_hist3 <-ess_hist3 + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 13)) + scale_y_continuous(limits = c(0,4000))
    

    ## Execute only after combined is done and matlab scripts executed to make simulation xmls
    # 
    # file.copy(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.log"), paste0("../FullDatesInference/simulation/posterior_rates/no_filter/run_",run,"/"), overwrite = TRUE )
    # # system(paste("/Applications/BEAST\\ 2.6.2/bin/logcombiner -b 0 -resample 100000 -log h5n1.log -o h5n1_comb.log"))
    # setwd(paste0("../FullDatesInference/simulation/posterior_rates/no_filter/run_",run))
    # system(paste0("java -jar /Users/jugne/Documents/Source/score_16072020.jar -seed random -overwrite ","h5n1_combined.sim.xml"))
    # setwd(this.dir)
    # 
    # file.copy(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.log"), paste0("../FullDatesInference/simulation/posterior_rates/root_filter_max/run_",run,"/h5n1_combined.log"), overwrite = TRUE )
    # # system(paste("/Applications/BEAST\\ 2.6.2/bin/logcombiner -b 0 -resample 100000 -log h5n1.log -o h5n1_comb.log"))
    # system(paste0("mkdir -p ../FullDatesInference/simulation/posterior_rates/root_filter_max/run_",run))
    # setwd(paste0("../FullDatesInference/simulation/posterior_rates/root_filter_max/run_",run))
    # system(paste0("java -jar /Users/jugne/Documents/Source/score_16072020.jar -seed random -overwrite ","h5n1_combined.rootFilter.sim.xml"))
    # setwd(this.dir)
    # 
    # file.copy(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max/h5n1_combined.rootFilter.log"), paste0("../FullDatesInference/simulation/posterior_rates/segment_root_filter_max/run_",run,"/h5n1_combined.log"), overwrite = TRUE )
    # # system(paste("/Applications/BEAST\\ 2.6.2/bin/logcombiner -b 0 -resample 100000 -log h5n1.log -o h5n1_comb.log"))
    # system(paste0("mkdir -p ../FullDatesInference/simulation/posterior_rates/segment_root_filter_max/run_",run))
    # setwd(paste0("../FullDatesInference/simulation/posterior_rates/segment_root_filter_max/run_",run))
    # system(paste0("java -jar /Users/jugne/Documents/Source/score_16072020.jar -seed random -overwrite ","h5n1_combined.rootFilter.sim.xml"))
    # setwd(this.dir)
    
    # l1 <- read.table(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.log"), header=TRUE, sep="\t")
    # str<- paste0("Maximum posterior sample nr: ",l1$Sample[which(l1$posterior==max(l1$posterior))])
    # write(str, file=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined/maxPosterior.txt"))
    # 
    # l1 <- read.table(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.log"), header=TRUE, sep="\t")
    # str<- paste0("Maximum posterior sample nr: ",l1$Sample[which(l1$posterior==max(l1$posterior))])
    # write(str, file=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/maxPosterior.txt"))
    # 
    # l1 <- read.table(paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max/h5n1_combined.rootFilter.log"), header=TRUE, sep="\t")
    # str<- paste0("Maximum posterior sample nr: ",l1$Sample[which(l1$posterior==max(l1$posterior))])
    # write(str, file=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max/maxPosterior.txt"))
}

# ggsave("posterior_ess_full.pdf",ess_hist, width = 15)
# ggsave("posterior_ess_root_filtering.pdf",ess_hist2, width = 15)
# ggsave("posterior_ess_segment_root_filtering.pdf",ess_hist3, width = 15)

