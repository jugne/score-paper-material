library(ggplot2)
library(scales)
library(ggthemes)
library("gridExtra")
library(wesanderson)
# library(wesanderson)

# clear workspace
rm(list = ls())
# fit <- 2
fit<- 4
directory <- "root_filter_max"
for(run in seq(1,10)){
  
  # Set the directory to the directory of the file
  # this.dir <- dirname(parent.frame(2)$ofile)
  this.dir <- "/Users/jugne/Documents/SCORE-paper/scripts"
  setwd(this.dir)
  system(paste0("mkdir -p",directory,"/run_",run))
  
  # window_size <- c(0.5, 1, 1.5,2,3,4,5,10)#, 50)
  window_size <- c(0.5, 1, 1.5,2,3,4,5)
  
  dat_full <- data.frame(value = double(), position = character(), window_size = integer());
  dat_fit <- data.frame(value = double(), position = character(), window_size = integer());
  dat_unfit <- data.frame(value = double(), position = character(), window_size = integer());
  dat_prior <- data.frame(value = double(), variable = character(), window_size = integer());
  # dat_mean_rates <- data.frame(value = double(), variable = character(), window_size = integer());
  # dat_mean_rates_corr <- data.frame(value = double(), variable = character(), window_size = integer());
  dat_post_rates <- data.frame(value = double(), variable = character(), window_size = integer());
  
  diff_prior <- data.frame(value = double(), variable = integer(), method = character());
  # diff_mean_rates <- data.frame(value = double(), variable = integer(), method = character());
  # diff_mean_rates_corr <- data.frame(value = double(), variable = integer(), method = character());
  diff_post_rates <- data.frame(value = double(), variable = integer(), method = character());
  diff_full <- data.frame(value = double(), variable = integer(), method = character());
  diff_fit <- data.frame(value = double(), variable = integer(), method = character());
  diff_unfit <- data.frame(value = double(), variable = integer(), method = character());
  
  dist_diff <- data.frame(value = double(), variable = integer(), method = character());
  
  for (i in window_size) {
    if (!file.exists())
    system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ", i,
                  " -burnin 0 -minTipDistance ",fit," ../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.trees root_filter_max/run_",run,"/reaMigAnalysis_", i,"_",fit,".txt "))
    
    data = read.table(paste0("root_filter_max/run_",run,"/reaMigAnalysis_",i, "_",fit,".txt"), header=TRUE)
    
    rea_window <- data.frame(value = data$n_reassortment_window/data$network_length_window, variable = "Full On Window")
    rea_outside_window <- data.frame(value = (data$n_reassortment-data$n_reassortment_window)/
                                       (data$network_length-data$network_length_window), 
                                     variable = "Full Off Window")
    
    dat_full<- rbind(dat_full,cbind(rbind(rea_window, rea_outside_window), data.frame(window_size = i)))
    diff_full<- rbind(diff_full, data.frame(value = rea_window$value - rea_outside_window$value, variable = i, method = "SCoRe_full"))
    
    
    ### fit 
    
    rea_fit_window <- data.frame(value = data$n_reassortment_window_fit/data$network_length_window_fit, variable = "Fit On Window")
    rea_fit_off_window <- data.frame(value = data$n_reassortment_off_window_fit/data$network_length_off_window_fit, variable = "Fit Off Window")
    
    dat_fit<- rbind(dat_fit,cbind(rbind(rea_fit_window, rea_fit_off_window), data.frame(window_size = i)))
    diff_fit<- rbind(diff_fit, data.frame(value = rea_fit_window$value - rea_fit_off_window$value, variable = i, method = "SCoRe_fit"))
    
    ### unfit 
    
    rea_unfit_window <- data.frame(value = data$n_reassortment_window_unfit/data$network_length_window_unfit, variable = "Unfit On Window")
    rea_unfit_off_window <- data.frame(value = data$n_reassortment_off_window_unfit/data$network_length_off_window_unfit, variable = "Unfit Off Window")
    
    dat_unfit<- rbind(dat_unfit,cbind(rbind(rea_unfit_window, rea_unfit_off_window), data.frame(window_size = i)))
    diff_unfit<- rbind(diff_unfit, data.frame(value = rea_unfit_window$value - rea_unfit_off_window$value, variable = i, method = "SCoRe_unfit"))
    
    
    ### Prior ###
    
    # system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ",
    #               i," -burnin 0 -minTipDistance 2 /Users/jugne/Documents/SCORE-paper/TipDateInference/h5n1_prior/reject/run_",run,"/h5n1_rep0.typed.network.trees run_",run,"/reaMigAnalysis_prior_reject_",i, ".txt"))
    # 
    # data_prior = read.table(paste0("run_",run,"/reaMigAnalysis_prior_reject_",i, ".txt"), header=TRUE)
    # 
    # rea_window_prior <- data.frame(value = data_prior$n_reassortment_window/data_prior$network_length_window, variable = "On Window")
    # full_outside_window_prior <- data.frame(value = (data_prior$n_reassortment-data_prior$n_reassortment_window)/
    #                                           (data_prior$network_length-data_prior$network_length_window),
    #                                         variable = "Off Window")
    # 
    # 
    # dat_prior <- rbind(dat_prior, cbind(rbind(rea_window_prior, full_outside_window_prior), data.frame(window_size = i)))
    # diff_prior <- rbind(diff_prior, data.frame(value = rea_window_prior$value - full_outside_window_prior$value, variable = i, method = "prior"))
    
    # ### Mean Rates, No Data Simulation ###
    # 
    # # system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ",
    # #               i," -burnin 0 -minTipDistance 2 /Users/jugne/Documents/SCORE-paper/TipDateInference/simulation/h5n1_run2_combined.sim.network.trees reaMigAnalysis_mean_rates_no_data_reject_",i, ".txt"))
    # 
    # data_mean_rates = read.table(paste0("reaMigAnalysis_mean_rates_no_data_reject_",i, ".txt"), header=TRUE)
    # 
    # rea_window_mean_rates <- data.frame(value = data_mean_rates$n_reassortment_window/data_mean_rates$network_length_window, variable = "On Window")
    # full_outside_window_mean_rates <- data.frame(value = (data_mean_rates$n_reassortment-data_mean_rates$n_reassortment_window)/
    #                                           (data_mean_rates$network_length-data_mean_rates$network_length_window),
    #                                         variable = "Off Window")
    # 
    # 
    # dat_mean_rates <- rbind(dat_mean_rates, cbind(rbind(rea_window_mean_rates, full_outside_window_mean_rates), data.frame(window_size = i)))
    # diff_mean_rates <- rbind(diff_mean_rates, data.frame(value = rea_window_mean_rates$value - full_outside_window_mean_rates$value, variable = i, method = "mean rates"))
    # 
    # ### Mean Rates Forced Correlation, No Data Simulation ###
    # 
    # # system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ",
    # #               i," -burnin 0 -minTipDistance 2 /Users/jugne/Documents/SCORE-paper/TipDateInference/simulation/h5n1_run2_combined_forced_corr.sim.network.trees reaMigAnalysis_mean_rates_corr_no_data_reject_",i, ".txt"))
    # 
    # data_mean_rates_corr = read.table(paste0("reaMigAnalysis_mean_rates_corr_no_data_reject_",i, ".txt"), header=TRUE)
    # 
    # rea_window_mean_rates_corr <- data.frame(value = data_mean_rates_corr$n_reassortment_window/data_mean_rates_corr$network_length_window, variable = "On Window")
    # full_outside_window_mean_rates_corr <- data.frame(value = (data_mean_rates_corr$n_reassortment-data_mean_rates_corr$n_reassortment_window)/
    #                                                (data_mean_rates_corr$network_length-data_mean_rates_corr$network_length_window),
    #                                              variable = "Off Window")
    # 
    # 
    # dat_mean_rates_corr <- rbind(dat_mean_rates_corr, cbind(rbind(rea_window_mean_rates_corr, full_outside_window_mean_rates_corr), data.frame(window_size = i)))
    # diff_mean_rates_corr <- rbind(diff_mean_rates_corr, data.frame(value = rea_window_mean_rates_corr$value - full_outside_window_mean_rates_corr$value, variable = i, method = "mean rates forced correlation"))
    
    
    
    ### Posterior Rates, No Data Simulation ### 
    
    system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ",
                  i," -burnin 0 -minTipDistance ",fit," ../FullDatesInference/simulation/posterior_rates/root_filter_max/run_",run,"/h5n1_combined.rootFilter.sim.network.trees root_filter_max/run_",run,"/reaMigAnalysis_post_rates_no_data_reject_",i,"_",fit,".txt"))
    
    data_post_rates = read.table(paste0("root_filter_max/run_",run,"/reaMigAnalysis_post_rates_no_data_reject_",i,"_",fit,".txt"), header=TRUE)
    
    rea_window_post_rates <- data.frame(value = data_post_rates$n_reassortment_window/data_post_rates$network_length_window, variable = "On Window")
    full_outside_window_post_rates <- data.frame(value = (data_post_rates$n_reassortment-data_post_rates$n_reassortment_window)/
                                                   (data_post_rates$network_length-data_post_rates$network_length_window),
                                                 variable = "Off Window")
    
    
    dat_post_rates <- rbind(dat_post_rates, cbind(rbind(rea_window_post_rates, full_outside_window_post_rates), data.frame(window_size = i)))
    diff_post_rates <- rbind(diff_post_rates, data.frame(value = rea_window_post_rates$value - full_outside_window_post_rates$value, variable = i, method = "posterior rates simulation"))
  }
  
  dat_full$window_size <- as.factor(dat_full$window_size)
  dat_fit$window_size <- as.factor(dat_fit$window_size)
  dat_unfit$window_size <- as.factor(dat_unfit$window_size)
  
  # dat_prior$window_size <- as.factor(dat_prior$window_size)
  # dat_mean_rates$window_size <- as.factor(dat_mean_rates$window_size)
  # dat_mean_rates_corr$window_size <- as.factor(dat_mean_rates_corr$window_size)
  dat_post_rates$window_size <- as.factor(dat_post_rates$window_size)
  
  # diff_data <- rbind(diff_prior, diff_mean_rates, diff_mean_rates_corr, diff_post_rates, diff_full, diff_fit, diff_unfit)
  # diff_data <- rbind(diff_prior, diff_post_rates, diff_full, diff_fit, diff_unfit)
  diff_data <- rbind(diff_post_rates, diff_full, diff_fit, diff_unfit)
  diff_data$variable <- as.factor(diff_data$variable)
  
  
  
  p_full<- ggplot(dat_full, aes(x = window_size, y = value, color = variable)) +
    geom_violin(scale = "width", adjust = 1, width = 0.5) +
    labs(##title = "Full Network", sub#title = "Analysis for full network and when network lineages are split in fit and unfit ",
      x="Window size in years", y="Reassortment rate") +
    scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
  
  p_fit<- ggplot(dat_fit, aes(x = window_size, y = value, color = variable)) +
    geom_violin(scale = "width", adjust = 1, width = 0.5) +
    labs(#title = "Fit Network Edges", sub#title = "Window size in years ",
      x="Window size in years", y="Reassortment rate") +
    scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
  
  p_unfit<- ggplot(dat_unfit, aes(x = window_size, y = value, color = variable)) +
    geom_violin(scale = "width", adjust = 1, width = 0.5) +
    labs(#title = "Unfit Network Edges", sub#title = "Window size in years ",
      x="Window size in years", y="Reassortment rate") +
    scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
  
  pdf(paste0("root_filter_max/run_",run,"/h5n1_reaMigAnalysis_reject_",fit,".pdf"))
  grid.arrange(p_full, p_fit, p_unfit)
  dev.off()
  
  
  p_combined<- ggplot(rbind(dat_full, dat_fit, dat_unfit), aes(x = window_size, y = value, color = variable)) +
    geom_violin(scale = "width", draw_quantiles = 0.5) +
    labs(#title = "Combined", sub#title = "Analysis for full network and when network lineages are split in fit and unfit",
      x="Window size in years", y="Reassortment rate") +
    scale_colour_colorblind() + scale_fill_colorblind() + theme_minimal() + theme(text = element_text(size=24))
  
  ggsave(plot=p_combined, paste0("root_filter_max/run_",run,"/h5n1_reject_combined_",fit,".pdf"), width=17)
  
  
  # p_mean_rates <- ggplot(dat_mean_rates, aes(x = window_size, y = value, color = variable)) +
  #   geom_violin(scale = "width", adjust = 1, width = 0.5) +
  #   labs(#title = "Reassortment and Migration Correlation",
  #        x="Window size in years ", y="Reassortment rate") +
  #   scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
  # 
  # ggsave(plot=p_mean_rates, "h5n1_reaMigAnalysis_mean_rates_reject.pdf")
  
  p_post_rates <- ggplot(dat_post_rates, aes(x = window_size, y = value, color = variable)) +
    geom_violin(scale = "width", adjust = 1, width = 0.5) +
    labs(#title = "Reassortment and Migration Correlation",
      x="Window size in years ", y="Reassortment rate") +
    scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
  
  ggsave(plot=p_post_rates, paste0("root_filter_max/run_",run,"/h5n1_reaMigAnalysis_post_rates_reject_",fit,".pdf"))
  
  # p_prior <- ggplot(dat_prior, aes(x = window_size, y = value, color = variable)) +
  #   geom_violin(scale = "width", adjust = 1, width = 0.5) +
  #   labs(#title = "Reassortment and Migration Correlation",
  #        x="Window size in years ", y="Reassortment rate") +
  #   scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
  # 
  # ggsave(plot=p_prior, paste0("run_",run,"/h5n1_reaMigAnalysis_prior_reject.pdf"))
  
  
  p_combined_diff<- ggplot(diff_data, aes(x = variable, y = value, color=method)) +
    geom_violin(scale = "width", trim=T, draw_quantiles = 0.5) +
    labs(#title = "Reassortment and Migration Correlation", sub#title = "Difference between on and off window reassortment rates when sampling from the prior and using data to inform posterior", 
      x="Window size in years",
      y="Reassortment rate difference between \n on and off window, before migration event") +
    scale_colour_colorblind() + scale_fill_colorblind() + theme_minimal() + theme(text = element_text(size=24)) + scale_y_continuous(limits = c(-0.4, 0.4))
  
  ggsave(plot=p_combined_diff, paste0("root_filter_max/run_",run,"/h5n1_combined_reject_diff_",fit,".pdf"), width=17)
  save.image(file=paste0("root_filter_max_run_",run,".RData"))
  
}
# score_full = wes_palette(n=4, name="Darjeeling1")[1]
# score_fit = wes_palette(n=4, name="Darjeeling1")[2]
# score_unfit = wes_palette(n=4, name="Darjeeling1")[3]
# mean_rates = wes_palette(n=4, name="GrandBudapest1")[3]
# 
# pdf("mean_rates_to_infernece_diff.pdf", width = 15, height = 8) 
# par(cex=1.5, mfrow=c(2,4))
# for (i in window_size) {
#   full_diff <- (diff_full$method=="score_full" & diff_full$variable==i)
#   fit_diff <- (diff_fit$method=="score_fit" & diff_fit$variable==i)
#   unfit_diff <- (diff_unfit$method=="score_unfit" & diff_fit$variable==i)
#   mean_diff <- (diff_mean_rates$method=="mean rates" & diff_mean_rates$variable==i)
#   
#   plot(ecdf(diff_fit$value[fit_diff]),col=score_fit, ylab="Fn(x)", xlab="Difference between on and off window reassortment rates", main=paste("Window size:",i))
#   lines(ecdf(diff_mean_rates$value[mean_diff]),col=mean_rates)
#   lines(ecdf(diff_full$value[full_diff]), col=score_full)
#   lines(ecdf(diff_unfit$value[unfit_diff]),col=score_unfit)
#   
#   ks_full<-ks.boot(diff_full$value[full_diff], diff_mean_rates$value[mean_diff] , nboots=500)
#   ks_fit<-ks.boot(diff_fit$value[fit_diff], diff_mean_rates$value[mean_diff] , nboots=500)
#   ks_unfit<-ks.boot(diff_unfit$value[unfit_diff], diff_mean_rates$value[mean_diff] , nboots=500)
#   
#   legend("bottomright", inset=0.05,
#          c(paste0("score full,D=",round(ks_full$ks$statistic, 3)),paste0("score fit,D=",round(ks_fit$ks$statistic, 3)),paste0("score unfit,D=",round(ks_unfit$ks$statistic, 3)), "sim. mean rates"),
#          lty=c(1,1), col=c(score_full,score_fit,score_unfit,mean_rates))
#   }
# dev.off()
# 
# 
# pdf("mean_rates_to_infernece_window.pdf", width = 15, height = 8) 
# par(cex=1.5, mfrow=c(2,4))
# for (i in window_size) {
#   full_window <- (dat_full$variable=="Full On Window" & dat_full$window_size==i)
#   fit_window <- (dat_fit$variable=="Fit On Window" & dat_fit$window_size==i)
#   unfit_window <- (dat_unfit$variable=="Unfit On Window" & dat_unfit$window_size==i)
#   mean_window <- (dat_mean_rates$variable=="On Window" & dat_mean_rates$window_size==i)
#   
#   plot(ecdf(dat_fit$value[fit_window]),col=score_fit, ylab="Fn(x)", xlab="Reassortment rates on window", main=paste("Window size:",i))
#   lines(ecdf(dat_mean_rates$value[mean_window]),col=mean_rates)
#   lines(ecdf(dat_full$value[full_window]), col=score_full)
#   lines(ecdf(dat_unfit$value[unfit_window]),col=score_unfit)
#   
#   ks_full<-ks.boot(dat_full$value[full_window], dat_mean_rates$value[mean_window] , nboots=500)
#   ks_fit<-ks.boot(dat_fit$value[fit_window], dat_mean_rates$value[mean_window] , nboots=500)
#   ks_unfit<-ks.boot(dat_unfit$value[unfit_window], dat_mean_rates$value[mean_window] , nboots=500)
#   
#   legend("bottomright", inset=0.05,
#          c(paste0("score full,D=",round(ks_full$ks$statistic, 3)),paste0("score fit,D=",round(ks_fit$ks$statistic, 3)),paste0("score unfit,D=",round(ks_unfit$ks$statistic, 3)), "sim. mean rates"),
#          lty=c(1,1), col=c(score_full,score_fit,score_unfit,mean_rates))
# }
# dev.off()
# 
# 
# pdf("mean_rates_to_infernece_off_window.pdf", width = 15, height = 8) 
# par(cex=1.5, mfrow=c(2,4))
# for (i in window_size) {
#   full_off_window <- (dat_full$variable=="Full Off Window" & dat_full$window_size==i)
#   fit_off_window <- (dat_fit$variable=="Fit Off Window" & dat_fit$window_size==i)
#   unfit_off_window <- (dat_unfit$variable=="Unfit Off Window" & dat_unfit$window_size==i)
#   mean_off_window <- (dat_mean_rates$variable=="Off Window" & dat_mean_rates$window_size==i)
#   
#   plot(ecdf(dat_fit$value[fit_off_window]),col=score_fit, ylab="Fn(x)", xlab="Reassortment rates off window", main=paste("Window size:",i))
#   lines(ecdf(dat_mean_rates$value[mean_off_window]),col=mean_rates)
#   lines(ecdf(dat_full$value[full_off_window]), col=score_full)
#   lines(ecdf(dat_unfit$value[unfit_off_window]),col=score_unfit)
#   
#   ks_full<-ks.boot(dat_full$value[full_off_window], dat_mean_rates$value[mean_off_window] , nboots=500)
#   ks_fit<-ks.boot(dat_fit$value[fit_off_window], dat_mean_rates$value[mean_off_window] , nboots=500)
#   ks_unfit<-ks.boot(dat_unfit$value[unfit_off_window], dat_mean_rates$value[mean_off_window] , nboots=500)
#   
#   legend("bottomright", inset=0.05,
#          c(paste0("score full,D=",round(ks_full$ks$statistic, 3)),paste0("score fit,D=",round(ks_fit$ks$statistic, 3)),paste0("score unfit,D=",round(ks_unfit$ks$statistic, 3)), "sim. mean rates"),
#          lty=c(1,1), col=c(score_full,score_fit,score_unfit,mean_rates))
# }
# dev.off()
