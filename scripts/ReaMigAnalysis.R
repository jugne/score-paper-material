library(ggplot2)
library(scales)
library(ggthemes)
library(gridExtra)
library(wesanderson)
library(xtable)
library(stringr)
# library(wesanderson)

# clear workspace
rm(list = ls())
fits<- c(1, 2, 4,6)
# fits <- c(2)
# directory <- "no_filter"
this.dir <- "/Users/jugne/Documents/SCORE-paper/scripts"
window_size <- c(0.5, 1.0, 1.5, 2.0,3.0,4.0,5.0)
# window_size <- c(2.0)
# directories<-c("no_filter","root_filter_max", "segment_root_filter_max")
# directories<-c("root_filter_max", "segment_root_filter_max")
directories<-c("no_filter", "segment_root_filter_max")
for (fit in fits){
  for(directory in directories){

all_runs <- data.frame()
summary_total<- data.frame()
dat_fit_unfit_diff <- data.frame(value = double(), run = integer())
diff_fit_unfit_window <- data.frame(value = double(), run = integer());
for(run in seq(1,10)){

# Set the directory to the directory of the file
setwd(this.dir)

# system(paste0("rm -r run_",run))
# system(paste0("mkdir -p ",directory,"/run_",run))

if (directory=="no_filter"){
  combined_trees_inf <- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.typed.network.trees")
  combined_trees_sim <- paste0("../FullDatesInference/simulation/posterior_rates/no_filter/run_",run,"/h5n1_combined.sim.network.trees")
} else if (directory=="root_filter_max"){
  combined_trees_inf <- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.trees")
  combined_trees_sim <- paste0("../FullDatesInference/simulation/posterior_rates/root_filter_max/run_",run,"/h5n1_combined.rootFilter.sim.network.trees")
} else if (directory=="segment_root_filter_max"){
  combined_trees_inf <- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max/h5n1_combined.rootFilter.trees")
  combined_trees_sim <- paste0("../FullDatesInference/simulation/posterior_rates/segment_root_filter_max/run_",run,"/h5n1_combined.rootFilter.sim.network.trees")
}

dat_full <- data.frame(value = double(), position = character(), window_size = integer());
dat_fit <- data.frame(value = double(), position = character(), window_size = integer());

dat_unfit <- data.frame(value = double(), position = character(), window_size = integer());
dat_prior <- data.frame(value = double(), window_size = character(), window_size = integer());
# dat_mean_rates <- data.frame(value = double(), network_split = character(), window_size = integer());
# dat_mean_rates_corr <- data.frame(value = double(), network_split = character(), window_size = integer());
dat_post_rates <- data.frame(value = double(), window_size = character(), window_size = integer());

diff_prior <- data.frame(value = double(), network_split = integer(), method = character());
# diff_mean_rates <- data.frame(value = double(), network_split = integer(), method = character());
# diff_mean_rates_corr <- data.frame(value = double(), network_split = integer(), method = character());
diff_post_rates <- data.frame(value = double(), window_size = integer(), method = character());
diff_full <- data.frame(value = double(), window_size = integer(), method = character());
diff_fit <- data.frame(value = double(), window_size = integer(), method = character());
diff_unfit <- data.frame(value = double(), window_size = integer(), method = character());


dist_diff <- data.frame(value = double(), window_size = integer(), method = character());

str<-""
summary_window <- data.frame()
summary_<- data.frame(row.names = c("Total reassortment node count", "On window reassortment node count",  "Off window reassortment node count", "Network length"))

for (i in window_size) {
  analysis_file_inf <- paste0(directory,"/new/fit_",fit,"/run_",run,"/reaMigAnalysis_",format(i, nsmall = 1),"_",fit,"_new.txt")
  analysis_file_sim <- paste0(directory,"/new/fit_",fit,"/run_",run,"/reaMigAnalysis_post_rates_no_data_reject_",format(i, nsmall = 1),"_",fit,"_new.txt")
  if (!file.exists(analysis_file_inf)){
      system(paste0("mkdir -p ",directory,"/new/fit_",fit,"/run_",run))
      system(paste0("java -jar ./ReaMigAnalysis3.jar -maxReaMigDistance ", i,
                  " -burnin 0 -minTipDistance ", fit, " ", combined_trees_inf," ", analysis_file_inf))
  }
      
  data = read.table(analysis_file_inf, header=TRUE)
  # round since there are values -4*10^(-14) which are obviously zero, but should not be negative
  data <- round(data, digits=10)
  
  rea_window <- data.frame(value = data$n_reassortment_window/data$network_length_window, fitness="Full", window_split = "On")
  rea_outside_window <- data.frame(value = data$n_reassortment_off_window/
                                      (data$network_length-data$network_length_window), fitness="Full",
                                   window_split = "Off")
  
  dat_full<- rbind(dat_full,cbind(rbind(rea_window, rea_outside_window), data.frame(window_size = i)))
  diff_full<- rbind(diff_full, data.frame(value = rea_window$value - rea_outside_window$value, window_size = i, method = "SCoRe_full"))
  
  
  ### fit 
  
  rea_fit_window <- data.frame(value = data$n_reassortment_window_fit/data$network_length_window_fit, fitness="Fit",
                               window_split = "On")
  rea_fit_off_window <- data.frame(value = data$n_reassortment_off_window_fit/data$network_length_off_window_fit, fitness="Fit",
                                   window_split = "Off")
  
  dat_fit<- rbind(dat_fit,cbind(rbind(rea_fit_window, rea_fit_off_window), data.frame(window_size = i)))
  diff_fit<- rbind(diff_fit, data.frame(value = rea_fit_window$value - rea_fit_off_window$value, window_size = i, method = "SCoRe_fit"))
  
  ### unfit 
  
  rea_unfit_window <- data.frame(value = data$n_reassortment_window_unfit/data$network_length_window_unfit, fitness="Unfit",
                                 window_split = "On")
  rea_unfit_off_window <- data.frame(value = data$n_reassortment_off_window_unfit/data$network_length_off_window_unfit,fitness="Unfit",
                                     window_split = "Off")
  
  dat_unfit<- rbind(dat_unfit,cbind(rbind(rea_unfit_window, rea_unfit_off_window), data.frame(window_size = i)))
  diff_unfit<- rbind(diff_unfit, data.frame(value = rea_unfit_window$value - rea_unfit_off_window$value, window_size = i, method = "SCoRe_unfit"))
  
  if (i == 2){
    diff_fit_unfit_window <- rbind(diff_fit_unfit_window, data.frame(value=data$n_reassortment_window_fit/data$network_length_fit - data$n_reassortment_window_unfit/data$network_length_unfit, run=run))
  }
  ### Prior ###
  
  # system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ",
  #               i," -burnin 0 -minTipDistance 2 /Users/jugne/Documents/SCORE-paper/TipDateInference/h5n1_prior/reject/run_",run,"/h5n1_rep0.typed.network.trees run_",run,"/reaMigAnalysis_prior_reject_",i, ".txt"))
  # 
  # data_prior = read.table(paste0("run_",run,"/reaMigAnalysis_prior_reject_",i, ".txt"), header=TRUE)
  # 
  # rea_window_prior <- data.frame(value = data_prior$n_reassortment_window/data_prior$network_length_window, network_split = "On Window")
  # full_outside_window_prior <- data.frame(value = (data_prior$n_reassortment-data_prior$n_reassortment_window)/
  #                                           (data_prior$network_length-data_prior$network_length_window),
  #                                         network_split = "Off Window")
  # 
  # 
  # dat_prior <- rbind(dat_prior, cbind(rbind(rea_window_prior, full_outside_window_prior), data.frame(window_size = i)))
  # diff_prior <- rbind(diff_prior, data.frame(value = rea_window_prior$value - full_outside_window_prior$value, network_split = i, method = "prior"))
  
  # ### Mean Rates, No Data Simulation ###
  # 
  # # system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ",
  # #               i," -burnin 0 -minTipDistance 2 /Users/jugne/Documents/SCORE-paper/TipDateInference/simulation/h5n1_run2_combined.sim.network.trees reaMigAnalysis_mean_rates_no_data_reject_",i, ".txt"))
  # 
  # data_mean_rates = read.table(paste0("reaMigAnalysis_mean_rates_no_data_reject_",i, ".txt"), header=TRUE)
  # 
  # rea_window_mean_rates <- data.frame(value = data_mean_rates$n_reassortment_window/data_mean_rates$network_length_window, network_split = "On Window")
  # full_outside_window_mean_rates <- data.frame(value = (data_mean_rates$n_reassortment-data_mean_rates$n_reassortment_window)/
  #                                           (data_mean_rates$network_length-data_mean_rates$network_length_window),
  #                                         network_split = "Off Window")
  # 
  # 
  # dat_mean_rates <- rbind(dat_mean_rates, cbind(rbind(rea_window_mean_rates, full_outside_window_mean_rates), data.frame(window_size = i)))
  # diff_mean_rates <- rbind(diff_mean_rates, data.frame(value = rea_window_mean_rates$value - full_outside_window_mean_rates$value, network_split = i, method = "mean rates"))
  # 
  # ### Mean Rates Forced Correlation, No Data Simulation ###
  # 
  # # system(paste0("java -jar ./ReaMigAnalysis.jar -maxReaMigDistance ",
  # #               i," -burnin 0 -minTipDistance 2 /Users/jugne/Documents/SCORE-paper/TipDateInference/simulation/h5n1_run2_combined_forced_corr.sim.network.trees reaMigAnalysis_mean_rates_corr_no_data_reject_",i, ".txt"))
  # 
  # data_mean_rates_corr = read.table(paste0("reaMigAnalysis_mean_rates_corr_no_data_reject_",i, ".txt"), header=TRUE)
  # 
  # rea_window_mean_rates_corr <- data.frame(value = data_mean_rates_corr$n_reassortment_window/data_mean_rates_corr$network_length_window, network_split = "On Window")
  # full_outside_window_mean_rates_corr <- data.frame(value = (data_mean_rates_corr$n_reassortment-data_mean_rates_corr$n_reassortment_window)/
  #                                                (data_mean_rates_corr$network_length-data_mean_rates_corr$network_length_window),
  #                                              network_split = "Off Window")
  # 
  # 
  # dat_mean_rates_corr <- rbind(dat_mean_rates_corr, cbind(rbind(rea_window_mean_rates_corr, full_outside_window_mean_rates_corr), data.frame(window_size = i)))
  # diff_mean_rates_corr <- rbind(diff_mean_rates_corr, data.frame(value = rea_window_mean_rates_corr$value - full_outside_window_mean_rates_corr$value, network_split = i, method = "mean rates forced correlation"))
  
  
  
  ### Posterior Rates, No Data Simulation ### 
  
  if (!file.exists(analysis_file_sim)){
      system(paste0("java -jar ./ReaMigAnalysis3.jar -maxReaMigDistance ", i,
                  " -burnin 0 -minTipDistance ", fit, " ", combined_trees_sim," ", analysis_file_sim))
  }
  
  data_post_rates = read.table(analysis_file_sim, header=TRUE)
  data_post_rates <- round(data_post_rates, digits=10)
  
  rea_window_post_rates <- data.frame(value = data_post_rates$n_reassortment_window/data_post_rates$network_length_window,
                                      window_split = "On")
  full_outside_window_post_rates <- data.frame(value = (data_post_rates$n_reassortment-data_post_rates$n_reassortment_window)/
                                                      (data_post_rates$network_length-data_post_rates$network_length_window),
                                               window_split = "Off")
  
  
  dat_post_rates <- rbind(dat_post_rates, cbind(rbind(rea_window_post_rates, full_outside_window_post_rates), data.frame(window_size = i)))
  diff_post_rates <- rbind(diff_post_rates, data.frame(value = rea_window_post_rates$value - full_outside_window_post_rates$value, window_size = i, method = "posterior rates simulation"))
  

  str <- paste(str, paste(paste0("\\hline \\hline \\multirow{6}{*}{",i,"} & \\multirow{2}{*}{On}"),  "Reass. node count",  round(mean(data$n_reassortment_window),2), round(mean(data$n_reassortment_window_fit),2),
                          round(mean(data$n_reassortment_window_unfit),2), round(mean(data_post_rates$n_reassortment_window),2), round(mean(data_post_rates$n_reassortment_window_fit),2), round(mean(data_post_rates$n_reassortment_window_unfit),2), sep=" & "),
               paste(" &    & Network length ", round(mean(data$network_length_window),2),  round(mean(data$network_length_window_fit),2), round(mean(data$network_length_window_unfit),2),
                     round(mean(data_post_rates$network_length_window),2), round(mean(data_post_rates$network_length_window_fit),2), round(mean(data_post_rates$network_length_window_unfit),2), sep=" & "),
               paste("  &   & \\textbf{Reass. rate} ", paste0("\\textbf{",  round(mean(data$n_reassortment_window/data$network_length_window),2), "}"), paste0("\\textbf{", round(mean(data$n_reassortment_window_fit/data$network_length_window_fit),2), "}"),
                     paste0("\\textbf{", round(mean(data$n_reassortment_window_unfit/data$network_length_window_unfit),2), "}"), paste0("\\textbf{", round(mean(data_post_rates$n_reassortment_window/data_post_rates$network_length_window),2), "}"),
                     paste0("\\textbf{", round(mean(data_post_rates$n_reassortment_window_fit/data_post_rates$network_length_window_fit),2), "}"), paste0("\\textbf{", round(mean(data_post_rates$n_reassortment_window_unfit/data_post_rates$network_length_window_unfit),2), "}"), sep=" &  "),
               paste( paste0(" \\cline{2-9} & \\multirow{2}{*}{Off} & Reass. node count"),  round(mean(data$n_reassortment_off_window),2), round(mean(data$n_reassortment_off_window_fit),2),  round(mean(data$n_reassortment_off_window_unfit),2),
                      round(mean(data_post_rates$n_reassortment_off_window),2), round(mean(data_post_rates$n_reassortment_off_window_fit),2), round(mean(data_post_rates$n_reassortment_off_window_unfit),2), sep=" & "),
               paste(" &    & Network length ", round(mean(data$network_length-data$network_length_window),2),  round(mean(data$network_length_off_window_fit),2), round(mean(data$network_length_off_window_unfit),2),
                     round(mean(data_post_rates$network_length-data_post_rates$network_length_window),2), round(mean(data_post_rates$network_length_off_window_fit),2), round(mean(data_post_rates$network_length_off_window_unfit),2), sep=" & "),
               paste("  &   & \\textbf{Reass. rate} ", paste0("\\textbf{",  round(mean(data$n_reassortment_off_window/(data$network_length-data$network_length_window)),2), "}"), paste0("\\textbf{", round(mean(data$n_reassortment_off_window_fit)/mean(data$network_length_off_window_fit),2), "}"),
                     paste0("\\textbf{", round(mean(data$n_reassortment_off_window_unfit/data$network_length_off_window_unfit),2), "}"), paste0("\\textbf{", round(mean(data_post_rates$n_reassortment_off_window)/mean(data_post_rates$network_length-data_post_rates$network_length_window),2), "}"),
                     paste0("\\textbf{", round(mean(data_post_rates$n_reassortment_off_window_fit)/mean(data_post_rates$network_length_off_window_fit),2), "}"), paste0("\\textbf{", round(mean(data_post_rates$n_reassortment_off_window_unfit/data_post_rates$network_length_off_window_unfit),2), "}"), sep=" &  "), sep=" \\\\ ")
  
  
  
  summary_table_score <- data.frame(full_score=c(mean(data$n_reassortment_window), mean(data$network_length_window), mean(data$n_reassortment_window)/mean(data$network_length_window),
                                           mean(data$n_reassortment_off_window), mean(data$network_length-data$network_length_window), mean(data$n_reassortment_off_window)/mean(data$network_length-data$network_length_window)),
                                    fit_score=c(mean(data$n_reassortment_window_fit), mean(data$network_length_window_fit), mean(data$n_reassortment_window_fit)/mean(data$network_length_window_fit),
                                          mean(data$n_reassortment_off_window_fit), mean(data$network_length_off_window_fit),mean(data$n_reassortment_off_window_fit)/mean(data$network_length_off_window_fit)),
                                    unfit_score=c(mean(data$n_reassortment_window_unfit), mean(data$network_length_window_unfit), mean(data$n_reassortment_window_unfit)/mean(data$network_length_window_unfit),
                                            mean(data$n_reassortment_off_window_unfit), mean(data$network_length_off_window_unfit),mean(data$n_reassortment_off_window_unfit)/mean(data$network_length_off_window_unfit)))
                                    # row.names = c("Reassortment node count", "Network length", "Reassortment rate",
                                    #               "Reassortment node count_", "Network length_", "Reassortment rate_"))

  
  summary_table_score <- round(summary_table_score, 3)
  
  summary_table_sim <- data.frame(full_sim=c(mean(data_post_rates$n_reassortment_window), mean(data_post_rates$network_length_window), mean(data_post_rates$n_reassortment_window)/mean(data_post_rates$network_length_window),
                                               mean(data_post_rates$n_reassortment_off_window), mean(data_post_rates$network_length-data_post_rates$network_length_window), mean(data_post_rates$n_reassortment_off_window)/mean(data_post_rates$network_length-data_post_rates$network_length_window)),
                                  fit_sim=c(mean(data_post_rates$n_reassortment_window_fit), mean(data_post_rates$network_length_window_fit), mean(data_post_rates$n_reassortment_window_fit)/mean(data_post_rates$network_length_window_fit),
                                              mean(data_post_rates$n_reassortment_off_window_fit), mean(data_post_rates$network_length_off_window_fit),mean(data_post_rates$n_reassortment_off_window_fit)/mean(data_post_rates$network_length_off_window_fit)),
                                  unfit_sim=c(mean(data_post_rates$n_reassortment_window_unfit), mean(data_post_rates$network_length_window_unfit), mean(data_post_rates$n_reassortment_window_unfit)/mean(data_post_rates$network_length_window_unfit),
                                                mean(data_post_rates$n_reassortment_off_window_unfit), mean(data_post_rates$network_length_off_window_unfit),mean(data_post_rates$n_reassortment_off_window_unfit)/mean(data_post_rates$network_length_off_window_unfit)))
                                  # row.names = c("Reassortment node count", "Network length", "Reassortment rate",
                                  #               "Reassortment node count_", "Network length_", "Reassortment rate_"))
  
  summary_table_sim <- round(summary_table_sim, 3)
  
  summary_window <- rbind(summary_window ,cbind(Window=c(rep("On",3), rep("Off", 3)),summary_table_score, summary_table_sim))
}

### total fit-unfit diff

dat_fit_unfit_diff <- rbind(dat_fit_unfit_diff, data.frame(value =   (data$n_reassortment_total_fit/data$network_length_fit) - (data$n_reassortment_total_unfit/data$network_length_unfit), run=run))


dat_full$window_size <- as.factor(dat_full$window_size)
dat_fit$window_size <- as.factor(dat_fit$window_size)
dat_unfit$window_size <- as.factor(dat_unfit$window_size)

# dat_prior$window_size <- as.factor(dat_prior$window_size)
# dat_mean_rates$window_size <- as.factor(dat_mean_rates$window_size)
# dat_mean_rates_corr$window_size <- as.factor(dat_mean_rates_corr$window_size)
dat_post_rates$window_size <- as.factor(dat_post_rates$window_size)
dat_post_rates$window_split <- as.factor(dat_post_rates$window_split)

# diff_data <- rbind(diff_prior, diff_mean_rates, diff_mean_rates_corr, diff_post_rates, diff_full, diff_fit, diff_unfit)
# diff_data <- rbind(diff_prior, diff_post_rates, diff_full, diff_fit, diff_unfit)
diff_data <- rbind(diff_post_rates, diff_full, diff_fit, diff_unfit)
diff_data$window_size <- as.factor(diff_data$window_size)


# max_lim <- max(boxplot.stats(dat_full$value)$stats,boxplot.stats(dat_fit$value)$stats, boxplot.stats(dat_unfit$value)$stats)*1.5
# 
# p_full<- ggplot(dat_full, aes(x = window_size, y = value, color = window_split)) +
#   geom_boxplot(outlier.shape = NA) +
#   coord_cartesian(ylim = c(0,max_lim))+
#   labs(title = "Full Network",
#        x="Window size in years",  y="Reassortment rate") +
#   scale_colour_colorblind(name = "Window split:") + theme_bw() +
#   ggtitle("Full Network") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#     text = element_text(size=10), legend.background = element_rect(fill="white", linetype="blank"),
#     legend.key.size = unit(1,"line"), legend.position=c(.16,.9), legend.direction = "horizontal")
# 
# p_fit<- ggplot(dat_fit, aes(x = window_size, y = value, color = window_split)) +
#   geom_boxplot(outlier.shape = NA) +
#   coord_cartesian(ylim = c(0,max_lim))+
#   labs(title = "Fit Network Edges",
#        x="Window size in years", y="Reassortment rate") +
#   scale_colour_colorblind(name = "Window split:") + theme_bw() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         text = element_text(size=10), legend.background = element_rect(fill="white", linetype="blank"),
#         legend.key.size = unit(1,"line"), legend.position=c(.16,.9), legend.direction = "horizontal")
# 
# p_unfit<- ggplot(dat_unfit, aes(x = window_size, y = value, color = window_split)) +
#   geom_boxplot(outlier.shape = NA) +
#   coord_cartesian(ylim = c(0,max_lim))+
#   labs(title = "Unfit Network Edges", x="Window size in years", y="Reassortment rate") +
#   scale_colour_colorblind(name = "Window split:") + theme_bw() +
#   theme(text = element_text(size=10), legend.background = element_rect(fill="white", linetype="blank"),
#           legend.key.size = unit(1,"line"), legend.position=c(.16,.88), legend.direction = "horizontal")
# 
# pdf(paste0(directory,"/fit_",fit,"/run_",run,"/h5n1_reaMigAnalysis_reject.pdf"))
# grid.arrange(p_full, p_fit, p_unfit)
# dev.off()
# 
# okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# dat_combined <- rbind(dat_full, dat_fit, dat_unfit)
# dat_combined$window_size <- factor(dat_combined$window_size)
# dat_combined$window_split <- factor(dat_combined$window_split)
# dat_combined$fitness <- factor(dat_combined$fitness)
# p_combined<- ggplot(dat_combined, aes(x = window_size, y = value, color = window_split, fill=fitness)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.4, lwd=1, fatten=1) +
#   coord_cartesian(ylim = c(0,max_lim))+
#   labs(#title = "Combined", sub#title = "Analysis for full network and when network lineages are split in fit and unfit",
#        x="Window size in years", y="Reassortment rate") + scale_fill_manual(values=c("#0072B2", "#D55E00", "#CC79A7"), name = "Fitness split") +
#   scale_colour_manual(values=c("#000000", "#009E73"), name = "Window split")+
#   theme_minimal() + theme(text = element_text(size=24), legend.background = element_rect(fill="white", linetype="blank"),
#                           legend.key.size = unit(3,"line"), legend.position=c(.3,.85), legend.direction = "horizontal")
# 
# ggsave(plot=p_combined, paste0(directory,"/fit_",fit,"/run_",run,"/h5n1_reject_combined.pdf"), width=17)
# 
# 
# # p_mean_rates <- ggplot(dat_mean_rates, aes(x = window_size, y = value, color = network_split)) +
# #   geom_violin(scale = "width", adjust = 1, width = 0.5) +
# #   labs(#title = "Reassortment and Migration Correlation",
# #        x="Window size in years ", y="Reassortment rate") +
# #   scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
# #
# # ggsave(plot=p_mean_rates, "h5n1_reaMigAnalysis_mean_rates_reject.pdf")
# 
# max_lim <- max(boxplot.stats(dat_post_rates$value)$stats)*1.3
# min_lim <- min(boxplot.stats(dat_post_rates$value)$stats)/2
# p_post_rates <- ggplot(dat_post_rates, aes(x = window_size, y = value, color = window_split)) +
#   geom_boxplot(outlier.shape = NA, fatten=1) +
#   # stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red", fill="red") +
#   coord_cartesian(ylim = c(min_lim,max_lim))+
#   labs(#title = "Reassortment and Migration Correlation",
#        x="Window size in years ", y="Reassortment rate") +
#   scale_color_colorblind(name = "Window split:") + theme_bw() +
#   theme(legend.key.size = unit(3,"line"), text = element_text(size=24), legend.background = element_rect(fill="white", linetype="blank"),
#         legend.position=c(.28,.93), legend.direction = "horizontal")
# 
# ggsave(plot=p_post_rates, paste0(directory,"/fit_",fit,"/run_",run,"/h5n1_reaMigAnalysis_post_rates_reject.pdf"))
# 
# # p_prior <- ggplot(dat_prior, aes(x = window_size, y = value, color = network_split)) +
# #   geom_violin(scale = "width", adjust = 1, width = 0.5) +
# #   labs(#title = "Reassortment and Migration Correlation",
# #        x="Window size in years ", y="Reassortment rate") +
# #   scale_colour_colorblind() + scale_fill_colorblind() + theme_bw()
# #
# # ggsave(plot=p_prior, paste0("run_",run,"/h5n1_reaMigAnalysis_prior_reject.pdf"))
# 
# 
# max_lim <- max(boxplot.stats(diff_data$value)$stats)*2
# p_combined_diff<- ggplot(diff_data, aes(x = window_size, y = value, color=method)) +
#     geom_boxplot(outlier.shape = NA,  lwd=1, fatten=1) +
#   coord_cartesian(ylim = c(-max_lim, max_lim))+
#     labs(#title = "Reassortment and Migration Correlation", sub#title = "Difference between on and off window reassortment rates when sampling from the prior and using data to inform posterior",
#          x="Window size in years",
#          y="Reassortment rate difference between \n on and off window, before migration event") +
#   scale_colour_colorblind(labels = c("Posterior rates simulation", "SCoRe. Full network", "SCoRe. Fit network edges", "SCoRe. Unfit network edges")) +
#   theme_minimal() + theme(legend.title = element_blank(), text = element_text(size=24), legend.background = element_rect(fill="white", linetype="blank"),
#                           legend.key.size = unit(3,"line"), legend.position=c(.28,.87), legend.direction = "horizontal") + scale_y_continuous(limits = c(-0.4, 0.4)) +guides(color=guide_legend(ncol=2))
# 
# ggsave(plot=p_combined_diff, paste0(directory,"/fit_",fit,"/run_",run,"/h5n1_combined_reject_diff.pdf"), width=17)

# summary_total_<- data.frame(full_score=c(mean(data$n_reassortment), mean(data$network_length), mean(data$n_reassortment)/mean(data$network_length)),
#                             fit_score=c(mean(data$n_reassortment_total_fit), mean(data$network_length_fit), mean(data$n_reassortment_total_fit)/mean(data$network_length_fit)),
#                             unfit_score=c(mean(data$n_reassortment_total_unfit), mean(data$network_length_unfit), mean(data$n_reassortment_total_unfit)/mean(data$network_length_unfit)),
#                             full_sim=c(mean(data_post_rates$n_reassortment), mean(data_post_rates$network_length), mean(data_post_rates$n_reassortment)/mean(data_post_rates$network_length)),
#                             fit_sim=c(mean(data_post_rates$n_reassortment_total_fit), mean(data_post_rates$network_length_fit), mean(data_post_rates$n_reassortment_total_fit)/mean(data_post_rates$network_length_fit)),
#                             unfit_sim=c(mean(data_post_rates$n_reassortment_total_unfit), mean(data_post_rates$network_length_unfit), mean(data_post_rates$n_reassortment_total_unfit)/mean(data_post_rates$network_length_unfit)),
#                             row.names = c("Total reassortment node count", "Network length", "\\frac{total reassortment node count}{network length}"))
# summary_total_<- round(summary_total_, 3)
# summary_total <- rbind(summary_total, summary_total_)
# 
# str<-paste0("\\begin{table}[ht]
# \\centering
# 	\\begin{tabular}{ccr|rrr|rrr}
# 	              Window &  Network    on/off   &  &  \\multicolumn{3}{c}{Score}  &  \\multicolumn{3}{c}{Simulation}  \\\\
# 	                size &  window  &  &  Full &   Fit &  Unfit &  Full & Fit & Unfit", str, " 
# \\end{tabular}
# \\caption[Subset ",run,". Fitness threshold ",fit," years.]{\\textbf{Subset ",run,". Fitness threshold ",fit," years.} Network length, reassortment node count and reassortment rate for full network, fit/unfit network edges, additionally partitioned to be either on or off window before the migration event. }
# \\end{table}")
# 
# write(str, file=paste0("/Users/jugne/Documents/score-text/tables/",directory,"_run_",run,"_fit_",fit,"summary_window.tex"))
# write(str, file=paste0(directory,"/fit_",fit,"/run_",run,"/run_",run,"_fit_",fit,"summary_window.tex"))
# # write(print(xtable(summary_window)), file=paste0(directory,"/fit_",fit,"/run_",run,"/summary_window.tex"))
# write(print(xtable(summary_total)), file=paste0(directory,"/fit_",fit,"/run_",run,"/run_",run,"_fit",fit,"_summary_total.tex"))
# write(print(xtable(summary_total)), file=paste0("/Users/jugne/Documents/score-text/tables/run_",run,"_fit_",fit,"summary_total.tex"))


window <- cbind(diff_data, data.frame(run=rep(run, length(diff_data$value))))
all_runs <- rbind(all_runs, window)
}
for (w in window_size){
  runs <- all_runs[which(all_runs$window_size==w),]
  max_lim <- max(runs$value)
  min_lim <- min(runs$value)
  runs$run <- factor(runs$run)
  
  p_combined_diff<- ggplot(runs, aes(x = run, y = value, color=method)) +
    # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(outlier.shape = NA,  lwd=1, fatten=1) +
    coord_cartesian(ylim = c(-0.21, 0.21))+
    labs(title = paste0("Window size: ",w),
      x="Subset",
      y="Reassortment rate difference between \n on and off window, before migration event") +
    scale_colour_colorblind(labels = c("Posterior rates simulation", "SCoRe. Full network", "SCoRe. Fit network edges", "SCoRe. Unfit network edges")) +
    theme_minimal() + theme(legend.title = element_blank(), text = element_text(size=24), legend.background = element_rect(fill="white", linetype="blank"),
                            legend.key.size = unit(3,"line"), legend.position="none", legend.direction = "horizontal") +guides(color=guide_legend(ncol=2))

  ggsave(plot=p_combined_diff, paste0(directory,"/new/fit_",fit,"/h5n1_combined_fit_",fit,"_window_",w,"_diff.pdf"), width=17)
}


# dat_fit_unfit_diff$run<-as.factor(dat_fit_unfit_diff$run)
# 
# p<- ggplot(dat_fit_unfit_diff, aes(x=run, y=value))+
#   geom_boxplot(outlier.shape = NA,  lwd=1, fatten=1, notch = T)
# 
# 
# diff_fit_unfit_window$run<-as.factor(diff_fit_unfit_window$run)
# pp<- ggplot(diff_fit_unfit_window, aes(x=run, y=value))+
#   geom_boxplot(outlier.shape = NA,  lwd=1, fatten=1, notch = T)

# score_full = wes_palette(n=4, name="Darjeeling1")[1]
# score_fit = wes_palette(n=4, name="Darjeeling1")[2]
# score_unfit = wes_palette(n=4, name="Darjeeling1")[3]
# mean_rates = wes_palette(n=4, name="GrandBudapest1")[3]
# 
# pdf("mean_rates_to_infernece_diff.pdf", width = 15, height = 8) 
# par(cex=1.5, mfrow=c(2,4))
# for (i in window_size) {
#   full_diff <- (diff_full$method=="score_full" & diff_full$network_split==i)
#   fit_diff <- (diff_fit$method=="score_fit" & diff_fit$network_split==i)
#   unfit_diff <- (diff_unfit$method=="score_unfit" & diff_fit$network_split==i)
#   mean_diff <- (diff_mean_rates$method=="mean rates" & diff_mean_rates$network_split==i)
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
#   full_window <- (dat_full$network_split=="Full On Window" & dat_full$window_size==i)
#   fit_window <- (dat_fit$network_split=="Fit On Window" & dat_fit$window_size==i)
#   unfit_window <- (dat_unfit$network_split=="Unfit On Window" & dat_unfit$window_size==i)
#   mean_window <- (dat_mean_rates$network_split=="On Window" & dat_mean_rates$window_size==i)
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
#   full_off_window <- (dat_full$network_split=="Full Off Window" & dat_full$window_size==i)
#   fit_off_window <- (dat_fit$network_split=="Fit Off Window" & dat_fit$window_size==i)
#   unfit_off_window <- (dat_unfit$network_split=="Unfit Off Window" & dat_unfit$window_size==i)
#   mean_off_window <- (dat_mean_rates$network_split=="Off Window" & dat_mean_rates$window_size==i)
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
    # addtorow <- list()
    # addtorow$pos <- list(0,0,0)
    # addtorow$command <- c("& \\multicolumn{21}{c}{Network split} \\\\\n",
    #                       "& \\multicolumn{6}{c}{Window size (in years): 0.5} & \\multicolumn{6}{c}{Window size (in years): 1.0} 
    #                     & \\multicolumn{6}{c}{Window size (in years): 1.5}  & \\multicolumn{6}{c}{Window size (in years): 2.0}  
    #                     & \\multicolumn{6}{c}{Window size (in years): 3.0}  & \\multicolumn{6}{c}{Window size (in years): 4.0}
    #                     & \\multicolumn{6}{c}{Window size (in years): 5.0} \\\\\n",
    #                       "& &\\multicolumn{3}{c|}{SCoRe} & \\multicolumn{3}{c}{Simulation} &\\multicolumn{3}{||c|}{SCoRe} & \\multicolumn{3}{c}{Simulation}
    #                     &\\multicolumn{3}{||c|}{SCoRe} & \\multicolumn{3}{c}{Simulation} &\\multicolumn{3}{||c|}{SCoRe} & \\multicolumn{3}{c}{Simulation}
    #                     &\\multicolumn{3}{||c|}{SCoRe} & \\multicolumn{3}{c}{Simulation} &\\multicolumn{3}{||c|}{SCoRe} & \\multicolumn{3}{c}{Simulation} 
    #                     &\\multicolumn{3}{||c|}{SCoRe} & \\multicolumn{3}{c}{Simulation} &\\multicolumn{3}{||c|}{SCoRe} & \\multicolumn{3}{c}{Simulation}\\\\\n,
    #                        & Full & Fit & Unfit & Full & Fit & Unfit \\\\\n")

  }
}

