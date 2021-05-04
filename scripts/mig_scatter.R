library(ggplot2)
library(bayesplot)
library(hexbin)
library(grid)
library(gridExtra) 
library(ggpubr)
library(directlabels)
library(ggExtra)
library(MASS)
library(reshape)


# clear workspace
rm(list = ls())

this.dir <- "/Users/jugne/Documents/SCORE-paper/scripts"
directories<-c("no_filter", "segment_root_filter_max")
# directories<-c("no_filter")
out.dir<-"scatter_plots"
dir.create(paste0(this.dir,"/",out.dir))

# color_scheme_set("teal")
color_scheme_set("darkgray")
legend<-F
for(directory in directories){

for(run in seq(1,10)){
  # Set the directory to the directory of the file
  setwd(this.dir)

  dir.create(paste0(out.dir, "/", directory))
  if (directory=="no_filter"){
    combined_log<- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.log")
    title<-"A"
  } else if (directory=="root_filter_max"){
    combined_log <- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/root_filter_max/h5n1_combined.rootFilter.log")
  } else if (directory=="segment_root_filter_max"){
    combined_log <- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output/segment_root_filter_max/h5n1_combined.rootFilter.log")
    title<-"B"
    legend<-T
  }

  data = read.table(combined_log, header=TRUE)
  if (run==1){
    data_combined <- data.frame(b_migrationRate.Anseriformes_to_Galliformes=data$b_migrationRate.Anseriformes_to_Galliformes, b_migrationRate.Galliformes_to_Anseriformes=data$b_migrationRate.Galliformes_to_Anseriformes)
  }
  else{
    data_combined <- rbind(data_combined,data.frame(b_migrationRate.Anseriformes_to_Galliformes=data$b_migrationRate.Anseriformes_to_Galliformes, b_migrationRate.Galliformes_to_Anseriformes=data$b_migrationRate.Galliformes_to_Anseriformes))
  }

#   # get the kde2d information:
#   mv.kde <- kde2d(data$b_migrationRate.Anseriformes_to_Galliformes, data$b_migrationRate.Galliformes_to_Anseriformes, n = 400)
#   dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
#   dy <- diff(mv.kde$y[1:2])
#   sz <- sort(mv.kde$z)
#   c1 <- cumsum(sz) * dx * dy
# 
#   # specify desired contour levels:
#   prob <- c(0.95)
# 
#   # plot:
#   dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
#   dc <- melt(mv.kde$z)
#   dc$prob <- approx(sz,1-c1,dc$value)$y
# 
#     p <-  mcmc_hex(mcmc(data), c("b_migrationRate.Anseriformes_to_Galliformes", "b_migrationRate.Galliformes_to_Anseriformes"), binwidth = 0.02) + #ggplot(dc,aes(x=Var1,y=Var2))+
#       geom_contour(data=dc,aes(x=Var1,y=Var2, z=prob),color="red",breaks=prob)+#scale_color_fermenter(palette = "Greens",labels=c("0.95"), breaks=prob)+
#       labs(x="Backward Migration: Ans -> Gal", y="Backward Migration: Gal -> Ans", title=paste0("Subset ", run))  + theme(legend.position = "none", text = element_text(size=30), plot.title = element_text(size=50)) + xlim(c(0,4.1)) + ylim(c(0,4.1))
# 
#   ggsave(plot=p,filename=paste0(out.dir,"/",directory,"/p",run,".pdf"), width=8, height=8)
}

  # get the kde2d information:
  mv.kde <- kde2d(data_combined$b_migrationRate.Anseriformes_to_Galliformes, data_combined$b_migrationRate.Galliformes_to_Anseriformes, n = 400)
  dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
  dy <- diff(mv.kde$y[1:2])
  sz <- sort(mv.kde$z)
  c1 <- cumsum(sz) * dx * dy

  # specify desired contour levels:
  prob <- c(0.95)

  # plot:
  dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
  dc <- melt(mv.kde$z)
  dc$prob <- approx(sz,1-c1,dc$value)$y
  if (legend){
    p <-  mcmc_hex(mcmc(data_combined), c("b_migrationRate.Anseriformes_to_Galliformes", "b_migrationRate.Galliformes_to_Anseriformes"), binwidth = 0.02) + #ggplot(dc,aes(x=Var1,y=Var2))+
      geom_contour(data=dc,aes(x=X1,y=X2, z=prob),color="red",breaks=prob)+#scale_color_fermenter(palette = "Greens",labels=c("0.95"), breaks=prob)+
    theme(legend.position=c(.9, .85), legend.text = element_text(size = 30), text = element_text(size=30)) + xlim(c(0,4.1)) + ylim(c(0,4.1))+
      labs(x="Backward Migration: Ans -> Gal", y="Backward Migration: Gal -> Ans", title = "SCoRe. Segment Root Filtering") 
  } else{
    p <-  mcmc_hex(mcmc(data_combined), c("b_migrationRate.Anseriformes_to_Galliformes", "b_migrationRate.Galliformes_to_Anseriformes"), binwidth = 0.02) + #ggplot(dc,aes(x=Var1,y=Var2))+
      geom_contour(data=dc,aes(x=X1,y=X2, z=prob),color="red",breaks=prob)+#scale_color_fermenter(palette = "Greens",labels=c("0.95"), breaks=prob)+
      theme(legend.position="none", text = element_text(size=30)) + xlim(c(0,4.1)) + ylim(c(0,4.1)) + labs(x="Backward Migration: Ans -> Gal", y="Backward Migration: Gal -> Ans",  
                                                                                                           title = "SCoRe. No filtering")
  }


  # p<-ggExtra::ggMarginal(p, type = "histogram")
  # print(direct.label(p, list("angled.boxes")))
  ggsave(plot=p,filename=paste0(out.dir,"/",directory,"/p_all.pdf"), width=8, height=8)
}




  dir.create(paste0(out.dir, "/mascot"))
  directories_m<- c("no_filter", "segment_root_filter_max")
for(directory in directories_m){
  for(run in seq(1,10)){
    # Set the directory to the directory of the file
    setwd(this.dir)

    dir.create(paste0(out.dir, "/mascot/", directory))
    if (directory=="no_filter"){
      combined_log<- paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/combined/h5n1_combined.log")
      title<-"A"
      legend<-F
    } else if (directory=="segment_root_filter_max"){
      combined_log <- paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/output/segment_filter_max/h5n1_combined.rootFilter.log")
      title<-"B"
      legend<-T
    }

    data = read.table(combined_log, header=TRUE)
      if (run==1){
        data_combined <- data.frame(b_migrationRate.Anseriformes_to_Galliformes=data$b_migrationRate.Anseriformes_to_Galliformes, b_migrationRate.Galliformes_to_Anseriformes=data$b_migrationRate.Galliformes_to_Anseriformes)
      }
      else{
        data_combined <- rbind(data_combined,data.frame(b_migrationRate.Anseriformes_to_Galliformes=data$b_migrationRate.Anseriformes_to_Galliformes, b_migrationRate.Galliformes_to_Anseriformes=data$b_migrationRate.Galliformes_to_Anseriformes))
      }
    
    # # get the kde2d information: 
    # mv.kde <- kde2d(data$b_migrationRate.Anseriformes_to_Galliformes, data$b_migrationRate.Galliformes_to_Anseriformes, n = 400)
    # dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
    # dy <- diff(mv.kde$y[1:2])
    # sz <- sort(mv.kde$z)
    # c1 <- cumsum(sz) * dx * dy
    # 
    # # specify desired contour levels:
    # prob <- c(0.95)
    # 
    # # plot:
    # dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
    # dc <- melt(mv.kde$z)
    # dc$prob <- approx(sz,1-c1,dc$value)$y
    # p <-  mcmc_hex(mcmc(data), c("b_migrationRate.Anseriformes_to_Galliformes", "b_migrationRate.Galliformes_to_Anseriformes"), binwidth = 0.02) + #ggplot(dc,aes(x=Var1,y=Var2))+
    #   geom_contour(data=dc,aes(x=Var1,y=Var2, z=prob),color="red",breaks=prob)+#scale_color_fermenter(palette = "Greens",labels=c("0.95"), breaks=prob)+
    #   theme(legend.position="none", text = element_text(size=30)) + xlim(c(0,4.1)) + ylim(c(0,4.1), title=paste0("Subset ", run)) + labs(x="Backward Migration: Ans -> Gal", y="Backward Migration: Gal -> Ans")
    # 
    # ggsave(plot=p,filename=paste0(out.dir,"/mascot/",directory,"/p",run,".pdf"), width=8, height=8) 
  }
  
  # get the kde2d information: 
  mv.kde <- kde2d(data_combined$b_migrationRate.Anseriformes_to_Galliformes, data_combined$b_migrationRate.Galliformes_to_Anseriformes, n = 400)
  dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
  dy <- diff(mv.kde$y[1:2])
  sz <- sort(mv.kde$z)
  c1 <- cumsum(sz) * dx * dy
  
  # specify desired contour levels:
  prob <- c(0.95)
  
  # plot:
  dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
  dc <- melt(mv.kde$z)
  dc$prob <- approx(sz,1-c1,dc$value)$y
  if (legend){
    p <-  mcmc_hex(mcmc(data_combined), c("b_migrationRate.Anseriformes_to_Galliformes", "b_migrationRate.Galliformes_to_Anseriformes"), binwidth = 0.02) + #ggplot(dc,aes(x=Var1,y=Var2))+
      geom_contour(data=dc,aes(x=X1,y=X2, z=prob),color="red",breaks=prob)+#scale_color_fermenter(palette = "Greens",labels=c("0.95"), breaks=prob)+
      theme(legend.position=c(.9, .85), legend.text = element_text(size = 30), text = element_text(size=30)) + xlim(c(0,4.1)) + ylim(c(0,4.1))+
      labs(x="Backward Migration: Ans -> Gal", y="Backward Migration: Gal -> Ans", title = "MASCOT. Segment Root Filtering") 
  }else{
    p <-  mcmc_hex(mcmc(data_combined), c("b_migrationRate.Anseriformes_to_Galliformes", "b_migrationRate.Galliformes_to_Anseriformes"), binwidth = 0.02) + #ggplot(dc,aes(x=Var1,y=Var2))+
      geom_contour(data=dc,aes(x=X1,y=X2, z=prob),color="red",breaks=prob)+#scale_color_fermenter(palette = "Greens",labels=c("0.95"), breaks=prob)+
      theme(legend.position="none", text = element_text(size=30)) + xlim(c(0,4.1)) + ylim(c(0,4.1))+
      labs(x="Backward Migration: Ans -> Gal", y="Backward Migration: Gal -> Ans",  title = "MASCOT. No filtering") 
  }

  
  
  # p<-ggExtra::ggMarginal(p, type = "histogram")
  # print(direct.label(p, list("angled.boxes")))
  ggsave(plot=p,filename=paste0(out.dir,"/mascot/",directory,"/p_all.pdf"), width=8, height=8)

}

  ## check percentage of networks for which both segment roots have identical heighs run_1
  combined_log<- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_1/combined/h5n1_combined.log")
  data = read.table(combined_log, header=TRUE)
  ha_root <- data$ha.tree.height
  na_root <- data$na.tree.height
  network_root <- data$network.obsHeight
  
  n_match <- length(which(ha_root==na_root))
  n_total <- length(ha_root)
  
  percent_match_score<-n_match*100/n_total
  
  combined_log<- paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_1/combined/h5n1_combined.log")
  data = read.table(combined_log, header=TRUE)
  ha_root <- data$ha.tree.height
  na_root <- data$na.tree.height
  network_root <- data$network.obsHeight
  
  n_match <- length(which(ha_root==na_root))
  n_total <- length(ha_root)
  
  percent_match_masc<-n_match*100/n_total
  
  write.csv(data.frame(score=percent_match_score, mascot=percent_match_masc), file=paste0(out.dir,"/run_1_matching_seg_roots.csv"))
  
  ## check percentage of networks for which both segment roots have identical heighs, run 5
  
  combined_log<- paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_8/combined/h5n1_combined.log")
  data = read.table(combined_log, header=TRUE)
  ha_root <- data$ha.tree.height
  na_root <- data$na.tree.height
  network_root <- data$network.obsHeight
  
  n_match <- length(which(ha_root==na_root))
  n_total <- length(ha_root)
  
  percent_match_score<-n_match*100/n_total
  
  combined_log<- paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_8/combined/h5n1_combined.log")
  data = read.table(combined_log, header=TRUE)
  ha_root <- data$ha.tree.height
  na_root <- data$na.tree.height
  network_root <- data$network.obsHeight
  
  n_match <- length(which(ha_root==na_root))
  n_total <- length(ha_root)
  
  percent_match_masc<-n_match*100/n_total
  
  write.csv(data.frame(score=percent_match_score, mascot=percent_match_masc), file=paste0(out.dir,"/run_8_matching_seg_roots.csv"))
  
  