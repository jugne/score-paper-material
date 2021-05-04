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
fit<- 1 # it is the same for all fitness values

this.dir <- "/Users/jugne/Documents/SCORE-paper/scripts"
window_size <- 0.5 # it is the same for ll window sizes
directories<-c("no_filter","root_filter_max", "segment_root_filter_max")


  state_lengths <- NA
  for(directory in directories){
    
    state_lengths1 <- data.frame()
    
    for(run in seq(1,10)){

      setwd(this.dir)
      
      analysis_file_inf <- paste0(directory,"/fit_",fit,"/run_",run,"/reaMigAnalysis_",window_size,"_",fit,".txt")
      analysis_file_sim <- paste0(directory,"/fit_",fit,"/run_",run,"/reaMigAnalysis_post_rates_no_data_reject_",window_size,"_",fit,".txt")

        
      data = read.table(analysis_file_inf, header=TRUE)
      # round since there are values -4*10^(-14) which are obviously zero, but should not be negative
      data <- round(data, digits=10)
      
      
      table_ <- data.frame(run=run,n_length_gal=mean(data$network_length_Galliformes), n_length_ans=mean(data$network_length_Anseriformes))  
      colnames(table_) <- c("run", paste0(directory, "_n_length_gal"), paste0(directory, "_n_length_ans"))
      state_lengths1<- rbind(state_lengths1, table_)

    }
    if (is.na(state_lengths))
      state_lengths <- state_lengths1
    else
      state_lengths <- cbind(state_lengths, state_lengths1)
  }
write(print(xtable(state_lengths)), file="state_lengths.tex")
