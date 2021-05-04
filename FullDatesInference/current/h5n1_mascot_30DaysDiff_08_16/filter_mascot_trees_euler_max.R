#!/usr/bin/env Rscript

library(treeio)
library(stringr)

rm(list = ls())
set.seed(354237364)

processFile = function(haPath, naPath, haPathOut, naPathPut, max_samples) {
  con1 = file(haPath, "r")
  con2 = file(naPath, "r")
  removeLines <- c()
  i<-1
  while ( TRUE ) {
    line1 = readLines(con1, n = 1)
    line2 = readLines(con2, n = 1)
    if ( length(line1) == 0 ) {
      break
    } 
    if (i>max_samples){
      removeLines <- append(removeLines, i)
      i<-i+1
      next
    }
      
    if (substr(line2, 1, 4) != "tree"){
      write(line1,file=haPathOut,append=TRUE)
      write(line2,file=naPathPut,append=TRUE)
    } else {
      sub1 <- strsplit(line1, "&")[[1]]
      sub2 <- strsplit(line2, "&")[[1]]
      sub1<- sub1[length(sub1)]
      sub2<- sub2[length(sub2)]
      
      p.ha_ans <- as.numeric(substr(sub1, 14, 18))
      p.na_ans <- as.numeric(substr(sub2, 14, 18))
      p.ha_gal <- as.numeric(substr(sub1, 32, 36))
      p.na_gal <- as.numeric(substr(sub2, 32, 36))
      
      # aa ag ga gg
      p <- c( p.ha_ans*p.na_ans, p.ha_ans*p.na_gal, p.ha_gal*p.na_ans,  p.ha_gal*p.na_gal)
      # p<- c(p.ha_ans, p.ha_gal)
      p <- p/sum(p)
      # s<- sample(c("a", "g"), size=1, prob=p)
      s <- sample(c("aa", "ag", "ga", "gg"),size=1, prob=p)
      if (s != "aa"){
      # if (s != "a"){
        removeLines <- append(removeLines, i)
      } else{
        write(line1,file=haPathOut,append=TRUE)
        write(line2,file=naPathPut,append=TRUE)
      }
      i<- i+1
    } 
  }
  close(con1)
  close(con2)
  return(removeLines)
}

# args = commandArgs(trailingOnly=TRUE)
# # run=args[1]
# # run=1
# path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_', run, "/")
# setwd(paste0(path,"output"))

for (run in 1:10){
  
  path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_', run, "/")
  setwd(paste0(path,"output"))
#   path = "/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16"
#   setwd(paste0(path, "/run_",run,"/output"))
  out <- "root_filter_max"
  
  dir.create(out)
  ha_files <- list.files(path=paste0(path, "output"), pattern="*.ha.str.trees")
  na_files <- list.files(path=paste0(path, "output"), pattern="*.na.str.trees")
  for (ha_file in ha_files){
    na_file <- str_replace(ha_file, ".ha.", ".na.")
    # resample log file to match logging rate of log file to the trees file if there is a mismatch
    # x <- paste0('"/Applications/BEAST 2.6.2/bin/logcombiner" -log ', str_replace(ha_file, ".ha.str.trees", ".log"), ' -o ',
    #             paste0(out,"/",str_replace(ha_file, ".ha.str.trees", ".subsampled.log")), ' -b 0 -resample 500000')	
    # system(x)

    
    # y <- paste0('echo "End;" >> ', ha_file)
    # yy <- paste0('echo "End;" >> ', na_file)
    # system(y)
    # system(yy)
    max_samples <- 3001
    removeLines <- processFile(ha_file, na_file, 
                               paste0(out, "/", str_replace(ha_file, ".trees", ".rootFilter.trees")), 
                                      paste0(out, "/", str_replace(na_file, ".trees", ".rootFilter.trees")), max_samples)

    
    # ha_trees_post <- readLines(file(ha_file))
    log_post <- readLines(file(paste0(str_replace(ha_file, ".ha.str.trees", ".log"))))
    # na_trees_post <- readLines(file(na_file))
    
    # offsetTree <- which(grepl("^tree STATE_0", ha_trees_post))-1
    offsetLog <- which(grepl("^0", log_post))-1
    # offsetTreeLines <- removeLines+offsetTree
    offsetLogLines <- removeLines + offsetLog
    
    # ha_trees_post <- ha_trees_post[-offsetTreeLines]
    log_post <- log_post[-offsetLogLines]
    # na_trees_post <- na_trees_post[-offsetTreeLines]

    
    # writeLines(ha_trees_post, file(paste0(out, "/", str_replace(ha_file, ".trees", ".rootFilter.trees"))))
    writeLines(log_post, file(paste0(out, "/", str_replace(ha_file, ".ha.str.trees", ".rootFilter.log"))))
    # writeLines(na_trees_post, file(paste0(out, "/", str_replace(na_file, ".trees", ".rootFilter.trees"))))
    
    # renumber otherwise tracer won't read it
    x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner" -log ', paste0(out, "/", str_replace(ha_file, ".ha.str.trees", ".rootFilter.log")), ' -o ',
                paste0(out, "/", str_replace(ha_file, ".ha.str.trees", ".reNum.rootFilter.log")), ' -b 0')
    system(x)
    
    unlink(paste0(out, "/", str_replace(ha_file, ".ha.str.trees", ".subsampled.rootFilter.log")))
  }

}

for (run in 1:10){
  path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_', run, "/")
  setwd(paste0(path,"output"))
  out <- "root_filter_max"
  ha_files <- list.files(path=paste0(path, "output"), pattern="*.ha.str.trees")
  s <- ""
  ss_ha <- ""
  ss_na <- ""
  for (file in ha_files){
    s <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".ha.str.trees", ".rootFilter.log")))
    ss_ha <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".ha.str.trees", ".ha.str.rootFilter.trees")))
    ss_na <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".na.str.trees", ".na.str.rootFilter.trees")))
  }
  x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', s,' -o ',
              paste0(out, "/", str_replace(ha_files[3], "rep2.ha.str.trees", "combined.rootFilter.log")), ' -b 10')
  system(x)

  x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', ss_ha,' -o ',
              paste0(out, "/", str_replace(ha_files[2], "rep1.ha.str.trees", "combined.rootFilter.trees")), ' -b 10')
  system(x)

}


# trees = read.beast("/Users/jugne/Documents/SCORE-paper/FullDatesInference/h5n1_mascot_30DaysDiff_08_16/run_1/output/h5n1_rep0.ha.str.trees");
