library(stringr)
rm(list = ls())
set.seed(354237364)
out <- "segment_root_filter_max"

processFile = function(networkPath, networkPathOut, max_samples) {
  con = file(networkPath, "r")
  log <- read.table(paste0(out, "/",str_replace(networkPath, ".typed.network.trees", ".networksToRemove.txt")), header=T)
  removeLines <- c()
  i<-1
  j<-1
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    } 
    if (i>max_samples){
      removeLines <- append(removeLines, i)
      i<-i+1
      next
    }
    
    if (substr(line, 1, 4) != "tree"){
      write(line,file=networkPathOut,append=TRUE)
    } else {
      if (length(log$n_sample)>=j && log$n_sample[j]==i){
        removeLines <- append(removeLines, i)
        j<- j+1
      } else{
        write(line,file=networkPathOut,append=TRUE)
      }
      i<- i+1
    } 
  }
  close(con)
  rm(log)
  return(removeLines)
}

# for (run in 1:1){
#   path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_', run, "/")
#   setwd(paste0(path,"output"))
#   out <- paste0(path,"output","/segment_root_filter_max")
#   dir.create(out)
#   files <- list.files(path=paste0(path, "output"), pattern="*.typed.network.trees")
#   for (file in files){
#     
#     x <-paste0("java -jar /Users/jugne/Documents/Source/segmentRootFilter.jar ", file, " ",
#                paste0(out, "/", str_replace(file, ".typed.network.trees", ".networksToRemove.txt")), " -burnin 0")
#     system(x)
#     
#     max_samples <- 3001
#     removeLines <- processFile(file,
#                                paste0(out, "/", str_replace(file, ".trees", ".rootFilter.trees")), max_samples)
#     log_post <- readLines(file(paste0(str_replace(file, ".typed.network.trees", ".log"))))
#     offsetLog <- which(grepl("^0", log_post))-1
#     offsetLogLines <- removeLines + offsetLog
#     log_post <- log_post[-offsetLogLines]
#     writeLines(log_post, file(paste0(out, "/", str_replace(file, ".typed.network.trees", ".rootFilter.log"))))
#     x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner" -log ', paste0(out, "/", str_replace(file, ".typed.network.trees", ".rootFilter.log")), ' -o ',
#                 paste0(out, "/", str_replace(file, ".typed.network.trees", ".reNum.rootFilter.log")), ' -b 0')
#     system(x)
#     
#     # unlink(paste0(out, "/", str_replace(file, ".ha.str.trees", ".subsampled.rootFilter.log")))
#     
#   }
# }

for (run in 1:1){

  path=paste0('/Users/jugne/Documents/SCORE-paper/FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_', run, "/")
  # setwd(paste0(path,"output"))
  out <- paste0(path,"output","/segment_root_filter_max")
  files <- list.files(path=paste0(path, "output"), pattern="*.typed.network.trees")
  s <- ""
  ss_ha <- ""
  ss_na <- ""
  for (file in files){
    s <- paste0(s, " -log ", paste0(out, "/", str_replace(file, ".typed.network.trees", ".rootFilter.log")))
    ss_ha <- paste0(ss_ha, " -log ", paste0(out, "/", str_replace(file, ".typed.network.trees", ".typed.network.rootFilter.trees")))
  }
  x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', s,' -o ',
              paste0(out, "/", str_replace(files[3], "rep2.typed.network.trees", "combined.rootFilter.log")), ' -b 10')
  # system(x)

  x <- paste0('"/Applications/BEAST 2.6.3/bin/logcombiner"', ss_ha,' -o ',
              paste0(out, "/", str_replace(files[1], "rep0.typed.network.trees", "combined.rootFilter.trees")), ' -b 10')
  # system(x)
}

