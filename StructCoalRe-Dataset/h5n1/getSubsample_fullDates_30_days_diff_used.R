library(ggplot2)
library(stringr)
library(lubridate)

# clear workspace
rm(list = ls())

set.seed(297957368)
# Set the directory to the directory of the file
# directory <- dirname(parent.frame(2)$ofile)
# setwd(directory)

directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory)

source(paste(directory, "/../datasetUtils_used.R", sep=""))

segments = c("ha", "na")
virusname = "h5n1"
used_order = c("Anseriformes", "Galliformes")
nr_samples = 200
min_length = c(800, 1000)
max_length = c(1800, 1500)
max_segment_length = c(1700,1350)
range=c(2008,2016)

for (run in 1:10){
  system(paste0("mkdir fasta_fullDates_30DaysDiff_08_16/run_",run))
# clean the dataset to only include sequences from birds and those that have all segments
if (file.exists(paste(directory, "/cleaned/metadata_",virusname,".tsv", sep=""))){
  print("metadata before alignment already exists, skip this step")
}else{
  t = list()
  f = list()
  for (i in seq(1,length(segments))){
    t[[i]] = read.table(paste(directory, "/data/metadata_", virusname, "_", segments[[i]], ".tsv", sep=""), header=T, sep="\t")
    f[[i]] = readFastaFile(paste(directory, "/data/sequences_", virusname, "_", segments[[i]], ".fasta", sep=""))
  }
  
  # define the order to be used, other absorbs everything else
  t.red = getUniqueMetaData(directory, virusname, segments, min_length, max_length, f, t)
  
  # get the "used" orders, with anything not specified as others
  t.red = getReducedOrder(t.red, used_order)
  # print t.red and sequences to new file for alignemnt
  printFastaAndMeta(directory, virusname, segments, f, t.red)
}

# read in sequences
if (file.exists(paste(directory, "/cleaned/metadata_",virusname,".cleaned.tsv", sep=""))){
  print("cleaned metadata already exists, skip this step")
}else{
  # read in metadata
  t.red = read.table(paste(directory, "/cleaned/metadata_", virusname, ".tsv", sep=""), header=T, sep="\t")
  getCleanedSeq(directory, virusname,segments, t.red)
}

# read in the cleaned data and build a quick phylogeny for all the segments
t.red = read.table(paste(directory, "/cleaned/metadata_", virusname, ".cleaned.tsv", sep=""), header=T, sep="\t")

# get the "used" orders, with anything not specified as others
t.red = getReducedOrder(t.red, used_order)

# remove rows with unknown regions
t.red = t.red[-which(t.red$region=="?"),]

# remove rowns with order of "other"
t.red = t.red[-which(t.red$used_order=="other"),]

# remove rowns with incomplete dates
t.red = t.red[-which(str_detect(t.red$date,".*XX")),]


# remove rows with dates below set range
t.red = t.red[-which(as.integer(format(as.Date(t.red$date, "%Y"), "%Y"))<range[1]),]
# remove rows with dates above set range
if (length(which(as.integer(format(as.Date(t.red$date, "%Y"), "%Y"))>range[2]))>0){
  t.red = t.red[-which(as.integer(format(as.Date(t.red$date, "%Y"), "%Y"))>range[2]),]
}

idxx = c();
keep = c();
t.red.copy <- t.red
for (ii in 1:length(t.red[,1])){
  date = as.Date((t.red$date[ii]), "%Y-%m-%d")
  region = as.character(t.red$region[ii])
  order = as.character(t.red$order[ii])
  idx = which(abs(decimal_date(as.Date(t.red$date, "%Y-%m-%d"))- decimal_date(date))<0.08333333 &
                as.character(t.red$region) == region &  
                as.character(t.red$order)==order)
  if (all(idx %in% idxx))
    next
  if (length(idx)>1){
    ss <- sample(idx,1)
    while (ss %in% keep | ss %in% idxx){
      ss <- sample(idx,1)
    }
    keep<-append(keep, ss)
    idxx <- unique(append(idxx, idx))
  }
}
t.keep <- t.red.copy[keep,]
t.red <- t.red[-idxx,]
t.red <- rbind(t.red, t.keep)


# for (ii in 1:length(t.red[,1])){
#   date = as.Date((t.red$date[ii]))
#   country = as.character(t.red$country[ii])
#   order = as.character(t.red$order[ii])
#   idx = which(abs(as.integer(as.Date(t.red$date[ii], "%Y-%m-%d")- date))<10 & 
#                 as.character(t.red$country) == country &  
#                 as.character(t.red$order)==order)
#   if (length(idx)>1){
#     keep = append(keep, sample(idx,1))
#     idx=idx[-which(idx==sample(idx,1))]
#     idxx = unique(append(idxx, idx))
#   }
#  
# }
# t.red = t.red[-idxx,]

# plot the distribution of birds
ggplot() + geom_histogram(data=t.red, aes(order), stat="count") + xlab("order")
ggplot() + geom_histogram(data=t.red, aes(used_order), stat="count") + xlab("order")



t.red.ans <- t.red[which(t.red$used_order=="Anseriformes"),]
t.red.gal <- t.red[which(t.red$used_order=="Galliformes"),]
# get a random subset that subsamples over time and space
use_samples_1 = getNrSamples(t.red.ans, nr_samples/2)
use_samples_2 = getNrSamples(t.red.gal, nr_samples/2)

use_samples = c(which(t.red$strain %in% t.red.ans[use_samples_1, "strain"]), which(t.red$strain %in% t.red.gal[use_samples_2, "strain"]))
# use_samples = which(t.red$strain %in% t.red.gal[use_samples_2, "strain"])

# use_samples = c(use_samples_1,use_samples_2)

# plot the number of samples by year
# plot.year = data.frame(year=getSamplingYear(t.red)[use_samples],region=t.red[use_samples, "region"],order=t.red[use_samples, "order"])
# p_year_ans = ggplot() + geom_histogram(data=plot.year[which(plot.year$order=="Anseriformes"),], aes(year), stat="count") + xlab("year")
# p_year_gal = ggplot() + geom_histogram(data=plot.year[which(plot.year$order=="Galliformes"),], aes(year), stat="count") + xlab("year")
# p_region_ans = ggplot() + geom_histogram(data=plot.year[which(plot.year$order=="Anseriformes"),], aes(region), stat="count") + xlab("region")
# p_region_gal = ggplot() + geom_histogram(data=plot.year[which(plot.year$order=="Galliformes"),], aes(region), stat="count") + xlab("region")
# 
# ggsave(paste0("fasta_fullDates_30DaysDiff_08_16/run_",run,"_year_ans.pdf"), p_year_ans)
# ggsave(paste0("fasta_fullDates_30DaysDiff_08_16/run_",run,"_year_gal.pdf"), p_year_gal)
# ggsave(paste0("fasta_fullDates_30DaysDiff_08_16/run_",run,"_region_ans.pdf"), p_region_ans)
# ggsave(paste0("fasta_fullDates_30DaysDiff_08_16/run_",run,"_region_gal.pdf"), p_region_gal)

# print the data to file
for (i in seq(1,length(segments))){
  fasta = readFastaFile(paste(directory, "/cleaned/sequences_", virusname, "_", segments[[i]], ".cleaned.afasta", sep=""))
  first = T
  for (j in seq(1,length(use_samples))){
    header = paste(t.red[use_samples[[j]], "strain"], t.red[use_samples[[j]], "region"], 
                   t.red[use_samples[[j]], "order"], t.red[use_samples[[j]], "date"], sep="|")
    seq.val = as.character(fasta[which(fasta$strain == t.red[use_samples[[j]], "strain"]), "sequence"])
    new.dat = data.frame(strain=header, sequence=seq.val)
    if (first){
      new.fasta = new.dat
      first = F
    }else{
      new.fasta = rbind(new.fasta, new.dat)
    }
  }
  # for (i in seq(1,length(segments))){
    con = file(paste(directory, "/fasta_fullDates_30DaysDiff_08_16/run_",run,"/sequences_",virusname, "_", segments[[i]], ".fasta", sep=""), "w")
    for (j in seq(1,length(new.fasta$strain))){
      write(paste(">",new.fasta[j, "strain"], "",sep=""), file=con, append=TRUE)
      write(paste(new.fasta[j, "sequence"], "",sep=""), file=con, append=TRUE)
    }
  # }
  close(con)
  
}
}


