library("seqinr")
library("FastaUtils")

write.table(t.red[use_samples,], "subsample_metadata.tsv", , col.names=NA, sep="\t");
subset.fasta(file ="data/sequences_h5n1_ha.fasta", subset = t.red[use_samples,]$strain, out = "sequences_h5n1_ha.subsample.fasta");
subset.fasta(file ="data/sequences_h5n1_na.fasta", subset = t.red[use_samples,]$strain, out = "sequences_h5n1_na.subsample.fasta");
subset.fasta<-function(file=NULL,subset=NULL,out=paste(file,".subset",sep="")){
  library(Biostrings)
  sequences<-readDNAStringSet(file)
  if (all(as.character(subset) %in% names(sequences))==FALSE) stop("There are names in 'subset' not present in the fasta file")
  pos<-match(as.character(subset),names(sequences))
  writeXStringSet(sequences[pos],filepath=out)}
