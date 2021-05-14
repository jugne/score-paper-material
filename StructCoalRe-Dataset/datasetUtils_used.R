
readFastaFile = function (file){
  i=1
  con = file(file, "r")
  first=T
  firstseq=T
  c = 0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      new.fasta = data.frame(strain=name, sequence=seq)
      if (firstseq){
        firstseq=F
        fasta=new.fasta
      }else{
        fasta = rbind(fasta, new.fasta)
      }
      break
    }
    
    if (startsWith(line, ">")){
      if (first){
        first = F
      }else{
        new.fasta = data.frame(strain=name, sequence=seq)
        if (firstseq){
          firstseq=F
          fasta=new.fasta
        }else{
          fasta = rbind(fasta, new.fasta)
        }
      }
      name = gsub(">","",line)
      seq = ""
    }else{
      seq = paste(seq,line, sep="")
    }
  }
  
  close(con)
  return (fasta)
}


getUniqueMetaData = function (directory, virusname, segments, min_length, max_length, f, t){
  # for each element in t[[1]], see if they are in all other metadata
  rem.vals = c()
  for (i in seq(1,length(t[[1]]$strain))){
    for (j in seq(2,length(segments))){
      if (!is.element(t[[1]][[i, "strain"]], t[[j]]$strain)){
        rem.vals = append(rem.vals, i)
      }
    }
  }
  
  t.return = na.omit(t[[1]][-unique(rem.vals),])
  rem.vals = c()
  for (i in seq(1,length(t.return$strain))){
    if (sum(tolower(t.return$strain)==tolower(t.return[i,"strain"]))!=1){
      rem.vals = append(rem.vals, i)
    }
  }
  
  t.return = na.omit(t.return[-unique(rem.vals),])
  
  # rem.vals everything without an order
  bird.order = read.table(paste(directory, "/../birdOrder.csv", sep=""), header=T, sep=",", quote="")
  
  # get the order of every bird
  t.return$order = as.character(NA)
  for (i in seq(1,length(t.return$strain))){
    if (grepl("/",as.character(t.return[i, "strain"]))){
      name = strsplit(as.character(t.return[i, "strain"]), split="/")[[1]][[2]]
      t.return[i,"order"] = getOrder(name, bird.order)
    }else{
      t.return[i,"order"] = NA
    }
  }
  
  t.return = na.omit(t.return)
  
  # rem.vals everything for which the sequences are not within the length limits
  rem.vals = c()
  for (i in seq(1,length(t.return$strain))){
    for (j in seq(1,length(segments))){
      seq.length = nchar(as.character((f[[j]][which(f[[j]]$strain==as.character(t.return[i,"strain"])), "sequence"])))
      print(seq.length)
      if (seq.length<min_length[[j]] || seq.length>max_length[[j]]){
        rem.vals = append(rem.vals, i)
      }
    }
  }
  
  t.return = na.omit(t.return[-unique(rem.vals),])
  
  return (t.return)
}

printFastaAndMeta = function (directory, virusname, segments, f, t.red) {
  for (i in seq(1,length(segments))){
    con = file(paste(directory, "/cleaned/sequences_",virusname, "_", segments[[i]], ".fasta", sep=""), "w")
    for (j in seq(1,length(t.red$strain))){
      ind = which(f[[i]]$strain==as.character(t.red[j,"strain"]))
      write(paste(">",t.red[j, "strain"], "",sep=""), file=con, append=TRUE)
      write(paste(f[[i]][ind, "sequence"], "",sep=""), file=con, append=TRUE)
    }
  }
  close(con)
  
  # print meta data
  metafile = paste(directory, "/cleaned/metadata_",virusname,".tsv", sep="")
  write.table(t.red, file=metafile, quote=FALSE, sep='\t', col.names = NA)
  
  
}

getOrder = function (name, bird.order){
  index = which(bird.order$name==tolower(name))
  if (length(index)==1){
    order = as.character(bird.order[index, "order"])
  }else{
    order = NA
  }
  return (order)
}

getReducedOrder = function(t, used_order){

  t$used_order = NA
  for (i in seq(1,length(t$strain))){
    ind = which(used_order==t[i, "order"])
    if (length(ind)>0){
      t[i,"used_order"] = used_order[[ind]]
    }
  }
  t[which(is.na(t$used_order)), "used_order"] = "other"
  return(t)
}

getSamplingYear = function(t.red){
  # get the weights for each sample based on region and year
  year = rep(0,length(t.red$date))

  # get all the years
  for (i in seq(1,length(year))){
    year[[i]] = as.numeric(strsplit(as.character(t.red[i, "date"]), split="-")[[1]][[1]])
  }
  return (year)
}
  
getNrSamples = function(t.red, nr_samples){
  # get the weights for each sample based on region and year
  year = rep(0,length(t.red$date))
  region = rep("",length(t.red$date))
  
  # get all the years
  for (i in seq(1,length(year))){
    year[[i]] = as.numeric(strsplit(as.character(t.red[i, "date"]), split="-")[[1]][[1]])
    region[[i]] = as.character(t.red[i, "region"])
  }
  # compute the sampling weights of a year as the inverse number of samples
  year_weight = rep(0,length(t.red$date))
  for (i in seq(1,length(year))){
    year_weight[[i]] = 1/length(which(year==year[[i]]))
  }
  
  year_weight[is.infinite(year_weight)]=0
  # get from which years to take samples
  sampling_years = year[sample(seq(1,length(year)), nr_samples, replace=F, prob=year_weight)]
  unique.sampling_year = unique(sampling_years)
  use_samples = c()
  for (i in seq(1,length(unique.sampling_year))){
    # get the number of samples from this year
    year.samples = length(which(sampling_years==unique.sampling_year[[i]]))
    # get all samples in this year
    in.year = which(year==unique.sampling_year[[i]])
    # compute the regional weights
    regional_weight = rep(0,length(in.year))
    for (j in seq(1,length(in.year))){
      regional_weight[j] = 1/length(which(region[in.year]==region[in.year[[j]]]))
    }
    if (year.samples>1){
      use_samples = append(use_samples, sample(in.year, year.samples, replace=F, prob=regional_weight))
    }else{
      use_samples = append(use_samples, in.year)
    }
  }
  return (use_samples)
}

getCleanedSeq = function(directory, virusname, segments, t.red){
  f = list()
  remove = c()
  remove_positions = c()
  for (i in seq(1,length(segments))){
    seg.fasta = readFastaFile(paste(directory, "/cleaned/sequences_", virusname, "_", segments[[i]], ".afasta", sep=""))
    # clean the segment fasta file by remove "bad" sequences and finding the initial start codon
    codonpositions = unlist(gregexpr( pattern = "ATG" , seg.fasta$sequence))
    uni.positions = sort(unique(codonpositions))
    star.position = 0
    # get the first start codon with more than 60% frequency (arbitrary)
    for (j in seq(1,length(uni.positions))){
      freq = length(which(codonpositions==uni.positions[[j]]))/length(seg.fasta$strain)
      if (freq>0.6){
        star.position = uni.positions[[j]]
        break
      }
    }
    seg.fasta$cleanseq=NA
    # remove all chars before the start.position and after the  max_segment length
    for (j in seq(1,length(seg.fasta$strain))){
      string.val = as.character(seg.fasta[j,"sequence"])
      substring.val = substr(string.val,star.position, nchar(string.val))
      #change leading and trailing gaps to N's
      k = 1
      while(substring(substring.val,k,k)=="-"){
        substring(substring.val,k,k) ="N"
        k=k+1
      }
      k=nchar(substring.val)
      while(substring(substring.val,k,k)=="-"){
        substring(substring.val,k,k) ="N"
        k=k-1
      }
      
      seg.fasta[j,"cleanseq"] = substring.val
    }
    
    # remove all sequence that cause a gap in more than 95% of other sequences
    gappositions.list = gregexpr( pattern = "-" , seg.fasta$cleanseq)
    
    gappositions = unlist(gappositions.list)
    uni.gaps = sort(unique(gappositions))
    pos = c()
    # get the first start codon with more than 60% frequency (arbitrary)
    freq=rep(0,0,0)
    for (j in seq(1,length(uni.gaps))){
      freq[[j]] = length(which(gappositions==uni.gaps[[j]]))/(length(seg.fasta$strain)-sum(substring(seg.fasta$cleanseq,uni.gaps[[j]],uni.gaps[[j]])=="N"))
      if (freq[[j]]>0.95){
        pos = append(pos, uni.gaps[[j]])
      }
    }
    remove_positions[[i]] = pos
    print(pos+star.position)
    # get all sequences that cause those gaps
    for (j in seq(1,length(gappositions.list))){
      gap.positions = unlist(gregexpr( pattern = "[-N]" , seg.fasta[j,"cleanseq"]))
      # check if any of the gaps is not! in pos
      if(sum(is.element(pos,gap.positions))!=length(pos)){
        remove = append(remove, as.character(seg.fasta[j, "strain"]))
      }
    }
    f[[i]] = seg.fasta
  }
  
  for (i in seq(1,length(f))){
    f[[i]] = f[[i]][-which(is.element(f[[i]]$strain, remove)),]
  }
  t.red = t.red[-which(is.element(t.red$strain, remove)),]
  
  # remove all gap positions from seq.fasta
  for (i in seq(1,length(f))){
    f[[i]]$nogaps=NA
    for (j in seq(1,length(f[[i]]$strain))){
      str.val = as.character(f[[i]][j, "cleanseq"])
      for (k in seq(length(remove_positions[[i]]), 1)){
        substring(str.val, remove_positions[[i]][[k]], remove_positions[[i]][[k]]) = ":"
      }
      f[[i]][j, "nogaps"] = gsub(":","",str.val)
    }
  }

  for (i in seq(1,length(segments))){
    con = file(paste(directory, "/cleaned/sequences_",virusname, "_", segments[[i]], ".cleaned.afasta", sep=""), "w")
    for (j in seq(1,length(t.red$strain))){
      ind = which(f[[i]]$strain==as.character(t.red[j,"strain"]))
      write(paste(">",t.red[j, "strain"], "",sep=""), file=con, append=TRUE)
      write(paste(f[[i]][ind, "nogaps"], "",sep=""), file=con, append=TRUE)
    }
  }
  close(con)
  
  # print meta data
  metafile = paste(directory, "/cleaned/metadata_",virusname,".cleaned.tsv", sep="")
  write.table(t.red, file=metafile, quote=FALSE, sep='\t', col.names = NA)
  
  
  return (f)
}