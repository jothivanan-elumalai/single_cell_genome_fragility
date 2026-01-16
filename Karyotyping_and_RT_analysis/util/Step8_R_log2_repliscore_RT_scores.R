#!/usr/bin/env Rscript
#scRepli-seq pipeline
#Computation of log2 percentile score following mappability correction

library("pracma")
library("AneuFinder")

args = commandArgs(TRUE)
fragment=args[1]
outdir=args[2]
ref_Rdata=args[3]
blacklist=args[4]
genome_file=args[5]
repliscore_df=args[6]
path=args[7]

source(normalizePath(paste0(path,"/util/Aneufinder_Optional_script.R")))

window=200000
sliding=40000

name=basename(fragment)
name=sub("_mapq10_blacklist_fragment.Rdata", "", name)

ptm <- startTimedMessage(paste0(name,"\t","log2repliscore_RT_calculation.....","\n"))

rep_df=read.table(repliscore_df,sep="\t",header=T)
rep_df$Repliscore_without_X_rev = (100 - rep_df$Repliscore_without_X)/100
percentile_value <- rep_df[rep_df[, 1] == name, 4]


percentile_scale = function (datalist) {
x=datalist[[1]]
c=elementMetadata(x)$counts
prtl_val=quantile(c[c!=0],probs=c(percentile_value)) # changed by JE
log2rep=log2(c/prtl_val)
out=data.frame(chr=seqnames(x),
                 starts=start(x)-1,
                 ends=end(x),
                 map_percentile_scale_log2=log2rep)
return(out)
}


##loading black list and genome Info##
genome_tmp <- read.table(genome_file,sep="\t") #
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

lcm=Lcm(window,sliding)
options(scipen=100)
print(paste("sliding size is",lcm/sliding,"for window:",window,"interval:",sliding))
sliding_window_data=list()
for (j in 1:(lcm/sliding)){
  tmp_out=NULL
  for (chr in chromosomes){
    chr_size=genome$UCSC_seqlength[genome$UCSC_seqlevel==chr]
    s=(j-1)*sliding+window
    if (s>=chr_size){print(paste("skip",chr));next}
    tmp=data.frame(chr=chr,
                   start=seq(s,chr_size,window) - window + 1,
                   end=seq(s,chr_size,window))
    tmp_out=rbind(tmp_out,tmp)
  }
  sliding_window_data[[j]]=tmp_out
}

x_gr=list()
for (i in 1:length(sliding_window_data)){
  x_gr[[i]]=makeGRangesFromDataFrame(sliding_window_data[[i]])
}

load(fragment)
test=binReads(raw_reads,
              assembly=genome,
              blacklist=blacklist,
              bins=x_gr)

ref_data = load(ref_Rdata)
control_100_S_G1_reads = get(ref_data)

test_map=list()
for (i in 1:(length(test)-1)){
test_map[[i]]=correctMappability(test[[i]],same.binsize = T,control_100_S_G1_reads,assembly=genome)
}

test_map_repliscore_log2=NULL
for (i in 1:length(test_map)){
  test_map_repliscore_log2=rbind(test_map_repliscore_log2,percentile_scale(test_map[[i]]))
}

out=test_map_repliscore_log2
a=out[order(out$chr,out$starts),]
a1=data.frame(chr=a[,1],s=round((a[,2]+a[,3])/2-sliding/2),e=round((a[,2]+a[,3])/2-sliding/2)+sliding,c=a$map_percentile_scale_log2)
a2=a1[!is.infinite(a1$c),]

write.table(a2,paste(outdir,"/",name,"_w200ks40k_map_count_repliscore_log2.bedGraph",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
