#!/usr/bin/env Rscript
#scRepli-seq pipeline

library(AneuFinder)

options(scipen=100)
args = commandArgs(TRUE)
in_dir=args[1]
out_dir=args[2]
out_name=args[3]

# in_dir="/home/jothivanan_elumalai/screpliseq_test/Aneu_analysis/bins"
files=list.files(in_dir,pattern=".Rdata",full.names=T)

out=NULL
for (file in files){
#file="/home/jothivanan_elumalai/screpliseq_test/Aneu_analysis/bins/Sample_P285_09_1_R1_100k_mapq10_blacklist_bin.Rdata"
load(file) #bin file

#200k bins data
data=data.frame(chr=seqnames(bins_reads$`binsize_2e+05`),
                starts=start(bins_reads$`binsize_2e+05`)-1,
                ends=end(bins_reads$`binsize_2e+05`),
                elementMetadata(bins_reads$`binsize_2e+05`)$counts)

rmZero=data[,4] !=0
data_rmZero=data[rmZero,]
med=median(data_rmZero[,4])
med_scaled=data_rmZero[,4]/med
MAD=mad(med_scaled)
MAD_log=mad(log2(med_scaled))

name=basename(file)

out=rbind(out,c(name,MAD_log))
}

colnames(out)=c("Name","MAD_score_log2")

write.table(out,paste0(out_dir,"/",out_name,"_MAD_scores_log2.txt"),sep="\t",row.names=F,col.names=T,quote=F)
