#!/usr/bin/env Rscript
# edivisive by 500kb for checking the karyotype

library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm9)

args <- commandArgs(TRUE)
binfile <- args[1]
out_dir <- args[2]
genome <- args[3] # "hg19" or "mm9"

name <- basename(binfile)
name <- sub("_mapq10_blacklist_bin.Rdata", "", name)

load(binfile)

bsgenome <- if (genome == "hg19") {
    BSgenome.Hsapiens.UCSC.hg19
} else if (genome == "mm9") {
    BSgenome.Mmusculus.UCSC.mm9
} else {
    stop("Unsupported genome: ", genome)
}

bins_500k_gc <- correctGC(bins_reads$`binsize_5e+05`,
    GC.BSgenome = bsgenome,
    method = "loess"
)

bins_500k_ediv <- findCNVs(bins_500k_gc,
    ID = name,
    method = "edivisive",
    R = 50
)

save(bins_500k_ediv, file = paste(out_dir, "/", name, "_500k_edivisive.Rdata", sep = ""))

pdf(paste(out_dir, "/", name, "_500k_edivisive_histogram.pdf", sep = ""))
print(plot(bins_500k_ediv, type = "histogram"))
dev.off()
