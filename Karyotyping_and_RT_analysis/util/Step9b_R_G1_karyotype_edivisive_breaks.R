#!/usr/bin/env Rscript

library(AneuFinder)

args <- commandArgs(TRUE)
ediv_models_dir <- args[1]
out_dir <- args[2]
granges <- args[3]


ediv_models_list <- list.files(ediv_models_dir,
    pattern = ".Rdata",
    full.names = TRUE
)

for (file in ediv_models_list) {
    name <- basename(file)
    name <- sub("_500kb_edivisive.Rdata", "", name)
    grange <- paste0(granges, "/", name, "_mapq10_blacklist_fragment.Rdata")
    fragment <- loadFromFiles(grange)
    hmm <- NULL
    hmm <- loadFromFiles(file)
    bp <- getBreakpoints(hmm, fragment, confint = 0.95)
    bpr <- refineBreakpoints(hmm, fragment, bp, confint = 0.95)
    save(bpr, file = paste0(out_dir, "/", name, "_500kb_CNV_breakpoints.Rdata"))
    exportCNVs(bpr,
        filename = paste0(out_dir, "/", name, "_500kb"),
        export.CNV = TRUE, export.breakpoints = TRUE, cluster = FALSE
    )
}
