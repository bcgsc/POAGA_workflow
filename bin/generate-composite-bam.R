#!/usr/bin/env Rscript
library(rtracklayer)
library(GenomicAlignments)

args <- commandArgs(T)
dir <- args[1]
output <- args[2]

load(paste0(dir,"/WW_composite_file.RData"))
load(paste0(dir,"/WC_composite_file.RData"))

export(WW_composite_file, BamFile(paste0(output, "WW_composite.bam")))
export(WC_composite_file, BamFile(paste0(output, "WC_composite.bam")))
