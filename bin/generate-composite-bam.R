#!/usr/bin/env Rscript
library(rtracklayer)
library(GenomicAlignments)

args <- commandArgs(T)
dir <- args[1]
ww_output <- args[2]
wc_output <- args[3]

load(paste0(dir,"/WW_composite_file.RData"))
load(paste0(dir,"/WC_composite_file.RData"))

export(WW_composite_file, BamFile(ww_output))
export(WC_composite_file, BamFile(wc_output))
