#!/usr/bin/env Rscript
library(breakpointR)
library(dplyr)
library(tidyr)

get_states <- function(bpr_object){
    df_obj <- as.data.frame(bpr_object$counts[width(bpr_object$counts) >= 10000000]) %>%
        select(seqnames, states) %>%
        group_by(seqnames) %>%
        summarize(state = ifelse(n_distinct(states) == 1, first(states), NA)) %>%
        pivot_wider(names_from = seqnames, values_from = state) %>%
        mutate(library = strsplit(bpr_object$ID, "\\.")[[1]][1]) %>%
        select(last_col(), sort(everything()))
        return(df_obj)
}

args <- commandArgs(T)
output_file <- args[1]

filelist <- list.files("BPR_output/data", pattern="RData$", full.names=TRUE)
rdata <- loadFromFiles(filelist)

Library <- unname(unlist(lapply(rdata, function(x) strsplit(x$ID, "\\.")[[1]][1])))
summary_states <- bind_rows(unname(lapply(rdata, get_states)))
write.table(summary_states, output_file, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")