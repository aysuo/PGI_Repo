#!/usr/bin/env Rscript
library(data.table)

args=commandArgs(trailingOnly=TRUE)

inpath <- args[1]
outpath <- args[2]

data <- fread(inpath, header=T, sep="\t")
data$Z <- qnorm(1 - data$P/2)*sign(data$BETA)
data$SE <- data$BETA/data$Z
data$Z <- NULL

write.table(data,outpath, row.names=F, quote=F, sep="\t")
