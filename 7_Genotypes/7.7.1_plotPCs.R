#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
args=commandArgs(trailingOnly=TRUE)
PCs=args[1]
out=args[2]

data <-fread(PCs,header=T)

pdf(file=out)  

PC12<-ggplot(data, aes(x=PC1, y=PC2, color=SUPERPOP)) + geom_point()
PC13<-ggplot(data, aes(x=PC1, y=PC3, color=SUPERPOP)) + geom_point()
PC14<-ggplot(data, aes(x=PC1, y=PC4, color=SUPERPOP)) + geom_point()
PC23<-ggplot(data, aes(x=PC2, y=PC3, color=SUPERPOP)) + geom_point()
PC24<-ggplot(data, aes(x=PC2, y=PC4, color=SUPERPOP)) + geom_point()
PC34<-ggplot(data, aes(x=PC3, y=PC4, color=SUPERPOP)) + geom_point()

print(PC12)
print(PC13)
print(PC14)
print(PC23)
print(PC24)
print(PC34)

dev.off()