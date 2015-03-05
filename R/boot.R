library("ggplot2")
library("reshape")
library("reshape2")
library("data.table")
#library("plyr")
library("igraph")
#for data normalization
library("preprocessCore")
library("DESeq")
library("vsn")
library("limma")

library("pheatmap")
library("NMF")
library("NbClust")
library(plyr); library(dplyr)
library("grid")
library("sva")
library("xlsx")

output_dir = "./R/objects"
figures_dir = "./figures"

dir.create(output_dir)
dir.create(figures_dir)
