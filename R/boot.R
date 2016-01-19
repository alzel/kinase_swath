#library("org.Sc.sgd.db")
library("ggplot2")
library("reshape")
library("reshape2")
#library("data.table")
library("igraph")
#for data normalization
#library("preprocessCore")
library("DESeq")
#library("vsn")
library("limma")
#install.packages(c("NMF","grid", "NbClust", "dplyr","plyr", "xlsx"))
library("pheatmap")
library("NMF")
library("NbClust")
library("plyr"); library("dplyr")
library("grid")
library("gridExtra")
library("sva")
#library("xlsx")
library("tidyr")

library("RColorBrewer")
#library("ConsensusClusterPlus")


output_dir = "./R/objects"
figures_dir = "./figures"

dir.create(output_dir)
dir.create(figures_dir)
