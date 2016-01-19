#!/usr/bin/env Rscript
# Metabolite data exploration and figures for paper

rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis2d"

load("./R/objects/dataPPP_AA.imputed.create_datasets.RData")
load("./R/objects/metabolite2iMM904._load_.RData")
load("./R/objects/gene.annotations._load_.RData")
load("./R/objects/exp_metadata._clean_.RData")


orf2name = droplevels(unique(gene.annotations[,c("V4", "V6")]))
orf2name$V4 = as.character(orf2name$V4)
orf2name$V6 = as.character(orf2name$V6)
orf2name$V6[orf2name$V6 == ""] = orf2name$V4[orf2name$V6 == ""]
names(orf2name) = c("ORF", "gene_name")

metabolite.data.imputed = dataPPP_AA.imputed$metabolites

metabolite.matrix = metabolite.data.imputed
metabolite.matrix = metabolite.matrix[,which(!(colnames(metabolite.matrix) %in% c("AMP", "ADP", "ATP")))]


#removing outliers
metabolite.matrix.cov = cov.rob(metabolite.matrix)
metabolite.matrix.md = mahalanobis(metabolite.matrix, center=metabolite.matrix.cov$center, cov=metabolite.matrix.cov$cov)
n = nrow(metabolite.matrix); p=ncol(metabolite.matrix)

plot(qchisq((1:n)/(n+1),p), sort(metabolite.matrix.md), 
     xlab = expression(paste(chi^2," quantiles")),
     ylab = "Sorted Machalanobis distances")
abline(0,1)

p = recordPlot()
plots.list = lappend(plots.list, p)



#decided to remove 3 points
metabolite.matrix.f = metabolite.matrix
metabolite.matrix.f = metabolite.matrix.f[!(rownames(metabolite.matrix.f) %in% names(sort(-metabolite.matrix.md)[1:3])),]

metabolite.matrix.f = metabolite.matrix.f[rownames(metabolite.matrix.f) %in% c("WT", as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),] #taking only kinase deletions

toPlot = scale(metabolite.matrix.f)
s1 = prcomp(toPlot)

# toPlot = toPlot[!rownames(toPlot) %in% names(c(which(s1$x[,1] == max(s1$x[,1])), which(s1$x[,2] == max(s1$x[,2])))),]
# s1 = prcomp(toPlot)
xPC = 1
yPC = 2
labels = orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)]
labels[is.na(labels)] = "WT"

biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], cex=0.66,
       xlabs=labels,
       ylabs=toupper(metabolite2iMM904$model_name[match(rownames(s1$rotation[,c(xPC,yPC)]), metabolite2iMM904$id)]),
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))

abline(h=0,v=0)
p = recordPlot()
plots.list = lappend(plots.list, p)


xPC = 2
yPC = 3
biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], cex=0.66,
       xlabs=labels,
       ylabs=toupper(metabolite2iMM904$model_name[match(rownames(s1$rotation[,c(xPC,yPC)]), metabolite2iMM904$id)]),
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
abline(h=0,v=0)

p = recordPlot()
plots.list = lappend(plots.list, p)


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")

# ### PCA plot for paper
# 
# S = prcomp(scale(metabolite.matrix.f))
# 
# toPlot = as.data.frame(S$x[,c(1:5)])
# 
# 
# toPlot$labels = orf2name$gene_name[match(rownames(toPlot), orf2name$ORF)]
# 
# round(S$sdev[1]^2/sum(S$sdev^2),2)
# round(S$sdev[2]^2/sum(S$sdev^2),2)
# 
# 
# p1 = ggplot(toPlot, aes(x=PC1, y=PC2, label=labels)) +
#   geom_point() +
#   geom_text(hjust=0, just=0, col="blue") +
#   xlab(paste("PC1,",round(S$sdev[1]^2/sum(S$sdev^2),2))) +
#   ylab(paste("PC2,",round(S$sdev[2]^2/sum(S$sdev^2),2))) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   theme(aspect.ratio = 1) +
#   panel_border()
# 
# 
# p2 = ggplot(toPlot, aes(x=PC2, y=PC3, label=labels)) +
#   geom_point() +
#   geom_text(hjust=0, just=0, col="blue") +
#   xlab(paste("PC2,",round(S$sdev[2]^2/sum(S$sdev^2),2))) +
#   ylab(paste("PC3,",round(S$sdev[3]^2/sum(S$sdev^2),2))) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   theme(aspect.ratio = 1) +
#   panel_border()
# 
# 
# toPlot = as.data.frame(S$rotation[,1:5])
# 
# p2 = ggplot(proteins.FC.f.stats, aes(x=gene_name, y=n_metabolic, fill=as.factor(part))) + 
#   geom_bar(stat="identity", width=.5) + 
#   coord_flip() +
#   theme_classic()
# 
# 
# plot_grid(p1, p2, labels = c("A", "B"))



