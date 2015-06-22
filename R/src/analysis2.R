#!/usr/bin/env Rscript
# exploration of metabolite data

rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis2"

load("./R/objects/peptides.matrix.combat.RData")
load("./R/objects/metabolites.data.RData")
load("./R/objects/metabolites.folds.RData")
load("./R/objects/protein_annotations.RData")
load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/kinase_classes._clean_.RData")


orf2name = unique(data.frame(ORF = protein_annotations$SystName, 
                             gene_name = protein_annotations$sgdName))





p = ggplot(metabolites.data, aes(x=batch, y=log(value,2))) +
           geom_boxplot()
plots.list = lappend(plots.list, p)

metabolites.folds = metabolites.folds %>% group_by(id) %>% mutate(median=median(value, na.rm=T))
metabolites.folds = metabolites.folds[with(metabolites.folds, order(median)),]
metabolites.folds$id = factor(metabolites.folds$id, levels=as.character(unique(metabolites.folds$id)))


p = ggplot(metabolites.folds, aes(x=id, y=log(value,2))) +
           geom_boxplot() + 
           coord_flip() +
           ylim(-2,4)
           
plots.list = lappend(plots.list, p)


metabolites.folds.matrix = dcast(metabolites.folds, formula=id~variable, value.var="value")
rownames(metabolites.folds.matrix) = metabolites.folds.matrix$id
metabolites.folds.matrix = as.matrix(metabolites.folds.matrix[,-1])

dissimilarity = as.dist(1 - cor(metabolites.folds.matrix,use="pairwise.complete.obs", method="spearman"), diag=T, upper=T)
pheatmap(dissimilarity)
p = recordPlot()
plots.list = lappend(plots.list, p)


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 



toPlot = log(metabolites.folds.imputed.matrix)
s1 = prcomp(toPlot)
toPlot = toPlot[!rownames(toPlot) %in% names(c(which(s1$x[,1] == max(s1$x[,1])), which(s1$x[,2] == max(s1$x[,2])))),]
s1 = prcomp(toPlot)

xPC = 1
yPC = 2
biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], 
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))

xPC = 1
yPC = 3
biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], 
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))

xPC = 2
yPC = 4
biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], 
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))


## -- PCA of metabolites ---- 
library("MASS")
library("Amelia")

#imputing missing values
#TODO: check validity of imputation

metabolites.folds.matrix.t = t(metabolites.folds.matrix)
metabolites.folds.imputed = amelia(metabolites.folds.matrix.t, logs=colnames(metabolites.folds.matrix.t), m=25)
metabolites.folds.imputed.matrix = Reduce("+",metabolites.folds.imputed$imputations)/length(metabolites.folds.imputed$imputations)

#removing outliers
metabolites = log(metabolites.folds.imputed.matrix)
metabolites.cov = cov.rob(metabolites)
md = mahalanobis(metabolites, center=metabolites.cov$center, cov=metabolites.cov$cov)
n = nrow(metabolites); p=ncol(metabolites)

plot(qchisq((1:n)/(n+1),p), sort(md), 
     xlab = expression(paste(chi^2," quantiles")),
     ylab = "Sorted Machalanobis distances")
abline(0,1)
p = recordPlot()
plots.list = lappend(plots.list, p)

#decided to remove 4 points
metabolites.f = metabolites
metabolites.f = metabolites[!(rownames(metabolites) %in% names(sort(-md)[1:4])),]
metabolites.f = metabolites.f[rownames(metabolites.f) %in% droplevels(exp_metadata$ORF[exp_metadata$type == "Kinase"]),] #taking only kinase deletions

toPlot = metabolites.f
s1 = prcomp(toPlot)
# toPlot = toPlot[!rownames(toPlot) %in% names(c(which(s1$x[,1] == max(s1$x[,1])), which(s1$x[,2] == max(s1$x[,2])))),]
# s1 = prcomp(toPlot)
par(mfrow=c(1,2))
xPC = 1
yPC = 2


biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)],
       xlabs=orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)],
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
abline(h=0,v=0)

xPC = 1
yPC = 3
biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], 
       xlabs=orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)],
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
abline(h=0,v=0)
p = recordPlot()
plots.list = lappend(plots.list, p)

# xPC = 2
# yPC = 3
# biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], 
#        xlabs=orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)],
#        xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
#        ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
# abline(h=0,v=0)

load("./R/objects/proteins.FC.RData")
load("./R/objects/peptides.ratios.RData")


proteins.FC.wide = dcast(proteins.FC,KO~ORF, value.var="logFC")
proteins.FC.matrix = as.matrix(proteins.FC.wide[,-1])
rownames(proteins.FC.matrix) = proteins.FC.wide$KO

proteins.FC.matrix.f = proteins.FC.matrix[match(rownames(metabolites.f), rownames(proteins.FC.matrix)),]

peptides.ratios.wide = dcast(peptides.ratios, KO~peptide, value.var="ratio.FC")
peptides.ratios.matrix = as.matrix(peptides.ratios.wide[,-1])
rownames(peptides.ratios.matrix) = peptides.ratios.wide$KO

peptides.ratios.matrix.f = peptides.ratios.matrix[match(rownames(metabolites.f), rownames(peptides.ratios.matrix)),]

hist(log(peptides.ratios.matrix.f))
prot.pca = prcomp(proteins.FC.matrix.f)
ratios.pca = prcomp(log(peptides.ratios.matrix.f))

plot(ratios.pca$x[,c(1,2)])



summary(lm(exp(metabolites.f)~prot.pca$x[,c(1:15)]))
prot.pca$x

cor(prot.pca$x[,1:10], metabolites.f)


a = hclust(dist(proteins.FC.matrix.f))
b = hclust(dist(metabolites.f))
c = hclust(dist(log(peptides.ratios.matrix.f)))

dendextend::cor_cophenetic(a,c)


