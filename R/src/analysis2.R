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
load("./R/objects/metabolite2iMM904._load_.RData")
load("./R/objects/metabolitesTCA_metadata._clean_.RData")
load("./R/objects/gene.annotations._load_.RData")

library(cowplot)

orf2name = droplevels(unique(gene.annotations[,c("V4", "V6")]))
orf2name$V4 = as.character(orf2name$V4)
orf2name$V6 = as.character(orf2name$V6)
orf2name$V6[orf2name$V6 == ""] = orf2name$V4[orf2name$V6 == ""]
names(orf2name) = c("ORF", "gene_name")


p = ggplot(metabolites.data, aes(x=batch, y=log(value,2))) +
           geom_boxplot()


plots.list = lappend(plots.list, p)

metabolites.folds = metabolites.folds %>% group_by(id) %>% mutate(median=median(value, na.rm=T))
metabolites.folds = metabolites.folds[with(metabolites.folds, order(median)),]

metabolites.folds$name = metabolite2iMM904$official_name[match(metabolites.folds$id, metabolite2iMM904$id)]
metabolites.folds$name = factor(metabolites.folds$name, levels=as.character(unique(metabolites.folds$name)))

metabolites.folds$inTCA_method = metabolites.folds$variable %in% metabolitesTCA_metadata$ORF

p = ggplot(metabolites.folds, aes(x=name, y=log(value,2))) +
           geom_boxplot() + 
           geom_point(data=metabolites.folds %>% filter(inTCA_method == T), aes(x=name, y=log(value,2)), color="red") +
           coord_flip()  +
           ylim(c(-2, 4)) +
           xlab("Metabolite") +
           ylab("Log2(fold-change)")



file_name = paste(fun_name,"metabolites_folds_PPP", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width = 8.27)



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
metabolites.folds.matrix.t = metabolites.folds.matrix.t[,which(!colnames(metabolites.folds.matrix.t) %in% c("ADP", "AMP", "ATP"))]
metabolites.folds.imputed = amelia(metabolites.folds.matrix.t, logs=colnames(metabolites.folds.matrix.t), m=25)
metabolites.folds.imputed.matrix = Reduce("+",metabolites.folds.imputed$imputations)/length(metabolites.folds.imputed$imputations)

#removing outliers
metabolites = log(metabolites.folds.imputed.matrix)
metabolites.cov = cov.rob(metabolites)
md = mahalanobis(metabolites, center=metabolites.cov$center, cov=metabolites.cov$cov)
n = nrow(metabolites); p=ncol(metabolites)

names(md) %in% metabolitesTCA_metadata$ORF

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

dev.off()
biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)],
       xlabs=orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)],
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))

points(s1$x[rownames(s1$x) %in% metabolitesTCA_metadata$ORF,c(xPC,yPC)], col="blue")

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

### PCA plot for paper

S = prcomp(metabolites.f)
metabolites.f
toPlot = as.data.frame(S$x[,c(1:5)])



toPlot$labels = orf2name$gene_name[match(rownames(toPlot), orf2name$ORF)]
  
round(S$sdev[1]^2/sum(S$sdev^2),2)
round(S$sdev[2]^2/sum(S$sdev^2),2)


p1 = ggplot(toPlot, aes(x=PC1, y=PC2, label=labels)) +
        geom_point() +
        geom_text(hjust=0, just=0, col="blue") +
        xlab(paste("PC1,",round(S$sdev[1]^2/sum(S$sdev^2),2))) +
        ylab(paste("PC2,",round(S$sdev[2]^2/sum(S$sdev^2),2))) +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        theme(aspect.ratio = 1) +
        panel_border()



p2 = ggplot(toPlot, aes(x=PC2, y=PC3, label=labels)) +
        geom_point() +
        geom_text(hjust=0, just=0, col="blue") +
        xlab(paste("PC2,",round(S$sdev[2]^2/sum(S$sdev^2),2))) +
        ylab(paste("PC3,",round(S$sdev[3]^2/sum(S$sdev^2),2))) +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        theme(aspect.ratio = 1) +
        panel_border()


toPlot = as.data.frame(S$rotation[,1:5])

p2 = ggplot(proteins.FC.f.stats, aes(x=gene_name, y=n_metabolic, fill=as.factor(part))) + 
  geom_bar(stat="identity", width=.5) + 
  coord_flip() +
  theme_classic()


plot_grid(p1, p2, labels = c("A", "B"))


       
       