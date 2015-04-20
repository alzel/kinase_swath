#!/usr/bin/env Rscript
# exploration of metabolite data

rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis2"

load("./R/objects/proteins.matrix.f.deseq.combat.RData")
load("./R/objects/peptides.matrix.f.RData")

load("./R/objects/metabolites.data.RData")
load("./R/objects/metabolites.folds.RData")



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
p           
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


