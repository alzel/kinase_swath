rm(list = ls())
library(caret)
library(relaimpo)
source("./R/boot.R")
library("cowplot")
library("scales")

#summarizes results from linear regression models

load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")
input_path = "./results/2015-08-03/linear_models"
load("./R/objects/proteins.matrix.combat.quant.RData")

load("./R/objects/gene.annotations._load_.RData")
orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]


fun_name = "linear_models_summary"

filesToProcess = dir(path=input_path, pattern = "[12].[01]+.linear_models.RData$", recursive=F)
filesToProcess = grep(pattern="imputed", invert=T, filesToProcess, value=T)

pattern.p = "data.(\\w+).(.*?).([0-9]+).([0-9]+).linear_models.RData$"

matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

read_models.caret = function(x) {
  z <<- x
  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))
  
  if (length(my_models$caret_models) > 1 ) {
    t = resamples(my_models$caret_models)
    t.summary = summary(t)
    fit = my_models$caret_models[[which.min(t.summary$statistics$RMSE[,"Mean"])]]$finalModel
  } else {
    fit  = my_models$caret_models[[1]]$finalModel
  }
  fit.s = summary(fit)
  
  if (length(fit$coefficients) > 2 ) {
    table = data.frame(feature = names(calc.relimp(fit, type="lmg")@lmg),
                       importance = calc.relimp(fit, type="lmg")@lmg,
                       row.names=NULL)
    table$betas = fit$coefficients[match(table$feature,  names(fit$coefficients))]
    
  } else if (length(fit$coefficients) == 2) {
    
    table = data.frame(feature = names(fit$coefficients)[2],
                       importance = fit.s$r.squared, row.names=NULL,
                       betas = fit$coefficients[2])
  } else {
    table = data.frame(feature = NA,
                       importance = NA,
                       betas = NA)
  } 
    
  table$r.squared = fit.s$r.squared
  table$adj.r.squared = fit.s$adj.r.squared
  table$dataset = factor(x[[2]])
  table$metabolite  = factor(x[[3]])
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = factor(x[[4]])
  table$ismetIncluded  = factor(x[[5]])
  table$file =  factor(x[[1]])
  return(table)
}


read_models.aic = function(x) {
  z <<- x
  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))

  table = my_models$sample_models
  table$dataset = factor(x[[2]])
  table$metabolite  = factor(x[[3]])
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = factor(x[[4]])
  table$ismetIncluded  = factor(x[[5]])
  table$file =  factor(x[[1]])
  return(table)
}



file.list = lapply(matches, FUN=read_models.caret)
all_final_models.caret = do.call(rbind.data.frame, file.list)

file.list = lapply(matches, FUN=read_models.aic)
all_final_models.aic = do.call(rbind.data.frame, file.list)

all_final_models.aic.summary = all_final_models.aic %>% group_by(formula, dataset, metabolite, isImputed, degree, ismetIncluded) %>% summarize(r.squared = mean(r.squared),
                                                                                                                                               adj.r.squared = mean(adj.r.squared),
                                                                                                                                               aic = median(aic))

formulas.coef = all_final_models.aic %>% dplyr::select(formula, coeficients) %>% distinct()
all_final_models.aic.merged = merge(all_final_models.aic.summary, formulas.coef, by="formula")
all_final_models.aic.merged = droplevels(all_final_models.aic.merged[!all_final_models.aic.merged$coeficients == "(Intercept)",])



all_final_models.caret$isMetabolite = ifelse(all_final_models.caret$feature %in% metabolite2iMM904$id,1,0)
all_final_models.caret = droplevels(filter(all_final_models.caret, dataset != "PPP_AA"))

      
## getting coverages
yeast.model = iMM904
yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters
edgelist = unique(droplevels(na.omit(subset(yeast.model, metabolite != "h"  & metabolite !="h2o" & metabolite != "nh4" & metabolite != "ppi" , select = c("metabolite", "gene")))))
B <- graph.data.frame(edgelist)
V(B)$type <- V(B)$name %in% edgelist$metabolite
stopifnot(is.bipartite(B))


measured.proteins = row.names(proteins.matrix.combat.quant)

coverage = ddply(all_final_models.caret, .(dataset, metabolite), 
      .fun=function(x) {
  
        dataset = as.character(unique(x$dataset))
        i = as.character(unique(x$metabolite))
        
        current_nodes = as.character(metabolite2iMM904$model_name[metabolite2iMM904$id == i])    
        measured.metabolites = as.character(unique(all_final_models.caret$metabolite[all_final_models.caret$dataset == dataset ]))
        measured.metabolites.model = as.character(unique(metabolite2iMM904$model_name[metabolite2iMM904$id %in% measured.metabolites]))
  
        if (length(current_nodes) == 0) {
          message(paste("No ", i, "found in network"))
          next()
        }
        
        SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=2, nodes=current_nodes))))    
        genes = unique(V(SUB)$name[V(SUB)$type == 0])
        metabolites = unique(V(SUB)$name[V(SUB)$type == 1])
        
        coverage.met = data.frame(metabolite = i, 
                                  dataset = dataset, 
                                  type = "metabolic",
                                  neighbour = metabolites,
                                  isMeasured = ifelse(metabolites %in% measured.metabolites.model, 1, 0))
        
        coverage.gene = data.frame(metabolite = i, 
                                  dataset = dataset, 
                                  type = "gene",
                                  neighbour = genes,
                                  isMeasured = ifelse(genes %in% measured.proteins, 1, 0))
        
        return(rbind(coverage.met, coverage.gene))
        
      })  
  


coverage.stats = coverage %>% group_by(metabolite, dataset, type) %>% summarise(coverage = sum(isMeasured)/length(isMeasured))
coverage.stats.wide = dcast(formula=metabolite+dataset~type, value.var="coverage", data=coverage.stats)
  
ggplot(coverage.stats.wide, aes(x=gene, y=metabolic)) + geom_abline(intercept=0, slope=1)+
      xlim(0,1)+
      ylim(0,1)+
      geom_point()



control.summary.wide = dcast(data=control.summary.long, formula=metabolite+dataset~variable, value.var="value")

control.stats = merge(control.summary.wide, coverage.stats.wide, by=c("metabolite", "dataset"))

ggplot(control.stats, aes(y=rela.metabolic, x=metabolic)) + geom_point()

control.stats$labels = metabolite2iMM904$model_name[match(control.stats$metabolite, metabolite2iMM904$id)]
p = ggplot(filter(control.stats, dataset=="TCA_AA"), aes(y=rela.gene, x=gene)) + 
    geom_point() + 
    geom_text(aes(y=rela.gene,size=1, hjust=0, vjust=0,  x=gene, label=labels))+
    ylab(expression(paste("Relative importance as a fraction of model's R"^"2", sep=""))) +
    xlab(expression("Fraction of measured metabolic enzymes")) +
    scale_y_continuous(labels = percent, limits = c(0,1)) +
    scale_x_continuous(labels = percent, limits = c(0,1)) +
    theme(legend.position = "none",
          axis.title = element_text(size = rel(1.5)), 
          aspect.ratio = 1) +
    background_grid(major = "xy", minor = "none")

file_name = paste(fun_name,"gene_coverage.vs.importance", "png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width = 8.27, dpi=300)



control.summary = all_final_models.caret %>% group_by(metabolite, dataset, isImputed, degree, ismetIncluded) %>% 
  summarize(metabolic = sum(importance[isMetabolite == 1]),
            gene = sum(importance[isMetabolite == 0]),
            rela.metabolic = metabolic/sum(importance),
            rela.gene = gene/sum(importance),
            r.squared = r.squared[1], 
            adj.r.squared = adj.r.squared[1])


control.summary.long = melt(as.data.frame(control.summary %>% ungroup %>% filter(isImputed == 0, degree==2, ismetIncluded ==1) %>% 
                                            dplyr::select(metabolite, dataset, rela.metabolic, rela.gene)), id.vars=c("metabolite", "dataset"))


control.summary.merged = merge(control.summary.long , control.summary %>% ungroup %>% filter(isImputed == 0, degree==2, ismetIncluded ==1)  %>% 
                              select(metabolite, dataset, r.squared, adj.r.squared),  by=c("metabolite","dataset"))

dtst = "TCA_AA"
toPlot = droplevels(control.summary.merged %>% filter(dataset != dtst))
toPlot = toPlot %>% ungroup() %>% arrange(variable, adj.r.squared)
toPlot$metabolite = factor(toPlot$metabolite, levels = unique(as.character(toPlot$metabolite)))
toPlot$label = metabolite2iMM904$official_name[match(toPlot$metabolite, metabolite2iMM904$id)]
toPlot$label = factor(toPlot$label, levels = unique(as.character(toPlot$label)))

#no-barplot
p = ggplot(toPlot, aes(x=label, y=value)) + 
  geom_bar(stat="identity",  aes(fill=variable)) + 
  geom_point(data = toPlot , aes(x=label, y= adj.r.squared)) +
  geom_linerange(data = toPlot , aes(x=label,ymin=0, ymax=adj.r.squared)) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  xlab("Metabolite model") +
  ylab("Relative importance of features") + 
  #scale_fill_manual("Features",values=c(rela.metabolic="lightgreen",
  #                                     rela.gene="lightblue") )+
  theme(axis.text.y = element_text(size = rel(0.75)),
        legend.position = "none") +
  coord_flip()


file_name = paste(fun_name,"needleplot.adj.r2", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width = 8.27)


p = ggplot(toPlot, aes(x=label, y=value)) + 
  geom_bar(stat="identity",  aes(fill=variable)) + 
  geom_point(data = toPlot , aes(x=label, y= adj.r.squared)) +
  geom_linerange(data = toPlot , aes(x=label,ymin=0, ymax=adj.r.squared)) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  xlab("Metabolite model") +
  ylab("Relative importance of features") + 
  #scale_fill_manual("Features",values=c(rela.metabolic="lightgreen",
  #                                      rela.gene="lightblue") )+
  theme(axis.text.y = element_text(size = rel(0.75)),
        legend.position = "none") +
  coord_flip()


file_name = paste(fun_name,"importance.needleplot.adj.r2", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width = 8.27)




toPlot2 = droplevels(dplyr::filter(control.stats, dataset != dtst))
toPlot2$metabolite = factor(toPlot2$metabolite, levels = unique(as.character(toPlot$metabolite)))

p = ggplot(toPlot2, aes(x=metabolite)) +
        geom_bar(aes(y=-metabolic, fill="Metabolic"), stat="identity" ) +
        geom_bar(aes(y=gene, fill="Enzymatic"), stat="identity") +
        ylab("Fraction of measured metabolic/enzyme neighbours") +
        scale_fill_manual("Feature",values=c(Metabolic="lightgreen",
                                              Enzymatic="lightblue") )+
        scale_y_continuous(labels = percent) +
        theme(axis.text.y = element_text(size = rel(0.75))) +
        coord_flip()

file_name = paste(fun_name,"coverage", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width = 8.27)



graph_dataset = all_final_models.caret %>% filter(dataset=="TCA_AA", ismetIncluded == 1, degree == 2)
graph_dataset$effector = ifelse(graph_dataset$betas > 0, 1, 0)
edgelist = graph_dataset %>% dplyr::select(feature, metabolite, importance, effector, betas)

B <- graph.data.frame(edgelist)
V(B)$type <- V(B)$name %in% edgelist$metabolite
V(B)$rsquared <- graph_dataset$r.squared[match(V(B)$name, graph_dataset$metabolite)]

name_idx = na.omit(match(V(B)$name, orf2name$ORF))
B_idx = which(!is.na(match(V(B)$name, orf2name$ORF)))

V(B)$label = 1:length(V(B)$name)
V(B)$label[B_idx] = as.character(orf2name$gene_name[name_idx])


name_idx = na.omit(match(V(B)$name, metabolite2iMM904$id))
B_idx = which(!is.na(match(V(B)$name, metabolite2iMM904$id)))
V(B)$label[B_idx] = as.character(metabolite2iMM904$model[name_idx])

write.graph(graph=B,file="test.graphml",format="graphml")




heatmap = read.delim("./data/2015-08-17/heatmap.txt", header=T)
heatmap$cluster = NULL
metInDataset = filter(all_final_models.caret, dataset == "TCA_AA") %>% dplyr::select(metabolite) %>% distinct()

yeast.model = iMM904
yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters
edgelist = unique(droplevels(na.omit(subset(yeast.model, metabolite != "h"  & metabolite !="h2o" , select = c("metabolite", "gene")))))

B <- graph.data.frame(edgelist)
V(B)$type <- V(B)$name %in% edgelist$metabolite



current_nodes = as.character(metabolite2iMM904$model_name[match(as.character(metInDataset$metabolite), metabolite2iMM904$id)])
all_final_models.caret$isMetabolite = ifelse(all_final_models.caret$feature %in% metabolite2iMM904$id,1,0)

SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=1, nodes=current_nodes)))) 

#enrichment of essential genes
load("./R/objects/essential_ORFs._load_.RData")
coverage$isEssential = ifelse(coverage$neighbour %in% unique(essential_ORFs$ORF_name),1,0)


coverage$isImportant = ifelse(coverage$neighbour %in% all_final_models.caret$feature,1,0)


genes = unique(filter(coverage, type == "gene" ) %>% select(-metabolite, -dataset))
test.matrix = matrix(c(nrow(filter(genes, isImportant==1, isEssential==1)),
                       nrow(filter(genes, isImportant==0, isEssential==1)),
                       nrow(filter(genes, isImportant==1, isEssential==0)),
                       nrow(filter(genes, isImportant==0, isEssential==0))),nrow = 2,
                       dimnames = list(inModel  = c("Y", "N"),
                                       isEssential = c("Y", "N")))
fisher.test(test.matrix)

heatmap.long = melt(heatmap, id.vars=c("ORF", "gene"))
coverage.heatmap.merged = merge(filter(coverage, type == "gene"), heatmap.long,   
                                       by.x = c("metabolite", "neighbour"), by.y = c("variable", "ORF"))


coverage.heatmap.merged$isImportant = coverage.heatmap.merged$neighbour %in% all_final_models.caret$feature

ggplot(droplevels(coverage.heatmap.merged), aes(x=metabolite, y=log(abs(value)), colour=isImportant, shape = ifelse(isMeasured==1, "a", "b"))) + 
       geom_point() + 
       facet_wrap(~dataset, scales="free")

ggplot(droplevels(coverage.heatmap.merged, dataset=="AA"), aes(y=(value), x=isImportant)) + 
        geom_boxplot(alpha=0.1) + ylim(-100,100)
        
dplyr::select(all_final_models.caret, feature, importance, betas, metabolite, dataset, isImputed, degree, ismetIncluded)
control.heatmap.merged = merge(dplyr::select(all_final_models.caret, feature, importance, betas, metabolite, dataset,
                                             isImputed, degree, ismetIncluded), heatmap.long,   
                                by.x = c("metabolite", "feature"), by.y = c("variable", "ORF"))


## F6P example ####

yeast.model = iMM904
yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters

edgelist = unique(droplevels(na.omit(subset(yeast.model, metabolite != "h"  & metabolite !="h2o" , select = c("metabolite", "gene")))))


measured.metabolites = as.character(unique(all_final_models.caret$metabolite[all_final_models.caret$dataset == "TCA" ]))
measured.metabolites.model = as.character(unique(metabolite2iMM904$model_name[metabolite2iMM904$id %in% measured.metabolites]))


edgelist = edgelist[edgelist$gene %in% measured.proteins,]
#edgelist = edgelist[edgelist$metabolite %in% measured.metabolites.model,]

B <- graph.data.frame(edgelist)
V(B)$type <- V(B)$name %in% edgelist$metabolite
current_nodes = as.character(metabolite2iMM904$model_name[metabolite2iMM904$id %in% "ADP"])

name_idx = na.omit(match(V(B)$name, orf2name$ORF))
B_idx = which(!is.na(match(V(B)$name, orf2name$ORF)))
V(B)$label = 1:length(V(B)$name)
V(B)$label[B_idx] = as.character(orf2name$gene_name[name_idx])

name_idx = na.omit(match(V(B)$name, metabolite2iMM904$model_name))
B_idx = which(!is.na(match(V(B)$name, metabolite2iMM904$model_name)))
V(B)$label[B_idx] = as.character(metabolite2iMM904$official_name[name_idx])
adp_model = c("YDL185W", "YER052C", "YGL062W","YGL253W", "YGR204W", "YGR240C", "YJL026W", "YJL130C", "YJR109C", "YJR121W", "YKL001C")

SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=1, nodes=current_nodes)))) 
V(SUB)$model = ifelse(V(SUB)$name %in% adp_model, 1, 0)

adp_model = c("YDL185W", "YER052C", "YGL062W","YGL253W", "YGR204W", "YGR240C", "YJL026W", "YJL130C", "YJR109C", "YJR121W", "YKL001C")

orf2name$gene_name[match(adp_model, orf2name$ORF)]


file_name = paste(fun_name,"adp_example", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
pdf(file=file_path, height=8.27, width = 8.27)
plot.igraph(x=SUB, vertex.color=V(SUB)$model)
dev.off()


edgelist = edgelist[edgelist$gene %in% measured.proteins,]
edgelist = edgelist[edgelist$metabolite %in% measured.metabolites.model,]

B <- graph.data.frame(edgelist, directed="F")
V(B)$type <- V(B)$name %in% edgelist$metabolite
current_nodes = as.character(metabolite2iMM904$model_name[metabolite2iMM904$id %in% "F6P"])

name_idx = na.omit(match(V(B)$name, orf2name$ORF))
B_idx = which(!is.na(match(V(B)$name, orf2name$ORF)))
V(B)$label = 1:length(V(B)$name)
V(B)$label[B_idx] = as.character(orf2name$gene_name[name_idx])

name_idx = na.omit(match(V(B)$name, metabolite2iMM904$model_name))
B_idx = which(!is.na(match(V(B)$name, metabolite2iMM904$model_name)))
V(B)$label[B_idx] = as.character(metabolite2iMM904$official_name[name_idx])

F6P_model = c("YGL253W", "YKL104C", "YPR074C")

SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=2, nodes=current_nodes)))) 
V(SUB)$model = ifelse(V(SUB)$name %in% F6P_model, 1, 0)

file_name = paste(fun_name,"f6p_example", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
pdf(file=file_path, height=8.27, width = 8.27)
plot.igraph(x=SUB, vertex.color=V(SUB)$model)
dev.off()





