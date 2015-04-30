rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis1"

## -- Load of data ----
load("./R/objects/proteins.FC.combat.RData")

load("./R/objects/pathway2orf.RData")
load("./R/objects/pathways.RData")
load("./R/objects/GO_slim.function.RData")
load("./R/objects/GO_slim.process.RData")
load("./R/objects/GO_slim.compartment.RData")

load("./R/objects/kinase_classes._clean_.RData")
load("./R/objects/exp_metadata._clean_.RData")

## -- clustering based on fold-changes ----

proteins.FC = proteins.FC.combat
reference = unique(as.character(proteins.FC$reference))

lb = -4
ub = 4
toPlot = proteins.FC
toPlot[toPlot$logFC < 0,]$logFC = ifelse(toPlot[toPlot$logFC < 0,]$logFC < lb, lb, toPlot[toPlot$logFC < 0,]$logFC)
toPlot[toPlot$logFC > 0,]$logFC = ifelse(toPlot[toPlot$logFC > 0,]$logFC > ub, ub, toPlot[toPlot$logFC > 0,]$logFC)

FC_thr = 1
pval_thr = 0.05
proteins.FC.stats = data.frame( ratio_sign = round(sum(proteins.FC$logFC < -FC_thr & proteins.FC$p.value_BH < pval_thr)/sum(proteins.FC$logFC > FC_thr & proteins.FC$p.value_BH < pval_thr),2),
                                ratio =  round(sum(proteins.FC$logFC < 0)/sum(proteins.FC$logFC > 0),2),
                                n_prot = length(unique(proteins.FC$ORF)),
                                n_sign = round(sum(proteins.FC$p.value_BH<0.05)/length(unique(proteins.FC$contrasts)),2),
                                x_min = -4,
                                y_max =max(-log10(proteins.FC$p.value_BH)))

toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)


p1 = ggplot(toPlot, aes(y=-log10(p.value_BH), x=logFC)) +
         geom_point(aes(color = sign), alpha=0.5) + 
         xlim(c(lb,ub)) +
         geom_hline(y=-log(pval_thr,10),linetype=3) +
         geom_vline(x=c(FC_thr,-FC_thr),linetype=3) +
         geom_text(data=proteins.FC.stats, aes(x=x_min, y=y_max, hjust=0, vjust=1,
                                             label=paste(paste0("#_mean_sign = ", n_sign),
                                                         paste0("up/down_ratio = ",ratio),
                                                         paste0("ratio_sign = ",   ratio_sign),
                                                         paste0("#_prot = ", n_prot), sep="\n"))) +
        xlab(paste("Log2(fold-change)")) +
        theme(aspect.ratio = 1, legend.position = "none", text=element_text(size=18))

p2 = ggplot(toPlot, aes(x=logFC)) +
            geom_histogram(colour = "white", fill = "black", binwidth = 0.15) +
            xlab(paste("Log2(fold-change)")) +
            theme(aspect.ratio = 1, text=element_text(size=18))

g = arrangeGrob(p1,p2 , ncol=2, main=paste("Reference:", reference))

plots.list = lappend(plots.list, g)

proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]

p_thr = 0.05
fc_thr = 1
top = 100
tmp.f = table(droplevels(filter(proteins.FC.f, p.value_BH < p_thr, (logFC < fc_thr | logFC > fc_thr))$ORF))
selected = names(sort(-tmp.f)[1:top])

proteins.FC.f.selected = droplevels(proteins.FC.f[proteins.FC.f$ORF %in% selected,])
proteins.FC.wide = droplevels(dcast(proteins.FC.f.selected, formula=ORF~KO, value.var="logFC"))


# cor_dist = as.dist(1 - abs(cor(proteins.FC.wide[,-1])))
# toClust = t(proteins.FC.wide[,-1])
# colnames(toClust) = proteins.FC.wide$ORF
# cl = NbClust(data = t(proteins.FC.wide[,-1]), diss=cor_dist, distance=NULL, method="ward.D2", index=c("ch"))


cl_res = ConsensusClusterPlus(as.matrix(proteins.FC.wide[,-1]), maxK=6,
                              reps=1000, pItem=0.8, pFeature=1,finalLinkage="complete",
                              clusterAlg="hc",distance="spearman",seed=123 )

# cl_res2 = ConsensusClusterPlus(as.matrix(t(proteins.FC.wide[,-1])), maxK=20,
#                               reps=1000, pItem=0.8, pFeature=1,title=title, finalLinkage="complete",
#                               clusterAlg="hc",distance="euclidean",seed=123 )



cl_V = 3
#cl_H = 2

annotation = data.frame(batch_date = exp_metadata$batch_date[match(colnames(proteins.FC.wide)[-1], exp_metadata$ORF)])
rownames(annotation) = colnames(proteins.FC.wide)[-1]
annotation$clusters = factor(cl_res[[cl_V]]$consensusClass)
#annotation$classes = droplevels(kinase_classes$class[match(KO_types$gene[match(colnames(proteins.FC.wide)[-1], KO_types$ORF)], kinase_classes$kinase)])

col_breaks = c(-4,-3,-2,-1,1,2,3,4)
#col_breaks = c(0)
file_name = paste("heatmap_proteins.p_thr.fc_thr.top", p_thr, fc_thr, top, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

aheatmap(proteins.FC.wide[,-1], annCol=annotation, 
         Colv=cl_res[[cl_V]]$consensusTree,
         color = paste("-RdBu", length(col_breaks)-1, sep=":"),
         breaks = col_breaks,
         filename=file_path)

proteins.FC.clusters = droplevels(merge(proteins.FC, data.frame(KO=rownames(annotation), cl = annotation$clusters), by="KO"))


## -- cluster enrichments ----

enrichments = list()
for (i in unique(proteins.FC.clusters$cl)) {
  for (j in c("up","down")) {
    tmp.selected = c()
    cluster = i
#     cluster = 1
#     j = "down"
    if (j == "up") {
      tmp.selected = table(droplevels(filter(proteins.FC.clusters, p.value_BH < p_thr, (logFC > fc_thr), cl == cluster)$ORF))  
    } else {
      tmp.selected = table(droplevels(filter(proteins.FC.clusters, p.value_BH < p_thr, (logFC < -fc_thr), cl == cluster)$ORF))  
    }
    
    important_size = round(max(table(droplevels(filter(proteins.FC.clusters, cl == cluster)$ORF)))/3)
    
    orf_thr = names(tmp.selected[tmp.selected >= important_size])
    orf_universe = unique(as.character(droplevels(filter(proteins.FC.clusters, cl == cluster)$ORF)))
    
    pathway.counts = pathway_enrichments(orf_thr=orf_thr, orf_universe=orf_universe, pathway2orf=GO_slim.process)
    pathway.counts$cl = factor(cluster)
    pathway.counts$direction = factor(j)
    enrichments[[i]] = pathway.counts
  }
}

enrichments.clusters = do.call(rbind.data.frame, enrichments)
View(enrichments.clusters[enrichments.clusters$p.value <0.05,])

## ---- protein variations ----
if (F) {
  
  proteins.long.ko_mean = proteins.long %>% 
    group_by(KO, ORF) %>%
    summarise(signal = mean(value))
  
  proteins.long.ko_mean = proteins.long.ko_mean[grep("none", proteins.long.ko_mean$KO, ignore.case=T, invert=T),]
  
  proteins.ko_mean = dcast(proteins.long.ko_mean, ORF~KO, value.var="signal")
  
  proteins.matrix.ko_mean = as.matrix(proteins.ko_mean[,-1])
  rownames(proteins.matrix.ko_mean) = proteins.ko_mean$ORF
  
  proteins.cv.variation = apply(proteins.matrix.ko_mean, 1, 
                                FUN=function(x) {
                                  sqrt(exp(sd(x)^2)-1)        
                                } )
  
  
  #proteins variation
  toPlot=data.frame(ORF=names(proteins.cv.variation), 
                    cv.variation = proteins.cv.variation) 
  p = ggplot(data=toPlot, aes(x=cv.variation)) +
    geom_histogram(binwidth=0.05, fill="white", colour="black") +
    xlab("Coefficient of variation of protein intensity across all KO") +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=14),
          aspect.ratio = 1)
  plots.list = lappend(plots.list, p)
  
  upper = 0.95
  lower = 0.05
  signal_upper = names(proteins.cv.variation[proteins.cv.variation > quantile(proteins.cv.variation,upper)])
  signal_lower = names(proteins.cv.variation[proteins.cv.variation < quantile(proteins.cv.variation,lower)])
  universe  = names(proteins.cv.variation)
  
  
  
  proteins.z_scores = apply(proteins.matrix.ko_mean, 1, 
                            FUN=function(x) {
                              unlog_x = exp(x)
                              sd.unlog_x = sd(unlog_x,na.rm=T)
                              mean.unlog_x = mean(unlog_x)
                              
                              sd.log = sqrt(log(1+sd.unlog_x^2/mean.unlog_x^2))
                              mean.log = log(mean.unlog_x)-0.5*sd.log^2
                              z.log = (x - mean.log)/sd.log
                              return(z.log)
                            } )
  
  
  
  #write.table(unique(colnames(proteins.z_scores)), file="proteins.txt", quote=F, row.names=F, col.names=F)
  mutants.z_sums = apply(proteins.z_scores, 1,
                         FUN=function(x) {
                           sum(abs(x))
                         })
  
  hist(mutants.z_sums, breaks=30, xlim=c(400,2000))
  points(mutants.z_sums[names(mutants.z_sums)=="WT"], 2, col="red", pch=20)
  points(mutants.z_sums[names(mutants.z_sums)=="YPL115C"], 2, col="red", pch=20)
  points(mutants.z_sums[names(mutants.z_sums)=="YPL185W"], 2, col="red", pch=20)
  
  text(mutants.z_sums[names(mutants.z_sums)=="WT"], 3, labels="WT")
}


## ---- comprarison with PPI network ----
if (F) {
  load("./R/objects/BIOGRID.ORGANISM.RData")
  
  BIOGRID.ORGANISM = tbl_df(BIOGRID.ORGANISM)
  BIOGRID.ORGANISM$Pubmed.ID = factor(BIOGRID.ORGANISM$Pubmed.ID)
  interactions = BIOGRID.ORGANISM %>% select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B,Experimental.System, Experimental.System.Type, Pubmed.ID)
  
  interactions.p = filter(interactions, Experimental.System.Type == "physical") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)
  interactions.g = filter(interactions, Experimental.System.Type == "genetic") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)
  G.phys = graph.data.frame(select(interactions.p, Systematic.Name.Interactor.A, Systematic.Name.Interactor.B), directed=F)
  G.gene = graph.data.frame(select(interactions.g, Systematic.Name.Interactor.A, Systematic.Name.Interactor.B), directed=F)
  
  proteins.FC.sums = proteins.FC.clusters %>% group_by(KO) %>% summarize(perturbation = sum(abs(logFC)))
  
  tmp = data.frame(KO = names(mutants.z_sums), 
                   perturb = proteins.FC.sums$perturbation[match(names(mutants.z_sums), proteins.FC.sums$KO)],
                   z_sums = mutants.z_sums, 
                   degree.phys = degree(G.phys)[match(names(mutants.z_sums), names(degree(G.phys)))],
                   degree.gene = degree(G.gene)[match(names(mutants.z_sums), names(degree(G.gene)))])
  
  
  proteins.z_scores.long = melt(proteins.z_scores, id.vars=row.names)
  names(proteins.z_scores.long) = c("KO", "ORF", "value")
  all_enrichments = ddply(proteins.z_scores.long, .(KO), 
                          .fun = function(z) {
                            #z = proteins.z_scores.long[proteins.z_scores.long$KO == "WT",]
                            signal = z$ORF[abs(z$value) >= 1.98]
                            universe =  z$ORF
                            
                            tmp.p = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.process)
                            tmp.p$type = "process"
                            tmp.f = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.function)
                            tmp.f$type = "function"
                            tmp.c = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.compartment)
                            tmp.c$type = "compartment"
                            tmp.path = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=pathway2orf)
                            tmp.path$type = "kegg_pathways"
                            
                            return(rbind(tmp.p, tmp.c, tmp.c, tmp.path))
                            
                          })
  
}


## -- comparison all vs all protein changes ----

all_KOs = unique(proteins.FC$KO)
kinases = as.character(unique(exp_metadata$ORF[exp_metadata$type == "Kinase"]))
all_KOs.kinases = all_KOs[all_KOs %in% kinases]


ind = 1:length(all_KOs.kinases)
tmp.ind  = t(combn(ind, 2))
tmp.names = matrix(all_KOs[as.vector(tmp.ind)], ncol=2, byrow=F) 
comparisons = tmp.names
#names(comparisons) = c("KO1", "KO2")

padj_thr = 0.05
int_abs = rep(0, nrow(comparisons))
int_up = rep(0, nrow(comparisons))
int_down = rep(0, nrow(comparisons))

union_abs = rep(0, nrow(comparisons))
union_up = rep(0, nrow(comparisons))
union_down = rep(0, nrow(comparisons))

for (i in 1:nrow(comparisons)) {
  pair = comparisons[i,]
  KO1 = pair[1]
  KO2 = pair[2]
  
  KO1_changed = filter(proteins.FC, KO == KO1, p.value_BH < padj_thr)
  KO2_changed = filter(proteins.FC, KO == KO2, p.value_BH < padj_thr)
  
  int_abs[i]   = length(intersect(KO1_changed$ORF, KO2_changed$ORF))
  int_up[i]    = length(intersect(KO1_changed$ORF[KO1_changed$logFC > 0 ], KO2_changed$ORF[KO2_changed$logFC > 0 ]))
  int_down[i]  = length(intersect(KO1_changed$ORF[KO1_changed$logFC < 0 ], KO2_changed$ORF[KO2_changed$logFC < 0 ]))
  union_abs[i] = length(union(KO1_changed$ORF, KO2_changed$ORF))
  union_up[i]  = length(union(KO1_changed$ORF[KO1_changed$logFC > 0 ], KO2_changed$ORF[KO2_changed$logFC > 0 ]))
  union_down[i] = length(union(KO1_changed$ORF[KO1_changed$logFC < 0 ], KO2_changed$ORF[KO2_changed$logFC < 0 ]))
}

comparisons.df = as.data.frame(cbind(comparisons, int_abs, int_up, int_down, union_abs, union_up, union_down))

#making matrix all vs all
#intersection
tmp = dcast(comparisons.df, formula=V1~V2, value.var="int_abs")
name.tmp = as.character(tmp$V1[1])
tmp.extraCol = cbind(NA,tmp[,-1])
colnames(tmp.extraCol)[1] = name.tmp
rownames(tmp.extraCol) = as.character(tmp$V1)

tmp.extraColRow = rbind(tmp.extraCol,NA)
rownames(tmp.extraColRow)[nrow(tmp.extraColRow)] = colnames(tmp.extraCol)[ncol(tmp.extraCol)] 

comparisons.matrix.tmp = as.matrix(tmp.extraColRow)
comparisons.matrix.tmp.t = t(comparisons.matrix.tmp)

comparisons.matrix.tmp[lower.tri(comparisons.matrix.tmp)] = comparisons.matrix.tmp.t[lower.tri(comparisons.matrix.tmp.t)]
comparisons.matrix.tmp[is.na(comparisons.matrix.tmp)] = 0 #diagonal

stopifnot(isSymmetric(comparisons.matrix.tmp))

comparisons.matrix = matrix(as.numeric(comparisons.matrix.tmp), ncol=ncol(comparisons.matrix.tmp), byrow=T)
colnames(comparisons.matrix) = colnames(comparisons.matrix.tmp)
rownames(comparisons.matrix) = rownames(comparisons.matrix.tmp)

# union matrix
# tmp = dcast(comparisons.df, formula=V1~V2, value.var="union_abs")
# name.tmp = as.character(tmp$V1[1])
# tmp.extraCol = cbind(NA,tmp[,-1])
# colnames(tmp.extraCol)[1] = name.tmp
# rownames(tmp.extraCol) = as.character(tmp$V1)
# 
# tmp.extraColRow = rbind(tmp.extraCol,NA)
# rownames(tmp.extraColRow)[nrow(tmp.extraColRow)] = colnames(tmp.extraCol)[ncol(tmp.extraCol)] 
# 
# union.matrix.tmp = as.matrix(tmp.extraColRow)
# union.matrix.tmp.t = t(union.matrix.tmp)
# 
# union.matrix.tmp[lower.tri(union.matrix.tmp)] = union.matrix.tmp.t[lower.tri(union.matrix.tmp.t)]
# union.matrix.tmp[is.na(union.matrix.tmp)] = 0 #diagonal
# 
# stopifnot(isSymmetric(union.matrix.tmp))
# 
# union.matrix = matrix(as.numeric(union.matrix.tmp), ncol=ncol(union.matrix.tmp), byrow=T)
# colnames(union.matrix) = colnames(union.matrix.tmp)
# rownames(union.matrix) = rownames(union.matrix.tmp)


#### heatmap of intersection
col_breaks = sort(c(seq(0,500,100),50))
toPlot = comparisons.matrix
annotation = data.frame(class = droplevels(kinase_classes$class[match(colnames(toPlot), kinase_classes$ORF)]),
                        batch_date = exp_metadata$batch_date[match(colnames(toPlot), exp_metadata$ORF)])
annotation$class = as.character(annotation$class)
annotation$class[is.na(annotation$class)] = "unclassified"
annotation$class = factor(annotation$class)


colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]
rownames(toPlot) = exp_metadata$gene[match(rownames(toPlot), exp_metadata$ORF)]

toPlot[is.infinite(toPlot)] = 0
toPlot[toPlot >= max(col_breaks)] = max(col_breaks)

annot_Rnames =  colnames(toPlot)
rownames(annotation) = annot_Rnames
#TODO: union set with all mutants
pheatmap(toPlot, breaks=col_breaks, annotation_col=annotation, color=brewer.pal(n=length(col_breaks)-1, name="Greys"))
    
proteins.FC$KO_gene = exp_metadata$gene[match(proteins.FC$KO, exp_metadata$ORF)]
proteins.FC$KO_type = exp_metadata$Type[match(proteins.FC$KO, exp_metadata$ORF)]
proteins.FC$KO_class = as.character(kinase_classes$class[match(proteins.FC$KO, kinase_classes$ORF)])
proteins.FC$KO_class[is.na(proteins.FC$KO_class)] = "unclassified"
proteins.FC$KO_class = factor(proteins.FC$KO_class)

## -- changes per kinase class ----
proteins.FC.classes.stats = proteins.FC %>% group_by(KO, KO_class) %>% summarise(changes = sum(p.value_BH <padj_thr))
p = ggplot(proteins.FC.classes.stats, aes(x=KO_class, y=changes)) +
      geom_boxplot()

tbl.pie = sort(table(proteins.FC.classes.stats$KO_class))
pie(tbl.pie, clockwise=F, init.angle=90, col=brewer.pal(length(tbl.pie),"Set3"))
p = recordPlot()
plots.list = lappend(plots.list, p)


# -- enrichments of changed genes ----- 

#orf2name = gene.annotations %>% distinct(V4, V6) %>% dplyr::select(V4, V6)

proteins.FC = tbl_df(proteins.FC)
proteins.FC$KO = factor(proteins.FC$KO)

all_enrichments = ddply(proteins.FC, .(KO), 
                        .fun = function(z) {
                          #z = proteins.z_scores.long[proteins.z_scores.long$KO == "WT",]
                          signal = z$ORF[abs(z$p.value_BH)< 0.05]
                          universe =  z$ORF              
                          tmp.p = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.process)
                          tmp.p$type = "process"
                          tmp.f = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.function)
                          tmp.f$type = "function"
                          tmp.c = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.compartment)
                          tmp.c$type = "compartment"
                          tmp.path = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=pathway2orf)
                          tmp.path$type = "kegg_pathways"
                        
                          return(rbind(tmp.p, tmp.c, tmp.c, tmp.path))
                          
                        })

enrich_thr = 0.05

#kegg
kegg.enrichments = dcast(droplevels(all_enrichments[all_enrichments$type == "kegg_pathways",]), pathway~KO, value.var="p.value")
kegg.enrichments.matrix = as.matrix(kegg.enrichments[,-1])
rownames(kegg.enrichments.matrix) = kegg.enrichments$pathway

kegg.enrichments.matrix[is.na(kegg.enrichments.matrix)] = 0
kegg.enrichments.matrix[kegg.enrichments.matrix > enrich_thr] = 0
kegg.enrichments.matrix[kegg.enrichments.matrix != 0] = 1
kegg.enrichments.matrix = kegg.enrichments.matrix[rowSums(kegg.enrichments.matrix) > 3,]

pathway2desription = pathway2orf %>% distinct(pathway, description) %>% dplyr::select(pathway, description)
pathway2desription$description = sub(pattern=" - Saccharomyces cerevisiae (budding yeast)", replacement="" , x=pathway2desription$description, fixed=T)


kegg.enrichments.matrix.kinases = kegg.enrichments.matrix[,colnames(kegg.enrichments.matrix) %in% kinases]
toPlot = kegg.enrichments.matrix.kinases
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]

pheatmap(toPlot, clustering_distance_cols="binary", cluster_rows=T, col=c("black","white"),annotation = annotation, legend=F,
         labels_row=pathway2desription$description[match(rownames(toPlot), pathway2desription$pathway)],)
         
p = recordPlot()
plots.list = lappend(plots.list, p)



#GO slim process
process.enrichments = dcast(droplevels(all_enrichments[all_enrichments$type == "process",]),pathway~KO, value.var="p.value", fun.aggregate=min)
process.enrichments.matrix = as.matrix(process.enrichments[,-1])
rownames(process.enrichments.matrix) = process.enrichments$pathway

process.enrichments.matrix[is.na(process.enrichments.matrix)] = 0
process.enrichments.matrix[process.enrichments.matrix > enrich_thr] = 0
process.enrichments.matrix[process.enrichments.matrix != 0] = 1
process.enrichments.matrix = process.enrichments.matrix[rowSums(process.enrichments.matrix) != 0,]

GOprocess2desription = GO_slim.process %>% distinct(pathway, description) %>% dplyr::select(pathway, description)

process.enrichments.matrix.kinases = process.enrichments.matrix[,colnames(process.enrichments.matrix) %in% kinases]
toPlot = process.enrichments.matrix.kinases

pheatmap(toPlot, clustering_distance_cols="binary", cluster_rows=F, 
         labels_row = GOprocess2desription$description[match(rownames(process.enrichments.matrix), GOprocess2desription$pathway)],col=c("black","white"),annotation = annotation, legend=F,
         labels_col = exp_metadata$gene[match(colnames(process.enrichments.matrix), exp_metadata$ORF)])
p = recordPlot()
plots.list = lappend(plots.list, p)


#PCA by pathways
proteins.FC.kegg = tbl_df(merge(proteins.FC, pathway2orf, by="ORF"))
proteins.FC.kegg.fractions = proteins.FC.kegg %>% group_by(KO, pathway) %>% summarise(pathway.fraction = sum(p.value_BH < 0.05)/length(p.value_BH))

kegg.fractions = dcast(proteins.FC.kegg.fractions, formula=pathway~KO, value.var="pathway.fraction")
kegg.fractions.matrix = as.matrix(kegg.fractions[,-1])
rownames(kegg.fractions.matrix) = kegg.fractions$pathway
kegg.fractions.matrix = kegg.fractions.matrix[rowSums(kegg.fractions.matrix) != 0,]


toPlot = kegg.fractions.matrix[,colnames(kegg.fractions.matrix) %in% kinases]

rownames(toPlot) = pathway2desription$description[match(rownames(toPlot), pathway2desription$pathway)]
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]
S = prcomp(t(toPlot))
biplot(S, cex=0.66)
p = recordPlot()
plots.list = lappend(plots.list, p)


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 
