rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis1"

## -- Load of data ----
load("./R/objects/proteins.matrix.combat.quant.FC.RData")
load("./R/objects/proteins.matrix.combat.quant.RData")


load("./R/objects/pathway2orf._load_.RData")
load("./R/objects/pathways._load_.RData")
load("./R/objects/GO_slim.function.RData")
load("./R/objects/GO_slim.process.RData")
load("./R/objects/GO_slim.compartment.RData")

load("./R/objects/kinase_classes._clean_.RData")
load("./R/objects/exp_metadata._clean_.RData")

library(cowplot)
library(scales)
load("./R/objects/pathway2orf._load_.RData")
load("./R/objects/pathways._load_.RData")
load("./R/objects/gene.annotations._load_.RData")
load("./R/objects/modules.RData")
load("./R/objects/kegg_categories._load_.RData")


load("./R/objects/gene.annotations._load_.RData")
orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]

load("./R/objects/yeastract._load_.RData")
load("./R/objects/GO_slim.raw._load_.RData")
load("./R/objects/STRING._load_.RData")
load("./R/objects/yeast.ppi._load_.RData")


getFC_thr = function(proteins.matrix = proteins.matrix.combat, pval_thr = 0.01) {
  
  #proteins.matrix = proteins.matrix.combat
  #pval_thr = 0.05
  
  #checking WT samples to define FC
  
  exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
  cl = pam(exp_metadata$aquisition_date.str, 7)
  exp_metadata$batch_kmeans = cl$clustering
    
  pheno_wt = as.data.frame(exp_metadata[match(colnames(proteins.matrix), exp_metadata$sample_name),])
  pheno_wt = pheno_wt[pheno_wt$type == "Standard Mix",]
  pheno_wt = pheno_wt[pheno_wt$batch_kmeans %in% names(table(pheno_wt$batch_kmeans))[table(pheno_wt$batch_kmeans)  >= 3],]
  
  #plot(exp_metadata$aquisition_date, exp_metadata$batch_kmeans)
  
  #points(pheno_wt$aquisition_date, pheno_wt$batch_kmeans, col="red")
  
  proteins.matrix.WT = proteins.matrix[,match(pheno_wt$sample_name, colnames(proteins.matrix))]
  
#   s = prcomp(t(proteins.matrix.WT))
#   plot(s$x[,c(1,2)], col=pheno_wt$batch_kmeans)
  
  pheno_wt$factor = factor(paste(pheno_wt$ORF, pheno_wt$batch_kmeans, sep="."))
  
  X = model.matrix(~factor + 0, data=pheno_wt)
  colnames(X) = levels(pheno_wt$factor)
  
  tbl.tmp = table(pheno_wt$factor)
  reference = names(tbl.tmp)[tbl.tmp == max(tbl.tmp)][1]
  
  matrix = proteins.matrix.WT
  
  lm.fit_model = lmFit(matrix, X)
  ph = unique(as.character(pheno_wt$factor))
  contrasts = paste0( ph[ph !=reference] ,"-", reference)  
  
  mc = makeContrasts(contrasts=contrasts, levels=X)    
  c.fit = contrasts.fit(lm.fit_model, mc)
  eb = eBayes(c.fit)
  
  folds = rowFolds(data=exp(matrix), groups=pheno_wt$factor, reference=reference)
  folds = log(folds, 2)
  
  folds_tmp = melt(eb$coefficients, id.vars="row.names")
  #folds_tmp$contrasts = factor(paste(folds_tmp$contrasts, "-", reference, sep=""))
  pvals_tmp = melt(eb$p.value, id.vars="row.names")
    
  names(folds_tmp) = c("ORF", "contrasts", "logFC")
  names(pvals_tmp) = c("ORF", "contrasts", "p.value")
  
  folds_tmp$contrasts = factor(folds_tmp$contrasts)
  pvals_tmp$contrasts = factor(pvals_tmp$contrasts)
    
  proteins.FC = merge(folds_tmp, pvals_tmp, all=T,
                      by=c("ORF", "contrasts"))
  
  ##multiple testing correction
  proteins.FC$p.value_BH = p.adjust(proteins.FC$p.value, method="BH")
  proteins.FC$p.value_bonferroni = p.adjust(proteins.FC$p.value, method="bonferroni")
  
  proteins.FC$KO = sub(x = proteins.FC$contrasts, pattern=paste("(.*?)-", reference, sep=""), replacement="\\1")
  proteins.FC$reference = reference
  
  data = abs(proteins.FC$logFC[proteins.FC$p.value_BH < pval_thr])
  
  
  

  file_name = paste(fun_name, "getFC_thr", "pdf", sep=".")
  file_path = paste(figures_dir, file_name, sep="/")
  
  pdf(file_path, paper="a4")
  par(pty="s")
  hist(data, breaks=50, main="Expected fold changes in Standatd Mix")
  fc_thr = max(data)

  abline(v=fc_thr, lty=2)
  legend("topleft", bg=NULL, bty="n", 
         legend=paste("Median FC=", round(fc_thr,2)))
  dev.off()
  
  return(abs(fc_thr))
}


#### VAR ####
protein.matrix = proteins.matrix.combat.quant
proteins.FC = proteins.matrix.combat.quant.FC

## -- coverage: basic stats ----
EC.genes = gene.annotations[gene.annotations$V3 == "EC number",]

KEGG.pathways = distinct(pathway2orf, pathway, ORF) 
KEGG.pathways$EC.number = EC.genes$V1[match(KEGG.pathways$ORF, EC.genes$V4)]
KEGG.pathways = droplevels(KEGG.pathways %>% filter( !is.na(EC.number)))
KEGG.pathways = KEGG.pathways %>% group_by(ORF) %>% mutate(n = n())

  
measured.ORFs = data.frame(ORF = unique(rownames(protein.matrix)))
measured.ORFs = droplevels(merge(measured.ORFs, EC.genes, by.x="ORF", by.y = "V4")[,c(1,2)])
names(measured.ORFs) = c("ORF", "EC.number")

KEGG.pathways.stats = tbl_df(KEGG.pathways) %>% group_by(pathway, description) %>% 
  summarize(inData.ORF = sum(unique(ORF) %in% measured.ORFs[,1]),
            inData.EC = sum(unique(EC.number) %in% measured.ORFs[,2]),
            nORF = length(unique(ORF)),
            nEC = length(unique(EC.number)),
            EC.coverage = inData.EC/nEC,
            ORF.coverage = inData.ORF/nORF)


KEGG.pathways.stats = merge(KEGG.pathways.stats, kegg_categories, by="pathway")
KEGG.pathways = merge(KEGG.pathways, kegg_categories, by="pathway")




selected = c("Amino acid metabolism",                       
             "Carbohydrate metabolism",                                          
             "Energy metabolism",                                      
             "Glycan biosynthesis and metabolism",                                    
             "Metabolism of cofactors and vitamins",       
             "Metabolism of other amino acids",                
             "Nucleotide metabolism",                                                         
             "Metabolism of terpenoids and polyketides",
             "Biosynthesis of other secondary metabolites",
             "Lipid metabolism")

KEGG.pathways.f = droplevels(KEGG.pathways[KEGG.pathways$B %in% selected,])

KEGG.pathways.stats.f = droplevels(KEGG.pathways.stats[KEGG.pathways.stats$B %in% selected,])
KEGG.pathways.stats.f = KEGG.pathways.stats.f %>% group_by(B) %>% distinct(pathway) %>% mutate(avgEC = mean(EC.coverage, na.rm=T)) %>%  arrange(EC.coverage) 
KEGG.pathways.stats.f$pathway = factor(KEGG.pathways.stats.f$pathway, levels = as.character(KEGG.pathways.stats.f$pathway))



p = ggplot(KEGG.pathways.stats.f, aes(x=pathway, y=EC.coverage, fill=B)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels = KEGG.pathways.stats.f$C ) +
  coord_flip() + theme_light()
plots.list = lappend(plots.list, p)

KEGG.pathways.stats.f = KEGG.pathways.stats.f %>% group_by(B) %>% distinct(pathway) %>% mutate(avgORF = mean(ORF.coverage, na.rm=T)) %>% arrange(ORF.coverage) 
KEGG.pathways.stats.f$pathway = factor(KEGG.pathways.stats.f$pathway, levels = as.character(KEGG.pathways.stats.f$pathway))


p = ggplot(KEGG.pathways.stats.f, aes(x=pathway, y=ORF.coverage, fill=B)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels = KEGG.pathways.stats.f$C ) +
  coord_flip() + theme_light()
plots.list = lappend(plots.list, p)


toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Carbohydrate metabolism"))
p1 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip() + 
  theme(aspect.ratio = 1)

toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Amino acid metabolism"))
p2 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip() + 
  theme(aspect.ratio = 1)

toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Energy metabolism"))
p3 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip() + 
  theme(aspect.ratio = 1)

toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Lipid metabolism"))
p4 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip() + 
  theme(aspect.ratio = 1)

# -- coverage: 4 together ----
select.tmp = data.frame(cat = c("Carbohydrate metabolism", "Amino acid metabolism", "Lipid metabolism", "Energy metabolism"),
                        order = c(1,2,3,4))

toPlot = droplevels(KEGG.pathways.stats.f[KEGG.pathways.stats.f$B %in% select.tmp$cat ,])
toPlot$B = factor(toPlot$B, levels = select.tmp$cat)
toPlot$order = select.tmp$order[match(toPlot$B, select.tmp$cat)]

toPlot =  toPlot[with(toPlot, order(-order, EC.coverage)),]
toPlot$pathway = factor(toPlot$pathway, levels = as.character(toPlot$pathway))

p = ggplot(toPlot, aes(x=pathway, y=EC.coverage, fill=B)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip() + theme_light()

plots.list = lappend(plots.list, p)




## -- vulcano plot ----

reference = unique(as.character(proteins.FC$reference))

lb = -4
ub = 4
toPlot = proteins.FC
toPlot[toPlot$logFC < 0,]$logFC = ifelse(toPlot[toPlot$logFC < 0,]$logFC < lb, lb, toPlot[toPlot$logFC < 0,]$logFC)
toPlot[toPlot$logFC > 0,]$logFC = ifelse(toPlot[toPlot$logFC > 0,]$logFC > ub, ub, toPlot[toPlot$logFC > 0,]$logFC)


pval_thr = 0.01
set.seed(123)
FC_thr = getFC_thr(proteins.matrix=proteins.matrix.combat.quant, pval_thr=pval_thr)

proteins.FC.stats = data.frame( ratio_sign = round(sum(proteins.FC$logFC < -FC_thr & proteins.FC$p.value_BH < pval_thr)/sum(proteins.FC$logFC > FC_thr & proteins.FC$p.value_BH < pval_thr),2),
                                ratio =  round(sum(proteins.FC$logFC < 0)/sum(proteins.FC$logFC > 0),2),
                                n_prot = length(unique(proteins.FC$ORF)),
                                n_sign = round(sum(proteins.FC$p.value_BH<pval_thr & abs(proteins.FC$logFC)>FC_thr )/length(unique(proteins.FC$contrasts)),2),
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


#file_name = paste("volcano_proteins","png", sep=".")
#file_path = paste(figures_dir, file_name, sep="/")
#png(file_path, res=150, width=11.69+0.1*11.69, height=8.27+0.1*8.27, units="in")
#dev.off()

## -- boxplots of kinases mutants ----
proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]
proteins.FC.f$isMetabolic = proteins.FC.f$ORF %in% unique(KEGG.pathways.f$ORF)
proteins.FC.f$isEnzyme = proteins.FC.f$ORF %in% unique(EC.genes$V4)

toPlot = proteins.FC.f
toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)

toPlot$label = factor(orf2name$gene_name[match(toPlot$KO, orf2name$ORF)])


p3 = ggplot(toPlot, aes(x=KO, y = logFC)) + geom_boxplot() + 
          geom_point(data=filter(toPlot,sign == 1), color="cyan", alpha=0.5) +
          scale_x_discrete(labels=toPlot$label)+
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
plots.list = lappend(plots.list, p3)

toPlot = filter(toPlot, isMetabolic==T)
p4 = ggplot(toPlot, aes(x=KO, y = logFC)) + geom_boxplot() + 
            geom_point(data=filter(toPlot,sign == 1), color="cyan", alpha=0.5) +
            scale_x_discrete(labels=toPlot$label)+
            ggtitle("Metabolic")
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
plots.list = lappend(plots.list, p4)




pathways.FC = merge(KEGG.pathways.f, proteins.FC.f, by="ORF")
toPlot = pathways.FC
p = ggplot(toPlot, aes(x=C, y=logFC, fill=B)) + 
                      geom_boxplot()
                      


if(F) {
  p_thr = 0.05
  fc_thr = FC_thr
  top = 1000
  
  tmp.f = table(droplevels(filter(proteins.FC.f, p.value_BH < p_thr, (logFC < fc_thr | logFC > fc_thr))$ORF))
  
  tmp.f = tmp.f[tmp.f > round(length(unique(proteins.FC.f$contrasts))/3)]
  
  selected = as.character(na.omit(names(sort(-tmp.f)[1:top])))
  proteins.FC.f.selected = droplevels(filter(proteins.FC.f[proteins.FC.f$ORF %in% selected,], isMetabolic == T))
  proteins.FC.wide = droplevels(dcast(proteins.FC.f.selected, formula=ORF~KO, value.var="logFC"))
  
  
  # cor_dist = as.dist(1 - abs(cor(proteins.FC.wide[,-1])))
  # toClust = t(proteins.FC.wide[,-1])
  # colnames(toClust) = proteins.FC.wide$ORF
  # cl = NbClust(data = t(proteins.FC.wide[,-1]), diss=cor_dist, distance=NULL, method="ward.D2", index=c("ch"))
  
  
  cl_res = ConsensusClusterPlus(as.matrix(proteins.FC.wide[,-1]), maxK=10,
                                reps=100, pItem=0.8, pFeature=1,finalLinkage="complete",
                                clusterAlg="hc",distance="spearman",seed=123 )
  
  #cl_H = 2
  
  annotation = data.frame(batch_date = exp_metadata$batch_date[match(colnames(proteins.FC.wide)[-1], exp_metadata$ORF)])
  rownames(annotation) = colnames(proteins.FC.wide)[-1]
  annotation$clusters = factor(cl_res[[cl_V]]$consensusClass)
  #annotation$classes = droplevels(kinase_classes$class[match(KO_types$gene[match(colnames(proteins.FC.wide)[-1], KO_types$ORF)], kinase_classes$kinase)])
  
  col_breaks = c(-2,-1,-0.5,0.5,1,2)
  #col_breaks = c(0)
  file_name = paste("heatmap_proteins.p_thr.fc_thr.top", p_thr, fc_thr, top, "pdf", sep=".")
  file_path = paste(figures_dir, file_name, sep="/")
  
  aheatmap(proteins.FC.wide[,-1], annCol=annotation, 
           Colv=cl_res[[cl_V]]$consensusTree,
           color = paste("-RdBu", length(col_breaks)-1, sep=":"),
           breaks = col_breaks)
  
  pheatmap(proteins.FC.wide[,-1])
  
  proteins.FC.clusters = droplevels(merge(proteins.FC.f, data.frame(KO=rownames(annotation), cl = annotation$clusters), by="KO"))
  
  
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
  
}

#View(enrichments.clusters[enrichments.clusters$p.value <0.05,])

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
  interactions = BIOGRID.ORGANISM %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B,Experimental.System, Experimental.System.Type, Pubmed.ID)
  
  interactions.p = filter(interactions, Experimental.System.Type == "physical") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)
  interactions.g = filter(interactions, Experimental.System.Type == "genetic") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)
  G.phys = graph.data.frame(dplyr::select(interactions.p, Systematic.Name.Interactor.A, Systematic.Name.Interactor.B), directed=F)
  G.gene = graph.data.frame(dplyr::select(interactions.g, Systematic.Name.Interactor.A, Systematic.Name.Interactor.B), directed=F)
  
  
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
tmp.names = matrix(all_KOs.kinases[as.vector(tmp.ind)], ncol=2, byrow=F) 
comparisons = tmp.names

# -- PPI similarities with intersection --

load("./R/objects/yeast.ppi._load_.RData")
yeast.ppi = tbl_df(yeast.ppi)

yeast.ppi.p = droplevels(yeast.ppi %>% filter(Author == "Breitkreutz A (2010)", Experimental.System.Type == "physical") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B))
yeast.ppi.g = droplevels(yeast.ppi %>% filter(Author == "Costanzo M (2010)", Experimental.System.Type == "genetic") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B))



G.phys = graph.data.frame(droplevels(yeast.ppi.p %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)
G.gene = graph.data.frame(droplevels(yeast.ppi.g %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)


#names(comparisons) = c("KO1", "KO2")

padj_thr = pval_thr

int_abs = rep(0, nrow(comparisons))
int_up = rep(0, nrow(comparisons))
int_down = rep(0, nrow(comparisons))

union_abs = rep(0, nrow(comparisons))
union_up = rep(0, nrow(comparisons))
union_down = rep(0, nrow(comparisons))


int_ppi = rep(0, nrow(comparisons))
for (i in 1:nrow(comparisons)) {
  pair = comparisons[i,]
  KO1 = pair[1]
  KO2 = pair[2]
  
  KO1_changed = filter(proteins.FC.f, isMetabolic == T, KO == KO1, p.value_BH < padj_thr)
  KO2_changed = filter(proteins.FC.f, KO == KO2, p.value_BH < padj_thr)
  
  int_abs[i]   = length(intersect(KO1_changed$ORF[abs(KO1_changed$logFC) > FC_thr ], KO2_changed$ORF[abs(KO2_changed$logFC) > FC_thr ]))
  int_up[i]    = length(intersect(KO1_changed$ORF[KO1_changed$logFC > FC_thr ], KO2_changed$ORF[KO2_changed$logFC > FC_thr ]))
  int_down[i]  = length(intersect(KO1_changed$ORF[KO1_changed$logFC < -FC_thr ], KO2_changed$ORF[KO2_changed$logFC < -FC_thr ]))
  union_abs[i] = length(union(KO1_changed$ORF, KO2_changed$ORF))
  union_up[i]  = length(union(KO1_changed$ORF[KO1_changed$logFC > FC_thr ], KO2_changed$ORF[KO2_changed$logFC > FC_thr ]))
  union_down[i] = length(union(KO1_changed$ORF[KO1_changed$logFC < -FC_thr ], KO2_changed$ORF[KO2_changed$logFC < -FC_thr ]))
  
#    if (KO1 %in% V(GRAPH)$name & KO2 %in% V(GRAPH)$name) {
#      int_ppi[i] = length(intersect(V(GRAPH)$name[neighbors(GRAPH,v=KO1)], V(GRAPH)$name[neighbors(GRAPH,v=KO2)]))
#    } else {
#      int_ppi[i] = NA  
#    } 
}

comparisons.df = as.data.frame(cbind(comparisons, int_abs, int_up, int_down, union_abs, union_up, union_down, int_ppi))

#making matrix all vs all

#intersection
tmp = dcast(comparisons.df, formula=V1~V2, value.var="int_abs")
#tmp = dcast(comparisons.df, formula=V1~V2, value.var="int_ppi")

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


#ppi

tmp = dcast(comparisons.df, formula=V1~V2, value.var="int_ppi")

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

comparisons.matrix.ppi = matrix(as.numeric(comparisons.matrix.tmp), ncol=ncol(comparisons.matrix.tmp), byrow=T)
colnames(comparisons.matrix.ppi) = colnames(comparisons.matrix.tmp)
rownames(comparisons.matrix.ppi) = rownames(comparisons.matrix.tmp)

library(ade4)
mantel.rtest(as.dist(comparisons.matrix.ppi), as.dist(comparisons.matrix))


cl_res = ConsensusClusterPlus(comparisons.matrix, maxK=6,
                              reps=1000, pItem=0.9, pFeature=1,
                              clusterAlg="hc",seed=123 )






tmp = melt(comparisons.matrix, formula=rownames(comparisons.matrix)~colnames(comparisons.matrix))

hcl = hclust(dist(comparisons.matrix))
tmp$X1 = factor(tmp$X1, levels = hcl$labels[hcl$order])
tmp$X2 = factor(tmp$X2, levels = hcl$labels[hcl$order])

ggplot(tmp,aes(x=X1, y=X2, fill=value)) + geom_tile()

##-- heatmap of intersection ----

col_breaks = sort(c(seq(0,100,10),5))
toPlot = comparisons.matrix

cl_V = 3



annotation = data.frame(class = droplevels(kinase_classes$class[match(colnames(toPlot), kinase_classes$ORF)]),
                        batch_date = exp_metadata$batch_date[match(colnames(toPlot), exp_metadata$ORF)],
                        clusters = factor(cl_res[[cl_V]]$consensusClass))
annotation$class = as.character(annotation$class)
annotation$class[is.na(annotation$class)] = "unclassified"
annotation$class = factor(annotation$class)


colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]
rownames(toPlot) = exp_metadata$gene[match(rownames(toPlot), exp_metadata$ORF)]

#toPlot[is.infinite(toPlot)] = 0
#toPlot[toPlot >= max(col_breaks)] = max(col_breaks)

annot_Rnames =  colnames(toPlot)
rownames(annotation) = annot_Rnames
#TODO: union set with all mutants
pheatmap(toPlot, breaks=col_breaks, annotation_col=annotation, color=brewer.pal(n=length(col_breaks)-1, name="Greys"))

parts = 3
S = pheatmap(toPlot,  annotation_col=annotation, cutree_cols=parts, colorRampPalette(c("white", "black"))(10))



#cutting to 3 parts
KO.stats = data.frame(KO = droplevels(exp_metadata$ORF[match(names(cutree(S$tree_col,parts)), exp_metadata$gene)]))
KO.stats$cluster = factor(cutree(S$tree_col,parts))

KO.stats$degree.p = degree(G.phys)[match(KO.stats$KO, names(degree(G.phys)))]
KO.stats$degree.g = degree(G.gene)[match(KO.stats$KO, names(degree(G.gene)))]

annotation$cl = KO.stats$cluster
annotation$clusters = NULL

pheatmap(toPlot,  annotation_col=annotation, 
             cutree_cols=parts)

p = recordPlot()
plots.list = lappend(plots.list, p)



p = ggplot(melt(KO.stats, id.vars=c("KO", "cluster")), aes(fill=variable, x=cluster, y=value)) +
          geom_boxplot()
plots.list = lappend(plots.list, p)


proteins.FC$KO_gene = exp_metadata$gene[match(proteins.FC$KO, exp_metadata$ORF)]
proteins.FC$KO_type = exp_metadata$Type[match(proteins.FC$KO, exp_metadata$ORF)]
proteins.FC$KO_class = as.character(kinase_classes$class[match(proteins.FC$KO, kinase_classes$ORF)])
proteins.FC$KO_class[is.na(proteins.FC$KO_class)] = "unclassified"
proteins.FC$KO_class = factor(proteins.FC$KO_class)

## -- changes per kinase class ----
proteins.FC.classes.stats = proteins.FC %>% group_by(KO, KO_class) %>% summarise(changes = sum(p.value_BH <padj_thr))
p = ggplot(proteins.FC.classes.stats, aes(x=KO_class, y=changes)) +
      geom_boxplot()

plots.list = lappend(plots.list, p)

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

pheatmap(toPlot, clustering_distance_cols="binary", clustering_distance_rows="binary", 
         cluster_rows=T, col=c("lightgrey","black"),annotation = annotation, legend=F,
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
process.enrichments.matrix = process.enrichments.matrix[rowSums(process.enrichments.matrix) >3,]

GOprocess2desription = GO_slim.process %>% distinct(pathway, description) %>% dplyr::select(pathway, description)

process.enrichments.matrix.kinases = process.enrichments.matrix[,colnames(process.enrichments.matrix) %in% kinases]
toPlot = process.enrichments.matrix.kinases
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]

pheatmap(toPlot, clustering_distance_cols="binary",  clustering_distance_rows="binary", cluster_rows=T,
         labels_row = GOprocess2desription$description[match(rownames(process.enrichments.matrix), GOprocess2desription$pathway)],col=c("lightgrey","black"),annotation = annotation)
       
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
s1 = prcomp(t(toPlot))
biplot(s1, cex=0.66)
abline(h=0,v=0)
l1 = sort(s1$rotation[,1],decreasing=T)
l2 = sort(S1$rotation[,2], decreasing=T)
p = recordPlot()
biplot.default

xPC = 1
yPC = 2
idx = match(unique(names(sort(abs(c(s1$rotation[,1], s1$rotation[,2])),decreasing=T)[1:20])), rownames(s1$rotation))
loads = s1$rotation[idx,c(xPC,yPC)]

biplot(x=s1$x[,c(xPC,yPC)],y=loads, cex=0.66, ylim = c(-3,2),expand=1.5,
     xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
     ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
abline(h=0,v=0)
p = recordPlot()
plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")

## -- perturbation vs shortest pathways ---- 

# mapkpathway <- parseKGML(file="./data/2015-05-29/ko04011.kgml")
# 
# mapkG2 <- KEGGpathway2Graph(mapkpathway, genesOnly=T, expandGenes=T)
# 
# G = igraph.from.graphNEL(mapkG2, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
# 
# V(G)$name = sub(pattern="sce:(.*?)", replacement="\\1", x=V(G)$name)
# V(G)$gene_name = orf2name$gene_name[match(V(G)$name, orf2name$ORF)]
# V(G)$ORF = V(G)$name
# V(G)$name = V(G)$gene_name
# 
# #first nucleous proteins
# nucleous = c("DIG1","DIG2","MCM1", "STE12", "FAR1", "RLM1", "MSN2", "MSN4", "SWI4", "SWI6")
# paths = shortest.paths(G, to=nucleous, mode="out" )
# degrees = degree(G, mode="all")
#
# paths[is.infinite(paths)] = 100
# tmp = apply(paths, 1, min) #minimal 
# 
# proteins.FC.f.stats.inG$distance = tmp[match(proteins.FC.f.stats.inG$gene_name , names(tmp))]
# proteins.FC.f.stats.inG$degree = degrees[match(proteins.FC.f.stats.inG$gene_name , names(degrees))]

yeast.ppi = tbl_df(yeast.ppi)
yeast.ppi.p.all = droplevels(yeast.ppi %>% filter(Experimental.System.Type == "physical") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B))
yeast.ppi.g.all = droplevels(yeast.ppi %>% filter(Experimental.System.Type == "genetic") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B))

yeast.ppi.p = droplevels(yeast.ppi %>% filter(Author == "Breitkreutz A (2010)", Experimental.System.Type == "physical") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B))
yeast.ppi.g = droplevels(yeast.ppi %>% filter(Author == "Costanzo M (2010)", Experimental.System.Type == "genetic") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B))

string.exp = STRING %>% filter(experimental > 900)
string.all = STRING %>% filter(combined_score > 900)


G.string.exp = graph.data.frame(droplevels(string.exp %>% dplyr::select(ORF1, ORF2)), directed=F)
G.string.all = graph.data.frame(droplevels(string.all %>% dplyr::select(ORF1, ORF2)), directed=F)

G.phys = graph.data.frame(droplevels(yeast.ppi.p %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)
G.gene = graph.data.frame(droplevels(yeast.ppi.g %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)
#G.phys.all  = graph.data.frame(droplevels(yeast.ppi.p.all %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)
#G.gene.all  = graph.data.frame(droplevels(yeast.ppi.g.all %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)


proteins.FC.f = droplevels(proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),])
proteins.FC.f$sign = ifelse(abs(proteins.FC.f$logFC) >= FC_thr & proteins.FC.f$p.value_BH < pval_thr, 1,0)
proteins.FC.f$isMetabolic = proteins.FC.f$ORF %in% unique(KEGG.pathways.f$ORF)
proteins.FC.f$isEnzyme = proteins.FC.f$ORF %in% unique(EC.genes$V4)

proteins.FC.f.stats = proteins.FC.f %>% filter(isMetabolic == T) %>% group_by(KO) %>% dplyr::summarise(changes = sum(sign==1),
                                                                          perturbation = sum(abs(logFC[p.value_BH < pval_thr])))



proteins.FC.f.stats$KO %in% unique(as.vector(GO_slim.raw$V1[grep(x=GO_slim.raw$V5, pattern="transcription factor activity")]))
proteins.FC.f.stats$KO %in% unique(yeastract)

yeastract$TF_ORF = orf2name$ORF[match(yeastract$TF, orf2name$gene_name)]

yeastract.stats = yeastract %>% group_by(TF, TF_ORF) %>% summarise(n = n())

TFs = unique(as.vector(GO_slim.raw$V1[grep(x=GO_slim.raw$V5, pattern="nucleic acid binding transcription factor activity")]))

GRAPH = G.string.exp
phys.TF_ORFs = TFs[TFs %in% V(GRAPH)$name]

kinases = unique(proteins.FC.f.stats$KO)[unique(proteins.FC.f.stats$KO) %in% V(GRAPH)$name]
paths = shortest.paths(GRAPH, v=kinases ,  to=phys.TF_ORFs)
paths.all = get.all.shortest.paths(GRAPH, from=kinases, to=phys.TF_ORFs)

min.paths = apply(paths, 1, min)
min.paths[is.infinite(min.paths)] = NA

paths.long = melt(paths) 
names(paths.long) = c("KO", "TF","value")
paths.long.stats = paths.long %>% group_by(KO) %>% summarise( mean.min = mean(yeastract.stats$n[match(TF[value %in% min(value)], yeastract.stats$TF_ORF)],na.rm=T),
                                                              sum.min = sum(yeastract.stats$n[match(TF[value %in% min(value)], yeastract.stats$TF_ORF)],na.rm=T),
                                                              n = length(yeastract.stats$n[match(TF[value %in% min(value)], yeastract.stats$TF_ORF)]))


proteins.FC.f.stats$gene_name = orf2name$gene_name[match(proteins.FC.f.stats$KO, orf2name$ORF)]
proteins.FC.f.stats$phys.degree = degree(G.phys)[match(proteins.FC.f.stats$KO, names(degree(G.phys)))]
proteins.FC.f.stats$betweenness = betweenness(GRAPH)[match(proteins.FC.f.stats$KO, names(degree(GRAPH)))]
proteins.FC.f.stats$degree = degree(GRAPH)[match(proteins.FC.f.stats$KO, names(degree(GRAPH)))]
#proteins.FC.f.stats$gene.all.degree = degree(G.gene.all)[match(proteins.FC.f.stats$KO, names(degree(G.gene.all)))]
proteins.FC.f.stats$gene.degree = degree(G.gene)[match(proteins.FC.f.stats$KO, names(degree(G.gene)))]
proteins.FC.f.stats$min.paths = min.paths[match(proteins.FC.f.stats$KO, names(min.paths))]
proteins.FC.f.stats$min.paths.comb = ifelse(min.paths[match(proteins.FC.f.stats$KO, names(min.paths))] <= 1, 0, 1)

proteins.FC.f.stats = merge(proteins.FC.f.stats, paths.long.stats, by="KO", all=T)
proteins.FC.f.stats$min.paths.comb[is.na(proteins.FC.f.stats$min.paths.comb)] = "none"
proteins.FC.f.stats$min.paths[is.na(proteins.FC.f.stats$min.paths)] = "none"

proteins.FC.f.stats.long = melt(as.data.frame(proteins.FC.f.stats), id.vars=c("KO", "gene_name", "min.paths", "min.paths.comb")) 


control.min.paths.comb = proteins.FC.f.stats.long %>% filter(min.paths.comb == 0)
control.min.paths = proteins.FC.f.stats.long %>% filter(min.paths == 1)


min.paths.comb.stats = proteins.FC.f.stats.long %>% filter(min.paths.comb != 0) %>% 
                             group_by(min.paths.comb, variable) %>% 
                             summarize( FC.mean = mean(value, na.rm=T)/mean(control.min.paths.comb[control.min.paths.comb[,5] == variable, "value"], na.rm=T),
                                        FC.median = median(value, na.rm=T)/median(control.min.paths.comb[control.min.paths.comb[,5] == variable, "value"], na.rm=T),
                                        p.value = wilcox.test(value,control.min.paths.comb[control.min.paths.comb[,5] == variable, "value"])$p.value)

min.paths.stats = proteins.FC.f.stats.long %>% filter(min.paths != 1) %>% 
                             group_by(min.paths, variable) %>% 
                             summarize( FC.mean = mean(value, na.rm=T)/mean(control.min.paths[control.min.paths[,5] == variable, "value"], na.rm=T),
                                        FC.median = median(value, na.rm=T)/median(control.min.paths[control.min.paths[,5] == variable, "value"], na.rm=T),
                                        p.value = wilcox.test(value, control.min.paths[control.min.paths[,5] == variable, "value"])$p.value)


min.paths.stats$symbol = ""
min.paths.stats$symbol = ifelse(min.paths.stats$p.value<0.05, "*", "")
#min.paths.stats$symbol = ifelse(min.paths.stats$p.value<0.001, "**", "")


min.paths.comb.stats$symbol = ""
min.paths.comb.stats$symbol = ifelse(min.paths.comb.stats$p.value<0.05, "*", "")
#min.paths.comb.stats$symbol = ifelse(min.paths.comb.stats$p.value<0.001, "**", "")


p1 = ggplot(proteins.FC.f.stats.long, aes(x=factor(min.paths.comb), y=value)) +
  geom_boxplot() +
  geom_text(data=min.paths.comb.stats, aes(x = factor(min.paths.comb) , y=1, label=symbol), col="red", size=10)+
  scale_x_discrete(labels=c("First", ">1", "None")) +
  xlab("Shortest distance to transcription factor") +
  facet_wrap(~variable, scales="free")        

p2 = ggplot(proteins.FC.f.stats.long, aes(x=factor(min.paths), y=value)) +
  geom_boxplot() +
  geom_text(data=min.paths.stats, aes(x = factor(min.paths) , y=1, label=symbol), col="red", size=10)+
  xlab("Shortest distance to transcription factor") +
  facet_wrap(~variable, scales="free") 

g = arrangeGrob(p1,p2 , ncol=2, main=paste("Reference:", reference))
plots.list = lappend(plots.list, g)

