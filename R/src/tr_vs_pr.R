rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")
library("cowplot")
plots.list = list()
fun_name = "tr_vs_pr"

load("./R/objects/proteins.matrix.combat.quant.FC.RData")
load("./R//objects/transcriptome.FC._clean_.RData")
load("./R/objects/orf2ko._load_.RData")

load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/pathway2orf._load_.RData")
load("./R/objects/pathways._load_.RData")
load("./R/objects/gene.annotations._load_.RData")
load("./R/objects/kegg_categories._load_.RData")
load("./R/objects/iMM904._load_.RData")

orf2name = droplevels(unique(gene.annotations[,c("V4", "V6")]))
orf2name$V4 = as.character(orf2name$V4)
orf2name$V6 = as.character(orf2name$V6)
orf2name$V6[orf2name$V6 == ""] = orf2name$V4[orf2name$V6 == ""]
names(orf2name) = c("ORF", "gene_name")
EC.genes = gene.annotations[gene.annotations$V3 == "EC number",]

KEGG.pathways = distinct(pathway2orf, pathway, ORF) 
#KEGG.pathways$EC.number = EC.genes$V1[match(KEGG.pathways$ORF, EC.genes$V4)]

KEGG.pathways = merge(KEGG.pathways, dplyr::select(EC.genes, V1,V4), by.x="ORF",by.y=c("V4"))
names(KEGG.pathways)[length(KEGG.pathways)] = "EC.number"

KEGG.pathways = droplevels(KEGG.pathways %>% filter( !is.na(EC.number)))
KEGG.pathways = KEGG.pathways %>% group_by(ORF) %>% mutate(n = n())

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


## -- transcriptome vs proteome

transcriptome.FC$KO = NULL
transcriptome.FC = transcriptome.FC[,c("ORF", "logFC", "p.value", "p.value_BH", "KO.ORF")]
names(transcriptome.FC)[length(transcriptome.FC)] = "KO"
transcriptome.FC$isMetabolic = transcriptome.FC$ORF %in% unique(KEGG.pathways.f$ORF)
transcriptome.FC$isiMM904 = transcriptome.FC$ORF %in% unique(as.character(iMM904$gene))
transcriptome.FC.f = transcriptome.FC[transcriptome.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]


tr.pr.FC = merge(transcriptome.FC, proteins.matrix.combat.quant.FC, by=c("KO", "ORF"), suffixes=c(".tr", ".pr"))

pval_thr = 0.05

tr.pr.cor.ORF = tr.pr.FC %>% group_by(ORF, isMetabolic) %>% 
  summarise(cor = cor(logFC.tr, logFC.pr, method="spearman"))

tr.pr.cor.KO = tr.pr.FC %>% group_by(KO, isMetabolic) %>% 
  summarise(cor = cor(logFC.tr, logFC.pr, method="spearman"))

p = ggplot(tr.pr.cor.KO, aes(x=cor, fill = isMetabolic)) + 
    geom_density(alpha=0.25)
plots.list = lappend(plots.list, p)


#abs(logFC.tr) < log2(1.7),
tr.pr.cor.ORF.f = tr.pr.FC %>% group_by(ORF) %>% filter(p.value.tr < pval_thr) %>%
  summarise(cor = cor(logFC.tr, logFC.pr, method="pearson"),
            n = n())

tr.pr.cor.KO.f = tr.pr.FC %>% group_by(KO) %>% filter(p.value.tr < pval_thr ) %>%
  summarise(cor = cor(logFC.tr, logFC.pr, method="pearson"),
            n = n())

toPlot = tr.pr.cor.KO.f[tr.pr.cor.KO.f$n>=10,] %>% arrange(cor)
toPlot$KO.name = orf2name$gene_name[match(toPlot$KO, orf2name$ORF)]
toPlot$KO.name = factor(toPlot$KO.name, levels=unique(toPlot$KO.name))


p = ggplot(toPlot, aes(x = KO.name, y = cor,  size=n)) + 
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.5,-0.25,0.25, 0.5), linetype = 3) +  
  ylim(-0.8,0.8) +
  coord_flip()
p
plots.list = lappend(plots.list, p)



## -- intersections ----

all_KOs = unique(sort(as.character(transcriptome.FC$KO)))
kinases = as.character(unique(exp_metadata$ORF[exp_metadata$type == "Kinase"]))
all_KOs.kinases = all_KOs[all_KOs %in% kinases]


ind = 1:length(all_KOs.kinases)
tmp.ind  = t(combn(ind, 2))
tmp.names = matrix(all_KOs.kinases[as.vector(tmp.ind)], ncol=2, byrow=F) 
comparisons = tmp.names



int_abs = rep(0, nrow(comparisons))
int_up = rep(0, nrow(comparisons))
int_down = rep(0, nrow(comparisons))

union_abs = rep(0, nrow(comparisons))
union_up = rep(0, nrow(comparisons))
union_down = rep(0, nrow(comparisons))

FC_thr = log2(1)

for (i in 1:nrow(comparisons)) {
  pair = comparisons[i,]
  KO1 = pair[1]
  KO2 = pair[2]
    
  KO1_changed = filter(transcriptome.FC.f,  KO == KO1, p.value_BH < pval_thr)
  KO2_changed = filter(transcriptome.FC.f,  KO == KO2, p.value_BH < pval_thr)
  
  
  int_abs[i]   = length(intersect(KO1_changed$ORF[abs(KO1_changed$logFC) > FC_thr ], KO2_changed$ORF[abs(KO2_changed$logFC) > FC_thr ]))
  #int_up[i]    = length(intersect(KO1_changed$ORF[KO1_changed$logFC > FC_thr ], KO2_changed$ORF[KO2_changed$logFC > FC_thr ]))
  #int_down[i]  = length(intersect(KO1_changed$ORF[KO1_changed$logFC < -FC_thr ], KO2_changed$ORF[KO2_changed$logFC < -FC_thr ]))
  #union_abs[i] = length(union(KO1_changed$ORF, KO2_changed$ORF))
  #union_up[i]  = length(union(KO1_changed$ORF[KO1_changed$logFC > FC_thr ], KO2_changed$ORF[KO2_changed$logFC > FC_thr ]))
  #union_down[i] = length(union(KO1_changed$ORF[KO1_changed$logFC < -FC_thr ], KO2_changed$ORF[KO2_changed$logFC < -FC_thr ]))
}

#comparisons.df = as.data.frame(cbind(comparisons, int_abs, int_up, int_down, union_abs, union_up, union_down))
comparisons.df = as.data.frame(cbind(comparisons, int_abs))



#making matrix all vs all
#intersection
tmp = reshape2::dcast(comparisons.df, formula=V1~V2, value.var="int_abs")

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

toPlot = comparisons.matrix

rownames(toPlot) = exp_metadata$gene[match(rownames(toPlot), exp_metadata$ORF)]
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]



transcriptome.FC.f.stats = transcriptome.FC.f %>% group_by(KO)%>% filter(p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% summarize(n=n(),
                                                                                                                         n_metabolic = sum(isMetabolic == T))
transcriptome.FC.f.stats$gene_name = exp_metadata$gene[match(transcriptome.FC.f.stats$KO, exp_metadata$ORF)]
transcriptome.FC.f.stats$metabolic.fraction = transcriptome.FC.f.stats$n_metabolic/transcriptome.FC.f.stats$n

tmp.comparisons = comparisons.matrix
diag(tmp.comparisons) = NA

transcriptome.FC.f.stats$mean_int = rowMeans(tmp.comparisons, na.rm=T)[match(transcriptome.FC.f.stats$KO, rownames(comparisons.matrix))]
transcriptome.FC.f.stats = transcriptome.FC.f.stats[order(transcriptome.FC.f.stats$n),]

library(splines)
library(cowplot)
p = ggplot(transcriptome.FC.f.stats, aes(x = n, y = mean_int )) + 
            geom_point() +
            stat_smooth(method = "lm", formula=y~ns(x,2), se=F) +
            ylab("Mean intersection size") +
            xlab("Number of total transcripts changed per mutant") +
            scale_size(range = c(5, 15)) + theme(aspect.ratio = 5/8)
            
plots.list = lappend(plots.list, p)


p = ggplot(transcriptome.FC.f.stats, aes(x = n, y = metabolic.fraction)) + 
            geom_point() +
            stat_smooth(method = "loess", se=F) +
            scale_size(range = c(5, 15)) + theme(aspect.ratio = 1) +
            xlab("Number of total transcripts changed per mutant") +
            ylab("Fraction of metabolic genes affected by kinase deletion")

plots.list = lappend(plots.list, p)


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")





