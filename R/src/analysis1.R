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

# write.table(x=droplevels(exp_metadata %>% filter(type == "Kinase" | type == "Standard Mix" | type == "Wild Type") %>% select(sample_name)),quote=F, row.names=F, col.names=F,
#             file="selected_samples.txt")

library(dplyr)
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

string.exp = STRING %>% filter(experimental > 900)
string.all = STRING %>% filter(combined_score > 900)

G.string.exp = graph.data.frame(droplevels(string.exp %>% dplyr::select(ORF1, ORF2)), directed=F)
G.string.all = graph.data.frame(droplevels(string.all %>% dplyr::select(ORF1, ORF2)), directed=F)
GRAPH = G.string.exp


#### VAR ####
protein.matrix = proteins.matrix.combat.quant
proteins.FC = proteins.matrix.combat.quant.FC


## -- coverage: basic stats ----
load("./R/objects/essential_ORFs._load_.RData")
load("./R/objects/flux_coupling._load_.RData")
load("./R/objects/model_reaction2ORF._load_.RData")

EC.genes = gene.annotations[gene.annotations$V3 == "EC number",]
KEGG.pathways = distinct(pathway2orf, pathway, ORF) 
KEGG.pathways = merge(KEGG.pathways, dplyr::select(EC.genes, V1,V4) %>% distinct(), by.x="ORF",by.y=c("V4"))
names(KEGG.pathways)[length(KEGG.pathways)] = "EC.number"

KEGG.pathways = KEGG.pathways %>% 
  filter(!is.na(EC.number)) %>%
  group_by(ORF) %>% mutate(n = n())

KEGG.pathways = droplevels(KEGG.pathways) 
  
measured.ORFs = data.frame(ORF = unique(rownames(protein.matrix)))
measured.ORFs = droplevels(merge(measured.ORFs, EC.genes, by.x="ORF", by.y = "V4")[,c(1,2)])
names(measured.ORFs) = c("ORF", "EC.number")

coupled_reactions = as.character(unique(flux_coupling$reactions[grep(pattern="DIR|FULLY", ignore.case=T, x=flux_coupling$biomass)]))
couplings = merge(na.omit(model_reaction2ORF), data.frame(coupled_reactions), by.x = "reaction", by.y = "coupled_reactions")
couplings = merge(couplings, dplyr::select(EC.genes, V1,V4), by.y=c("V4"), by.x="ORF")
names(couplings) = c("ORF", "reaction", "EC.number")

couplings.ec = unique(as.character(couplings$EC.number))

KEGG.pathways = tbl_df(KEGG.pathways) %>% 
  mutate(is.inData    = ifelse(ORF %in% measured.ORFs$ORF, 1, 0),
         is.Coupled   = ifelse(ORF %in% couplings$ORF, 1, 0),
         is.Essential = ifelse(ORF %in% essential_ORFs$ORF_name, 1, 0))

KEGG.pathways = merge(KEGG.pathways, kegg_categories, by="pathway")

  
file_name = paste("KEGG.pathways", fun_name, "RData", sep=".")
file_path = paste(output_dir, file_name, sep="/")
save(KEGG.pathways, file=file_path)

KEGG.pathways.stats = KEGG.pathways %>%
  group_by(pathway) %>% 
  summarize(EC.coverage  = length(unique(EC.number[is.inData == 1]))/length(unique(EC.number)),
            ORF.coverage = length(unique(ORF[is.inData == 1]))/length(unique(ORF)),
            EC.active  = length(unique(EC.number[is.Coupled == 1 & is.inData == 1 ]))/length(unique(EC.number)),
            ORF.active = length(unique(ORF[is.Coupled == 1 & is.inData == 1]))/length(unique(ORF)))

KEGG.pathways.stats = merge(KEGG.pathways.stats, kegg_categories, by="pathway")
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

p = ggplot(KEGG.pathways.stats.f, aes(x=pathway, y=EC.active, fill=B)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels = KEGG.pathways.stats.f$C ) +
  coord_flip()

plots.list = lappend(plots.list, p)

KEGG.pathways.stats.f = KEGG.pathways.stats.f %>% group_by(B) %>% distinct(pathway) %>% mutate(avgORF = mean(ORF.coverage, na.rm=T)) %>% arrange(ORF.coverage) 
KEGG.pathways.stats.f$pathway = factor(KEGG.pathways.stats.f$pathway, levels = as.character(KEGG.pathways.stats.f$pathway))

file_name = paste("KEGG.pathways.stats.f", fun_name, "RData", sep=".")
file_path = paste(output_dir, file_name, sep="/")
save(KEGG.pathways.stats.f, file=file_path)


p = ggplot(KEGG.pathways.stats.f, aes(x=pathway, y=ORF.coverage, fill=B)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels = KEGG.pathways.stats.f$C ) +
  coord_flip() + theme_light()

toPlot = KEGG.pathways.stats.f %>% group_by(B) %>% 
  summarise(avgEC = mean(EC.coverage),
            avgORF = mean(ORF.coverage),
            avgActive = mean(EC.active))

p = ggplot(toPlot, aes(x=B, y=avgEC, fill=B)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels = toPlot$B ) +
  coord_flip() + theme_light()


plots.list = lappend(plots.list, p)

toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Carbohydrate metabolism"))
p1 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip()
plots.list = lappend(plots.list, p1)

toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Amino acid metabolism"))
p2 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip()
plots.list = lappend(plots.list, p2)

toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Energy metabolism"))
p3 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip()
plots.list = lappend(plots.list, p3)  

toPlot = droplevels(filter(KEGG.pathways.stats.f, B == "Lipid metabolism"))
p4 = ggplot(toPlot, aes(x=pathway, y=EC.coverage)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=percent, limits = c(0,1)) +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip() 

plots.list = lappend(plots.list, p4)



# -- coverage: 4 together ----
select.tmp = data.frame(cat = c("Carbohydrate metabolism", "Amino acid metabolism", "Lipid metabolism", "Energy metabolism"),
                        order = c(1,2,3,4))

toPlot = KEGG.pathways.stats.f[KEGG.pathways.stats.f$B %in% select.tmp$cat ,]
toPlot$B = factor(toPlot$B, levels = select.tmp$cat)
toPlot$order = select.tmp$order[match(toPlot$B, select.tmp$cat)]

toPlot =  toPlot[with(toPlot, order(-order, EC.coverage)),]
toPlot$pathway = factor(toPlot$pathway, levels = as.character(toPlot$pathway))

p = ggplot(toPlot, aes(x=pathway, y=EC.coverage, fill=B)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels = toPlot$C) +
  coord_flip() + theme_light()
p
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

#write.table(x=proteins.FC.f %>% filter(isMetabolic == T) %>% select(ORF) %>% distinct(ORF), quote=F, row.names=F, col.names=F, file="metabolic.enzymes.txt")
#write.table(x=unique(KEGG.pathways.f$ORF), quote=F, row.names=F, col.names=F, file="all_metabolic.enzymes.txt")

toPlot = proteins.FC.f
toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)
toPlot$label = factor(orf2name$gene_name[match(toPlot$KO, orf2name$ORF)])

toPlot = toPlot %>% group_by(KO) %>% mutate(perturbation = abs(min(logFC)) + max(logFC)) %>% ungroup() %>% arrange(-perturbation)
toPlot$label = factor(toPlot$label, levels=unique(toPlot$label))


p3 = ggplot(toPlot, aes(x=label, y = logFC)) + geom_boxplot() + 
          geom_point(data=filter(toPlot,sign == 1), color="cyan", alpha=0.5) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          panel_border()
plots.list = lappend(plots.list, p3)


toPlot = filter(toPlot, isMetabolic==T)
p4 = ggplot(toPlot, aes(x=KO, y = logFC)) + geom_boxplot() + 
            geom_point(data=filter(toPlot,sign == 1), color="cyan", alpha=0.5) +
            scale_x_discrete(labels=toPlot$label)+
            ggtitle("Metabolic") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

t.test(abs(toPlot$logFC[toPlot$isMetabolic == T]), abs(toPlot$logFC[toPlot$isMetabolic == F]))
plots.list = lappend(plots.list, p4)



# pathways.FC = merge(KEGG.pathways.f, proteins.FC.f, by="ORF")
# toPlot = pathways.FC
# p = ggplot(toPlot, aes(x=C, y=logFC, fill=B)) + 
#                       geom_boxplot()
                      
## -- coverage of essential, flux coupling, sentinels ####

load("./R/objects/essential_ORFs._load_.RData")
load("./R/objects/flux_coupling._load_.RData")
load("./R/objects/model_reaction2ORF._load_.RData")

essential_all = unique(as.character(essential_ORFs$ORF_name[essential_ORFs$ORF_name %in% KEGG.pathways.f$ORF]))
EC.essential_all = merge(data.frame(essential_all), dplyr::select(EC.genes, V1,V4), by.y=c("V4"), by.x="essential_all")
names(EC.essential_all) = c("ORF", "EC.number")


coupled_reactions = as.character(unique(flux_coupling$reactions[grep(pattern="DIR|FULLY", ignore.case=T, x=flux_coupling$biomass)]))
couplings = merge(na.omit(model_reaction2ORF), data.frame(coupled_reactions), by.x = "reaction", by.y = "coupled_reactions")

couplings = merge(couplings, dplyr::select(EC.genes, V1,V4), by.y=c("V4"), by.x="ORF")
names(couplings) = c("ORF", "reaction", "EC.number")


measured.ORFs = merge(measured.ORFs, model_reaction2ORF, all.x=T)

coverage.stats = data.frame(essential_ORF = sum(unique(rownames(protein.matrix)) %in% essential_all),
                            total_essental_ORF = length(unique(essential_all)),
                            essential_EC = sum(unique(measured.ORFs$EC.number) %in% unique(EC.essential_all$EC.number)),
                            total_essential_EC = length(unique(EC.essential_all$EC.number)),
                            
                            coupled_ORF = sum(unique(rownames(protein.matrix)) %in% unique(couplings$ORF)),
                            total_coupled_ORF = length(unique(couplings$ORF)),
                            
                            coupled_EC = sum(unique(measured.ORFs$EC.number) %in% unique(couplings$EC.number)),
                            total_coupled_EC = length(unique(couplings$EC.number)),
                            
                            coupled_reactions = sum(unique(measured.ORFs$EC.nu) %in% unique(couplings$EC.number)),
                            total_coupled_EC = length(unique(couplings$EC.number))  )


coverage.stats$fraction_essential_ORF = with(coverage.stats, essential_ORF/total_essental_ORF)
coverage.stats$fraction_essential_EC  = with(coverage.stats, essential_EC/total_essential_EC)
coverage.stats$fraction_coupled_ORF  = with(coverage.stats, coupled_ORF/total_coupled_ORF)
coverage.stats$fraction_coupled_EC  = with(coverage.stats, coupled_EC/total_coupled_EC)

tmp.df = data.frame(fraction = t(coverage.stats %>% dplyr::select(fraction_essential_ORF, fraction_essential_EC, fraction_coupled_ORF, fraction_coupled_EC)),
                    type="present")
tmp.df$id = rownames(tmp.df)

coverage.fractions = rbind(tmp.df, data.frame(fraction = (1 - tmp.df[,1]), type = "absent", id = rownames(tmp.df)))

p = ggplot(coverage.fractions, aes(x = "" , y = fraction, fill = type)) +
           geom_bar(width=1,stat = "identity") +
           scale_fill_manual(values = c("red", "yellow")) +
           scale_y_continuous(labels = percent_format()) +
           coord_polar("y") + facet_wrap(~id, ncol=2) + theme_gray() 

p = p + geom_text(data=coverage.fractions, 
              aes(y = fraction, label = percent(fraction)))


plots.list = lappend(plots.list, p)


KEGG.pathways.f

## -- comparison all vs all protein changes ----

all_KOs = unique(proteins.FC$KO)
kinases = as.character(unique(exp_metadata$ORF[exp_metadata$type == "Kinase"]))
all_KOs.kinases = all_KOs[all_KOs %in% kinases]


ind = 1:length(all_KOs.kinases)
tmp.ind  = t(combn(ind, 2))
tmp.names = matrix(all_KOs.kinases[as.vector(tmp.ind)], ncol=2, byrow=F) 
comparisons = tmp.names

# -- PPI similarities with intersection --

string.exp = STRING %>% filter(experimental > 900)
string.all = STRING %>% filter(combined_score > 900)

G.string.exp = graph.data.frame(droplevels(string.exp %>% dplyr::select(ORF1, ORF2)), directed=F)
G.string.all = graph.data.frame(droplevels(string.all %>% dplyr::select(ORF1, ORF2)), directed=F)
GRAPH = G.string.exp
GRAPH.all = G.string.all

#names(comparisons) = c("KO1", "KO2")


int_abs = rep(0, nrow(comparisons))
int_abs_enzymes = rep(0, nrow(comparisons))
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
  
  KO1_changed = filter(proteins.FC.f,  KO == KO1, p.value_BH < pval_thr)
  KO2_changed = filter(proteins.FC.f,  KO == KO2, p.value_BH < pval_thr)
  
  KO1_changed.enzymes = filter(KO1_changed, isMetabolic == T)
  KO2_changed.enzymes = filter(KO2_changed, isMetabolic == T)
  
  int_abs[i]   = length(intersect(KO1_changed$ORF[abs(KO1_changed$logFC) > FC_thr ], KO2_changed$ORF[abs(KO2_changed$logFC) > FC_thr ]))
  
  int_abs_enzymes[i]   = length(intersect(KO1_changed.enzymes$ORF[abs(KO1_changed.enzymes$logFC) > FC_thr ], KO2_changed.enzymes$ORF[abs(KO2_changed.enzymes$logFC) > FC_thr ]))
  
  int_up[i]    = length(intersect(KO1_changed$ORF[KO1_changed$logFC > FC_thr ], KO2_changed$ORF[KO2_changed$logFC > FC_thr ]))
  int_down[i]  = length(intersect(KO1_changed$ORF[KO1_changed$logFC < -FC_thr ], KO2_changed$ORF[KO2_changed$logFC < -FC_thr ]))
  union_abs[i] = length(union(KO1_changed$ORF, KO2_changed$ORF))
  union_up[i]  = length(union(KO1_changed$ORF[KO1_changed$logFC > FC_thr ], KO2_changed$ORF[KO2_changed$logFC > FC_thr ]))
  union_down[i] = length(union(KO1_changed$ORF[KO1_changed$logFC < -FC_thr ], KO2_changed$ORF[KO2_changed$logFC < -FC_thr ]))
  
  if (KO1 %in% V(GRAPH)$name & KO2 %in% V(GRAPH)$name) {
      int_ppi[i] = length(intersect(V(GRAPH)$name[neighbors(GRAPH,v=KO1)], V(GRAPH)$name[neighbors(GRAPH,v=KO2)]))
  } else {
    int_ppi[i] = NA  
  } 
}

comparisons.df = as.data.frame(cbind(comparisons, int_abs, int_abs_enzymes, int_up, int_down, union_abs, union_up, union_down, int_ppi))


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

#intersection of enzymes
tmp = dcast(comparisons.df, formula=V1~V2, value.var="int_abs_enzymes")
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

comparisons.matrix.enzymes = matrix(as.numeric(comparisons.matrix.tmp), ncol=ncol(comparisons.matrix.tmp), byrow=T)
colnames(comparisons.matrix.enzymes) = colnames(comparisons.matrix.tmp)
rownames(comparisons.matrix.enzymes) = rownames(comparisons.matrix.tmp)

#intersection of up 
tmp = dcast(comparisons.df, formula=V1~V2, value.var="int_up")
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

comparisons.matrix.up = matrix(as.numeric(comparisons.matrix.tmp), ncol=ncol(comparisons.matrix.tmp), byrow=T)
colnames(comparisons.matrix.up) = colnames(comparisons.matrix.tmp)
rownames(comparisons.matrix.up) = rownames(comparisons.matrix.tmp)

#intersection of down 
tmp = dcast(comparisons.df, formula=V1~V2, value.var="int_down")
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

comparisons.matrix.down = matrix(as.numeric(comparisons.matrix.tmp), ncol=ncol(comparisons.matrix.tmp), byrow=T)
colnames(comparisons.matrix.down) = colnames(comparisons.matrix.tmp)
rownames(comparisons.matrix.down) = rownames(comparisons.matrix.tmp)



##-- heatmap of intersection ----

toPlot = comparisons.matrix
annotation = data.frame(class = droplevels(kinase_classes$class[match(colnames(toPlot), kinase_classes$ORF)]),
                        batch_date = exp_metadata$batch_date[match(colnames(toPlot), exp_metadata$ORF)])

annotation$class = as.character(annotation$class)
annotation$class[is.na(annotation$class)] = "unclassified"
annotation$class = factor(annotation$class)

rownames(toPlot) = exp_metadata$gene[match(rownames(toPlot), exp_metadata$ORF)]
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]

file_name = paste("all_intersections", "heatmap.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

parts = 5
S = pheatmap(toPlot,  annotation_col=annotation, cutree_cols=parts, cutree_rows=parts, colorRampPalette(c("white", "black"))(10),
             cellwidth = 8, cellheight = 8, filename=file_path)
             #cellwidth = 8, cellheight = 8)


a = hclust(dist(comparisons.matrix.down))
b = hclust(dist(comparisons.matrix.up))

proteins.FC.f.stats = proteins.FC.f %>% group_by(KO)%>% filter(p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% summarize(n=n(),
                                                                                                                         n_metabolic = sum(isMetabolic == T),
                                                                                                                         n_up = sum(logFC > FC_thr),
                                                                                                                         n_down = sum(logFC < FC_thr))

proteins.FC.f.stats$gene_name = orf2name$gene_name[match(proteins.FC.f.stats$KO, orf2name$ORF)]
#proteins.FC.f.stats$gene_name = factor(proteins.FC.f.stats$gene_name, levels = S$tree_row$labels[S$tree_col$order])

proteins.FC.f.stats$part = cutree(S$tree_row, parts)[match(proteins.FC.f.stats$gene_name, S$tree_row$labels)]
proteins.FC.f.stats$metabolic.fraction = proteins.FC.f.stats$n_metabolic/proteins.FC.f.stats$n
proteins.FC.f.stats = proteins.FC.f.stats[with(proteins.FC.f.stats, order(-n)),]

toPlot = comparisons.matrix
rownames(toPlot) = exp_metadata$gene[match(rownames(toPlot), exp_metadata$ORF)]
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]

dend = as.dendrogram(hclust(dist(toPlot)))
dent.ordered = rotate(dend, proteins.FC.f.stats$gene_name)
plot(dent.ordered)

p = recordPlot()
plots.list = lappend(plots.list, p)


p = ggplot(proteins.FC.f.stats, aes(x=gene_name, y=n, fill=as.factor(part))) + 
  geom_bar(stat="identity", width=.5) + 
  coord_flip() +
  theme_classic()

file_name = paste("all_protein_sets", "barplot.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(p, filename=file_path, width=8.27, height=11.69)





## -- up/down heatmap ----
proteins.FC.f.stats = proteins.FC.f.stats[match(labels(dent.ordered), proteins.FC.f.stats$gene_name),]
proteins.FC.f.stats$gene_name = factor(proteins.FC.f.stats$gene_name, levels = proteins.FC.f.stats$gene_name)

comparisons.matrix.down = comparisons.matrix.down[match(proteins.FC.f.stats$KO, rownames(comparisons.matrix.down)),]
comparisons.matrix.down = comparisons.matrix.down[,match(proteins.FC.f.stats$KO, colnames(comparisons.matrix.down))]

comparisons.matrix.up = comparisons.matrix.up[match(proteins.FC.f.stats$KO, rownames(comparisons.matrix.up)),]
comparisons.matrix.up = comparisons.matrix.up[,match(proteins.FC.f.stats$KO, colnames(comparisons.matrix.up))]

zeros.matrix = matrix(data=0, nrow=nrow(comparisons.matrix.down), 
                       ncol=ncol(comparisons.matrix.down))

colnames(zeros.matrix) = colnames(comparisons.matrix.up)
rownames(zeros.matrix) = rownames(comparisons.matrix.up)

zeros.matrix[upper.tri(comparisons.matrix.up)] = comparisons.matrix.up[upper.tri(comparisons.matrix.up)]
zeros.matrix[lower.tri(comparisons.matrix.down)] = -comparisons.matrix.down[lower.tri(comparisons.matrix.down)]

toPlot = zeros.matrix
toPlot = toPlot/nrow(protein.matrix)
toPlot = ifelse(toPlot > abs(min(toPlot)), abs(min(toPlot)),toPlot)

annotation = data.frame(class = droplevels(kinase_classes$class[match(colnames(toPlot), kinase_classes$ORF)]),
                        batch_date = exp_metadata$batch_date[match(colnames(toPlot), exp_metadata$ORF)])

annotation$class = as.character(annotation$class)
annotation$class[is.na(annotation$class)] = "unclassified"
annotation$class = factor(annotation$class)

rownames(toPlot) = exp_metadata$gene[match(rownames(toPlot), exp_metadata$ORF)]
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]
rownames(annotation) = rownames(toPlot)

file_name = paste("up_down_intersections", fun_name,"pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")


pheatmap(toPlot, cluster_rows=F, cluster_cols=F, cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(7), filename=file_path)

library(scales)
toPlot = melt(as.data.frame(proteins.FC.f.stats %>% select(KO, gene_name, n_up,n_down)), id.vars=c("KO", "gene_name"))
toPlot$value_fraction = toPlot$value/nrow(protein.matrix)
toPlot$gene_name = factor(toPlot$gene_name, levels=rev(unique(toPlot$gene_name)))

p = ggplot(toPlot, aes(x=gene_name, y=value_fraction, fill=variable) )+ 
  geom_bar(stat="identity", width=.5) + 
  coord_flip() +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,0.7,0.1), labels=percent)

file_name = paste("up_down_sets", fun_name, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(p, filename=file_path, width=8.27, height=11.69)



##-- heatmap of intersection of metabolic genes ----

toPlot = comparisons.matrix.enzymes

annotation = data.frame(class = droplevels(kinase_classes$class[match(colnames(toPlot), kinase_classes$ORF)]),
                        batch_date = exp_metadata$batch_date[match(colnames(toPlot), exp_metadata$ORF)])

annotation$class = as.character(annotation$class)
annotation$class[is.na(annotation$class)] = "unclassified"
annotation$class = factor(annotation$class)

rownames(toPlot) = exp_metadata$gene[match(rownames(toPlot), exp_metadata$ORF)]
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]

file_name = paste("metabolic_intersections", "heatmap.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

parts = 4
S = pheatmap(toPlot,  annotation_col=annotation, cutree_cols=parts, cutree_rows=parts, colorRampPalette(c("white", "black"))(10), 
             cellwidth = 8, cellheight = 8, filename=file_path)



proteins.FC.f.stats = proteins.FC.f %>% group_by(KO)%>% filter(p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% summarize(n=n(),
                                                                                                                         n_metabolic = sum(isMetabolic == T))
proteins.FC.f.stats$gene_name = exp_metadata$gene[match(proteins.FC.f.stats$KO, exp_metadata$ORF)]
proteins.FC.f.stats$gene_name = factor(proteins.FC.f.stats$gene_name, levels = S$tree_row$labels[S$tree_col$order])
proteins.FC.f.stats$part = cutree(S$tree_row, parts)[match(proteins.FC.f.stats$gene_name, S$tree_row$labels)]
proteins.FC.f.stats$metabolic.fraction = proteins.FC.f.stats$n_metabolic/proteins.FC.f.stats$n




p2 = ggplot(proteins.FC.f.stats, aes(x=gene_name, y=n_metabolic, fill=as.factor(part))) + 
      geom_bar(stat="identity", width=.5) + 
      coord_flip() +
      theme_classic()

file_name = paste("metabolic_sets", "barplot.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(p2, filename=file_path, width=8.27, height=11.69)

proteins.FC.f.stats$part = cutree(S$tree_row, parts)[match(proteins.FC.f.stats$gene_name, S$tree_row$labels)]
proteins.FC.f.stats$metabolic.fraction = proteins.FC.f.stats$n_metabolic/proteins.FC.f.stats$n

tmp.comparisons = comparisons.matrix

diag(tmp.comparisons) = NA

proteins.FC.f.stats$mean_int = rowMeans(tmp.comparisons, na.rm=T)[match(proteins.FC.f.stats$KO, rownames(comparisons.matrix))]
proteins.FC.f.stats = proteins.FC.f.stats[order(proteins.FC.f.stats$n),]
proteins.FC.f.stats$degree = degree(GRAPH)[match(proteins.FC.f.stats$KO, names(degree(GRAPH)))]





# -- enrichments of changed genes ----- 

#orf2name = gene.annotations %>% distinct(V4, V6) %>% dplyr::select(V4, V6)
enrich_thr = 0.05
proteins.FC = tbl_df(proteins.FC)
proteins.FC$KO = factor(proteins.FC$KO)

all_enrichments = ddply(filter(proteins.FC.f, isMetabolic==T), .(KO), 
                        .fun = function(z) {
                          #z = proteins.z_scores.long[proteins.z_scores.long$KO == "WT",]
                          signal = z$ORF[z$p.value_BH < pval_thr & abs(z$logFC) > FC_thr]
                          universe =  z$ORF              
                          #                           tmp.p = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.process)
                          #                           tmp.p$type = "process"
                          #                           tmp.f = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.function)
                          #                           tmp.f$type = "function"
                          #                           tmp.c = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.compartment)
                          #tmp.c$type = "compartment"
                          tmp.path = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=pathway2orf[pathway2orf$pathway %in% KEGG.pathways.f$pathway,])
                          tmp.path$type = "kegg_pathways"
                          
                          return(tmp.path)
                          
                        })

#kegg
kegg.enrichments = dcast(droplevels(all_enrichments[all_enrichments$type == "kegg_pathways",]), pathway~KO, value.var="p.value")
kegg.enrichments.matrix = as.matrix(kegg.enrichments[,-1])
rownames(kegg.enrichments.matrix) = kegg.enrichments$pathway

kegg.enrichments.matrix[is.na(kegg.enrichments.matrix)] = 0
kegg.enrichments.matrix[kegg.enrichments.matrix > enrich_thr] = 0
kegg.enrichments.matrix[kegg.enrichments.matrix != 0] = 1
proteins.FC.f.stats$enrichment = as.vector(colSums(kegg.enrichments.matrix)[match(proteins.FC.f.stats$KO, colnames(kegg.enrichments.matrix))])

library(splines)

min.stats = proteins.FC.f.stats[na.omit(which(proteins.FC.f.stats$enrichment == 0)),]
min.stats = min.stats[order(-min.stats$n)[1:3],]
max.stats = proteins.FC.f.stats[order(-proteins.FC.f.stats$enrichment)[1:3],]

p = ggplot(proteins.FC.f.stats, aes(x = n, y = metabolic.fraction)) + 
            geom_point(aes(size= enrichment,colour=factor(part))) +
            geom_text(data=max.stats, aes(x = n, y=metabolic.fraction, label=gene_name)) +
            geom_text(data=min.stats, aes(x = n, y=metabolic.fraction, label=gene_name)) +
            stat_smooth(method = "loess", se=T) +
            scale_size(range = c(5, 15)) + theme(aspect.ratio = 1)

file_name = paste(fun_name,"metabolic.fraction", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27+0.1*8.27, width = 11.69+0.1*11.69)


file_name = paste(fun_name, "mean_int", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

p = ggplot(proteins.FC.f.stats, aes(x = n, y = mean_int )) + 
          geom_point(aes(size= degree,colour=factor(part))) +
          stat_smooth(method = "lm", formula=y~ns(x,2), se=F) +
          scale_size(range = c(5, 15)) + theme(aspect.ratio = 5/8)
ggsave(filename=file_path, plot=p, height=8.27+0.1*8.27, width = 11.69+0.1*11.69)

file_name = paste(fun_name, "mean_similarity", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

proteins.FC.f.stats = proteins.FC.f.stats[order(proteins.FC.f.stats$n),]
p = ggplot(proteins.FC.f.stats, aes(x = 1:nrow(proteins.FC.f.stats), y = mean_int/n )) + 
  geom_point(aes(size=8, colour=factor(part))) +
  stat_smooth(method = "loess") + theme(aspect.ratio = 5/8)
ggsave(filename=file_path, plot=p, height=8.27+0.1*8.27, width = 11.69+0.1*11.69)


    #stat_smooth(method = "lm", formula=y~ns(x,2), se=F)

#calculate number of measured ORFs per KEGG metabolic pathway
KEGG.pathways.f.summary = KEGG.pathways.f[KEGG.pathways.f$ORF %in% rownames(protein.matrix),] %>% group_by(pathway, description) %>% summarize(n=n())
selected.pathways = as.character(KEGG.pathways.f.summary$pathway[KEGG.pathways.f.summary$n >= 5])

kegg.enrichments.matrix.f = kegg.enrichments.matrix[rownames(kegg.enrichments.matrix) %in% selected.pathways,]
kegg.enrichments.matrix.f = kegg.enrichments.matrix.f[rowSums(kegg.enrichments.matrix.f) >=3,]
kegg.enrichments.matrix.f = kegg.enrichments.matrix.f[,colSums(kegg.enrichments.matrix.f) >=3]

pathway2desription = pathway2orf %>% distinct(pathway, description) %>% dplyr::select(pathway, description)
pathway2desription$description = sub(pattern=" - Saccharomyces cerevisiae (budding yeast)", replacement="" , x=pathway2desription$description, fixed=T)


kegg.enrichments.matrix.kinases = kegg.enrichments.matrix.f[,colnames(kegg.enrichments.matrix.f) %in% kinases]
toPlot = kegg.enrichments.matrix.kinases
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]

file_name = paste("kegg_enrichments", "heatmap.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

pheatmap(toPlot, clustering_distance_cols="binary", clustering_distance_rows="binary", 
         cluster_rows=T, col=c("lightgrey","black"),legend=F, cellwidth = 8, cellheight = 8,
         labels_row=pathway2desription$description[match(rownames(toPlot), pathway2desription$pathway)], filename=file_path)
         #labels_row=pathway2desription$description[match(rownames(toPlot), pathway2desription$pathway)])


# all enrichments, for supplementary
kegg.enrichments.matrix.f = kegg.enrichments.matrix[rownames(kegg.enrichments.matrix) %in% selected.pathways,]
kegg.enrichments.matrix.f = kegg.enrichments.matrix.f[rowSums(kegg.enrichments.matrix.f) >0,]
kegg.enrichments.matrix.f = kegg.enrichments.matrix.f[,colSums(kegg.enrichments.matrix.f) >0]

pathway2desription = pathway2orf %>% distinct(pathway, description) %>% dplyr::select(pathway, description)
pathway2desription$description = sub(pattern=" - Saccharomyces cerevisiae (budding yeast)", replacement="" , x=pathway2desription$description, fixed=T)


kegg.enrichments.matrix.kinases = kegg.enrichments.matrix.f[,colnames(kegg.enrichments.matrix.f) %in% kinases]
toPlot = kegg.enrichments.matrix.kinases
colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]

file_name = paste("kegg_enrichments_all", "heatmap.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

pheatmap(toPlot, clustering_distance_cols="binary", clustering_distance_rows="binary", 
         cluster_rows=T, col=c("lightgrey","black"),legend=F, cellwidth = 8, cellheight = 8,
         labels_row=pathway2desription$description[match(rownames(toPlot), pathway2desription$pathway)], filename=file_path)
         #labels_row=pathway2desription$description[match(rownames(toPlot), pathway2desription$pathway)])


#cutting to 4 parts
KO.stats = data.frame(KO = droplevels(exp_metadata$ORF[match(names(cutree(S$tree_col,parts)), exp_metadata$gene)]))
KO.stats$cluster = factor(cutree(S$tree_col,parts))
KO.stats$degree = degree(GRAPH)[match(KO.stats$KO, names(degree(GRAPH)))]
KO.stats$degree.all = degree(GRAPH.all)[match(KO.stats$KO, names(degree(GRAPH.all)))]


p = ggplot(melt(KO.stats, id.vars=c("KO", "cluster")), aes(fill=variable, x=cluster, y=value)) +
          geom_boxplot() 
plots.list = lappend(plots.list, p)

## -- changes per kinase class ----
proteins.FC.f$KO_gene = exp_metadata$gene[match(proteins.FC.f$KO, exp_metadata$ORF)]
proteins.FC.f$degree = degree(GRAPH)[match(proteins.FC.f$KO, names(degree(GRAPH)))]
proteins.FC.f$degree.all = degree(GRAPH.all)[match(proteins.FC.f$KO, names(degree(GRAPH.all)))]
proteins.FC.f$cluster = factor(cutree(S$tree_col,parts)[match(proteins.FC.f$KO_gene, names(cutree(S$tree_col,parts)))])
proteins.FC.f$KO_type = exp_metadata$Type[match(proteins.FC.f$KO, exp_metadata$ORF)]
proteins.FC.f$KO_class = as.character(kinase_classes$class[match(proteins.FC.f$KO, kinase_classes$ORF)])
proteins.FC.f$KO_class[is.na(proteins.FC.f$KO_class)] = "unclassified"
proteins.FC.f$KO_class = factor(proteins.FC.f$KO_class)


proteins.FC.classes.stats = proteins.FC.f %>% filter(isMetabolic == T, p.value_BH < pval_thr, abs(logFC) > FC_thr ) %>% 
                            group_by(KO, KO_class) %>% summarise(changes = n(),
                                                                 degree = degree[1])
p = ggplot(proteins.FC.classes.stats, aes(x=KO_class, y=degree)) +
      geom_boxplot() +
      geom_point(data=proteins.FC.classes.stats, aes(x=jitter(as.numeric(KO_class))))
plots.list = lappend(plots.list, p)

p = ggplot(proteins.FC.classes.stats, aes(x=KO_class, y=changes)) +
      geom_boxplot() +
      geom_point(data=proteins.FC.classes.stats, aes(x=jitter(as.numeric(KO_class))))
plots.list = lappend(plots.list, p)

proteins.FC.clusters.stats = proteins.FC.f %>% filter(isMetabolic == T, p.value_BH < pval_thr, abs(logFC) > FC_thr ) %>% 
                             group_by(KO, cluster) %>% summarise(changes = n(),
                                                                 degree = degree[1])

p = ggplot(proteins.FC.clusters.stats, aes(x=cluster, y=changes)) +
          geom_boxplot() +
          geom_point(data=proteins.FC.clusters.stats, aes(x=jitter(as.numeric(cluster))))

plots.list = lappend(plots.list, p)

p = ggplot(proteins.FC.clusters.stats, aes(x=cluster, y=degree)) +
  geom_boxplot() +
  geom_point(data=proteins.FC.clusters.stats, aes(x=jitter(as.numeric(cluster))))
plots.list = lappend(plots.list, p)

plots.list = lappend(plots.list, p)

tbl.pie = sort(table(proteins.FC.classes.stats$KO_class))
pie(tbl.pie, clockwise=F, init.angle=90, col=brewer.pal(length(tbl.pie),"Set3"))
p = recordPlot()
plots.list = lappend(plots.list, p)


         

# #GO slim process
# process.enrichments = dcast(droplevels(all_enrichments[all_enrichments$type == "process",]),pathway~KO, value.var="p.value", fun.aggregate=min)
# process.enrichments.matrix = as.matrix(process.enrichments[,-1])
# rownames(process.enrichments.matrix) = process.enrichments$pathway
# 
# process.enrichments.matrix[is.na(process.enrichments.matrix)] = 0
# process.enrichments.matrix[process.enrichments.matrix > enrich_thr] = 0
# process.enrichments.matrix[process.enrichments.matrix != 0] = 1
# process.enrichments.matrix = process.enrichments.matrix[rowSums(process.enrichments.matrix) >3,]
# 
# GOprocess2desription = GO_slim.process %>% distinct(pathway, description) %>% dplyr::select(pathway, description)
# 
# process.enrichments.matrix.kinases = process.enrichments.matrix[,colnames(process.enrichments.matrix) %in% kinases]
# toPlot = process.enrichments.matrix.kinases
# colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]
# 
# pheatmap(toPlot, clustering_distance_cols="binary",  clustering_distance_rows="binary", cluster_rows=T,
#          labels_row = GOprocess2desription$description[match(rownames(process.enrichments.matrix), GOprocess2desription$pathway)],col=c("lightgrey","black"),annotation = annotation)
#        
# p = recordPlot()
# plots.list = lappend(plots.list, p)


#PCA by pathways
# proteins.FC.kegg = tbl_df(merge(proteins.FC, pathway2orf, by="ORF"))
# proteins.FC.kegg.fractions = proteins.FC.kegg %>% group_by(KO, pathway) %>% summarise(pathway.fraction = sum(p.value_BH < 0.05)/length(p.value_BH))
# 
# kegg.fractions = dcast(proteins.FC.kegg.fractions, formula=pathway~KO, value.var="pathway.fraction")
# kegg.fractions.matrix = as.matrix(kegg.fractions[,-1])
# rownames(kegg.fractions.matrix) = kegg.fractions$pathway
# kegg.fractions.matrix = kegg.fractions.matrix[rowSums(kegg.fractions.matrix) != 0,]
# 
# toPlot = kegg.fractions.matrix[,colnames(kegg.fractions.matrix) %in% kinases]
# 
# rownames(toPlot) = pathway2desription$description[match(rownames(toPlot), pathway2desription$pathway)]
# colnames(toPlot) = exp_metadata$gene[match(colnames(toPlot), exp_metadata$ORF)]
# s1 = prcomp(t(toPlot))
# biplot(s1, cex=0.66)
# abline(h=0,v=0)
# l1 = sort(s1$rotation[,1],decreasing=T)
# l2 = sort(S1$rotation[,2], decreasing=T)
# p = recordPlot()
# biplot.default
# 
# xPC = 1
# yPC = 2
# idx = match(unique(names(sort(abs(c(s1$rotation[,1], s1$rotation[,2])),decreasing=T)[1:20])), rownames(s1$rotation))
# loads = s1$rotation[idx,c(xPC,yPC)]
# 
# biplot(x=s1$x[,c(xPC,yPC)],y=loads, cex=0.66, ylim = c(-3,2),expand=1.5,
#      xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
#      ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
# abline(h=0,v=0)
# p = recordPlot()
# plots.list = lappend(plots.list, p)
# 
# file_name = paste(fun_name, "report.pdf", sep=".")
# file_path = paste(figures_dir, file_name, sep="/")
# save_plots(plots.list, filename=file_path, type="l")

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


G.phys = graph.data.frame(droplevels(yeast.ppi.p %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)
G.gene = graph.data.frame(droplevels(yeast.ppi.g %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)
G.phys.all  = graph.data.frame(droplevels(yeast.ppi.p.all %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)
G.gene.all  = graph.data.frame(droplevels(yeast.ppi.g.all %>% dplyr::select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)), directed=F)


proteins.FC.f = droplevels(proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]) #only kinases
proteins.FC.f$sign = ifelse(abs(proteins.FC.f$logFC) >= FC_thr & proteins.FC.f$p.value_BH < pval_thr, 1,0)
proteins.FC.f$isMetabolic = proteins.FC.f$ORF %in% unique(KEGG.pathways.f$ORF)
proteins.FC.f$isEnzyme = proteins.FC.f$ORF %in% unique(EC.genes$V4)
#proteins.FC.f$cluster = factor(cutree(S$tree_col,parts)[match(proteins.FC.f$KO_gene, names(cutree(S$tree_col,parts)))])


proteins.FC.f.stats = proteins.FC.f %>% filter(p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% 
                      group_by(KO) %>% dplyr::summarise(changes = n(),
                                                        perturbation = sum(abs(logFC)))


#proteins.FC.f.stats$KO %in% unique(as.vector(GO_slim.raw$V1[grep(x=GO_slim.raw$V5, pattern="transcription factor activity")]))
#proteins.FC.f.stats$KO %in% unique(yeastract)

yeastract$TF_ORF = orf2name$ORF[match(yeastract$TF, orf2name$gene_name)]

yeastract.stats = yeastract %>% group_by(TF, TF_ORF) %>% summarise(n = n()) #number of genes yeastract connected to

TFs = unique(as.vector(GO_slim.raw$V1[grep(x=GO_slim.raw$V5, pattern="nucleic acid binding transcription factor activity")]))

#write.table(x=unique(rownames(protein.matrix)), file="proteins.txt",  quote=F, row.names=F, col.names=F)

phys.TF_ORFs = as.character(TFs[TFs %in% V(GRAPH)$name])
kinases = as.character(unique(proteins.FC.f.stats$KO)[unique(proteins.FC.f.stats$KO) %in% V(GRAPH)$name])

paths = igraph::shortest.paths(graph=GRAPH, v=kinases,  to=phys.TF_ORFs, algorithm="unweighted")

paths.all = get.all.shortest.paths(GRAPH, from=kinases, to=phys.TF_ORFs)



min.paths = apply(paths, 1, min)
min.paths[is.infinite(min.paths)] = NA

paths.long = melt(paths) 
names(paths.long) = c("KO", "TF","value")
paths.long$value[is.infinite(paths.long$value)] = NA

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

p3 = ggplot(proteins.FC.f.stats.long, aes(x=factor(cluster), y=min.paths)) +
            geom_boxplot()


g = arrangeGrob(p1,p2 , ncol=2, main=paste("Reference:", reference))
plots.list = lappend(plots.list, g)

p = ggplot(proteins.FC.f.stats.long %>% filter(min.paths.comb != "none") , aes(x=factor(min.paths.comb), y=value)) +
  geom_boxplot() +
  geom_text(data=min.paths.comb.stats %>% filter(min.paths.comb != "none"), aes(x = factor(min.paths.comb) , y=1, label=symbol), col="red", size=10)+
  scale_x_discrete(labels=c("First", ">1")) +
  xlab("Shortest distance to transcription factor") +
  facet_wrap(~variable, scales="free")        
plots.list = lappend(plots.list, p)


proteins.FC.f.stats.long2 = melt(proteins.FC.f.stats, id.vars=c("KO", "changes", "perturbation", "gene_name"))


proteins.FC.f.stats.long2$value = as.numeric(proteins.FC.f.stats.long2$value)
proteins.FC.f.stats.long2.stats = proteins.FC.f.stats.long2 %>% group_by(variable) %>% summarize(cor.changes = cor.test(changes, as.numeric(value))$estimate,
                                                                                                 pval.changes = cor.test(changes, as.numeric(value))$p.value,
                                                                                                 x_min = min(changes,na.rm=T) + min(changes,na.rm=T)*0.3,
                                                                                                 y_max = max(value,na.rm=T) - max(value,na.rm=T)*0.1,
                                                                                                 cor.label = paste(paste("r =",round(cor.changes,2)),
                                                                                                                   paste("p = ",round(pval.changes,2)),sep="\n"))


p = ggplot(proteins.FC.f.stats.long2  , aes(x=perturbation, y=value)) +
    geom_point() +
    geom_text(data=proteins.FC.f.stats.long2.stats, aes(x = x_min, y=y_max, label=cor.label)) +
    facet_wrap(~variable, scales="free") +
    theme(aspect.ratio = 1)
plots.list = lappend(plots.list, p)


## -- sentinels coverage  ----

load("./R/objects/sentinels.table._clean_.RData")
load("./R/objects/sentinelsSRM._clean_.RData")
SRM.sentinels <- read.delim("./data/2015-10-20/Sentinels_matrix_summary_tripl.tsv")


sentinels.table$isMeasured = ifelse(sentinels.table$ORF.ID %in% rownames(protein.matrix), 1, 0)
sentinelsSRM$isSRM = ifelse(sentinelsSRM$ID.ORF %in% rownames(protein.matrix), 1, 0)

#THIS NUMBER, SRM sentinels detected in our DATASET 
sum(unique(sentinelsSRM$ID.ORF) %in% SRM.sentinels$ORF.ID)/(length(unique(sentinelsSRM$ID.ORF))-1)


sum(unique(sentinels.table$ORF.ID[sentinels.table$Sentinel.Grade == "A"])  %in% rownames(protein.matrix))/length(unique(sentinels.table$ORF.ID[sentinels.table$Sentinel.Grade == "A"]))
sum(unique(sentinels.table$ORF.ID[sentinels.table$Sentinel.Grade == "B"])  %in% rownames(protein.matrix))/length(unique(sentinels.table$ORF.ID[sentinels.table$Sentinel.Grade == "B"]))

toPlot =  dplyr::select(sentinels.table, ORF.ID, Sentinel.Grade, isMeasured, isSRM) %>% distinct(ORF.ID, Sentinel.Grade, isMeasured)
p = ggplot(toPlot , aes(fill=factor(isMeasured), x = Sentinel.Grade,
                   y = (..count..)/sum(..count..))) + 
       geom_bar() + 
       stat_bin(geom = "text",   aes(label = paste(round((..count..)/sum(..count..)*100), "%"))) +
       scale_y_continuous(labels = percent)

plots.list = lappend(plots.list, p)



### examples #####
#HOG
example1 = proteins.FC.f[proteins.FC.f$KO %in% orf2name[orf2name$gene_name %in% c("SLT2", "HOG1"),]$ORF,]
example1$gene_name = orf2name$gene_name[match(example1$ORF, orf2name$ORF)]
example1$sign = ifelse(abs(example1$logFC) >= FC_thr & example1$p.value_BH < pval_thr, 1,0)
example1.df = dcast(example1 %>% filter(sign == 1), KO~gene_name, value.var="logFC")

example1.matrix = as.matrix(example1.df[,-1])
example1.df$gene_name = orf2name$gene_name[match(example1.df$KO, orf2name$ORF)]
rownames(example1.matrix) = example1.df[,"gene_name"]

file_name = paste("HOG1_SLT2", "example.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
pheatmap(t(na.omit(t(example1.matrix))), cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(51), filename=file_path)

#STE7/STE20
example1 = proteins.FC.f[proteins.FC.f$KO %in% orf2name[orf2name$gene_name %in% c("STE7", "STE20"),]$ORF,]
example1$gene_name = orf2name$gene_name[match(example1$ORF, orf2name$ORF)]
example1$sign = ifelse(abs(example1$logFC) >= FC_thr & example1$p.value_BH < pval_thr, 1,0)
example1.df = dcast(example1 %>% filter(sign == 1), KO~gene_name, value.var="logFC")

example1.matrix = as.matrix(example1.df[,-1])
example1.df$gene_name = orf2name$gene_name[match(example1.df$KO, orf2name$ORF)]
rownames(example1.matrix) = example1.df[,"gene_name"]

file_name = paste("STE7_STE20", "example.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
pheatmap(t(na.omit(t(example1.matrix))), cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(51), filename=file_path)

#GCN2/SNF1
example1 = proteins.FC.f[proteins.FC.f$KO %in% orf2name[orf2name$gene_name %in% c("GCN2", "SNF1"),]$ORF,]
example1$gene_name = orf2name$gene_name[match(example1$ORF, orf2name$ORF)]
example1$sign = ifelse(abs(example1$logFC) >= FC_thr & example1$p.value_BH < pval_thr, 1,0)
example1.df = dcast(example1 %>% filter(sign == 1), KO~gene_name, value.var="logFC")

example1.matrix = as.matrix(example1.df[,-1])
example1.df$gene_name = orf2name$gene_name[match(example1.df$KO, orf2name$ORF)]
rownames(example1.matrix) = example1.df[,"gene_name"]


file_name = paste("GCN3_SNF1", "example.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
pheatmap(t(na.omit(t(example1.matrix))), cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(51), filename=file_path)

#c("FMP48", "FUS3", "MEK1", "PRK1", "TPK3")
example1 = proteins.FC.f[proteins.FC.f$KO %in% orf2name[orf2name$gene_name %in% c("FMP48", "FUS3", "MEK1", "PRK1", "TPK3"),]$ORF,]
example1$gene_name = orf2name$gene_name[match(example1$ORF, orf2name$ORF)]
example1$sign = ifelse(abs(example1$logFC) >= FC_thr & example1$p.value_BH < pval_thr, 1,0)
example1.df = dcast(example1 %>% filter(sign == 1), KO~gene_name, value.var="logFC")

example1.matrix = as.matrix(example1.df[,-1])
example1.df$gene_name = orf2name$gene_name[match(example1.df$KO, orf2name$ORF)]
rownames(example1.matrix) = example1.df[,"gene_name"]


file_name = paste("FMP48", "example.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
pheatmap(t(na.omit(t(example1.matrix))), cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(51), filename=file_path)



# "FUS3",KSS1
example1 = proteins.FC.f[proteins.FC.f$KO %in% orf2name[orf2name$gene_name %in% c("FUS3", "KSS1"),]$ORF,]
example1$gene_name = orf2name$gene_name[match(example1$ORF, orf2name$ORF)]
example1$sign = ifelse(abs(example1$logFC) >= FC_hr & example1$p.value_BH < pval_thr, 1,0)
example1.df = dcast(example1 %>% filter(sign == 1), KO~gene_name, value.var="logFC")
example1.df = dcast(example1, KO~gene_name, value.var="logFC")

example1.df = example1.df[,c(1,which(colnames(example1.df) %in% (example1 %>% filter(sign == 1) %>% select(gene_name))$gene_name))]

example1.matrix = as.matrix(example1.df[,-1])
example1.df$gene_name = orf2name$gene_name[match(example1.df$KO, orf2name$ORF)]
rownames(example1.matrix) = example1.df[,"gene_name"]


file_name = paste("KSS1_FUS3", "example.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
pheatmap(t(na.omit(t(example1.matrix))), cellheight=8, cellwidth=2,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(51), filename=file_path)
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(51))


### network example


current_nodes = orf2name$ORF[match(c("FMP48", "FUS3", "MEK1", "PRK1", "TPK3"), orf2name$gene_name)]

SUB = induced.subgraph(G.string.all, unique(unlist(neighborhood(G.string.all, order=1, nodes=current_nodes))))
V(SUB)$label  = orf2name$gene_name[match(V(SUB)$name, orf2name$ORF)]

SUB = simplify(as.undirected(SUB, mode="collapse", edge.attr.comb="random"))

plot.igraph(SUB, vertex.label=V(SUB)$label, vertex.size=10,
                 layout=layout.fruchterman.reingold)


file_name = paste("example", fun_name, "graphml", sep=".")
file_path = paste(output_dir, file_name, sep="/")

write.graph(SUB, file=file_path, format="graphml")

### kinase network

d = degree(GRAPH, mode="all")
if (sum(d == 0) > 0) {
  d = d + 1
} 

fit2 <- power.law.fit(d,  implementation="plfit")
dd = degree.distribution(GRAPH, cumulative=F, mode="all")




file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")


