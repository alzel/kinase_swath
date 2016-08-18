rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "Figure1"

library(scales)
#library(cowplot)

load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/gene.annotations._load_.RData")

orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             sgd = gene.annotations$V5,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]

# KEGG pathway Coverage -------------------------------------

load("./R/objects/KEGG.pathways.analysis1.RData")
load("./R/objects/kegg_categories._load_.RData")

KEGG.pathways.stats <- KEGG.pathways %>%
  group_by(pathway) %>% 
  summarize(EC.coverage  = length(unique(EC.number[is.inData == 1]))/length(unique(EC.number)),
            ORF.coverage = length(unique(ORF[is.inData == 1]))/length(unique(ORF)),
            EC.active  = length(unique(EC.number[is.Coupled == 1 ]))/length(unique(EC.number)),
            ORF.active = length(unique(ORF[is.Coupled == 1 ]))/length(unique(ORF)),
            EC.active.inData  = length(unique(EC.number[is.Coupled == 1 & is.inData == 1 ]))/length(unique(EC.number[is.Coupled == 1])))


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

KEGG.pathways.stats.f  = droplevels(KEGG.pathways.stats[KEGG.pathways.stats$B %in% selected,])
toPlot = KEGG.pathways.stats.f

KEGG.pathways.stats.f.long = KEGG.pathways.stats.f %>% 
  gather(stats, coverage, -pathway, -A, -B, -C)



toPlot = droplevels(KEGG.pathways.stats.f.long %>% filter(stats == "EC.active.inData" | stats == "EC.coverage"))
p = ggplot(toPlot, aes(x=pathway, y=coverage, fill=stats)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_x_discrete(labels = KEGG.pathways.stats.f$C ) +
  coord_flip()



toPlot = KEGG.pathways.stats.f %>% group_by(B) %>% 
  summarise(avgEC = mean(EC.coverage, na.rm = T),
            avgORF = mean(ORF.coverage, na.rm = T),
            avgActive = mean(EC.active, na.rm = T),
            avgActive.inData = mean(EC.active.inData, na.rm = T)) %>%
  ungroup() %>%
  arrange(avgEC)

toPlot = KEGG.pathways.stats.f.long %>% 
  filter(stats == "EC.active.inData" | stats == "EC.coverage") %>%
  group_by(B, stats) %>%
  summarise(avg = mean(coverage, na.rm=T)) %>%
  ungroup() %>%
  arrange(avg)

toPlot = droplevels(toPlot)
toPlot$avg[is.nan(toPlot$avg)] <- NA
arrange_vector = as.character((toPlot %>% filter(stats == "EC.coverage") %>% arrange(avg))$B)
tmp.vector = as.character((toPlot %>% arrange(avg))$B)

toPlot$B = factor(toPlot$B, levels = arrange_vector)
toPlot$stats = factor(toPlot$stats, levels = sort(unique(toPlot$stats), decreasing = T))

toPlot %>% ungroup() %>% arrange(stats)

p.coverage <- ggplot(toPlot, aes(x=B, y=avg, fill=stats)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_y_continuous(breaks = seq(0,1,0.2),labels = percent_format(), limits = c(0,1)) +
  xlab("") +
  ylab("Average coverage by metabolism part") +
  coord_flip() +
  scale_fill_grey(start = 0.8, 
                  end = 0.2,
                  name = "",
                  breaks = c("EC.coverage", "EC.active.inData"),
                  labels = c("All reactions", "Active reactions")) +
  cowplot::background_grid(major = "x", minor = "none") +
  theme_bw() +
  theme(legend.position = c(.66, 0.1)) 

coverage_thr <- 0.1
toPlot <- KEGG.pathways.stats.f %>% 
  group_by(A) %>%
    summarise(pathways_coverage = sum(EC.coverage > coverage_thr)/length(EC.coverage))
toPlot <- rbind.data.frame(toPlot, data.frame(A = "not_covered", pathways_coverage = 1-toPlot$pathways_coverage[1]))

p.pack_man <- ggplot(toPlot, aes(x="", y=pathways_coverage, fill=A)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c("black", "lightgrey")) +
  coord_polar("y", start = pi/2) +
  scale_y_continuous(breaks = seq(0,1,0.2),labels = percent_format(), limits = c(0,1)) +
  labs(title = "Metabolic pathways\ncovered") +
  theme(line = element_blank(),
        legend.justification=c(1,0), 
        legend.position="none",
        legend.key.size	=  unit(0.25, "cm"))+
  xlab("") +
  ylab("")

toPlot <- KEGG.pathways.stats.f %>% 
  dplyr::select(A, EC.coverage, EC.active.inData) %>%
  gather(stats, coverage, -A) %>%
  group_by(stats) %>%
    summarise(MEAN_coverage = mean(coverage, na.rm = T))

toPlot$stats = factor(toPlot$stats, levels = c("EC.active.inData", "EC.coverage"))
p.mean_cov <- ggplot(toPlot, aes(x=stats, y=MEAN_coverage)) +
  geom_bar(width=0.85, stat="identity", position = "dodge") +
  scale_y_continuous(breaks = seq(0,1,0.2),labels = percent_format(), limits = c(0,1)) +
  ylab("Average coverage per pathway") +
  xlab("Mapped reactions") +
  scale_x_discrete(breaks = c("EC.coverage", "EC.active.inData"),
                  labels = c("Overall", "Active")) +
  theme_bw() +
  theme(legend.position = "none")


toPlot <- KEGG.pathways.stats.f %>% 
  dplyr::select(A, EC.coverage, EC.active.inData) %>%
  gather(stats, coverage, -A) %>%
  group_by(stats) %>%
  summarise(MEAN_coverage = mean(coverage, na.rm = T))

toPlot$stats = factor(toPlot$stats, levels = c("EC.active.inData", "EC.coverage"))

toPlot2 <- KEGG.pathways.stats.f %>% 
  group_by(A) %>%
  summarise(pathways_coverage = sum(EC.coverage > coverage_thr)/length(EC.coverage))

toPlot <- rbind(toPlot, rename(toPlot2, stats = A, MEAN_coverage = pathways_coverage) )
toPlot$non_coverered <- 1 - toPlot$MEAN_coverage
toPlot <- reshape::melt(as.data.frame(toPlot), id.vars = c("stats"))

toPlot$stats <- factor(toPlot$stats, levels = as.character((toPlot %>% filter(variable  == "MEAN_coverage") %>% arrange(value) %>% dplyr::select(stats))[,1]))

p.pack_man2 <- ggplot(toPlot, aes(x=stats, y=value, fill = variable)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y", start = pi/2 ) +
  scale_y_continuous(breaks = seq(0,1,0.25),labels = percent_format(), limits = c(0,1)) +
  theme_bw() +
  theme(legend.justification=c(1,0), 
        legend.position="none",
        legend.key.size	=  unit(0.25, "cm"))+
  xlab("") +
  ylab("")


toSave <- KEGG.pathways.stats
toSave$label <- sub(x = toSave$pathway, pattern = "path:sce", replacement = "")

toSave <- toSave %>% filter(EC.coverage > 0.3) %>% dplyr::select(label )

file_name = paste("KEGG_pathways_highlight", fun_name, "csv", sep = ".")
file_path = paste(suppl_dir, file_name, sep="/")
write.csv(x = toSave, file = file_path, row.names = F, quote = F, eol = "\n", col.names = NA)


toSave <- KEGG.pathways %>% filter(is.inData == T)
toSave$label <- paste("EC", as.character(toSave$EC.number), sep = "")
toSave <- toSave %>% dplyr::select(label ) %>% distinct()


file_name = paste("KEGG_EC_highlight", fun_name, "csv", sep = ".")
file_path = paste(suppl_dir, file_name, sep="/")
write.table(x = toSave, file = file_path, row.names = F, quote = F, eol = "\n", col.names = F)

## -- CVs of mix samples -----

load("./R/objects/proteins.matrix.sva.0.5.1.RData")
#load("./R/objects/fragments.matrix.quant.combat.RData")
#load("./R/objects/peptides.matrix.quant.combat.RData")

data.raw = proteins.matrix.sva.0.5.1
data.raw.f = data.raw[,grepl(pattern = "mix", ignore.case = T, x = colnames(data.raw) )]
data.raw.f.df = tbl_df(data.frame(ORF = row.names(data.raw.f), as.data.frame(data.raw.f)))

data.raw.f.proteins = data.raw.f.df %>%
  gather(sample, Tsignal, -ORF)
data.raw.f.proteins$type = "proteins"

# data.raw = peptides.matrix.quant.combat
# data.raw.f = data.raw[,grepl(pattern = "mix", ignore.case = T, x = colnames(data.raw) )]
# data.raw.f.df = tbl_df(data.frame(ORF = row.names(data.raw.f), as.data.frame(data.raw.f)))
# 
# data.raw.f.peptides = data.raw.f.df %>%
#   gather(sample, Tsignal, -ORF)
# data.raw.f.peptides$type = "peptides"
# 
# data.raw = fragments.matrix.quant.combat
# data.raw.f = data.raw[,grepl(pattern = "mix", ignore.case = T, x = colnames(data.raw) )]
# data.raw.f.df = tbl_df(data.frame(ORF = row.names(data.raw.f), as.data.frame(data.raw.f)))
# 
# data.raw.f.fragments = data.raw.f.df %>%
#   gather(sample, Tsignal, -ORF)
# data.raw.f.fragments$type = "fragments"
# 
 #data.all = rbind.data.frame(data.raw.f.fragments, data.raw.f.peptides, data.raw.f.proteins)
 data.all = rbind.data.frame(data.raw.f.proteins)
 data.all$signal = exp(data.all$Tsignal)
 data.all$type = factor(data.all$type, levels = c("fragments", "peptides", "proteins"))

data.all.summary = data.all %>% group_by(ORF, type) %>% summarise(CV = sd(signal, na.rm=T)/mean(signal, na.rm=T))


toPlot <- data.all.summary

p = ggplot(data = toPlot, aes(x=type, y=CV*100)) + 
  geom_violin() + 
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_y_log10(limits=c(5,100), breaks=c(5, 10, 15,20, 50,  100)) +
  annotation_logticks(sides="l") +
  ylab("Signal variation for the duration of experiment, CV")
p 

toPlot <- data.all.summary %>% filter(type == "proteins")
p.CV = ggplot(data = toPlot, aes(x=CV)) + 
  geom_density() + 
  scale_x_log10(limits=c(5,100)/100, breaks=c(5, 10, 15,20, 50,  100)/100, labels = percent_format())+
  annotation_logticks(sides="b") +
  xlab("Signal variation during the duration of experiment, CV%") +
  geom_vline(xintercept = genefilter::shorth(toPlot$CV), linetype=5) +
  theme_bw() + 
  theme(aspect.ratio = 1)
  #background_grid(major = "x", minor="x")

plots.list <- lappend(plots.list, p.CV)

## -- Volcano plots ---- 
# load("./R/objects/proteins.matrix.combat.quant.FC.RData")
# load("./R/objects/proteins.matrix.combat.quant.RData")

load("./R/objects/proteins.matrix.sva.0.5.1.RData")
load("./R/objects/proteins.matrix.sva.0.5.1.FC.RData")
load("./R/objects/iMM904._load_.RData")

protein.matrix = proteins.matrix.sva.0.5.1
proteins.FC = proteins.matrix.sva.0.5.1.FC

reference = unique(as.character(proteins.FC$reference))

lb = -4
ub = 4

proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]

#saving perturbations 
toSave <- proteins.FC.f
toSave$isMetabolic = proteins.FC.f$ORF %in% unique(as.character(iMM904$gene))
toSave$ORF.name <- orf2name$gene_name[match(toSave$ORF, orf2name$ORF)]
toSave$KO.name <- orf2name$gene_name[match(toSave$KO, orf2name$ORF)]

toSave <- toSave %>% 
  dplyr::select(ORF, ORF.name, KO, KO.name, logFC, p.value, p.value_BH, p.value_bonferroni, isMetabolic, reference) %>%
  arrange(KO)

file_name = paste("supplementary_file1", fun_name, "csv", sep = ".")
file_path = paste(suppl_dir, file_name, sep="/")
write.csv(x = toSave, file = file_path, row.names = F,  eol = "\n")

proteins.FC.f$isMetabolic = proteins.FC.f$ORF %in% unique(KEGG.pathways.f$ORF)
proteins.FC.f$isiMM904 = proteins.FC.f$ORF %in% unique(as.character(iMM904$gene))

toPlot = proteins.FC.f %>% filter(isiMM904 == T)
toPlot[toPlot$logFC < 0,]$logFC = ifelse(toPlot[toPlot$logFC < 0,]$logFC < lb, lb, toPlot[toPlot$logFC < 0,]$logFC)
toPlot[toPlot$logFC > 0,]$logFC = ifelse(toPlot[toPlot$logFC > 0,]$logFC > ub, ub, toPlot[toPlot$logFC > 0,]$logFC)

pval_thr = 0.01
set.seed(123)
FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)

proteins.FC.stats = data.frame( ratio_sign = round(sum(proteins.FC$logFC < -FC_thr & proteins.FC$p.value_BH < pval_thr)/sum(proteins.FC$logFC > FC_thr & proteins.FC$p.value_BH < pval_thr),2),
                                ratio =  round(sum(proteins.FC$logFC < 0)/sum(proteins.FC$logFC > 0),2),
                                n_prot = length(unique(proteins.FC$ORF)),
                                n_sign = round(sum(proteins.FC$p.value_BH<pval_thr & abs(proteins.FC$logFC)>FC_thr )/length(unique(proteins.FC$contrasts)),2),
                                x_min = -4,
                                y_max =max(-log10(proteins.FC$p.value_BH)))

toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)



p.volcano1 = ggplot(toPlot, aes(y=-log10(p.value_BH), x=logFC)) +
  geom_point(aes(color = sign), alpha=0.5) + 
  xlim(c(lb,ub)) +
  geom_hline(yintercept = -log(pval_thr,10),linetype=3) +
  geom_vline(xintercept = c(FC_thr,-FC_thr),linetype=3) +
  geom_text(data=proteins.FC.stats, aes(x=x_min, y=y_max, hjust=0, vjust=1,
                                        label=paste(paste0("On average perturbed = ", round(n_sign/n_prot*100, 2), "%"),
                                                    paste0("up/down_ratio = ",ratio), sep="\n"))) +
                                                    #paste0("ratio_sign = ",   ratio_sign),
                                                    #paste0("#_prot = ", n_prot), sep="\n"))) +
  xlab("Protein expression change in mutant, log2(fold-change)") +
  theme_bw() +
  theme(aspect.ratio = 5/8, legend.position = "none")

p.volcano_rect = ggplot(toPlot, aes(y=-log10(p.value_BH), x=logFC)) +
  geom_point(aes(color = sign), alpha=0.5) + 
  xlim(c(lb,ub)) +
  geom_hline(yintercept = -log(pval_thr,10),linetype=3) +
  geom_vline(xintercept = c(FC_thr,-FC_thr),linetype=3) +
#   geom_text(data=proteins.FC.stats, aes(x=x_min, y=y_max, hjust=0, vjust=1,
#                                         label=paste(paste0("On average perturbed = ", round(n_sign/n_prot*100, 2), "%"),
#                                                     paste0("up/down_ratio = ",ratio), sep="\n"))) +
  #paste0("ratio_sign = ",   ratio_sign),
  #paste0("#_prot = ", n_prot), sep="\n"))) +
  xlab("Protein expression change in mutant, log2(fold-change)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        legend.position = "none")



p.volcano_inset = ggplot(toPlot, aes(x=logFC)) +
  geom_histogram(colour = "white", fill = "black", binwidth = 0.25) +
  #geom_rect(aes(ymin=0, ymax=Inf, xmin=-FC_thr, xmax=FC_thr), inherit.aes=F, fill="grey", alpha = 0.01) +
  xlab(paste("Log2(fold-change)")) +
  #geom_vline(xintercept = c(FC_thr,-FC_thr),linetype=2) +
  xlim(c(-2,2)) +
  theme_bw() +
  theme(aspect.ratio = 1)

p.volcano_combined <- p.volcano1 + annotation_custom(grob = ggplotGrob(p.volcano_inset), xmin = 0, xmax = 5, ymin = 20 , ymax = 120)



##  ---- fold-change boxplots for supplementary ------

toPlot = proteins.FC.f

toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)
toPlot$label = orf2name$gene_name[match(toPlot$KO, orf2name$ORF)]

toPlot = toPlot %>% group_by(KO) %>% mutate(perturbation = abs(min(logFC)) + max(logFC)) %>% ungroup() %>% arrange(-perturbation)
#label.genes <- unique(tolower(toPlot$label))
#my_labels = lapply(label.genes, FUN = function(x) expression(paste(italic(get(x), Delta))))
toPlot$label <- factor(toPlot$label, levels = unique(toPlot$label))

s.kinase_boxplot_all <- ggplot(toPlot, aes(x=label, y = logFC)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_point(data=filter(toPlot,sign == 1), color="cyan", alpha=0.5) +
  geom_boxplot(outlier.shape = NA) + 
  xlab("Kinase mutant") +
  ylab("Log2(fold-change)") +
  ylim(-4, 4) +
  #scale_x_discrete(labels = as.vector(my_labels)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

toPlot.stats <- toPlot %>% group_by(KO, isiMM904) %>%
  summarize(diff_range = diff(range(logFC, na.rm=T)),
            diff_IQR = IQR(logFC, na.rm=T),
            n = sum(sign)) %>% ungroup() %>% arrange(-diff_IQR) 
  
#plots.list <- lappend(plots.list, s.kinase_boxplot_all)
toPlot.stats.res <- t.test(toPlot.stats$diff_IQR[toPlot.stats$isiMM904 == T], toPlot.stats$diff_IQR[toPlot.stats$isiMM904 == F])
toPlot$KO <- factor(toPlot$KO, levels = unique(as.character(toPlot.stats$KO)))


toPlot.stats <- toPlot %>%
  filter(isiMM904) %>% 
  group_by(label) %>%
    summarize(n = sum(sign)) %>% 
  ungroup() %>% 
  arrange(-n)
toPlot.stats$label <- factor(toPlot.stats$label, levels = toPlot.stats$label)
toPlot$label <- factor(toPlot$label, levels = toPlot.stats$label)

s.kinase_boxplot_metabolic_only <- ggplot(toPlot %>% filter(isiMM904), aes(x = label, y = logFC)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  #geom_point(data=filter(toPlot,sign == 1), color="cyan", alpha=0.5) +
  geom_boxplot() + 
  #annotate("text", x = 80, y = 2, label = round(toPlot.stats.res$p.value,2)) +
  #geom_text(data = NA, aes(x=1, y= -2, label = round(toPlot.stats.res$p.value,2))) +
  xlab("Kinase mutant") +
  ylab(expression(paste("Protein expression, ", log[2], "(mutant/wild-type)", sep=""))) +
  ylim(-3, 3) +
  #scale_x_discrete(labels = as.vector(my_labels)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"), 
        legend.position=c(0.8, 0.2), aspect.ratio = 5/8)

s.barplot_metabolic_counts <- ggplot(toPlot.stats, aes(x=label, y=n)) + 
  geom_bar(stat="identity", width=.5) +
  ylab("Number of perturbed metabolic enzymes") +
  xlab("Kinase mutant") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"), aspect.ratio = 5/8)
        
s.kinase_boxplot_metabolic_only$toScale <- T
plots.list <- lappend(plots.list, s.kinase_boxplot_metabolic_only)

s.barplot_metabolic_counts$toScale <- T
plots.list <- lappend(plots.list, s.barplot_metabolic_counts)



toPlot.stats <- toPlot %>% group_by(KO) %>%
  summarize(diff_range = diff(range(logFC, na.rm=T)),
            sign = sum(sign),
            diff_IQR = IQR(logFC, na.rm=T)) %>% ungroup() %>% arrange(-diff_IQR) 

toPlot$KO <- factor(toPlot$KO, levels = unique(as.character(toPlot.stats$KO)))
s.kinase_boxplot_all <- ggplot(toPlot, aes(x=KO, y = logFC)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_point(data=filter(toPlot,sign == 1), color="cyan", alpha=0.5) +
  geom_boxplot(outlier.shape = NA) + 
  #annotate("text", x = 80, y = 2, label = round(toPlot.stats.res$p.value,2)) +
  #geom_text(data = NA, aes(x=1, y= -2, label = round(toPlot.stats.res$p.value,2))) +
  scale_fill_grey(name=paste("is metabolic enzyme?",  "p-value", " = ", round(toPlot.stats.res$p.value,2), sep = ""),
                  breaks=c(T, F),
                  labels=c("TRUE", "FALSE")) +
  xlab("Kinase mutant") +
  ylab(expression(paste("Protein expression, ", log[2], "(mutant/wild-type)", sep=""))) +
  ylim(-4, 4) +
  #scale_x_discrete(labels = as.vector(my_labels)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position=c(0.8, 0.2), aspect.ratio = 5/8)


toPlot.stats$label <- as.character(exp_metadata$gene[match(toPlot.stats$KO, exp_metadata$ORF)])
toSave <- toPlot.stats
file_name = paste("kinase_mutant_stats", fun_name, "csv", sep = ".")
file_path = paste(suppl_dir, file_name, sep="/")
write.table(x = toSave, file = file_path, row.names = F, quote = F, eol = "\n", col.names = T, sep = "\t")

s.kinase_boxplot_all$toScale <- T
plots.list <- lappend(plots.list, s.kinase_boxplot_all)



### -- perturbation histogram -------------

toPlot = proteins.FC.f %>% filter(isiMM904 == T)
toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)
all_measured_enzymes <- as.vector((proteins.FC.f %>% filter(isiMM904 ==T ) %>% dplyr::select(ORF) %>% distinct())$ORF)

toPlot <- toPlot %>% 
  group_by(KO) %>%
  summarize(value = sum(sign),
            fraction = value/length(all_measured_enzymes)) %>% 
  arrange(value)

toPlot$label <- as.character(exp_metadata$gene[match(toPlot$KO, exp_metadata$ORF)])
toPlot$label <- factor(toPlot$label, levels = unique(toPlot$label))

p.enzyme_perturbation <- ggplot(toPlot, aes(x = fraction)) +
  geom_histogram(colour = "white", fill = "black", binwidth = 0.05) +
  xlab(paste("Fraction of perturbed metabolic network per kinase mutant")) +
  theme_bw() +
  theme(aspect.ratio = 1)

p.volcano_combined2 <- p.volcano1 + annotation_custom(grob = ggplotGrob(p.enzyme_perturbation), xmin = 0, xmax = 5, ymin = 20 , ymax = 120)

# IRT chromatogram stability ------------------


load("./R/objects/dataset_irt_openswath._clean_.RData")
toPlot <- dataset_irt_openswath %>% 
  dplyr::select(Sequence, filename, RT, ProteinName, m_score) %>% 
  filter(m_score < 0.01, ProteinName == "1/Biognosys_iRT") %>% distinct()

toPlot$sample_name <- sub(pattern = ".mzML", x = basename(toPlot$filename), replacement = "")
selected.samples  = unique(as.character(exp_metadata$sample_name[exp_metadata$type == "Kinase" | exp_metadata$type == "Wild Type" | exp_metadata$type == "Standard Mix"]))

toPlot <- toPlot %>% filter(sample_name %in% selected.samples)
toPlot$aquisition_date <- exp_metadata$aquisition_date[match(toPlot$sample_name, exp_metadata$sample_name)]
toPlot$batch.exp <- exp_metadata$batch.exp[match(toPlot$sample_name, exp_metadata$sample_name)]
toPlot$aquisition_date.str = as.Date(strptime(toPlot$aquisition_date, format="%Y-%m-%d %H:%M:%S"))

tmp.selected <- toPlot %>% 
  group_by(Sequence, batch.exp) %>% 
    summarise(SD = sd(RT, na.rm=T)) %>% 
  group_by(Sequence) %>%
    summarise(MEAN.SD = mean(SD, na.rm=T)) %>%
  ungroup() %>%
    mutate(toRemove = ifelse(MEAN.SD == max(MEAN.SD),1,0))

toPlot <- toPlot %>% filter(Sequence %in% tmp.selected$Sequence[tmp.selected$toRemove == 0])
library(ggthemes)
library(scales)

p.irt <- ggplot(toPlot, aes(x=aquisition_date.str, y=RT/60)) +
  geom_point(aes(colour=Sequence)) +
  geom_vline(data = toPlot[grepl("mix", toPlot$filename, ignore.case = T),], 
             aes(xintercept=as.numeric(aquisition_date.str)), linetype=1, alpha=0.75, colour="darkgrey") +
  xlab("Acquisition date") +
  ylab("Retention time, min") +
  scale_colour_tableau() +
  theme(legend.position="none",
        aspect.ratio = 5/8)

toPlot <- toPlot %>% 
  group_by(filename) %>%
  mutate(iRT = (RT - sort(RT)[1])/(max(RT, na.rm = T) - sort(RT)[1]))

# p.irt.normalized <- ggplot(toPlot, aes(x=aquisition_date.str, y=iRT)) +
#   geom_point(aes(colour=Sequence)) +
#   geom_vline(data = toPlot[grepl("mix", toPlot$filename, ignore.case = T),], 
#              aes(xintercept=as.numeric(aquisition_date.str)), linetype=1, alpha=0.75, colour="darkgrey") +
#   xlab("Acquisition date") +
#   ylab("Retention time, min") +
#   scale_colour_tableau() +
#   theme(legend.position="none",
#         aspect.ratio = 5/8)

forModels <- toPlot %>% 
  group_by(Sequence) %>% 
  mutate(Z_iRT = (iRT - mean(iRT, na.rm = T))/sd(iRT, na.rm = T),
         toRemove  = ifelse(abs(Z_iRT) > 2, 1, 0)) %>%
  filter(toRemove != 1) %>%
  group_by(Sequence) %>%
  summarize(iRT_model = mean(iRT, na.rm=T))

toPlot <- ddply(toPlot, .(filename), 
      .fun = function(x) {
        z <<- x
        #x = toPlot[toPlot$filename == "./data.raw/KL_St_Mix_11.mzML",]
        tmp = data.frame(RT  = x$RT, iRT_model = forModels$iRT_model[match(x$Sequence, forModels$Sequence)])
        
        lm.fit = lm(iRT_model~RT, data = na.omit(tmp))
        tmp.dataset = data.frame(RT = x$RT)
        x$iRT_recalibrated = predict(lm.fit, tmp.dataset)
        return(x)
      })  


p.irt.recalibrated <- ggplot(toPlot, aes(x=aquisition_date.str, y=iRT_recalibrated)) +
  geom_point(aes(colour=Sequence)) +
  geom_vline(data = toPlot[grepl("mix", toPlot$filename, ignore.case = T),], 
             aes(xintercept=as.numeric(aquisition_date.str)), linetype=1, alpha=0.75, colour="darkgrey") +
  xlab("Acquisition date") +
  ylab("Relative retention time") +
  #scale_y_date() +
  scale_colour_tableau() +
  theme_bw() +
  theme(legend.position="none",
        aspect.ratio = 5/8)


gradient_length = 2400
toPlot.stats <- toPlot %>% group_by(Sequence) %>% summarise(median_irt = median(iRT_recalibrated),
                                            sd_rt = sd(RT, na.rm = T),
                                            mean_rt = mean(RT, na.rm = T),
                                            cv_rt = sd_rt/mean_rt,
                                            mean_irt = mean(iRT_recalibrated, na.rm = T),
                                            sd_irt = sd(iRT_recalibrated, na.rm = T),
                                            cv_irt = sd_irt/mean_irt,
                                            scaled_sd = sd_rt/gradient_length)

mean(toPlot.stats$scaled_sd)

p.irt.recalibrated2 <- ggplot() +
  geom_point(data = toPlot, aes(x=aquisition_date.str, y=iRT_recalibrated,colour=Sequence)) +
  geom_point(data = toPlot.stats,
               aes(x= as.Date("2014-08-22", format="%Y-%m-%d"), 
                   y = mean_irt, colour=Sequence)) +
  geom_errorbar(data = toPlot.stats, width=2,
               aes(x = as.Date("2014-08-22", format="%Y-%m-%d"), 
                   ymax = mean_irt + sd_irt,  
                   ymin = mean_irt - sd_irt, colour=Sequence)) +
  geom_text(data = toPlot.stats, 
            aes(x = as.Date("2014-08-22", format="%Y-%m-%d"), 
                y = mean_irt,
                label = percent(scaled_sd),
                colour=Sequence)) +
  geom_vline(data = toPlot[grepl("mix", toPlot$filename, ignore.case = T),], 
             aes(xintercept=as.numeric(aquisition_date.str)), linetype=1, alpha=0.75, colour="darkgrey") +
  xlab("Acquisition date") +
  ylab("Relative retention time") +
  #scale_y_date() +
  scale_colour_tableau() +
  theme_bw() +
  theme(legend.position="none",
        aspect.ratio = 5/8)


## -- transcriptome vs proteome ----------------

# load("./R/objects/proteins.matrix.combat.quant.FC.RData")
load("./R/objects/proteins.matrix.sva.0.5.1.FC.RData")
load("./R//objects/transcriptome.FC._clean_.RData")
load("./R/objects/orf2ko._load_.RData")
load("./R/objects/gene.annotations._load_.RData")

orf2name = droplevels(unique(gene.annotations[,c("V4", "V6")]))
orf2name$V4 = as.character(orf2name$V4)
orf2name$V6 = as.character(orf2name$V6)
orf2name$V6[orf2name$V6 == ""] = orf2name$V4[orf2name$V6 == ""]
names(orf2name) = c("ORF", "gene_name")

transcriptome.FC$KO = NULL
transcriptome.FC = transcriptome.FC[,c("ORF", "logFC", "p.value", "p.value_BH", "KO.ORF")]
names(transcriptome.FC)[length(transcriptome.FC)] = "KO"
transcriptome.FC$isMetabolic = transcriptome.FC$ORF %in% unique(iMM904$gene)
transcriptome.FC$isiMM904 = transcriptome.FC$ORF %in% unique(as.character(iMM904$gene))
transcriptome.FC.f = transcriptome.FC[transcriptome.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]


proteins.FC = proteins.matrix.sva.0.5.1.FC
tr.pr.FC = merge(transcriptome.FC, proteins.FC, by=c("KO", "ORF"), suffixes=c(".tr", ".pr"))

pval_thr = 0.01

tr.pr.cor.ORF = tr.pr.FC %>% group_by(ORF, isMetabolic) %>% 
  summarise(cor = cor(logFC.tr, logFC.pr, method="spearman"))

tr.pr.cor.KO = tr.pr.FC %>% group_by(KO, isMetabolic) %>% 
  summarise(cor = cor(logFC.tr, logFC.pr, method="spearman"))

p.metabolic_correlations <- ggplot(tr.pr.cor.KO, aes(x=cor, fill = isMetabolic)) +  #to Supplementary
  geom_density(alpha=0.5) +
  xlab(expression(paste("Spearman's correlation coefficient, ", rho))) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.8), aspect.ratio = 1) 
  
plots.list <- lappend(plots.list, p.metabolic_correlations)

tr.pr.cor.ORF.f = tr.pr.FC %>% group_by(ORF) %>% filter(p.value.tr < pval_thr) %>%
  summarise(cor = cor(logFC.tr, logFC.pr, method="pearson"),
            n = n())

tr.pr.cor.KO.f = tr.pr.FC %>% group_by(KO) %>% filter(isiMM904) %>%
  summarise(#cor = cor(logFC.tr, logFC.pr, method="pearson"),
            cor = cor.test(logFC.tr, logFC.pr)$estimate,
            p.value = cor.test(logFC.tr, logFC.pr)$p.value,
            n = n())


toPlot = tr.pr.cor.KO.f %>% arrange(cor)
toPlot$KO.name = orf2name$gene_name[match(toPlot$KO, orf2name$ORF)]
#toPlot$KO.name = paste(Hmisc::capitalize(tolower(toPlot$KO.name)), "p", sep="")
toPlot$KO.name = factor(toPlot$KO.name, levels=unique(toPlot$KO.name))


p.tr_vs_pr_cor = ggplot(toPlot, aes(x = KO.name, y = cor,  size=n)) + 
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.5,-0.25,0.25, 0.5), linetype = 3) +  
  ylab(expression(paste("Pearson correlation between\ngene and protein expression changes, ", r))) +
  xlab("Kinase mutant") +
  ylim(-0.8,0.8) +
  #guides(guide = guide_legend(tittle = "Deferentially expressed transcripts")) +
  scale_size(name="Number of common transcripts\nand proteins changed per mutant") +
  #coord_flip() + 
  theme_bw() +
  theme(legend.position = c(0.25, 0.3),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))


p.tr_vs_pr_cor.v = ggplot(toPlot, aes(x = KO.name, y = cor)) + 
  geom_point() +
  geom_point(data=toPlot[toPlot$p.value >0.01,], aes(x=KO.name, y=cor), colour = "grey") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.5,-0.25,0.25, 0.5), linetype = 3) +  
  ylab(expression(paste("Pearson correlation between gene and protein expression log fold-changes, ", r))) +
  xlab("Kinase mutant") +
  ylim(-0.8,0.8) +
  #guides(guide = guide_legend(tittle = "Deferentially expressed transcripts")) +
  scale_size(name="Number of common transcripts\nand proteins changed per mutant") +
  coord_flip() + 
  theme_bw() +
  theme(legend.position = c(0.25, 0.3),
        axis.text.y = element_text(face = "italic"))

plots.list <- lappend(plots.list, p.tr_vs_pr_cor.v)

toPlot = tr.pr.cor.KO.f
toPlot$KO.name = orf2name$gene_name[match(toPlot$KO, orf2name$ORF)]
#toPlot$KO.name = paste(Hmisc::capitalize(tolower(toPlot$KO.name)), "p", sep="")
toPlot$KO.name = factor(toPlot$KO.name, levels=unique(toPlot$KO.name))

p.tr_vs_pr_cor.hist <- ggplot(toPlot, aes(x = cor)) +
  geom_histogram(colour = "white", fill = "black", binwidth = 0.05) +
  geom_vline(aes(xintercept = median(cor)), col="red", linetype = 2) +
  xlab(expression(paste("Pearson correlation between gene and protein expression changes, ", r))) +
  theme_bw() +
  theme(aspect.ratio = 1)
  

# --- saturation plots for supplementary ---------------

pval_thr = 0.01
FC_thr  = log2(1)

transcriptome.FC.f.metabolic <- tbl_df(transcriptome.FC.f) %>% 
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr, isiMM904 == T)

transcriptome.FC.f.stats <- transcriptome.FC.f %>% 
  group_by(KO) %>%
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr) %>% 
  summarize(n_metabolic = sum(isMetabolic == T),
            n_yeast  = sum(isiMM904),
            n_total = n())

combinations <- combn(unique(as.character(transcriptome.FC.f.metabolic$KO)), 2)
combinations.df <- as.data.frame(t(combinations))
KO.genes <- transcriptome.FC.f.metabolic %>% dplyr::select(KO, ORF) %>% distinct()
KO.genes$KO <- as.character(KO.genes$KO)
KO.genes$ORF <- as.character(KO.genes$ORF)

overlaps <- ddply(combinations.df, .(V1, V2),
                  
                  .fun = function(x) {
                    df1 <- KO.genes %>% filter(KO == x$V1)
                    df2 <- KO.genes %>% filter(KO == x$V2)
                    
                    result <- bind_rows(df1, df2) %>% 
                      group_by(ORF) %>%
                      summarise(gene_n = n()) %>% 
                      summarize(intersection = sum(gene_n == 2),
                                union = n(),
                                overlap = intersection/union)
                    
                    
                    return(result)
                  }  )

all_measured_enzymes <- as.vector((transcriptome.FC.f %>% filter(isiMM904 ==T ) %>% dplyr::select(ORF) %>% distinct())$ORF)
transcriptome.FC.f.metabolic$KO.gene <- orf2name$gene_name[match(transcriptome.FC.f.metabolic$KO, orf2name$ORF)]

FC.f.metabolic.stats <-transcriptome.FC.f.metabolic %>% 
  group_by(KO.gene) %>% 
  summarise(n = n(),
            n_pos = sum(logFC > 0)/length(all_measured_enzymes),
            n_neg = sum(logFC < 0)/length(all_measured_enzymes)) %>% 
  ungroup() %>% arrange(n)


overlaps.stats <- overlaps %>% group_by(V1) %>%
  summarize(mean.intersection = mean(intersection, na.rm=T)) %>%
  mutate(gene_name = orf2name$gene_name[match(V1, orf2name$ORF)]) %>%
  left_join(FC.f.metabolic.stats, by = c("gene_name" = "KO.gene"))

toPlot <- overlaps.stats %>% 
  ungroup() %>% 
  mutate(n_fraction = n/length(all_measured_enzymes),
         mean.intersection_fraction = mean.intersection/length(all_measured_enzymes))
library(splines)

p.saturation <- ggplot(toPlot, aes(x = n_fraction, y = mean.intersection_fraction)) +
  geom_point() +
  stat_smooth(method = "lm", formula=y~ns(x,2), se=F) +
  scale_size(range = c(1, 5)) +
  theme_bw() +
  xlab("Fraction of perturbed metabolic network, %") +
  theme(legend.position = c(0.7, 0.5))


transcriptome.FC.f.stats$n_yeast_fraction <- transcriptome.FC.f.stats$n_metabolic/transcriptome.FC.f.stats$n_total
toPlot <- transcriptome.FC.f.stats

p.constant_fraction <- ggplot(toPlot %>% filter(n_total > 0 ), aes(x = n_total, y = n_yeast_fraction)) +
  geom_point() +
  geom_smooth(method = "loess", se=F) +
  xlab("Number of proteins changes per mutant") +
  ylab("Fraction of metabolic enzymes affected by kinase") +
  theme_bw() + theme(aspect.ratio = 5/8)

plots.list <- lappend(plots.list, p.constant_fraction)



# Batch Effects ---------------------
#load("./R/objects/peptides.matrix.RData")
#load("./R/objects/peptides.matrix.quant.combat.RData")

load("./R/objects/peptides.matrix.RData")
load("./R/objects/peptides.matrix.sva.0.5.1.sva_batch_effects.RData")



before = peptides.matrix
after  = peptides.matrix.sva.0.5.1

stopifnot(dim(before) == dim(after))

set.seed(1234)
exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))

pca = prcomp(t(before), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = data.frame(x_var, y_var, type="before")

scores = as.data.frame(pca$x[,1:5])
scores$type = "before"
scores$sample.id = rownames(scores)

pca = prcomp(t(after), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = rbind(annot,data.frame(x_var, y_var, type="after"))

scores = rbind(scores, data.frame(sample.id = rownames(scores), pca$x[,1:5], type = "after"))
scores$batch = factor(exp_metadata$batch.exp.n[match(scores$sample.id, exp_metadata$sample_name)])
scores.mix = scores[grepl(pattern="mix", ignore.case=T, x=rownames(scores)),]
annot$text = paste(annot$x_var, annot$y_var)
scores$type = factor(scores$type, levels=c("before", "after"))


library(ggthemes)

p.batch <- ggplot(scores, aes(x=PC1, y=PC2)) + 
  geom_point(size=2, aes(col=batch))+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(data=scores[grepl(pattern="mix", ignore.case=T, x=scores$sample.id),], 
             aes(x=PC1, y=PC2), size=3, col="black", shape=17) +
  geom_text(data = annot, aes(x=0, y=0, label=text)) +
  facet_wrap(~type, scales="fixed") + 
  #scale_colour_tableau() +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "none")

plots.list <- lappend(plots.list, p.batch)


## ---- cophenetic correlation -----------
library(dendextend)
load("./R/objects/dataPPP_AA.imputed.create_datasets.RData")


left_data <- scale(dataPPP_AA.imputed$proteins.log.quant)
right_data <- scale(dataPPP_AA.imputed$metabolites)

rownames(left_data) <- exp_metadata$gene[match(rownames(left_data), exp_metadata$ORF)]
rownames(right_data) <- exp_metadata$gene[match(rownames(right_data), exp_metadata$ORF)]

d_left <- left_data %>% dist %>% hclust %>% as.dendrogram
d_right <- right_data %>% dist %>% hclust %>% as.dendrogram
d_left_right <- dendlist(proteome = d_left, metabolome = d_right)

file_name = paste("coephenetic_correlation", fun_name, "pdf", sep = ".")
file_path = paste(figures_dir, file_name, sep="/")

pdf(file = file_path, width = 210*0.039 , height = 297*0.039)
  plot(d_left_right, main_left = names(d_left_right)[1], 
     main_right = names(d_left_right)[2], 
     sub = paste("Cophenetic correlation =", format(cor_cophenetic(d_left_right), digits = 2, scientific = T)), 
     cex_main_left = 0.85, cex_main_right = 0.85,
     cex_sub = 0.75)
dev.off()

p = recordPlot()
plots.list = lappend(plots.list, p)


### ---- sentinels ------
#load("./R/objects/sentinels.proteins.matrix.quant.combat.FC.RData")
#load("./R/objects/sentinels.proteins.matrix.quant.combat_pseudo.FC.PSEUDO.RData")
#load("./R/objects/sentinels.proteins.matrix.quant.combat.RData")
load("./R/objects/sentinels.table._clean_.RData")

#load("./R/objects/proteins.matrix.combat.quant.FC.RData")
#load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/sentinels.proteins.matrix.sva.0.5.1.RData")
load("./R/objects/sentinels.proteins.matrix.sva.0.5.1.FC.WT.RData")
load("./R/objects/iMM904._load_.RData")

protein.matrix = sentinels.proteins.matrix.sva.0.5.1
proteins.FC = sentinels.proteins.matrix.sva.0.5.1.FC

proteins.FC$gene_name <- orf2name$gene_name[match(proteins.FC$ORF, orf2name$ORF)]
reference = unique(as.character(proteins.FC$reference))
proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]

sentinels.df <- dcast(proteins.FC.f, formula = "KO~ORF", value.var = "logFC")
sentinels.matrix <- as.matrix(sentinels.df[,-c(1,2)])
rownames(sentinels.matrix) <- sentinels.df[,1]


pval_thr = 0.01
set.seed(123)
FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)


sentinels.table$short_info <- sub(pattern = "Marker for ", replacement = "", x = sentinels.table$Sentinel)
sentinels.table$short_info <- trimws(sentinels.table$short_info)
sentinels.table$ORF.ID <- trim(sentinels.table$ORF.ID)

proteins.FC.f$KO.gene <- exp_metadata$gene[match(proteins.FC.f$KO, exp_metadata$ORF)]

toPlot <- proteins.FC.f %>% 
  filter(ORF != "iRT", abs(logFC) > FC_thr, p.value_BH < pval_thr) %>%
  group_by(KO) %>%
  mutate(sign = n()) %>%
  group_by(ORF) %>%
  mutate(popularity = n())
#%>%
#  filter(popularity, sign >= 15)

sentinels_desc_collapsed = sentinels.table %>% dplyr::select(ORF.ID, short_info) %>% distinct() %>%
  filter(ORF.ID %in% toPlot$ORF) %>%
  group_by(short_info) %>%
  mutate(info_n = n()) %>%
  group_by(ORF.ID) %>%
  summarise(description_full = paste(sort(short_info),  collapse = "|"),
            description = ifelse(length(short_info) > 2, "multiple", paste(sort(short_info),  collapse = "|")),
            popular_info = paste(max(info_n), short_info[which.max(info_n)], collapse = ":" )) %>% 
  ungroup %>% 
  arrange(popular_info)

sentinels_desc_single = sentinels.table %>% dplyr::select(ORF.ID, short_info) %>%  filter(ORF.ID %in% toPlot$ORF) %>% distinct()


lb = -2.5
ub = 2.5

my_breaks <- seq(lb, ub, 0.5)
my_breaks[5] <- FC_thr
my_breaks[7] <- -FC_thr
my_levels = levels(cut(toPlot$logFC, breaks = my_breaks))

my_colours1 <- brewer.pal(name = "Reds", n = 4)
my_colours2 <- rev(brewer.pal(name = "Blues", n = 4))
my_colours = c(my_colours2, my_colours1)

toPlot$my_fill <- as.character(cut(toPlot$logFC, breaks = my_breaks))
toPlot$my_fill[which(is.na(toPlot$my_fill))] <- ifelse(toPlot$logFC[which(is.na(toPlot$my_fill))] <= lb, lb, ub)
toPlot$my_fill <- factor(toPlot$my_fill, levels = my_levels)

toPlot$info <- sentinels_desc_collapsed$description[match(toPlot$ORF, sentinels_desc_collapsed$ORF.ID)]
toPlot$sentinel <- paste(toPlot$info, Hmisc::capitalize(tolower(toPlot$gene_name)), sep = "::")

toPlot$popular_info <- sentinels_desc_collapsed$popular_info[match(toPlot$ORF, sentinels_desc_collapsed$ORF.ID)]
toPlot$sentinel_popular <- paste(toPlot$popular_info, paste(Hmisc::capitalize(tolower(toPlot$gene_name)), "p", sep=""), sep = "::")

# sorted = sort(unique(toPlot$sentinel))
# sorted <- c(sorted[grep("multiple", x = sorted, invert = T)], sorted[grep("multiple", x = sorted)])
#toPlot$sentinel <- factor(toPlot$sentinel, levels = rev(sorted))
sorted <- sort(unique(toPlot$sentinel_popular))
#sorted <- sorted[c(grep(12, sorted,invert =T), grep(12, sorted))]
toPlot$sentinel_popular <- factor(toPlot$sentinel_popular, levels = sorted)

p.heatmap_sentinels <- ggplot(toPlot) +  
  geom_tile(aes(x = KO.gene, y = sentinel_popular, fill = my_fill ), colour="grey") +
  #   scale_fill_gradient2(low="#1F78B4",high="#E31A1C",mid ="white",
  #                        breaks = seq(-0.75, 0.75, 0.25),
  #                        midpoint=0)  +
  scale_fill_manual(values = my_colours, 
                    name = "log2(mutant/WT)",
                    labels = c("<-2", "-1.5", "-1", "-0.41", "0.41", "1", "1.5", ">2")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="italic"),  
        legend.position = "left", aspect.ratio = 1, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x="Kinase mutant", y = "Marker for biological process")

plots.list <- lappend(plots.list, p.heatmap_sentinels)

p.heatmap_sentinels_slide <- ggplot(toPlot) +  
  geom_tile(aes(x = KO.gene, y = sentinel_popular, fill = my_fill ), colour="grey") +
  #   scale_fill_gradient2(low="#1F78B4",high="#E31A1C",mid ="white",
  #                        breaks = seq(-0.75, 0.75, 0.25),
  #                        midpoint=0)  +
  scale_fill_manual(values = my_colours, 
                    name = "log2(mutant/WT)",
                    labels = c("<-2", "-1.5", "-1", "-0.41", "0.41", "1", "1.5", ">2")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="italic"),  
        legend.position = "left", aspect.ratio = 1, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x="Kinase mutant", y = "Marker for biological process")

p.heatmap_sentinels_slide

file_name = "sentinels_slide.pdf"
file_path = paste(figures_dir, file_name, sep="/")

save_plot(plot = p.heatmap_sentinels_slide, filename = file_path, base_height = 8.27*1.5, base_width = 11.69*1.5 )

### ---- sentinel coverage ---- ####

load("./R/objects/sentinelsSRM._clean_.RData")

sum(unique(sentinelsSRM$ORF) %in% rownames(sentinels.proteins.matrix.sva.0.5.1))/length(unique(sentinelsSRM$ORF))

#sum((sentinels.table %>% filter(Sentinel.Grade == "A") %>% select(ORF.ID) %>% distinct())$ORF.ID %in% rownames(proteins.matrix.combat.quant))


# --- constant fraction & saturation -----
#load("./R/objects/proteins.matrix.combat.quant.FC.RData")
#load("./R/objects/proteins.matrix.combat.quant.RData")

load("./R/objects/proteins.matrix.sva.0.5.1.FC.RData")
load("./R/objects/proteins.matrix.sva.0.5.1.RData")

load("./R/objects/iMM904._load_.RData")

EC.genes = gene.annotations[gene.annotations$V3 == "EC number",]
protein.matrix <- proteins.matrix.sva.0.5.1
proteins.FC <- proteins.matrix.sva.0.5.1.FC

pval_thr = 0.01
set.seed(123)
FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)

proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]
#proteins.FC.f$isMetabolic = proteins.FC.f$ORF %in% unique(KEGG.pathways.f$ORF)
proteins.FC.f$isEnzyme = proteins.FC.f$ORF %in% unique(EC.genes$V4)
proteins.FC.f$isiMM904 = proteins.FC.f$ORF %in% unique(as.character(iMM904$gene))

proteins.FC.f.stats <- proteins.FC.f %>% 
  group_by(KO) %>%
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr) %>% 
  summarize(#n_metabolic = sum(isMetabolic == T),
            n_yeast  = sum(isiMM904),
            n_total = n())
proteins.FC.f.stats$gene_name <- orf2name$gene_name[match(proteins.FC.f.stats$KO, orf2name$ORF)]

filter(proteins.FC.f.stats, gene_name %in% c("FMP48", "NPR1", "NNK1", "VPS15"))

proteins.FC.f.stats$n_yeast_fraction <- proteins.FC.f.stats$n_yeast/proteins.FC.f.stats$n_total
toPlot <- proteins.FC.f.stats 
p.constant_fraction <- ggplot(toPlot, aes(x = n_total, y = n_yeast_fraction)) +
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Number of proteins changes per mutant") +
  ylab("Fraction of metabolic enzymes affected by kinase") +
  theme_bw()

all_proteins <- as.character(unique(proteins.FC.f$ORF))
all_measured_enzymes <- as.vector((proteins.FC.f %>% filter(isiMM904 ==T ) %>% dplyr::select(ORF) %>% distinct())$ORF)

proteins.FC.f.metabolic <- tbl_df(proteins.FC.f) %>% 
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr, isiMM904 == T)

proteins.FC.f.metabolic$KO.gene <- orf2name$gene_name[match(proteins.FC.f.metabolic$KO, orf2name$ORF)]
stopifnot(!any(is.na(proteins.FC.f.metabolic$KO.gene)))

FC.f.metabolic.stats <- proteins.FC.f.metabolic %>% 
  group_by(KO.gene) %>% 
  summarise(n = n(),
            n_pos = sum(logFC > 0)/length(all_measured_enzymes),
            n_neg = sum(logFC < 0)/length(all_measured_enzymes)) %>% 
  ungroup() %>% arrange(n)

combinations <- combn(unique(proteins.FC.f.metabolic$KO), 2)
combinations.df <- as.data.frame(t(combinations))

KO.genes <- proteins.FC.f.metabolic %>% dplyr::select(KO, ORF) %>% distinct()


overlaps <- ddply(combinations.df, .(V1, V2),
                  
                  .fun = function(x) {
                    df1 <- KO.genes %>% filter(KO == x$V1)
                    df2 <- KO.genes %>% filter(KO == x$V2)
                    
                    result <- bind_rows(df1, df2) %>% 
                      group_by(ORF) %>%
                      summarise(gene_n = n()) %>% 
                      summarize(intersection = sum(gene_n == 2),
                                union = n(),
                                overlap = intersection/union)
                    
                    
                    return(result)
                  }  )

overlaps.stats <- overlaps %>% group_by(V1) %>%
  summarize(mean.intersection = mean(intersection, na.rm=T)) %>%
  mutate(gene_name = orf2name$gene_name[match(V1, orf2name$ORF)]) %>%
  left_join(FC.f.metabolic.stats, by = c("gene_name" = "KO.gene"))

overlaps.stats$degree <- igraph::degree(GRAPH)[match(overlaps.stats$V1, names(igraph::degree(GRAPH)))]

toPlot <- overlaps.stats %>% 
  ungroup() %>% 
  mutate(n_fraction = n/length(all_measured_enzymes),
         mean.intersection_fraction = mean.intersection/length(all_measured_enzymes))
library(splines)

p.saturation <- ggplot(toPlot, aes(x = n_fraction, y = mean.intersection_fraction)) +
  geom_point(aes(size =degree)) +
  stat_smooth(method = "lm", formula=y~ns(x,2), se=F) +
  scale_size(range = c(1, 5)) +
  theme_bw() +
  xlab("Fraction of perturbed metabolic network, %") +
  theme(legend.position = c(0.7, 0.5))

# -- PPI and other networks ----- 
load("./R/objects/yeastract._load_.RData")
load("./R/objects/GO_slim.raw._load_.RData")
load("./R/objects/STRING._load_.RData")
load("./R/objects/yeast.ppi._load_.RData")

string.exp = STRING %>% filter(experimental > 900)
string.all = STRING %>% filter(combined_score > 900)

G.string.exp = graph.data.frame(droplevels(string.exp %>% dplyr::select(ORF1, ORF2)), directed=F)
G.string.all = graph.data.frame(droplevels(string.all %>% dplyr::select(ORF1, ORF2)), directed=F)
GRAPH = G.string.exp

# ---- pathways and TFs -------------

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
proteins.FC.f$isiMM904 = proteins.FC.f$ORF %in% unique(as.character(iMM904$gene))


proteins.FC.f.stats = proteins.FC.f %>% filter(isiMM904, p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% 
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


min.paths.comb.stats = tbl_df(proteins.FC.f.stats.long %>% filter(min.paths.comb != 0)) %>% 
  group_by(min.paths.comb, variable) %>% 
  summarize( FC.mean = mean(value, na.rm=T)/mean(control.min.paths.comb[control.min.paths.comb[,5] == variable, "value"], na.rm=T),
             FC.median = median(value, na.rm=T)/median(control.min.paths.comb[control.min.paths.comb[,5] == variable, "value"], na.rm=T),
             p.value = wilcox.test(value,control.min.paths.comb[control.min.paths.comb[,5] == variable, "value"])$'p.value')


min.paths.stats = proteins.FC.f.stats.long %>% filter(min.paths != 1) %>% 
  group_by(min.paths, variable) %>% 
  summarize( FC.mean = mean(value, na.rm=T)/mean(control.min.paths[control.min.paths[,5] == variable, "value"], na.rm=T),
             FC.median = median(value, na.rm=T)/median(control.min.paths[control.min.paths[,5] == variable, "value"], na.rm=T),
             p.value = wilcox.test(value, control.min.paths[control.min.paths[,5] == variable, "value"])$'p.value')


min.paths.stats$symbol = ""
min.paths.stats$symbol = ifelse(min.paths.stats$p.value<0.05, "*", "")
#min.paths.stats$symbol = ifelse(min.paths.stats$p.value<0.001, "**", "")


min.paths.comb.stats$symbol = ""
min.paths.comb.stats$symbol = ifelse(min.paths.comb.stats$p.value<0.05, "*", "")
#min.paths.comb.stats$symbol = ifelse(min.paths.comb.stats$p.value<0.001, "**", "")

toPlot <- proteins.FC.f.stats.long %>% filter(min.paths.comb != "none", variable == "changes")
toText <- min.paths.comb.stats  %>% filter(min.paths.comb != "none", variable == "changes")

p.dist_changes = ggplot(toPlot , aes(x=factor(min.paths.comb), y=value)) +
  stat_boxplot(geom ='errorbar', width = 0.33) +
  geom_boxplot(width = 0.5) + 
  geom_text(data=toText, aes(x = factor(min.paths.comb) , y=1, label=symbol), col="red", size=10)+
  scale_x_discrete(labels=c("First", ">1")) +
  xlab("Length of the shortest path in PPI network from kinase to transcription factor") +
  ylab("Number of perturbed proteins") +
  theme_bw()


toPlot <- proteins.FC.f.stats.long %>% filter(min.paths.comb != "none", variable == "betweenness")
toText <- min.paths.comb.stats  %>% filter(min.paths.comb != "none", variable == "betweenness")
library(scales)
p.dist_between <- ggplot(toPlot , aes(x=factor(min.paths.comb), y=value)) +
  stat_boxplot(geom ='errorbar', width = 0.33) +
  geom_boxplot(width = 0.5) + 
  geom_text(data=toText, 
            aes(x = factor(min.paths.comb) , y=1, label=symbol), col="red", size=10)+
  scale_x_discrete(labels=c("First", ">1")) +
  scale_y_continuous(name="Number of shortest paths through kinase, node betweeness", labels = scientific) +
  xlab("") +
  theme_bw()

toPlot <- proteins.FC.f.stats.long %>% filter(min.paths.comb != "none", variable == "degree")
toText <- min.paths.comb.stats  %>% filter(min.paths.comb != "none", variable == "degree")
p.dist_degree <- ggplot(toPlot , aes(x=factor(min.paths.comb), y=value)) +
  stat_boxplot(geom ='errorbar', width = 0.33) +
  geom_boxplot(width = 0.5) + 
  geom_text(data=toText, 
            aes(x = factor(min.paths.comb) , y=1, label=symbol), col="red", size=10)+
  scale_x_discrete(labels=c("First", ">1")) +
  scale_y_continuous(name="Number of neigbours of kinasse in PPI, node degree") +
  xlab("") +
  theme_bw()


toPlot <- proteins.FC.f.stats.long %>% 
  filter(variable == "betweenness" | variable == "degree"| variable == "changes") %>%
  dcast(formula = "KO+gene_name~variable", value.var = "value") %>%
  melt(id.vars = c("KO", "gene_name", "changes"))

toPlot.stats <- toPlot %>% group_by(variable) %>%
  summarise(x = max(changes, na.rm = T) - 0.2*max(changes,  na.rm = T),
            y = max(value, na.rm = T) - 0.2*max(value, na.rm = T),
            cor = cor.test(changes, value)$estimate,
            p.value = cor.test(changes, value)$p.value)
toPlot.stats$variable = factor(toPlot.stats$variable, levels = c("degree", "betweenness"))
toPlot$variable = factor(toPlot$variable, levels = c("degree", "betweenness"))

p.network_vs_changes <- ggplot(toPlot, aes(x = changes, y = value)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE) +
  geom_text(data = toPlot.stats, 
            aes(x = x,  y = y,
                label = paste("r = ", round(cor, 2), "\n",
                              "p-value = ", round(p.value, 2), sep = ""))) +
  facet_wrap(~variable, scales = "free_y") +
  xlab("Number of differentially expressed enzymes per kinase mutant") +
  theme_bw() + 
  theme(aspect.ratio = 1) 

plots.list <- lappend(plots.list, p.network_vs_changes)

## -- growth analysis correlation of growth with gene targets ----

combinedData_raw <- read.csv("./data/2016-08-12/combinedData_raw.csv")
growth.data <- combinedData_raw %>% dplyr::select(X, max.mass, lag, max.slope, integral, Plateau.best, max.slope.time, dM.dT.best)

#Manually exclude outliers from Stephan's  ./data/2016-08-12/dataPrep.py 
#negativeLag = combinedData.index[combinedData['lag']<0].tolist() #Result of this line pasted below
didntGrow = c('YDR127W', 'YHR183W', 'YLR362W', 'YPR060C', 'YPR074C', 'YDR127W', 'YNL241C', 'YPR060C', 'YPR074C', 'YGL148W', 'YGL148W', 'YGR262C')
veryDiffReplicates = c('YOL045W', 'YLR362W', 'YOL045W', 'YLR362W')
growth.data.f <- growth.data %>% filter(!(X %in% didntGrow) | !(X %in% veryDiffReplicates), lag >0)

load("./R/objects/iMM904._load_.RData")

protein.matrix <- proteins.matrix.sva.0.5.1
proteins.FC <- proteins.matrix.sva.0.5.1.FC

reference <- unique(as.character(proteins.FC$reference))
pval_thr = 0.01
set.seed(123)
FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)

proteins.FC.f <- proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]
proteins.FC.f$isiMM904 <- proteins.FC.f$ORF %in% unique(as.character(iMM904$gene))
proteins.FC.f$sign  <- ifelse( abs(proteins.FC.f$logFC) > FC_thr &  proteins.FC.f$p.value_BH < pval_thr, 1, 0 )

proteins.FC.f.stats <- proteins.FC.f %>% 
  filter(sign == 1, isiMM904) %>% 
  group_by(KO) %>%
  summarise(n = n())

toPlot <- inner_join(proteins.FC.f.stats, growth.data.f, by = c("KO" = "X"))
toPlot.stats <- data.frame( cor = cor.test(toPlot$n, toPlot$max.slope)$estimate,
                            p.value = cor.test(toPlot$n, toPlot$max.slope)$p.value)

p.growth <- ggplot(toPlot, aes(x = n, y = max.slope )) +
  geom_point() +
  geom_text(data = toPlot.stats, size = 5,
            aes(x = 40, y = 0.35, 
                label = paste("r = ", round(toPlot.stats$cor, 2), "\n",
                              "p-value = ", round(toPlot.stats$p.value, 2), sep = ""))) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Number of differentially expressed enzymes per kinase mutant comparing to wild type strain") +
  ylab("Exponential growth rate, slope of growth curve") +
  theme_bw() +
  theme(aspect.ratio = 1)

plots.list <- lappend(plots.list, p.growth)

## ---- stats of kinase perturbation, quantified proteins etc ----

load("./R/objects/proteins.matrix.sva.0.5.1.RData")
load("./R/objects/proteins.matrix.sva.0.5.1.FC.RData")
load("./R/objects/iMM904._load_.RData")

protein.matrix = proteins.matrix.sva.0.5.1
proteins.FC = proteins.matrix.sva.0.5.1.FC

pval_thr = 0.01
set.seed(123)
FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)


proteins.FC.f = droplevels(proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]) #only kinases
proteins.FC.f$sign = ifelse(abs(proteins.FC.f$logFC) >= FC_thr & proteins.FC.f$p.value_BH < pval_thr, 1,0)
proteins.FC.f$isiMM904 = proteins.FC.f$ORF %in% unique(as.character(iMM904$gene))


proteins.FC.f.stats.all = proteins.FC.f %>% filter(p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% 
  group_by(KO) %>% dplyr::summarise(changes = n()) %>%
  ungroup() %>% summarise(min = min(changes),
                          max = max(changes),
                          median = median(changes),
                          average = mean(changes))


proteins.FC.f.stats.metabolic = proteins.FC.f %>% filter(isiMM904, p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% 
  group_by(KO) %>% dplyr::summarise(changes = n()) %>%
  ungroup() %>% summarise(min = min(changes),
                          max = max(changes),
                          median = median(changes),
                          average = mean(changes))

stats_table <- data.frame(stats_name = character(),
                          value = character(), 
                          rel_value = character(),
                          comment=character(),
                          stringsAsFactors=FALSE)

#total expressed proteins in yeast 4,517 from 
# Sina Ghaemmaghami et al "Global analysis of protein expression in yeast" Nature 425, 737-741 (16 October 2003) | doi:10.1038/nature02046;
total_proteins = 4517
stats_tmp <- data.frame(stats = "yeast_expressed_ORFs", 
                        value = total_proteins, 
                       comment = "Total number of yeast expressed ORFs from Ghaemmaghami et. al Nature 2003")
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)


#cummulative proteins affected by kinase deletion
stats_tmp <- data.frame(stats = "cummulatite_all", 
                        value = nrow(proteins.FC.f %>% filter(p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% dplyr::select(ORF) %>% distinct()), 
                        comment = "Cummulative changes of all proteins")
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)

#cummulative enzymes proteins affected by kinase deletion
stats_tmp <- data.frame(stats = "cummulatite_enzymes", 
                        value = nrow(proteins.FC.f %>% filter(isiMM904, p.value_BH < pval_thr, abs(logFC) > FC_thr) %>% dplyr::select(ORF) %>% distinct()), 
                        comment = "Cummulative changes of enzymes")
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)

#total quantified proteins in every samples
stats_tmp <- data.frame(stats = "across_samples_all", 
                        value = nrow(proteins.FC.f %>% dplyr::select(ORF) %>% distinct()), 
                        comment = "Total quantified proteins in every sample")
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)

#total quantified enzymes in every samples
stats_tmp <- data.frame(stats = "across_samples_enzymes", 
                        value = nrow(proteins.FC.f %>% filter(isiMM904) %>% dplyr::select(ORF) %>% distinct()), 
                        comment = "Total quantified enzymes in every sample")
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)



#perturbed enzymes per kinase mutant sumaries
tmp <- t(proteins.FC.f.stats.metabolic)
stats_tmp <- data.frame(stats = paste0("enzymes_", rownames(tmp)),
                        value = tmp[,1],
                        comment = paste0("Perturbed enzymes ", rownames(tmp)))
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)

#perturbed all proteins per kinase mutant sumaries
tmp <- t(proteins.FC.f.stats.all)
stats_tmp <- data.frame(stats = paste0("all_", rownames(tmp)),
                        value = tmp[,1],
                        comment = paste0("Perturbed all proteins ", rownames(tmp)))
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)



# adding data about all samples
load("./R/objects/protein_annotations_trypsin._clean_.RData")
load("./R/objects/fragments.data.RData")

protein_annotations.unique <- tbl_df(protein_annotations) %>% 
  ungroup() %>% 
  dplyr::select(strippedSequence, SystName) %>% 
  distinct()


protein_annotations.all <- protein_annotations %>% 
  dplyr::select(strippedSequence, SystName) 

dataset.raw.f.good <- fragments.data %>% left_join(protein_annotations.all, by = c("EG.StrippedSequence" = "strippedSequence"))
dataset.raw.f.good <- dataset.raw.f.good %>% rename(ProteinName = SystName)

fdr_thr1 = 0.01

total_samples = length(unique(dataset.raw.f.good$R.Label))

selected <- unique(as.character(exp_metadata$sample_name[exp_metadata$type == "Kinase" | exp_metadata$type == "Wild Type" | exp_metadata$type == "Standard Mix"]))


proteins.detected <- dataset.raw.f.good %>% filter(R.Label %in% selected) %>%
  ungroup() %>%
  select(ProteinName, R.Label, EG.Qvalue) %>%
  gather(stats, value, -R.Label, -EG.Qvalue) %>%
  group_by(stats, value, R.Label) %>% 
    summarise(min.EG.Qvalue = min(EG.Qvalue, na.rm=T)) %>%
    mutate(isiMM904 = value %in% unique(as.character(iMM904$gene)))

proteins.detected.sumaries <- proteins.detected %>%
  ungroup() %>%
  group_by(R.Label) %>% 
    summarize(n = sum(min.EG.Qvalue < fdr_thr1),
              n_enzymes = sum(min.EG.Qvalue < fdr_thr1 & isiMM904)) %>%
  ungroup() %>%
    summarise(min_all = min(n),
              max_all = max(n),
              median_all = median(n),
              mean_all = mean(n),
              min_enzymes = min(n_enzymes),
              max_enzymes = max(n_enzymes),
              median_enzymes = median(n_enzymes),
              mean_enzymes = mean(n_enzymes))



#Total detected/quantified proteins in whole dataset
stats_tmp <- data.frame(stats = "total_detected_all", 
                        value = nrow(proteins.detected %>% ungroup() %>% filter(min.EG.Qvalue < fdr_thr1) %>%  dplyr::select(value) %>% distinct()) , 
                        comment = paste0("total detected all proteins FDR < ", fdr_thr1))
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)

#Total detected/quantified enzymes in whole dataset
stats_tmp <- data.frame(stats = "total_detected_enzymes", 
                        value = nrow(proteins.detected %>% ungroup() %>% filter(isiMM904, min.EG.Qvalue < fdr_thr1) %>%  dplyr::select(value) %>% distinct()) , 
                        comment = paste0("total detected enzymes FDR < ", fdr_thr1))
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)

#Total detected/quantified enzymes in whole dataset summaries
tmp <- t(proteins.detected.sumaries)
stats_tmp <- data.frame(stats =  rownames(tmp),
                        value = tmp[,1],
                        comment = paste0("Detected in whole dataset ", rownames(tmp)))
stats_tmp$rel_value <- with(stats_tmp, value/total_proteins)
stats_table <- rbind(stats_table, stats_tmp)

library("gridExtra")
p <- tableGrob(stats_table)
p$landscape = T
plots.list = lappend(plots.list, p)


## ----- figure version 1 --------


plot_figure_v1 <- function () {
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(90, 60)))
  
  grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.irt, vp = viewport(layout.pos.row = 1:20, layout.pos.col = 1:30)) #experiment timeline with RT stability
  grid.text("B", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 31),gp=gpar(fontsize=20, col="black"))
  print(p.volcano_combined, vp = viewport(layout.pos.row = 1:20, layout.pos.col = 31:60)) #Volcano plot
  grid.text("C", just=c("left", "centre"), vp = viewport(layout.pos.row = 21, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.coverage, vp = viewport(layout.pos.row = 21:40, layout.pos.col = 1:25))  #pathway coverage detailed
  print(p.pack_man, vp = viewport(layout.pos.row = 21:30, layout.pos.col = 26:35)) #pathway pie chart
  print(p.mean_cov, vp = viewport(layout.pos.row = 31:40, layout.pos.col = 26:35)) #pathway average
  grid.text("D", just=c("left", "centre"), vp = viewport(layout.pos.row = 21, layout.pos.col = 36),gp=gpar(fontsize=20, col="black"))
  #print(p.CV, vp = viewport(layout.pos.row = 21:30, layout.pos.col = 36:60)) #pathway average
  #grid.text("E", just=c("left", "centre"), vp = viewport(layout.pos.row = 21, layout.pos.col = 36),gp=gpar(fontsize=20, col="black"))
  print(p.batch, vp = viewport(layout.pos.row = 21:40, layout.pos.col = 36:60)) #pathway average
}

file_name = "Figure1_v01_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
  plot_figure_v1() 
dev.off()

file_name = "Figure1_v01_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
  plot_figure_v1() 
dev.off()

# --- figure version 2 ---- 

plot_figure_v2 <- function() {
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(135, 100)))
  
  grid.text("a", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.irt.recalibrated, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 1:40)) #experiment timeline with RT stability
  grid.text("b", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 41),gp=gpar(fontsize=20, col="black"))
  print(p.coverage , vp = viewport(layout.pos.row = 1:30, layout.pos.col = 41:80))  #pathway coverage detailed
  print(p.pack_man, vp = viewport(layout.pos.row = 1:15, layout.pos.col = 81:100)) #pathway pie chart
  print(p.mean_cov , vp = viewport(layout.pos.row = 16:30, layout.pos.col = 81:100)) #pathway average
  grid.text("c", just=c("left", "centre"), vp = viewport(layout.pos.row = 31, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.volcano_combined, vp = viewport(layout.pos.row = 31:60, layout.pos.col = 1:40)) #Volcano plot
  grid.text("d", just=c("left", "centre"), vp = viewport(layout.pos.row = 31, layout.pos.col = 42),gp=gpar(fontsize=20, col="black"))
  print(p.tr_vs_pr_cor, vp = viewport(layout.pos.row = 31:60, layout.pos.col = 45:100)) #Volcano plot
  
}

file_name = "Figure1_v02_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")
pdf(file_path, height=247/25.4*2, width=183/25.4*2)
plot_figure_v2()
dev.off()

file_name = "Figure1_v02_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
plot_figure_v2() 
dev.off()

plot_figure_v3 <- function() {
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(135, 100)))
  
  grid.text("a", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.irt.recalibrated, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 1:40)) #experiment timeline with RT stability
  grid.text("b", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 41),gp=gpar(fontsize=20, col="black"))
  print(p.heatmap_sentinels , vp = viewport(layout.pos.row = 1:60, layout.pos.col = 20:100))  #sentinels
  
  
  grid.text("c", just=c("left", "centre"), vp = viewport(layout.pos.row = 31, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.volcano_rect, vp = viewport(layout.pos.row = 31:45, layout.pos.col = 1:15)) #Volcano plot
  print(p.enzyme_perturbation, vp = viewport(layout.pos.row = 31:45, layout.pos.col = 16:30))
  print(p.constant_fraction, vp = viewport(layout.pos.row = 45:60, layout.pos.col = 1:40))
  print(p.coverage , vp = viewport(layout.pos.row = 61:90, layout.pos.col = 1:27))  #pathway coverage detailed
  print(p.pack_man, vp = viewport(layout.pos.row = 61:75, layout.pos.col = 27:40)) #pathway pie chart
  print(p.mean_cov , vp = viewport(layout.pos.row = 76:90, layout.pos.col = 27:40)) #pathway average
  
  print(p.dist_changes, vp = viewport(layout.pos.row = 96:115, layout.pos.col = 41:51))
  grid.text("g", just=c("left", "centre"), vp = viewport(layout.pos.row = 100, layout.pos.col = 67),gp=gpar(fontsize=20, col="black"))
  print(p.dist_degree, vp = viewport(layout.pos.row = 96:115, layout.pos.col = 52:63))
  print(p.dist_between, vp = viewport(layout.pos.row = 96:115, layout.pos.col = 64:75))
  
  grid.text("d", just=c("left", "centre"), vp = viewport(layout.pos.row = 31, layout.pos.col = 42),gp=gpar(fontsize=20, col="black"))
  #print(p.tr_vs_pr_cor, vp = viewport(layout.pos.row = 31:60, layout.pos.col = 45:100)) #Volcano plot
  
}


file_name = "Figure1_v03_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")
pdf(file_path, height=247/25.4*2, width=183/25.4*2)
plot_figure_v3()
dev.off()


file_name = "Figure1_v03_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
plot_figure_v3() 
dev.off()


plot_figure_v4<- function() {
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(135, 100)))
  
  grid.text("a", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.irt.recalibrated2, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 1:40)) #experiment timeline with RT stability

  grid.text("c", just=c("left", "centre"), vp = viewport(layout.pos.row = 31, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.volcano_rect, vp = viewport(layout.pos.row = 31:45, layout.pos.col = 1:15)) #Volcano plot
  print(p.enzyme_perturbation, vp = viewport(layout.pos.row = 31:45, layout.pos.col = 16:30))
  print(p.constant_fraction, vp = viewport(layout.pos.row = 45:60, layout.pos.col = 1:40))
  print(p.coverage , vp = viewport(layout.pos.row = 61:90, layout.pos.col = 1:27))  #pathway coverage detailed
  print(p.pack_man, vp = viewport(layout.pos.row = 61:75, layout.pos.col = 27:40)) #pathway pie chart
  print(p.mean_cov , vp = viewport(layout.pos.row = 76:90, layout.pos.col = 27:40)) #pathway average
  
  print(p.dist_changes, vp = viewport(layout.pos.row = 96:115, layout.pos.col = 41:51))
  grid.text("g", just=c("left", "centre"), vp = viewport(layout.pos.row = 100, layout.pos.col = 67),gp=gpar(fontsize=20, col="black"))
  print(p.dist_degree, vp = viewport(layout.pos.row = 96:115, layout.pos.col = 52:63))
  print(p.dist_between, vp = viewport(layout.pos.row = 96:115, layout.pos.col = 64:75))
  print(p.tr_vs_pr_cor.v, vp = viewport(layout.pos.row = 1:60, layout.pos.col = 41:80)) #Volcano plot
  print(p.volcano_inset, vp = viewport(layout.pos.row = 96:115, layout.pos.col = 1:30)) #Volcano plot
  print(p.tr_vs_pr_cor.hist, vp = viewport(layout.pos.row = 116:135, layout.pos.col = 1:30)) #Volcano plot
  
  
}


file_name = "Figure1_v04_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")
pdf(file_path, height=247/25.4*2, width=183/25.4*2)
plot_figure_v4()
dev.off()


file_name = "Figure1_v04_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
plot_figure_v4() 
dev.off()






#----- Supplementary figures ----- 

file_name = paste("supplementary", fun_name, sep = ".")
file_path = paste(figures_dir, file_name, sep="/")

lapply(seq_along(plots.list) , 
       function(x) {
        tryCatch({
          p <- plots.list[[x]]
           
          width = 210 
          height = 297
          
          scale = 1
          if (length(p$toScale) != 0 && p$toScale == T){
             scale = 2
          }
          
          if (length(p$landscape) != 0 && p$landscape == T){
            width = 297
            height = 210
          }
          
          if(any(grep("^g",class(p)))) {
             ggplot2::ggsave(filename = paste(file_path, x , "pdf", sep = "."), device = NULL,
                             plot = p, width = width , height = height, units = "mm", scale = scale)
            ggplot2::ggsave(filename = paste(file_path, x , "png", sep = "."), device = NULL, dpi = 150,
                            plot = p, width = width , height = height, units = "mm", scale = scale)  
          } else {
              pdf(file = paste(file_path, x , "pdf", sep = "."), width = width*0.039 , height = height*0.039)
              par(cex = 0.8)
              print(p)
              
              png(file = paste(file_path, x , "pdf", sep = "."), width = width*0.039 , height = height*0.039, res = 150)
              par(cex = 0.8)
              print(p)
            dev.off()
          }
         }, error = function(e) {
           message(paste("Plot", x, "sucks!" ))
           return(NULL)
         }, finally = {
           message(paste("processed plot", x))
         })
       })
