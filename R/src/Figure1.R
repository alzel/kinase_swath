rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "Figure1"

library(scales)
library(cowplot)
load("./R/objects/exp_metadata._clean_.RData")


# KEGG pathway Coverage -------------------------------------

load("./R/objects/KEGG.pathways.analysis1.RData")
load("./R/objects/kegg_categories._load_.RData")

KEGG.pathways.stats <- KEGG.pathways %>%
  group_by(pathway) %>% 
  summarize(EC.coverage  = length(unique(EC.number[is.inData == 1]))/length(unique(EC.number)),
            ORF.coverage = length(unique(ORF[is.inData == 1]))/length(unique(ORF)),
            EC.active  = length(unique(EC.number[is.Coupled == 1 ]))/length(unique(EC.number)),
            ORF.active = length(unique(ORF[is.Coupled == 1 ]))/length(unique(ORF)))

KEGG.pathways.stats2 <- KEGG.pathways %>% 
  filter(is.Coupled == 1) %>%
  group_by(pathway) %>% 
    summarize(EC.active.inData  = length(unique(EC.number[is.Coupled == 1 & is.inData == 1 ]))/length(unique(EC.number[is.Coupled == 1])))


KEGG.pathways.stats <- merge(KEGG.pathways.stats, KEGG.pathways.stats2, all = T)
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
  background_grid(major = "x", minor = "none") +
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

toPlot <- KEGG.pathways.stats.f %>% filter(EC.coverage > 0) %>%
  select(A, EC.coverage, EC.active.inData) %>%
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
  theme(legend.position = "none")

## -- CVs of mix samples -----

load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/fragments.matrix.quant.combat.RData")
load("./R/objects/peptides.matrix.quant.combat.RData")

data.raw = proteins.matrix.combat.quant
data.raw.f = data.raw[,grepl(pattern = "mix", ignore.case = T, x = colnames(data.raw) )]
data.raw.f.df = tbl_df(data.frame(ORF = row.names(data.raw.f), as.data.frame(data.raw.f)))

data.raw.f.proteins = data.raw.f.df %>%
  gather(sample, Tsignal, -ORF)
data.raw.f.proteins$type = "proteins"

data.raw = peptides.matrix.quant.combat
data.raw.f = data.raw[,grepl(pattern = "mix", ignore.case = T, x = colnames(data.raw) )]
data.raw.f.df = tbl_df(data.frame(ORF = row.names(data.raw.f), as.data.frame(data.raw.f)))

data.raw.f.peptides = data.raw.f.df %>%
  gather(sample, Tsignal, -ORF)
data.raw.f.peptides$type = "peptides"

data.raw = fragments.matrix.quant.combat
data.raw.f = data.raw[,grepl(pattern = "mix", ignore.case = T, x = colnames(data.raw) )]
data.raw.f.df = tbl_df(data.frame(ORF = row.names(data.raw.f), as.data.frame(data.raw.f)))

data.raw.f.fragments = data.raw.f.df %>%
  gather(sample, Tsignal, -ORF)
data.raw.f.fragments$type = "fragments"

data.all = rbind.data.frame(data.raw.f.fragments, data.raw.f.peptides, data.raw.f.proteins)
data.all$signal = exp(data.all$Tsignal)
data.all$type = factor(data.all$type, levels = c("fragments", "peptides", "proteins"))

data.all.summary = data.all %>% group_by(ORF, type) %>% summarise(CV = sd(signal, na.rm=T)/mean(signal, na.rm=T))


toPlot <- data.all.summary

p = ggplot(data = toPlot, aes(x=type, y=CV*100)) + 
  geom_violin() + 
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_y_log10(limits=c(5,100), breaks=c(5, 10, 15,20, 50,  100))+
  annotation_logticks(sides="l") +
  ylab("Signal variation for the duration of experiment, CV") +
  background_grid(major = "y", minor="y")


toPlot <- data.all.summary %>% filter(type == "proteins")
p.CV = ggplot(data = toPlot, aes(x=CV)) + 
  geom_density() + 
  scale_x_log10(limits=c(5,100)/100, breaks=c(5, 10, 15,20, 50,  100)/100, labels = percent_format())+
  annotation_logticks(sides="b") +
  xlab("Signal variation during the duration of experiment, CV%") +
  geom_vline(xintercept = genefilter::shorth(toPlot$CV), linetype=5) +
  background_grid(major = "x", minor="x")

## -- Volcano plots ---- 
load("./R/objects/proteins.matrix.combat.quant.FC.RData")
load("./R/objects/proteins.matrix.combat.quant.RData")

protein.matrix = proteins.matrix.combat.quant
proteins.FC = proteins.matrix.combat.quant.FC

reference = unique(as.character(proteins.FC$reference))

lb = -4
ub = 4
toPlot = proteins.FC
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


p1 = ggplot(toPlot, aes(y=-log10(p.value_BH), x=logFC)) +
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
  theme(aspect.ratio = 5/8, legend.position = "none")

p2 = ggplot(toPlot, aes(x=logFC)) +
  geom_histogram(colour = "white", fill = "black", binwidth = 0.25) +
  #geom_rect(aes(ymin=0, ymax=Inf, xmin=-FC_thr, xmax=FC_thr), inherit.aes=F, fill="grey", alpha = 0.01) +
  xlab(paste("Log2(fold-change)")) +
  geom_vline(xintercept = c(FC_thr,-FC_thr),linetype=3) +
  xlim(c(lb,ub)) +
  theme(aspect.ratio = 1)

p.volcano <- p1 + annotation_custom(grob = ggplotGrob(p2), xmin = 0, xmax = 5, ymin = 20 , ymax = 120)


# IRT chromatogram stability ------------------

load("./R/objects/dataset_irt_openswath._clean_.RData")
toPlot <- dataset_irt_openswath %>% 
  select(Sequence, filename, RT, ProteinName, m_score) %>% 
  filter(m_score < 0.01,  ProteinName == "1/Biognosys_iRT") %>% distinct()

toPlot$sample_name <- sub(pattern = ".mzML", x = basename(toPlot$filename), replacement = "")

toPlot <- toPlot %>% filter(sample_name %in% exp_metadata$sample_name)
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

View(toPlot)
p.irt.recalibrated <- ggplot(toPlot, aes(x=aquisition_date.str, y=iRT_recalibrated)) +
  geom_point(aes(colour=Sequence)) +
  geom_vline(data = toPlot[grepl("mix", toPlot$filename, ignore.case = T),], 
             aes(xintercept=as.numeric(aquisition_date.str)), linetype=1, alpha=0.75, colour="darkgrey") +
  xlab("Acquisition date") +
  ylab("Relative retention time") +
  #scale_y_date() +
  scale_colour_tableau() +
  theme(legend.position="none",
        aspect.ratio = 5/8)


## -- transcriptome vs proteome ----------------

load("./R/objects/proteins.matrix.combat.quant.FC.RData")
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
transcriptome.FC$isMetabolic = transcriptome.FC$ORF %in% unique(KEGG.pathways.f$ORF)
transcriptome.FC.f = transcriptome.FC[transcriptome.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]


tr.pr.FC = merge(transcriptome.FC, proteins.matrix.combat.quant.FC, by=c("KO", "ORF"), suffixes=c(".tr", ".pr"))

pval_thr = 0.05

tr.pr.cor.ORF = tr.pr.FC %>% group_by(ORF, isMetabolic) %>% 
  summarise(cor = cor(logFC.tr, logFC.pr, method="spearman"))

tr.pr.cor.KO = tr.pr.FC %>% group_by(KO, isMetabolic) %>% 
  summarise(cor = cor(logFC.tr, logFC.pr, method="spearman"))

p.metabolic_correlations <- ggplot(tr.pr.cor.KO, aes(x=cor, fill = isMetabolic)) +  #to Supplementary
  geom_density(alpha=0.5) +
  xlab(expression(paste("Spearman's correlation coefficient, ", rho))) +
  theme(legend.position = c(0.1, 0.8), aspect.ratio = 1)


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


p.tr_vs_pr_cor = ggplot(toPlot, aes(x = rev(KO.name), y = cor,  size=n)) + 
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.5,-0.25,0.25, 0.5), linetype = 3) +  
  ylab(expression(paste("Pearson correlation between\ngene and protein expression changes, ", r))) +
  xlab("Kinase mutant") +
  ylim(-0.8,0.8) +
  #guides(guide = guide_legend(tittle = "Deferentially expressed transcripts")) +
  scale_size(name="Number of common transcripts\nand proteins changed per mutant") +
  #coord_flip() + 
  theme(legend.position = c(0.25, 0.3),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

  
# Batch Effects ---------------------
load("./R/objects/peptides.matrix.RData")
load("./R/objects/peptides.matrix.quant.combat.RData")

before = peptides.matrix
after  = peptides.matrix.quant.combat

stopifnot(dim(before) == dim(after))

set.seed(1234)
exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
exp_metadata$batch_kmeans = pam(exp_metadata$aquisition_date.str, 7)$clustering

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
scores$batch_kmeans = factor(exp_metadata$batch_kmeans[match(scores$sample.id, exp_metadata$sample_name)])
scores$batch = factor(exp_metadata$batch.exp.n[match(scores$sample.id, exp_metadata$sample_name)])
scores.mix = scores[grepl(pattern="mix", ignore.case=T, x=rownames(scores)),]
annot$text = paste(annot$x_var, annot$y_var)
scores$type = factor(scores$type, levels=c("before", "after"))


library(cowplot)

p.batch <- ggplot(scores, aes(x=PC1, y=PC2)) + 
  geom_point(size=2, aes(col=batch_kmeans))+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(data=scores[grepl(pattern="mix", ignore.case=T, x=scores$sample.id),], 
             aes(x=PC1, y=PC2), size=3, col="black", shape=17) +
  geom_text(data = annot, aes(x=-50, y=-50, label=text)) +
  facet_wrap(~type, scales="fixed") + 
  scale_colour_tableau() +
  theme(aspect.ratio = 1,
        legend.position = "bottom")

file_name = "Figure1_v01_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(90, 60)))
  
  grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.irt, vp = viewport(layout.pos.row = 1:20, layout.pos.col = 1:30)) #experiment timeline with RT stability
  grid.text("B", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 31),gp=gpar(fontsize=20, col="black"))
  print(p.volcano, vp = viewport(layout.pos.row = 1:20, layout.pos.col = 31:60)) #Volcano plot
  grid.text("C", just=c("left", "centre"), vp = viewport(layout.pos.row = 21, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.coverage, vp = viewport(layout.pos.row = 21:40, layout.pos.col = 1:25))  #pathway coverage detailed
  print(p.pack_man, vp = viewport(layout.pos.row = 21:30, layout.pos.col = 26:35)) #pathway pie chart
  print(p.mean_cov, vp = viewport(layout.pos.row = 31:40, layout.pos.col = 26:35)) #pathway average
  grid.text("D", just=c("left", "centre"), vp = viewport(layout.pos.row = 21, layout.pos.col = 36),gp=gpar(fontsize=20, col="black"))
  #print(p.CV, vp = viewport(layout.pos.row = 21:30, layout.pos.col = 36:60)) #pathway average
  #grid.text("E", just=c("left", "centre"), vp = viewport(layout.pos.row = 21, layout.pos.col = 36),gp=gpar(fontsize=20, col="black"))
  print(p.batch, vp = viewport(layout.pos.row = 21:40, layout.pos.col = 36:60)) #pathway average
dev.off()


file_name = "Figure1_v02_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
grid.newpage() 
pushViewport(viewport(layout = grid.layout(135, 100)))

grid.text("a", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
print(p.irt.recalibrated, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 1:40)) #experiment timeline with RT stability
grid.text("b", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 41),gp=gpar(fontsize=20, col="black"))
print(p.coverage, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 41:80))  #pathway coverage detailed
print(p.pack_man, vp = viewport(layout.pos.row = 1:15, layout.pos.col = 81:100)) #pathway pie chart
print(p.mean_cov, vp = viewport(layout.pos.row = 16:30, layout.pos.col = 81:100)) #pathway average
grid.text("c", just=c("left", "centre"), vp = viewport(layout.pos.row = 31, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
print(p.volcano, vp = viewport(layout.pos.row = 31:60, layout.pos.col = 1:40)) #Volcano plot
grid.text("d", just=c("left", "centre"), vp = viewport(layout.pos.row = 31, layout.pos.col = 42),gp=gpar(fontsize=20, col="black"))
print(p.tr_vs_pr_cor, vp = viewport(layout.pos.row = 31:60, layout.pos.col = 45:100)) #Volcano plot

dev.off()





