rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "Figure2"

library(gridExtra)


#library(scales)
#library(cowplot)
load("./R/objects/pathway2orf._load_.RData")
load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/proteins.matrix.combat.quant.FC.RData")
load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/gene.annotations._load_.RData")
load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")


orf2name = droplevels(unique(gene.annotations[,c("V4", "V6")]))
orf2name$V4 = as.character(orf2name$V4)
orf2name$V6 = as.character(orf2name$V6)
orf2name$V6[orf2name$V6 == ""] = orf2name$V4[orf2name$V6 == ""]
names(orf2name) = c("ORF", "gene_name")

protein.matrix = proteins.matrix.combat.quant
proteins.FC = proteins.matrix.combat.quant.FC

reference = unique(as.character(proteins.FC$reference))

pval_thr = 0.01

FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)

load("./R/objects/KEGG.pathways.analysis1.RData")
EC.genes = gene.annotations[gene.annotations$V3 == "EC number",]

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

proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]
proteins.FC.f$isMetabolic = proteins.FC.f$ORF %in% unique(KEGG.pathways.f$ORF)
proteins.FC.f$isEnzyme = proteins.FC.f$ORF %in% unique(EC.genes$V4)
proteins.FC.f$isiMM904 = proteins.FC.f$ORF %in% unique(as.character(iMM904$gene))


all_proteins <- as.character(unique(proteins.FC.f$ORF))
all_measured_enzymes <- as.vector((proteins.FC.f %>% filter(isiMM904 ==T ) %>% dplyr::select(ORF) %>% distinct())$ORF)


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


# -- kinase specificity ----------

load("./R/objects/all_final_models.importance.models_summary2.RData")
load("./R/objects/all_final_models.models_summary2.RData")
load("./R/objects/file.list.models_summary2.RData")


metabolites.models.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  dplyr::select(model, RMSE, Rsquared, normalization, dataset, metabolite, degree, preprocessing) %>% 
  distinct() %>%
  group_by(metabolite, normalization, degree, preprocessing) %>%
  #group_by(metabolite) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite) %>% filter(degree <= 5) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))


predictors.dataset <- left_join(metabolites.models.long, all_final_models.importance) 

metabolite.predictors <- predictors.dataset %>% filter(Overall >= 80) %>%
  dplyr::select(metabolite, degree, gene ) %>% group_by(metabolite, degree, gene) %>% distinct()


proteins.FC.f.stats <- proteins.FC.f %>% 
  group_by(KO) %>%
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr) %>% 
  summarize(n_metabolic = sum(isMetabolic == T),
            n_yeast  = sum(isiMM904),
            n_total = n())

changed.genes <- proteins.FC.f %>% 
  group_by(KO) %>% arrange() %>%
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr) %>%
  dplyr::select(ORF, KO)

metabolite.predictors <- metabolite.predictors %>%
  group_by(metabolite, degree) %>%
  mutate(n_predictors = n())

changed.predictors <- inner_join(metabolite.predictors, changed.genes, by = c("gene" = "ORF" ))


kinase.distances <- ddply(changed.predictors, .(metabolite, degree, n_predictors), 
                          .fun = function(x) {
                            tmp.df <- dcast(data = x, formula = "KO ~ gene", fun.aggregate = length)
                            tmp.matrix <- as.matrix(tmp.df[,-1])
                            rownames(tmp.matrix) = tmp.df[,1]
                            d = dist(tmp.matrix, method = "binary")
                            d.matrix <- as.matrix(d)
                            d.matrix[upper.tri(d.matrix)] <- NA
                            diag(d.matrix) <- NA
                            return(melt(d.matrix))
                            
                          } )

metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]

toPlot <- kinase.distances %>% filter(metabolite %in% metabolite.order$metabolite, n_predictors >= 5)
toPlot$metabolite.label <- toupper(metabolite2iMM904$model_name[match(toPlot$metabolite, metabolite2iMM904$id)])
toPlot.stats <- toPlot %>% 
  group_by(metabolite, metabolite.label) %>%
  summarise(median = median(1 - value, na.rm = T)) %>% 
  ungroup() %>%
  arrange(median)



toPlot$metabolite <- factor(toPlot$metabolite, levels = rev(toPlot.stats$metabolite))

toPlot$metabolite.label <- factor(toPlot$metabolite.label, levels = rev(unique(toPlot.stats$metabolite.label)))

p.specificity <- ggplot(toPlot, aes(x = 1 - value)) +
  geom_density(fill="black") +
  facet_grid(metabolite.label ~ . ) +
  xlab("Predictors response similarity following kinase deletion, Jaccard index") +
  #theme_bw() +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        strip.text.y = element_text(angle=0),
        strip.background = element_rect(colour = NULL), aspect.ratio = 1/10) 

p.specificity$toScale <- F

plots.list <- lappend(plots.list, p.specificity)

p.specificity.h <- ggplot(toPlot, aes(x = 1 - value)) +
  geom_density(fill="black") +
  facet_wrap(~ metabolite, scales = "free_y") +
  xlab("Predictors response similarity following kinase deletion, Jaccard index") +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        #strip.text.y = element_text(angle=0),
        strip.background = element_rect(colour = NULL))


toPlot.stats$metabolite.label <- metabolite2iMM904$official_name[match(toPlot.stats$metabolite, metabolite2iMM904$id)] 

toPlot.stats$y = seq(from = 2, to = 0.3, length.out = nrow(toPlot.stats))

p.specificity_hist <- ggplot(toPlot, aes(x = 1 - value)) +
  geom_density(fill="black") +
  geom_segment(data=toPlot.stats, aes(x = median, y = y, xend=median , yend=0), 
               colour="red") +
  geom_text(data=toPlot.stats, aes(x = median, y = y, label = metabolite.label), 
            colour="red", check_overlap = TRUE, hjust=-0.1, vjust=0, size=2) +
  geom_point(data=toPlot.stats, aes(x = median, y = y), colour="red", size=1) +
  xlab("Predictors response similarity following kinase deletion, Jaccard index") +
  theme_bw()
  #scale_x_log10(breaks=seq(0,1,0.1))+
  #annotation_logticks(sides="b") +



## ------ enzyme overlaps --------

proteins.FC.f.metabolic <- tbl_df(proteins.FC.f) %>% 
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr, isiMM904 == T)

# all enzymes
proteins.FC.f.metabolic$KO.gene <- orf2name$gene_name[match(proteins.FC.f.metabolic$KO, orf2name$ORF)]
stopifnot(!any(is.na(proteins.FC.f.metabolic$KO.gene)))

x = proteins.FC.f.metabolic
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.all <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))


#upregulated
x = proteins.FC.f.metabolic %>% filter(logFC > 0)
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.up <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))

#downregulated
x = proteins.FC.f.metabolic %>% filter(logFC < 0)
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.down <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))


FC.f.metabolic.stats <- proteins.FC.f.metabolic %>% 
  group_by(KO.gene) %>% 
  summarise(n = n(),
            n_pos = sum(logFC > 0)/length(all_measured_enzymes),
            n_neg = sum(logFC < 0)/length(all_measured_enzymes)) %>% 
  ungroup() %>% arrange(n)

cl = hclust(dist(d.matrix.all))
cl <- dendextend::rotate(cl, order = as.character(FC.f.metabolic.stats$KO.gene))
d.matrix.all <- d.matrix.all[cl$order,cl$order]


d.matrix.up <- d.matrix.up[rownames(d.matrix.up)[match(rownames(d.matrix.all),rownames(d.matrix.up))], 
                           colnames(d.matrix.up)[match(colnames(d.matrix.all),colnames(d.matrix.up))]]
d.matrix.down <- d.matrix.down[rownames(d.matrix.down)[match(rownames(d.matrix.all),rownames(d.matrix.down))], 
                           colnames(d.matrix.down)[match(colnames(d.matrix.all),colnames(d.matrix.down))]]


zeros.matrix <- matrix(data=0, ncol = ncol(d.matrix.down), nrow = nrow(d.matrix.up))
zeros.matrix[upper.tri(zeros.matrix)] <- d.matrix.up[upper.tri(d.matrix.up)]
zeros.matrix[lower.tri(zeros.matrix)] <- d.matrix.down[lower.tri(d.matrix.down)]*-1

similarities <- list()
tmp <- data.frame(similarity = d.matrix.up[upper.tri(d.matrix.up)],
                  type = "up",
                  dataset = "metabolic")
similarities <- lappend(similarities, tmp)

tmp <- data.frame(similarity = d.matrix.down[upper.tri(d.matrix.down)],
                  type = "down",
                  dataset = "metabolic")
similarities <- lappend(similarities, tmp)

tmp <- data.frame(similarity = d.matrix.all[upper.tri(d.matrix.all)],
                  type = "all",
                  dataset = "metabolic")
similarities <- lappend(similarities, tmp)



toPlot <- melt(zeros.matrix)
toPlot$x.name <- factor(rownames(d.matrix.all)[toPlot$X1], levels = rownames(d.matrix.all))
toPlot$y.name <- factor(colnames(d.matrix.all)[toPlot$X2], levels = colnames(d.matrix.all))

my_breaks <- seq(1, -1, -0.25)
my_colours <- rev(brewer.pal(name = "RdBu", n = length(my_breaks) - 1))
my_colours[c(4,5)] <- "white"
#pheatmap(d.matrix.all, cluster_rows = F, cluster_cols = F)

p.heatmap <- ggplot(toPlot) +  
  geom_tile(aes(x = x.name, y = y.name, fill = cut(value, breaks = my_breaks)), colour="grey") +
#   scale_fill_gradient2(low="#1F78B4",high="#E31A1C",mid ="white",
#                        breaks = seq(-0.75, 0.75, 0.25),
#                        midpoint=0)  +
  scale_fill_manual(values = my_colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        aspect.ratio = 1, legend.position = c(0.2, 0.8) ) +
  labs(x="", y = "")

toPlot <- FC.f.metabolic.stats %>% 
  dplyr::select(KO.gene, n_pos, n_neg) %>% as.data.frame()%>%
  melt(id.vars = "KO.gene") 
toPlot$KO.gene <- factor(toPlot$KO.gene, levels = rownames(d.matrix.all)) 

library(ggthemes)
p.barplot = ggplot(toPlot, aes(x=KO.gene, y=value, fill=variable)) + 
  geom_bar(stat="identity", width=.5) + 
  labs(x = "", y = "Fraction of perturbed metabolic network") +
  coord_flip() + 
  scale_fill_manual(values = my_colours[c(length(my_colours),1)]) +
  theme_few() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.7, 0.2))


p.barplot.h = ggplot(toPlot, aes(x=KO.gene, y=value, fill=variable)) + 
  geom_bar(stat="identity", width=.5) + 
  labs(x = "", y = "Fraction of perturbed metabolic network") +
  #coord_flip() + 
  scale_fill_manual(values = my_colours[c(length(my_colours),1)]) +
  theme_few() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.2, 0.7))


#### ---- overlaps up/down for supplementary --------------------
zeros.matrix <- matrix(data=0, ncol = ncol(d.matrix.all), nrow = nrow(d.matrix.all))
zeros.matrix[upper.tri(zeros.matrix)] <- d.matrix.all[upper.tri(d.matrix.all)]
zeros.matrix[lower.tri(zeros.matrix)] <- d.matrix.all[lower.tri(d.matrix.all)]
diag(zeros.matrix) <- max(zeros.matrix)
toPlot <- melt(zeros.matrix)

toPlot$x.name <- factor(rownames(d.matrix.all)[toPlot$X1], levels = rownames(d.matrix.all))
toPlot$y.name <- factor(colnames(d.matrix.all)[toPlot$X2], levels = colnames(d.matrix.all))

p.heatmap_all <- ggplot(toPlot) +  
  geom_tile(aes(x = x.name, y = y.name, fill = value), colour="grey")+
  scale_fill_gradient2(low="white", high = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        aspect.ratio = 1, legend.position = c(0.2, 0.8) ) +
  labs(x="", y = "")

cl$labels = ""
p.dendro = ggdendro::ggdendrogram(cl, rotate=T)

gp1<-ggplotGrob(p.dendro)
gp2<-ggplotGrob(p.heatmap_all)
 
maxWidth = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxWidth)
gp2$heights[2:5] <- as.list(maxWidth)
p.combined <- grid.arrange(gp2, gp1, ncol=2,widths=c(7/10,3/10))

p.combined$toScale <- T
plots.list <- lappend(plots.list, p.combined)

### ----- heatmap of all genes ---------------


proteins.FC.f.non_metabolic <- tbl_df(proteins.FC.f) %>% 
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr, isiMM904 == F)

# all enzymes
proteins.FC.f.non_metabolic$KO.gene <- orf2name$gene_name[match(proteins.FC.f.non_metabolic$KO, orf2name$ORF)]
stopifnot(!any(is.na(proteins.FC.f.non_metabolic$KO.gene)))

x = proteins.FC.f.non_metabolic
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.all <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))


#upregulated
x = proteins.FC.f.non_metabolic %>% filter(logFC > 0)
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.up <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))

#downregulated
x = proteins.FC.f.non_metabolic %>% filter(logFC < 0)
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.down <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))

non_metabolic_proteins <- all_proteins[!all_proteins %in% all_measured_enzymes]
FC.f.non_metabolic.stats <- proteins.FC.f.non_metabolic %>% 
  group_by(KO.gene) %>% 
  summarise(n = n(),
            n_pos = sum(logFC > 0)/length(non_metabolic_proteins),
            n_neg = sum(logFC < 0)/length(non_metabolic_proteins)) %>% 
  ungroup() %>% arrange(n)

cl = hclust(dist(d.matrix.all))
cl <- dendextend::rotate(cl, order = as.character(FC.f.non_metabolic.stats$KO.gene))
d.matrix.all <- d.matrix.all[cl$order,cl$order]


d.matrix.up <- d.matrix.up[rownames(d.matrix.up)[match(rownames(d.matrix.all),rownames(d.matrix.up))], 
                           colnames(d.matrix.up)[match(colnames(d.matrix.all),colnames(d.matrix.up))]]
d.matrix.down <- d.matrix.down[rownames(d.matrix.down)[match(rownames(d.matrix.all),rownames(d.matrix.down))], 
                               colnames(d.matrix.down)[match(colnames(d.matrix.all),colnames(d.matrix.down))]]


zeros.matrix <- matrix(data=0, ncol = ncol(d.matrix.down), nrow = nrow(d.matrix.up))
zeros.matrix[upper.tri(zeros.matrix)] <- d.matrix.up[upper.tri(d.matrix.up)]
zeros.matrix[lower.tri(zeros.matrix)] <- d.matrix.down[lower.tri(d.matrix.down)]*-1

tmp <- data.frame(similarity = d.matrix.up[upper.tri(d.matrix.up)],
                  type = "up",
                  dataset = "non_metabolic")
similarities <- lappend(similarities, tmp)

tmp <- data.frame(similarity = d.matrix.down[upper.tri(d.matrix.down)],
                  type = "down",
                  dataset = "non_metabolic")
similarities <- lappend(similarities, tmp)

tmp <- data.frame(similarity = d.matrix.all[upper.tri(d.matrix.all)],
                  type = "all",
                  dataset = "non_metabolic")
similarities <- lappend(similarities, tmp)



toPlot <- melt(zeros.matrix)
toPlot$x.name <- factor(rownames(d.matrix.all)[toPlot$X1], levels = rownames(d.matrix.all))
toPlot$y.name <- factor(colnames(d.matrix.all)[toPlot$X2], levels = colnames(d.matrix.all))

my_breaks <- seq(1, -1, -0.25)

my_colours <- rev(brewer.pal(name = "RdBu", n = length(my_breaks) - 1))
my_colours[c(4,5)] <- "white"

p.heatmap_non_metabolic <- ggplot(toPlot) +  
  geom_tile(aes(x = x.name, y = y.name, fill = cut(value, breaks = my_breaks)), colour="grey") +
  #   scale_fill_gradient2(low="#1F78B4",high="#E31A1C",mid ="white",
  #                        breaks = seq(-0.75, 0.75, 0.25),
  #                        midpoint=0)  +
  scale_fill_manual(values = my_colours, name = "Perturbation overlap") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = c(0.1, 0.8), 
        legend.background = element_rect(colour = NA),
        aspect.ratio = 1) +
  labs(x="", y = "")


toPlot <- FC.f.non_metabolic.stats %>% 
  dplyr::select(KO.gene, n_pos, n_neg) %>% as.data.frame()%>%
  melt(id.vars = "KO.gene") 
toPlot$KO.gene <- factor(toPlot$KO.gene, levels = rownames(d.matrix.all)) 


p.barplot_non_metabolic <- ggplot(toPlot, aes(x=KO.gene, y=value, fill=variable)) + 
  geom_bar(stat="identity", width=.5) + 
  labs(x = "", y = "Perturbed fraction of non metabolic proteins") +
  coord_flip() + 
  scale_fill_manual(values = my_colours[c(length(my_colours),1)], name = "") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.7, 0.1),
        legend.background = element_rect(colour = NA)) 



pEmpty <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_blank() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank())

p.combined <- grid.arrange(pEmpty + theme(aspect.ratio = 4/1), pEmpty + theme(aspect.ratio = 1/4), p.heatmap_non_metabolic, p.barplot_non_metabolic, ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
p.combined$toScale <- T
plots.list <- lappend(plots.list, p.combined)

# --- similarities of metabolic vs non metabolic -------------------

similarities.df <- bind_rows(similarities)
toPlot <- similarities.df %>% filter(type == "all")
library(ggthemes)
p.value <- (wilcox.test(toPlot$similarity[toPlot$dataset == "metabolic"], toPlot$similarity[toPlot$dataset != "metabolic"]))$p.value

s.similarities <- ggplot(toPlot, aes(y = similarity, x = dataset)) +
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_boxplot(width = 0.5) +
  scale_x_discrete(labels = c("metabolic" = "Metabolic enzymes", "non_metabolic" = "The rest proteins")) + 
  xlab("Subset of dataset") +
  annotate("text", x = 1.5, y = 0.7, label = paste("P-value = ", base::format(p.value, scientific = T, digits = 2 ) )) +
  theme_bw() + theme(aspect.ratio = 1)

plots.list <- lappend(plots.list, s.similarities)


# --- saturation plots ---------------
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

overlaps.stats$degree <- degree(GRAPH)[match(overlaps.stats$V1, names(degree(GRAPH)))]

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
  
  
toText <- cor.test(toPlot$degree, toPlot$n_fraction, use = "pairwise.complete.obs" )

p.degree_vs_metabolic <- ggplot(toPlot, aes(x=degree, y = n_fraction)) +
  geom_point() + 
  annotate(geom = "text", x = 100, y=0.4, 
            label = paste("r = ", round(toText$estimate,2), " p-value = ", round(toText$p.value,2), sep="")) +
  theme_bw() + theme(aspect.ratio = 1) +
  xlab("Number of protein-protein interactions kinase involved") +
  ylab("Perturbed fraction of metabolic network")

plots.list <- lappend(plots.list, p.degree_vs_metabolic)


proteins.FC.f.stats$n_yeast_fraction <- proteins.FC.f.stats$n_metabolic/proteins.FC.f.stats$n_total
toPlot <- proteins.FC.f.stats 
p.constant_fraction <- ggplot(toPlot, aes(x = n_total, y = n_yeast_fraction)) +
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Number of proteins changes per mutant") +
  ylab("Fraction of metabolic enzymes affected by kinase") +
  theme_bw()




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




#  --- pathway enrichment ---- 
enrich_thr = 0.05
proteins.FC = tbl_df(proteins.FC)
proteins.FC$KO = factor(proteins.FC$KO)

all_enrichments = ddply(filter(proteins.FC.f, isMetabolic==T), .(KO), 
                        .fun = function(z) {
                          #z = proteins.z_scores.long[proteins.z_scores.long$KO == "WT",]
                          signal = z$ORF[z$p.value_BH < pval_thr & abs(z$logFC) > FC_thr]
                          universe =  z$ORF              
                   
                          tmp.path = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=pathway2orf[pathway2orf$pathway %in% KEGG.pathways.f$pathway,])
                          tmp.path$type = "kegg_pathways"
                          
                          return(tmp.path)
                          
                        })


kegg.enrichments = dcast(droplevels(all_enrichments[all_enrichments$type == "kegg_pathways",]), pathway~KO, value.var="p.value")
kegg.enrichments.matrix = as.matrix(kegg.enrichments[,-1])
rownames(kegg.enrichments.matrix) = kegg.enrichments$pathway

kegg.enrichments.matrix[is.na(kegg.enrichments.matrix)] = 0
kegg.enrichments.matrix[kegg.enrichments.matrix > enrich_thr] = 0
kegg.enrichments.matrix[kegg.enrichments.matrix != 0] = 1
proteins.FC.f.stats$enrichment = as.vector(colSums(kegg.enrichments.matrix)[match(proteins.FC.f.stats$KO, colnames(kegg.enrichments.matrix))])

#calculate number of measured ORFs per KEGG metabolic pathway
KEGG.pathways.f.summary = KEGG.pathways.f[KEGG.pathways.f$ORF %in% rownames(protein.matrix),] %>% group_by(pathway, description) %>% summarize(n=n())
selected.pathways = as.character(KEGG.pathways.f.summary$pathway[KEGG.pathways.f.summary$n >= 5])

kegg.enrichments.matrix.f = kegg.enrichments.matrix[rownames(kegg.enrichments.matrix) %in% selected.pathways,]
kegg.enrichments.matrix.f = kegg.enrichments.matrix.f[rowSums(kegg.enrichments.matrix.f) >=2,]
kegg.enrichments.matrix.f = kegg.enrichments.matrix.f[,colSums(kegg.enrichments.matrix.f) >=2]

pathway2desription = pathway2orf %>% distinct(pathway, description) %>% dplyr::select(pathway, description)
pathway2desription$description = sub(pattern=" - Saccharomyces cerevisiae (budding yeast)", replacement="" , x=pathway2desription$description, fixed=T)
pathway2desription$description = sub(pattern="metabolism", replacement="m." , x=pathway2desription$description, fixed=T)
pathway2desription$description = sub(pattern="interconversions", replacement="int." , x=pathway2desription$description, fixed=T)

kegg.enrichments.matrix.kinases = kegg.enrichments.matrix.f[,colnames(kegg.enrichments.matrix.f) %in% kinases]
toPlot = melt(kegg.enrichments.matrix.kinases)

toPlot$KO.name = exp_metadata$gene[match(toPlot$X2, exp_metadata$ORF)]
toPlot$description = pathway2desription$description[match(toPlot$X1, pathway2desription$pathway)]

toPlot <- toPlot %>% group_by(X1) %>% mutate(n = sum(value)) %>% ungroup()  %>% arrange(n)
toPlot$description <- factor(toPlot$description, levels = as.character(unique(toPlot$description)) )
p.kegg_enrich <- ggplot(toPlot, aes(x = KO.name, y = description, fill=factor(value))) +
  geom_tile(colour="grey")  +
  scale_fill_manual(values = c("lightgrey", "black")) +
  theme_bw() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")



# --- slopes -----

my_means <- function(proteins.matrix) {
  
  proteins.long = melt(proteins.matrix, id.vars="rownames")
  names(proteins.long) = c("EG.StrippedSequence", "R.Label", "signal")
  proteins.long$ORF = exp_metadata$ORF[match(proteins.long$R.Label, exp_metadata$sample_name)]
  proteins.long.mean = tbl_df(proteins.long) %>% group_by(EG.StrippedSequence, ORF) %>% summarize(mean = mean(signal))
  proteins.mean.df = dcast(proteins.long.mean, formula=EG.StrippedSequence~ORF, value.var="mean")
  
  proteins.mean.matrix = as.matrix(proteins.mean.df[,-1])
  rownames(proteins.mean.matrix) = as.matrix(proteins.mean.df$EG.StrippedSequence)
  return(proteins.mean.matrix)  
}


protein.matrix.mean = my_means(exp(proteins.matrix.combat.quant))

all.kinases <- c("WT", unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])))
combinations <- combn(all.kinases, 2)
combinations.df <- as.data.frame(t(combinations))


slopes <- ddply(combinations.df, .(V1, V2),
                  .fun = function(x) {
                    var1 <- protein.matrix.mean[,x$V1]
                    var2 <- protein.matrix.mean[,x$V2]
                    tmp <- data.frame(var1, var2)
                    fit <- lm(var1~var2, data=tmp)
                    
                    result <- as.vector(fit$coefficients[2])
                    cor.tmp <- cor(tmp$var2, tmp$var1)
                      return(data.frame(slope = result, 
                                        cor = cor.tmp))
                  }  )




toPlot <- slopes %>% filter(V1 == "WT")
toPlot$label <- orf2name$gene_name[match(toPlot$V2, orf2name$ORF)]


toPlot$n_yeast <- proteins.FC.f.stats$changes[match(toPlot$V2, proteins.FC.f.stats$KO)]

toPlot.stats = data.frame(cor = (cor.test(toPlot$cor, toPlot$slope, use = "pairwise.complete.obs"))$estimate,
                          p.value = (cor.test(toPlot$cor, toPlot$slope, use = "pairwise.complete.obs"))$p.value)

toPlot$label <- orf2name$gene_name[match(toPlot$V2, orf2name$ORF)]
s.slopes_vs_cor <- ggplot(toPlot, aes(x = slope, y = cor)) +
  geom_point(aes(size = n_yeast), colour = "darkgrey") +
  geom_text(aes(label = label), check_overlap = T) +
  #geom_smooth(method = lm, se = F) +
  #annotate("text", x = 1, y = 0.9, label = paste( "r = ", round(toPlot.stats$cor, 2), "p-value = ", format(toPlot.stats$p.value,digits = 2, scientific = T))) +
  scale_size_continuous(name="Total metabolic\nenzymes perturbed") +
  theme_bw() +
  theme(aspect.ratio = 1, 
        legend.position = c(0.2, 0.7), 
        legend.background = element_rect(colour = NA)) +
  xlab("Linear regression slope between WT and kinase mutants proteome profiles") +
  ylab("Proteome profile correalations between WT and kinase mutants, Pearson's r")


toPlot <- slopes
s.slopes <- ggplot(toPlot, aes(x = slope)) + 
       geom_histogram(colour = "white") +
       theme_bw() + 
       annotate("text", x=0.6, y = 50, label = round(min(toPlot$slope),2)) +
       annotate("text", x=1.2, y = 50, label = round(max(toPlot$slope),2)) +
       theme(aspect.ratio = 1) +
       xlab("Linear regression slope between kinase mutants proteome profiles")
        

s.cors <- ggplot(toPlot, aes(x = cor)) + 
  geom_histogram(colour = "white") +
  theme_bw() + 
  #annotate("text", x=0.6, y = 50, label = round(min(toPlot$cor),2)) +
  #annotate("text", x=1.2, y = 50, label = round(max(toPlot$cor),2)) +
  theme(aspect.ratio = 1) +
  xlab("Proteome profile correalations between WT and kinase mutants, Pearson's r")

s.combined.slopes <- arrangeGrob(s.slopes, s.cors, s.slopes_vs_cor, nrow = 1)
s.combined.slopes$toScale <- T

plots.list <- lappend(plots.list, s.combined.slopes)



# -- Figure 2 -------

plot_figure_v2 <- function() {
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(135, 100)))
  
  
  print(p.heatmap, vp = viewport(layout.pos.row = 1:80, layout.pos.col = 1:80)) #all overlaps of perturbations
  grid.text("a", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.barplot, vp = viewport(layout.pos.row = 1:80, layout.pos.col = 81:100))
  grid.text("b", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 81),gp=gpar(fontsize=20, col="black"))
  print(p.kegg_enrich, vp = viewport(layout.pos.row = 80:100, layout.pos.col = 1:50))
  grid.text("c", just=c("left", "centre"), vp = viewport(layout.pos.row = 80, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.saturation, vp = viewport(layout.pos.row = 80:95, layout.pos.col = 51:75))
  grid.text("d", just=c("left", "centre"), vp = viewport(layout.pos.row = 80, layout.pos.col = 51),gp=gpar(fontsize=20, col="black"))
  print(p.constant_fraction, vp = viewport(layout.pos.row = 80:95, layout.pos.col = 76:100))
  grid.text("e", just=c("left", "centre"), vp = viewport(layout.pos.row = 80, layout.pos.col = 76),gp=gpar(fontsize=20, col="black"))
  print(p.specificity_hist, vp = viewport(layout.pos.row = 100:135, layout.pos.col = 1:66))
  grid.text("f", just=c("left", "centre"), vp = viewport(layout.pos.row = 100, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  
  print(p.dist_changes, vp = viewport(layout.pos.row = 105:130, layout.pos.col = 67:77))
  grid.text("g", just=c("left", "centre"), vp = viewport(layout.pos.row = 100, layout.pos.col = 67),gp=gpar(fontsize=20, col="black"))
  print(p.dist_degree, vp = viewport(layout.pos.row = 105:130, layout.pos.col = 78:88))
  print(p.dist_between, vp = viewport(layout.pos.row = 105:130, layout.pos.col = 89:99))
  
}

file_name = "Figure2_v02_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")
pdf(file_path, height=247/25.4*2, width=183/25.4*2)
  plot_figure_v2()
dev.off()


file_name = "Figure2_v02_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
  plot_figure_v2() 
dev.off()


file_name = paste("supplementary", fun_name, sep = ".")
file_path = paste(figures_dir, file_name, sep="/")

lapply(seq_along(plots.list) , 
       function(x) {
         
         tryCatch({
           p <- plots.list[[x]]
           scale = 1
           if (length(p$toScale) != 0 && p$toScale == T  ){
             scale = 2
           }
           ggplot2::ggsave(filename = paste(file_path, x , "pdf", sep = "."), device = NULL,
                           plot = p, width = 210 , height = 297, units = "mm", scale = scale)
         }, error = function(e) {
           message(paste("Plot", "x", "sucks!" ))
           return(NULL)
         }, finally = {
           message(paste("processed plot", x))
         })
         
         
       })


