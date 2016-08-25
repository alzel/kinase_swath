library(plyr);library(dplyr)
library(reshape2)
library(ggplot2)
library(sva)
library(tidyr)

rm(list = ls())
output_dir = "./R/objects"
figures_dir = "./figures"
fun_name = "signaling_similarity"
dir.create(output_dir)
dir.create(figures_dir)


# -- load pathways ------ 
ReactomePathways <- read.delim("./data/2016-04-29/ReactomePathways.txt", header=FALSE, stringsAsFactors = F)
names(ReactomePathways) <- c("reactome_id", "description", "species")
ReactomePathways <- ReactomePathways %>% filter(species == "Saccharomyces cerevisiae")


UniProt2Reactome <- read.delim("./data/2016-04-29/UniProt2Reactome.txt", header=FALSE, stringsAsFactors = F)
names(UniProt2Reactome) <- c("uniprot_id", "reactome_id", "url", "description", "confidence", "species")
UniProt2Reactome <- tbl_df(UniProt2Reactome) %>% 
  filter(species == "Saccharomyces cerevisiae") %>% 
  group_by(reactome_id, uniprot_id) %>% distinct()


ReactomePathwaysRelation <- read.delim("./data/2016-04-29/ReactomePathwaysRelation.txt", header=FALSE, stringsAsFactors = F)
names(ReactomePathwaysRelation) <- c("parent", "child")


# --- analysis ----

load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/gene.annotations._load_.RData")
load("./R/objects/similarities.long.Figure2.RData")
load("./R/objects/iMM904._load_.RData")

load("./R/objects/proteins.matrix.sva.0.5.1.RData")
load("./R/objects/proteins.matrix.sva.0.5.1.FC.RData")

proteins.matrix <- proteins.matrix.sva.0.5.1[rownames(proteins.matrix.sva.0.5.1) %in% unique(as.character(iMM904$gene)),]
proteins.FC <- proteins.matrix.sva.0.5.1.FC[proteins.matrix.sva.0.5.1.FC$ORF %in% unique(as.character(iMM904$gene)),]

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

protein.matrix.mean = my_means(exp(proteins.matrix))


kinase_orfs <- unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"]))

uniprot2orf.kinases <- gene.annotations %>%  filter(V3 ==  "UniProt/Swiss-Prot ID", V4 %in% kinase_orfs) %>%
  dplyr::select(V1, V4, V6) %>% distinct()
uniprot2orf.kinases[] <- lapply(uniprot2orf.kinases, as.character)
names(uniprot2orf.kinases) <- c("uniprot_id", "ORF", "gene_name")
uniprot2orf.kinases$gene_name[uniprot2orf.kinases$gene_name == ""] <- uniprot2orf.kinases$ORF[uniprot2orf.kinases$gene_name == ""]


#overlap similarities upregulated
similarities.long.df.up <- droplevels(bind_rows(similarities.long) %>% filter(type == "up", dataset == "metabolic"))

similarities.long.df.up$X1_uniprot <- uniprot2orf.kinases$uniprot_id[match(similarities.long.df.up$X1, uniprot2orf.kinases$gene_name)]
similarities.long.df.up$X2_uniprot<- uniprot2orf.kinases$uniprot_id[match(similarities.long.df.up$X2, uniprot2orf.kinases$gene_name)]
similarities.long.df.up$X1 <- uniprot2orf.kinases$ORF[match(similarities.long.df.up$X1, uniprot2orf.kinases$gene_name)]
similarities.long.df.up$X2 <- uniprot2orf.kinases$ORF[match(similarities.long.df.up$X2, uniprot2orf.kinases$gene_name)]

similarities.long.df.up$sim_type  = "overlap.up"
similarities.long.df.up$value <- with(similarities.long.df.up, ifelse(X1 == X2, NA, value))
similarities.long.df.up <- droplevels(similarities.long.df.up[!is.na(similarities.long.df.up$value),])


#overlap similarities downregulated
similarities.long.df.down <- droplevels(bind_rows(similarities.long) %>% filter(type == "down", dataset == "metabolic"))

similarities.long.df.down$X1_uniprot <- uniprot2orf.kinases$uniprot_id[match(similarities.long.df.down$X1, uniprot2orf.kinases$gene_name)]
similarities.long.df.down$X2_uniprot<- uniprot2orf.kinases$uniprot_id[match(similarities.long.df.down$X2, uniprot2orf.kinases$gene_name)]
similarities.long.df.down$X1 <- uniprot2orf.kinases$ORF[match(similarities.long.df.down$X1, uniprot2orf.kinases$gene_name)]
similarities.long.df.down$X2 <- uniprot2orf.kinases$ORF[match(similarities.long.df.down$X2, uniprot2orf.kinases$gene_name)]

similarities.long.df.down$sim_type  = "overlap.down"
similarities.long.df.down$value <- with(similarities.long.df.down, ifelse(X1 == X2, NA, value))
similarities.long.df.down <- droplevels(similarities.long.df.down[!is.na(similarities.long.df.down$value),])


#overlap similarities all
similarities.long.df.all <- droplevels(bind_rows(similarities.long) %>% filter(type == "all", dataset == "metabolic"))
similarities.long.df.all$X1_uniprot <- uniprot2orf.kinases$uniprot_id[match(similarities.long.df.all$X1, uniprot2orf.kinases$gene_name)]
similarities.long.df.all$X2_uniprot<- uniprot2orf.kinases$uniprot_id[match(similarities.long.df.all$X2, uniprot2orf.kinases$gene_name)]
similarities.long.df.all$X1 <- uniprot2orf.kinases$ORF[match(similarities.long.df.all$X1, uniprot2orf.kinases$gene_name)]
similarities.long.df.all$X2 <- uniprot2orf.kinases$ORF[match(similarities.long.df.all$X2, uniprot2orf.kinases$gene_name)]

similarities.long.df.all$sim_type  = "overlap.all"
similarities.long.df.all$value <- with(similarities.long.df.all, ifelse(X1 == X2, NA, value))
similarities.long.df.all <- droplevels(similarities.long.df.all[!is.na(similarities.long.df.all$value),])



#protein correlation based

proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]
proteins.FC.f.wide <- dcast(proteins.FC.f, formula = "ORF~KO", value.var = "logFC")
proteins.FC.f.matrix <- as.matrix(proteins.FC.f.wide[,-1])
rownames(proteins.FC.f.matrix) <- proteins.FC.f.wide$ORF

tmp.cor <- cor(proteins.FC.f.matrix)
diag(tmp.cor) = NA
proteins.FC.f.matrix.cor.long <- melt(tmp.cor, varnames = c("X1", "X2"))
proteins.FC.f.matrix.cor.long$X1_uniprot <- uniprot2orf.kinases$uniprot_id[match(proteins.FC.f.matrix.cor.long$X1, uniprot2orf.kinases$ORF)]
proteins.FC.f.matrix.cor.long$X2_uniprot <- uniprot2orf.kinases$uniprot_id[match(proteins.FC.f.matrix.cor.long$X2, uniprot2orf.kinases$ORF)]
proteins.FC.f.matrix.cor.long$sim_type  = "pearson_fc"
proteins.FC.f.matrix.cor.long <- proteins.FC.f.matrix.cor.long[!is.na(proteins.FC.f.matrix.cor.long$value),]



toSelect = colnames(protein.matrix.mean)[colnames(protein.matrix.mean) %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"]))]
tmp.cor <- cor(protein.matrix.mean[,toSelect])
tmp.cor <- cor(protein.matrix.mean)

diag(tmp.cor) = NA
protein.matrix.mean.cor.long <- melt(tmp.cor, varnames = c("X1", "X2"))
protein.matrix.mean.cor.long$X1_uniprot <- uniprot2orf.kinases$uniprot_id[match(protein.matrix.mean.cor.long$X1, uniprot2orf.kinases$ORF)]
protein.matrix.mean.cor.long$X2_uniprot <- uniprot2orf.kinases$uniprot_id[match(protein.matrix.mean.cor.long$X2, uniprot2orf.kinases$ORF)]
protein.matrix.mean.cor.long$sim_type  = "pearson"
protein.matrix.mean.cor.long <- protein.matrix.mean.cor.long[!is.na(protein.matrix.mean.cor.long$value),]


UniProt2Reactome.kinases <- tbl_df(UniProt2Reactome) %>% 
  dplyr::filter(uniprot_id %in% uniprot2orf.kinases$uniprot_id) %>% 
  group_by(reactome_id) %>%
  mutate(n = n()) %>% 
  filter(n >1) %>% arrange(reactome_id)

UniProt2Reactome.kinases <- rename(UniProt2Reactome.kinases, pathway = reactome_id)


load("./R/objects/pathway2orf._load_.RData")
pathway2orf$uniprot_id = uniprot2orf.kinases$uniprot_id[match(pathway2orf$ORF, uniprot2orf.kinases$ORF)]



orf2kegg.kinases <- pathway2orf %>% 
  dplyr::filter(!is.na(uniprot_id)) %>%
  group_by(pathway) %>%
  mutate(n = n()) %>% 
  filter(n >1) %>% arrange(pathway) %>% 
  dplyr::select(pathway, uniprot_id, n)



orf2kegg.kinases$pathway_base = "kegg"
UniProt2Reactome.kinases$pathway_base = "reactome"

pathway_kinases <- bind_rows(orf2kegg.kinases, UniProt2Reactome.kinases)


pathway_similarities <- ddply(pathway_kinases, .(pathway, pathway_base), 
                               .fun = function(x) {
                                 return.list = list()
                                 
                                 
                                 tmp <- data.frame(t(combn(unique(x$uniprot_id), 2)))
                                 tmp$n = unique(x$n)
                                 return.list$similarities.up <- left_join(tmp, similarities.long.df.up, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
                                 return.list$similarities.down <- left_join(tmp, similarities.long.df.down, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
                                 return.list$similarities.all <- left_join(tmp, similarities.long.df.all, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
                                 
                                 return.list$cor <- left_join(tmp, protein.matrix.mean.cor.long, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
                                 return.list$cor_fc <- left_join(tmp, proteins.FC.f.matrix.cor.long, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
                                 
                                 return(bind_rows(return.list))
                               })

pathway_similarities$sample_type = "signal" 


all_similarities <- bind_rows(similarities.long.df.all, similarities.long.df.down, similarities.long.df.up, 
                              proteins.FC.f.matrix.cor.long, protein.matrix.mean.cor.long) %>% filter(!is.na(value))


pathway_similarities.grouped <- bind_rows(pathway_similarities %>% filter(n >= 2) %>% mutate(n_cut = ">=2"),
                                          pathway_similarities %>% filter(n >= 3) %>% mutate(n_cut = ">=3"),
                                          pathway_similarities %>% filter(n >= 5) %>% mutate(n_cut = ">=5"))

set.seed(1234)
random_similarities.grouped <- ddply(pathway_similarities.grouped %>% filter(n_cut == ">=3"), 
      .(pathway_base, sim_type, n_cut), 
      .fun = function(x) {
        #x = pathway_similarities.grouped %>% filter(pathway_base == "kegg", sim_type == "overlap.up", n_cut == ">=3")
        universe = unique(c(x$X1,x$X2))
        ret = all_similarities %>% filter(X1_uniprot %in% universe, 
                                          X2_uniprot %in% universe,
                                          sim_type == unique(x$sim_type)) %>% 
                                    sample_n(1000, replace = T) %>%
                                    mutate(pathway = "none", sample_type = "random")
        return(ret)
      })
  
random_similarities.grouped$sample_type = "random" 

pathway_similarities.grouped.dataset <- bind_rows(pathway_similarities.grouped, random_similarities.grouped)

file_name = paste("pathway_similarities.grouped.dataset", fun_name, "RData", sep=".")
file_path = paste(output_dir, file_name, sep="/")
save(pathway_similarities.grouped.dataset, file=file_path)



toPlot <- pathway_similarities.grouped.dataset %>% filter(sim_type != "pearson_fc",
                                                          sim_type != "overlap.all",
                                                          n_cut == ">=3")
toPlot$sample_type <- factor(toPlot$sample_type, levels = c("signal", "random"))
toPlot$sim_type <- as.factor(toPlot$sim_type)

toPlot.stats <- toPlot %>% group_by(sim_type, pathway_base, n_cut) %>% 
  summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"])$'p.value'))

toPlot.stats.medians <- toPlot %>% 
  group_by(sim_type, pathway_base, n_cut, sample_type) %>% 
  summarise(median_value = median(value, na.rm = T))

toPlot.stats$padj <- p.adjust(toPlot.stats$pval, method = "BH")

ggplot() +
  geom_density(data = toPlot, aes(x = value, fill = sample_type), alpha = 0.5) +
  facet_grid(pathway_base+n_cut~sim_type, scales = "free") +
  geom_text(data = toPlot.stats, aes(x=0.5, y = 5, label= paste("p-value=", format(padj, digits=2, scientific=T)))) +
  geom_vline(data = toPlot.stats.medians, aes(xintercept = median_value, colour = sample_type), linetype = 2) +
  theme_bw() + 
  theme(legend.position = c(0.1, 0.5))

# 
# 
# 
# reactome_similarities <- ddply(UniProt2Reactome.kinases, .(reactome_id), 
#       .fun = function(x) {
#         return.list = list()
#         
#         x <- UniProt2Reactome.kinases[UniProt2Reactome.kinases$pathway == "R-SCE-2559580",]
#         
#         tmp <- data.frame(t(combn(unique(x$uniprot_id), 2)))
#         
#         return.list$similarities.up <- left_join(tmp, similarities.long.df.up, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
#         return.list$similarities.down <- left_join(tmp, similarities.long.df.down, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
#         return.list$similarities.all <- left_join(tmp, similarities.long.df.all, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
#         
#         return.list$cor <- left_join(tmp, protein.matrix.mean.cor.long, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
#         return.list$cor_fc <- left_join(tmp, proteins.FC.f.matrix.cor.long, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
#         
#         return(bind_rows(return.list))
#       })
# 
# reactome_similarities <- tbl_df(reactome_similarities)
# #reactome_similarities <- reactome_similarities %>% group_by(sim_type) %>% distinct(X1, X2) %>% ungroup()
# reactome_similarities$sample_type = "signal" 
# 
# 
# 
# KEGG_similarities <- ddply(orf2kegg.kinases, .(pathway), 
#                                .fun = function(x) {
#                                  return.list = list()
#                                  
#                                  x <- orf2kegg.kinases[orf2kegg.kinases$pathway == "path:sce04011",]
#                                  
#                                  tmp <- data.frame(t(combn(unique(as.character(x$ORF)), 2)))
#                                  
#                                  return.list$similarities.up <- left_join(tmp, similarities.long.df.up, by = c("X1" = "X1", "X2" = "X2"))
#                                  return.list$similarities.down <- left_join(tmp, similarities.long.df.down, by = c("X1" = "X1", "X2" = "X2"))
#                                  return.list$similarities.all <- left_join(tmp, similarities.long.df.all, by = c("X1" = "X1", "X2" = "X2"))
#                                  
#                                  return.list$cor <- left_join(tmp, protein.matrix.mean.cor.long, by = c("X1" = "X1", "X2" = "X2"))
#                                  return.list$cor_fc <- left_join(tmp, proteins.FC.f.matrix.cor.long, by = c("X1" = "X1", "X2" = "X2"))
#                                  
#                                  return(bind_rows(return.list))
#                                })
# 
# #KEGG_similarities <- KEGG_similarities %>% group_by(sim_type) %>% distinct(X1, X2) %>% ungroup()
# KEGG_similarities$sample_type = "signal" 
# 
# mapped_kinases.reactome <- as.vector(na.omit(unique(c(reactome_similarities$X1 , reactome_similarities$X2 ))))
# 
# set.seed(123)
# random_similarities.reactome <- bind_rows(similarities.long.df.all, similarities.long.df.down, 
#                                           similarities.long.df.up, proteins.FC.f.matrix.cor.long, protein.matrix.mean.cor.long ) %>% 
#   filter(!is.na(value)) %>% filter(X1_uniprot %in% mapped_kinases.reactome, X2_uniprot %in% mapped_kinases.reactome) %>%
#   group_by(dataset, sim_type, type) %>% 
#   sample_n(1000, replace=T) %>% mutate(reactome_id = "none", sample_type = "random")
# 
# reactome_similarities <- rename(reactome_similarities, X1_uniprot = X1, X2_uniprot = X2, X1 =  X1, X2 = X2) 
# 
# 
# similarities.dataset.reactome <- bind_rows(reactome_similarities, random_similarities.reactome) %>% filter(!is.na(sim_type))
# 
# file_name = paste("similarities.dataset.reactome", fun_name, "RData", sep=".")
# file_path = paste(output_dir, file_name, sep="/")
# save(similarities.dataset.reactome, file=file_path)
# 
# 
# 
# load("./R/objects/similarities.dataset.reactome.signaling_similarity.RData")
# toPlot <- similarities.dataset.reactome
# toPlot$sample_type <- factor(toPlot$sample_type, levels = c("signal", "random"))
# 
# toPlot$sim_type <- as.factor(toPlot$sim_type)
# 
# toPlot.stats <- toPlot %>% group_by(sim_type) %>% 
#   summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"], alternative = "greater")$'p.value'))
# 
# ggplot(toPlot , aes(x = value, fill = sample_type)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(breaks = seq(0,1, by = 0.25)) +
#   facet_wrap(~sim_type, scales = "free") +
#   theme_bw() + 
#   theme(legend.position = c(0.1, 0.5))
# 
# ###
# 
# mapped_kinases.kegg <- as.vector(na.omit(unique(c(KEGG_similarities$X1 , KEGG_similarities$X2 ))))
# set.seed(123)
# random_similarities.kegg <- bind_rows(similarities.long.df, proteins.FC.f.matrix.cor.long, protein.matrix.mean.cor.long ) %>% 
#   filter(!is.na(value)) %>% filter(X1 %in% mapped_kinases.kegg, X2 %in% mapped_kinases.kegg) %>%
#   group_by(dataset, sim_type, type) %>% 
#   sample_n(1000, replace=T) %>% mutate(pathway = "none", sample_type = "random")
# 
# #reactome_similarities <- rename(reactome_similarities, X1_uniprot = X1, X2_uniprot = X2, X1 =  X1.y, X2 = X2.y) 
# 
# similarities.dataset.kegg <- bind_rows(KEGG_similarities, random_similarities.kegg) %>% filter(!is.na(sim_type))
# toPlot <- similarities.dataset.kegg
# toPlot$sample_type <- factor(toPlot$sample_type, levels = c("signal", "random"))
# 
# toPlot$sim_type <- as.factor(toPlot$sim_type)
# 
# toPlot.stats <- toPlot %>% group_by(sim_type) %>% 
#   summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"], alternative = "greater")$'p.value'))
# 
# ggplot(toPlot , aes(x = value, fill = sample_type)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(breaks = seq(0,1, by = 0.25)) +
#   facet_wrap(~sim_type, scales = "free") +
#   theme_bw() + 
#   theme(legend.position = c(0.1, 0.5))
# 
#   
# str(similarities.dataset.kegg)
# similarities.dataset.reactome <- rename(similarities.dataset.reactome, pathway = reactome_id )
#                                  
# similarities.dataset.reactome$pathway_base <- "reactome"                                 
# similarities.dataset.kegg$pathway_base <- "kegg"                                 
# 
# similarities.dataset <- bind_rows(similarities.dataset.reactome, similarities.dataset.kegg)
# 
# file_name = paste("similarities.dataset", fun_name, "RData", sep=".")
# file_path = paste(output_dir, file_name, sep="/")
# save(similarities.dataset, file=file_path)
# 
