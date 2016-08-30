rm(list=ls())
source("./R/boot.R")
source("./R/functions.R")
library(cowplot)


output_dir = "./R/objects"
figures_dir = "./figures"
fun_name = "signaling_similarity"
dir.create(output_dir)
dir.create(figures_dir)
plots.list = list()

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
  filter(n >2) %>% arrange(reactome_id)

UniProt2Reactome.kinases <- dplyr::rename(UniProt2Reactome.kinases, pathway = reactome_id)


load("./R/objects/pathway2orf._load_.RData")
pathway2orf$uniprot_id = uniprot2orf.kinases$uniprot_id[match(pathway2orf$ORF, uniprot2orf.kinases$ORF)]

View(pathway2orf)

orf2kegg.kinases <- pathway2orf %>% 
  dplyr::filter(!is.na(uniprot_id)) %>%
  group_by(pathway) %>%
  mutate(n = n()) %>% 
  filter(n >2) %>% arrange(pathway) %>% 
  dplyr::select(pathway, uniprot_id, n)



orf2kegg.kinases$pathway_base = "kegg"
UniProt2Reactome.kinases$pathway_base = "reactome"

pathway_kinases <- bind_rows(orf2kegg.kinases, UniProt2Reactome.kinases)



set.seed(123)
pathway_kinases_random <- ddply(pathway_kinases, 
                                .(pathway_base),
                                .fun = function(x) {
                                  # x = pathway_kinases %>% filter(pathway_base == "kegg")
                                  uniprots = unique((pathway_kinases %>% filter(pathway_base == unique(x$pathway_base)) %>% dplyr::select(uniprot_id))$uniprot_id)
                                  z <<- ddply(x, 
                                        .(pathway), 
                                        .fun = function(y) {
                          #                 y = pathway_kinases %>% filter(pathway == "path:sce04011",
                          #                                                pathway_base == "kegg")
                                          
                                          n = unique(y$n)
                                          tmp = as.matrix(replicate(1000, sample(uniprots, n)))
                                          colnames(tmp) = paste0("rep_", 1:ncol(tmp), sep="")
                                          tmp.long <- melt(tmp, varnames = c("X1", "X2"))
                                          tmp.long$n <- n
                                          return(tmp.long)
                                        }
                                        )
                                  
                                  return(data.frame(pathway = paste0(z$pathway,"_",  z$X2),
                                                   uniprot_id = z$value,
                                                   n = z$n))
                                })


getPathwaySimilarity <- function(pathways){
  ret = ddply(pathways, 
        .(pathway, pathway_base), 
        .fun = function(x) {
          return.list = list()
#           z <<- x
#           x <- z
          tmp <<- data.frame(t(combn(unique(as.character(x$uniprot_id)), 2)))
          tmp$n = x$n[1]
          return.list$similarities.up <- left_join(tmp, similarities.long.df.up, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
          return.list$similarities.down <- left_join(tmp, similarities.long.df.down, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
          return.list$similarities.all <- left_join(tmp, similarities.long.df.all, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
          
          return.list$cor <- left_join(tmp, protein.matrix.mean.cor.long, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
          return.list$cor_fc <- left_join(tmp, proteins.FC.f.matrix.cor.long, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
          
          return(bind_rows(return.list))
        })
  return(ret)
}

getRandomRepeatStats <- function(pathway_similarities, random_similarities) {
  random_repeats.stats <- ddply(random_similarities,
                         .(pathway_base, sim_type, rep),
                         .fun = function(x) {
                           pathway_similarities %>% 
                             filter(sim_type ==  x$sim_type[1], pathway_base == x$pathway_base[1]) %>%
                             summarise(pval = wilcox.test(value, x$value)$'p.value',
                                       rep = x$rep[1],
                                       diff_median = abs(median(value, na.rm = T) - median(x$value, na.rm = T)))
                         })
  
  random_repeats.stats$padj <- p.adjust(random_repeats.stats$pval,  method = "BH")
  return(random_repeats.stats)
}
 
pathway_similarities <- getPathwaySimilarity(pathways = pathway_kinases)
pathway_similarities$sample_type = "signal" 
pathway_similarities_random <- getPathwaySimilarity(pathways = pathway_kinases_random)
pathway_similarities_random$sample_type = "random"

pathways_dataset <- bind_rows(pathway_similarities, pathway_similarities_random)

pathway_similarities_random.repeats <- pathway_similarities_random %>% separate(pathway, into=c("pahtway_id", "rep"), sep="_rep_" , remove = F)

makeReport <- function(random_repeats.stats, pathway_similarities, pathway_similarities_random.repeats, sim_type_) {
  pvalx=0.5
  pvaly= 5
  #produces plots double panel 
  pval_thr = 0.05
  toPlot <- random_repeats.stats
  
  toPlot.long <- toPlot %>%
    dplyr::select(-pval) %>%
    melt(id.vars = c("pathway_base", "sim_type", "rep"))
  
  toPlot.stats <- toPlot %>% filter(sim_type != "pearson_fc") %>%
    group_by(pathway_base, sim_type) %>%
    summarise(padj = sum(padj<pval_thr)/length(padj),
              pval = sum(pval<pval_thr)/length(pval))
  
  toPlot.stats$diff_median = NA

 #text for pvalues 
 toPlot.long.stats <- toPlot.long %>% filter(grepl(pattern = "overlap", sim_type)) %>%
    group_by(pathway_base, sim_type) %>%
    summarise(diff_median = median(value[variable == "diff_median"], na.rm = T),
              padj = pval_thr) %>% 
    gather(variable, value, -pathway_base, -sim_type)
    
  toPlot.stats.text <- toPlot.stats %>% dplyr::select(-pval) %>% 
    filter(grepl(pattern = "overlap", sim_type)) %>%
    gather(variable, value, -pathway_base, -sim_type)
    
  histogram_data <- toPlot.long %>% filter(grepl(pattern = "overlap", sim_type))
  
  if (sim_type_ == "pearson") {
    # pearson pvalues
    toPlot.long.stats <- toPlot.long %>%
      group_by(pathway_base, sim_type) %>%
      summarise(diff_median = median(value[variable == "diff_median"], na.rm = T),
                padj = pval_thr) %>%
      gather(variable, value, -pathway_base, -sim_type) %>%
      filter(variable != "diff_median", sim_type == "pearson")  
    toPlot.stats.text <- toPlot.stats %>% dplyr::select(-pval) %>%
      gather(variable, value, -pathway_base, -sim_type) %>%
      filter(variable != "diff_median", sim_type == "pearson")
    histogram_data <- toPlot.long %>% filter(sim_type == "pearson", variable != "diff_median")
    
  }
  
  pvalues_panel <- ggplot() +
    geom_histogram(data = histogram_data, aes(x= value), bins = 15, colour = "white", fill = "black") +
    geom_vline(data = toPlot.long.stats, aes(xintercept=value), linetype = 2, colour = "red") +
    geom_text(data = toPlot.stats.text, aes(x = pval_thr, y = 100, label = value)) +
    facet_grid(pathway_base+sim_type~variable, scales = "free") +
    theme_bw() +
    theme(aspect.ratio = 5/8)
  
  # single random example of overlaps
  set.seed(123)
  toPlot.single <- bind_rows(pathway_similarities, pathway_similarities_random.repeats %>% group_by(rep) %>% sample_n(1)) %>% 
    filter(sim_type != "overlap.all", grepl("overlap", sim_type) )
  
  if (sim_type_ == "pearson") {
    toPlot.single <- bind_rows(pathway_similarities, pathway_similarities_random.repeats %>% group_by(rep) %>% sample_n(1)) %>% 
      filter(sim_type == "pearson")
    pvalx=0.8
    pvaly= 5
  }
  
  toPlot.single$sample_type <- factor(toPlot.single$sample_type, levels = c("signal", "random"))
  toPlot.single$sim_type <- as.factor(toPlot.single$sim_type)
  
  toPlot.stats <- toPlot.single %>% group_by(sim_type, pathway_base) %>% 
    summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"])$'p.value'))
  
  toPlot.stats$padj <- p.adjust(toPlot.stats$pval, method = "BH")
  
  toPlot.stats.medians <- toPlot.single %>% 
    group_by(sim_type, pathway_base, sample_type) %>% 
    summarise(median_value = median(value, na.rm = T))
  
  single_panel <- ggplot() +
    geom_density(data = toPlot.single, aes(x = value, fill = sample_type), alpha = 0.5) +
    facet_grid(pathway_base+sim_type~., scales = "free") +
    geom_text(data = toPlot.stats, aes(x=pvalx, y = pvaly, label= paste("p-value=", format(pval, digits=2, scientific=T)))) +
    geom_vline(data = toPlot.stats.medians, aes(xintercept = median_value, colour = sample_type), linetype = 2) +
    theme_bw() + 
    theme(legend.position = c(0.1, 0.5), aspect.ratio = 5/8)
  return(plot_grid(single_panel, pvalues_panel, labels = c("A", "B")))
}



## --- co-occurrence of kinases in pathways ----

pathway_similarities$sorted_pair <- with(pathway_similarities, paste(sort(c(X1, X2)), collapse = "|"))
pathway_similarities$sorted_pair <- apply(pathway_similarities[ ,c("X1", "X2")] , 1 , 
                                          FUN = function(x) {
                                            #x = pathway_similarities[1 ,c("X1", "X2")]
                                            paste(sort(c(x[1], x[2])), collapse = "|")
                                          })

pathway_similarities %>% 
  group_by(pathway_base, sim_type, sample_type, sorted_pair) %>%
  mutate(member_of_n = n()) %>% 
  ungroup() %>% filter(sim_type == "overlap.up") %>%
  arrange(-member_of_n)

# "P12688|P18961" pair found in multiple pathways in reactome


## making tests of individual random shuffles repeats
pathway_similarities_random.repeats.stats <- getRandomRepeatStats(pathway_similarities = pathway_similarities %>% 
                                                                    filter(sim_type != "overlap.all"),
                                                                  random_similarities = pathway_similarities_random.repeats %>% 
                                                                    filter(sim_type != "overlap.all"))

# removing pair from intdividual tests
pathway_similarities_random.repeats.stats_nopair <- getRandomRepeatStats(pathway_similarities = pathway_similarities %>%
                                                                           filter(sim_type != "overlap.all" ,
                                                                                  !(X1 %in% c("P12688", "P18961")), 
                                                                                  !(X2 %in% c("P12688", "P18961"))),
                                                                         random_similarities =  pathway_similarities_random.repeats %>% 
                                                                           filter(sim_type != "overlap.all" ,
                                                                                  !(X1 %in% c("P12688", "P18961")), 
                                                                                  !(X2 %in% c("P12688", "P18961"))))

## -- saving reports ----
p <- makeReport(random_repeats.stats = pathway_similarities_random.repeats.stats, 
           pathway_similarities = pathway_similarities,
           pathway_similarities_random.repeats = pathway_similarities_random.repeats,
           sim_type_ =  "overlaps" )
plots.list <- lappend(plots.list, p)

p <- makeReport(random_repeats.stats = pathway_similarities_random.repeats.stats, 
           pathway_similarities = pathway_similarities,
           pathway_similarities_random.repeats = pathway_similarities_random.repeats,
           sim_type = "pearson" )

plots.list <- lappend(plots.list, p)


# without frequent pair
p <- makeReport(random_repeats.stats = pathway_similarities_random.repeats.stats_nopair, 
           pathway_similarities = pathway_similarities %>% 
             filter(sim_type != "overlap.all" ,
                    !(X1 %in% c("P12688", "P18961")), 
                    !(X2 %in% c("P12688", "P18961"))),
           pathway_similarities_random.repeats = pathway_similarities_random.repeats  %>% 
             filter(sim_type != "overlap.all",
                    !(X1 %in% c("P12688", "P18961")), 
                    !(X2 %in% c("P12688", "P18961"))),
                    sim_type = "overlaps")
plots.list <- lappend(plots.list, p)

p <- makeReport(random_repeats.stats = pathway_similarities_random.repeats.stats_nopair, 
           pathway_similarities = pathway_similarities %>% 
             filter(sim_type != "overlap.all" ,
                    !(X1 %in% c("P12688", "P18961")), 
                    !(X2 %in% c("P12688", "P18961"))),
           pathway_similarities_random.repeats = pathway_similarities_random.repeats  %>% 
             filter(sim_type != "overlap.all",
                    !(X1 %in% c("P12688", "P18961")), 
                    !(X2 %in% c("P12688", "P18961"))),
           sim_type = "pearson")
plots.list <- lappend(plots.list, p)

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
           
           ggplot2::ggsave(filename = paste(file_path, x , "png", sep = "."), device = NULL,
                           plot = p, width = 210 , height = 297, dpi = 150, units = "mm", scale = scale)
           
         }, error = function(e) {
           message(paste("Plot", "x", "sucks!" ))
           return(NULL)
         }, finally = {
           message(paste("processed plot", x))
         })
         
         
       })




# 
# # single random example of overlaps
# set.seed(123)
# toPlot <- bind_rows(pathway_similarities, pathway_similarities_random.repeats %>% group_by(rep) %>% sample_n(1)) %>% 
#   filter(sim_type != "overlap.all", grepl("overlap", sim_type) )
# toPlot$sample_type <- factor(toPlot$sample_type, levels = c("signal", "random"))
# toPlot$sim_type <- as.factor(toPlot$sim_type)
# 
# toPlot.stats <- toPlot %>% group_by(sim_type, pathway_base) %>% 
#   summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"])$'p.value'))
# 
# toPlot.stats.medians <- toPlot %>% 
#   group_by(sim_type, pathway_base, sample_type) %>% 
#   summarise(median_value = median(value, na.rm = T))
# 
# toPlot.stats$padj <- p.adjust(toPlot.stats$pval, method = "BH")
# 
# p.overlaps_single <- ggplot() +
#   geom_density(data = toPlot, aes(x = value, fill = sample_type), alpha = 0.5) +
#   facet_grid(pathway_base+sim_type~., scales = "free") +
#   geom_text(data = toPlot.stats, aes(x=0.5, y = 5, label= paste("p-value=", format(pval, digits=2, scientific=T)))) +
#   geom_vline(data = toPlot.stats.medians, aes(xintercept = median_value, colour = sample_type), linetype = 2) +
#   theme_bw() + 
#   theme(legend.position = c(0.1, 0.5))
# 
# p.grid.overlaps <- plot_grid(p.overlaps_single, p.double_overlap, labels = c("A", "B"))
# 
# #single random example of pearson
# set.seed(123)
# toPlot <- bind_rows(pathway_similarities, pathway_similarities_random.repeats %>% group_by(rep) %>% sample_n(1)) %>% 
#   filter(sim_type == "pearson")
# toPlot$sample_type <- factor(toPlot$sample_type, levels = c("signal", "random"))
# toPlot$sim_type <- as.factor(toPlot$sim_type)
# 
# toPlot.stats <- toPlot %>% group_by(sim_type, pathway_base) %>% 
#   summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"])$'p.value'))
# 
# toPlot.stats.medians <- toPlot %>% 
#   group_by(sim_type, pathway_base, sample_type) %>% 
#   summarise(median_value = median(value, na.rm = T))
# 
# toPlot.stats$padj <- p.adjust(toPlot.stats$pval, method = "BH")
# 
# p.pearson_single <- ggplot() +
#   geom_density(data = toPlot, aes(x = value, fill = sample_type), alpha = 0.5) +
#   facet_grid(pathway_base+sim_type~., scales = "free") +
#   geom_text(data = toPlot.stats, aes(x=0.9, y = 5, label= paste("p-value=", format(pval, digits=2, scientific=T)))) +
#   geom_vline(data = toPlot.stats.medians, aes(xintercept = median_value, colour = sample_type), linetype = 2) +
#   theme_bw() + 
#   theme(legend.position = c(0.1, 0.5))
# p.grid.pearson <- plot_grid(p.pearson_single, p.double_pearson, labels = c("A", "B"))
# 
# 
# #
# 
# pathway_kinases %>% group_by(pathway_base) %>% distinct(pathway) %>% summarise(n = n())
# 
# unique(pathway_kinases[pathway_kinases$pathway_base == "kegg",]$pathway)
# 
# # all_similarities <- bind_rows(similarities.long.df.all, similarities.long.df.down, similarities.long.df.up, 
# #                               proteins.FC.f.matrix.cor.long, protein.matrix.mean.cor.long) %>% filter(!is.na(value))
# # 
# # pathway_similarities.grouped <- bind_rows(pathway_similarities %>% filter(n >= 2) %>% mutate(n_cut = ">=2"),
# #                                           pathway_similarities %>% filter(n >= 3) %>% mutate(n_cut = ">=3"),
# #                                           pathway_similarities %>% filter(n >= 5) %>% mutate(n_cut = ">=5"))
# # 
# # set.seed(1234)
# # random_similarities.grouped <- ddply(pathway_similarities.grouped %>% filter(n_cut == ">=3"), 
# #       .(pathway_base, sim_type, n_cut), 
# #       .fun = function(x) {
# #         #x = pathway_similarities.grouped %>% filter(pathway_base == "kegg", sim_type == "overlap.up", n_cut == ">=3")
# #         universe = unique(c(x$X1,x$X2))
# #         ret = all_similarities %>% filter(X1_uniprot %in% universe, 
# #                                           X2_uniprot %in% universe,
# #                                           sim_type == unique(x$sim_type)) %>% 
# #                                     sample_n(1000, replace = T) %>%
# #                                     mutate(pathway = "none", sample_type = "random")
# #         return(ret)
# #       })
# #   
# # random_similarities.grouped$sample_type = "random" 
# # pathway_similarities.grouped.dataset <- bind_rows(pathway_similarities.grouped, random_similarities.grouped)
# # 
# # file_name = paste("pathway_similarities.grouped.dataset", fun_name, "RData", sep=".")
# # file_path = paste(output_dir, file_name, sep="/")
# # save(pathway_similarities.grouped.dataset, file=file_path)
# # 
# # 
# # 
# # 
# # 
# # toPlot <- pathway_similarities.grouped.dataset %>% filter(sim_type != "pearson_fc",
# #                                                           sim_type != "overlap.all",
# #                                                           n_cut == ">=3")
# # toPlot$sample_type <- factor(toPlot$sample_type, levels = c("signal", "random"))
# # toPlot$sim_type <- as.factor(toPlot$sim_type)
# # 
# # toPlot.stats <- toPlot %>% group_by(sim_type, pathway_base, n_cut) %>% 
# #   summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"])$'p.value'))
# # 
# # toPlot.stats.medians <- toPlot %>% 
# #   group_by(sim_type, pathway_base, n_cut, sample_type) %>% 
# #   summarise(median_value = median(value, na.rm = T))
# # 
# # toPlot.stats$padj <- p.adjust(toPlot.stats$pval, method = "BH")
# # 
# # ggplot() +
# #   geom_density(data = toPlot, aes(x = value, fill = sample_type), alpha = 0.5) +
# #   facet_grid(pathway_base+n_cut~sim_type, scales = "free") +
# #   geom_text(data = toPlot.stats, aes(x=0.5, y = 5, label= paste("p-value=", format(padj, digits=2, scientific=T)))) +
# #   geom_vline(data = toPlot.stats.medians, aes(xintercept = median_value, colour = sample_type), linetype = 2) +
# #   theme_bw() + 
# #   theme(legend.position = c(0.1, 0.5))
