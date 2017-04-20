rm(list = ls())
library(tidyverse)
library(forcats)
library(viridis)
library(ggthemes)
library(forcats)
library(stringr)
library(grid)

fun_name = "MCA"
figures_dir = "./figures"

lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

plots.list = list()
load("./R/objects/exp_metadata._clean_.RData")

## -- misc settings ------
#don't format the following vectors!!!
parameter_order <- c('HXK2
HXK1
GLK1
PGM1
PGI1
PFK1
PFK2
FBA1
GPD1
GPD2
HOR2
PGK1
GPM1
ENO1
ENO2
CDC19
PDC1
PDC5
PDC6
ADH1
ADH5') %>% str_split("\n") %>% unlist()

reaction_order <- c('HXT
HXK_HXK2
HXK_HXK1
HXK_GLK1
PGM
TPS
TPP
PGI
PFK
FBA
TPI
GPD
GPP
TDH_TDH1
TDH_TDH3
PGK
GPM
ENO_ENO1
ENO_ENO2
PYK_CDC19
ATPase
PDC_PDC1
PDC_PDC5
PDC_PDC6
ADH_ADH1
acetate_branch
UGP
udp_to_utp
AK') %>% str_split("\n") %>% unlist()

metabolite_order <- c('GLC
G6P
F6P
F16bP
GAP
DHAP
BPG
P3G
P2G
PEP
PYR
AcAld
G1P
G3P
UDP
UTP
T6P
ADP
ATP
NAD') %>% str_split("\n") %>% unlist()


## -- Steady-concentration fluxes ----
dataset.ss <- read_delim("./results/2017-02-22/ss_fluxes_concentrations_in_kinase_KOs.csv", delim = "\t")
WT.ss <- dataset.ss %>% filter(KO == "WT")

dataset.ss <- left_join(dataset.ss, WT.ss, by = c("type" = "type", "variable" = "variable"))  %>% 
  mutate(value_change = (value.x - value.y)/value.y) %>% dplyr::select(-matches(".y$")) 
names(dataset.ss) <- sub(x = names(dataset.ss), pattern = ".x$", replacement = "")


dataset.ss <- dataset.ss %>% group_by(variable, type) %>% 
  filter(ss_status < 10-6) %>% mutate(z_value = (value - mean(value, na.rm=T))/sd(value, na.rm=T))


p.ss <- ggplot(data = dataset.ss, aes(x = variable, y = z_value)) +
  geom_boxplot() +
  facet_wrap(~type, scales = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 5/8)

library(scales)
p.ss_changes <- ggplot(data = dataset.ss, aes(x = variable, y = value_change)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-3,3), labels=percent) +
  labs(y="Value change compared to WT, %") +
  facet_wrap(~type, scales = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 5/8)
plots.list = lappend(plots.list, p.ss_changes)

dataset.ss <- dataset.ss %>% mutate(abs_value_change = ifelse(abs(value_change)>3,3, abs(value_change)))

p.ss_changes_density <- ggplot(data = dataset.ss, aes(x = abs_value_change)) +
  geom_density(fill = "black") +
  scale_x_continuous(limits = c(0,3), labels=percent) +
  facet_wrap(~type, scales = "free") +
  labs(x = "Absolute value change compared to WT, %") +
  theme_bw(base_family = "Helvetica") +
  theme(aspect.ratio = 5/8) +
  theme(panel.grid = element_blank())
plots.list = lappend(plots.list, p.ss_changes_density)

p.ss_density <- ggplot(data = dataset.ss %>% filter(type == "flux"), aes(x = value)) +
  geom_density() +
  facet_wrap(~variable, scales = "free") +
  theme_classic()

p.ss_density$landscape = T
p.ss_density$toScale = T
plots.list = lappend(plots.list, p.ss_density)


p.ss_density_conc <- ggplot(data = dataset.ss %>% filter(type == "conc"), aes(x = value)) +
  geom_density() +
  facet_wrap(~variable, scales = "free") +
  theme_classic()

p.ss_density_conc$landscape = T
p.ss_density_conc$toScale = T
plots.list = lappend(plots.list, p.ss_density_conc)



#small helper to tidy data for pca
tidy_pca = function(x) {
  x.matrix <- x[,-1] %>% as.matrix()
  rownames(x.matrix) <- as.data.frame(x)[,1]
  return(x.matrix)
}

flux_matrix <- dataset.ss %>% ungroup() %>%
  filter(type == "flux", ss_status < 10-6, is.nan(z_value) == F | is.na(z_value) == F | is.infinite(z_value) == F) %>% 
  dplyr::select(KO, variable, z_value) %>% filter(!is.nan(z_value), !is.na(z_value)) %>%
  spread(variable, z_value) %>% tidy_pca 

fluxes_pca <- prcomp(flux_matrix)

conc_matrix <- dataset.ss %>% ungroup() %>%
  filter(type == "conc", is.nan(z_value) == F ) %>% 
  dplyr::select(KO, variable, z_value) %>% filter(!is.nan(z_value), !is.na(z_value)) %>%
  spread(variable, z_value) %>% tidy_pca()
  
conc_pca <- conc_matrix %>%  prcomp()

conc_loadings <- conc_pca %>% .$rotation %>% reshape2::melt()

toPlot <- conc_pca
xlabel <- paste("PC1", round(toPlot$sdev[1]/sum(toPlot$sdev),2))
ylabel <- paste("PC2", round(toPlot$sdev[2]/sum(toPlot$sdev),2))

p.scores_conc <- toPlot$x %>% as_tibble() %>% mutate(gene_name = exp_metadata$gene[match(rownames(toPlot$x), exp_metadata$ORF)]) %>%
  ggplot(aes(x = PC1 , y = PC2)) +
    geom_point() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_text(aes(label = gene_name)) +
    theme(aspect.ratio = 1) +
    xlab(label = xlabel) +
    ylab(label = ylabel)

toPlot <- conc_pca %>% .$rotation %>% reshape2::melt()
names(toPlot) = c("variable", "component", "value")
toPlot <- toPlot %>% filter(component %in% c("PC1", "PC2")) %>% mutate(variable = fct_reorder(variable, value))
p.loadings_conc <-  ggplot(toPlot, aes(x = variable, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Steady-state concentrations")

p.grid_conc <- plot_grid(p.scores_conc, p.loadings_conc, labels = c("A", "B"), align = "h" )
p.grid_conc$landscape = T
plots.list = lappend(plots.list, p.grid_conc)

toPlot <- fluxes_pca
xlabel <- paste("PC1", round(toPlot$sdev[1]/sum(toPlot$sdev),2))
ylabel <- paste("PC2", round(toPlot$sdev[2]/sum(toPlot$sdev),2))

p.scores_fluxes <- toPlot$x %>% as_tibble() %>% mutate(gene_name = exp_metadata$gene[match(rownames(toPlot$x), exp_metadata$ORF)]) %>%
  ggplot(aes(x = PC1 , y = PC2)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = gene_name)) +
  theme(aspect.ratio = 1) +
  xlab(label = xlabel) +
  ylab(label = ylabel)

#loadings
toPlot = fluxes_pca %>% .$rotation %>% reshape2::melt()
names(toPlot) = c("variable", "component", "value")

p.loadings_fluxes <- toPlot %>% filter(component %in% c("PC1", "PC2")) %>% 
  mutate(variable = fct_reorder(variable, value)) %>%
  ggplot(aes(x = variable, y = value, fill = component)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Steady-state fluxes")

p.grid_fluxes <- plot_grid(p.scores_fluxes, p.loadings_fluxes, labels = c("A", "B"), align = "h")
p.grid_fluxes$landscape = T
plots.list = lappend(plots.list, p.grid_fluxes)


## -- MCA for fluxes/concentrations ----

flux_mca <- read_delim("./results/2017-02-22/mca_fluxes_in_kinase_KOs.csv", delim = "\t")
conc_mca <- read_delim("./results/2017-02-22/mca_conc_in_kinase_KOs.csv", delim = "\t")

dataset.mca <- bind_rows(flux_mca, conc_mca) %>% dplyr::rename(VAR = variable) %>% group_by(var_type, KO, VAR, parameter) %>% distinct(.keep_all = T) %>% ungroup() ## unclear why multiple metabolites???? 
dataset.mca$KO_name <- as.character(exp_metadata$gene[match(dataset.mca$KO, exp_metadata$ORF)])


toPlot <- dataset.mca %>% group_by(var_type, VAR, KO) %>% 
  summarise(min = min(CC, na.rm = T),
            min_plot = ifelse(min < -3, -3, min),
            min_parameter = parameter[which.min(CC)][1],
            max = max(CC, na.rm = T),
            max_plot = ifelse(max > 3, 3, max),
            max_parameter = parameter[which.max(CC)][1],
            max_abs = max(abs(CC), na.rm = T),
            max_abs_parameter = parameter[which.max(abs(CC))][1]) %>% 
  dplyr::select(var_type, VAR, KO, max_plot, min_plot ) %>% as.data.frame() %>%
  reshape2::melt(id.vars = c("var_type", "VAR", "KO"))

toPlot %>% 
p.min_max <- ggplot(aes(x = VAR, y = value)) +
    ylim(c(-3,3)) +
    geom_violin() +
    geom_jitter(aes(colour = variable), alpha = 0.3) +
    geom_point(data = toPlot %>% filter( KO == "WT")) +
    facet_wrap(~var_type, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
plots.list = lappend(plots.list, p.min_max)

selected_var = "ADH_ADH1"
toPlot <- dataset.mca %>% ungroup() %>%
  filter(VAR == selected_var, var_type == "flux") %>%
  mutate(parameter = factor(parameter, levels = parameter_order))
p.example <- ggplot(toPlot,  aes(x = parameter, y = CC)) +
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot() +
    #geom_jitter() +
    ylab("Flux control coeffiecient") +
    geom_point(data = toPlot %>% filter(KO == "WT"), aes(colour = "red")) +
    theme_bw() +
    theme(legend.position = "none")

plots.list = lappend(plots.list, p.example)

p.example_zoom <- ggplot(toPlot,  aes(x = parameter, y = CC)) +
  ylab(paste("Flux control coeficients of", selected_var)) +
  geom_jitter(alpha = 0.25) +
  ylim(0, 1) +
  geom_point(data = toPlot %>% filter(KO == "WT"), aes(colour = "red")) +
  #ggrepel::geom_text_repel(data = toPlot %>% filter(parameter %in% c("HXK2", "PFK1")), aes(label = KO_name, x = parameter, y = CC), segment.alpha = 0.25, segment.size = 0.25) +
  theme_bw() +
  theme(legend.position="none")


plots.list = lappend(plots.list, p.example_zoom)
 
#overall flux 
# toPlot <- dataset.mca %>% group_by(var_type, parameter, KO, KO_name) %>% filter(var_type == "flux", VAR != "AK") %>% summarise(overal_CC = sqrt(sum(CC^2, na.rm = T)))
# toPlot.examples <- toPlot %>% filter(KO_name %in% c("WT", "VPS15"))
# ggplot(toPlot, aes(y = overal_CC, x = parameter)) +
#   geom_jitter(alpha = 0.25) +
#   ylim(0,5) +
#   geom_point(data = toPlot.examples, aes(y = overal_CC, x = parameter, colour = KO_name)) +
#   geom_line(data = toPlot.examples , aes(y = overal_CC, x = as.numeric(parameter), colour = KO_name))
#   #ggrepel::geom_text_repel(data = toPlot %>% filter(parameter == "HXK2"), aes(label = KO_name, x = parameter, y = overal_CC), segment.alpha = 0.25, segment.size = 0.25)


## -- pca FCC ----
CCflux_matrix <- dataset.mca %>% filter(var_type == "flux", VAR != "AK", !is.na(CC)) %>% 
  dplyr::select(KO, parameter, VAR, CC) %>%  unite(parameter_VAR, parameter, VAR) %>%
  spread(parameter_VAR, CC) %>% tidy_pca()


CCfluxes_pca <- prcomp(CCflux_matrix)
cl <- cluster::pam(CCflux_matrix, 4)

toPlot <- CCfluxes_pca
xlabel <- paste("PC1", round(toPlot$sdev[1]/sum(toPlot$sdev),2))
ylabel <- paste("PC2", round(toPlot$sdev[2]/sum(toPlot$sdev),2))


p.scores_CCfluxes <- toPlot$x %>% as_tibble() %>% mutate(gene_name = exp_metadata$gene[match(rownames(toPlot$x), exp_metadata$ORF)],
                                                         cluster = cl$clustering[match(rownames(toPlot$x), names(cl$clustering))]) %>% 
  ggplot(aes(x = PC1 , y = PC2)) +
  geom_point(aes(colour = factor(cluster))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggrepel::geom_text_repel(aes(label = gene_name, x = PC1, y = PC2), size = 3,  segment.alpha = 0.25, segment.size = 0.25) +
  theme(aspect.ratio = 1) +
  xlab(label = xlabel) +
  ylab(label = ylabel) +
  theme_bw(base_family = "Helvetica") +
  theme(legend.position = "none", panel.grid = element_blank())
  

toPlot <- CCfluxes_pca %>% .$rotation %>% reshape2::melt()
names(toPlot) = c("variable", "component", "value")
toPlot <- toPlot %>% filter(component %in% c("PC1", "PC2")) %>% 
  mutate(variable = fct_reorder(variable, value)) %>%  
  arrange(component, -abs(value)) %>% 
  group_by(component) %>% 
  dplyr::filter(row_number() <= 30)

p.loadings_FCC <-  ggplot(toPlot, aes(x = variable, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Steady-state FCC") + 
  theme_bw(base_family = "Helvetica") +
  theme(legend.position = "none", panel.grid = element_blank())

p.grid_FCC <- plot_grid(p.scores_CCfluxes, p.loadings_FCC, labels = c("A", "B"), align = "h" )
p.grid_conc$landscape = T
plots.list = lappend(plots.list, p.grid_conc)
plots.list = lappend(plots.list, p.scores_CCfluxes)

sel <- toPlot %>% ungroup %>% arrange(component, desc(abs(value))) %>% filter(row_number() == 1) %>% dplyr::select(variable)

p.first_loading <- dataset.mca %>% filter(var_type == "flux", VAR != "AK", !is.na(CC)) %>% 
  dplyr::select(KO, parameter, VAR, CC) %>% 
  unite(parameter_VAR, parameter, VAR) %>% 
  filter(parameter_VAR %in% sel$variable) %>%
  mutate(cluster = cl$clustering[match(KO, names(cl$clustering))]) %>%
  ggplot(aes(fill = factor(cluster), x = CC)) +
    geom_density(alpha = 0.5) +
    xlab("Control coefficent of HXK2 over HXK_GLK1") +
    theme_bw(base_family = "Helvetica") +
    theme(legend.position = "none", panel.grid = element_blank())
    

### -- example HXK2 GLK ---
load("./R/objects/proteins.matrix.sva.0.5.1.RData")
load("./R/objects/orf2name._clean_.RData")

my_means <- function(proteins.matrix) {
  
  proteins.long = reshape2::melt(proteins.matrix, id.vars="rownames")
  names(proteins.long) = c("EG.StrippedSequence", "R.Label", "signal")
  proteins.long$ORF = exp_metadata$ORF[match(proteins.long$R.Label, exp_metadata$sample_name)]
  proteins.long.mean = tbl_df(proteins.long) %>% group_by(EG.StrippedSequence, ORF) %>% summarize(mean = mean(signal))
  proteins.mean.df = reshape2::dcast(proteins.long.mean, formula=EG.StrippedSequence~ORF, value.var="mean")
  
  proteins.mean.matrix = as.matrix(proteins.mean.df[,-1])
  rownames(proteins.mean.matrix) = as.matrix(proteins.mean.df$EG.StrippedSequence)
  return(proteins.mean.matrix)  
}

protein.matrix.mean = my_means(exp(proteins.matrix.sva.0.5.1))

proteins <- protein.matrix.mean %>% as_tibble() %>% 
  mutate(ORF = rownames(protein.matrix.mean)) %>% 
  dplyr::select(ORF, everything()) %>%
  gather(KO, value, none:YPR111W) %>%
  mutate(KO_name = exp_metadata$gene[match(KO, exp_metadata$ORF)],
         gene_name = orf2name$gene_name[match(ORF, orf2name$ORF)]) %>% filter(KO != "none")

ratio_sel <- c("GLK1", "HXK2")
toPlot <- proteins %>% filter(gene_name %in% ratio_sel) %>% group_by(ORF) %>% 
  mutate(norm_value=(value-min(value))/(max(value)-min(value))) %>% 
  arrange(KO, gene_name) %>% 
  group_by(KO, KO_name) %>% 
    summarise(ratio = norm_value[1]/norm_value[2]) %>% 
    mutate(cluster = cl$clustering[match(KO, names(cl$clustering))])

p.cluster_ratio <- toPlot %>% 
  ggplot(aes(y = ratio, x = factor(cluster) )) +
    stat_boxplot(geom ='errorbar', width = 0.25) +
    geom_boxplot(aes(fill = factor(cluster)), width = 0.75) + 
    ylab(paste("Proteins expression ratio,  ", ratio_sel[1], "/" ,ratio_sel[2], sep = "" )) +
    xlab("Cluster") +
    ylim(0, 5) +
    theme_bw(base_family = "Helvetica") +
    theme(legend.position = "none", panel.grid = element_blank()) 

p.cor_CC_ratio <- inner_join(toPlot, p.first_loading$data, by = c("KO", "cluster")) %>% 
  ggplot(aes(x = ratio, y = CC)) + 
    geom_point() +
    xlab(paste("Proteins expression ratio,  ", ratio_sel[1], "/" ,ratio_sel[2], sep = "" )) +
    ylab("Control coefficent of HXK2 over HXK_GLK1") +
    theme_bw(base_family = "Helvetica") +
    theme(legend.position = "none", panel.grid = element_blank()) 

# example of heatmaps from MAPKK
thr = 1
selected <- c("STE20", "STE11", "STE7", "FUS3")
p.FCC_heatmap <- dataset.mca %>% filter(var_type == "flux",  KO_name %in% selected , !is.na(CC)) %>% 
  mutate(CC_plot = ifelse(CC > thr, thr, CC ),
         CC_plot = ifelse(CC_plot < -thr, -thr, CC_plot),
         KO_name = factor(KO_name, levels = selected),
         VAR = factor(VAR, levels = rev(reaction_order)),
         parameter = factor(parameter, levels = parameter_order)) %>% 
  ggplot(aes(y = VAR, x = parameter, fill = CC_plot)) + 
  geom_tile(color="white", size=0.1) +
  #scale_fill_viridis(breaks = seq(-thr, thr, 0.1), limits = c(-thr,thr)) +
  scale_fill_gradient2(breaks = seq(-thr, thr, 0.25), limits = c(-thr,thr)) +
  facet_wrap(~KO_name, ncol = 4) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  theme_tufte(base_family="Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.ticks=element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.25)) +
  theme(strip.text.x = element_text(hjust = 0)) +
  theme(axis.text=element_text(size=7)) +
  ylab("Fluxes, J") +
  xlab("Enzymes, E")

#p.FCC_heatmap <- p.FCC_heatmap + theme(legend.position = c(.95, .95), legend.justification = c("right", "top"),legend.box.just = "right")

thr = 1
selected_var = "ADH_ADH1"

toPlot <- dataset.mca %>% ungroup() %>% filter(var_type == "flux",  !is.na(CC)) %>% filter(VAR == selected_var) %>%
  mutate(CC_plot = ifelse(CC > thr, thr, CC ),
         CC_plot = ifelse(CC_plot < -thr, -thr, CC_plot),
         VAR = factor(VAR, levels = rev(reaction_order))) %>% 
  group_by(KO_name) %>% mutate(rank = dense_rank(desc(CC)))

rank_order <- (toPlot %>% filter(KO == "WT") %>% ungroup %>% arrange(rank) %>% dplyr::select(parameter))$parameter %>% as.character()
toPlot <- toPlot %>% mutate(parameter = factor(parameter, levels = rank_order))

p.FCC_rank_heatmap <- toPlot %>%
  ggplot(aes(y = KO_name, x = parameter)) +
    geom_tile(color="white", size=0.1, aes(fill = factor(rank))) +
    #geom_point(aes(size = CC_plot))
    #geom_text(aes(label = rank), size = 3)+
    ggtitle(paste("Control variation of flux through", selected_var)) +
    scale_fill_grey(start = 0.05, end = 1, name = "Control rank") +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    coord_equal() +
    theme_tufte(base_family="Helvetica") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.ticks = element_blank()) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.25)) +
    theme(strip.text.x = element_text(hjust = 0)) +
    theme(axis.text = element_text(size=7)) + 
    theme(plot.title = element_text(size=8)) +
    ylab("Kinase mutants") +
    xlab("Enzymes, E")
    
plots.list = lappend(plots.list, p.FCC_rank_heatmap)

selected_var = "ATP"
toPlot <- dataset.mca %>% filter(var_type == "conc",  !is.na(CC)) %>% filter(VAR == selected_var) %>%
  mutate(CC_plot = ifelse(CC > thr, thr, CC ),
         CC_plot = ifelse(CC_plot < -thr, -thr, CC_plot)) %>% 
  group_by(KO_name) %>% mutate(rank = dense_rank(desc(CC)))

rank_order <- (toPlot %>% filter(KO == "WT") %>% ungroup %>% arrange(rank) %>% select(parameter))$parameter %>% as.character()
toPlot <- toPlot %>% mutate(parameter = factor(parameter, levels = rank_order))

p.CCC_rank_heatmap <- toPlot %>%
  ggplot(aes(y = KO_name, x = parameter)) +
  geom_tile(color="white", size=0.1, aes(fill = factor(rank))) +
  #geom_point(aes(size = CC_plot))
  ggtitle(paste("Control variation of ", selected_var, "concentration")) +
  scale_fill_grey(start = 0.05, end = 1, name = "Control rank") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  theme_tufte(base_family="Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.ticks = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.25)) +
  theme(strip.text.x = element_text(hjust = 0)) +
  theme(axis.text = element_text(size=7)) + 
  theme(plot.title = element_text(size=8)) +
  ylab("Kinase mutants") +
  xlab("Enzymes, E")

plots.list = lappend(plots.list, p.CCC_rank_heatmap)


## ---- MCA stats  ----

stats_table <- data.frame(stats_name = character(),
                          value = character(), 
                          rel_value = character(),
                          comment=character(),
                          stringsAsFactors=FALSE)

#number of strains with considerable (thr) changes comparing to WT
thr = 1
total_kinases = nrow(dataset.ss %>% ungroup() %>% distinct(KO)) - 1
summary_stats <- dataset.ss %>% ungroup %>% filter(abs_value_change > thr) %>% dplyr::select(KO, type) %>% distinct() %>% group_by(type) %>% summarise(n = n()/total_kinases)

stats_tmp <- data.frame(stats = "flux_changed_kinase_fraction", 
                        value = round((summary_stats %>% filter(type == "flux"))$n, 2),
                        comment = paste("Fraction of kinase mutants with >", thr, "changes"))
stats_table <- rbind(stats_table, stats_tmp)


stats_tmp <- data.frame(stats = "conc_changed_kinase_fraction", 
                        value = round((summary_stats %>% filter(type == "conc"))$n, 2),
                        comment = paste("Fraction of kinase mutants with >", thr, "changes"))
stats_table <- rbind(stats_table, stats_tmp)


#differences in FCC comparing to WT strain
"We observed that for XX % of enzymes the overall control over the fluxes (Methods) in XX% of kinase mutants 
for at is more than 2-fold different in absolute terms from the wild-type strain. "




## ----- figure version 1 --------

plot_figure_v1 <- function () {
  grid.newpage() 
  
  pushViewport(viewport(layout = grid.layout(161, 60)))
  grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=14, col="black"))
  print(p.ss_changes_density, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 1:60)) #fluxes and concentration at steady-state
  print(p.example_zoom, vp = viewport(layout.pos.row = 31:60, layout.pos.col = 1:60)) # example 
  print(p.example, vp = viewport(layout.pos.row = 33:52, layout.pos.col = 31:58)) #inset in example
  print(p.scores_CCfluxes + theme(aspect.ratio = 1), vp = viewport(layout.pos.row = 61:120, layout.pos.col = 1:40))
  print(p.first_loading + theme(aspect.ratio = 5/8), vp = viewport(layout.pos.row = 61:80, layout.pos.col = 41:60))
  print(p.cluster_ratio + theme(aspect.ratio = 5/8), vp = viewport(layout.pos.row = 61:80, layout.pos.col = 21:40))
  print(p.cor_CC_ratio + theme(aspect.ratio = 1), vp = viewport(layout.pos.row = 81:100, layout.pos.col = 41:60))
  print(p.FCC_heatmap, vp = viewport(layout.pos.row = 120:160, layout.pos.col = 1:60)) 
}

file_name = "MCA_figure_v1.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
plot_figure_v1() 
dev.off()

file_name = "MCA_figure_v1.png"
file_path = paste(figures_dir, file_name, sep="/")
png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
plot_figure_v1() 
dev.off()




file_name = paste("results", fun_name, sep = ".")
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


