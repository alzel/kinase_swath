rm(list = ls())
library(tidyverse)
library(forcats)
library(cowplot)


fun_name = "MCA"
figures_dir = "./figures"

lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

plots.list = list()


load("./R/objects/exp_metadata._clean_.RData")

## -- Steady-concentration fluxes ----
dataset.ss <- read_delim("./results/2017-02-22/ss_fluxes_concentrations_in_kinase_KOs.csv", delim = "\t")
dataset.ss <- dataset.ss %>% group_by(variable, type) %>% mutate(z_value = (value - mean(value, na.rm=T))/sd(value, na.rm=T))

p.ss <- ggplot(data = dataset.ss, aes(x = variable, y = z_value)) +
  geom_boxplot() +
  facet_wrap(~type, scales = "free") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

plots.list = lappend(plots.list, p.ss)

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
  filter(type == "flux", is.nan(z_value) == F ) %>% 
  dplyr::select(KO, variable, z_value) %>% 
  spread(variable, z_value) %>% tidy_pca 

fluxes_pca <- flux_matrix  %>% prcomp()

conc_matrix <- dataset.ss %>% ungroup() %>%
  filter(type == "conc", is.nan(z_value) == F ) %>% 
  dplyr::select(KO, variable, z_value) %>% 
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

dataset.mca <- bind_rows(flux_mca, conc_mca) %>% dplyr::rename(VAR = variable)


dataset.mca %>% group_by(var_type, VAR, KO) %>% 
  summarise(min = min(CC, na.rm = T),
            min_plot = ifelse(min < -3, -3, min),
            min_parameter = parameter[which.min(CC)][1],
            max = max(CC, na.rm = T),
            max_plot = ifelse(max > 3, 3, max),
            max_parameter = parameter[which.max(CC)][1],
            max_abs = max(abs(CC), na.rm = T),
            max_abs_parameter = parameter[which.max(abs(CC))][1]) %>% 
  dplyr::select(var_type, VAR, KO, max_plot, min_plot ) %>% reshape2::melt(id.vars = c("var_type", "VAR", "KO")) %>% 
ggplot(aes(x = VAR, y = value)) +
    ylim(c(-3,3)) +
    geom_violin() +
    geom_jitter(aes(colour = variable), alpha = 0.3) +
    facet_wrap(~var_type, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

toPlot <- dataset.mca %>% 
  filter(VAR == "HXT", var_type == "flux") 
p.example <-  ggplot(toPlot,  aes(x = parameter, y = CC)) +
    ggtitle("Control coefficients of flux through HXT") +
    geom_boxplot() +
    #geom_jitter() +
    geom_point(data = toPlot %>% filter(KO == "WT"), aes(colour = "red")) +
    theme_bw() + theme(aspect.ratio = 5/8) 
plots.list = lappend(plots.list, p.example)

p.example_zoom <- ggplot(toPlot,  aes(x = parameter, y = CC)) +
  ggtitle("Control coefficients of flux through HXT") +
  geom_jitter(alpha = 0.25) +
  geom_point(data = toPlot %>% filter(KO == "WT"), aes(colour = "red")) +
  ylim(c(0, 0.5)) +
  theme_bw() + 
  theme(legend.position="none", aspect.ratio = 5/8)
plots.list = lappend(plots.list, p.example_zoom)

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


