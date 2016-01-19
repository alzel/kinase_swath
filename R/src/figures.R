#!/usr/bin/env Rscript
# Paper figures

rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "figures"

# Metabolite models results figure
load("R/objects/model_results.RData")
load("./R/objects/metabolite2iMM904._load_.RData")


detach("package:cowplot", unload=T)

## -- screen results PPP ----
model_results$dataset = factor(model_results$dataset)
screen_results = droplevels(model_results %>% filter(dataset == "screenPPP_AA"))
selected_PPP = c("PEP", "Glu", "X1.6.FP.", "X2...3.PG", "Pyr", "X6PGlu", "DHAP", "F6P", "R5P", "G3P", "G6P...F6P", "X5P...Rl5P", "E4P", "S7P")

screen_results.selected = droplevels(screen_results[screen_results$metabolite %in% selected_PPP,])

toPlot = screen_results.selected
p = ggplot(toPlot, aes(x = metabolite, y = r.squared)) + 
       stat_boxplot(geom ='errorbar') + 
       geom_boxplot(width=0.66, fill="black", color = "white") + 
       scale_x_discrete(labels = metabolite2iMM904$official_name[match(levels(toPlot$metabolite), metabolite2iMM904$id)]) +
       coord_flip()

plots.list = lappend(plots.list, p)

## -- seleceted samples TCA ----

model_results$dataset = factor(model_results$dataset)
selectedTCA = droplevels(model_results %>% filter(dataset == "selectedTCA"))

selectedTCA$metabolite = factor(selectedTCA$metabolite, 
                         levels = unique((selectedTCA %>% arrange(median.r.squared))$metabolite))
toPlot = selectedTCA

p = ggplot(toPlot, aes(x = metabolite, y = r.squared)) + 
          stat_boxplot(geom ='errorbar') + 
          geom_boxplot(width=0.66, fill="black", color = "white") + 
          scale_x_discrete(labels = metabolite2iMM904$official_name[match(levels(toPlot$metabolite), metabolite2iMM904$id)]) +
          coord_flip()

plots.list = lappend(plots.list, p)

## -- seleceted samples AA ----

model_results$dataset = factor(model_results$dataset)
selectedAA = tbl_df(droplevels(model_results %>% filter(dataset == "selectedAA")))

selectedAA$metabolite = factor(selectedAA$metabolite, 
                                levels = unique((selectedAA %>% arrange(median.r.squared))$metabolite))
toPlot = selectedAA
p = ggplot(toPlot, aes(x = metabolite, y = r.squared)) + 
            stat_boxplot(geom ='errorbar') + 
            geom_boxplot(width=0.66, fill="black", color = "white") + 
            scale_x_discrete(labels = metabolite2iMM904$official_name[match(levels(toPlot$metabolite), metabolite2iMM904$id)]) +
            coord_flip()
plots.list = lappend(plots.list, p)
load("./R/objects/KEGG.pathways.stats.f.analysis1.RData")


KEGG.pathways.summary = KEGG.pathways.stats.f %>% group_by(A) %>% summarise(avgEC = mean(EC.coverage))
KEGG.pathways.summary = KEGG.pathways.stats.f %>% group_by(B) %>% summarise(avgEC = avgEC[1])
KEGG.pathways.summary = KEGG.pathways.summary %>% filter(B != "Biosynthesis of other secondary metabolites")
KEGG.pathways.summary$B = factor(KEGG.pathways.summary$B, levels=KEGG.pathways.summary$B[order(KEGG.pathways.summary$avgEC)])

toPlot = KEGG.pathways.summary
library(scales)
#slide for embo
p = ggplot(toPlot, aes(x=B, y=avgEC)) +
  geom_bar(stat="identity", width=0.66) +
  scale_y_continuous(breaks = seq(0,0.9,0.1),labels = percent_format()) +
  ylab("Average coverage of reactions by measured proteins") +
  xlab("KEGG pathways category") +
  coord_flip()

file_name = paste(fun_name,"pathway_coverage", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width = 8.27)


plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")

# ggplot(KEGG.pathways.stats.f, aes(x = B, y = EC.coverage)) +
#        geom_boxplot() +
#        geom_point(data=KEGG.pathways.stats.f, aes(x=jitter(as.numeric(B)))) + coord_flip()

# 
