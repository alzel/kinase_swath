rm(list = ls())
library(caret)
library(relaimpo)
source("./R/boot.R")
source("./R/functions.R")
library("cowplot")
library("scales")

### summarizes results from linear regression models
### Analysis of BRENDA enzymes
### produces examples of glutamate
### Makes figures of energy metabolites


load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")
input_path = "./results/2015-09-29/linear_models"
load("./R/objects/proteins.matrix.combat.quant.RData")

load("./R/objects/gene.annotations._load_.RData")
orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]


fun_name = "lm_models2"

#filesToProcess = dir(path=input_path, pattern = "[12].[01]+.linear_models.RData$", recursive=F)
filesToProcess = dir(path=input_path, pattern = "[123].[01]+.linear_models.RData$", recursive=F)
filesToProcess = grep(pattern="imputed", invert=T, filesToProcess, value=T)

pattern.p = "data.(\\w+).(.*?).([0-9]+).([0-9]+).linear_models.RData$"

matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

# read_models.caret = function(x) {
#   z <<- x
#   file_name = paste(input_path,x[[1]], sep="/") 
#   my_models = get(load(file_name))
#   
#   if (length(my_models$caret_models) > 1 ) {
#     t = resamples(my_models$caret_models)
#     t.summary = summary(t)
#     fit = my_models$caret_models[[which.min(t.summary$statistics$RMSE[,"Mean"])]]$finalModel
#   } else {
#     fit  = my_models$caret_models[[1]]$finalModel
#   }
#   fit.s = summary(fit)
#   
#   if (length(fit$coefficients) > 2 ) {
#     table = data.frame(feature = names(calc.relimp(fit, type="lmg")@lmg),
#                        importance = calc.relimp(fit, type="lmg")@lmg,
#                        row.names=NULL)
#     table$betas = fit$coefficients[match(table$feature,  names(fit$coefficients))]
#     
#   } else if (length(fit$coefficients) == 2) {
#     
#     table = data.frame(feature = names(fit$coefficients)[2],
#                        importance = fit.s$r.squared, row.names=NULL,
#                        betas = fit$coefficients[2])
#   } else {
#     table = data.frame(feature = NA,
#                        importance = NA,
#                        betas = NA)
#   } 
#     
#   table$r.squared = fit.s$r.squared
#   table$adj.r.squared = fit.s$adj.r.squared
#   table$dataset = factor(x[[2]])
#   table$metabolite  = factor(x[[3]])
#   table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
#   table$degree  = factor(x[[4]])
#   table$ismetIncluded  = factor(x[[5]])
#   table$file =  factor(x[[1]])
#   return(table)
# }

read_models.aic = function(x) {
  z <<- x
  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))
  
  table = my_models$summaries
  
  table$dataset = factor(x[[2]])
  table$species  = factor(x[[3]])
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = factor(x[[4]])
  table$ismetIncluded  = factor(x[[5]])
  table$file =  factor(x[[1]])
  return(table)
}

# file.list = lapply(matches, FUN=read_models.caret)
# all_final_models.caret = do.call(rbind.data.frame, file.list)

file.list = lapply(matches, FUN=read_models.aic)
all_final_models.aic = do.call(rbind.data.frame, file.list)
all_final_models.aic$metabolite = all_final_models.aic$species
all_final_models.aic$metabolite = sub(x=all_final_models.aic$metabolite, pattern="log.quant.(.*)", replacement="\\1")
all_final_models.aic$metabolite = sub(x=all_final_models.aic$metabolite, pattern="log.(.*)", replacement="\\1")

all_final_models.aic$normalization = "bc"
all_final_models.aic$normalization[grep(pattern="log", x=all_final_models.aic$species)] = "log"
all_final_models.aic$normalization[grep(pattern="log.quant", x=all_final_models.aic$species)] = "log.quant"
all_final_models.aic$normalization = factor(all_final_models.aic$normalization)

auto_thr = 0.05
auto_thr.bonferonni = auto_thr/nrow(all_final_models.aic)
all_final_models.aic$isAutocorrelation.bonferoni = ifelse(all_final_models.aic$bg.p.value < auto_thr.bonferonni, 1, 0)
all_final_models.aic$isAutocorrelation = ifelse(all_final_models.aic$bg.p.value < auto_thr, 1, 0)

#selecting best representative model based on adj R2
models.filtered = all_final_models.aic %>% filter(type == "after", 
                                                  isImputed == 0 , isAutocorrelation.bonferoni == 0) %>%
                                                  group_by(dataset, metabolite, ismetIncluded, degree) %>% 
                                                  summarise(formula = formula[which.max(adj.r.squared)],
                                                            median.cv.r2 = median.cv.r2[which.max(adj.r.squared)],
                                                            adj.r.squared = adj.r.squared[which.max(adj.r.squared)],
                                                            r.squared = r.squared[which.max(adj.r.squared)],
                                                            model = model[which.max(adj.r.squared)],
                                                            type = type[which.max(adj.r.squared)],
                                                            datafile = datafile[which.max(adj.r.squared)])


plots.list = list()

## all adjusted R2 plots
toPlot = subset(all_final_models.aic, ismetIncluded == 0 & degree == 1 )
toPlot = all_final_models.aic
p = ggplot(toPlot, aes(x = metabolite, y = adj.r.squared, colour=type, shape = normalization)) +
           geom_point() +
           geom_point(data=toPlot, aes(y = median.cv.r2 , colour = "black")) +
           facet_wrap(ismetIncluded~dataset+degree, scales = "free") + 
           theme(axis.text.x = element_text(angle = 90, hjust = 1))

plots.list = lappend(plots.list, p)


## all median.cv.r2 plots
toPlot = subset(all_final_models.aic, type == "before")
p = ggplot(toPlot, aes(x = metabolite, y = median.cv.r2,  shape = normalization)) +
  geom_point() +
  #geom_point(data=toPlot, aes(y = median.cv.r2 , colour = 1)) +
  facet_wrap(ismetIncluded~dataset+degree, scales = "free") + 
  #facet_wrap(~dataset, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plots.list = lappend(plots.list, p)

toPlot = models.filtered %>% ungroup() %>% arrange(dataset, adj.r.squared)

toPlot$metabolite = factor(toPlot$metabolite, levels = unique(toPlot$metabolite))
p = ggplot(toPlot, aes(x = metabolite, y = adj.r.squared)) +
          geom_point() +
          geom_linerange(data = toPlot , aes(x=metabolite,ymin=0, ymax=adj.r.squared))+
          facet_wrap(ismetIncluded~dataset+degree, scales = "free") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
plots.list = lappend(plots.list, p)

# 
# all_final_models.aic.summary = all_final_models.aic %>% 
#   group_by(formula, dataset, metabolite, isImputed, degree, ismetIncluded) %>% 
#   summarize(r.squared = mean(r.squared),
#             adj.r.squared = mean(adj.r.squared),
#             aic = median(aic))
# 
# formulas.coef = all_final_models.aic %>% dplyr::select(formula, coeficients) %>% distinct()
# all_final_models.aic.merged = merge(all_final_models.aic.summary, formulas.coef, by="formula")
# all_final_models.aic.merged = droplevels(all_final_models.aic.merged[!all_final_models.aic.merged$coeficients == "(Intercept)",])
# 
# 
# 
# all_final_models.caret$isMetabolite = ifelse(all_final_models.caret$feature %in% metabolite2iMM904$id,1,0)
# all_final_models.caret = droplevels(filter(all_final_models.caret, dataset != "PPP_AA"))

      
## getting coverages
yeast.model = iMM904
yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters
edgelist = unique(droplevels(na.omit(subset(yeast.model, metabolite != "h"  & metabolite !="h2o" & metabolite != "nh4" & metabolite != "ppi" , select = c("metabolite", "gene")))))
B <- graph.data.frame(edgelist)
V(B)$type <- V(B)$name %in% edgelist$metabolite
stopifnot(is.bipartite(B))


measured.proteins = row.names(proteins.matrix.combat.quant)

coverage = ddply(all_final_models.aic, .(dataset, metabolite), 
      .fun=function(x) {
  
        dataset = as.character(unique(x$dataset))
        i = as.character(unique(x$metabolite))
        
        current_nodes = as.character(metabolite2iMM904$model_name[metabolite2iMM904$id == i])    
        measured.metabolites = as.character(unique(all_final_models.aic$metabolite[all_final_models.aic$dataset == dataset ]))
        measured.metabolites.model = as.character(unique(metabolite2iMM904$model_name[metabolite2iMM904$id %in% measured.metabolites]))
  
        if (length(current_nodes) == 0) {
          message(paste("No ", i, "found in network"))
          next()
        }
        
        SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=2, nodes=current_nodes))))    
        genes = unique(V(SUB)$name[V(SUB)$type == 0])
        metabolites = unique(V(SUB)$name[V(SUB)$type == 1])
        
        coverage.met = data.frame(metabolite = i, 
                                  dataset = dataset, 
                                  type = "metabolic",
                                  neighbour = metabolites,
                                  isMeasured = ifelse(metabolites %in% measured.metabolites.model, 1, 0))
        
        coverage.gene = data.frame(metabolite = i, 
                                  dataset = dataset, 
                                  type = "gene",
                                  neighbour = genes,
                                  isMeasured = ifelse(genes %in% measured.proteins, 1, 0))
        
        return(rbind(coverage.met, coverage.gene))
        
      })  
  

coverage.stats = coverage %>% group_by(metabolite, dataset, type) %>% summarise(coverage = sum(isMeasured)/length(isMeasured))
coverage.stats.wide = dcast(formula=metabolite+dataset~type, value.var="coverage", data=coverage.stats)

#R2 is independent of coverage 
toPlot = merge(models.filtered, coverage.stats.wide)

res = cor(toPlot$r.squared, toPlot$gene)
res = data.frame(res)
p = ggplot(toPlot, aes(x=gene, y=r.squared)) +
           geom_text(data = res, aes(x=0.25, y=0.5, label=paste("rho =", round(res,3)))) +
           geom_point()

plots.list = lappend(plots.list, p)

## -- model parameters ----
read_models.models = function(x) {
  z <<- x

  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))
  
  models = my_models$models
  
  tmp = data.frame()
  for (i in 1:length(models)) {
    coefficients1 = models[[i]]$before$coefficients[-1]
    tmp = rbind(tmp, data.frame(model = i, type = "before", coefficients = coefficients1, variables =  names(coefficients1)))
    tmp = rbind(tmp,data.frame(model = i, type = "before", coefficients = NA, variables = "stats"))
    
    coefficients2 = models[[i]]$after$coefficients[-1]
    tmp = rbind(tmp, data.frame(model = i, type = "after", coefficients = coefficients2, variables =  names(coefficients2)))
    tmp = rbind(tmp,data.frame(model = i, type = "after", coefficients = NA, variables = "stats"))
   
  }
  
  table = tmp
  
  table$dataset = factor(x[[2]])
  table$species  = factor(x[[3]])
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = factor(x[[4]])
  table$ismetIncluded  = factor(x[[5]])
  table$file =  factor(x[[1]])
  table = merge(table,  my_models$summaries, by = c("model", "type"))
  
  annotations = my_models$summaries %>% dplyr::select(model, type, adj.r.squared)
  colnames(annotations)[3] = "stats_text"

  annotations$variables  = "stats"
  table = merge(table, annotations, by = c("model", "type", "variables"), all=T)
  return(table)
}

file.list = lapply(matches, FUN=read_models.models)
all_final_models.models = do.call(rbind.data.frame, file.list)
all_final_models.models$metabolite = all_final_models.models$species
all_final_models.models$metabolite = sub(x=all_final_models.models$metabolite, pattern="log.quant.(.*)", replacement="\\1")
all_final_models.models$metabolite = sub(x=all_final_models.models$metabolite, pattern="log.(.*)", replacement="\\1")

all_final_models.models$normalization = "bc"
all_final_models.models$normalization[grep(pattern="log", x=all_final_models.models$species)] = "log"
all_final_models.models$normalization[grep(pattern="log.quant", x=all_final_models.models$species)] = "log.quant"
all_final_models.models$normalization = factor(all_final_models.models$normalization)

auto_thr = 0.05
auto_thr.bonferonni = auto_thr/nrow(all_final_models.models)
all_final_models.models$isAutocorrelation.bonferoni = ifelse(all_final_models.models$bg.p.value < auto_thr.bonferonni, 1, 0)
all_final_models.models$isAutocorrelation = ifelse(all_final_models.models$bg.p.value < auto_thr, 1, 0)


#dirty to get ec numbers for brenda extraction
# ec.numbers = gene.annotations %>% filter(gene.annotations$V3 == "EC number")
# ec.only = grep(pattern="-", unique(ec.numbers$V1), invert=T, value=T)
# write.table(x = ec.only, file="./results/2015-10-01/ec.numbers.list",row.names=F, col.names=F, quote=F) 


all_final_models.models$mode = factor(ifelse(all_final_models.models$coefficients > 0,1,0))

#selecting best representative model based on adj R2 out of all normalization methods
#all_final_models.models = all_final_models.models %>% filter(type == "after", ismetIncluded == 0, degree == 1,
all_final_models.models = all_final_models.models %>% filter(type == "after",  
                                   isImputed == 0 , isAutocorrelation.bonferoni == 0) %>%
                                   group_by(dataset, model, metabolite, ismetIncluded, degree) %>% 
                                   mutate(the_best = adj.r.squared == max(adj.r.squared)) %>% # best among normalization methods
                                   group_by(dataset, metabolite, ismetIncluded, degree) %>% 
                                   mutate(the_super_best = adj.r.squared == max(adj.r.squared)) #the best model

all_final_models.models$varname = orf2name$gene_name[match(all_final_models.models$variables, orf2name$ORF)]
all_final_models.models$varname[which(is.na(all_final_models.models$varname))] = as.character(all_final_models.models$variables[is.na(all_final_models.models$varname)])
tmp.lev = unique(sort(all_final_models.models$varname))
all_final_models.models$varname = factor(all_final_models.models$varname, levels=c("stats", tmp.lev[tmp.lev != "stats"]))



## ---- variable plots ----
for ( i in c("AA", "PPP_AA", "TCA", "TCA_AA") ) {
    
    toPlot  = dplyr::filter(all_final_models.models, type == "after", adj.r.squared > 0.1, isImputed == 0 , isAutocorrelation.bonferoni == 0, dataset == i, the_best == T)
        
    for( d in unique(as.character(toPlot$degree))) {
      toPlot.tmp = toPlot %>% ungroup() %>% filter(degree == d)
      variables = unique(toPlot.tmp$metabolite)
                    
      noVars = length(variables)
      noPlots = 12
      
      plotSequence <- c(seq(0, length(variables)-1, by = noPlots), noVars)
          
      for(ii in 2:length(plotSequence)){
        # select start and end of variables to plot
        start <- plotSequence[ii-1] + 1
        end <- plotSequence[ii]  
        
        tmp = droplevels(toPlot.tmp[toPlot.tmp$metabolite %in% variables[start:end],])
        
        p = ggplot() +
          geom_text(data=tmp, aes(x=factor(model), y = varname, label=round(stats_text,2))) +
          geom_point(data=tmp, aes(x=factor(model), y = varname, 
                                   size=abs(coefficients), color = mode)) +
          facet_wrap(dataset~metabolite+degree, scales="free") +
          xlab("Candidate model") + 
          ylab("Enzyme predictors") +
          theme_light() + 
          scale_size_continuous(name="Effect size",
                                breaks=c(0.25, 0.5, 1),
                                labels = c("low", "medium", "strong")) +
          scale_color_discrete(name="Predictor's effect",
                               breaks = c(0, 1),
                               labels = c("negative", "positive") )
        
        plots.list = lappend(plots.list, p)      

    }
  }
}



# BRENDA vs concentrations

brenda <- read.delim("./data/2015-10-07/brenda.txt")
load("./R/objects/proteins.matrix.combat.RData")

ec.gene = unique(gene.annotations[gene.annotations$V3 == "EC number",c(1,4)])

#coverage of measured EC 
coverage.ec = droplevels(coverage[coverage$type == "gene",])
coverage.ec = merge(coverage.ec, ec.gene, by.x = "neighbour", by.y = "V4")
coverage.ec$kegg_id = metabolite2iMM904$kegg_id[match(coverage.ec$metabolite, metabolite2iMM904$id)]

brenda.f = brenda[!brenda$KEGGID == "",]
brenda.f = brenda.f[grep(pattern="mutant|recombinant", x= brenda.f$commentary, invert=T),]




#######
load("./R/objects/dataTCA.create_datasets.RData")
load("./R/objects/dataAA.create_datasets.RData")

metabolitesTCA.long = melt(dataTCA$metabolites)
metabolitesTCA.long$dataset = "TCA"


metabolitesAA.long = melt(dataAA$metabolites)
metabolitesAA.long$dataset = "AA"


# adjust concentration with volume and OD from this paper: doi:10.1016/j.funbio.2009.11.002
my.vol = c(median = 45.54, sd = 0.9) * 1e-15 # cell vol
my.cells = 3.2 * 10^7 * 1.5*5 # median of spectrophotometre data
ex.vol = 100*1e-6

metabolitesTCA.long$concentration = metabolitesTCA.long$value*ex.vol/(my.cells*my.vol[1])/1000 # mM

#different dilution used fot AA protocol
ex.vol = 200*1e-6

metabolitesAA.long$concentration = metabolitesAA.long$value*ex.vol/(my.cells*my.vol[1])/1000 # mM


metabolites.long = rbind(metabolitesAA.long, metabolitesTCA.long)
metabolites.long = merge(metabolites.long, unique(droplevels(subset(metabolite2iMM904,select=c("id", "kegg_id")))), by.x="X2", by.y="id")
metabolites.long = merge(brenda.f, metabolites.long, by.x="KEGGID", by.y="kegg_id")
metabolites.long = metabolites.long[metabolites.long$kmValue > 0,]
metabolites.long$label = metabolite2iMM904$official_name[match(metabolites.long$KEGGID, metabolite2iMM904$kegg_id )]

models.summary = all_final_models.models %>% filter(degree==1, ismetIncluded == 0,  the_super_best == T) %>% 
  dplyr::select(metabolite, degree, dataset,variables, adj.r.squared) 

models.summary = models.summary[models.summary$variables != "stats",]
models.summary = merge(models.summary, ec.gene, by.x = "variables", by.y = "V4")
models.summary$kegg_id = metabolite2iMM904$kegg_id[match(models.summary$metabolite, metabolite2iMM904$id)]
models.summary = models.summary%>% arrange(metabolite, dataset)



ec.presence  = unique(ddply(metabolites.long, .variables=.(KEGGID, dataset), 
                            .fun=function(x) {
                              met_id = as.character(unique(x$KEGGID))
                              z = data.frame (ecNumber = x$ecNumber,
                                              inNetwork = x$ecNumber %in% coverage.ec$V1[coverage.ec$kegg_id == met_id],
                                              isPredictor = x$ecNumber %in% models.summary$V1[models.summary$kegg_id == met_id] )
                              return(z)
                            }))

metabolites.long = merge(metabolites.long, ec.presence, by = c("KEGGID", "dataset", "ecNumber"))
toPlot = metabolites.long %>% filter(inNetwork == T)

p = ggplot(toPlot, aes(x = concentration, y = kmValue)) +
  geom_point() +
  facet_wrap(~label, scales="free")
plots.list = lappend(plots.list, p)

# checking for saturation
metabolites.long = metabolites.long %>% mutate(ratio = concentration/kmValue)
toPlot = metabolites.long %>% filter(inNetwork == T)
points = metabolites.long %>% filter(isPredictor == T)


p = ggplot(toPlot, aes(y=log(ratio), x = label)) +  
  geom_boxplot() +
  geom_point(data=points, aes(y=log(ratio), x = label), col="red") +
  geom_hline(yintercept = 0) +
  facet_wrap(~dataset, scales = "free", ncol=1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plots.list = lappend(plots.list, p)

points = metabolites.long %>% filter(isPredictor == T)
points.AA = points %>% filter(dataset == "AA")
points.TCA = points %>% filter(dataset == "TCA")
sum(log(metabolites.long$ratio) > 0, na.rm=T) / nrow(metabolites.long[!is.na(metabolites.long$ratio),]) # total
sum(log(points.AA$ratio) > 0, na.rm=T) / nrow(points.AA[!is.na(points.AA$ratio),]) 


saturation.stats = metabolites.long %>% group_by(dataset, isPredictor) %>% summarize (p.value = t.test(log(ratio), mu=0, alternative="greater")$p.value)
saturation.stats$x = 5
saturation.stats$y = 300
saturation.stats$y[saturation.stats$isPredictor == T] = 250

toPlot = metabolites.long %>% filter(inNetwork == T)
p = ggplot(toPlot, aes(x=log(ratio))) +  
      geom_histogram()+
      geom_histogram(data=points, aes(x=log(ratio), fill=isPredictor)) +
      geom_text(data=saturation.stats, aes(x=x, y=y, color = isPredictor, label = paste("P-value =", round(p.value,2)))) +
      geom_vline(xintercept = 0) +
      facet_wrap(~dataset, scales = "free_y", ncol=1)

plots.list = lappend(plots.list, p)


### ratio whether it is predictor of not
a = (metabolites.long %>% filter(inNetwork == T, isPredictor == F) %>% dplyr::select(ratio))$ratio
b = (metabolites.long %>% filter(inNetwork == T, isPredictor == T) %>% dplyr::select(ratio))$ratio

toPlot = metabolites.long %>% filter(inNetwork == T)
stats = data.frame(label_text = c(median(b,na.rm=T)/median(a, na.rm=T),
                                wilcox.test(b,a)$p.value),
                   x = 1,
                   y = c(5,4))

p = ggplot(toPlot, aes(x = isPredictor, y = log(ratio))) + 
       geom_boxplot(width=0.2)+
       #geom_violin(alpha=0)+
       geom_text(data=stats, aes(x=x, y=y, label = label_text)) +
       panel_border() +
       theme(aspect.ratio = 8/3)
       #geom_point(data=toPlot, aes(x=jitter(as.numeric(isPredictor)) + 1 , y = log(ratio)), alpha = 0.1)
plots.list = lappend(plots.list, p)


# predictors AA vs TCA
toPlot = metabolites.long %>% filter(inNetwork == T, isPredictor == T)
#toPlot = metabolites.long
stats = data.frame(label_text = c(log(median(toPlot$ratio[toPlot$dataset == "AA"], na.rm=T)/median(toPlot$ratio[toPlot$dataset == "TCA"], na.rm=T)),
                             wilcox.test(log(toPlot$ratio[toPlot$dataset == "AA"]), log(toPlot$ratio[toPlot$dataset == "TCA"]))$p.value,
                             sum(log(toPlot$ratio) > 0, na.rm=T)/length(toPlot$ratio)), #above Km
                      
                   x = -7,
                   y = c(0.15, 0.13, 0.1))

toPlot$dataset = factor(toPlot$dataset)
  
p = ggplot() +  
      geom_density(data=toPlot, aes(x=log(ratio), fill=dataset), alpha = 0.5) +
      geom_text(data=stats, aes(x=x+2, y=y, label = label_text)) +
      panel_border() + 
      theme(aspect.ratio = 0.625)

plots.list = lappend(plots.list, p)


## -- glutamate/glutamine examples ----
measured.proteins = row.names(proteins.matrix.combat.quant)
yeast.model = iMM904
yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters
yeast.model = yeast.model[yeast.model$gene %in% measured.proteins, ]

yeast.model = yeast.model %>% group_by(metabolite, gene) %>% dplyr::mutate( from = ifelse(side == "substrate", as.character(metabolite), as.character(gene)),
                                                                            to   = ifelse(side == "substrate", as.character(gene), as.character(metabolite)),
                                                                            from.ec = ifelse(side == "substrate", as.character(metabolite), as.character(ec_number)),
                                                                            to.ec   = ifelse(side == "substrate", as.character(ec_number), as.character(metabolite)))

met = c("glutamate") 


current_nodes = as.character(metabolite2iMM904$model_name[metabolite2iMM904$id %in% met])
edgelist = data.frame(yeast.model %>% ungroup() %>% filter(metabolite %in% current_nodes))


tmp.edge = data.frame(edgelist)

# for (i in 1:nrow(edgelist)) {
#   if ( edgelist[i, "directionality"] == "<->" ) {
#     tmp.edge = rbind(tmp.edge, data.frame(metabolite = edgelist[i,"metabolite"],
#                                           gene = edgelist[i, "gene"],
#                                           from = edgelist[i,"to"], to = edgelist[i,"from"], directionality = "<->" ))
#   }
# }

model.edges = tmp.edge

selected.model = all_final_models.models %>% ungroup() %>% filter(metabolite %in% met, the_super_best == T)
selected.model$metabolite_name = as.character(metabolite2iMM904$model_name[match(selected.model$metabolite, metabolite2iMM904$id)])


model.edges = merge(model.edges, subset(selected.model, variables!= "stats", select = c("metabolite_name", "variables", "coefficients")), all.x=T,
                          by.x = c("metabolite","gene"),
                          by.y = c("metabolite_name", "variables"))
    
names(model.edges)[length(names(model.edges))] = "weights"
# 
# idx = match(model.edges$gene, selected.model$variables)
# model.edges$weights = 0
# model.edges$weights = selected.model$coefficients[idx]

#model.edges$weights[is.na(model.edges$weights)] = 0

model.edges = model.edges %>% arrange(gene, ec_number) %>% 
                               group_by(metabolite, gene) %>% 
                               mutate(isCorrect = ifelse(side == "substrate" && weights < 0 && directionality == "->", T, 
                                                        ifelse(side == "product" && weights > 0 && directionality == "->", T, 
                                                               ifelse(directionality == "<->", T, F))))

model.edges = model.edges %>% arrange(gene, ec_number) 
model.edges$gene_name = orf2name$gene_name[match(model.edges$gene, orf2name$ORF)]


graph.edges = unique(model.edges[,c("from", "to", "reaction", "ec_number", "gene", "directionality", "weights", "isCorrect", "gene_name")])
graph.edges$abs_weight= abs(graph.edges$weights)
graph.edges$abs_weight[graph.edges$abs_weight >= 1] = 1
graph.edges$effect = ifelse(graph.edges$weights > 0, 1, 0)
graph.edges$effect.label = cut(abs(graph.edges$abs_weight), breaks=3, labels=c("low", "medium", "strong"))

B <- graph.data.frame(graph.edges, directed=T)
#B = as.directed(B, mode = "arbitrary")

V(B)$type <- V(B)$name %in% unique(yeast.model$metabolite)


name_idx = na.omit(match(V(B)$name, orf2name$ORF))
B_idx = which(!is.na(match(V(B)$name, orf2name$ORF)))
V(B)$label = 1:length(V(B)$name)
V(B)$label[B_idx] = as.character(orf2name$gene_name[name_idx])

name_idx = na.omit(match(V(B)$name, metabolite2iMM904$model_name))
B_idx = which(!is.na(match(V(B)$name, metabolite2iMM904$model_name)))
V(B)$label[B_idx] = as.character(metabolite2iMM904$official_name[name_idx])


V(B)$model = ifelse(V(B)$name %in% as.character(selected.model$variables), 1, 0)
plot.igraph(B, edge.arrow.size=0.5, edge.arrow.width=0.5)

file_name = paste(met, fun_name, "graphml", sep=".")
file_path = paste(output_dir, file_name, sep="/")

write.graph(B, file=file_path, format="graphml")

## -- Energy metabolites example ----

example.models = all_final_models.models %>% ungroup() %>% filter(dataset == "PPP_AA", degree == 1, ismetIncluded == 0, the_super_best == T, metabolite %in% c("ATP", "ADP", "AMP")) %>% 
                                                           group_by(metabolite) %>% summarize(model = model[1],
                                                                                              type = type[1],
                                                                                              file = file[1],
                                                                                              adj.r.squared = adj.r.squared[1],
                                                                                              median.cv.r2 =  median.cv.r2[1])


tmp.list = list()
for(i in 1:nrow(example.models)) {
  tmp.list[[i]] = matrix(t(example.models[i,]), nrow=1)
}


read_models = function(x) {
  file_name = paste(input_path,x[[4]], sep="/") 
  my_models = get(load(file_name))

  models = my_models$models
  fit = models[[as.numeric(x[[2]])]][[x[[3]]]]
  yhat = predict(fit, fit$model[,-1])
  table = data.frame(metabolite = x[[1]],
                     model = x[[2]],
                     type = x[[3]],
                     file = x[[4]],
                     y = fit$model[,1],
                     yhat = yhat,
                     adj.r.squared = x[[5]],
                     median.cv.r2 =  x[[6]])
  return(table)
}



prediction.list = lapply(tmp.list, FUN=read_models)
prediction.models = do.call(rbind.data.frame, prediction.list)
toPlot = prediction.models

stats.text = prediction.models %>% group_by(metabolite, model) %>% summarise(adj.r.squared = as.numeric(as.character(adj.r.squared[1])),
                                                                             median.cv.r2 = as.numeric(as.character(median.cv.r2[1])))
stats.text$x = -2
stats.text$y = seq(2,1,length.out=nrow(stats.text))

p = ggplot(toPlot, aes(x = yhat, y=y, color = factor(metabolite))) +
      geom_point(size = 4) +
      geom_text(data=stats.text, aes(x=x,y=y,color=metabolite, label=round(adj.r.squared,2)))+
      geom_abline(slope = 1) +
      ylim(c(-3.5,3.5)) +
      xlim(c(-3.5,3.5)) + 
      xlab("Predicted metabolite levels, standartized value") +
      ylab("Observed metabolite levels, standartized value") +
      theme(aspect.ratio = 1)
plots.list = lappend(plots.list, p)

toPlot$metabolite = factor(toPlot$metabolite, levels = c("ATP", "ADP", "AMP"))
p = ggplot(toPlot, aes(x = yhat, y=y))+
      geom_point() +
      geom_text(data=stats.text, aes(x=x,y=y,label=round(adj.r.squared,2)))+
      facet_wrap(~metabolite, ncol = 1) +
      ylim(c(-3.5,3.5)) +
      xlim(c(-3.5,3.5)) + 
      xlab("Predicted metabolite levels, standartized value") +
      ylab("Observed metabolite levels, standartized value") +
      geom_smooth(method=lm,   # Add linear regression lines
                se=T,    # Don't add shaded confidence region
                fullrange=TRUE) +
      panel_border() +
      theme(aspect.ratio = 1)
plots.list = lappend(plots.list, p)





toPlot = all_final_models.models %>% ungroup() %>% filter(dataset != "PPP_AA", degree== 1, ismetIncluded == 0, the_super_best == T) %>% 
                                                          group_by(dataset, metabolite, the_super_best) %>% 
                                                          summarise(formula = unique(formula),
                                                                    median.cv.r2 = unique(median.cv.r2),
                                                                    adj.r.squared = unique(adj.r.squared),
                                                                    r.squared = unique(r.squared),
                                                                    model = unique(as.character(model)),
                                                                    type = unique(as.character(type)),
                                                                    datafile = unique(as.character(datafile)))

#metabolite order

metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]

toPlot = droplevels(toPlot)
toPlot$metabolite = factor(toPlot$metabolite, levels=as.character(metabolite.order$metabolite))

toPlot$met_name = metabolite2iMM904$official_name[match(toPlot$metabolite, metabolite2iMM904$id)]
toPlot$met_name = factor(toPlot$met_name, levels=rev(as.character(metabolite.order$met_name)))
toPlot$pathway = metabolite.order$pathway[match(toPlot$metabolite, metabolite.order$metabolite)]
toPlot$pathway = factor(toPlot$pathway, levels = as.character(metabolite.order$pathway))
p = ggplot(toPlot, aes(x = met_name, color=pathway)) +
    geom_point(data=toPlot, aes( y = adj.r.squared), colour = "black", size=5) +
    geom_point(data=toPlot, aes( y = median.cv.r2), colour="black", shape=17, size=5) +
    geom_linerange(data = toPlot , aes(ymin=0, ymax=adj.r.squared)) +
    coord_flip() + 
    background_grid(major = "x", minor = "none") +
    panel_border()
    
plots.list = lappend(plots.list, p)

p = ggplot(toPlot, aes(x = met_name, color=pathway)) +
    geom_point(data=toPlot, aes(y = adj.r.squared), colour = "black", size=5) +
    #geom_point(data=toPlot, aes( y = median.cv.r2), colour="black", shape=17, size=5) +
    geom_linerange(data = toPlot , aes(ymin=0, ymax=adj.r.squared)) +
    coord_flip() + 
    background_grid(major = "x", minor = "none") +
    panel_border()
plots.list = lappend(plots.list, p)

toPlot = toPlot %>% ungroup() %>% arrange(dataset, adj.r.squared)
toPlot$met_name = factor(toPlot$met_name, levels = toPlot$met_name)


p = ggplot(toPlot, aes(x = met_name)) +
  geom_point(data=toPlot, aes( y = adj.r.squared), colour = "black", size=3) +
  #geom_point(data=toPlot, aes( y = median.cv.r2), colour="black", shape=17, size=5) +
  geom_linerange(data = toPlot , aes(ymin=0, ymax=adj.r.squared)) +
  coord_flip() + 
  background_grid(major = "x", minor = "none") +
  panel_border()
plots.list = lappend(plots.list, p)


# saturation vs adj r2
# 
# brenda.summary = metabolites.long %>% group_by(dataset, X2, isPredictor) %>% summarise(median.ratio = median(log(ratio)))
# 
# models.summary = all_final_models.models %>% ungroup() %>% filter(dataset != "PPP_AA", the_super_best == T) %>% 
#                   group_by(dataset, metabolite, the_super_best) %>% 
#                   summarise(formula = unique(formula),
#                             median.cv.r2 = unique(median.cv.r2),
#                             adj.r.squared = unique(adj.r.squared),
#                             r.squared = unique(r.squared),
#                             model = unique(as.character(model)),
#                             type = unique(as.character(type)),
#                             datafile = unique(as.character(datafile)))
# 
# tmp = merge(brenda.summary, models.summary, by.x = c("dataset", "X2"), 
#                                             by.y = c("dataset", "metabolite" ))
# 
# tmp$isSaturated = ifelse(tmp$median.ratio > 0, 1, 0)
# 
# p = ggplot(tmp , aes(x = median.ratio, y = adj.r.squared, colour=isPredictor)) +
#             geom_point()
# 
# tmp %>% group_by(dataset, isPredictor) %>% summarise(cor(median.ratio,adj.r.squared, use="pairwise.complete" ))
# 
# pheatmap(cor(dataAA$metabolites, method="spearman"))
# cor.test(dataAA$metabolites[,"glutamate"],dataAA$metabolites[,"glutamine"], method="spearman")


### -- degree comparison ----
unique(all_final_models.models$dataset)

toPlot = all_final_models.models %>% filter(dataset %in%  c("TCA", "AA"), type=="after",  ismetIncluded == 0, the_super_best == T) %>% 
                                    group_by(dataset, metabolite, degree,  the_super_best) %>% 
                                    summarise(formula = unique(formula),
                                              median.cv.r2 = unique(median.cv.r2),
                                              adj.r.squared = unique(adj.r.squared),
                                              r.squared = unique(r.squared),
                                              model = unique(as.character(model)),
                                              type = unique(as.character(type)),
                                              datafile = unique(as.character(datafile)))
p = ggplot(toPlot, aes(x = adj.r.squared, fill = factor(degree))) + 
          geom_density(alpha = 0.5) +
          theme(aspect.ratio = 0.625) +
          panel_border()

plots.list = lappend(plots.list, p)


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")














