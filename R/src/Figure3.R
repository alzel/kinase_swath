rm(list = ls())
library(caret)
library(relaimpo)
source("./R/boot.R")
source("./R/functions.R")
library("cowplot")
library("scales")

### summarizes results from linear regression models


load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")
input_path = "./results/2016-02-24/linear_models"
load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/proteins.matrix.combat.RData")
load("./R/objects/exp_metadata._clean_.RData")

load("./R/objects/GO_slim.process._load_.RData")
load("./R/objects/GO.raw._load_.RData")
load("./R/objects/gene.annotations._load_.RData")

orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             sgd = gene.annotations$V5,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]


plots.list = list()

fun_name = "Figure3"

# lm models results overview ------------------

#filesToProcess = dir(path=input_path, pattern = "[12].[01]+.linear_models.RData$", recursive=F)
filesToProcess = dir(path=input_path, pattern = "[123].[01]+.linear_models.RData$", recursive=F)
filesToProcess = grep(pattern="imputed", invert=T, filesToProcess, value=T)
filesToProcess = grep(pattern="([1,3]+).([0-9]+).linear_models.RData", filesToProcess, value=T)

pattern.p = "data.(\\w+).(.*?).([0-9]+).([0-9]+).linear_models.RData$"

matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

read_models.aic = function(x) {
  z <<- x
  #x =  matches[[1]]
  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))
  
  table = my_models$summaries
  
  table$dataset = factor(x[[2]])
  table$species  = factor(x[[3]])
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = factor(x[[4]])
  table$ismetIncluded  = factor(x[[5]])
  table$file =  factor(x[[1]])
  table$total_predictors  = dim(my_models$input_data)[2]-1
  return(table)
}

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
auto_thr.bonferonni = auto_thr/length(unique(all_final_models.aic$formula))
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



toPlot = models.filtered %>% ungroup() %>% arrange(dataset, adj.r.squared)

toPlot$metabolite = factor(toPlot$metabolite, levels = unique(toPlot$metabolite))
p = ggplot(toPlot, aes(x = metabolite, y = adj.r.squared)) +
  geom_point() +
  geom_linerange(data = toPlot , aes(x=metabolite,ymin=0, ymax=adj.r.squared))+
  facet_wrap(ismetIncluded~dataset+degree, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




## summarizing protein table ----------------------
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

proteins.log.quant = t(my_means(proteins.matrix.combat.quant))
proteins.raw = exp(t(my_means(proteins.matrix.combat)))
proteins.log = t(my_means(proteins.matrix.combat))

read_models.aic.predict = function(x) {
  z <<- x
  x = z
  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))
  
  models = my_models$models
  input_data = my_models$input_data[,-1]
  trans.x <- my_models$trans.x
  new_data = NULL
  
  if (grepl(pattern = "log.quant", x = x[[1]])) {
    message("log.quant")
    new_data <- proteins.log.quant[,colnames(input_data)]
  } else if (grepl(pattern = "log", x = x[[1]])) {
    message("log")
    new_data <- proteins.log[,colnames(input_data)]
  } else {
    message("raw") #if raw data, applying boxcox transformation model as in input_data
    new_data <- proteins.raw[,colnames(input_data)]
    tmp.matrix = new_data
    trans.BC = preProcess(x = input_data , method = c("BoxCox"))
    new_data = predict(trans.BC, tmp.matrix)
  }

  #getting minimum and maximums for every proteins
  mins = apply(scale(proteins.raw, center = T, scale = T), 2, which.min)  
  maxes = apply(scale(proteins.raw, center = T, scale = T), 2, which.max)
  kinase_deletions = rownames(proteins.raw)
  
  if (!is.null(trans.x)) {
    new_data <- new_data[,names(trans.x$mean)]
    new_data <- scale(new_data, center = trans.x$mean, scale = trans.x$std) %*% trans.x$rotation #all data projected into pca space
  } else {
    new_data <- scale(new_data, center = T, scale = T)
  }
 
  models.list = list()
  for (i in 1:length(models)) {
    for (j in c("before", "after")) {
      top = 3
      current_model <- models[[i]][[j]]
      predictors <- coefficients(current_model)[-1]
      
      predictions <- predict(current_model, as.data.frame(new_data))
      expanded.model = NULL
      
      if (!is.null(trans.x) & any(grepl(x = names(predictors), pattern = "PC\\d+", perl=T))) {
        
        tmp.long <- melt(trans.x$rotation[,names(predictors[-1])], id.varr="rownames")
        #tmp.long$type <- j
        tmp.long$weight <- predictors[match(tmp.long$X2, names(predictors))]
        
        expanded.model <- tmp.long %>%  # selecting per component top loadings
          arrange(X2) %>%
          group_by(X2) %>%
          arrange(desc(value)) %>% 
          mutate(toSelect  = ifelse(value %in% head(value,top), 1, 
                                    ifelse(value %in% tail(value,top), 1, 0))) %>%
          filter(toSelect == 1) %>% 
          dplyr::select(-toSelect) %>%
          group_by(X1) %>%
          mutate(kinase.max = kinase_deletions[maxes[match(X1, names(maxes))]],
                 kinase.min = kinase_deletions[mins[match(X1, names(mins))]],
                 kinase.max.effect = predictions[match(kinase.max, names(predictions))],
                 kinase.min.effect = predictions[match(kinase.min, names(predictions))]) %>% ungroup()
        
        expanded.model <- expanded.model[gtools::mixedorder(as.character(expanded.model$X2)),] 
        expanded.model <- expanded.model %>% group_by(X1) %>% filter(row_number() == 1) #selecting highest loading based on component rank
        
      } else {
        expanded.model <- data.frame(X1 = names(predictors), X2 = NA, value = NA) %>%
        #expanded.model <- data.frame(X1 = names(predictors), X2 = NA, value = NA, type = j) %>% 
                   mutate(weight = predictors,
                          kinase.max = kinase_deletions[maxes[match(X1, names(maxes))]],
                          kinase.min = kinase_deletions[mins[match(X1, names(mins))]],
                          kinase.max.effect = predictions[match(kinase.max, names(predictions))],
                          kinase.min.effect = predictions[match(kinase.min, names(predictions))])
      }
      
      tmp <- cbind.data.frame(expanded.model, my_models$summaries  %>% filter(model == i, type == j))
      models.list <- lappend(models.list, tmp)
    }
  }
  
  table = do.call(rbind.data.frame, models.list)
  table$dataset = factor(x[[2]])
  table$species  = factor(x[[3]])
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = factor(x[[4]])
  table$ismetIncluded  = factor(x[[5]])
  table$file =  factor(x[[1]])
  table$total_predictors  = dim(my_models$input_data)[2]-1
  return(table)
}

file.list = lapply(matches, FUN=read_models.aic.predict)
all_final_models.models = do.call(rbind.data.frame, file.list) #stores MLR models including min/max kinases

all_final_models.models$metabolite = all_final_models.models$species
all_final_models.models$metabolite = sub(x=all_final_models.models$metabolite, pattern="log.quant.(.*)", replacement="\\1")
all_final_models.models$metabolite = sub(x=all_final_models.models$metabolite, pattern="log.(.*)", replacement="\\1")

all_final_models.models$normalization = "bc"
all_final_models.models$normalization[grep(pattern="log", x=all_final_models.models$species)] = "log"
all_final_models.models$normalization[grep(pattern="log.quant", x=all_final_models.models$species)] = "log.quant"
all_final_models.models$normalization = factor(all_final_models.models$normalization)

auto_thr = 0.05
auto_thr.bonferonni = auto_thr/length(unique(all_final_models.models$formula))
all_final_models.models$isAutocorrelation.bonferoni = ifelse(all_final_models.models$bg.p.value < auto_thr.bonferonni, 1, 0)
all_final_models.models$isAutocorrelation = ifelse(all_final_models.models$bg.p.value < auto_thr, 1, 0)

all_linear_models <- all_final_models.models %>% 
  filter(type == "after",  
        isImputed == 0 , isAutocorrelation.bonferoni == 0) %>%
      group_by(metabolite, ismetIncluded, degree) %>% 
        mutate(the_super_best = adj.r.squared == max(adj.r.squared)) #the best model

metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]

toPlot = all_linear_models %>% ungroup() %>% 
  filter(metabolite %in% metabolite.order$metabolite, 
         degree== 1, ismetIncluded == 0, the_super_best == T, isImputed == 0 , isAutocorrelation.bonferoni == 0) %>% 
  group_by(metabolite, the_super_best) %>% 
  summarise(formula = unique(formula),
            median.cv.r2 = unique(median.cv.r2),
            adj.r.squared = unique(adj.r.squared),
            r.squared = unique(r.squared),
            model = unique(as.character(model)),
            type = unique(as.character(type)),
            datafile = unique(as.character(datafile)),
            dataset = unique(as.character(dataset)),
            total_predictors = unique(total_predictors))


toPlot = droplevels(toPlot)
toPlot$metabolite = factor(toPlot$metabolite, levels=as.character(metabolite.order$metabolite))

toPlot$met_name = metabolite.order$met_name[match(toPlot$metabolite, metabolite.order$metabolite)]
toPlot$met_name = factor(toPlot$met_name, levels=rev(as.character(metabolite.order$met_name)))
toPlot$pathway = metabolite.order$pathway[match(toPlot$metabolite, metabolite.order$metabolite)]
toPlot$pathway = factor(toPlot$pathway, levels = as.character(metabolite.order$pathway))
aic.models <- toPlot

library(ggthemes)

p.all_met = ggplot(toPlot, aes(x = met_name, color=pathway)) +
  #geom_point(data=toPlot, aes( y = adj.r.squared), colour = "black", size=2) +
  #geom_point(data=toPlot, aes( y = median.cv.r2), colour="black", shape=17, size=5) +
  geom_linerange(data = toPlot , aes(ymin=0, ymax=adj.r.squared),size=3) +
  coord_flip() + 
  scale_colour_tableau("tableau20") +
  background_grid(major = "x", minor = "none") +
  panel_border() +
  theme(legend.justification=c(1,0), legend.position=c(-0.1,0)) +
  ylab("Explained variance of metabolite concentrations\nusing proteome data, adj. R2")

message(paste("Total number of models", format(sum(2^(toPlot %>% filter(metabolite != "ATP"))$total_predictors))))

# total theoretical number of models
sum(sapply((toPlot %>% filter(metabolite != "ATP"))$total_predictors, function(x)
  sum(choose(x, 1:x))))



# plotting graph -----------------
library(GGally)
library(network)
library(sna)


graph_dataset <- all_linear_models %>% 
  filter(ismetIncluded == 0, degree == 1, the_super_best == T, adj.r.squared > 0.3, p.value <0.01, weight > 0.2) %>%
  ungroup() %>%
  dplyr::select(X1, weight, kinase.max, kinase.min, kinase.max.effect,kinase.min.effect, metabolite, adj.r.squared)

met.graph_dataset <- graph_dataset %>% dplyr::select(X1, metabolite, weight)
met.graph_dataset$type <- "enz->met"
names(met.graph_dataset) <- c("from", "to", "weight", "type")

kinase.graph_dataset.min <- graph_dataset %>% dplyr::select(kinase.min, X1, kinase.min.effect)
names(kinase.graph_dataset.min) <- c("from", "to", "weight")
kinase.graph_dataset.min$type <- "kinase.min->enz"

kinase.graph_dataset.max <- graph_dataset %>% dplyr::select(kinase.max, X1, kinase.max.effect)
names(kinase.graph_dataset.max) <- c("from", "to", "weight")
kinase.graph_dataset.max$type <- "kinase.max->enz"


enzymes.tmp <- data.frame(id = graph_dataset$X1, 
                          adj.r.squared = NA, type = "enzyme", 
                          node_name = orf2name$gene_name[match(graph_dataset$X1, orf2name$ORF)]) %>% distinct()

kinases.tmp <- data.frame(id = c(graph_dataset$kinase.max, graph_dataset$kinase.min), 
                          adj.r.squared = NA, 
                          type = "kinase",
                          node_name = orf2name$gene_name[match(c(graph_dataset$kinase.max, graph_dataset$kinase.min), orf2name$ORF)]) %>% distinct()

metabolites.tmp <- data.frame(id = graph_dataset$metabolite, 
                              adj.r.squared = graph_dataset$adj.r.squared, 
                              type = "metabolite",
                              node_name = metabolite2iMM904$official_name[match(graph_dataset$metabolite, metabolite2iMM904$id)]) %>% distinct()


nodes <- rbind.data.frame(metabolites.tmp, enzymes.tmp, kinases.tmp)
edges <- rbind.data.frame(met.graph_dataset, kinase.graph_dataset.max, kinase.graph_dataset.min)
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T) 
colrs <- c("gray50", "tomato", "gold")

V(net)$color <- colrs[as.numeric(as.factor(V(net)$type))]
V(net)$degree.out <- igraph::degree(net, v=V(net),mode="out")

V(net)$size[V(net)$type == "metabolite"] <- as.numeric(cut(V(net)$adj.r.squared[V(net)$type == "metabolite"],5))
V(net)$size[V(net)$type == "kinase"] <- as.numeric(cut(V(net)$degree.out[V(net)$type == "kinase"],5))
V(net)$size[V(net)$type == "enzyme"] <- 1

V(net)$layer <- as.numeric(factor(V(net)$type, levels = c("kinase", "enzyme", "metabolite")))
V(net)$name <- V(net)$node_name


layout.k_partite <- function(g) {
  
  l <- layout.sugiyama(g)$layout[,2:1]
  l[,1] <- V(g)$layer
  l[,2] <- - l[,2] + 1 + max(l[,2])
  l
}
node_sizes <- data.frame(node = V(net)$name,
                         size = V(net)$size)

my_layout <-  layout.k_partite(net)
new_layout <- sweep(my_layout, MARGIN = 2, STATS = apply(my_layout,2, max), FUN = "/")

B = network(igraph::get.adjacency(net), directed = TRUE)

B %v% "x" <- new_layout[,1]
B %v% "y" <- new_layout[,2]
B %v% "color" = V(net)$color

B.nodes = network.vertex.names(B)
B.sizes = node_sizes$size[match(B.nodes, node_sizes$node)]
B %v% "size" = B.sizes
B %e% "effect" <- ifelse(E(net)$weight>0, "red", "green") 

p.graph = ggnet2(B, mode = c("x", "y"), label.size = 4, edge.alpha = 0.7,
           size = "size", 
           label = T, 
           color = "color", 
           edge.color = "effect") + guides(size = FALSE)




# plot(net, vertex.label.dist=0, edge.arrow.size=.08,
#      edge.curved=0, vertex.size=5,
#      vertex.color=V(net)$color, vertex.frame.color="white",
#      vertex.label=V(net)$node_name, vertex.label.color="black",
#      vertex.label.cex=.7, layout = my_layout)
# 
# file_name = paste("all_linear_models", fun_name, "graphml", sep=".")
# file_path = paste(output_dir, file_name, sep="/")
# write.graph(net, file=file_path, format="graphml")

# ----- glutamine example ------------------------

met = "glutamine"

stopifnot(length(met)==1)
example.models = all_linear_models %>% ungroup() %>% filter(degree == 1, ismetIncluded == 0, the_super_best == T, metabolite %in% met) %>% 
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
  #x = tmp.list[[1]]
  file_name = paste(input_path,x[[4]], sep="/") 
  my_models = get(load(file_name))
  
  models = my_models$models
  
  
  fit = models[[as.numeric(x[[2]])]][[x[[3]]]]
  yhat = predict(fit, fit$model[,-1])
  
  
  response.trans  = my_models$trans
  predictions.untransformed <- yhat
  y.untransformed <- fit$model[,1]
  
  if(!is.null(response.trans$std)) {
    predictions.untransformed <- predictions.untransformed*response.trans$std[1]
    y.untransformed <- y.untransformed*response.trans$std[1]
  }
  
  if(!is.null(response.trans$mean)) {
    predictions.untransformed <- predictions.untransformed + response.trans$mean[1]
    y.untransformed <- y.untransformed + response.trans$mean[1]
  }
  
  if(!is.null(response.trans$bc)) {
    predictions.untransformed <- inverse.BoxCoxTrans(response.trans$bc[[1]], predictions.untransformed)
    y.untransformed <- inverse.BoxCoxTrans(response.trans$bc[[1]], y.untransformed)
  }
  
  
  table = data.frame(metabolite = x[[1]],
                     model = x[[2]],
                     type = x[[3]],
                     file = x[[4]],
                     y = fit$model[,1],
                     yhat = yhat,
                     adj.r.squared = x[[5]],
                     median.cv.r2 =  x[[6]],
                     predictors = paste(names(coefficients(fit)[-1]), collapse = ":"),
                     coefficients = paste(coefficients(fit)[-1], collapse = ":"),
                     sample_name = rownames(fit$model),
                     predictions.untransformed = predictions.untransformed, 
                     y.untransformed = y.untransformed)
  return(table)
}


prediction.list = lapply(tmp.list, FUN=read_models)
prediction.models = do.call(rbind.data.frame, prediction.list)

prediction.models <- prediction.models %>% group_by(metabolite) %>% arrange(y)
prediction.models$gene = exp_metadata$gene[match(prediction.models$sample_name, exp_metadata$sample_name)]

prediction.models <- prediction.models %>% 
  group_by(metabolite) %>% 
    mutate(ntile20 = ntile(yhat, 20))
  
toMark <- prediction.models %>%
  group_by(ntile20) %>% filter(ntile20 %in% c(2,10, 19), row_number() == 1) %>% 
  bind_rows(prediction.models %>% 
              group_by(metabolite, gene) %>% filter(gene == "WT", row_number() == 1) ) 


toPlot = prediction.models
stats.text = prediction.models %>% group_by(metabolite, model) %>% summarise(adj.r.squared = as.numeric(as.character(adj.r.squared[1])),
                                                                             median.cv.r2 = as.numeric(as.character(median.cv.r2[1])))
stats.text$x = -2
stats.text$y = seq(2,1,length.out=nrow(stats.text))
#toPlot$metabolite = factor(toPlot$metabolite, levels = c("ATP", "ADP", "AMP"))
p.met = ggplot(toPlot, aes(x = yhat, y=y))+
  geom_point() +
  geom_text(data=stats.text, aes(x=x,y=y,label=round(adj.r.squared,2))) +
  geom_point(data = toMark, aes(x = yhat, y = y), colour="red") +
  geom_text(data = toMark, aes(x = yhat, y = y, label=gene), check_overlap = TRUE) +
  ylim(c(-2.2,2.2)) +
  xlim(c(-2.2,2.2)) + 
  xlab("Predicted metabolite levels,\nstandartized value") +
  ylab("Observed metabolite levels,\nstandartized value") +
  geom_smooth(method=lm,   # Add linear regression lines
              se=F,    # Don't add shaded confidence region
              fullrange=F) +
  panel_border() +
  theme(aspect.ratio = 1)

stats.text$x = sort(prediction.models$predictions.untransformed)[2]
stats.text$y = sort(prediction.models$predictions.untransformed, decreasing = T)[2]

p.met.untransformed = ggplot(toPlot, aes(x = predictions.untransformed, y = y.untransformed))+
  geom_point() +
  geom_text(data=stats.text, aes(x = x,y = y,label=round(adj.r.squared,2))) +
  geom_point(data = toMark, aes(x = predictions.untransformed, y = y.untransformed), colour="red") +
  geom_text(data = toMark, aes(x = predictions.untransformed, y = y.untransformed, label=gene), check_overlap = TRUE) +
  xlab("Predicted metabolite concentration, µM") +
  ylab("Observed metabolite concentration, µM") +
  geom_smooth(method=lm,   # Add linear regression lines
              se=F,    # Don't add shaded confidence region
              fullrange=F) +
  panel_border() +
  theme(aspect.ratio = 1)

#file_name = paste(input_path,x[[4]], sep="/") 
my_models = get(load(file_name))
models = my_models$models


proteins.log.quant = t(proteins.matrix.combat.quant)
proteins.raw = exp(t(proteins.matrix.combat))
proteins.log = t(proteins.matrix.combat)

#returns data for example plot
toMark_proteins <- ddply(toMark, .variables = .(metabolite, file, sample_name, gene), 
      .fun = function(x) {
        
        z <<- predictors <- unlist(strsplit(unique(as.character(x$predictors)), split=":"))

        if (grepl(pattern = "log.quant", x = unique(as.character(x$file)))) {
          message("log.quant")
          new_data <- proteins.log.quant[,predictors]
          tmp.matrix = new_data
          trans.BC = preProcess(x = new_data, method = c("center", "scale"))
          new_data = predict(trans.BC, tmp.matrix)
        } else if (grepl(pattern = "log", x = unique(as.character(x$file)))) {
          message("log")
          new_data <- proteins.log[,predictors]
          tmp.matrix = new_data
          trans.BC = preProcess(x = new_data, method = c("center", "scale"))
          new_data = predict(trans.BC, tmp.matrix)
        } else {
          message("raw") #if raw data, applying boxcox transformation model as in input_data
    
          new_data <- proteins.raw[,predictors]
          tmp.matrix = new_data
          trans.BC = preProcess(x = new_data, method = c("BoxCox", "center", "scale"))
          new_data = predict(trans.BC, tmp.matrix)
        }
        return(new_data[rownames(new_data) %in% x$sample_name,])
      }
)



tmp.data <- all_linear_models %>% filter(ismetIncluded == 0, degree == 1, the_super_best == T, metabolite == met)


first_row <- toMark_proteins[1,]
first_row[, 5:length(first_row)] <- 0
first_row$gene = "effect"
first_row$sample_name = "unknown"


toPlot <- melt(rbind(toMark_proteins,first_row), id.vars=c("metabolite", "file", "sample_name", "gene"))

toPlot_coefs<- data.frame(variable = unique(as.character(toPlot$variable)))
toPlot_coefs$weight <- tmp.data$weight[match(toPlot_coefs$variable, tmp.data$X1)]
toPlot_coefs$effector <- ifelse(toPlot_coefs$weight >0, "positive", "negative")

toPlot_coefs <- toPlot_coefs %>% arrange(desc(abs(weight)))
toPlot_coefs$variable <- factor(toPlot_coefs$variable, levels = as.character(toPlot_coefs$variable))
toPlot_coefs$variable_name <- orf2name$gene_name[match(toPlot_coefs$variable, orf2name$ORF)]
toPlot_coefs$variable_name <- factor(toPlot_coefs$variable_name, levels = as.character(toPlot_coefs$variable_name))

toPlot$variable = factor(toPlot$variable, levels = rev(as.character(toPlot_coefs$variable )))
toPlot$variable_name <- orf2name$gene_name[match(toPlot$variable, orf2name$ORF)]
toPlot$variable_name = factor(toPlot$variable_name, levels = rev(as.character(toPlot_coefs$variable_name )))

library(cowplot)

p.enzymes <- ggplot() + 
  geom_tile(data = toPlot, aes(y=variable_name, x=gene, fill=value)) +
  scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                       midpoint = 0, limit = c(-2,2), space = "Lab", 
                       name="Protein expression") +
  ylab("Effector enzyme") +
  xlab("Kinase knockout") +
  geom_point(data=toPlot_coefs, aes(y=variable_name, x=5, size=abs(weight), colour=effector)) +
  theme(legend.position="top", aspect.ratio = 1) 

p.enzymes <- ggdraw(switch_axis_position(p.enzymes , axis = 'x'))

# data_file = sub(x = unique(as.character(p.met$data$file)), pattern = "linear_models.", replacement = "")
# metabolite_data <- get(load(paste("./results/2016-02-24/data.AA/", data_file, sep = "")))
# colnames(metabolite_data) = c(met, orf2name$gene_name[match(colnames(metabolite_data), orf2name$ORF)][-1])

#GLN1 example

x <- tmp.list
file_name = paste(input_path,x[[1]][4], sep="/") 
my_models = get(load(file_name))
models = my_models$models
fit = models[[as.numeric(x[[1]][2])]][x[[1]][3]]
used.samples <- rownames(fit$after$model)
data.tmp = my_models$input_data
data.tmp <- data.tmp[rownames(data.tmp) %in% used.samples,]

toPlot <- data.tmp[, c(met, orf2name$ORF[orf2name$gene_name == "GLN1"])]
names(toPlot) <- c("y", "x")

p.gln1 <- ggplot(toPlot, aes(x =x , y = y) ) +
  geom_point() +
  xlab("Predicted metabolite concentration, µM") +
  ylab("Observed metabolite concentration, µM") +
  geom_smooth(method=lm,   # Add linear regression lines
              se=F,    # Don't add shaded confidence region
              fullrange=F) +
  panel_border() +
  theme(aspect.ratio = 1)
plots.list <- lappend(plots.list, p.gln1)

# -- Energy metabolite example ---------------

met = c("ATP", "ADP", "AMP")

example.models = all_linear_models %>% ungroup() %>% filter(degree == 1, ismetIncluded == 0, the_super_best == T, metabolite %in% met) %>% 
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
  #x = tmp.list[[1]]
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
                     median.cv.r2 =  x[[6]],
                     predictors = paste(names(coefficients(fit)[-1]), collapse = ":"),
                     coefficients = paste(coefficients(fit)[-1], collapse = ":"),
                     sample_name = rownames(fit$model))
  return(table)
}

read_predictions = function(x) {
  #x = tmp.list[[2]]
  file_name = paste(input_path,x[[4]], sep="/") 
  my_models = get(load(file_name))
  
  models = my_models$models
  
  fit = models[[as.numeric(x[[2]])]][[x[[3]]]]
  
  M = cor(fit$model[,-1])
  L = chol(M)
  nvars = dim(L)[1]
  n_obs = 1000
  
  r = t(L) %*% matrix(rnorm(nvars*n_obs, mean = 0, sd = 1), nrow=nvars, ncol=n_obs)
  r = as.data.frame(t(r))
  names(r) <- names(fit$model[,-1])
  
  #yhat_real = predict(fit, fit$model[,-1])
  
  yhat = predict(fit, r)
  conf = predict(fit, r, interval = "confidence", level = 0.95)
  pred = predict(fit, r, interval = "prediction", level = 0.95)
  
  pred.int <- as.data.frame(cbind(yhat,pred))
  conf.int <- as.data.frame(cbind(yhat,conf))
  
#    g.pred <- ggplot(pred.int, aes(x = yhat, y = fit)) +
#     theme_bw() +
#     ggtitle("Prediction interval for future observations from predict()") +
#     geom_point(data = xy, aes(x = x, y = y)) +
#     geom_smooth(data = pred.int, aes(ymin = lwr, ymax = upr), stat = "identity") 
#     #geom_smooth(data = conf.int, aes(ymin = lwr, ymax = upr), stat = "identity", fill="red")

  
   table = data.frame(metabolite = x[[1]],
                     model = x[[2]],
                     type = x[[3]],
                     file = x[[4]],
                     yhat = yhat,
                     conf.lwr = conf.int$lwr,
                     conf.upr = conf.int$upr,
                     pred.lwr = pred.int$lwr,
                     pred.upr = pred.int$upr,
                     adj.r.squared = x[[5]],
                     median.cv.r2 =  x[[6]] )
  return(table)
}


prediction.list = lapply(tmp.list, FUN=read_models)
prediction.models = do.call(rbind.data.frame, prediction.list)

prediction.models <- prediction.models %>% group_by(metabolite) %>% arrange(y)
prediction.models$gene = exp_metadata$gene[match(prediction.models$sample_name, exp_metadata$ORF)]


prediction.models <- prediction.models %>% 
  group_by(metabolite) %>% 
  mutate(ntile20 = ntile(yhat, 20))

simulations.list = lapply(tmp.list, FUN=read_predictions)
prediction.intervals = do.call(rbind.data.frame, simulations.list)


toPlot = prediction.models
stats.text = prediction.models %>% group_by(metabolite, model) %>% summarise(adj.r.squared = as.numeric(as.character(adj.r.squared[1])),
                                                                             median.cv.r2 = as.numeric(as.character(median.cv.r2[1])))
stats.text$x = -1
stats.text$y = seq(2,1,length.out=nrow(stats.text))
#toPlot$metabolite = factor(toPlot$metabolite, levels = c("ATP", "ADP", "AMP"))

toPlot$metabolite = factor(toPlot$metabolite, levels = c("ATP", "ADP", "AMP"))
prediction.intervals$metabolite = factor(prediction.intervals$metabolite, levels = c("ATP", "ADP", "AMP"))

p.energy = ggplot(prediction.intervals)+
  geom_point(data = toPlot, aes(x = yhat, y = y) ) +
  geom_text(data=stats.text, aes(x=x,y=y,label=round(adj.r.squared,2)))+
  facet_wrap(~metabolite, ncol = 1, scale="free") +
  ylim(c(-5.5,5.5)) +
  xlim(c(-2.5,2.5)) + 
  xlab("Predicted metabolite levels,\nstandartized value") +
  ylab("Observed metabolite levels, standartized value") +
  geom_smooth(data = prediction.intervals, aes(x = yhat, y = yhat, ymin = pred.lwr, ymax = pred.upr), stat = "identity") + 
  panel_border() +
  theme(aspect.ratio = 1)

# -- tRNAs ------
yeast.model = iMM904
AA.linear_models = all_linear_models %>% ungroup() %>% 
  filter(metabolite %in% metabolite.order$metabolite, 
         degree== 1, ismetIncluded == 0, isImputed == 0,
         isAutocorrelation.bonferoni == 0, dataset == "AA")


proteogenic_AA = metabolite2iMM904[grep(x = metabolite2iMM904$model_name, pattern = "-L|gly"),]  %>% 
  filter(!(id %in% c("Malate", "citrulline","homo.cysteine", "homo.serine"))) 

AA.linear_models$isProteogenic <- ifelse(AA.linear_models$metabolite %in% proteogenic_AA$id, 1, 0)

tRNAs <- right_join(orf2name, GO.raw[grep(pattern = "tRNA", x = GO.raw$V10),] %>%  dplyr::select(V2) %>% distinct(), by=c("sgd" = "V2"))

AA.linear_models$istRNA <- ifelse(AA.linear_models$X1 %in% tRNAs$ORF, 1, 0)

protAA = AA.linear_models %>% filter(isProteogenic == 1) %>% dplyr::select(metabolite) %>% distinct()

protAA_model_name <- metabolite2iMM904$model_name[match(protAA$metabolite, metabolite2iMM904$id)]

model.measured.AA <- yeast.model[yeast.model$metabolite %in% protAA_model_name,]
model.measured.AA$istRNA <- ifelse(model.measured.AA$gene %in% tRNAs$ORF, 1, 0)


model.measured.AA %>% 
  group_by(metabolite) %>%
  summarize(sum(istRNA))

#metaboblites with tRNA as predictors
AA.linear_models[AA.linear_models$X1 %in% tRNAs$ORF,] %>% filter(adj.r.squared > 0.25) %>% dplyr::select(metabolite) %>% distinct()

# --- BRENDA metabolites ---------------

yeast.model = iMM904
yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters
edgelist = unique(droplevels(na.omit(subset(yeast.model, metabolite != "h"  & metabolite !="h2o" & metabolite != "nh4" & metabolite != "ppi" , select = c("metabolite", "gene")))))
B <- graph.data.frame(edgelist)
V(B)$type <- V(B)$name %in% edgelist$metabolite
stopifnot(igraph::is.bipartite(B))

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
                   
                   SUB = igraph::induced.subgraph(B, unique(unlist(igraph::neighborhood(B, order=2, nodes=current_nodes))))    
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


## -- R2 is independent of coverage ----
toPlot = merge(models.filtered %>% filter(degree == 1), coverage.stats.wide)
res = cor(toPlot$adj.r.squared, toPlot$gene)
res = data.frame(res)
p.coverage_r2 = ggplot(toPlot, aes(x=gene, y=adj.r.squared)) +
  geom_text(data = res, aes(x=0.25, y=0.5, label=paste("rho =", round(res,3)))) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  xlab("Coverage of measured metabolizing enzymes, x100%") +
  ylab("Proportion of metabolite concentration variance explained by proteome data, adjusted R2") +
  theme(aspect.ratio = 1)

plots.list <- lappend(plots.list, p.coverage_r2)


brenda <- read.delim("./data/2015-10-07/brenda.txt")
load("./R/objects/proteins.matrix.combat.RData")

ec.gene = unique(gene.annotations[gene.annotations$V3 == "EC number",c(1,4)])



#coverage of measured EC 
coverage.ec = droplevels(coverage[coverage$type == "gene",])
coverage.ec = merge(coverage.ec, ec.gene, by.x = "neighbour", by.y = "V4")
coverage.ec$kegg_id = metabolite2iMM904$kegg_id[match(coverage.ec$metabolite, metabolite2iMM904$id)]
brenda[grep(pattern = "pro", x = brenda$substrate, ignore.case = T),] %>% View()

brenda.f = brenda[!brenda$KEGGID == "",]
brenda.f = brenda.f[grep(pattern="mutant|recombinant", x= brenda.f$commentary, invert=T),]

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



models.summary = all_linear_models %>% filter(degree==1, ismetIncluded == 0,  the_super_best)

models.summary = models.summary[models.summary$X1 != "stats",]
models.summary = merge(models.summary, ec.gene, by.x = "X1", by.y = "V4")
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

tRNA_predictors <- AA.linear_models %>% 
  filter(X1 %in% tRNAs$ORF) %>% 
  dplyr::select(metabolite, X1) %>% distinct() %>% 
  left_join(ec.gene, by = c("X1" = "V4"))

names(tRNA_predictors) <- c("metabolite", "gene_name", "ec")
  

# checking for saturation


metabolites.long = metabolites.long %>% mutate(ratio = concentration/kmValue)
toPlot = metabolites.long %>% filter(inNetwork == T)
points = metabolites.long %>% filter(isPredictor == T)

points2 = inner_join(points, tRNA_predictors, by = c("X2" = "metabolite", "ecNumber" = "ec"))

p = ggplot(toPlot, aes(y=log(ratio), x = label)) +  
  geom_boxplot() +
  geom_point(data=points, aes(y=log(ratio), x = label), col="red") +
  geom_hline(yintercept = 0) +
  facet_wrap(~dataset, scales = "free", ncol=1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toPlot <- points2
p.trna <- ggplot(toPlot, aes(x=log(ratio), fill = label)) +  
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0) +
  scale_fill_discrete(name = "") +
  xlab(expression(paste("Measured intracellular concentraions divided by ", K[M], " of tRNA charging enzymes, ln(ratio)"))) +
  theme_bw() +
  theme(legend.position = c(0.6, 0.7), legend.background = element_rect(color = NULL), legend.direction = "horizontal")

ordering <- toPlot %>% group_by(label) %>% 
  mutate(median_ratio = median(log(ratio), na.rm=T)) %>%  ungroup() %>%
  arrange(median_ratio) %>% dplyr::select(label) %>% distinct()

toPlot$label <- factor(toPlot$label, levels = rev(ordering$label))
p.trna_separate <- ggplot(toPlot, aes(x=log(ratio))) +  
  geom_density(fill = "black") +
  geom_vline(xintercept = 0) +
  facet_wrap(~label, scales = "free_y", ncol = 1) +
  scale_fill_discrete(name = "") +
  xlab(expression(paste("Measured intracellular concentraions divided by ", K[M], " of tRNA charging enzymes, ln(ratio)"))) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



# predictors AA vs TCA
toPlot = metabolites.long %>% filter(inNetwork == T, isPredictor == T)
#toPlot = metabolites.long
stats = data.frame(label_text = c(median(toPlot$ratio[toPlot$dataset == "AA"], na.rm=T)/median(toPlot$ratio[toPlot$dataset == "TCA"], na.rm=T),
                                  wilcox.test(log(toPlot$ratio[toPlot$dataset == "AA"]), log(toPlot$ratio[toPlot$dataset == "TCA"]))$p.value,
                                  sum(log(toPlot$ratio,2) > 0, na.rm=T)/length(toPlot$ratio)), #associated enzymes above Km
                   
                   x = -7,
                   y = c(0.15, 0.13, 0.1))

toPlot$dataset = factor(toPlot$dataset)

super_saturated <- toPlot %>% group_by(label) %>%
  summarise(super_saturated = sum(log(ratio,10)>1)) 

sum(super_saturated$super_saturated>0)/length(unique(super_saturated$label))


p.aatca <- ggplot() +  
  geom_density(data=toPlot, aes(x=log(ratio,2), fill=dataset), alpha = 0.5) +
  geom_text(data=stats, aes(x=x+2, y=y, label = label_text)) +
  xlab(expression(paste("Measured intracellular concentraions divided by ", K[M], " of predictor enzymes, ln(ratio)"))) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  theme( legend.position = c(0.7,0.7))
# 
# ggplot(toPlot, aes(x = concentration)) +
#   geom_density() +
#   facet_wrap(~label, scales = "free") +
#   geom_point(data = toPlot, aes(x = kmValue, y=0.1))

# ### ratio whether it is predictor of not
a = (metabolites.long %>% filter(inNetwork == T, isPredictor == F) %>% dplyr::select(ratio))$ratio
b = (metabolites.long %>% filter(inNetwork == T, isPredictor == T) %>% dplyr::select(ratio))$ratio
# 
toPlot = metabolites.long %>% filter(inNetwork == T)
stats = data.frame(label_text = c(median(b,na.rm=T)/median(a, na.rm=T),
                                   wilcox.test(b,a)$p.value),
                    x = 1,
                    y = c(5,4))
 
p.predictors = ggplot(toPlot, aes(x = isPredictor, y = log(ratio))) + 
 geom_boxplot(width=0.2)+
 geom_violin(alpha=0)+
 geom_text(data=stats, aes(x=x, y=y, label = label_text)) +
 panel_border() +
 theme(aspect.ratio = 8/3)



#geom_point(data=toPlot, aes(x=jitter(as.numeric(isPredictor)) + 1 , y = log(ratio)), alpha = 0.1)


# ---- comparing network degrees ------ 

load("./R/objects/all_final_models.models_summary2.RData")
load("./R/objects/file.list.models_summary2.RData")

metabolites.models.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  dplyr::select(model, RMSE, Rsquared, normalization, dataset, metabolite, degree, preprocessing) %>% 
  distinct() %>%
  group_by(metabolite, normalization, degree, preprocessing) %>%
  #group_by(metabolite) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite, degree) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))

length(unique(all_final_models$model))

metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]


toPlot <- metabolites.models.long %>% filter(metabolite %in% metabolite.order$metabolite)

toPlot$metabolite = factor(toPlot$metabolite, levels=as.character(metabolite.order$metabolite))
toPlot$met_name = metabolite.order$met_name[match(toPlot$metabolite, metabolite.order$metabolite)]
toPlot$met_name = factor(toPlot$met_name, levels=as.character(metabolite.order$met_name))
toPlot$pathway = metabolite.order$pathway[match(toPlot$metabolite, metabolite.order$metabolite)]
toPlot$pathway = factor(toPlot$pathway, levels = as.character(metabolite.order$pathway))

p.cv.all <- ggplot(toPlot, aes(x = met_name, y = Rsquared, fill=degree)) +
  geom_bar(position="dodge", stat = "identity") +
  ylab(expression(paste("Cross-validated ", R^2, sep=""))) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.2, 0.7)) +
  background_grid(major = "y", minor = "y") +
  panel_border() 

p.cv.density <- ggplot(toPlot, aes(x = Rsquared, fill=degree)) +
  geom_density(aes(alpha = 0.5))

tests.pvals <- toPlot %>% 
  group_by(degree) %>% 
  summarise(p.value1 = wilcox.test(toPlot$Rsquared[toPlot$degree == 1], Rsquared)$p.value,
            p.value3 = wilcox.test(toPlot$Rsquared[toPlot$degree == 3], Rsquared)$p.value) %>%
  as.data.frame() %>% 
  melt(id.vars = "degree") %>% 
  filter(value != 1) %>%
  group_by(degree, variable) %>%
  mutate(Rsquared = jitter(0.7,0.2))

p.cv.boxplots <- ggplot(toPlot, aes(x=degree, y = Rsquared, fill = degree)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot() +
  ylab(expression(paste("Cross-validated ", R^2, sep=""))) +
  xlab("Network distance") +
  geom_text(data=tests.pvals, aes(label=format(value, scientific=T))) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid = element_blank())

toPlot.points <- toPlot %>% 
  group_by(metabolite) %>%
  filter(degree == 1) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))

toPlot <- toPlot %>% 
  group_by(metabolite) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))

p.cv <- ggplot(toPlot, aes(x = met_name, y = Rsquared, fill=degree)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_point(data = toPlot.points, aes(x = met_name, y = Rsquared)) +
  ylab(expression(paste("Cross-validated ", R^2, sep=""))) +
  xlab("") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)), limits = c(0,0.8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.2, 0.7)) +
  background_grid(major = "y", minor = "y") +
  panel_border() 
mean(toPlot$Rsquared)
table(toPlot$model)

setdiff(p.cv$data$met_name, aic.models$met_name)

toPlot$Rsquared <- round(toPlot$Rsquared,2)
toPlot$RMSE <- round(toPlot$RMSE,2)
toPlot$preprocessing <- NULL
toPlot$metabolite <- NULL

degree2radius = data.frame(degree = c(1,3,5),
                           radius = c(1,2,3))

normalization = data.frame(normalization = c("log", "log.quant", "bc"),
                           label = c("log", "log Quantile", "Box-Cox"))
toPlot$radius <- degree2radius$radius[match(toPlot$degree, degree2radius$degree)]
toPlot$degree <- NULL
toPlot$normalization <- normalization$label[match(toPlot$normalization, normalization$normalization)]

names(toPlot) <- c("Algorithm", "RMSE", "RsquaredCV", "Data transformation", "Dataset", "Metabolite", "Pathway", "Network radius")
p <- tableGrob(toPlot)
p$toScale <- T
plots.list = lappend(plots.list, p)

# ---- networks part ------

predicted.metabolites.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  group_by(metabolite, normalization, degree, preprocessing) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))


tmp = file.list[[2]] #some of sample names are not annotated, some algorithms did not keep rownames in their model matrix
tmp <- tmp[tmp$model == "glmStepAICModel",] %>% 
  dplyr::select(knockout) %>%
  mutate(n = 1:length(knockout))

predicted.metabolites.long$knockout_name <- tmp$knockout[match(predicted.metabolites.long$knockout, tmp$n)]
predicted.metabolites.long$knockout_name[which(is.na(predicted.metabolites.long$knockout_name))] <- predicted.metabolites.long$knockout[which(is.na(predicted.metabolites.long$knockout_name))]
predicted.metabolites.long$knockout <- predicted.metabolites.long$knockout_name


predicted.metabolites <- predicted.metabolites.long %>%
  dplyr::select(pred, knockout, metabolite) %>% 
  ungroup() %>%
  dcast(formula = "metabolite ~ knockout", value.var = "pred")


#dcast(predicted.metabolites, formula = "knockout ~ metabolite", value.var = "pred")
predicted.metabolite.matrix <- as.matrix(predicted.metabolites[,-1])
rownames(predicted.metabolite.matrix) <- predicted.metabolites$metabolite

bks <- seq(-4, 4, 0.25)
pheatmap(predicted.metabolite.matrix, breaks = bks, ,
         color =  colorRampPalette(c("red4","white","blue4"))(length(bks)-1), filename = "test.pdf", 
         cellwidth = 10, cellheight = 10) 

predicted.metabolites.long %>% 
  filter(!(knockout %in% c("none", "WT"))) %>%
  group_by(knockout) %>%
  summarize(min_metabolite = metabolite[which.min(pred)],
            max_metabolite = metabolite[which.max(pred)])

metabolites.cl <- hclust(as.dist((1 - cor(t(predicted.metabolite.matrix)))/2))
genes.cl <- hclust(as.dist((1 - cor(predicted.metabolite.matrix))/2))

predicted.metabolites.long$knockout = factor(predicted.metabolites.long$knockout, levels = genes.cl$labels[genes.cl$order])
predicted.metabolites.long$metabolite = factor(predicted.metabolites.long$metabolite, levels = metabolites.cl$labels[metabolites.cl$order])

cl.metabolites <- hclust(dist(predicted.metabolite.matrix))
cl.genes <- hclust(dist(t(predicted.metabolite.matrix)))
plot(cl.metabolites)

predicted.metabolites.long$metabolite = factor(predicted.metabolites.long$metabolite, levels = rownames(predicted.metabolite.matrix)[cl.metabolites$order])
predicted.metabolites.long$knockout = factor(predicted.metabolites.long$knockout, levels = colnames(predicted.metabolite.matrix)[cl.genes$order])


ggplot() +
  geom_tile(data = predicted.metabolites.long, aes(x = metabolite, y = knockout, fill=pred)) +
  scale_fill_gradient2(mid = "white", low = brewer.pal(11, "RdBu")[1], high = brewer.pal(11, "RdBu")[11],
                       midpoint = 0, limit = c(-4,4),  
                       name="predicted Z-score")



cor.matrix <- cor(t(predicted.metabolite.matrix))
cor.matrix[lower.tri(cor.matrix)] <- NA
diag(cor.matrix) <- NA


library(GGally)
library(network)
library(sna)

met2met.dataset <- melt(cor.matrix)
names(met2met.dataset) <- c("from", "to", "weight")
met2met.dataset$type <- "met->met"
cor_thr = 0.5
met.graph_dataset <- na.omit(met2met.dataset[abs(met2met.dataset$weight) > cor_thr,])

#selected_metabolites <- unique(as.vector(t(na.omit(met2met.dataset[met2met.dataset$weight > 0.25, c(1,2)]))))
predicted.metabolite.matrix.f <- predicted.metabolite.matrix[,!(colnames(predicted.metabolite.matrix) %in% c("none", "WT"))]

min.max.kinases <- melt(predicted.metabolite.matrix.f) %>% 
  group_by(X2) %>%
  summarize(min.met = X1[which.min(value)],
            z_score.min = value[which.min(value)],
            max.met = X1[which.max(value)],
            z_score.max = value[which.max(value)])

min.max.metabolites <- melt(predicted.metabolite.matrix.f) %>% 
  group_by(X1) %>%
  summarize(min.kinase = X2[which.min(value)],
            z_score.min = value[which.min(value)],
            max.kinase = X2[which.max(value)],
            z_score.max = value[which.max(value)])

kin_thr = 1.64

kinase.graph_dataset.min <- min.max.kinases %>% dplyr::select(X2, min.met, z_score.min)
names(kinase.graph_dataset.min) <- c("from", "to", "weight")
kinase.graph_dataset.min$type <- "kinase->met.min"

kinase.graph_dataset.max <- min.max.kinases %>% dplyr::select(X2, max.met, z_score.max)
names(kinase.graph_dataset.max) <- c("from", "to", "weight")
kinase.graph_dataset.max$type <- "kinase->met.max"

metabolite.graph_dataset.min <- min.max.metabolites %>% dplyr::select(min.kinase, X1, z_score.min)
names(metabolite.graph_dataset.min) <- c("from", "to", "weight")
metabolite.graph_dataset.min$type <- "kinase.min->met"

metabolite.graph_dataset.max <- min.max.metabolites %>% dplyr::select(max.kinase, X1, z_score.max)
names(metabolite.graph_dataset.max) <- c("from", "to", "weight")
metabolite.graph_dataset.max$type <- "kinase.max->met"


kinase.graph_dataset <- rbind.data.frame(kinase.graph_dataset.max, kinase.graph_dataset.min)
kinase.graph_dataset.f <- kinase.graph_dataset %>% filter(abs(weight) > kin_thr)
#kinase.graph_dataset.f <- kinase.graph_dataset

metabolite.graph_dataset <- rbind.data.frame(metabolite.graph_dataset.max, metabolite.graph_dataset.min)
metabolite.graph_dataset.f <- metabolite.graph_dataset %>% filter(abs(weight) > kin_thr) 
#metabolite.graph_dataset.f <- metabolite.graph_dataset


#storing different network pictures 
metabolites.tmp <- data.frame(id = unique(c(unique(as.vector(t(met.graph_dataset[,c(1,2)]))), as.character(kinase.graph_dataset.f$to))))
metabolites.tmp$type = "metabolite"
metabolites.tmp$node_name = sub(pattern = "-L", 
                                x = toupper(metabolite2iMM904$model_name[match(metabolites.tmp$id, metabolite2iMM904$id)]),
                                replacement = "")

kinases.tmp <- data.frame(id = unique(kinase.graph_dataset.f$from),
                          type = "kinase",
                          node_name = exp_metadata$gene[match(unique(kinase.graph_dataset.f$from), exp_metadata$ORF)])

metabolites2.tmp <- data.frame(id = unique(c(unique(as.vector(t(met.graph_dataset[,c(1,2)]))), as.character(metabolite.graph_dataset.f$to))))
metabolites2.tmp$type = "metabolite"
metabolites2.tmp$node_name = sub(pattern = "-L", 
                                 x = toupper(metabolite2iMM904$model_name[match(metabolites2.tmp$id, metabolite2iMM904$id)]),
                                 replacement = "") 

kinases2.tmp <- data.frame(id = unique(metabolite.graph_dataset.f$from),
                           type = "kinase",
                           node_name = exp_metadata$gene[match(unique(metabolite.graph_dataset.f$from), exp_metadata$ORF)])



nodes2 <- rbind.data.frame(metabolites2.tmp, kinases2.tmp)
nodes2 <- nodes2 %>% filter(!(id %in% c("G6P...F6P", "X2...3.PG")))
#write.table(x = nodes2, file = "./paper/figures/nodes2_annotation.tsv", sep = "\t", quote = F, row.names = F, col.names = F)

edges2 <- rbind.data.frame(met.graph_dataset, metabolite.graph_dataset.f)
edges2 <- edges2 %>% filter(!(from  %in% c("G6P...F6P", "X2...3.PG" )), !(to  %in% c("G6P...F6P", "X2...3.PG")))


net2 <- graph_from_data_frame(vertices = nodes2, d=edges2, directed=T) 

nodes1 <- rbind.data.frame(metabolites.tmp, kinases.tmp)
nodes1 <- nodes1 %>% filter(!(id %in% c("G6P...F6P", "X2...3.PG")))
edges1 <- rbind.data.frame(met.graph_dataset, kinase.graph_dataset.f) 
edges1 <- edges1 %>% filter(!(from  %in% c("G6P...F6P", "X2...3.PG" )), !(to  %in% c("G6P...F6P", "X2...3.PG")))

net1 <- graph_from_data_frame(vertices = nodes1, d=edges1, directed=T) 

# #net <- graph_from_data_frame( d=edges, directed=T) 
# colrs <- c( "tomato", "gold")
# 
# V(net)$color <- colrs[as.numeric(as.factor(V(net)$type))]
# V(net)$shape <- c("square", "circle")[as.numeric(as.factor(V(net)$type))]
# 
# E(net)$lsize = 1
# E(net)$lsize[E(net)$type == "met->met"] <- as.numeric(cut_number(abs(E(net)$weight[E(net)$type == "met->met"]),5))
# E(net)$lsize[E(net)$type != "met->met"] <- as.numeric(cut_number(abs(E(net)$weight[E(net)$type != "met->met"]),5))
# 
# 
# file_name = paste("predictions_regression1", fun_name, "graphml", sep=".")
# file_path = paste(output_dir, file_name, sep="/")
# write.graph(net, file=file_path, format="graphml")

colrs <- c( "tomato", "gold")

kinase.effects <- edges2 %>% 
  filter(type != "met->met") %>% 
  group_by(from) %>%
  summarise(node_score=sum(abs(weight)))

perturbation.effects <- edges1 %>% 
  filter(type != "met->met") %>% 
  group_by(to) %>% 
  summarise(node_score=sum(abs(weight)),
            degree = n())

#V(net2)$color <- colrs[as.numeric(as.factor(V(net2)$type))]
#V(net2)$shape <- c("square", "circle")[as.numeric(as.factor(V(net2)$type))]
#V(net2)$effect <- kinase.effects$node_score[match(V(net2)$name, kinase.effects$from)]
#V(net2)$effect[is.na(V(net2)$effect)] <- 0
#V(net2)$perturbation_effect <- perturbation.effects$node_score[match(V(net2)$name, perturbation.effects$to)]

# 
# E(net2)$lsize = 1
# E(net2)$lsize[E(net2)$type == "met->met"] <- as.numeric(cut(abs(E(net2)$weight[E(net2)$type == "met->met"]),
#                                                             breaks = c(0.3, 0.45, 0.75,1)))
# 
# E(net2)$lsize[E(net2)$type != "met->met"] <- as.numeric(cut(abs(E(net2)$weight[E(net2)$type != "met->met"]), 
#                                                         breaks = qnorm(c(0.5, 1-0.05/2, 1-0.01/2, 1))))
# E(net2)$regulation <- NA
# E(net2)$regulation[E(net2)$type != "met->met"] = ifelse(E(net2)$weight[E(net2)$type != "met->met"] > 0 , "pos", "neg")
# 
# E(net2)$correlation <- NA
# E(net2)$correlation[E(net2)$type == "met->met"] = ifelse(E(net2)$weight[E(net2)$type == "met->met"] > 0 , "pos", "neg")
# 
# 
# file_name = paste("predictions_regression2", fun_name, "graphml", sep=".")
# file_path = paste(output_dir, file_name, sep="/")
# write.graph(net2, file=file_path, format="graphml")

# ----- metabolite picture -----
node_properties <- data.frame(node = V(net2)$name,
                              type = V(net2)$type,
                              node_name = V(net2)$node_name)

node_properties$node_score  = kinase.effects$node_score[match(node_properties$node, kinase.effects$from)]
node_properties$node_size = as.numeric(cut(node_properties$node_score, breaks = c(0,3,10,30)))
node_properties$node_size[is.na(node_properties$node_size)] <- 1
node_properties$node_size <- 2^(node_properties$node_size)

#B <- network(get.adjacency(net2, attr="weight"))

B <- network(edges2,matrix.type='edgelist', ignore.eval = FALSE)

#node attributes
B %v% "type" <- as.character(node_properties$type[match(B %v% "vertex.names", node_properties$node)])
colors <- c("cyan", "orange")
set.vertex.attribute(B, "color", colors[as.integer(as.factor(B %v% "type"))])
set.vertex.attribute(B, "label", as.character(node_properties$node_name[match(B %v% "vertex.names", node_properties$node)]))
set.vertex.attribute(B, "perturbation", node_properties$node_score[match(B %v% "vertex.names", node_properties$node)])
set.vertex.attribute(B, "node_size", node_properties$node_size[match(B %v% "vertex.names", node_properties$node)])

#edge attributes
set.edge.attribute(B, "effect", ifelse(B %e% "weight" < 0, "red", "green"))
set.edge.attribute(B, "line_type", ifelse((as.factor(B %e% "type") == "met->met"), 2, 1))

tmp.sizes <- B %e% "weight"

tmp.sizes[B %e% "type" != "met->met"] <- as.numeric(cut(abs(tmp.sizes[B %e% "type" != "met->met"]),breaks = qnorm(c(0.5, 1-0.05, 1-0.01, 1))))
tmp.sizes[B %e% "type" == "met->met"] <- 1
set.edge.attribute(B, "edge_size", (2^tmp.sizes)/5)
set.seed(123)


p.metabolite_picture <- ggnet2(B, label.size = 3, edge.alpha = 1,
                           label = T, 
                           color = "color",
                           edge.size = "edge_size",
                           node.size = "node_size",
                           edge.color = "effect",
                           node.label = "label",
                           edge.lty ="line_type" )


## ---  kinase picture ---- 

node_properties1 <- data.frame(node = V(net1)$name,
                               type = V(net1)$type,
                               node_name = V(net1)$node_name)

node_properties1$node_score  = perturbation.effects$degree[match(node_properties1$node, perturbation.effects$to )]
node_properties1$node_size = as.numeric(cut(node_properties1$node_score,  breaks = c(0,3,10,30)))
node_properties1$node_size[is.na(node_properties1$node_size)] <- 1
node_properties1$node_size <- 2^(node_properties1$node_size)

#B <- network(get.adjacency(net2, attr="weight"))
B1 <- network(edges1,matrix.type='edgelist', ignore.eval = FALSE)

#node attributes
B1 %v% "type" <- as.character(node_properties1$type[match(B1 %v% "vertex.names", node_properties1$node)])
colors <- c("cyan", "orange")
set.vertex.attribute(B1, "color", colors[as.integer(as.factor(B1 %v% "type"))])
set.vertex.attribute(B1, "label", as.character(node_properties1$node_name[match(B1 %v% "vertex.names", node_properties1$node)]))
set.vertex.attribute(B1, "perturbation", node_properties1$node_score[match(B1 %v% "vertex.names", node_properties1$node)])
set.vertex.attribute(B1, "node_size", node_properties1$node_size[match(B1 %v% "vertex.names", node_properties1$node)])

#edge attributes
set.edge.attribute(B1, "effect", ifelse(B1 %e% "weight" < 0, "red", "green"))
set.edge.attribute(B1, "line_type", ifelse((as.factor(B1 %e% "type") == "met->met"), 2, 1))

tmp.sizes <- B1 %e% "weight"

tmp.sizes[B1 %e% "type" != "met->met"] <- as.numeric(cut(abs(tmp.sizes[B1 %e% "type" != "met->met"]),breaks = qnorm(c(0.5, 1-0.05, 1-0.01, 1))))
tmp.sizes[B1 %e% "type" == "met->met"] <- 1
set.edge.attribute(B1, "edge_size", (2^tmp.sizes)/5)
set.seed(123)

p.kinase_picture <- ggnet2(B1, label.size = 3, edge.alpha = 1, arrow.size = 1, 
                               label = T, 
                               color = "color",
                               edge.size = "edge_size",
                               node.size = "node_size",
                               edge.color = "effect",
                               node.label = "label",
                               edge.lty ="line_type" )

# ----- gene overlaps -----

load("./R/objects/overlaps.gene_overlaps.RData")
load("./R/objects/dataset.genes.gene_overlaps.RData")
overlaps.f <- overlaps %>% filter(metabolite1 %in% metabolite.order$metabolite, metabolite2 %in% metabolite.order$metabolite) %>% droplevels()

overlap.plots = list()
distance <- unique(overlaps.f$degree)

for (i in 1:length(distance)) {
  
  overlaps.f$metabolite1.lab <- toupper(metabolite2iMM904$model_name[match(overlaps.f$metabolite1, metabolite2iMM904$id)])
  overlaps.f$metabolite2.lab <- toupper(metabolite2iMM904$model_name[match(overlaps.f$metabolite2, metabolite2iMM904$id)])
  overlaps.wide1 <- dcast(overlaps.f %>% filter(degree == distance[i]), formula = "metabolite1.lab~metabolite2.lab", value.var = "overlap")
  
  overlaps.wide1.matrix <- as.matrix(overlaps.wide1[,-1])
  rownames(overlaps.wide1.matrix) <- overlaps.wide1[,1]
  
  #diag(overlaps.wide1.matrix) <- 1
  
  cl = hclust(dist(overlaps.wide1.matrix))
  my.order = colnames(overlaps.wide1.matrix)[cl$order]
  melted_cormat <- melt(overlaps.wide1.matrix)
  
  melted_cormat$X1 = factor(melted_cormat$X1,levels = my.order)
  melted_cormat$X2 = factor(melted_cormat$X2,levels = my.order)
  
  toPlot = data.frame(Var1 = melted_cormat$X1,
                      Var2 = melted_cormat$X2,
                      value = ifelse(!as.numeric(melted_cormat$X1)>=as.numeric(melted_cormat$X2),melted_cormat$value,NA),
                      label = ifelse(as.numeric(melted_cormat$X1)>=as.numeric(melted_cormat$X2),melted_cormat$value,NA))
  toPlot$label[toPlot$Var1== toPlot$Var2] = ""
  
  text = data.frame(x = my.order,
                    y = my.order,
                    my_label = c("", my.order[-length(my.order)]),
                    my_label2 = my.order, 
                    stringsAsFactors = F)
  
  p = ggplot(toPlot, aes(x = Var1, y = Var2)) +
    geom_tile(fill="white") +
    geom_point(aes(size = abs(value), color=value),pch=20) +
    scale_colour_gradient(limits=c(0, 1), low="gray90", high="black") +
    geom_text(data = text , aes(x=x, y=y, label=my_label2)) +
    xlab("") +
    ylab("") +
    theme(axis.ticks = element_blank(), 
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          legend.position = c(0.66,0.33),
          aspect.ratio = 1)
  plots.list <- lappend(plots.list, p)
  overlap.plots[[i]] <- p
}



toPlot <- bind_rows(lapply(seq_along(overlap.plots), 
                           function (x) {
                             tmp <- overlap.plots[[x]]$data
                             tmp$degree <- distance[x]
                             return(tmp)}))


p.overlaps <-ggplot(toPlot, aes(x=value, fill=degree)) +
  geom_density(alpha = 0.5) +
  xlab("Overlap of metabolite neighbours") +
  theme_bw() +
  theme(legend.position = c(0.5,0.5),
        panel.grid = element_blank())
  
  
# overlaps with total network
total <- length(unique(dataset.genes$genes))
toPlot <- dataset.genes %>% 
  group_by(metabolite, degree) %>%
  distinct(genes) %>%
  summarize(fraction = n()/total)

p.network_overlaps <- ggplot(toPlot, aes(x=fraction, fill=degree)) +
  geom_density(alpha = 0.5) +
  xlab("Metabolite neighbourhood as a fraction metabolic network") +
  theme_bw() +
  theme(legend.position = c(0.5,0.5),
        panel.grid = element_blank())

# counting measured gene neighbours
metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]
dataset.genes %>% filter(degree == 1, metabolite %in% metabolite.order$metabolite) %>% dplyr::select(genes) %>% distinct() %>% summarise(n())





# ---- predictor concentrations ------------
load("../Michael_AA/data/2016-03-30/01_Workspace.Rdata")
load("./R/objects/all_final_models.importance.models_summary2.RData")

metabolites.models.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  dplyr::select(model, RMSE, Rsquared, normalization, dataset, metabolite, degree, preprocessing) %>% 
  distinct() %>%
  group_by(metabolite, normalization, degree, preprocessing) %>%
  #group_by(metabolite) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite) %>% filter(degree <= 5) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))

metabolites.models.long <- metabolites.models.long %>% dplyr::rename(algorithm = model)

#metabolites.models.long %>% group_by(metabolite, normalization, degree, preprocessing) %>% summarise(metabolite_n = n()) %>% View
predictors.dataset <- left_join(metabolites.models.long, all_final_models.importance) 

metabolite.predictors <- predictors.dataset %>% 
  group_by(metabolite,id) %>% 
  arrange(abs(loading)) %>% 
  filter(row_number()<=10) %>% filter(Overall >= 50) %>%
  ungroup() %>%
  dplyr::select(metabolite, degree, gene ) %>% group_by(metabolite, degree, gene) %>% distinct()

load("./R/objects/overlaps.gene_overlaps.RData")
load("./R/objects/dataset.genes.gene_overlaps.RData")

metabolite.predictors$gene_name <- orf2name$gene_name[match(metabolite.predictors$gene, orf2name$ORF)]
dataset.genes <- dataset.genes %>% dplyr::rename(gene = genes)

dataset.predictors <- left_join(dataset.genes, metabolite.predictors) %>% 
  mutate(is.Predictor = ifelse(is.na(gene_name), 0, 1)) %>%
  group_by(metabolite, degree) %>%
  filter(sum(is.Predictor)>0)

data.long <- tbl_df(melt(data, id.vars = c("ORF", "gene")))
dataset.predictors <- left_join(dataset.predictors, data.long, by = c("metabolite" = "variable", "gene" = "ORF") ) %>% filter(!is.na(value))


### background only non-predictors
load("./R/objects/dataset.genes.gene_overlaps.RData")

dataset.genes <- dataset.genes %>% dplyr::rename(gene = genes)

dataset.predictors <- left_join(dataset.genes, metabolite.predictors) %>% 
  mutate(is.Predictor = ifelse(is.na(gene_name), 0, 1)) %>% filter(is.Predictor == 1) %>%
  group_by(gene) %>% mutate(n_gene = n())


bg <- left_join(dataset.genes %>% filter(!(gene %in% metabolite.predictors$gene)), metabolite.predictors) %>% 
  mutate(is.Predictor = 0)
bg <- bg %>% group_by(metabolite, gene) %>% distinct()

dataset.predictors <- bind_rows(dataset.predictors, bg)           
#dataset.predictors <- bind_rows(dataset.predictors %>% filter(n_gene <10), bg) 

data.long <- tbl_df(melt(data, id.vars = c("ORF", "gene")))
dataset.predictors <- left_join(dataset.predictors, data.long, by = c("metabolite" = "variable", "gene" = "ORF") ) %>% filter(!is.na(value))


toPlot <- dataset.predictors %>% filter(!is.na(value))

toPlot <- toPlot %>% group_by(metabolite) %>%
  mutate(valueZ = (value - mean(value, na.rm = T))/sd(value, na.rm = T))

# adjust concentration with volume and OD from this paper: doi:10.1016/j.funbio.2009.11.002
my.vol = c(median = 45.54, sd = 0.9) * 1e-15 # cell vol
my.cells = 3.2 * 10^7 * 1.5*5 # median of spectrophotometre data
#different dilution used fot AA protocol
ex.vol = 200*1e-6

toPlot$value = toPlot$value*ex.vol/(my.cells*my.vol[1])/1000 # mM



toPlot.stats <- ddply(toPlot, .(metabolite),
                      .fun = function(x) {
                        Z <<- x
                        x <- Z
                        x = x[!is.na(x$value),]
                        n_pred = sum(x$is.Predictor == 1)
                        n_bg = sum(x$is.Predictor == 0)
                        if (n_pred >2 && n_bg >2 ) {
                          #p.value = leveneTest(value~as.factor(is.Predictor), data = x)$`Pr(>F)`
                          p.value_bartlett = bartlett.test(value~as.factor(is.Predictor), data = x)$p.value
                          p.value_wilcox  = wilcox.test(x$value[x$is.Predictor == 1], x$value[x$is.Predictor == 0])$p.value
                          p.value_var = var.test(x$value[x$is.Predictor == 1], x$value[x$is.Predictor == 0], alternative = "greater")$p.value
                          return(data.frame(ymax = max(x$value), p.value_bartlett, p.value_wilcox, p.value_var))
                        }
                      })

toPlot.stats$label = ifelse(toPlot.stats$p.value_var <0.05, "*", "")
toPlot.stats$label = ifelse(toPlot.stats$p.value_var <0.01, "**", toPlot.stats$label)



p.enzyme_predictor <- ggplot(toPlot, aes(y = value, x = is.Predictor)) +
  geom_boxplot(aes(colour = factor(is.Predictor)))  +
  geom_jitter(data = toPlot, width = .3, aes(y = value, x = is.Predictor, colour = factor(is.Predictor))) +
  geom_text(data = toPlot.stats, size = 10, aes(y = ymax-0.20*ymax, x = 0.5, label = label ) ) +
  facet_wrap(~metabolite, scales = "free", nrow = 1) + 
  theme_bw() +
  xlab("") +
  ylab("Intracellular metabolite concentration, mM") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

#ranges
toPlot.summary <- toPlot %>% group_by(metabolite, is.Predictor) %>%  
  summarise(IQR = IQR(valueZ),
            rangeZ = sum(abs(range(valueZ))))

toPlot.summary.stats <- toPlot.summary %>% ungroup() %>%
  summarise(p.value = wilcox.test(rangeZ[is.Predictor == 1], rangeZ[is.Predictor == 0])$p.value)

toPlot.summary.medians <- toPlot.summary %>% group_by(is.Predictor) %>%
  summarise(medianZ = median(rangeZ, na.rm = T))

p.rangez <- ggplot(toPlot.summary, aes(x = rangeZ, fill= factor(is.Predictor))) + 
  geom_density(alpha =.5) +
  geom_vline(data = toPlot.summary.medians, linetype=2, aes(xintercept = medianZ, colour = factor(is.Predictor))) +
  annotate(geom = "text", x = 6, y = 0.3, 
           label = paste("P-value =", format(toPlot.summary.stats$p.value, digits = 2, scientific = T), sep = "")) +
  xlab("Sum of absolute extreme values of metabolite z-scores") +
  theme_bw() + 
  theme(legend.position = c(0.7, 0.2),
        aspect.ratio = 1)



#### ----- Example of few metabolites  vs WT -------

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
metabolites.long$label = metabolite2iMM904$official_name[match(metabolites.long$kegg_id, metabolite2iMM904$kegg_id )]

metabolites.long$gene_name <- exp_metadata$gene[match(metabolites.long$X1, exp_metadata$sample_name)] 
metabolites.long$ORF <- exp_metadata$ORF[match(metabolites.long$X1, exp_metadata$sample_name)] 

# metabolites.long %>% filter(gene_name %in% c("VPS15", "WT"), X2 %in% c("citrate", "histidine")) %>% View()
# metabolites.long %>% filter(gene_name %in% c("MKK1", "VPS15","AKL1", "WT"), X2 %in% c("citrate", "histidine")) %>% View()

load("./R/objects/all_final_models.models_summary2.RData")
load("./R/objects/file.list.models_summary2.RData")
load("./R/objects/proteins.matrix.quant.combat.RData")

proteins.matrix = proteins.matrix.quant.combat

met = c("aspartate", "tyrosine", "threonine", "tryptophan", "lysine", "homo.serine", "glycine")
metabolites.long.selected <- metabolites.long %>% filter(gene_name %in% c("VPS15", "WT"), X2 %in% met, dataset == "AA")


metabolite.predictors <- predictors.dataset %>% 
  group_by(metabolite,id) %>% 
  arrange(abs(loading)) %>% 
  filter(row_number()<=5) %>% filter(Overall >= 80) %>%
  ungroup() %>%
  dplyr::select(metabolite, degree, gene ) %>% 
  group_by(metabolite, degree, gene) %>% distinct() %>% 
  group_by(metabolite) %>%
  filter(row_number()<= 5) 

metabolite.predictors.selected <- metabolite.predictors %>% filter(metabolite %in% met)

#metabolite example plot for VPS15
toPlot <- metabolites.long.selected %>% 
  group_by(X2, gene_name) %>%
  summarize(mean_concentration = mean(concentration, na.rm = T),
            sd_concentration = sd(concentration, na.rm = T))

toPlot$gene_name = factor(toPlot$gene_name, levels = c("WT", "VPS15"))

p.aa_conc <- ggplot(toPlot, aes(y = mean_concentration, x = gene_name, fill = gene_name)) +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~X2, ncol = 1, scales = "free") +
  geom_errorbar(aes(ymax = mean_concentration + sd_concentration, ymin = mean_concentration - sd_concentration), 
                position = "dodge", width=0.25) +
  ylab("Intracellular metabolite concentraion, mM") +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.position = c(0.5, 0.5))

protein.matrix.long <- tbl_df(melt(proteins.matrix))
protein.matrix.long.f <- protein.matrix.long %>% 
  filter( X2 %in% as.character(unique(metabolites.long.selected$X1)), 
          X1 %in% unique(metabolite.predictors.selected$gene))

protein.matrix.long.f$gene_name <- exp_metadata$gene[match(protein.matrix.long.f$X2, exp_metadata$sample_name)]

metabolite.predictors.selected.proteins <- left_join(metabolite.predictors.selected, protein.matrix.long.f, by = c("gene" = "X1"))

toPlot = metabolite.predictors.selected.proteins %>% 
  group_by(metabolite, gene) %>%
  mutate(value_rel = exp(value)/max(exp(value))) %>%
  group_by(metabolite, gene, gene_name) %>%
  summarize(mean_value = mean(value_rel, na.rm = T),
            sd_value = sd(value_rel, na.rm = T))

toPlot$predictor_label <- orf2name$gene_name[match(toPlot$gene, orf2name$ORF)]
toPlot$gene_name = factor(toPlot$gene_name, levels = c("WT", "VPS15"))

p.aa_predictors <- ggplot(toPlot, aes(y = mean_value, x = predictor_label, fill = gene_name )) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(data = toPlot, aes(ymax = mean_value + sd_value, ymin = mean_value - sd_value), 
                position = position_dodge(width = 0.9), width = 0.5 ) +
  facet_wrap(~metabolite, ncol = 1, scales = "free") +
  scale_y_continuous(limits = c(0,1.25), breaks = c(0, 0.5, 1) ) +
  ylab("Relative protein abundance") +
  theme_bw() + 
  theme(legend.position = c(0.5,0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


# ------ predictors specificity ----------

load("./R/objects/all_final_models.importance.models_summary2.RData")
load("./R/objects/all_final_models.models_summary2.RData")
load("./R/objects/file.list.models_summary2.RData")
load("./R/objects/protens.F")
load("./R/objects/proteins.matrix.combat.quant.RData")

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

metabolites.models.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  dplyr::select(model, RMSE, Rsquared, normalization, dataset, metabolite, degree, preprocessing) %>% 
  distinct() %>%
  group_by(metabolite, normalization, degree, preprocessing) %>%
  #group_by(metabolite) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite) %>% filter(degree <= 5) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))

metabolites.models.long <- metabolites.models.long %>% rename(algorithm = model)

predictors.dataset <- left_join(metabolites.models.long, all_final_models.importance) 
# 
# metabolite.predictors <- predictors.dataset %>% filter(Overall >= 50) %>%
#   dplyr::select(metabolite, degree, gene ) %>% group_by(metabolite, degree, gene) %>% distinct()

metabolite.predictors <- predictors.dataset %>% 
  group_by(metabolite,id) %>% 
  arrange(abs(loading)) %>% 
  filter(row_number()<=10) %>% filter(Overall >= 50) %>%
  ungroup() %>%
  dplyr::select(metabolite, degree, gene ) %>% group_by(metabolite, degree, gene) %>% distinct()

# proteins.FC.f.stats <- proteins.FC.f %>% 
#   group_by(KO) %>%
#   filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr) %>% 
#   summarize(n_metabolic = sum(isMetabolic == T),
#             n_yeast  = sum(isiMM904),
#             n_total = n())
# 
# changed.genes <- proteins.FC.f %>% 
#   group_by(KO) %>% arrange() %>%
#   filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr) %>% 
#   dplyr::select(ORF, KO)
#

changed.genes <- melt(protein.matrix.mean)
 names(changed.genes) <- c("ORF", "KO", "value")
 changed.genes <- changed.genes %>% filter(KO != "none")

# changed.genes <- proteins.FC.f %>% 
#    group_by(KO) %>% arrange() %>%
#    dplyr::select(ORF, KO, logFC)


metabolite.predictors <- metabolite.predictors %>%
  group_by(metabolite, degree) %>%
  mutate(n_predictors = n())

# cosineDist <- function(x){
#   as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
# }

changed.predictors <- inner_join(metabolite.predictors, changed.genes, by = c("gene" = "ORF" ))

changed.predictors %>% group_by(metabolite, degree, n_predictors) %>% mutate(n_ko = length(unique(KO))) 

kinase.distances <- ddply(changed.predictors, .(metabolite, degree, n_predictors), 
                          .fun = function(x) {
                            #Z <<-x
                            #x <- Z
                            #tmp.df <- dcast(data = x, formula = "KO ~ gene", fun.aggregate = length)
                            
                            tmp.df <- reshape2::dcast(data = x, formula = "KO~gene", value.var = "value")
                            tmp.matrix <- as.matrix(tmp.df[,-1])
                            rownames(tmp.matrix) = tmp.df[,1]
                            d = stats::dist(scale(tmp.matrix), method = "euclidean" )
                            #d = stats::dist(tmp.matrix, method = "binary" )
                            
                            d.matrix <- as.matrix(d)
                            d.matrix[upper.tri(d.matrix)] <- NA
                            diag(d.matrix) <- NA
                            return(melt(d.matrix))
                            
                          } )

metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]

toPlot <- kinase.distances 
toPlot$metabolite.label <- toupper(metabolite2iMM904$model_name[match(toPlot$metabolite, metabolite2iMM904$id)])

toPlot.stats <- toPlot %>% 
  ungroup() %>% mutate(value100 = value/max(value, na.rm = T)) %>%
  group_by(metabolite, metabolite.label) %>%
  summarise(median = median(value100, na.rm = T)) %>% 
  ungroup() %>%
  arrange(median) %>% mutate(row_n = row_number())


toPlot$metabolite <- factor(toPlot$metabolite, levels = rev(toPlot.stats$metabolite))
toPlot$metabolite.label <- factor(toPlot$metabolite.label, levels = rev(unique(toPlot.stats$metabolite.label)))
toPlot <- toPlot %>% ungroup() %>% mutate(value100 = value/max(value, na.rm = T))

selected_rows <- seq(1, nrow(toPlot.stats), 5)
toPlot.stats.selected <- toPlot.stats %>%  filter(row_n  %in% selected_rows)


p.specificity <- ggplot(toPlot, aes(x = value100)) +
  geom_density(fill="black") +
  facet_grid(metabolite.label ~ . ) +
  xlab("Predictors response distance following kinase deletion") +
  xlim(c(0,1)) +
  scale_x_continuous(labels=percent) +
  #theme_bw() +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        strip.text.y = element_text(angle=0),
        strip.background = element_rect(colour = NULL), aspect.ratio = 1/10) 
p.specificity$toScale <- F
plots.list <- lappend(plots.list, p.specificity)

toPlot.selected <- filter(toPlot, metabolite.label %in% toPlot.stats.selected$metabolite.label)

toPlot.stats.selected$metabolite.label <- factor(toPlot.stats.selected$metabolite.label, levels = toPlot.selected$metabolite.label)
p.specificity.sel <- ggplot() +
  geom_density(data = toPlot.selected, aes(x = value100), fill="black") +
  geom_vline(data = toPlot.stats.selected, aes(xintercept = median), colour="red", linetype="longdash") +
  facet_grid(metabolite.label ~ . ) +
  xlab("Predictors response distance following kinase deletion") +
  xlim(c(0,1)) +
  scale_x_continuous(labels=percent) +
  #theme_bw() +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        strip.text.y = element_text(angle=0),
        strip.background = element_rect(colour = NULL), aspect.ratio = 1/10) 



p.specificity_violin <- ggplot(toPlot, aes(x = metabolite.label, y = value100)) +
  geom_violin(fill="black") +
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  ylab("Predictors response distance following kinase deletion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


toPlot <- predicted.metabolites.long %>% 
  group_by(metabolite) %>%
  summarise(RSD = sd(pred.untransformed, na.rm = T)/mean(pred.untransformed, na.rm = T))

toPlot$metabolite <- factor(toPlot$metabolite, levels = rev(toPlot.stats$metabolite))
toPlot$metabolite.label <- toupper(metabolite2iMM904$model_name[match(toPlot$metabolite, metabolite2iMM904$id)])
toPlot$metabolite.label <- factor(toPlot$metabolite.label, levels = rev(unique(toPlot.stats$metabolite.label)))

p.metabolites_rsd <- ggplot(toPlot, aes(x = metabolite.label, y = RSD )) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(labels=percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

toPlot$metabolite.label <- factor(toPlot$metabolite.label, levels = unique(toPlot.stats$metabolite.label))
p.metabolites_rsd.flipped <- ggplot(toPlot, aes(x = metabolite.label, y = RSD )) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels=percent) +
  coord_flip( ylim=c(0,1)) +
  theme_bw()

p.specificity_rsd <- grid.arrange(p.specificity, p.metabolites_rsd.flipped, ncol=2)
plots.list <- lappend(plots.list, p.specificity_rsd)

toPlot.selected <- filter(toPlot, metabolite.label %in% toPlot.stats.selected$metabolite.label )
p.metabolites_rsd.sel <- ggplot(toPlot.selected, aes(x = metabolite.label, y = RSD )) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels=percent) +
  coord_flip(ylim = c(0,1)) +
  theme_bw()
  

## -- toy example of saturation ----

mm_kinetics <- function(S, Km, Vmax) {
  V <- (Vmax*S)/(Km + S)
  ret <- data.frame(S, Km, Vmax, V)
  return(ret)
}

mmS_kinetics <- function(V, Km, Vmax) {
  S <- Km/(Vmax-V)
  ret <- data.frame(S, Km, Vmax, V)
  return(ret)
}


Vmax = 10
Km = 0.05
S = seq(0, 5, 0.001)

# simple MM kinetics
data.list <- lapply(seq(10), FUN = function(x) {
  mm_kinetics(S, Km, x)})
toPlot <- do.call(rbind.data.frame, data.list)

p.mm_kinetics <- ggplot(data = toPlot , aes(y = log(S/Km), x = V/Vmax)) +
  geom_point() +
  theme_bw()
plots.list <- lappend(plots.list, p.mm_kinetics)

# substare responsivness vs enzyme abundance
Vmax = seq(1, 10, 0.01)
Km = 0.05
data.list <- lapply(seq(0.001, 5,length.out = 10), FUN = function(x) {
  mmS_kinetics(x, Km, Vmax)})

toPlot <- do.call(rbind.data.frame, data.list)

toPlot <- toPlot %>% mutate(delta_S = log(S/Km) - log(lag(S/Km)),
                            delta_Vmax = Vmax - lag(Vmax),
                            delta_Vratio = V/Vmax - lag(V/Vmax) )

p.S_vs_Vmax <- ggplot(data = toPlot , aes(y = log(S/Km), x = Vmax, colour=factor(V))) +
  geom_point() + 
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = c(0.7, 0.5))
plots.list <- lappend(plots.list, p.S_vs_Vmax)


p.responsivness1 <- ggplot(data = toPlot , aes(y = delta_S/delta_Vmax, x = log(S/Km), colour=factor(V))) +
  geom_point() +
  ylim(-20, 10) + 
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = c(0.7, 0.5))
plots.list <- lappend(plots.list, p.responsivness1)

p.responsivness2 <- ggplot(data = toPlot , aes(y = delta_S/delta_Vratio, x = log(S/Km), colour=factor(V))) +
  geom_point() +
  ylim(-10, 100) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = c(0.7, 0.5))
plots.list <- lappend(plots.list, p.responsivness2)


## -- MLR models predictors plots ----


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

all_final_models.models$metabolite.label <- metabolite2iMM904$official_name[match(all_final_models.models$metabolite, metabolite2iMM904$id)]


toPlot  = dplyr::filter(all_final_models.models, type == "after", 
                        metabolite %in% metabolite.order$metabolite, 
                        adj.r.squared > 0.25, isImputed == 0 , isAutocorrelation.bonferoni == 0, the_best == T, degree == 1)
p.aic_models = ggplot() +
  geom_text(data=toPlot, aes(x=factor(model), y = varname, label=round(stats_text,2))) +
  geom_point(data=toPlot, aes(x=factor(model), y = varname, 
                              size=abs(coefficients), color = mode)) +
  facet_wrap(dataset~metabolite.label, scales="free") +
  xlab("Candidate model") + 
  ylab("Enzyme predictors") +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  scale_size_continuous(name="Effect size",
                        breaks=c(0.25, 0.5, 1),
                        labels = c("low", "medium", "strong")) +
  scale_color_discrete(name="Predictor's effect",
                       breaks = c(0, 1),
                       labels = c("negative", "positive") )


file_name = paste("supplementary_models.aic", fun_name, "pdf", sep = ".")
file_path = paste(figures_dir, file_name, sep="/")

ggsave(filename = file_path, plot = p.aic_models,  width = 210*2.1 , height = 297*2.1, units = "mm")


## -- dataset statistics ---- 

load("./R/objects/dataAA.create_datasets.RData")
load("./R/objects/dataTCA.create_datasets.RData")
load("./R/objects/dataPPP_AA.create_datasets.RData")


dataAA.long <- melt(dataAA$metabolites, varnames = "rownames")
names(dataAA.long) <- c("sample_name", "metabolite", "value")
dataAA.long$genotype <- as.character(exp_metadata$ORF[match(dataAA.long$sample_name, exp_metadata$sample_name)])
dataAA.long$dataset = "AA"

dataTCA.long <- melt(dataTCA$metabolites, varnames = "rownames")
names(dataTCA.long) <- c("sample_name", "metabolite", "value")
dataTCA.long$genotype <- as.character(exp_metadata$ORF[match(dataTCA.long$sample_name, exp_metadata$sample_name)])
dataTCA.long$dataset = "TCA"

dataPPP_AA.long <- melt(dataPPP_AA$metabolites, varnames = "rownames")
names(dataPPP_AA.long) <- c("sample_name", "metabolite", "value")
dataPPP_AA.long$genotype <- dataPPP_AA.long$sample_name
dataPPP_AA.long$dataset = "PPP_AA"


metabolites.long <- do.call(rbind.data.frame, list(dataTCA.long, dataAA.long, dataPPP_AA.long))

metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]

metabolites.long$metabolite_name <- toupper(metabolite2iMM904$model_name[match(metabolites.long$metabolite, metabolite2iMM904$id)])

toPlot <- metabolites.long %>% 
  filter(metabolite %in% metabolite.order$metabolite) %>% 
  group_by(genotype, metabolite_name, dataset) %>% summarise(n = sum(!is.na(value))) %>% 
  dcast(formula = "genotype~metabolite_name+dataset", value.var = "n") %>% 
  melt(id.vars = "genotype")

toPlot$value[toPlot$value == 0] = NA
toPlot$value[toPlot$value == 4] = 3

toPlot <- toPlot %>% group_by(genotype) %>% mutate(rowsum = sum(!is.na(value)))
toPlot <- toPlot %>% group_by(variable) %>% mutate(colsum = sum(!is.na(value)))

toPlot$genotype_name <- as.character(exp_metadata$gene[match(toPlot$genotype, exp_metadata$ORF)])
toPlot$genotype_name <- factor(toPlot$genotype_name, levels = unique(toPlot$genotype_name[order(toPlot$rowsum)]))

toPlot$genotype <- factor(toPlot$genotype, levels = unique(toPlot$genotype[order(toPlot$rowsum)]))
toPlot$variable <- factor(toPlot$variable, levels = unique(toPlot$variable[order(-toPlot$colsum)]))

toPlot <- toPlot %>% separate(variable, c("metabolite_name","dataset"), sep = "_", extra = "merge", remove = F)
toPlot$metabolite_name <- factor(toPlot$metabolite_name, levels = unique(toPlot$metabolite_name[order(-toPlot$colsum)]))

toPlot.stats <- toPlot %>% 
  group_by(variable) %>% 
  summarise(n = sum(value, na.rm = T))

p <- ggplot() +
  geom_point(data = toPlot, aes(x = variable, y = genotype_name, shape = factor(value), colour = dataset)) +
  geom_text(data = toPlot.stats, aes(y = 1, x = variable, label = n)) +
  xlab("Measured metabolite") +
  ylab("Genotype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        aspect.ratio = 1, legend.position = c(0.75,0.25)) +
  labs(shape = "Number of samples", colour = "Protocol")
p <- cowplot::ggdraw(cowplot::switch_axis_position(p , axis = 'x'))

file_name = paste("supplementary_data_description", fun_name, "pdf", sep = ".")
file_path = paste(figures_dir, file_name, sep="/")

ggsave(filename = file_path, plot = p,  width = 210*1.5 , height = 297*1.5, units = "mm")

### --- PCA of metabolite data ----

load("./R/objects/dataPPP_AA.imputed.create_datasets.RData")
load("./R/objects/dataPPP_AA.create_datasets.RData")

metabolite.data.imputed = dataPPP_AA.imputed$metabolites
metabolite.matrix = metabolite.data.imputed
metabolite.matrix = metabolite.matrix[,which(!(colnames(metabolite.matrix) %in% c("AMP", "ADP", "ATP")))]

#removing outliers
metabolite.matrix.cov = cov.rob(metabolite.matrix)
metabolite.matrix.md = mahalanobis(metabolite.matrix, center=metabolite.matrix.cov$center, cov=metabolite.matrix.cov$cov)
n = nrow(metabolite.matrix); p=ncol(metabolite.matrix)

plot(qchisq((1:n)/(n+1),p), sort(metabolite.matrix.md), 
     xlab = expression(paste(chi^2," quantiles")),
     ylab = "Sorted Machalanobis distances")
abline(0,1)

p = recordPlot()
plots.list = lappend(plots.list, p)


#decided to remove 3 points
metabolite.matrix.f = metabolite.matrix
metabolite.matrix.f = metabolite.matrix.f[!(rownames(metabolite.matrix.f) %in% names(sort(-metabolite.matrix.md)[1:3])),]

metabolite.matrix.f = metabolite.matrix.f[rownames(metabolite.matrix.f) %in% c("WT", as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),] #taking only kinase deletions

toPlot = scale(metabolite.matrix.f)
s1 = prcomp(toPlot)

# toPlot = toPlot[!rownames(toPlot) %in% names(c(which(s1$x[,1] == max(s1$x[,1])), which(s1$x[,2] == max(s1$x[,2])))),]
# s1 = prcomp(toPlot)
xPC = 1
yPC = 2
labels = orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)]
labels[is.na(labels)] = "WT"

biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)], cex=0.66,
       xlabs=labels,
       ylabs=toupper(metabolite2iMM904$model_name[match(rownames(s1$rotation[,c(xPC,yPC)]), metabolite2iMM904$id)]),
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))

abline(h=0,v=0)
p = recordPlot()
plots.list = lappend(plots.list, p)

metabolite.data = dataPPP_AA$metabolites
metabolite.data.scaled = scale(metabolite.data)

#counting how many mutants has at least one metabolite with -+2sigma perturbation
perturbation.stats <- apply(metabolite.data.scaled, 1, function(x) {
  if(any(abs(na.omit(x)) > 2)) {
    return(TRUE)
  }
  return(FALSE)
})
sum(perturbation.stats)/length(perturbation.stats)

## correlation of metabolite data with metabolic enzymes -----
load("./R/objects/dataAA.create_datasets.RData")
load("./R/objects/dataTCA.create_datasets.RData")
load("./R/objects/dataPPP_AA.create_datasets.RData")

measured.enzymes = measured.proteins[measured.proteins %in%  iMM904$gene]

corPPP_AAA.long <- melt(cor(dataPPP_AA$metabolite, dataPPP_AA$proteins[,measured.enzymes], use="pairwise.complete.obs"), id.vars = rownames )
corTCA.long <- melt(cor(dataTCA$metabolite, dataTCA$proteins[,measured.enzymes], use="pairwise.complete.obs"), id.vars = rownames)
corAA.long <- melt(cor(dataAA$metabolite, dataAA$proteins[,measured.enzymes], use="pairwise.complete.obs"), id.vars = rownames)


metabolites.cor  <- rbind(corPPP_AAA.long, corTCA.long, corAA.long)
names(metabolites.cor) = c("metabolite", "enzyme", "cor")
metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]

metabolites.cor$metabolite_name <- toupper(metabolite2iMM904$model_name[match(metabolites.cor$metabolite, metabolite2iMM904$id)])
metabolites.cor$metabolite_name <- metabolite.order$met_name[match(metabolites.cor$metabolite, metabolite.order$metabolite)]

toPlot <- metabolites.cor %>% 
  filter(metabolite %in% metabolite.order$metabolite)

toPlot$metabolite_name <- factor(toPlot$metabolite_name, levels = as.character(metabolite.order$met_name))

toPlot.stats <- toPlot %>% 
  group_by(metabolite_name) %>%
  summarise(median.cor = median(cor, na.rm=T))

p.cor <- ggplot(toPlot, aes(x = cor)) +
  geom_density() +
  geom_vline(data = toPlot.stats, aes(xintercept = median.cor), col="red", linetype = 3) +
  facet_wrap(~metabolite_name, scales="free") +
  xlab("Pearson correlation between enzyme and metabolite levels")
p.cor$toScale = T
plots.list = lappend(plots.list, p.cor)


# ---- Figure3 --------------

plot_figure3_v1 <- function() {
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(90, 60)))
  grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.energy, vp = viewport(layout.pos.row = 1:45, layout.pos.col = 1:10)) #energy metabolites R
  
  grid.text("B", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 21),gp=gpar(fontsize=20, col="black"))
  print(p.all_met, vp = viewport(layout.pos.row = 1:45, layout.pos.col = 11:30)) #all R
  
  grid.text("C", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 41),gp=gpar(fontsize=20, col="black"))
  print(p.enzymes, vp = viewport(layout.pos.row = 24:45, layout.pos.col = 31:45)) #example heatmap of glutamine
  print(p.met, vp = viewport(layout.pos.row = 24:45, layout.pos.col = 46:60)) #same
}

file_name = "Figure3_v01_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
  plot_figure3_v1()
dev.off()


file_name = "Figure3_v01_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
  plot_figure3_v1()
dev.off()


plot_figure3_v2 <- function() {
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(90, 60)))
  grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
  print(p.energy, vp = viewport(layout.pos.row = 1:45, layout.pos.col = 1:10)) #energy metabolites R
  grid.text("B", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 11),gp=gpar(fontsize=20, col="black"))
  print(p.all_met, vp = viewport(layout.pos.row = 1:45, layout.pos.col = 11:30)) #all R
  
  grid.text("C", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 31),gp=gpar(fontsize=20, col="black"))
  print(p.enzymes, vp = viewport(layout.pos.row = 1:15, layout.pos.col = 46:60)) #example heatmap of glutamine
  print(p.met.untransformed, vp = viewport(layout.pos.row = 16:30, layout.pos.col = 46:60)) #same
  print(p.gln1, vp = viewport(layout.pos.row = 16:30, layout.pos.col = 31:46))
  print(p.trna_separate, vp = viewport(layout.pos.row = 31:45, layout.pos.col = 31:46)) #example heatmap of glutamine
  print(p.aatca, vp = viewport(layout.pos.row = 31:45, layout.pos.col = 46:60)) #same
  print(p.predictors, vp = viewport(layout.pos.row = 45:70, layout.pos.col = 46:60)) #same
  
    
  
}

file_name = "Figure3_v02_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
plot_figure3_v2()
dev.off()


file_name = "Figure3_v02_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
plot_figure3_v2()
dev.off()



# ---- Figure4 --------------


plot_figure4_v1 <- function() {
  
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(130, 60)))
  grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=14, col="black"))
  print(p.cv, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 1:40)) #all metaoblites R
  
  grid.text("B", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 41),gp=gpar(fontsize=14, col="black"))
  print(p.cv.boxplots, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 41:50)) #boxplots of Rs vs degree
  
  grid.text("C", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 51),gp=gpar(fontsize=14, col="black"))
  print(p.overlaps, vp = viewport(layout.pos.row = 1:15, layout.pos.col = 51:60)) #boxplots of Rs vs degree
  
  grid.text("D", just=c("left", "centre"), vp = viewport(layout.pos.row = 16, layout.pos.col = 51),gp=gpar(fontsize=14, col="black"))
  print(p.network_overlaps, vp = viewport(layout.pos.row = 15:30, layout.pos.col = 51:60)) #boxplots of Rs vs degree
  
  
  print(p.metabolite_picture, vp = viewport(layout.pos.row = 21:60, layout.pos.col = 1:30))
  print(p.kinase_picture, vp = viewport(layout.pos.row = 21:60, layout.pos.col = 31:60))
  
  print(p.aa_conc, vp = viewport(layout.pos.row = 55:90, layout.pos.col = 1:5))
  print(p.aa_predictors, vp = viewport(layout.pos.row = 55:90, layout.pos.col = 6:15))
  print(p.enzyme_predictor, vp = viewport(layout.pos.row = 55:65, layout.pos.col = 1:60))
  print(p.rangez, vp = viewport(layout.pos.row = 66:90, layout.pos.col = 16:30))
  print(p.specificity_violin, vp = viewport(layout.pos.row = 66:90, layout.pos.col = 31:60))
  print(p.metabolites_rsd, vp = viewport(layout.pos.row = 91:110, layout.pos.col = 31:60))
  
}

file_name = "Figure4_v01_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
  plot_figure4_v1()
dev.off()


file_name = "Figure4_v01_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
  plot_figure4_v1()
dev.off()



plot_figure4_v2 <- function() {
  
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(130, 60)))
  grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=14, col="black"))
  print(p.cv, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 1:40)) #all metaoblites R
  
  grid.text("B", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 41),gp=gpar(fontsize=14, col="black"))
  print(p.cv.boxplots, vp = viewport(layout.pos.row = 1:30, layout.pos.col = 41:50)) #boxplots of Rs vs degree
  
  grid.text("C", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 51),gp=gpar(fontsize=14, col="black"))
  print(p.overlaps, vp = viewport(layout.pos.row = 1:15, layout.pos.col = 51:60)) #boxplots of Rs vs degree
  
  grid.text("D", just=c("left", "centre"), vp = viewport(layout.pos.row = 16, layout.pos.col = 51),gp=gpar(fontsize=14, col="black"))
  print(p.network_overlaps, vp = viewport(layout.pos.row = 15:30, layout.pos.col = 51:60)) #boxplots of Rs vs degree
  
  
  print(p.metabolite_picture, vp = viewport(layout.pos.row = 21:60, layout.pos.col = 1:30))
  print(p.kinase_picture, vp = viewport(layout.pos.row = 21:60, layout.pos.col = 31:60))
  
  print(p.aa_conc, vp = viewport(layout.pos.row = 55:90, layout.pos.col = 1:5))
  print(p.aa_predictors, vp = viewport(layout.pos.row = 55:90, layout.pos.col = 6:15))
  print(p.enzyme_predictor, vp = viewport(layout.pos.row = 55:65, layout.pos.col = 1:60))
  print(p.rangez, vp = viewport(layout.pos.row = 66:90, layout.pos.col = 16:30))
  #summary(p.metabolites_rsd$data$RSD)
  print(p.specificity.sel, vp = viewport(layout.pos.row = 66:90, layout.pos.col = 31:45))
  print(p.metabolites_rsd.sel, vp = viewport(layout.pos.row = 66:90, layout.pos.col = 46:60)) 
}



file_name = "Figure4_v02_scripted.png"
file_path = paste(figures_dir, file_name, sep="/")

png(file_path, height=247/25.4*2, width=183/25.4*2, units = "in", res = 150)
plot_figure4_v2()
dev.off()


file_name = "Figure4_v02_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
plot_figure4_v2()
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
           
           ggplot2::ggsave(filename = paste(file_path, x , "png", sep = "."), device = NULL,
                           plot = p, width = 210 , height = 297, dpi = 150, units = "mm", scale = scale)
           
         }, error = function(e) {
           message(paste("Plot", "x", "sucks!" ))
           return(NULL)
         }, finally = {
           message(paste("processed plot", x))
         })
       })



