rm(list = ls())
library(caret)
library(relaimpo)
source("./R/boot.R")
source("./R/functions.R")
library("cowplot")
library("scales")

### summarizes results from linear regression models
### makes 
### produces examples of glutamate
### Makes figures of energy metabolites


load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")
input_path = "./results/2016-02-24/linear_models"
load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/proteins.matrix.combat.RData")
load("./R/objects/exp_metadata._clean_.RData")


load("./R/objects/gene.annotations._load_.RData")
orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]


fun_name = "lm_models3"

# lm models results overview ------------------

#filesToProcess = dir(path=input_path, pattern = "[12].[01]+.linear_models.RData$", recursive=F)
filesToProcess = dir(path=input_path, pattern = "[123].[01]+.linear_models.RData$", recursive=F)
filesToProcess = grep(pattern="imputed", invert=T, filesToProcess, value=T)
filesToProcess = grep(pattern="([1,3]+).([0-9]+).linear_models.RData", filesToProcess, value=T)

pattern.p = "data.(\\w+).(.*?).([0-9]+).([0-9]+).linear_models.RData$"

matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

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
          select(-toSelect) %>%
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
  return(table)
}

file.list = lapply(matches, FUN=read_models.aic.predict)
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

all_linear_models <- all_final_models.models %>% 
  filter(type == "after",  
        isImputed == 0 , isAutocorrelation.bonferoni == 0) %>%
      group_by(metabolite, ismetIncluded, degree) %>% 
        mutate(the_super_best = adj.r.squared == max(adj.r.squared)) #the best model


metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]


toPlot = all_linear_models %>% ungroup() %>% filter(metabolite %in% metabolite.order$metabolite, degree== 1, ismetIncluded == 0, the_super_best == T, isImputed == 0 , isAutocorrelation.bonferoni == 0) %>% 
  group_by(metabolite, the_super_best) %>% 
  summarise(formula = unique(formula),
            median.cv.r2 = unique(median.cv.r2),
            adj.r.squared = unique(adj.r.squared),
            r.squared = unique(r.squared),
            model = unique(as.character(model)),
            type = unique(as.character(type)),
            datafile = unique(as.character(datafile)))

toPlot = droplevels(toPlot)
toPlot$metabolite = factor(toPlot$metabolite, levels=as.character(metabolite.order$metabolite))

toPlot$met_name = metabolite2iMM904$official_name[match(toPlot$metabolite, metabolite2iMM904$id)]
toPlot$met_name = factor(toPlot$met_name, levels=rev(as.character(metabolite.order$met_name)))
toPlot$pathway = metabolite.order$pathway[match(toPlot$metabolite, metabolite.order$metabolite)]
toPlot$pathway = factor(toPlot$pathway, levels = as.character(metabolite.order$pathway))
library(ggthemes)
p.all_met = ggplot(toPlot, aes(x = met_name, color=pathway)) +
  #geom_point(data=toPlot, aes( y = adj.r.squared), colour = "black", size=2) +
  #geom_point(data=toPlot, aes( y = median.cv.r2), colour="black", shape=17, size=5) +
  geom_linerange(data = toPlot , aes(ymin=0, ymax=adj.r.squared),size=3) +
  coord_flip() + 
  scale_colour_tableau("tableau20") +
  background_grid(major = "x", minor = "none") +
  panel_border() +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  ylab("Explained variance of metabolite concentrations\nusing proteome data, adj. R2")


# plotting graph -----------------
library(GGally)
library(network)
library(sna)

graph_dataset <- all_linear_models %>% 
  filter(ismetIncluded == 0, degree == 1, the_super_best == T, adj.r.squared > 0.3, p.value <0.01, weight > 0.2) %>%
  ungroup() %>%
  select(X1, weight, kinase.max, kinase.min, kinase.max.effect,kinase.min.effect, metabolite, adj.r.squared)

met.graph_dataset <- graph_dataset %>% select(X1, metabolite, weight)
met.graph_dataset$type <- "enz->met"
names(met.graph_dataset) <- c("from", "to", "weight", "type")

kinase.graph_dataset.min <- graph_dataset %>% select(kinase.min, X1, kinase.min.effect)
names(kinase.graph_dataset.min) <- c("from", "to", "weight")
kinase.graph_dataset.min$type <- "kinase.min->enz"

kinase.graph_dataset.max <- graph_dataset %>% select(kinase.max, X1, kinase.max.effect)
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

B = network(get.adjacency(net), directed = TRUE)

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

# glutamate example ------------------------
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
              se=T,    # Don't add shaded confidence region
              fullrange=F) +
  panel_border() +
  theme(aspect.ratio = 1)

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
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-2,2), space = "Lab", 
                       name="Protein expression") +
  ylab("Effector enzyme") +
  xlab("Kinase knockout") +
  geom_point(data=toPlot_coefs, aes(y=variable_name, x=5, size=abs(weight), colour=effector)) +
  theme(legend.position="bottom", aspect.ratio = 1) 
p.enzymes <- ggdraw(switch_axis_position(p.enzymes , axis = 'x'))
  

#Figure3 

file_name = "Figure3_v01_scripted.pdf"
file_path = paste(figures_dir, file_name, sep="/")

pdf(file_path, height=247/25.4*2, width=183/25.4*2)
grid.newpage() 
pushViewport(viewport(layout = grid.layout(90, 60)))

grid.text("A", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=20, col="black"))
print(p.all_met, vp = viewport(layout.pos.row = 1:45, layout.pos.col = 1:20)) #all metaoblites R

grid.text("B", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 21),gp=gpar(fontsize=20, col="black"))
print(p.graph, vp = viewport(layout.pos.row = 1:45, layout.pos.col = 21:30)) #graph with kinases

grid.text("C", just=c("left", "centre"), vp = viewport(layout.pos.row = 1, layout.pos.col = 41),gp=gpar(fontsize=20, col="black"))
print(p.enzymes, vp = viewport(layout.pos.row = 24:45, layout.pos.col = 31:45)) #graph with kinases
print(p.met, vp = viewport(layout.pos.row = 24:45, layout.pos.col = 46:60)) #graph with kinases

dev.off()

