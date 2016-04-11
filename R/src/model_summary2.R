rm(list = ls())
library(caret)
source("./R/boot.R")
source("./R/functions.R")

#summarizes results from regression models from 2016-02-24
# chooses best prediction models from different ML algorithm
# makes graph file for cytoscape

load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")
input_path = "./results/2016-02-24/regression_models"

load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/gene.annotations._load_.RData")
orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]

fun_name = "models_summary2"

pattern.p = "data.(\\w+).(.*?).([1,3,5]+).([0-9]+).([pca.p]+).models.RData$"
filesToProcess = dir(path=input_path, pattern = pattern.p, recursive=F)
matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)



read_models = function(x) {
  #x = matches[[1]]
  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))
  my_models$models
  t = resamples(my_models$models)
  t.long = melt(t$values[,grep("Rsquared", names(t$values))],)
  t.long.stats = t.long %>% group_by(variable) %>% summarise(meanR2=mean(value, na.rm=T)) %>% arrange(meanR2)
  t.long.stats$variable = t.long.stats$variable
  table = t.long.stats
  
  table$dataset = x[[2]]
  table$metabolite  = x[[3]]
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = x[[4]]
  table$ismetIncluded  = x[[5]]
  table$preprocessing = x[[6]]
  table$file =  x[[1]]
  return(table)
}



file.list = lapply(matches, FUN=read_models)

all_final_models = do.call(rbind.data.frame, file.list)
all_final_models$species <- all_final_models$metabolite
all_final_models$metabolite = sub(x=all_final_models$metabolite, pattern="log.quant.(.*)", replacement="\\1")
all_final_models$metabolite = sub(x=all_final_models$metabolite, pattern="log.(.*)", replacement="\\1")

all_final_models$normalization = "bc"
all_final_models$normalization[grep(pattern="log", x=all_final_models$species)] = "log"
all_final_models$normalization[grep(pattern="log.quant", x=all_final_models$species)] = "log.quant"
all_final_models$normalization = factor(all_final_models$normalization)


## screen results
annotation = data.frame(variable =                   
                          c("earthModel~Rsquared" , "rfProfile~Rsquared", "treebagModel~Rsquared",  
                            "rfModel~Rsquared",   "svmRModel~Rsquared", "gbmModel~Rsquared",
                            "lmProfile~Rsquared", "rpartModel~Rsquared", "fobaModel~Rsquared", 
                            "rlmProfile~Rsquared", "plsModel~Rsquared", "lassoModel~Rsquared", 
                            "enetModel~Rsquared",  "nnetModel~Rsquared", "glmStepAICModel~Rsquared"),
                        name = 
                          c("Multivariate Adaptive Regression Spline", "Random Forest with variable selection", "Bagged CART",
                            "Random Forest", "Support Vector Machines", "Stochastic Gradient Boosting",
                            "Linear regression with Recursive elimination", "CART", "Ridge Regression with Variable Selection",
                            "Robust regression with Recursive elimination", "Partial Least Squares", "The lasso",
                            "Elasticnet", "Model Averaged Neural Network", "Generalized Linear Model with Feature Selection"))


ddply(all_final_models, .(dataset, degree, isImputed, ismetIncluded, normalization, preprocessing), .fun=function(x) {
  
  
  model.df = dcast(formula=metabolite+file~variable, value.var="meanR2",  data=subset(x, metabolite != "ATP" & 
                                                                                         metabolite != "AMP" & 
                                                                                         metabolite != "ADP"))
  model.matrix = as.matrix(model.df[,-c(1,2)])
  rownames(model.matrix) = model.df$metabolite
  model.matrix = model.matrix[,!colnames(model.matrix) %in% c("mtModel~Rsquared", "cbModel~Rsquared","ctreeModel~Rsquared")]
  file_name = paste(fun_name,"cross-validated", "heatmap", 
                    unique(x$dataset), unique(x$isImputed), unique(x$degree), unique(x$ismetIncluded), 
                    unique(x$normalization), unique(x$preprocessing), "pdf", sep=".")
  file_path = paste(figures_dir, file_name, sep="/")
  toPlot = t(model.matrix)                       
  if (unique(x$isImputed) == 1) {
    colnames(toPlot) = sub(pattern="imputed.(*.?)", x=colnames(toPlot), replacement="\\1")
  }
    
  pheatmap(toPlot, color=brewer.pal(name="Greys", n=5), filename=file_path, width=11.7, height=8.27,
           labels_row = annotation$name[match(rownames(toPlot), annotation$variable)],
           labels_col = metabolite2iMM904$official_name[match(colnames(toPlot), metabolite2iMM904$id)],
           cellwidth=12, cellheight=12)
  return(model.df)
})
         



# reading models with predictions -------

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

load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/proteins.matrix.combat.RData")

proteins.log.quant = t(my_means(proteins.matrix.combat.quant))
proteins.raw = exp(t(my_means(proteins.matrix.combat)))
proteins.log = t(my_means(proteins.matrix.combat))


read_models.caret.predict = function(x) {
  z <<- x
  x = z
  
  file_name = paste(input_path, x[[1]], sep="/") 
  my_models = get(load(file_name))
  
  
  bad_idx <- which(names(my_models$models) == "rfModel") #removing random forest, because of predict method
  
  models = my_models$models[-bad_idx]
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
  
  
  if (!is.null(trans.x)) {
    new_data <- new_data[,names(trans.x$mean)]
    new_data <- scale(new_data, center = trans.x$mean, scale = trans.x$std) %*% trans.x$rotation #all data projected into pca space
  } else {
    new_data <- scale(new_data, center = T, scale = T)
  }
  
  models.list = list()
  
  train_idx = which(unlist(lapply(lapply(models, class), function(x) x == "train")) == TRUE)
  models = models[train_idx]
  #getting model predictions
  predictions <- predict(models, as.data.frame(new_data))
  tmp.m = matrix(unlist(predictions), ncol = length(predictions))
  colnames(tmp.m) <-  names(predictions)
  rownames(tmp.m) <-  names(predictions[[1]])
  tmp.m.long <- reshape2::melt(tmp.m)
  names(tmp.m.long) <- c("knockout", "model", "pred")
  
  #getting models stats
  t <- resamples(models)
  t.long = reshape2::melt(t$values)
  t.long.summary.wide <- t.long %>% 
    separate(variable, c("model","stats"), sep="~")  %>% 
    group_by(model, stats) %>% 
      summarise(mean_stats = mean(value, na.rm=T),
                median_stats = median(value, na.rm=T)) %>%
      dplyr::select(model, stats, mean_stats) %>%
      spread(stats, mean_stats)
  
  t.long.summary.wide$is.minRMSE <- ifelse(t.long.summary.wide$RMSE == min(t.long.summary.wide$RMSE),1,0)
  t.long.summary.wide$is.maxRsquared <- ifelse(t.long.summary.wide$Rsquared == max(t.long.summary.wide$Rsquared),1,0)
  
  table <- full_join(tmp.m.long, t.long.summary.wide, by = "model")
  
  table$dataset = x[[2]]
  table$metabolite  = x[[3]]
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = x[[4]]
  table$ismetIncluded  = x[[5]]
  table$preprocessing = x[[6]]
  table$file =  x[[1]]
  return(table)
}

file.list = lapply(matches, FUN=read_models.caret.predict)

all_final_models = do.call(rbind.data.frame, file.list)
all_final_models$species <- all_final_models$metabolite
all_final_models$metabolite = sub(x=all_final_models$metabolite, pattern="log.quant.(.*)", replacement="\\1")
all_final_models$metabolite = sub(x=all_final_models$metabolite, pattern="log.(.*)", replacement="\\1")

all_final_models$normalization = "bc"
all_final_models$normalization[grep(pattern="log", x=all_final_models$species)] = "log"
all_final_models$normalization[grep(pattern="log.quant", x=all_final_models$species)] = "log.quant"
all_final_models$normalization = factor(all_final_models$normalization)


file_name = paste("all_final_models", fun_name, "RData", sep=".")
file_path = paste(output_dir, file_name, sep="/")
save(all_final_models, file=file_path)


tmp = file.list[[2]]
tmp <- tmp[tmp$model == "glmStepAICModel",] %>% 
  dplyr::select(knockout) %>%
  mutate(n = 1:length(knockout))


predicted.metabolites.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  group_by(metabolite, normalization, degree, preprocessing) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))



predicted.metabolites.long$knockout_name <- tmp$knockout[match(predicted.metabolites.long$knockout, tmp$n)]
predicted.metabolites.long$knockout_name[which(is.na(predicted.metabolites.long$knockout_name))] <- predicted.metabolites.long$knockout[which(is.na(predicted.metabolites.long$knockout_name))]

predicted.metabolites.long$knockout <- predicted.metabolites.long$knockout_name



predicted.metabolites <- predicted.metabolites.long %>%
  dplyr::select(pred, knockout, metabolite) %>% 
  ungroup() %>%
  spread(knockout, pred)


#dcast(predicted.metabolites, formula = "knockout ~ metabolite", value.var = "pred")
predicted.metabolite.matrix <- as.matrix(predicted.metabolites[,-1])
rownames(predicted.metabolite.matrix) <- predicted.metabolites$metabolite

pheatmap(predicted.metabolite.matrix,  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 
         
predicted.metabolites.long %>% 
  filter(!(knockout %in% c("none", "WT"))) %>%
  group_by(knockout) %>%
  summarize(min_metabolite = metabolite[which.min(pred)],
            max_metabolite = metabolite[which.max(pred)])

metabolites.cl <- hclust(as.dist((1 - cor(t(predicted.metabolite.matrix)))/2))
genes.cl <- hclust(as.dist((1 - cor(predicted.metabolite.matrix))/2))

predicted.metabolites.long$knockout = factor(predicted.metabolites.long$knockout, levels = genes.cl$labels[genes.cl$order])
predicted.metabolites.long$metabolite = factor(predicted.metabolites.long$metabolite, levels = metabolites.cl$labels[metabolites.cl$order])

ggplot() +
  geom_tile(data = predicted.metabolites.long, aes(x = metabolite, y = knockout, fill=pred)) +
  scale_fill_gradient2(mid = "white", low = brewer.pal(11, "RdBu")[1], high = brewer.pal(11, "RdBu")[11],
                       midpoint = 0, limit = c(-4,5),  
                       name="predicted Z-score")

cor.matrix <- cor(t(predicted.metabolite.matrix))
cor.matrix[lower.tri(cor.matrix)] <- NA
diag(cor.matrix) <- NA



toPlot <- melt(cor.matrix)
toPlot$X1 <-  factor(toPlot$X1, levels =  metabolites.cl$labels[metabolites.cl$order])
toPlot$X2 <-  factor(toPlot$X2, levels =  rev(metabolites.cl$labels[metabolites.cl$order]))

ggplot(toPlot, aes(x=X1, y=X2, size=abs(value))) +
  geom_point()



###plotting as network -------------------

# plotting graph -----------------
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
#kinase.graph_dataset.f <- kinase.graph_dataset %>% filter(abs(weight) > kin_thr)

metabolite.graph_dataset <- rbind.data.frame(metabolite.graph_dataset.max, metabolite.graph_dataset.min)
#metabolite.graph_dataset.f <- metabolite.graph_dataset %>% filter(abs(weight) > kin_thr) 
metabolite.graph_dataset.f <- metabolite.graph_dataset

metabolites.tmp <- data.frame(id = unique(c(unique(as.vector(t(met.graph_dataset[,c(1,2)]))), as.character(kinase.graph_dataset.f$to))))
metabolites.tmp$type = "metabolite"
metabolites.tmp$node_name = metabolite2iMM904$official_name[match(metabolites.tmp$id, metabolite2iMM904$id)]

# 
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

perturbation.effects <- edges2 %>% 
  filter(type != "met->met") %>% 
  group_by(to) %>% 
  summarise(node_score=sum(abs(weight)))

V(net2)$color <- colrs[as.numeric(as.factor(V(net2)$type))]
V(net2)$shape <- c("square", "circle")[as.numeric(as.factor(V(net2)$type))]
#V(net2)$effect <- kinase.effects$node_score[match(V(net2)$name, kinase.effects$from)]
#V(net2)$effect[is.na(V(net2)$effect)] <- 0
#V(net2)$perturbation_effect <- perturbation.effects$node_score[match(V(net2)$name, perturbation.effects$to)]


E(net2)$lsize = 1
E(net2)$lsize[E(net2)$type == "met->met"] <- as.numeric(cut(abs(E(net2)$weight[E(net2)$type == "met->met"]),
                                                            breaks = c(0.3, 0.45, 0.75,1)))

E(net2)$lsize[E(net2)$type != "met->met"] <- as.numeric(cut(abs(E(net2)$weight[E(net2)$type != "met->met"]), 
                                                        breaks = qnorm(c(0.5, 1-0.05/2, 1-0.01/2, 1))))
E(net2)$regulation <- NA
E(net2)$regulation[E(net2)$type != "met->met"] = ifelse(E(net2)$weight[E(net2)$type != "met->met"] > 0 , "pos", "neg")

E(net2)$correlation <- NA
E(net2)$correlation[E(net2)$type == "met->met"] = ifelse(E(net2)$weight[E(net2)$type == "met->met"] > 0 , "pos", "neg")


file_name = paste("predictions_regression2", fun_name, "graphml", sep=".")
file_path = paste(output_dir, file_name, sep="/")
write.graph(net2, file=file_path, format="graphml")



node_properties <- data.frame(node = V(net2)$name,
                              type = V(net2)$type,
                              node_name = V(net2)$node_name)



node_properties$node_score  = kinase.effects$node_score[match(node_properties$node, kinase.effects$from)]
node_properties$node_size = as.numeric(cut(node_properties$node_score, breaks = c(0, 2, 5, 10, 30)))
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


p.kinase_picture <- ggnet2(B, label.size = 3, edge.alpha = 1,
       label = T, 
       color = "color",
       node.size = "node_size",
       edge.color = "effect",
       node.label = "label",
       edge.lty ="line_type" )
       


# load("./R/objects/KEGG.pathways.analysis1.RData")
# 
# selected = c("Amino acid metabolism",                       
#              "Carbohydrate metabolism",                                          
#              "Energy metabolism",                                      
#              "Glycan biosynthesis and metabolism",                                    
#              "Metabolism of cofactors and vitamins",       
#              "Metabolism of other amino acids",                
#              "Nucleotide metabolism",                                                         
#              "Metabolism of terpenoids and polyketides",
#              "Biosynthesis of other secondary metabolites",
#              "Lipid metabolism")
# 
# KEGG.pathways.f = droplevels(KEGG.pathways[KEGG.pathways$B %in% selected,])
# 
# my_matrix <- t(proteins.raw)
# i = 1 #index group 1
# j = 1 #index group 2
# results.tmp <- matrix(rep(0, ncol(my_matrix)^2), ncol = ncol(my_matrix))
# while (i <= ncol(my_matrix)){
#   j = i + 1
#   while (j <= ncol(my_matrix)) {
#     results.tmp[i,j] <- coef(lm(my_matrix[,i]~my_matrix[,j]))[2]
#     j = j + 1
#   }
#   i = i + 1
# }
# hist(results.tmp[upper.tri(results.tmp)])
# pheatmap(results.tmp)






