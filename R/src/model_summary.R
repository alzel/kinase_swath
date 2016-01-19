library(caret)
#summarizes results from regression models

load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")
input_path = "./results/2015-08-03/models/"

fun_name = "models_summary"

filesToProcess = dir(path=input_path, pattern = ".RData$", recursive=F)
pattern.p = "data.(\\w+).(.*?).([0-9]+).([0-9]+).models.RData$"

matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

read_models = function(x) {
  file_name = paste(input_path,x[[1]], sep="/") 
  my_models = get(load(file_name))

  t = resamples(my_models)
    
  t.long = melt(t$values[,grep("Rsquared", names(t$values))],)
  t.long.stats = t.long %>% group_by(variable) %>% summarise(meanR2=mean(value, na.rm=T)) %>% arrange(meanR2)
  t.long.stats$variable = factor(t.long.stats$variable, levels=t.long.stats$variable)
  table = t.long.stats
  table$dataset = factor(x[[2]])
  table$metabolite  = factor(x[[3]])
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = factor(x[[4]])
  table$ismetIncluded  = factor(x[[5]])
  table$file =  factor(x[[1]])
  return(table)
}



file.list = lapply(matches, FUN=read_models)
all_final_models = do.call(rbind.data.frame, file.list)


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


ddply(all_final_models, .(dataset, degree, isImputed, ismetIncluded), .fun=function(x) {
  
  model.df = dcast(formula=metabolite+file~variable, value.var="meanR2",  data=subset(x, metabolite != "ATP" & 
                                                                                         metabolite != "AMP" & 
                                                                                         metabolite != "ADP"))
  model.matrix = as.matrix(model.df[,-c(1,2)])
  rownames(model.matrix) = model.df$metabolite
  model.matrix = model.matrix[,!colnames(model.matrix) %in% c("mtModel~Rsquared", "cbModel~Rsquared","ctreeModel~Rsquared")]
  file_name = paste(fun_name,"cross-validated", "heatmap", 
                    unique(x$dataset), unique(x$isImputed), unique(x$degree), unique(x$ismetIncluded), "pdf", sep=".")
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
         


