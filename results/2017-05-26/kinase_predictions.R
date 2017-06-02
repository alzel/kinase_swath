rm(list = ls())
library(tidyverse)


load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")

input_path = "./results/2017-05-26/regression_models"

load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/orf2name._clean_.RData")


# check if all models are present
load("./R/objects/all_final_models.models_summary2.RData")

metabolites.models.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  dplyr::select(file, model, RMSE, Rsquared, normalization, dataset, metabolite, degree, preprocessing) %>% 
  distinct() %>%
  group_by(file, metabolite, normalization, degree, preprocessing) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite) %>% filter(degree <= 5) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))





untransform <- function(trans, data.transformed) {
  
  stopifnot(class(trans) == "preProcess")
  ignored = c()
  if (!is.null(trans$method)) {
    methods <- unlist(as.list(my_models$trans$call$method)[-1])  
  } else {
    methods <- names(trans$method)
    ignored <- trans$method$ignore
  }
  
  new_data <- data.transformed[, !colnames(data.transformed) %in% ignored]
  
  if (any(methods == "pca")) {
    new_data = as.matrix(new_data) %*% t(trans$rotation)
  }
  
  if (any(methods == "scale")) {
    new_data <- scale(new_data, center = FALSE , scale=1/trans$std)
  }
  
  if (any(methods == "center")) {
    new_data <- scale(new_data, center = -1 * trans$mean, scale=FALSE)
  }
  
  if (any(methods == "BoxCox")) {
    new_data <- sapply(names(trans$bc), 
                       FUN = function(i) {
                         inverse.BoxCoxTrans(object = trans$bc[[i]], newdata = new_data[,i])
                       })  
  } 
  
  new_data <- as.data.frame(new_data)
  return (new_data)
}

inverse.BoxCoxTrans <- function(object, newdata) {
  if(!is.vector(newdata) || !is.numeric(newdata)) stop("newdata should be a numeric vector")
  if(is.na(object$lambda)) return(newdata) 
  
  lambda <- object$lambda
  if(object$lambda < object$fudge & object$lambda > -object$fudge)
    lambda <- 0
  else if(object$lambda < 1 + object$fudge & object$lambda > 1 - object$fudge) {
    #lambda <- 1
    warning(paste("No transformation applied, lambda within the fudge", "lambda:", lambda, "tolerance:", object$fudge))
    return(newdata)
  }
  if(lambda == 0) exp(newdata) else (lambda*newdata + 1)^(1/lambda) 
}




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


load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/proteins.matrix.combat.RData")

# proteins.log.quant = t(my_means(proteins.matrix.combat.quant))
# proteins.raw = exp(t(my_means(proteins.matrix.combat)))
# proteins.log = t(my_means(proteins.matrix.combat))


proteins.log.quant = t(proteins.matrix.combat.quant)
proteins.raw = exp(t(proteins.matrix.combat))
proteins.log = t(proteins.matrix.combat)



pattern.p = "data.(\\w+).(.*?).([1,3,5]+).([0-9]+).(\\w+).(\\w+).([pca.p]+).models.RData$"
filesToProcess = dir(path=input_path, pattern = pattern.p, recursive=F)
matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)


# tmp_df <- bind_rows(lapply(matches, FUN = function(x) {as.tibble(x)}))
# names(tmp_df) = c("file", "dataset", "species", "degree", "ismetincluded", "knockout", "model", "preprocessing")

# tmp_df <- tmp_df %>% ungroup() %>% 
#   dplyr::mutate(metabolite = stringr::str_replace(string = species, pattern="log.quant.(.*)", replacement="\\1")) %>%
#   dplyr::mutate(metabolite = stringr::str_replace(string=metabolite, pattern="log.(.*)", replacement="\\1"),
#                 normalization = "bc") 

#tmp_df$normalization[grep(pattern="log", x=tmp_df$species)] = "log"
#tmp_df$normalization[grep(pattern="log.quant", x=tmp_df$species)] = "log.quant"
#tmp_df %>% group_by(knockout) %>% summarise(n = n())



read_models.caret.predict = function(x) {
  
  # X <<- x
  # x <- X
  # x <- matches[[1]]
  
  file_name = paste(input_path, x[[1]], sep="/") 
  
  my_models = get(load(file_name))
  
  
  message(paste("processing:", file_name))
  
  input_data = my_models$input_data[,-1]
  models = my_models$models
  
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
    
    #trans.BC = preProcess(x = input_data , method = c("BoxCox"))
    
    # new_data2 = predict(trans.BC, tmp.matrix)
    #transfroming all data using BoxCox estimated on training data    
    new_data <- sapply(names(trans.x$bc), 
                       FUN = function(i) {
                         predict(object = trans.x$bc[[i]], newdata = new_data[,i])
                       }) 
  }
  
  # scaling and projecting new data to PCA space of training data  
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
  predictions <- predict(models, newdata = as.data.frame(new_data))
  
  names(predictions[[1]]) =  rownames(as.data.frame(new_data))
  
  ### untransforming predictors
  response.trans  = my_models$trans
  predictions.untransformed <- predictions
  knockout_names <- names(predictions.untransformed[[which(any(!sapply(sapply(predictions.untransformed, names), is.null)))]])

  predictions.untransformed <- lapply(predictions.untransformed,
                                      FUN = function(x) {
                                        tmp = as.vector(x)
                                        names(tmp) = knockout_names
                                        return(tmp)
                                      })

  if(!is.null(response.trans$std)) {
    predictions.untransformed <- lapply(predictions.untransformed, FUN = function(x) {
      x*response.trans$std[1]
    })
  }

  if(!is.null(response.trans$mean)) {
    predictions.untransformed <- lapply(predictions.untransformed, FUN = function(x) {
      x + response.trans$mean[1]
    })
  }

  if(!is.null(response.trans$bc)) {
    predictions.untransformed <- lapply(predictions.untransformed, FUN = function(x) {
      inverse.BoxCoxTrans(response.trans$bc[[1]], x)
    })
  }
  
  tmp.m = matrix(unlist(predictions), ncol = length(predictions))
  colnames(tmp.m) <-  names(predictions)
  rownames(tmp.m) <-  names(predictions[[1]])
  
  tmp.m2  = matrix(unlist(predictions.untransformed), ncol = length(predictions.untransformed))
  colnames(tmp.m2) <-  names(predictions.untransformed)
  rownames(tmp.m2) <-  names(predictions.untransformed[[1]])

  tmp.m.long <- reshape2::melt(tmp.m)
  tmp.m.long$knockout <- exp_metadata$ORF[match(tmp.m.long$Var1, exp_metadata$sample_name)]
  tmp.m.long <- tmp.m.long %>% dplyr::select(knockout, everything())
  
  names(tmp.m.long) <- c("knockout", "sample_name", "model", "pred")

  tmp.m.long2 <- reshape2::melt(tmp.m2)
  tmp.m.long2$knockout <- exp_metadata$ORF[match(tmp.m.long2$Var1, exp_metadata$sample_name)]
  tmp.m.long2 <- tmp.m.long2 %>% dplyr::select(knockout, everything())
  
  names(tmp.m.long2) <- c("knockout", "sample_name", "model", "pred.untransformed")
  table  <- full_join(tmp.m.long, tmp.m.long2) %>% filter(knockout == x[6])

  #table  <- tmp.m.long %>% filter(knockout == x[6])
  table$dataset = x[[2]]
  table$species  = x[[3]]
  table$isImputed = ifelse(length(grep(pattern="imputed", x=x[[3]])) == 0, 0, 1)
  table$degree  = x[[4]]
  table$ismetIncluded  = x[[5]]
  table$preprocessing = x[[8]]
  table$file =  x[[1]]
  
  return(table)
}




file.list = lapply(matches, FUN=read_models.caret.predict)

model_predictions <- bind_rows(file.list)

model_predictions <- model_predictions %>% group_by(species) %>% 
  dplyr::mutate(metabolite = stringr::str_replace(string = species, pattern="log.quant.(.*)", replacement="\\1")) %>%
  dplyr::mutate(metabolite = stringr::str_replace(string = metabolite, pattern="log.(.*)", replacement="\\1"),
         normalization = "bc") 


model_predictions$normalization[grep(pattern="log", x=model_predictions$species)] = "log"
model_predictions$normalization[grep(pattern="log.quant", x=model_predictions$species)] = "log.quant"


# metabolite data
load("./R/objects/dataAA.create_datasets.RData")
load("./R/objects/dataTCA.create_datasets.RData")
load("./R/objects/dataPPP_AA.create_datasets.RData")

dataAA.long <- reshape2::melt(dataAA$metabolites, varnames = "rownames")
names(dataAA.long) <- c("sample_name", "metabolite", "value")
dataAA.long$genotype <- as.character(exp_metadata$ORF[match(dataAA.long$sample_name, exp_metadata$sample_name)])
dataAA.long$dataset = "AA"

dataTCA.long <- reshape2::melt(dataTCA$metabolites, varnames = "rownames")
names(dataTCA.long) <- c("sample_name", "metabolite", "value")
dataTCA.long$genotype <- as.character(exp_metadata$ORF[match(dataTCA.long$sample_name, exp_metadata$sample_name)])
dataTCA.long$dataset = "TCA"

dataPPP_AA.long <- reshape2::melt(dataPPP_AA$metabolites, varnames = "rownames")
names(dataPPP_AA.long) <- c("sample_name", "metabolite", "value")
dataPPP_AA.long$genotype <- dataPPP_AA.long$sample_name
dataPPP_AA.long$dataset = "PPP_AA"

metabolites.long <- do.call(rbind.data.frame, list(dataTCA.long, dataAA.long, dataPPP_AA.long))

metabolite.order <- read.delim("./data/2015-10-16/metabolites.txt")
metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]

metabolites.long$metabolite_name <- toupper(metabolite2iMM904$model_name[match(metabolites.long$metabolite, metabolite2iMM904$id)])

metabolites.stats <- metabolites.long %>% filter(!is.na(value)) %>%
  dplyr::filter(metabolite %in% metabolite.order$metabolite) %>% ungroup() %>%
  dplyr::distinct(genotype, dataset, metabolite) %>% 
  dplyr::group_by(genotype) %>% 
  dplyr::summarise(n = n()) %>% arrange(desc(n)) %>% dplyr::mutate(rank = row_number(desc(n)))

# metabolites.sumary <- metabolites.long %>% filter(!is.na(value)) %>%
#   dplyr::group_by(dataset, metabolite) %>% 
#   dplyr::mutate(lambda = forecast::BoxCox.lambda(value)) %>% 
#   dplyr::mutate(value_BoxCox = forecast::BoxCox(value, lambda[1])) %>% 
#   dplyr::group_by(dataset, genotype, metabolite) %>% 
#   dplyr::summarise(mean_value = mean(value, na.rm = T),
#                    mean_value_BoxCox = mean(value_BoxCox, na.rm = T)) %>% 
#   dplyr::group_by(dataset, metabolite) %>% arrange(metabolite) %>% 
#   dplyr::mutate(z_mean_value = (mean_value - mean(mean_value, na.rm = T))/sd(mean_value, na.rm = T),
#                 z_mean_value_BoxCox = (mean_value_BoxCox - mean(mean_value_BoxCox, na.rm = T))/sd(mean_value_BoxCox, na.rm = T)) 
# 

metabolites.sumary <- metabolites.long %>% filter(!is.na(value)) %>%
  dplyr::group_by(dataset, metabolite) %>%
  dplyr::mutate(lambda = forecast::BoxCox.lambda(value)) %>%
  dplyr::mutate(value_BoxCox = forecast::BoxCox(value, lambda[1])) %>%
  # dplyr::group_by(dataset, genotype, metabolite) %>%
  # dplyr::summarise(mean_value = mean(value, na.rm = T),
  #                  mean_value_BoxCox = mean(value_BoxCox, na.rm = T)) %>%
  dplyr::group_by(dataset, metabolite) %>% arrange(metabolite) %>%
  # dplyr::mutate(z_mean_value = (mean_value - mean(mean_value, na.rm = T))/sd(mean_value, na.rm = T),
  #               z_mean_value_BoxCox = (mean_value_BoxCox - mean(mean_value_BoxCox, na.rm = T))/sd(mean_value_BoxCox, na.rm = T))
  dplyr::mutate(z_mean_value = (value - mean(value, na.rm = T))/sd(value, na.rm = T),
                z_mean_value_BoxCox = (value_BoxCox - mean(value_BoxCox, na.rm = T))/sd(value_BoxCox, na.rm = T))




metabolites.sumary <- metabolites.sumary %>% dplyr::rename(knockout = genotype)

load("./R/objects/metabolites.data._clean_.RData")
load("./R/objects/metabolite_metadata._clean_.RData")

metabolites.sumary_PPP_AA <- metabolites.sumary %>% filter(dataset == "PPP_AA") %>% dplyr::select(-sample_name)

model_predictions_PPP_AA <- model_predictions %>% filter(dataset == "PPP_AA") %>% 
  group_by(knockout, model, dataset, species, isImputed, degree, ismetIncluded, preprocessing, file, metabolite, normalization) %>%
  dplyr::summarise(pred = mean(pred, na.rm = T), 
                   pred.untransformed = mean(pred.untransformed, na.rm = T))
    
predicted_measured <- left_join(model_predictions %>% ungroup() %>% filter(dataset != "PPP_AA"), metabolites.sumary %>% ungroup() %>% filter(dataset != "PPP_AA"),
                                by = c("sample_name" = "sample_name", "metabolite" = "metabolite", "knockout" = "knockout", "dataset" = "dataset") ) %>% ungroup()


predicted_measured_PPP_AA <- left_join(metabolites.sumary_PPP_AA %>% ungroup(), model_predictions_PPP_AA %>% ungroup(),
          by = c("metabolite" = "metabolite", "knockout" = "knockout", "dataset" = "dataset") ) %>% ungroup()


predicted_measured <- bind_rows(predicted_measured, predicted_measured_PPP_AA)

predicted_measured <- predicted_measured %>% 
  dplyr::group_by(metabolite, dataset) %>% arrange(knockout) %>% 
  dplyr::mutate(z_pred.untransformed = (pred.untransformed - mean(pred.untransformed, na.rm = T))/sd(pred.untransformed,na.rm = T)) 


by_knockout <- predicted_measured %>% group_by(knockout) %>% 
  filter(!is.na(z_mean_value)) %>% 
  nest()

prediction_model <- function(df) {
  lm(z_pred.untransformed ~ z_mean_value, data = df, na.action = "na.omit")
}  

modify_df <- function(data, glance) {
  glance$max1 = max(data$pred, na.rm = T)
  glance$max2 = min(data$z_mean_value_BoxCox, na.rm = T) + 1 
  return(glance)
}  


model.stats <- by_knockout %>%
  dplyr::mutate(model = map(data, prediction_model)) %>% 
  dplyr::mutate(glance = map(model, broom::glance)) %>%
  dplyr::mutate(glance = map2(data, glance, modify_df)) %>%
  tidyr::unnest(glance, .drop = T) %>% arrange(desc(r.squared))


#selected_genotypes <- c("YBL088C","YBR097W","YDR122W","YDR477W","YDR490C","YGR188C","YHR079C","YHR082C","YIL042C")
selected_genotypes <- (metabolites.stats %>% filter(rank <=100, genotype != "WT"))$genotype
selected_genotypes <- (model.stats %>% filter(knockout != "WT", p.value < 0.05) %>% arrange(desc(r.squared)) %>% filter(row_number() <= 10))$knockout


toPlot <- predicted_measured %>% filter(knockout %in% selected_genotypes) %>% arrange(knockout)

toPlot.stats <- model.stats %>% filter(knockout %in% selected_genotypes) %>% arrange(desc(r.squared))
toPlot.stats$gene_label <- factor(orf2name$gene_name[match(toPlot.stats$knockout, orf2name$ORF)], levels = orf2name$gene_name[match(toPlot.stats$knockout, orf2name$ORF)])

toPlot$met_label <- metabolite2iMM904$model_name[match(toPlot$metabolite, metabolite2iMM904$id)]
toPlot$gene_label <- factor(orf2name$gene_name[match(toPlot$knockout, orf2name$ORF)], levels = toPlot.stats$gene_label)




p.pheno_geno <- ggplot(toPlot) +
  geom_point(aes(x = z_pred.untransformed, y = z_mean_value)) +
  geom_text(aes(label = paste("R2 =", format(r.squared,  digits = 2)), x = max2, y = max1), 
             data = toPlot.stats , vjust = "top", hjust = "right", parse = F) +
  geom_abline(intercept = 0, slope = 1) +
  # geom_text(data = toPlot, aes(label = met_label, x = z_pred.untransformed, y = z_mean_value),  
  #            check_overlap = T) +
  # ggrepel::geom_text_repel(data = toPlot, aes(label = met_label, x = pred, y = z_mean_value_BoxCox), 
  #                  segment.alpha = 0.25, segment.size = 0.25) +
  facet_wrap(~ gene_label , scales = "free") + 
  theme_bw() + 
  labs(x = "Predicted, standartised value",
       y = "Observed, standartised value")

p.pheno_geno


predicted_measured %>% filter(knockout == "YKL025C") %>% 
  ggplot(aes(z_mean_value, z_pred.untransformed)) +
    geom_point()
