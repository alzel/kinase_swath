#!/usr/bin/env Rscript
# linear models of metabolism for TCA metabolites from selected strains

rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis2c"


load("./R/objects/metabolite_metadata._clean_.RData")
load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/metabolites.data._clean_.RData")
load("./R/objects/protein_annotations._load_.RData")
load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/proteins.matrix.combat.RData")
load("./R/objects/exp_metadata._clean_.RData")

load("./R/objects/aa_metadata._clean_.RData")
load("./R/objects/aa.data._clean_.RData")

library(leaps)
#load("./R/objects/yeast.model._load_.RData")

library("faraway")
library("car")
library("cowplot")
library("zoo")
library("scales")
library("bootstrap")

library("MASS")
library("Amelia")


load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")


orf2name = unique(data.frame(ORF = protein_annotations$SystName, 
                             gene_name = protein_annotations$sgdName))

getResPred = function(response.matrix, predictors.matrix, i,  B, order=NULL, include.metabolites = F) {
  #B - bipartite graph with type == 1 of metabolites 
  
  stopifnot(!any(is.null(predictors.matrix)| is.null(response.matrix) | !is.bipartite(B)))
  
  if(is.null(order)) {
    order = 1  
  }
  
  current_nodes = as.character(metabolite2iMM904$model_name[metabolite2iMM904$id == i])
  stopifnot(i %in% colnames(response.matrix))
  
  if (length(current_nodes) == 0) {
    message(paste("No ", i, "found in network"))
    return(-1)
  }
  
  
  SUB = induced.subgraph(B, base::unique(unlist(neighborhood(B, order=order, nodes=current_nodes))))  
  #SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=2, nodes=current_nodes))))  
    
  genes.pred = unique(V(SUB)$name[V(SUB)$type == 0])
  genes.pred = colnames(predictors.matrix)[colnames(predictors.matrix) %in% genes.pred]
  
  if (length(genes.pred) == 0) {
    return(-1)
  }
  additional.present = integer(0)
  
  if (order > 1 & include.metabolites) {
    metabolites.pred = V(SUB)$name[V(SUB)$type == 1]
    metabolites.pred = metabolites.pred[!(metabolites.pred %in% current_nodes)] #removing response
    additional.pred = as.vector(metabolite2iMM904$id[metabolite2iMM904$model_name %in% metabolites.pred])
    additional.present = colnames(response.matrix)[which(colnames(response.matrix) %in% additional.pred)]
  }
  
  if (length(additional.present) == 0) {
    tmp.df = as.data.frame(cbind(response.matrix[,i], predictors.matrix[,genes.pred]))
    colnames(tmp.df) = c(i, genes.pred)
    return(tmp.df)
  } else {
    tmp.df = as.data.frame(cbind(response.matrix[,i], response.matrix[,additional.present], predictors.matrix[,genes.pred]))
    colnames(tmp.df) = c(i,additional.present, genes.pred)
    return(tmp.df)
  }
}


metabolite_models = function(predictors.matrix = NULL, response.matrix = NULL, models_dir = NULL, 
                             suffix = "scaled", scale.var=T, before=F, after=F,
                             order = 1, include.metabolites=F) {

   before = F
   after = F
   predictors.matrix=dataPPP_AA$proteins
   response.matrix= dataPPP_AA$metabolites
   suffix = "scaled"
   scale.var = T
   include.metabolites = F
   order = 1
   
  stopifnot(!any(is.null(predictors.matrix)| is.null(response.matrix) | is.null(models_dir) | is.null(fun_name)))
  
  yeast.model = iMM904
  yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters
  edgelist = unique(droplevels(na.omit(subset(yeast.model, metabolite != "h"  & metabolite !="h2o" , select = c("metabolite", "gene")))))
    
  B <- graph.data.frame(edgelist)
  V(B)$type <- V(B)$name %in% edgelist$metabolite
  stopifnot(is.bipartite(B))
  
  models.list = list()
  fit.list = list()
  res = data.frame()
  
  for(ji in 1:length(colnames(response.matrix))) {
    
    i <<- colnames(response.matrix)[ji]
    i = "ATP"
    #tmp.i = droplevels(subset(yeast.model.merged, id == i))
    #tmp.df = as.data.frame(cbind(response.matrix[,i], predictors.matrix[,colnames(predictors.matrix) %in% tmp.i$gene]))
    
    tmp.df = getResPred(response.matrix=response.matrix, 
                        predictors.matrix=predictors.matrix, 
                        B=B, i=i, order=order, include.metabolites=include.metabolites)
    
    if (length(tmp.df) == 1 && tmp.df  == -1) {
      next
    }
    
    if (ncol(tmp.df) == 2 ) {
      next
    }
    
    tmp.df.scaled = na.omit(as.data.frame(scale(tmp.df, center=T, scale=T)))
        
    if (scale.var == F) {
      tmp.df.scaled = na.omit(as.data.frame(scale(tmp.df, center=F, scale=F)))
      suffix = "unscaled"
    }
        
    #PCA if metabolite has more connections than samples
    if ( sum(!is.na(tmp.df.scaled[,1])) - ceiling(sum(!is.na(tmp.df.scaled[,1]))/10) < (ncol(tmp.df.scaled) - 1) ) {
      tmp.df.scaled = as.data.frame(scale(tmp.df, center=T, scale=T))
      s = prcomp(tmp.df.scaled[,-1])
      tmp.s = summary(s)
      include_components = sum(tmp.s$importance[3,] < 0.95)
      tmp.df.scaled = as.data.frame(cbind(tmp.df.scaled[,1],s$x[,1:include_components]))
      
      colnames(tmp.df.scaled)[1] = i
    }
    
    sample_size = round(nrow(tmp.df.scaled)*0.9)
    total_samples = nrow(tmp.df.scaled)
    
    k = 100
    sample.matrix = matrix(rep(0, k*sample_size), nrow=k)
    
    set.seed(123)
    for (ik in 1:k) {
      selected = sample(1:total_samples, size=sample_size)  
      sample.matrix[ik,] = selected
    }
    
    formulas = c() #storing models
    
    #choosing best model using exhaustive approach
    
    if (ncol(tmp.df.scaled)-1 < 40) {
      for(j in 1:nrow(sample.matrix)) {
        
        sample.data = tmp.df.scaled[sample.matrix[j,],]
        
        NVMAX = ncol(tmp.df.scaled) - 1
        if ( NVMAX >= ceiling(nrow(tmp.df.scaled)/2) ) {
          #NVMAX = round(nrow(sample.data)/10)*9 - 1 #to ensure 10-fold cross validation validity
          NVMAX = ceiling(nrow(tmp.df.scaled)/2)
        }
        
        b <- regsubsets(formula(paste0(i," ~ ", ".")), data=sample.data, nbest=1, nvmax=NVMAX)
        rs = summary(b)
        
        n_points = nrow(na.omit(sample.data))
        k_params = apply(rs$which, 1, sum)
        tmp.best = data.frame(n_params = apply(rs$which, 1, sum),
                              cp = rs$cp,
                              adjr2 = rs$adjr2,
                              aic = n_points*log(rs$rss/n_points) + 2*k_params,
                              bic = rs$bic)
        
        best_idx = which(tmp.best$aic == min(tmp.best$aic))
        #best_idx = which(tmp.best$adjr2 == max(tmp.best$adjr2))
        #best_idx = which(tmp.best$cp == min(tmp.best$cp))
        #best_idx = which(tmp.best$bic == min(tmp.best$bic))
        F1 = paste0(i," ~ ", paste(sort(names(which(rs$which[best_idx,-1] == TRUE))), collapse=" + "))
        formulas = c(formulas, F1)
      }
      
    } else { #using step approach
      
      for(j in 1:nrow(sample.matrix)) {
        #j = 1
        sample.data <- na.omit(tmp.df.scaled[sample.matrix[j,],])
        
        #null.lm <- lm(formula=paste0(i," ~ ", "1"), data=sample.data)
        #full.lm <- lm(formula=paste0(i," ~ ", "."), data=sample.data)
        
        null.lm <- do.call("lm", list(paste0(i," ~ ", "1"),
                                             data = sample.data))
        
        full.lm <- do.call("lm", list(paste0(i," ~ ", "."),
                                             data = sample.data))
                
        result = step(null.lm, scope=list(lower=null.lm, upper=full.lm), direction="both", steps=100000)
        tmp_vars = paste(sort(names(result$coefficients)[-1]), collapse=" + ")
        
        if ( tmp_vars == "" ) {
          tmp_vars = 1
        }
        F1 = paste0(i, " ~ ", tmp_vars)
        formulas = c(formulas, F1)
      }
    }
        
    sumarries.df = data.frame()
    fits.list = list()
    best.models = na.omit(names(sort(-table(formulas))[1:5]))
    
    for(j in 1:nrow(sample.matrix)) {
      sample.data = tmp.df.scaled[sample.matrix[j,],]
      m.list = list()
      for (m in 1:length(best.models)) {
        current.model = best.models[m]
        fit = lm(formula=as.formula(current.model), sample.data)
        fit.s = summary(fit)
        m.list[[m]] = fit
        
        tmp = data.frame(model = m, 
                         data = j,
                         formula = best.models[m],
                         coeficients = rownames(fit.s$coefficients),
                         values = fit.s$coefficients[,1],
                         r.squared = fit.s$r.squared,
                         adj.r.squared = fit.s$adj.r.squared,                 
                         p.value = ifelse(is.null(fit.s$fstatistic[1]), NA, (1 - pf(fit.s$fstatistic[1],fit.s$fstatistic[2],fit.s$fstatistic[3]))),
                         aic = AIC(fit))
        
        sumarries.df = rbind(sumarries.df, tmp)  
      }
      fits.list[[j]] = m.list
    } 
    
    sumarries.df.stats = sumarries.df %>% group_by(model,data) %>% summarize( r.squared =  r.squared[1],
                                                                              adj.r.squared = adj.r.squared[1],
                                                                              aic = aic[1],
                                                                              p.value = p.value[1])
    
    
    sumarries.df.stats = sumarries.df.stats %>% group_by(model) %>% 
      mutate(median.adj.r.squared = median(adj.r.squared, na.rm=T),
             mean.adj.r.squared = mean(adj.r.squared, na.rm=T),
             mean.r.squared = mean(r.squared, na.rm=T),
             median.r.squared = median(r.squared, na.rm=T),
             median.aic = median(aic, na.rm=T)) %>% ungroup() %>% 
      #arrange(median.aic)
      arrange(desc(median.adj.r.squared))
    
    all.data = na.omit(tmp.df.scaled)
    F1.best = as.formula(best.models[sumarries.df.stats$model[1]])
    fit = lm(formula=F1.best, data=all.data)
   
    
    outliers = outlierTest(fit)
    cooks_thr = 4/(nrow(na.omit(tmp.df.scaled)) - length(fit$coefficients)-2)
    
    clean.data = droplevels(all.data[!(rownames(all.data) %in% names(outliers$bonf.p < 0.05)),])
    #cooks_thr = 4/nrow(na.omit(tmp.df.scaled))
    
    if ( (nrow(all.data) - sum(cooks.distance(fit) > cooks_thr) - ceiling(nrow(all.data)/10)) < length(fit$coefficients)) {
      cooks_thr  = sort(cooks.distance(fit), decreasing=T)[3]
    }
    
    #all.data = droplevels(all.data[!(rownames(all.data) %in% names(which(cooks.distance(fit) == max(cooks.distance(fit))))),])
    clean.data = droplevels(clean.data[!(rownames(clean.data) %in% names(which((cooks.distance(fit) > cooks_thr) == TRUE))),])
        
    beta.changes = data.frame(((fit$coefficients + lm.influence(fit)$coef) - fit$coefficients )/fit$coefficients)
    beta.changes$sample = rownames(beta.changes)
    beta.changes.long = melt(beta.changes[,-1], id.vars="sample")
    
    n_thr = ceiling(length(unique(beta.changes.long$variable))/3)
    beta.changes.long$sample[abs(beta.changes.long$value) > 1]
    fit.after = lm(formula=F1.best, data=clean.data)
    
    #cross-valitation of all data
    theta.fit <- function(x,y){lsfit(x,y)}
    theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef} 
    
    # matrix of predictors
    X <- as.matrix(all.data[,names(fit$coefficients)[-1]])
    # vector of predicted values
    y <- all.data[,1]
    
    cv.results <- crossval(x=X,y=y,theta.fit=theta.fit,theta.predict=theta.predict, ngroup=10)
    cv.r.squared.all = cor(y,cv.results$cv.fit)**2
    
    theta.fit.after <- function(x,y){lsfit(x,y)}
    theta.predict.after <- function(fit.after,x){cbind(1,x)%*%fit.after$coef} 
    
    #cross-validation of clean data
    X <- as.matrix(clean.data[,names(fit.after$coefficients)[-1]])
    # vector of predicted values
    y <- clean.data[,1]
    cv.results.after <- crossval(x=X,y=y,theta.fit=theta.fit.after,theta.predict=theta.predict.after, ngroup=10)
    cv.r.squared.after = cor(y,cv.results.after$cv.fit)**2 
    
    if (before) {
      diag.plots1 = my_diagnostics(fit) 
      p_before = plot_grid(diag.plots1[[1]],
                           diag.plots1[[2]],
                           diag.plots1[[3]],
                           diag.plots1[[4]],
                           diag.plots1[[5]],
                           diag.plots1[[6]],
                           diag.plots1[[7]],
                           diag.plots1[[8]], labels = c(paste(i,"before")), ncol=2)  
      
      file_name = paste(fun_name,i,suffix, "diag_before.pdf", sep=".")
      file_path = paste(models_dir, file_name, sep="/")
      save_plot(filename=file_path, plot=p_before, base_height=15, base_aspect_ratio=1)
    }
    
    if (after) {
      diag.plots2 = my_diagnostics(fit.after, influence=F)
      p_after = plot_grid(diag.plots2[[1]],
                          diag.plots2[[2]],
                          diag.plots2[[3]],
                          diag.plots2[[4]],
                          diag.plots2[[5]],
                          diag.plots2[[6]],
                          diag.plots2[[7]],
                          diag.plots2[[8]], labels = c(paste(i,"after")), ncol=2)
      
      file_name = paste(fun_name,i,suffix, "diag_after.pdf", sep=".")
      file_path = paste(models_dir, file_name, sep="/")
      save_plot(filename=file_path, plot=p_after, base_height=15, base_aspect_ratio=1)
    }
    
    
    p1 = ggplot(sumarries.df.stats, aes(x=factor(model),y=adj.r.squared)) + 
                geom_boxplot()
    
    p2 = ggplot(sumarries.df, aes(x=coeficients,y=values, fill=factor(model))) +
      geom_boxplot() +
      geom_hline(yintercept = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    sum = summary(fit.after)
    p3 = tableGrob(round(sum$coefficients,5))
    model.summary = data.frame(r.sq = round(sum$r.squared,3),
                               adj.r.sq = round(sum$adj.r.squared,3),
                               fstats = round(ifelse(is.null(sum$fstatistic[1]),NA,sum$fstatistic[1]),3),
                               numdf = ifelse(is.null(sum$fstatistic[2]),NA, sum$fstatistic[2]),
                               dendf = ifelse(is.null(sum$fstatistic[3]),NA, sum$fstatistic[3]),
                               pval = ifelse(is.null(sum$fstatistic[1]),NA, (1 - pf(sum$fstatistic[1],sum$fstatistic[2],sum$fstatistic[3]))))
    
    p4 = tableGrob(model.summary)
    
    #heatmap of correlations
    pheat = cor_heatmap(clean.data=clean.data)
    
    #cross-validation
    X <- fit.after$model[,-1]
    y <- fit.after$model[,1]
    theta.fit <- function(x,y){lsfit(x,y)}
    theta.predict <- function(fit.after,x){cbind(1,x)%*%fit.after$coef}
    
    set.seed(123) 
    CVs = c()
    for (tmp.i in 1:1000) {
      cv.results <- crossval(x=X,y=y,theta.fit=theta.fit,theta.predict=theta.predict, ngroup=10)
      cv.r.squared = cor(y,cv.results$cv.fit)**2
      CVs = c(CVs, cv.r.squared)
    }
    CVs = data.frame(CVs)   

    res.tmp = data.frame(metabolite = i, subset(sumarries.df.stats, model == sumarries.df.stats$model[1]))
    tmp.model.summary = model.summary
    names(tmp.model.summary) = paste("clean", names(tmp.model.summary), sep=".")
    
    tmp.model.summary$cv.r.squared.all = cv.r.squared.all
    tmp.model.summary$cv.r.squared.after = cv.r.squared.after
    tmp.model.summary$median.cvs = median(CVs$CVs, na.rm=T)
    
    res.tmp = cbind(res.tmp, tmp.model.summary)
    res = rbind(res, res.tmp)
    
    p9 = ggplot(data=CVs, aes(x=CVs)) + geom_density()
    
    g = arrangeGrob(p1,p2,p3,p4,p9,pheat, main=textGrob(i))
       
    file_name = paste(fun_name,i, suffix, "summary.pdf", sep=".")
    file_path = paste(models_dir, file_name, sep="/")
    ggsave(filename=file_path, width=11.69 + 0.2*(11.69), height=8.27 + 0.2*(8.27), plot=g)
    break
  } 
  res$metabolite = factor(res$metabolite, 
                          levels = unique((res %>% arrange(median.r.squared))$metabolite))

  return(res)
}

make_models = function(proteins, metabolites, models_dir) {
  stopifnot(!any(is.null(proteins)| is.null(metabolites) | is.null(models_dir) ))
  
  res = metabolite_models(predictors.matrix=proteins, response.matrix=metabolites, before = F, after = T, models_dir = models_dir, scale.var=T)
  
  res$metabolite = factor(res$metabolite, 
                          levels = unique((res %>% arrange(median.r.squared))$metabolite))
  return(res)
}

plotReport = function(toPlot, fun_name, models_dir) {
  clean = toPlot %>% group_by(metabolite) %>% summarize(clean.r.sq = clean.r.sq[1],
                                                            cv.r.squared.all = cv.r.squared.all[1],
                                                            cv.r.squared.after = cv.r.squared.after[1],
                                                            median.cvs = median.cvs[1])
  
  p = ggplot(toPlot, aes(x = metabolite, y = r.squared)) + 
    geom_boxplot() + 
    geom_point(data=clean, aes(x = metabolite, y = clean.r.sq), col="red") + #R2 with cleaned data
    geom_point(data=clean, aes(x = metabolite, y = cv.r.squared.all), shape=20, col="blue") + #single CV R2 with all data
    geom_point(data=clean, aes(x = metabolite, y = cv.r.squared.after), shape=19, col="green") + #single CV R2 with cleaned data
    geom_point(data=clean, aes(x = metabolite, y = median.cvs), shape=19, col="pink") + #median of cross validataion
    coord_flip()
  
  
  file_name = paste(fun_name, "report.pdf", sep=".")
  file_path = paste(models_dir, file_name, sep="/")
  ggsave(p, filename=file_path, height=8.27+0.1*8.27, width = 11.69+0.1*11.69) 
  
}


##-- models ----

## -- selected kinases TCA ----
clean_data_TCA = function(models_dir = models_dir) {
  
  stopifnot(!is.null(models_dir))
  plots.list = list()
  load("./R/objects/metabolitesTCA_metadata._clean_.RData")
  load("./R/objects/metabolitesTCA.data._clean_.RData")

  ## -- TCA metabolites batch effects ----
  metabolitesTCA.data$date = metabolitesTCA_metadata$measure_date[match(metabolitesTCA.data$sample_id, metabolitesTCA_metadata$sample_id)]
  metabolitesTCA.data$batch = metabolitesTCA_metadata$measure_batch[match(metabolitesTCA.data$sample_id, metabolitesTCA_metadata$sample_id)]
  
  metabolitesTCA.data = droplevels(metabolitesTCA.data[grep(pattern="_[123]+$", x=metabolitesTCA.data$sample_id, perl=T),])
  metabolitesTCA.data = droplevels(filter(metabolitesTCA.data, variable != "Glu"))
  
  wt_points = metabolitesTCA.data[grep(x=metabolitesTCA.data$sample_id, pattern="WT", ignore.case=T),]
  
  library(scales)
  p1 = ggplot(metabolitesTCA.data, aes(x=batch, y=value)) +
    geom_point() +
    geom_point(data=wt_points, aes(x=batch, y=value), col="red") +
    #geom_text(data=wt_points, hjust=1, vjust = 0, aes(x=batch, y=value, label=sample_id), col="red") +
    ggtitle("Before correction") +
    facet_wrap(~variable, scales="free")  
  
  metabolitesTCA.df = dcast(metabolitesTCA.data, formula=variable~sample_id)
  metabolitesTCA.matrix = as.matrix(metabolitesTCA.df[,-1])
  rownames(metabolitesTCA.matrix) = metabolitesTCA.df$variable
  
  
  metabolitesTCA.matrix.t = t(metabolitesTCA.matrix)
    
  set.seed(123)
  metabolitesTCA.imputed = amelia(metabolitesTCA.matrix.t, logs=colnames(metabolitesTCA.matrix.t), m=5)
  metabolitesTCA.imputed.matrix = Reduce("+",metabolitesTCA.imputed$imputations)/length(metabolitesTCA.imputed$imputations)
  
    
  pheno = metabolitesTCA_metadata[metabolitesTCA_metadata$sample_id %in% rownames(metabolitesTCA.imputed.matrix),]
  pheno = as.data.frame(droplevels(pheno[match(rownames(metabolitesTCA.imputed.matrix), pheno$sample_id),]))
  pheno$measure_batch = factor(pheno$measure_batch)
  
  mod = model.matrix(~as.factor(ORF), data=pheno)
  metabolitesTCA.matrix.combat = exp(ComBat(log(t(metabolitesTCA.imputed.matrix)), batch=pheno$measure_batch, mod=mod, par.prior=T))
  metabolitesTCA.matrix.combat.long = melt(t(metabolitesTCA.matrix.combat), id.vars="rownames")
  names(metabolitesTCA.matrix.combat.long) = c("sample_id","variable", "value")
  
  metabolitesTCA.matrix.combat.long$batch = metabolitesTCA_metadata$measure_batch[match(metabolitesTCA.matrix.combat.long$sample_id, metabolitesTCA_metadata$sample_id)]
  
  wt_points.combat = metabolitesTCA.matrix.combat.long[grep(x=metabolitesTCA.matrix.combat.long$sample_id, pattern="WT", ignore.case=T),]
      
  p2 = ggplot(metabolitesTCA.matrix.combat.long, aes(x=batch, y=value)) +
    geom_point() +
    geom_point(data=wt_points.combat, aes(x=batch, y=value), col="red") +
    #geom_text(data=wt_points, hjust=1, vjust = 0, aes(x=batch, y=value, label=sample_id), col="red") +
    ggtitle("After correction") +
    facet_wrap(~variable, scales="free")
  g = arrangeGrob(p1,p2, ncol=2)
  plots.list = lappend(plots.list, g)
  
  ## -- cleaning TCA data for linear models ----
  metabolitesTCA.long = merge(metabolitesTCA.matrix.combat.long, metabolitesTCA.data, by=c("sample_id", "variable", "batch"), suffixes=c(".combat", ".raw"))
  metabolitesTCA.long$value.combat[is.na(metabolitesTCA.long$value.raw)] = NA
  
  proteins.matrix = proteins.matrix.combat.quant
    
  #normalizeQuantiles(proteins.matrix.combat)
  
  proteins.long = melt(proteins.matrix, id.vars="rownames")
  names(proteins.long) = c("ORF", "R.Label", "signal")
  proteins.long$KO = exp_metadata$ORF[match(proteins.long$R.Label, exp_metadata$sample_name)]
   
  set.seed(123)
  metabolitesTCA_metadata.wt = metabolitesTCA_metadata %>% filter(ORF== "WT") %>% group_by(ORF) %>% sample_n(3)
  metabolitesTCA_metadata.NOwt = metabolitesTCA_metadata %>% filter(ORF != "WT") %>% group_by(sample, replicate) %>% distinct()
  
  TCAmetadata = rbind(metabolitesTCA_metadata.wt, metabolitesTCA_metadata.NOwt) %>% 
    group_by(measure_batch, ORF) %>% 
    mutate(prot.sample_id = exp_metadata$sample_name[base::sample(x=which(exp_metadata$ORF %in% ORF),replace=F)[1:length(ORF)]],
           aa.sample_id = pheno$sample_id[base::sample(x=which(pheno$ORF %in% ORF),replace=F)[1:length(ORF)]])
  TCAmetadata = droplevels(TCAmetadata)
  
  metabolitesTCA.final.df = dcast(metabolitesTCA.long, formula=sample_id~variable, value.var="value.combat")
  metabolitesTCA.final.matrix = as.matrix(metabolitesTCA.final.df[,-1])
  rownames(metabolitesTCA.final.matrix) = metabolitesTCA.final.df$sample_id
  metabolitesTCA.final.matrix = metabolitesTCA.final.matrix[rownames(metabolitesTCA.final.matrix) %in% TCAmetadata$sample_id,]
  rownames(metabolitesTCA.final.matrix) = TCAmetadata$prot.sample_id[match(rownames(metabolitesTCA.final.matrix), as.vector(TCAmetadata$sample_id))]
  
  #TCA proteins
  proteins.matrix.combat.t = t(proteins.matrix)
  proteins.matrix.combat.t.f = proteins.matrix.combat.t[rownames(proteins.matrix.combat.t) %in% as.vector(TCAmetadata$prot.sample_id),]
  
  both.present = intersect(rownames(proteins.matrix.combat.t.f), rownames(metabolitesTCA.final.matrix))
  
  proteins.matrix.combat.t.f = proteins.matrix.combat.t.f[rownames(proteins.matrix.combat.t.f) %in% both.present,]
  metabolitesTCA.final.matrix = metabolitesTCA.final.matrix[rownames(metabolitesTCA.final.matrix) %in% both.present,]
  
  proteinsTCA.present = proteins.matrix.combat.t.f[match(both.present, rownames(proteins.matrix.combat.t.f)),]
  metabolitesTCA.present = metabolitesTCA.final.matrix[match(both.present, rownames(metabolitesTCA.final.matrix)),]
      
  file_name = paste(fun_name, "clean.pdf", sep=".")
  file_path = paste(models_dir, file_name, sep="/")
  save_plots(plots.list, filename=file_path, type="l") 
  
  
  return(list(proteins = proteinsTCA.present,
              metabolites = metabolitesTCA.present))
}


## -- selected TCA -- ####
## -- 1 degree selected TCA -- ####

models_dir = "./figures/models/selected_kinases_TCA"
dataTCA = clean_data_TCA(models_dir = models_dir)

unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

#resultsTCA = make_models(proteins=dataTCA$proteins, dataTCA$metabolites, models_dir=models_dir)
resultsTCA = metabolite_models(predictors.matrix=dataTCA$proteins, 
                               response.matrix=dataTCA$metabolites, 
                               before = F, after = T, models_dir = models_dir, scale.var=T, order=1, include.metabolites=F)

plotReport(toPlot=resultsTCA, fun_name=fun_name, models_dir=models_dir)


## -- selected TCA quantiles -- ####
## -- 1 degree selected TCA quantiles -- ####

models_dir = "./figures/models/selected_kinases_TCA_quantiles"
dataTCA = clean_data_TCA(models_dir = models_dir)

unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

#resultsTCA = make_models(proteins=dataTCA$proteins, dataTCA$metabolites, models_dir=models_dir)
resultsTCA = metabolite_models(predictors.matrix=dataTCA$proteins, 
                               response.matrix=dataTCA$metabolites, 
                               before = F, after = T, models_dir = models_dir, scale.var=T, order=1, include.metabolites=F)

plotReport(toPlot=resultsTCA, fun_name=fun_name, models_dir=models_dir)




## -- 2 degree selected TCA with metabolites -- ####
models_dir = "./figures/models/selected_kinases_TCA.met"
unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

resultsTCA.met = metabolite_models(predictors.matrix=dataTCA$proteins, 
                                   response.matrix=dataTCA$metabolites, 
                                   before = F, after = T, models_dir = models_dir, scale.var=T, order=2, include.metabolites=T)

plotReport(toPlot=resultsTCA.met, fun_name=fun_name, models_dir=models_dir)




clean_data_AA = function(models_dir=NULL) {
  
  stopifnot(!is.null(models_dir))

  aa.data$date = aa_metadata$date[match(aa.data$sample_id, aa_metadata$sample_id)]
  aa.data$batch = aa_metadata$Batch[match(aa.data$sample_id, aa_metadata$sample_id)]
  
  aa_metadata.f = droplevels(aa_metadata[grep(pattern="^A_|cryoStar|cryo\\+|C_WT_[ABC]+.*?", x=aa_metadata$sample_id, invert=T),])
  
  aa.data = aa.data[aa.data$sample_id %in% aa_metadata.f$sample_id,]
  qc_points = aa.data[grep(x=aa.data$sample_id, pattern="WT", ignore.case=T),]
  
    
  p1 = ggplot(aa.data, aes(x=batch, y=value)) +
    geom_point() +
    geom_point(data=qc_points, aes(x=batch, y=value), col="red") +
    geom_text(data=qc_points, aes(x=batch, y=value, label=sample_id), col="red") +
    ggtitle("Before correction")+
    #scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d")) +
    facet_wrap(~variable, scales="free")
  
  pheno = aa_metadata.f
  
  aa.data.df = dcast(aa.data, sample_id~variable, value.var="value")
  aa.data.matrix = as.matrix(aa.data.df[,-1])
  rownames(aa.data.matrix) = aa.data.df$sample_id
  
  
  #correcting for batch effects AA data
  aa.data.matrix = aa.data.matrix[complete.cases(aa.data.matrix),]
  pheno = droplevels(pheno[match(rownames(aa.data.matrix), pheno$sample_id),])
  aa.data.matrix.t = t(aa.data.matrix)
  
  mod = model.matrix(~ORF, data=pheno)
  stopifnot(length(pheno$Batch) == ncol(aa.data.matrix.t))
  
  aa.data.matrix.combat = ComBat(log(aa.data.matrix.t), batch=pheno$Batch, mod=mod, par.prior=T)
  aa.data.combat.long = melt(t(aa.data.matrix.combat), id.vars="rownames")
  names(aa.data.combat.long) = c("sample_id","variable","value" )
  
  
  aa.data.combat.long$batch = aa_metadata.f$Batch[match(aa.data.combat.long$sample_id, aa_metadata.f$sample_id)]
  
  qc_points.combat = aa.data.combat.long[grep(x=aa.data.combat.long$sample_id, pattern="WT", ignore.case=T),]
  
  p2 = ggplot(aa.data.combat.long, aes(x=batch, y=value)) +
    geom_point() +
    geom_point(data=qc_points.combat, aes(x=batch, y=value), col="red") +
    ggtitle("After correction") +
    facet_wrap(~variable, scales="free")
  
  g = arrangeGrob(p1,p2, ncol=1)
  plots.list = lappend(plots.list, g)
  
  set.seed(123)
  aa_metadata.wt = aa_metadata.f %>% filter(ORF== "WT") %>% group_by(ORF) %>% sample_n(3)
  aa_metadata.NOwt = aa_metadata.f %>% filter(ORF != "WT") %>% group_by(ORF, Replicate) %>% distinct()
  
  AAmetadata = rbind(aa_metadata.wt, aa_metadata.NOwt) %>% 
    group_by(Batch, ORF) %>% 
    mutate(prot.sample_id = exp_metadata$sample_name[base::sample(x=which(exp_metadata$ORF %in% ORF),replace=F)[1:length(ORF)]])
  
  #aa.sample_id = pheno$sample_id[base::sample(x=which(pheno$ORF %in% ORF),replace=F)[1:length(ORF)]])
    
  aa.data.combat.long$ORF = aa_metadata$ORF[match(aa.data.combat.long$sample_id, aa_metadata$sample_id)]
  aa.data.combat.final.df = dcast(aa.data.combat.long, formula = sample_id~variable)
  aa.data.combat.final.df = aa.data.combat.final.df[aa.data.combat.final.df$sample_id %in% AAmetadata$sample_id,]
  aa.data.combat.final.matrix = as.matrix(aa.data.combat.final.df[,-1])
  
  rownames(aa.data.combat.final.matrix) = AAmetadata$prot.sample_id[match(aa.data.combat.final.df[,1], AAmetadata$sample_id)]
  aa.data.combat.final.matrix = aa.data.combat.final.matrix[!is.na(rownames(aa.data.combat.final.matrix)),]
  
  ##AA proteins
  proteins.matrix = proteins.matrix.combat
  
  
  proteins.matrix.combat.t = t(proteins.matrix)
  proteins.matrix.combat.t.f = proteins.matrix.combat.t[rownames(proteins.matrix.combat.t) %in% as.vector(AAmetadata$prot.sample_id),]
  
  both.present = intersect(rownames(proteins.matrix.combat.t.f), rownames(aa.data.combat.final.matrix))
  
  proteins.matrix.combat.t.f = proteins.matrix.combat.t.f[rownames(proteins.matrix.combat.t.f) %in% both.present,]
  aa.data.combat.final.matrix = aa.data.combat.final.matrix[rownames(aa.data.combat.final.matrix) %in% both.present,]
  
  proteinsAA.present = proteins.matrix.combat.t.f[match(both.present, rownames(proteins.matrix.combat.t.f)),]
  metabolitesAA.present = aa.data.combat.final.matrix[match(both.present, rownames(aa.data.combat.final.matrix)),]  
  
  metabolitesAA.present = metabolitesAA.present[,-which(colnames(metabolitesAA.present) == "isoleucine")]
    
  file_name = paste(fun_name, "clean.pdf", sep=".")
  file_path = paste(models_dir, file_name, sep="/")
  save_plots(plots.list, filename=file_path, type="l") 
  
  return(list(proteins = proteinsAA.present,
              metabolites = metabolitesAA.present))
    
}
## -- selected AA ----



dataAA = clean_data_AA(models_dir = models_dir)


## -- 1 degree selected AA -- ####
models_dir = "./figures/models/selected_kinases_AA"
unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

resultsAA = metabolite_models(predictors.matrix=dataAA$proteins, 
                               response.matrix=dataAA$metabolites, 
                               before = F, after = T, models_dir = models_dir, scale.var=T, order=1, include.metabolites=F)

plotReport(toPlot=resultsAA, fun_name=fun_name, models_dir=models_dir)


## -- 2 degree selected AA with metabolites -- ####
models_dir = "./figures/models/selected_kinases_AA.met"
unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

resultsAA = metabolite_models(predictors.matrix=dataAA$proteins, 
                              response.matrix=dataAA$metabolites, 
                              before = F, after = T, models_dir = models_dir, scale.var=T, order=2, include.metabolites=T)

plotReport(toPlot=resultsAA, fun_name=fun_name, models_dir=models_dir)



## -- selected AA quantiles  ----
models_dir = "./figures/models/selected_kinases_AA.quantiles"
unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

dataAA = clean_data_AA(models_dir = models_dir)
resultsAA = make_models(proteins=normalizeQuantiles(dataAA$proteins), dataAA$metabolites, models_dir=models_dir)

toPlot = resultsAA
clean = resultsAA %>% group_by(metabolite) %>% summarize(clean.r.sq = clean.r.sq[1],
                                                         cv.r.squared.all = cv.r.squared.all[1],
                                                         cv.r.squared.after = cv.r.squared.after[1],
                                                         median.cvs = median.cvs[1])

p = ggplot(toPlot, aes(x = metabolite, y = r.squared)) + 
  geom_boxplot() + 
  geom_point(data=clean, aes(x = metabolite, y = clean.r.sq), col="red") + #R2 with cleaned data
  geom_point(data=clean, aes(x = metabolite, y = cv.r.squared.all), shape=20, col="blue") + #single CV R2 with all data
  geom_point(data=clean, aes(x = metabolite, y = cv.r.squared.after), shape=19, col="green") + #single CV R2 with cleaned data
  geom_point(data=clean, aes(x = metabolite, y = median.cvs), shape=19, col="pink") + #median of cross validataion
  coord_flip()

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(models_dir, file_name, sep="/")
ggsave(p, filename=file_path, height=8.27+0.1*8.27, width = 11.69+0.1*11.69) 



clean_data_PPP_AA = function(models_dir = NULL) {
  
  stopifnot(!is.null(models_dir))
  
  unlink(models_dir, recursive = T, force = FALSE)
  dir.create(models_dir)
  
  metabolites.data$batch = factor(droplevels(metabolite_metadata$batch[match(metabolites.data$sample_id, metabolite_metadata$sample_id)]))
  
  WTs = as.character(droplevels(metabolites.data$sample_id[grep(pattern="wt", x=metabolite_metadata$ORF, ignore.case=T)]))
  WTs.subset = metabolites.data[metabolites.data$sample_id %in% WTs,]
  
  p = ggplot(metabolites.data, aes(x=batch, y=value)) +
    geom_boxplot() +
    geom_point(data=WTs.subset, aes(x=batch, y=value), col="red") +
    facet_wrap(~variable, scales="free")
  plots.list = lappend(plots.list, p)
 
  
  #imputing missing values
  #TODO: check validity of imputation
  
  metabolites.df = dcast(metabolites.data, formula=variable~sample_id)
  metabolites.matrix = as.matrix(metabolites.df[,-1])
  rownames(metabolites.matrix) = metabolites.df$variable
  
  
  metabolites.matrix.t = t(metabolites.matrix)
  set.seed(123)
  metabolites.imputed = amelia(metabolites.matrix.t, logs=colnames(metabolites.matrix.t), m=5)
  metabolites.imputed.matrix = Reduce("+",metabolites.imputed$imputations)/length(metabolites.imputed$imputations)
  
  
  metabolites.imputed.long = melt(metabolites.imputed.matrix, id.vars="rownames")
  names(metabolites.imputed.long) = c("sample_id","variable","value" )
  
  metabolites.imputed.long$batch = factor(droplevels(metabolite_metadata$batch[match(metabolites.imputed.long$sample_id, metabolite_metadata$sample_id)]))
  
    
  # -- batch correction for metabolite data ----
  
  WTs.batch = unique(as.character(metabolite_metadata$batch[grep(pattern="wt", x=metabolite_metadata$ORF, ignore.case=T)]))
  WTs.imputed.subset = metabolites.imputed.long[metabolites.imputed.long$sample_id %in% WTs,]
  #WTs.batch = WTs.batch[WTs.batch != 8]
  
  
  p1 = ggplot(metabolites.imputed.long, aes(x=batch, y=value)) +
    geom_boxplot() +
    geom_point(data=WTs.imputed.subset, aes(x=batch, y=value), col="red") +
    ggtitle("Before correction") +
    facet_wrap(~variable, scales="free")
  
  
  pheno = droplevels(metabolite_metadata[metabolite_metadata$batch %in% WTs.batch,])
  metabolites.imputed.matrix.t = t(metabolites.imputed.matrix)
  metabolites.imputed.matrix.t.f = metabolites.imputed.matrix.t[,match(pheno$sample_id, colnames(metabolites.imputed.matrix.t))]
  
  mod = model.matrix(~as.factor(ORF), data=pheno)
  
  metabolites.matrix.combat = ComBat(log(metabolites.imputed.matrix.t.f), batch=pheno$batch, mod=mod, par.prior=F)
  
  metabolites.matrix.combat.long = melt(t(metabolites.matrix.combat), id.vars="rownames")
  names(metabolites.matrix.combat.long) = c("sample_id","variable", "value")
  
  metabolites.matrix.combat.long$batch = droplevels(metabolite_metadata$batch[match(metabolites.matrix.combat.long$sample_id, metabolite_metadata$sample_id)])
  WTs.combat.subset = droplevels(metabolites.matrix.combat.long[metabolites.matrix.combat.long$sample_id %in% WTs,])
  
  p2 = ggplot(metabolites.matrix.combat.long, aes(x=batch, y=value)) +
    geom_boxplot() +
    geom_point(data=WTs.combat.subset, aes(x=batch, y=value), col="red") +
    ggtitle("After correction") +
    facet_wrap(~variable, scales="free")
  
  g = arrangeGrob(p1,p2, ncol=1)
  plots.list = lappend(plots.list, g)
  
  metabolites.matrix.combat.long$ORF = pheno$ORF[match(metabolites.matrix.combat.long$sample_id, pheno$sample_id)]
  metabolites.long.mean = tbl_df(metabolites.matrix.combat.long) %>% group_by(variable, ORF) %>% summarize(mean = mean(value))
  
  metabolites.long.merged = merge(metabolites.matrix.combat.long, metabolites.data, by=c("sample_id", "variable", "batch"), suffixes=c(".combat", ".raw"))
  metabolites.long.merged$value.combat[is.na(metabolites.long.merged$value.raw)] = NA
  metabolites.long.mean.models = tbl_df(metabolites.long.merged) %>% group_by(variable, ORF) %>% summarize(mean = mean(value.combat))
  
    
  proteins.matrix = proteins.matrix.combat.quant
  
  proteins.long = melt(proteins.matrix, id.vars="rownames")
  names(proteins.long) = c("EG.StrippedSequence", "R.Label", "signal")
  proteins.long$ORF = exp_metadata$ORF[match(proteins.long$R.Label, exp_metadata$sample_name)]
  
  proteins.long.mean = tbl_df(proteins.long) %>% group_by(EG.StrippedSequence, ORF) %>% summarize(mean = mean(signal))
  
  
  proteins.mean.df = dcast(proteins.long.mean, formula=EG.StrippedSequence~ORF, value.var="mean")
  proteins.mean.matrix = as.matrix(proteins.mean.df[,-1])
  rownames(proteins.mean.matrix) = as.matrix(proteins.mean.df$EG.StrippedSequence)
  
    
  metabolites.mean.models.df = dcast(metabolites.long.mean.models, formula=variable~ORF, value.var="mean")
  metabolites.mean.matrix = as.matrix(metabolites.mean.models.df[,-1])
  rownames(metabolites.mean.matrix) = metabolites.mean.models.df$variable
  
#   ## -- amino acids from michael ----
#   
#   load("./R/objects/aa_michael.metadata._clean_.RData")
#   load("./R/objects/aa_michael.data._clean_.RData")
#   
#   aa_michael.data$date = aa_michael.metadata$date[match(aa_michael.data$sample_id, aa_michael.metadata$sample_id)]
#   aa_michael.data$batch = aa_michael.metadata$batch[match(aa_michael.data$sample_id, aa_michael.metadata$sample_id)]
#   
#   qc_points = aa_michael.data[grep(x=aa_michael.data$sample_id, pattern="QC", ignore.case=T),]
#   
#   
#   library(scales)
#   p1 = ggplot(aa_michael.data, aes(x=batch, y=value)) +
#     geom_point() +
#     geom_point(data=qc_points, aes(x=batch, y=value), col="red") +
#     ggtitle("Before correction")+
#     #scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d")) +
#     facet_wrap(~variable, scales="free")
#   
#     
#   pheno = aa_michael.metadata
#   
#   aa_michael.data.df = dcast(aa_michael.data, sample_id~variable, value.var="value")
#   aa_michael.data.matrix = as.matrix(aa_michael.data.df[,-1])
#   rownames(aa_michael.data.matrix) = aa_michael.data.df$sample_id
#   
#   
#   #correcting for batch effects Michael's data
#   aa_michael.data.matrix = aa_michael.data.matrix[complete.cases(aa_michael.data.matrix),]
#   pheno = droplevels(pheno[match(rownames(aa_michael.data.matrix), pheno$sample_id),])
#   aa_michael.data.matrix.t = t(aa_michael.data.matrix)
#   
#   mod = model.matrix(~ORF, data=pheno)
#   stopifnot(length(pheno$batch) == ncol(aa_michael.data.matrix.t))
#   
#   aa_michael.data.matrix.combat = ComBat(log(aa_michael.data.matrix.t), batch=pheno$batch, mod=mod, par.prior=T)
#   aa_michael.data.combat.long = melt(t(aa_michael.data.matrix.combat), id.vars="rownames")
#   names(aa_michael.data.combat.long) = c("sample_id","variable","value" )
#   
#   aa_michael.data.combat.long$batch = aa_michael.metadata$batch[match(aa_michael.data.combat.long$sample_id, aa_michael.metadata$sample_id)]
#   qc_points.combat = aa_michael.data.combat.long[grep(x=aa_michael.data.combat.long$sample_id, pattern="QC", ignore.case=T),]
#   
#   
#   p2 = ggplot(aa_michael.data.combat.long, aes(x=batch, y=value)) +
#     geom_point() +
#     geom_point(data=qc_points.combat, aes(x=batch, y=value), col="red") +
#     ggtitle("After correction") +
#     facet_wrap(~variable, scales="free")
#   
#   g = arrangeGrob(p1,p2, ncol=1)
#   plots.list = lappend(plots.list, g)
#   
#   aa_michael.data.combat.long$ORF = aa_michael.metadata$ORF[match(aa_michael.data.combat.long$sample_id, aa_michael.metadata$sample_id)]
#   aa_michael.combat.long.mean = tbl_df(aa_michael.data.combat.long) %>% group_by(variable, ORF) %>% summarize(mean = mean(value))
#   
#   
#   aa_michael.mean.df = dcast(aa_michael.combat.long.mean, formula=variable~ORF, value.var="mean")
#   aa_michael.mean.matrix = as.matrix(aa_michael.mean.df[,-1])
#   rownames(aa_michael.mean.matrix) = aa_michael.mean.df$variable
  
  
#   metabolites.binded = rbind(metabolites.long.mean.models, aa_michael.combat.long.mean)
#   metabolites.all.mean.df = dcast(metabolites.binded, formula=variable~ORF, value.var="mean")
  
  metabolites.binded = metabolites.long.mean.models
  metabolites.all.mean.df = dcast(metabolites.binded, formula=variable~ORF, value.var="mean")
  
  p = ggplot(metabolites.binded, aes(x=exp(mean))) +
    geom_density() +
    facet_wrap(~variable, scales="free")
  plots.list = lappend(plots.list, p)
  
  metabolites.all.mean.matrix = as.matrix(metabolites.all.mean.df[,-1])
  rownames(metabolites.all.mean.matrix) = metabolites.all.mean.df$variable
  
  both.present = intersect(colnames(proteins.mean.matrix), colnames(metabolites.all.mean.matrix))
  
  proteins.mean.matrix.present = t(proteins.mean.matrix[,match(both.present, colnames(proteins.mean.matrix))])
  metabolites.all.mean.matrix.present = t(metabolites.all.mean.matrix[,match(both.present, colnames(metabolites.all.mean.matrix))])
    
  plots.list = lappend(plots.list, p)
  
  file_name = paste(fun_name, "clean.pdf", sep=".")
  file_path = paste(models_dir, file_name, sep="/")
  save_plots(plots.list, filename=file_path, type="l")  
  
  return(list(proteins = proteins.mean.matrix.present,
              metabolites = metabolites.all.mean.matrix.present))
}

models_dir = "./figures/models/screen_kinases_PPP_AA"
unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

dataPPP_AA = clean_data_PPP_AA(models_dir = models_dir)
resultsPPP_AA = make_models(proteins=dataPPP_AA$proteins, dataPPP_AA$metabolites, models_dir=models_dir)

plotReport(toPlot=resultsPPP_AA, fun_name=fun_name, models_dir=models_dir)
str(resultsTCA)

resultsAA$dataset = "selectedAA"
resultsPPP_AA$dataset = "screenPPP_AA"
resultsTCA$dataset = "selectedTCA"

model_results = rbind(resultsAA, resultsPPP_AA, resultsTCA)
file_name = paste("model_results", "RData", sep=".")
file_path = paste(output_dir, file_name, sep="/")
save(model_results,file=file_path)






