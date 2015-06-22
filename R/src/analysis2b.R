#!/usr/bin/env Rscript
# exploration of metabolite data, with batch corrections using Combat approach

rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis2a"

load("./R/objects/metabolite_metadata._clean_.RData")
load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/metabolites.data._clean_.RData")
load("./R/objects/protein_annotations._load_.RData")

orf2name = unique(data.frame(ORF = protein_annotations$SystName, 
                             gene_name = protein_annotations$sgdName))


metabolites.data$batch = factor(droplevels(metabolite_metadata$batch[match(metabolites.data$sample_id, metabolite_metadata$sample_id)]))

WTs = as.character(droplevels(metabolites.data$sample_id[grep(pattern="wt", x=metabolite_metadata$ORF, ignore.case=T)]))
WTs.subset = metabolites.data[metabolites.data$sample_id %in% WTs,]

p = ggplot(metabolites.data, aes(x=batch, y=value)) +
  geom_boxplot() +
  geom_point(data=WTs.subset, aes(x=batch, y=value), col="red") +
  facet_wrap(~variable, scales="free")
plots.list = lappend(plots.list, p)

           
library("MASS")
library("Amelia")

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


p = ggplot(metabolites.imputed.long, aes(x=batch, y=value)) +
  geom_boxplot() +
  geom_point(data=WTs.imputed.subset, aes(x=batch, y=value), col="red") +
  ggtitle("Before correction") +
  facet_wrap(~variable, scales="free")
plots.list = lappend(plots.list, p)


pheno = droplevels(metabolite_metadata[metabolite_metadata$batch %in% WTs.batch,])
metabolites.imputed.matrix.t = t(metabolites.imputed.matrix)
metabolites.imputed.matrix.t.f = metabolites.imputed.matrix.t[,match(pheno$sample_id, colnames(metabolites.imputed.matrix.t))]

mod = model.matrix(~as.factor(ORF), data=pheno)

metabolites.matrix.combat = ComBat(log(metabolites.imputed.matrix.t.f), batch=pheno$batch, mod=mod, par.prior=F)

metabolites.matrix.combat.long = melt(t(metabolites.matrix.combat), id.vars="rownames")
names(metabolites.matrix.combat.long) = c("sample_id","variable", "value")

metabolites.matrix.combat.long$batch = droplevels(metabolite_metadata$batch[match(metabolites.matrix.combat.long$sample_id, metabolite_metadata$sample_id)])
WTs.combat.subset = droplevels(metabolites.matrix.combat.long[metabolites.matrix.combat.long$sample_id %in% WTs,])

p = ggplot(metabolites.matrix.combat.long, aes(x=batch, y=value)) +
  geom_boxplot() +
  geom_point(data=WTs.combat.subset, aes(x=batch, y=value), col="red") +
  ggtitle("After correction") +
  facet_wrap(~variable, scales="free")
plots.list = lappend(plots.list, p)

metabolites.matrix.combat.long$ORF = pheno$ORF[match(metabolites.matrix.combat.long$sample_id, pheno$sample_id)]
metabolites.long.mean = tbl_df(metabolites.matrix.combat.long) %>% group_by(variable, ORF) %>% summarize(mean = mean(value))


metabolites.long.merged = merge(metabolites.matrix.combat.long, metabolites.data, by=c("sample_id", "variable", "batch"), suffixes=c(".combat", ".raw"))
metabolites.long.merged$value.combat[is.na(metabolites.long.merged$value.raw)] = NA

metabolites.long.mean.models = tbl_df(metabolites.long.merged) %>% group_by(variable, ORF) %>% summarize(mean = mean(value.combat))


load("./R/objects/proteins.matrix.combat.RData")
load("./R/objects/exp_metadata._clean_.RData")

proteins.matrix = proteins.matrix.combat

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


## -- amino acids from michael ----

load("./R/objects/aa_michael.metadata._clean_.RData")
load("./R/objects/aa_michael.data._clean_.RData")

aa_michael.data$date = aa_michael.metadata$date[match(aa_michael.data$sample_id, aa_michael.metadata$sample_id)]
aa_michael.data$batch = aa_michael.metadata$batch[match(aa_michael.data$sample_id, aa_michael.metadata$sample_id)]

qc_points = aa_michael.data[grep(x=aa_michael.data$sample_id, pattern="QC", ignore.case=T),]


library(scales)
p = ggplot(aa_michael.data, aes(x=batch, y=value)) +
  geom_point() +
  geom_point(data=qc_points, aes(x=batch, y=value), col="red") +
  ggtitle("Before correction")+
  #scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d")) +
  facet_wrap(~variable, scales="free")

plots.list = lappend(plots.list, p)

pheno = aa_michael.metadata

aa_michael.data.df = dcast(aa_michael.data, sample_id~variable, value.var="value")
aa_michael.data.matrix = as.matrix(aa_michael.data.df[,-1])
rownames(aa_michael.data.matrix) = aa_michael.data.df$sample_id


#correcting for batch effects Michael's data
aa_michael.data.matrix = aa_michael.data.matrix[complete.cases(aa_michael.data.matrix),]
pheno = droplevels(pheno[match(rownames(aa_michael.data.matrix), pheno$sample_id),])
aa_michael.data.matrix.t = t(aa_michael.data.matrix)

mod = model.matrix(~ORF, data=pheno)
stopifnot(length(pheno$batch) == ncol(aa_michael.data.matrix.t))

aa_michael.data.matrix.combat = ComBat(log(aa_michael.data.matrix.t), batch=pheno$batch, mod=mod, par.prior=T)
aa_michael.data.combat.long = melt(t(aa_michael.data.matrix.combat), id.vars="rownames")
names(aa_michael.data.combat.long) = c("sample_id","variable","value" )

aa_michael.data.combat.long$batch = aa_michael.metadata$batch[match(aa_michael.data.combat.long$sample_id, aa_michael.metadata$sample_id)]
qc_points.combat = aa_michael.data.combat.long[grep(x=aa_michael.data.combat.long$sample_id, pattern="QC", ignore.case=T),]


p = ggplot(aa_michael.data.combat.long, aes(x=batch, y=value)) +
  geom_point() +
  geom_point(data=qc_points.combat, aes(x=batch, y=value), col="red") +
  ggtitle("After correction") +
  facet_wrap(~variable, scales="free")

plots.list = lappend(plots.list, p)

aa_michael.data.combat.long$ORF = aa_michael.metadata$ORF[match(aa_michael.data.combat.long$sample_id, aa_michael.metadata$sample_id)]
aa_michael.combat.long.mean = tbl_df(aa_michael.data.combat.long) %>% group_by(variable, ORF) %>% summarize(mean = mean(value))


aa_michael.mean.df = dcast(aa_michael.combat.long.mean, formula=variable~ORF, value.var="mean")
aa_michael.mean.matrix = as.matrix(aa_michael.mean.df[,-1])
rownames(aa_michael.mean.matrix) = aa_michael.mean.df$variable


## -- making linear model of metabolism ----
library(leaps)
load("./R/objects/yeast.model._load_.RData")
load("./R/objects/iMM904._load_.RData")

metabolites.binded = rbind(metabolites.long.mean.models, aa_michael.combat.long.mean)
metabolites.all.mean.df = dcast(metabolites.binded, formula=variable~ORF, value.var="mean")

p = ggplot(metabolites.binded, aes(x=exp(mean))) +
  geom_density() +
  facet_wrap(~variable, scales="free")
plots.list = lappend(plots.list, p)


metabolites.all.mean.matrix = as.matrix(metabolites.all.mean.df[,-1])
rownames(metabolites.all.mean.matrix) = metabolites.all.mean.df$variable

# metabolites.all.mean.matrix.t = t(metabolites.all.mean.matrix)
# z.matrix = scale(exp(metabolites.all.mean.matrix.t), center=T, scale=T)
# metabolites.all.mean.matrix.t[abs(z.matrix) > 4] = NA
# metabolites.all.mean.matrix = t(metabolites.all.mean.matrix.t)



both.present = intersect(colnames(proteins.mean.matrix), colnames(metabolites.all.mean.matrix))

proteins.mean.matrix.present = proteins.mean.matrix[,match(both.present, colnames(proteins.mean.matrix))]
metabolites.all.mean.matrix.present = metabolites.all.mean.matrix[,match(both.present, colnames(metabolites.all.mean.matrix))]

load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")


yeast.model = iMM904
head(yeast.model)
head(metabolite2iMM904)

yeast.model.merged = droplevels(merge(yeast.model, metabolite2iMM904, by.x="metabolite", by.y="model_name"))

str(yeast.model.merged)

table(yeast.model.merged$metabolite)

library("faraway")
library("car")
library("cowplot")

models.list = list()
fit.list = list()

yeast.model.merged = yeast.model.merged[yeast.model.merged$gene %in% rownames(proteins.mean.matrix.present),]
yeast.model.merged$gene_name = orf2name$gene_name[match(yeast.model.merged$gene, orf2name$ORF)]

models_dir = "./figures/models/screen"
unlink(models_dir, recursive = T, force = FALSE)
dir.create(models_dir)

for(j in 1:length(levels(yeast.model.merged$id))) {
    
  i = levels(yeast.model.merged$id)[j]
  
  tmp.i = droplevels(subset(yeast.model.merged, id == i))
  tmp.df = as.data.frame(cbind(metabolites.all.mean.matrix.present[i,], t(proteins.mean.matrix.present[rownames(proteins.mean.matrix.present) %in% tmp.i$gene,])))
  
  matched.network = droplevels(tmp.i[match(colnames(tmp.df[,-1]), tmp.i$gene),])
  
  cols_orfs = colnames(tmp.df)
  idx_na = which(is.na(orf2name$gene_name[match(colnames(tmp.df), orf2name$ORF)]))
  colnames(tmp.df) = orf2name$gene_name[match(colnames(tmp.df), orf2name$ORF)]
  colnames(tmp.df)[idx_na] = cols_orfs[idx_na]
  colnames(tmp.df)[1] = i
  tmp.df.scaled = na.omit(as.data.frame(scale(tmp.df, center=T, scale=T)))
  
  #PCA if metabolite has more connections than samples
  if ( sum(!is.na(tmp.df.scaled[,1])) < ncol(tmp.df.scaled[,-1]) ) {
    tmp.df.scaled = as.data.frame(scale(tmp.df, center=T, scale=T))
    s = prcomp(tmp.df.scaled[,-1])
    tmp.s = summary(s)
    include_components = sum(tmp.s$importance[3,] < 0.99)
    tmp.df.scaled = as.data.frame(cbind(tmp.df.scaled[,1],s$x[,1:include_components]))
    colnames(tmp.df.scaled)[1] = i
  }
  
  
  sample_size = round(nrow(tmp.df.scaled)*2/3)
  sample_size = round(nrow(tmp.df.scaled)*0.8)
  
  total_samples = nrow(tmp.df.scaled)
  k = 100
  sample.matrix = matrix(rep(0,k*sample_size), nrow=k)
  
  set.seed(123)
  for (ik in 1:k) {
    selected = sample(1:total_samples, size=sample_size)  
    sample.matrix[ik,] = selected
  }

  formulas = c() #storing models
  
  #choosing best model using exhaustive approach
  
  if (ncol(tmp.df.scaled)-1 < 30) {
    for(j in 1:nrow(sample.matrix)) {
    
      sample.data = tmp.df.scaled[sample.matrix[j,],]
      b = regsubsets(formula(paste0(i," ~ ", ".")), data=sample.data, nbest=1, nvmax=(ncol(sample.data)-1))
      rs = summary(b)
      n_points = nrow(na.omit(sample.data))
      
      tmp.best = data.frame(n_params = apply(rs$which, 1, sum),
                            adjr2 = rs$adjr2,
                            aic = n_points*log(rs$rss/n_points) + 2*apply(rs$which, 1, sum))
      best_idx = which(tmp.best$aic == min(tmp.best$aic))
      F1 = paste0(i," ~ ", paste(sort(names(which(rs$which[best_idx,-1] == TRUE))), collapse=" + "))
      formulas = c(formulas, F1)
    }
        
  } else { #using step approach
    
    for(j in 1:nrow(sample.matrix)) {
      #j = 1
      
      sample.data = tmp.df.scaled[sample.matrix[j,],]
      
      null.lm = lm(formula=paste0(i," ~ ", "1"), data=sample.data)
      full.lm = lm(formula=paste0(i," ~ ", "."), data=sample.data)
      
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
  
  best.models = na.omit(names(sort(-table(formulas))[1:3]))
  
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
                       p.value = ifelse(is.null(fit.s$fstatistic[1]),NA, (1 - pf(fit.s$fstatistic[1],fit.s$fstatistic[2],fit.s$fstatistic[3]))))
      
      sumarries.df = rbind(sumarries.df, tmp)  
    }
    fits.list[[j]] = m.list
  } 
  
  sumarries.df.stats = sumarries.df %>% group_by(model,data) %>% summarize( r.squared =  r.squared[1],
                                                                            adj.r.squared = adj.r.squared[1],
                                                                            p.value = p.value[1])
  sumarries.df.stats = sumarries.df.stats %>% group_by(model) %>% 
    mutate(median.adj.r.squared = median(adj.r.squared, na.rm=T)) %>% ungroup() %>% 
    arrange(desc(median.adj.r.squared))
  
  all.data = na.omit(tmp.df.scaled)
  F1.best = as.formula(best.models[sumarries.df.stats$model[1]])
  fit = lm(formula=F1.best, data=all.data)
  diag.plots1 = my_diagnostics(fit)  
  
  outliers = outlierTest(fit)
  cooks_thr = 4/(nrow(na.omit(tmp.df.scaled)) - length(fit$coefficients)-2)
  
  all.data = droplevels(all.data[!(rownames(all.data) %in% names(outliers$bonf.p < 0.05)),])
  #all.data = droplevels(all.data[!(rownames(all.data) %in% names(which(cooks.distance(fit) == max(cooks.distance(fit))))),])
  clean.data = droplevels(all.data[!(rownames(all.data) %in% names(which((cooks.distance(fit) > cooks_thr) == TRUE))),])
  
  beta.changes = data.frame(((fit$coefficients + lm.influence(fit)$coef) - fit$coefficients )/fit$coefficients)
  beta.changes$sample = rownames(beta.changes)
  beta.changes.long = melt(beta.changes[,-1], id.vars="sample")
  hist(beta.changes.long$value)
  n_thr = ceiling(length(unique(beta.changes.long$variable))/3)
  beta.changes.long$sample[abs(beta.changes.long$value) > 1]
  filter(beta.changes.long, abs(value) > 1)
  
  fit.after = lm(formula=F1.best, data=clean.data)
  diag.plots2 = my_diagnostics(fit.after, influence=F)    
  
  p_before = plot_grid(diag.plots1[[1]],
                       diag.plots1[[2]],
                       diag.plots1[[3]],
                       diag.plots1[[4]],
                       diag.plots1[[5]],
                       diag.plots1[[6]],
                       diag.plots1[[7]],
                       diag.plots1[[8]], labels = c(paste(i,"before")), ncol=2)  
  
  p_after = plot_grid(diag.plots2[[1]],
                      diag.plots2[[2]],
                      diag.plots2[[3]],
                      diag.plots2[[4]],
                      diag.plots2[[5]],
                      diag.plots2[[6]],
                      diag.plots2[[7]],
                      diag.plots2[[8]], labels = c(paste(i,"after")), ncol=2)  
  
    
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
  #g = grid.arrange(p1,p2,p3,p4)
  
  g = arrangeGrob(p1,p2,p3,p4, main=i)
    
  file_name = paste(fun_name,i,"summary.pdf", sep=".")
  file_path = paste(models_dir, file_name, sep="/")
  ggsave(filename=file_path, width=11.69 + 0.2*(11.69), height=8.27 + 0.2*(8.27), plot=g)
  
      
  file_name = paste(fun_name,i,"diag_before.pdf", sep=".")
  file_path = paste(models_dir, file_name, sep="/")
  save_plot(filename=file_path, plot=p_before, base_height=15, base_aspect_ratio=1)
  
  file_name = paste(fun_name,i,"diag_after.pdf", sep=".")
  file_path = paste(models_dir, file_name, sep="/")
  save_plot(filename=file_path, plot=p_after, base_height=15, base_aspect_ratio=1)
  
  #####
#   p1 = tableGrob(stats)
#   
#   tmp.pie = c(stats$product, stats$substrate)
#   tmp.lbs = paste(paste("product", stats$product , round(tmp.pie/sum(tmp.pie),2)[1]), 
#                   paste("substate", stats$substrate , round(tmp.pie/sum(tmp.pie),2)[2]), consolidate=T)
#   
#   p2 = ggplot(melt(stats %>% dplyr::select(product,substrate)), aes(x = "", y = value, fill = variable)) +
#     geom_bar(width = 1, stat = "identity") +
#     scale_fill_manual(values = c("red", "yellow")) +
#     coord_polar("y") +
#     ggtitle(tmp.lbs)
#   cook = cooks.distance(best.model)
#   best.model = lm(formula=as.formula(paste(i," ~ ", paste(names(result$model)[-1],   collapse=" + "))), data=result$model, subset = (cook<max(cook)))
#   robust.model = rlm(formula=as.formula(paste(i," ~ ", paste(names(result$model)[-1],   collapse=" + "))), data=result$model)
#   plot(robust.model)
#   tmp.i[match(names(result$coefficients)[-1], tmp.i$gene),]
    
}


## -- PCA metabolites ---- 

metabolites.binded = rbind(metabolites.long.mean, aa_michael.combat.long.mean)
metabolites.all.mean.df = dcast(metabolites.binded, formula=variable~ORF, value.var="mean")



metabolites.all.mean.matrix = as.matrix(metabolites.all.mean.df[,-1])
rownames(metabolites.all.mean.matrix) = metabolites.all.mean.df$variable
metabolites.all.mean.matrix.t = t(metabolites.all.mean.matrix)
metabolites.all.mean.matrix = t(metabolites.all.mean.matrix.t[complete.cases(metabolites.all.mean.matrix.t),])


both.present = intersect(colnames(proteins.mean.matrix), colnames(metabolites.all.mean.matrix))
proteins.mean.matrix.present = proteins.mean.matrix[,match(both.present, colnames(proteins.mean.matrix))]
metabolites.all.mean.matrix.present = metabolites.all.mean.matrix[,match(both.present, colnames(metabolites.all.mean.matrix))]


a = hclust(dist(t(proteins.mean.matrix.present)))
b = hclust(dist(t(metabolites.all.mean.matrix.present)))

dendextend::cor_cophenetic(a,b)

## -- PCA of metabolites ##### 
metabolites.mean.imputed.df = dcast(metabolites.long.mean, formula=variable~ORF, value.var="mean")
metabolites.mean.imputed.matrix = as.matrix(metabolites.mean.imputed.df[,-1])
rownames(metabolites.mean.imputed.matrix) = metabolites.mean.imputed.df$variable

#removing outliers
metabolites =  t(metabolites.all.mean.matrix)
metabolites.cov = cov.rob(metabolites)
md = mahalanobis(metabolites, center=metabolites.cov$center, cov=metabolites.cov$cov)
n = nrow(metabolites); p=ncol(metabolites)

plot(qchisq((1:n)/(n+1),p), sort(md), 
     xlab = expression(paste(chi^2," quantiles")),
     ylab = "Sorted Machalanobis distances")
abline(0,1)

p = recordPlot()
plots.list = lappend(plots.list, p)

#decided to remove 5 points
metabolites.f = metabolites
metabolites.f = metabolites[!(rownames(metabolites) %in% names(sort(-md)[1:14])),]
metabolites.f = metabolites.f[rownames(metabolites.f) %in% droplevels(exp_metadata$ORF[exp_metadata$type == "Kinase"]),] #taking only kinase deletions

toPlot = t(metabolites.all.mean.matrix)
s1 = prcomp(toPlot, scale=T)
# toPlot = toPlot[!rownames(toPlot) %in% names(c(which(s1$x[,1] == max(s1$x[,1])), which(s1$x[,2] == max(s1$x[,2])))),]
# s1 = prcomp(toPlot)
par(mfrow=c(1,2))
xPC = 1
yPC = 2

biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)],
       xlabs=orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)],
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
abline(h=0,v=0)

xPC = 1
yPC = 3

biplot(s1$x[,c(xPC,yPC)],s1$rotation[,c(xPC,yPC)],
       xlabs=orf2name$gene_name[match(rownames(s1$x), orf2name$ORF)],
       xlab=paste(paste0("PC",xPC), round(s1$sdev[xPC]^2/sum(s1$sdev^2),2)),
       ylab=paste(paste0("PC",yPC), round(s1$sdev[yPC]^2/sum(s1$sdev^2),2)))
abline(h=0,v=0)

p = recordPlot()
plots.list = lappend(plots.list, p)


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 


