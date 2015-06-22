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
names(metabolites.matrix.combat.long) = c("sample_id","variable","value" )

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

load("./R/objects/yeast.model._load_.RData")
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



metabolite2name = data.frame(name = c("ADP","ADP", "AMP", "AMP", "ATP", "ATP", "DHAP", "E4P", "F6P", 
                      "G3P", "G6P...F6P","G6P...F6P", "G6P...F6P", "Glu", "Glu",      
                      "PEP", "Pyr", "R5P", "S7P", "X1.6.FP.", 
                      "X2...3.PG", "X2...3.PG", "X5P...Rl5P", "X5P...Rl5P", "X6PGlu"),
       model_name = c("ADP","ADPM", "AMP", "AMPM", "ATP","ATPM", "Glycerone phosphate", "D-Erythrose 4-phosphate", "beta-D-Fructose 6-phosphate",
                      "D-Glyceraldehyde 3-phosphate", "beta-D-Glucose 6-phosphate", "alpha-D-Glucose 6-phosphate","beta-D-Fructose 6-phosphate", "beta-D-Glucose", "aplha-D-Glucose",
                      "Phosphoenolpyruvate", "Pyruvate", "D-Ribose 5-phosphate", "Sedoheptulose 7-phosphate", "beta-D-Fructose 1,6-bisphosphate",
                      "2-Phospho-D-glycerate", "3-Phospho-D-glycerate", "D-Xylose-5-phosphate", "D-Ribulose 5-phosphate", "6-Phospho-D-gluconate"))


aa2name = data.frame(name = c("alanine", "alanine", "aspartate", "aspartate", "glutamate", "glutamate", 
                              "phenylalanine", "glycine", "glycine", "histidine", 
                              "isoleucine", "isoleucine", "lysine","lysine", "leucine","leucine", 
                              "methionine", "asparagine", "proline", "proline",
                              "glutamine", "arginine", "serine", "serine",
                              "threonine", "threonine", "valine", "tryptophan", "tryptophan", "tyrosine"),
                    model_name = c("L-Alanine","L-AlanineM", "L-Aspartate", "L-AspartateM", "L-Glutamate", "GlutamateM", 
                                   "L-Phenylalanine", "Glycine","GlycineM", "L-Histidine", 
                                   "L-Isoleucine", "L-IsoleucineM", "L-Lysine","L-LysineM", "L-Leucine","L-LeucineM",
                                   "L-Methionine", "L-Asparagine", "L-Proline", "L-ProlineM", 
                                   "L-Glutamine", "L-Arginine", "L-Serine",  "L-SerineM",
                                   "L-Threonine", "L-ThreonineM", "L-Valine", "L-Tryptophan", "L-TryptophanM", "L-Tyrosine"))


yeast.model.merged = droplevels(merge(yeast.model, rbind(metabolite2name,aa2name), by.x="metabolite", by.y="model_name"))


library("faraway")
library("car")

models.list = list()
fit.list = list()

for(j in 1:length(levels(yeast.model.merged$name))) {
    
  i = levels(yeast.model.merged$name)[j]
  i = "ATP"
  tmp.i = droplevels(subset(yeast.model.merged, name == i))
  tmp.df = as.data.frame(cbind(metabolites.all.mean.matrix.present[i,], t(proteins.mean.matrix.present[rownames(proteins.mean.matrix.present) %in% tmp.i$gene,])))
  matched.network = droplevels(tmp.i[match(colnames(tmp.df[,-1]), tmp.i$gene),])
  
    
  stats = data.frame(metabolite=i,
             n_enzymes = ncol(tmp.df[,-1]),
             coverage = round(ncol(tmp.df[,-1])/length(grep("^U\\d+", na.omit(unique(tmp.i$gene)), ignore.case=T, invert=T, value=T)),2),
             n_reactions = length(unique(matched.network$reaction)),
             product = as.vector(table(matched.network$side)[1]),
             substrate = as.vector(table(matched.network$side)[2]),
             ps.ratio  = round(as.vector(table(tmp.i$side)[1]/table(tmp.i$side)[2]),2))
  
  if (stats$n_enzymes == 0) {
    p1 = tableGrob(stats)
    grid.arrange(p1)
    p = recordPlot()
    dev.off()
    
    models.list = lappend(models.list, p)
    next
  }
  
  colnames(tmp.df)[1] = i
  tmp.df.scaled = as.data.frame(scale(tmp.df, center=F, scale=F))
              
  null.lm = lm(formula=paste0(i," ~ ", "1"), data=tmp.df.scaled)
  full.lm = lm(formula=paste0(i," ~ ", "."), data=tmp.df.scaled)
  
  result = step(null.lm, scope=list(lower=null.lm, upper=full.lm), direction="both", steps=100000)
  
  fit = eval(result$call)
  outliers = outlierTest(fit)
  
  
  #removing outliers
  while(outliers$signif) {
        
    tmp.df.scaled = tmp.df.scaled[!(rownames(tmp.df.scaled) %in% names(outliers$bonf.p < 0.05)),]
    
    null.lm = lm(formula=paste0(i," ~ ", "1"), data=tmp.df.scaled)
    full.lm = lm(formula=paste0(i," ~ ", "."), data=tmp.df.scaled)
    
    result = step(null.lm, scope=list(lower=null.lm, upper=full.lm), direction="both", steps=100000)
    fit = eval(result$call)
    outliers = outlierTest(fit)
  }
  sum = summary(fit)
  
  par(mfrow=c(2,2))
  plot(fit)
  
  
  fit.list[[j]] = fit

  # final model stats
  p1 = tableGrob(stats)
  
  tmp.pie = c(stats$product, stats$substrate)
  tmp.lbs = paste(paste("product", stats$product , round(tmp.pie/sum(tmp.pie),2)[1]), 
                  paste("substate", stats$substrate , round(tmp.pie/sum(tmp.pie),2)[2]), consolidate=T)
  
  p2 = ggplot(melt(stats %>% dplyr::select(product,substrate)), aes(x = "", y = value, fill = variable)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = c("red", "yellow")) +
    coord_polar("y") +
    ggtitle(tmp.lbs)
#   cook = cooks.distance(best.model)
#   best.model = lm(formula=as.formula(paste(i," ~ ", paste(names(result$model)[-1],   collapse=" + "))), data=result$model, subset = (cook<max(cook)))
#   robust.model = rlm(formula=as.formula(paste(i," ~ ", paste(names(result$model)[-1],   collapse=" + "))), data=result$model)
#   plot(robust.model)
#   tmp.i[match(names(result$coefficients)[-1], tmp.i$gene),]
  
  p3 = tableGrob(round(sum$coefficients,5))
  model.summary = data.frame(r.sq = sum$r.squared,
                             adj.r.sq = sum$adj.r.squared,
                             fstats = ifelse(is.null(sum$fstatistic[1]),NA,sum$fstatistic[1]),
                             numdf = ifelse(is.null(sum$fstatistic[2]),NA, sum$fstatistic[2]),
                             dendf = ifelse(is.null(sum$fstatistic[3]),NA, sum$fstatistic[3]),
                             pval = anova(fit)$'Pr(>F)'[1])

  p4 = tableGrob(model.summary)
  
  grid.arrange(p1,p3,p4, ncol=1)

  p = recordPlot()
  
  dev.off()
  
  models.list = lappend(models.list, p)
  
  par(pty="s", mfrow=c(3,2))
  plot(fit)
  title(main=i, outer=T)
  
  if (ncol(fit$model) > 1) {
    y_star = stats::predict(object=fit, tmp.df.scaled)
    #y_star = stats::predict(object=robust.model, tmp.df.scaled)
    plot(y_star,tmp.df.scaled[,1], xlab="Predicted", ylab=paste("Measured",i, "signal intensity"))
    legend("topleft", paste( paste("R^2=",round(sum$r.squared,2))),bty="n")
    
    p = recordPlot()
    dev.off()  
    
  }
models.list = lappend(models.list, p)
  
}

file_name = paste(fun_name, "models.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(models.list, filename=file_path, type="l")



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



