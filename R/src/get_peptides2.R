#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "get_peptides2"

## ---- data_load ----
load("./R/objects/peptides.data.RData") 
load("./R/objects/exp_metadata._clean_.RData")

## ---- selecting peptides based on spectronaut Q-value
peptides.data$EG.StrippedSequence = factor(peptides.data$EG.StrippedSequence)
peptides.data$R.Label = factor(peptides.data$R.Label)

peptides.data = tbl_df(peptides.data[peptides.data$R.Label %in% exp_metadata$sample_name,])

peptides.peak_sums <- group_by(peptides.data, batch_date, batch, batch.exp.n, R.Label, EG.StrippedSequence) %>%
  dplyr::summarise(count = n(),
                   signal = FG.TotalPeakArea[1],
                   EG.Qvalue = EG.Qvalue[1]) %>% group_by(R.Label, EG.StrippedSequence) %>% distinct(R.Label, EG.StrippedSequence) #TODO: in future batch variable (spectronaut batch) has to be removed



peptides.peak_sums$aquisition_date = exp_metadata$aquisition_date[match(peptides.peak_sums$R.Label, exp_metadata$sample_name)]


peptides.peak_sums.stats = peptides.peak_sums %>% group_by(R.Label, aquisition_date) %>% summarize(sum = sum(signal),
                                                                                        median = median(signal),
                                                                                        shorth = genefilter::shorth(signal))
peptides.peak_sums.stats$aquisition_date.str = as.Date(strptime(peptides.peak_sums.stats$aquisition_date, format="%Y-%m-%d %H:%M:%S"))

peptides.peak_sums.stats.mix = peptides.peak_sums.stats[grep(x=peptides.peak_sums.stats$R.Label, pattern="mix", ignore.case=T),]
peptides.peak_sums.stats.wt  = peptides.peak_sums.stats[peptides.peak_sums.stats$R.Label %in% exp_metadata$sample_name[exp_metadata$ORF == "WT"],]

library(scales)
p = ggplot(peptides.peak_sums.stats, aes(x=aquisition_date.str, y=median)) + 
       geom_point()+
       geom_point(data=peptides.peak_sums.stats.mix,aes(x=aquisition_date.str, y=median),col="red") + #MIX
       geom_point(data=peptides.peak_sums.stats.wt,aes(x=aquisition_date.str, y=median),col="blue") + #WT   
       scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%B-%d"))

plots.list = lappend(plots.list, p)



tmp.peptides = unique(as.character(peptides.peak_sums$EG.StrippedSequence))

#plotting random proteins
set.seed(123)
tmp.selected = sample(tmp.peptides,3)
toPlot = peptides.peak_sums[peptides.peak_sums$EG.StrippedSequence %in% tmp.selected,]
toPlot$aquisition_date.str = as.Date(strptime(toPlot$aquisition_date, format="%Y-%m-%d %H:%M:%S"))


toPlot.mix = toPlot[grep(x=toPlot$R.Label, pattern="mix", ignore.case=T),]
toPlot.wt  = toPlot[toPlot$R.Label %in% exp_metadata$sample_name[exp_metadata$ORF == "WT"],]

p = ggplot(toPlot, aes(x=aquisition_date.str, y=signal, col=batch.exp.n)) + 
          geom_point() +
          geom_point(data=toPlot.mix,aes(x=aquisition_date.str, y=signal),col="red") + #MIX
          geom_point(data=toPlot.wt, aes(x=aquisition_date.str, y=signal),col="blue") + #WT   
          scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d"))+
          facet_wrap(~EG.StrippedSequence, scales="free")


## -- batch clustering ----

exp_clusters = tbl_df(exp_metadata)
exp_clusters$aquisition_date.str = as.POSIXct(strftime(exp_clusters$aquisition_date, format="%Y-%m-%d %H:%M:%S"))

exp_clusters$kmeans_5 = kmeans(exp_clusters$aquisition_date.str, 5)$cl
exp_clusters$kmeans_6 = kmeans(exp_clusters$aquisition_date.str, 6)$cl
exp_clusters$kmeans_7 = kmeans(exp_clusters$aquisition_date.str, 7)$cl
exp_clusters$kmeans_8 = kmeans(exp_clusters$aquisition_date.str, 8)$cl

tmp.clusters = exp_clusters %>% select(sample_name, kmeans_5, kmeans_6, kmeans_7, kmeans_8)
tmp.clusters.long = melt(as.data.frame(tmp.clusters), id.vars="sample_name")


set.seed(123)
tmp.selected = sample(tmp.peptides,2)
toPlot = droplevels(peptides.peak_sums[peptides.peak_sums$EG.StrippedSequence %in% tmp.selected,])
toPlot$aquisition_date.str = as.Date(strptime(toPlot$aquisition_date, format="%Y-%m-%d %H:%M:%S"))


toPlot.merged = merge(toPlot, tmp.clusters.long, by.x="R.Label", by.y="sample_name")
toPlot.mix = toPlot.merged[grep(x=toPlot.merged$R.Label, pattern="mix", ignore.case=T),]
toPlot.wt  = toPlot.merged[toPlot.merged$R.Label %in% exp_metadata$sample_name[exp_metadata$ORF == "WT"],]
toPlot.merged$value = factor(toPlot.merged$value)

p = ggplot(toPlot.merged, aes(x=aquisition_date.str, y=signal, col=value)) + 
  geom_point() +
  geom_point(data=toPlot.mix,aes(x=aquisition_date.str, y=signal),col="red") + #MIX
  geom_point(data=toPlot.wt, aes(x=aquisition_date.str, y=signal),col="blue") + #WT   
  scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d"))+
  facet_wrap(variable~EG.StrippedSequence, scales="free")
plots.list = lappend(plots.list, p)


## 7 batch clusters decided ##

exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_clusters$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
exp_metadata$batch_kmeans = kmeans(exp_metadata$aquisition_date.str, 7)$cl



## -- filtering based on Q-value ----
qvalues.stats <- peptides.peak_sums %>% group_by(EG.StrippedSequence) %>%  dplyr::summarise(qvalue.median = median(EG.Qvalue))
p = ggplot(qvalues.stats, aes(x = log(qvalue.median))) + 
  geom_density() +
  geom_vline(xintercept = log(0.01))
plots.list = lappend(plots.list, p)

peptides.peak_sums = merge(peptides.peak_sums, qvalues.stats, by = c("EG.StrippedSequence"))
peptides.peak_sums = tbl_df(peptides.peak_sums)
peptides.peak_sums.f = filter(peptides.peak_sums, qvalue.median <= 0.01)

peptides.peak_sums.f$T_signal = with(peptides.peak_sums.f, log(signal))

thr_remove = 0 #removing thr_remove/2% from each side of data
peptides.peak_sums.trimmed = peptides.peak_sums.f

if (thr_remove > 0) {
  message(paste("Removing ", thr_remove," fraction of dataset", sep=""))
  peptides.peak_sums.trimmed = ddply(peptides.peak_sums, .(R.Label), 
                                     .fun = function(z) {
                                       v<<-z
                                       tmp.data = z[z$T_signal < quantile(z$T_signal, 1 - thr_remove/2) & 
                                                      z$T_signal > quantile(z$T_signal, 0 + thr_remove/2),]
                                       
                                       return(droplevels(tmp.data))
                                     })  
}

file_name = "peptides.peak_sums.trimmed.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.peak_sums.trimmed, file=file_path) 



## -- adjusting for batch effects ----
peptides.df = dcast(peptides.peak_sums.trimmed, formula=EG.StrippedSequence~R.Label, value.var="T_signal")
peptides.matrix = as.matrix(peptides.df[,-1])
rownames(peptides.matrix) = peptides.df$EG.StrippedSequence

pheno = exp_metadata[match(colnames(peptides.matrix), exp_metadata$sample_name),]

mod = model.matrix(~as.factor(ORF), data=pheno)

peptides.matrix.combat = ComBat(peptides.matrix, batch=pheno$batch_kmeans, mod=mod, par.prior=T)
peptides.matrix.quant.combat = ComBat(normalizeQuantiles(peptides.matrix), batch=pheno$batch_kmeans, mod=mod, par.prior=T)


file_name = "peptides.matrix.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix.combat,file=file_path) 

file_name = "peptides.matrix.quant.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix.quant.combat,file=file_path) 


# X = model.matrix(~ORF + batch_kmeans, data=pheno)
# fit <- lmFit(peptides.matrix, X)
# 
# Batch<- fit$coef[,ncol(X)]%*%t(X[,ncol(X)])
# Signal<- fit$coef[,1:(ncol(X)-1)]%*%t(X[,1:(ncol(X)-1)])
# after  = peptides.matrix - Batch


before = peptides.matrix
after = peptides.matrix.combat

# pca
message("plotting PCA results")

file_name = "PCA_batch_effects.pdf"
file_path = paste(figures_dir, file_name, sep="/")
pdf(file_path, width=11.7+0.1*11.7, height=8.27+0.1*8.27)
par(pty="s", mfrow=c(1,2))



pca = prcomp(t(before), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], cex=1.5, cex.lab=1.5, col=pheno$batch_date, pch=16, main="Before adjustments for batch effects", 
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
pca.mix = pca$x[grep(x=rownames(pca$x), pattern="mix", ignore.case=T),]
pca.wt = pca$x[match(pheno$sample_name[pheno$ORF=="WT"], rownames(pca$x)),]
points(pca.mix[,1], pca.mix[,2], pch=8, col="black", cex=3)
points(pca.wt[,1], pca.wt[,2], pch=2, col="black", cex=3)
text(pca$x[,1], pca$x[,2], labels=as.numeric(pheno$batch.exp), cex=0.5)

pca = prcomp(t(after), scale.=T)
#pca = prcomp(t(proteins.matrix.f., scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], cex=1.5, cex.lab=1.5, col=pheno$batch_date, pch=16, main="After adjustments for batch effects",
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
pca.mix = pca$x[grep(x=rownames(pca$x), pattern="mix", ignore.case=T),]
pca.wt = pca$x[match(pheno$sample_name[pheno$ORF=="WT"], rownames(pca$x)),]
points(pca.mix[,1], pca.mix[,2], pch=8, col="black", cex=3)
points(pca.wt[,1], pca.wt[,2], pch=2, col="black", cex=3)
text(pca$x[,1], pca$x[,2], labels=as.numeric(pheno$batch.exp), cex=0.5)


p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()



#individual examples
before = peptides.matrix
after =  peptides.matrix.combat


before.long = tbl_df(melt(before, id.vars="rownames"))
names(before.long) = c("EG.StrippedSequence", "R.Label", "signal")
before.long$aquisition_date =exp_metadata$aquisition_date[match(before.long$R.Label, exp_metadata$sample_name)]
before.long$batch_kmeans = factor(exp_metadata$batch_kmeans[match(before.long$R.Label, exp_metadata$sample_name)])
before.long$category = "before"

after.long = tbl_df(melt(after, id.vars="rownames"))
names(after.long) = c("EG.StrippedSequence", "R.Label", "signal")
after.long$aquisition_date =exp_metadata$aquisition_date[match(after.long$R.Label, exp_metadata$sample_name)]
after.long$batch_kmeans = factor(exp_metadata$batch_kmeans[match(after.long$R.Label, exp_metadata$sample_name)])
after.long$category = "after"


toPlot = rbind(before.long, after.long)
toPlot$category = factor(toPlot$category, levels=c("before", "after"))

set.seed(123)
tmp.selected = sample(tmp.peptides,5)
toPlot = toPlot[toPlot$EG.StrippedSequence %in% tmp.selected,]

toPlot$aquisition_date.str = as.Date(strptime(toPlot$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
toPlot.mix = toPlot[grep(x=toPlot$R.Label, pattern="mix", ignore.case=T),]
toPlot.mix.stats = toPlot.mix %>% group_by(EG.StrippedSequence, batch_kmeans, category) %>% summarize(mean = mean(exp(signal)))

toPlot.wt  = toPlot[toPlot$R.Label %in% exp_metadata$sample_name[exp_metadata$ORF == "WT"],]

p = ggplot(toPlot, aes(x=aquisition_date.str, y=exp(signal), col=batch_kmeans)) + 
  geom_point(size=3) +
  geom_point(data=toPlot.mix,aes(x=aquisition_date.str, y=exp(signal)),col="red") + #MIX
  geom_point(data=toPlot.wt, aes(x=aquisition_date.str, y=exp(signal)),col="blue") + #WT   
  scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d"))+
  facet_grid(EG.StrippedSequence~category, scales="free")


plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 
