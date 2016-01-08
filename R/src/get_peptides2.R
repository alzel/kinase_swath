#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "get_peptides2"

## ---- data_load ----
load("./R/objects/peptides.data._clean_.RData") 
load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/data.frame.peptides.raw._load_.RData")


peptides.data = tbl_df(peptides.data[peptides.data$R.Label %in% exp_metadata$sample_name,])

peptides.data = peptides.data %>% distinct(R.Label, FG.Id, fragment)

stopifnot(any(as.vector((peptides.data %>% 
                           select(R.Label, batch) %>% 
                           distinct(R.Label, batch) %>% 
                           group_by(R.Label) %>% 
                           mutate(n = n()) %>% select(n))[,2]) == 1))


# Intructions from Oliver Bernhardt
#for each precursor (modified sequence + charge) you do:
#  1) select all fragments that have "PossibleInterference== false" in ALL runs to use for quantification
#  2) If the set obtained in step 1 is < 3 then select all remaining fragments ("PossibleInterference== true" in AT LEAST ONE run) for further ranking
#  3) calculate the average interference score per fragment across all runs (for the set from step 2)
#  4) order these fragments based on their cross run average interference score (from low to high).
#  5) add fragments from this ordered list to the set you already obtained in step 1 till you have AT LEAST 3 fragments selected for quantification.



quant_thr = 3
peptides.data.tmp = peptides.data %>% 
  group_by(FG.Id, fragment) %>% 
  mutate(F.GoodInAll = ifelse(any(F.PossibleInterference == "True"),0,1),
         mean.F.InterferenceScore = mean(F.InterferenceScore, na.rm=T),
         median.F.InterferenceScore = median(F.InterferenceScore, na.rm=T)) %>%
  group_by(R.Label, FG.Id) %>%
  mutate(count.F.GoodInAll = sum(F.GoodInAll)) %>%
  group_by(R.Label, FG.Id) %>%
  arrange(median.F.InterferenceScore) %>%
  group_by(R.Label, FG.Id) %>%
  mutate(miss.count = count.F.GoodInAll - quant_thr,
         topquant = ifelse(miss.count < 0, abs(miss.count), 0),
         quanty_all = T,
         F.IsQuantified = ifelse(topquant != 0, 
                                 which(quanty_all) %in% c(which(F.GoodInAll == 1), which(F.GoodInAll != 1)[1:topquant[1]]), 
                                 which(quanty_all) %in% which(F.GoodInAll == 1)))

fragments.data = peptides.data.tmp %>% group_by(R.Label, FG.Id) %>%
  mutate(FG.TotalPeakArea_new = sum(F.PeakArea[F.IsQuantified])) %>% group_by(R.Label, EG.StrippedSequence, FG.Id ) %>% arrange()

fragments.data$miss.count = NULL
fragments.data$topquant = NULL
fragments.data$quanty_all = NULL


fragments.data = fragments.data %>% group_by(FG.Id) %>%  mutate(qvalue.median = median(unique(EG.Qvalue)))

file_name = "fragments.data.RData"
file_path = paste(output_dir, file_name, sep="/")
save(fragments.data, file=file_path) 

p = ggplot(fragments.data, aes(x = log(qvalue.median))) + 
  geom_density() +
  geom_vline(xintercept = log(0.01))

plots.list = lappend(plots.list, p)

## ---- selecting peptides based on spectronaut Q-value

fragments.data.f = filter(fragments.data, qvalue.median <= 0.01, F.IsQuantified == T)

# peptides.data$EG.StrippedSequence = factor(peptides.data$EG.StrippedSequence)
# peptides.data$R.Label = factor(peptides.data$R.Label)
# 
# peptides.data = tbl_df(peptides.data[peptides.data$R.Label %in% exp_metadata$sample_name,])
# 
# peptides.peak_stats <- peptides.data %>% group_by(R.Label, EG.StrippedSequence, FG.Id, batch) %>%
#                                          summarize(count = n(),
#                                                    count.PossibleInterference = length(F.PeakArea[F.PossibleInterference != "True"]),
#                                                    sum.F.PeakArea = sum(F.PeakArea[F.PossibleInterference != "True"]),
#                                                    sum.F.PeakArea.all = sum(F.PeakArea),
#                                                    signal = FG.TotalPeakArea[1])
# 
# 
# toPlot = droplevels(peptides.peak_stats[peptides.peak_stats$FG.Id %in% sample(x=levels(peptides.peak_stats$FG.Id), size=50),])
# p = ggplot(toPlot, aes(x=jitter(as.numeric(FG.Id)), y=jitter(count.PossibleInterference))) +
#       geom_point() +
#       scale_x_continuous(breaks=as.numeric(factor(levels(toPlot$FG.Id))), labels=levels(toPlot$FG.Id))+
#       xlab("FG.Id") +
#       theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
# plots.list = lappend(plots.list, p)





# peptides.peak_sums <- group_by(peptides.data, batch_date, batch, batch.exp.n, R.Label, EG.StrippedSequence) %>%
#   dplyr::summarise(count = n(),
#                    signal = FG.TotalPeakArea[1],
#                    EG.Qvalue = EG.Qvalue[1]) %>% group_by(R.Label, EG.StrippedSequence) %>% distinct(R.Label, EG.StrippedSequence) #TODO: in future batch variable (spectronaut batch) has to be removed

peptides.peak_sums <- fragments.data.f %>% 
                      group_by(batch,R.Label, EG.StrippedSequence, FG.Id) %>%
                      dplyr::summarise(count = n(),  signal = FG.TotalPeakArea_new[1]) %>%
                      group_by(batch,R.Label, EG.StrippedSequence) %>%
                      dplyr::summarise(count = sum(count),  
                                       signal = sum(signal))


peptides.peak_sums$aquisition_date = exp_metadata$aquisition_date[match(peptides.peak_sums$R.Label, exp_metadata$sample_name)]
peptides.peak_sums$batch.exp.n = exp_metadata$batch.exp.n[match(peptides.peak_sums$R.Label, exp_metadata$sample_name)]


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

set.seed(1234)
exp_clusters$kmeans_5 = pam(exp_clusters$aquisition_date.str, 5)$clustering
exp_clusters$kmeans_6 = pam(exp_clusters$aquisition_date.str, 6)$clustering
exp_clusters$kmeans_7 = pam(exp_clusters$aquisition_date.str, 7)$clustering
exp_clusters$kmeans_8 = pam(exp_clusters$aquisition_date.str, 8)$clustering

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
  geom_point(data=toPlot.mix,aes(x=aquisition_date.str, y=jitter(signal,50)),col="red") + #MIX
  geom_point(data=toPlot.wt, aes(x=aquisition_date.str, y=signal),col="blue") + #WT   
  scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d"))+
  facet_wrap(variable~EG.StrippedSequence, scales="free")

plots.list = lappend(plots.list, p)


## 7 batch clusters decided ##

exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
exp_metadata$batch_kmeans = pam(exp_metadata$aquisition_date.str, 7)$clustering
#exp_metadata$batch_kmeans = kmeans(exp_metadata$aquisition_date.str, 7)$cl
                                    

## -- filtering based on Q-value ----
# qvalues.stats <- peptides.peak_sums %>% group_by(EG.StrippedSequence) %>%  dplyr::summarise(qvalue.median = median(EG.Qvalue))
# p = ggplot(qvalues.stats, aes(x = log(qvalue.median))) + 
#   geom_density() +
#   geom_vline(xintercept = log(0.01))
# plots.list = lappend(plots.list, p)
# 
# peptides.peak_sums = merge(peptides.peak_sums, qvalues.stats, by = c("EG.StrippedSequence"))
# peptides.peak_sums = tbl_df(peptides.peak_sums)
#peptides.peak_sums.f = filter(peptides.peak_sums, qvalue.median <= 0.01)


peptides.peak_sums$T_signal = with(peptides.peak_sums, log(signal))

thr_remove = 0 #removing thr_remove/2% from each side of data
peptides.peak_sums.trimmed = peptides.peak_sums

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
#peptides.matrix.combat.vsn = normalizeVSN(ComBat(exp(peptides.matrix), batch=pheno$batch_kmeans, mod=mod, par.prior=T))
#peptides.matrix.vsn.combat = ComBat(normalizeVSN(exp(peptides.matrix)), batch=pheno$batch_kmeans, mod=mod, par.prior=T)
peptides.matrix.combat.quant = normalizeQuantiles(ComBat(peptides.matrix, batch=pheno$batch_kmeans, mod=mod, par.prior=T))


file_name = "peptides.matrix.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix.combat,file=file_path) 

file_name = "peptides.matrix.combat.quant.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix.combat.quant,file=file_path) 



# X = model.matrix(~ORF + batch_kmeans, data=pheno)
# fit <- lmFit(peptides.matrix, X)
# 
# Batch<- fit$coef[,ncol(X)]%*%t(X[,ncol(X)])
# Signal<- fit$coef[,1:(ncol(X)-1)]%*%t(X[,1:(ncol(X)-1)])
# after  = peptides.matrix - Batch


before = exp(peptides.matrix)
after = peptides.matrix.combat.quant



pca = prcomp(t(before), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = data.frame(x_var, y_var, type="before")

scores = as.data.frame(pca$x[,1:5])
scores$type = "before"
scores$sample.id = rownames(scores)

pca = prcomp(t(after), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = rbind(annot,data.frame(x_var, y_var, type="after"))

scores = rbind(scores, data.frame(sample.id = rownames(scores), pca$x[,1:5], type = "after"))
scores$batch_kmeans = factor(pheno$batch_kmeans[match(scores$sample.id, pheno$sample_name)])
scores$batch = factor(pheno$batch.exp.n[match(scores$sample.id, pheno$sample_name)])
scores.mix = scores[grepl(pattern="mix", ignore.case=T, x=rownames(scores)),]
scores$type = factor(scores$type, levels=c("before", "after"))

annot$text = paste(annot$x_var, annot$y_var)

library(cowplot)
p = ggplot(scores, aes(x=PC1, y=PC2)) + 
  geom_point(size=3, aes(col=batch) )+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(data=scores.mix, aes(x=PC1, y=PC2),size=3,col="black", shape=17) +
  geom_text(data = annot, aes(x=-50, y=-50, label=text)) +
  facet_wrap(~type, scales="fixed") + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = rel(1.5)))
        
file_name = paste(fun_name,"batch_effects", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width = 2*8.27)



# pca
# message("plotting PCA results")
# 
# file_name = "PCA_batch_effects.pdf"
# file_path = paste(figures_dir, file_name, sep="/")
# pdf(file_path, width=11.7+0.1*11.7, height=8.27+0.1*8.27)
# par(pty="s", mfrow=c(1,2))
# 
# 
# 
# pca = prcomp(t(before), scale.=T)
# x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
# y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
# 
# plot(pca$x[,1], pca$x[,2], cex=1.5, cex.lab=1.5, col=pheno$batch_date, pch=16, main="Before adjustments for batch effects", 
#      xlab=paste("PC1,", x_var), 
#      ylab=paste("PC2,", y_var))
# pca.mix = pca$x[grep(x=rownames(pca$x), pattern="mix", ignore.case=T),]
# pca.wt = pca$x[match(pheno$sample_name[pheno$ORF=="WT"], rownames(pca$x)),]
# points(pca.mix[,1], pca.mix[,2], pch=8, col="black", cex=3)
# points(pca.wt[,1], pca.wt[,2], pch=2, col="black", cex=3)
# text(pca$x[,1], pca$x[,2], labels=as.numeric(pheno$batch.exp), cex=0.5)
# 
# pca = prcomp(t(after), scale.=T)
# #pca = prcomp(t(proteins.matrix.f., scale.=T)
# x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
# y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
# 
# plot(pca$x[,1], pca$x[,2], cex=1.5, cex.lab=1.5, col=pheno$batch_date, pch=16, main="After adjustments for batch effects",
#      xlab=paste("PC1,", x_var), 
#      ylab=paste("PC2,", y_var))
# pca.mix = pca$x[grep(x=rownames(pca$x), pattern="mix", ignore.case=T),]
# pca.wt = pca$x[match(pheno$sample_name[pheno$ORF=="WT"], rownames(pca$x)),]
# points(pca.mix[,1], pca.mix[,2], pch=8, col="black", cex=3)
# points(pca.wt[,1], pca.wt[,2], pch=2, col="black", cex=3)
# text(pca$x[,1], pca$x[,2], labels=as.numeric(pheno$batch.exp), cex=0.5)
# 
# 
# p = recordPlot()
# plots.list = lappend(plots.list, p)
# dev.off()



#individual examples
#before = peptides.matrix
#after =  peptides.matrix.combat


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
p
plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 
