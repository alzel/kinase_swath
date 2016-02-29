#!/usr/bin/env Rscript
rm(list=ls())
source("./R/boot.R")
source("./R/functions.R")


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

## 7 batch clusters decided ##
exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
exp_metadata$batch_kmeans = pam(exp_metadata$aquisition_date.str, 7)$clustering

#removing outliers based on total Peak_Area
summary.stats <- peptides.data %>%
  select(R.Label, F.PeakArea) %>%
  group_by(R.Label) %>% 
  summarise(signal = sum(F.PeakArea)) %>%
  gather(stats, value, -R.Label)

summary.stats$batch_kmeans <- exp_metadata$batch_kmeans[match(summary.stats$R.Label, exp_metadata$sample_name)]

# Decide which samples should be removed --------------------

summary.stats = summary.stats %>% 
  group_by(stats, batch_kmeans) %>% 
  mutate(z_value = (value - mean(value, na.rm=T))/sd(value, na.rm=T),
         toRemove = ifelse(abs(z_value) > 3, 1, 0))

summary.stats$batch_kmeans = exp_metadata$batch_kmeans[match(summary.stats$R.Label, exp_metadata$sample_name)]
summary.stats$aquisition_date.str = exp_metadata$aquisition_date.str[match(summary.stats$R.Label, exp_metadata$sample_name)]
toPlot = summary.stats
toPlot$isMix = ifelse(grepl(toPlot$R.Label, pattern = "mix", ignore.case = T),1,0)

library(scales)
p <- ggplot(toPlot, aes(x=aquisition_date.str,y=value, colour=factor(batch_kmeans))) +
  geom_point() +
  geom_point(data = toPlot %>% filter(toRemove == 1), shape=4, colour="red", size=5) +
  geom_point(data = toPlot %>% filter(isMix == 1),  colour="black") +
  geom_text(data=toPlot %>% filter(toRemove == 1), 
            aes(x=aquisition_date.str, y=value, label=R.Label), vjust = 0.05, hjust = -0.1) +
  facet_wrap(~stats, scales="free", ncol = 1) +
  ylab("Total sum of peak areas") +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

plots.list = lappend(plots.list, p)

#selecting relevant for kinases paper analysis
summary.stats$isInteresting <- ifelse(summary.stats$R.Label %in% exp_metadata$sample_name[grepl(exp_metadata$type,ignore.case = T, pattern = "Kinase|Wild Type|Mix")],1,0)

file_name = paste("summary.stats", fun_name, "RData", sep = ".")
file_path = paste(output_dir, file_name, sep="/")
save(summary.stats, file=file_path) 


# Peptide summaries ----------------
selected = summary.stats %>% filter(toRemove == 0, isInteresting == 1) %>% ungroup %>% select(R.Label) %>% distinct()

# Intructions from Oliver Bernhardt
#for each precursor (modified sequence + charge) you do:
#  1) select all fragments that have "PossibleInterference== false" in ALL runs to use for quantification
#  2) If the set obtained in step 1 is < 3 then select all remaining fragments ("PossibleInterference== true" in AT LEAST ONE run) for further ranking
#  3) calculate the average interference score per fragment across all runs (for the set from step 2)
#  4) order these fragments based on their cross run average interference score (from low to high).
#  5) add fragments from this ordered list to the set you already obtained in step 1 till you have AT LEAST 3 fragments selected for quantification.

peptides.data.f <- peptides.data %>% filter(R.Label %in% selected$R.Label)

quant_thr = 3
peptides.data.tmp = peptides.data.f %>%
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

file_name = paste("fragments.data", fun_name, "RData", sep = ".")
file_path = paste(output_dir, file_name, sep="/")
save(fragments.data, file=file_path) 

p = ggplot(fragments.data, aes(x = log(qvalue.median))) + 
  geom_density() +
  geom_vline(xintercept = log(0.01))

plots.list = lappend(plots.list, p)

## ---- selecting peptides based on spectronaut Q-value
fragments.data.f = filter(fragments.data, qvalue.median <= 0.01, F.IsQuantified == T)

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
       scale_x_date(breaks = date_breaks("1 week"), minor_breaks = date_breaks("1 day"), labels=date_format("%m-%d"))
       #scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%B-%d"))
plots.list = lappend(plots.list, p)


tmp.peptides = unique(as.character(peptides.peak_sums$EG.StrippedSequence))

#plotting random proteins
set.seed(123)
tmp.selected = sample(tmp.peptides,3)
toPlot = peptides.peak_sums[peptides.peak_sums$EG.StrippedSequence %in% tmp.selected,]
toPlot$aquisition_date.str = as.Date(strptime(toPlot$aquisition_date, format="%Y-%m-%d %H:%M:%S"))


toPlot.mix = toPlot[grep(x=toPlot$R.Label, pattern="mix", ignore.case=T),]
toPlot.wt  = toPlot[toPlot$R.Label %in% exp_metadata$sample_name[exp_metadata$ORF == "WT"],]
library(scales)
p = ggplot(toPlot, aes(x=aquisition_date.str, y=signal, col=batch.exp.n)) + 
          geom_point() +
          geom_point(data=toPlot.mix,aes(x=aquisition_date.str, y=signal),col="red") + #MIX
          geom_point(data=toPlot.wt, aes(x=aquisition_date.str, y=signal),col="blue") + #WT   
          #scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d"))+
          scale_x_date(breaks = date_breaks("1 week"), minor_breaks = date_breaks("1 day"), labels=date_format("%m-%d"))+
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
  #scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d"))+
  scale_x_date(breaks = date_breaks("1 week"), minor_breaks = date_breaks("1 day"), labels=date_format("%m-%d"))+
  facet_wrap(variable~EG.StrippedSequence, scales="free")

plots.list = lappend(plots.list, p)


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
file_name = "peptides.matrix.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix,file=file_path) 


pheno = droplevels(exp_metadata[match(colnames(peptides.matrix), exp_metadata$sample_name),])

mod = model.matrix(~as.factor(ORF), data=pheno)

peptides.matrix.combat = ComBat(peptides.matrix, batch=pheno$batch_kmeans, mod=mod, par.prior=T)
file_name = "peptides.matrix.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix.combat,file=file_path) 


#peptides.matrix.combat.vsn = normalizeVSN(ComBat(exp(peptides.matrix), batch=pheno$batch_kmeans, mod=mod, par.prior=T))
#peptides.matrix.vsn.combat = ComBat(normalizeVSN(exp(peptides.matrix)), batch=pheno$batch_kmeans, mod=mod, par.prior=T)

peptides.matrix.combat.quant = normalizeQuantiles(ComBat(peptides.matrix, batch=pheno$batch_kmeans, mod=mod, par.prior=T))
file_name = "peptides.matrix.combat.quant.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix.combat.quant,file=file_path) 

peptides.matrix.quant.combat = ComBat(normalizeQuantiles(peptides.matrix), batch=pheno$batch_kmeans, mod=mod, par.prior=T)
file_name = "peptides.matrix.quant.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.matrix.quant.combat,file=file_path) 


## adjusting fragmetns for batch effects
fragments.df = dcast(fragments.data.f, formula=FG.Id+fragment~R.Label, value.var="F.PeakArea")
fragments.matrix = as.matrix(fragments.df[,-c(1,2)])
rownames(fragments.matrix) = paste(fragments.df$FG.Id, fragments.df$fragment, sep="")

pheno = exp_metadata[match(colnames(fragments.matrix), exp_metadata$sample_name),]

mod = model.matrix(~as.factor(ORF), data=pheno)

fragments.matrix.quant.combat = ComBat(normalizeQuantiles(log(fragments.matrix)), batch=pheno$batch_kmeans, mod=mod, par.prior=T)

file_name = "fragments.matrix.quant.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(fragments.matrix.quant.combat,file=file_path) 



# X = model.matrix(~ORF + batch_kmeans, data=pheno)
# fit <- lmFit(peptides.matrix, X)
# 
# Batch<- fit$coef[,ncol(X)]%*%t(X[,ncol(X)])
# Signal<- fit$coef[,1:(ncol(X)-1)]%*%t(X[,1:(ncol(X)-1)])
# after  = peptides.matrix - Batch


before = peptides.matrix
after  = peptides.matrix.quant.combat


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
  geom_point(size=3, aes(col=batch_kmeans) )+
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
  scale_x_date(breaks = date_breaks("1 week"), minor_breaks = date_breaks("1 day"), labels=date_format("%m-%d"))+
  facet_grid(EG.StrippedSequence~category, scales="free")

plots.list = lappend(plots.list, p)


# Proteins per sample -------------------------
dataset.tmp <- peptides.data.f %>% 
  select(R.Label, EG.StrippedSequence, EG.Qvalue ) %>%
  dcast(formula = "EG.StrippedSequence ~ R.Label", value.var = "EG.Qvalue", fun.aggregate = min)

dataset.tmp[,-1][is.infinite(as.matrix(dataset.tmp[,-1]))] <- NA
dataset.tmp[,-1][dataset.tmp[,-1] > 0.01] <- NA

load("./R/objects/protein_annotations_trypsin._clean_.RData")
protein_annotations.unique <- tbl_df(protein_annotations) %>% 
  ungroup() %>% 
  dplyr::select(strippedSequence, SystName) %>% 
  distinct()

measured_peptides <- unique(peptides.data.f$EG.StrippedSequence)

protein_annotations.all <- protein_annotations %>% 
  dplyr::select(strippedSequence, SystName) 


dataset.raw.f.good <- peptides.data.f %>% left_join(protein_annotations.all, by = c("EG.StrippedSequence" = "strippedSequence"))
dataset.raw.f.good <- dataset.raw.f.good %>% rename(ProteinName = SystName)

fdr_thr1 = 0.01
fdr_thr5 = 0.05

total_samples = length(unique(dataset.raw.f.good$R.Label))
entity_summaries <- dataset.raw.f.good %>%  
  ungroup() %>%
  #select(Fragment_Annotation, transition_group_id, FullPeptideName, ProteinName, filename, m_score) %>% str()
  select(ProteinName, R.Label, EG.Qvalue) %>%
  gather(stats, value, -R.Label, -EG.Qvalue) %>% 
  group_by(stats, value, R.Label) %>% 
  summarise(min.EG.Qvalue = min(EG.Qvalue, na.rm=T)) %>%
  group_by(R.Label) %>% filter(min.EG.Qvalue < fdr_thr1 ) %>%
  summarise(n = n())

p <- ggplot(entity_summaries, aes(x=n)) +
  geom_histogram() +
  xlab("Proteins per samples")
plots.list = lappend(plots.list, p)


## Statistics ######
stats_table <- data.frame(stats_name = character(),
                          value = character(), 
                          comment=character(),
                          stringsAsFactors=FALSE)

#Number of all aquired samples including outliers
stats_tmp <- data.frame(stats = "uniq_R.Label", 
                        value = length(unique(peptides.data$R.Label)), 
                        comment = "Number of all aquired samples including outliers")
stats_table <- rbind(stats_table, stats_tmp)

#number of analysed kinases
tmp.selected <- exp_metadata[match(unique(colnames(peptides.matrix)), exp_metadata$sample_name),]

stats_tmp <- data.frame(stats = "uniq_kinases", 
                        value = length(unique(tmp.selected$ORF[tmp.selected$type == "Kinase"])), 
                        comment = "Number of kinase mutants analysed")
stats_table <- rbind(stats_table, stats_tmp)

#number of samples with mix, kinases or wild_type
tmp.selected <- exp_metadata[match(unique(peptides.data.f$R.Label), exp_metadata$sample_name),]
stats_tmp <- data.frame(stats = "uniq_Interesting_samples", 
                        value = length(unique(colnames(peptides.matrix)[grepl(tmp.selected$type,ignore.case = T, pattern = "Kinase|Wild Type|Mix")])), 
                        comment = "Number of samples with mix, kinases or wild_type")
stats_table <- rbind(stats_table, stats_tmp)

#mean number of proteins detected in samples (irrespective of peptide non-uniq match)
stats_tmp <- data.frame(stats = "mean_proteins", 
                        value = mean(entity_summaries$n), 
                        comment = "#mean number of proteins detected in samples (irrespective of peptide non-uniq match)")
stats_table <- rbind(stats_table, stats_tmp)

#mean number of proteins detected in samples (irrespective of peptide non-uniq match)
stats_tmp <- data.frame(stats = "sd_proteins", 
                        value = sd(entity_summaries$n), 
                        comment = "#sd number of proteins detected in samples (irrespective of peptide non-uniq match)")
stats_table <- rbind(stats_table, stats_tmp)

library("gridExtra")
p <- tableGrob(stats_table)
plots.list = lappend(plots.list, p)


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 

