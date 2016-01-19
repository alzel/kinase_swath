#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "batch_effects"

## ---- data_load ----
load("./R/objects/peptides.data.RData") # clean.R::prepare_peptides()
load("./R/objects/experiment.map.RData")# load.R::load_batch_map()
load("./R/objects/dates_map.RData") # load.R::load_dates_map()
load("./R/objects/sample.map.RData") # load.R::load_sample_map()
load("./R/objects/sample_exp.map.RData") # load.R::load_sample_map()
load("./R/objects/peptides2orfs.RData")

#peptides.data$batch.exp.n = factor(as.numeric(peptides.data$batch.exp)) 


## ---- selecting peptides based on spectronaut Q-value
peptides.data$EG.StrippedSequence = factor(peptides.data$EG.StrippedSequence)
peptides.data$R.Label = factor(peptides.data$R.Label)

peptides.peak_sums <- group_by(peptides.data, batch_date, batch, batch.exp.n, R.Label,sample, replicate, EG.StrippedSequence) %>%
                      dplyr::summarise(count = n(),
                                       signal = FG.TotalPeakArea[1],
                                       EG.Qvalue = EG.Qvalue[1],
                                       FG.PrecursorMz = FG.PrecursorMz[1])


qvalues.stats.batch <- peptides.peak_sums %>% group_by(EG.StrippedSequence, batch.exp.n) %>%  dplyr::summarise(qvalue.median = median(EG.Qvalue))
p = ggplot(qvalues.stats.batch, aes(x=batch.exp.n, y=log(qvalue.median))) +
           geom_boxplot() +
           geom_hline(yintercept = log(0.01))
plots.list = lappend(plots.list, p)

qvalues.stats <- peptides.peak_sums %>% group_by(EG.StrippedSequence) %>%  dplyr::summarise(qvalue.median = median(EG.Qvalue))
p = ggplot(qvalues.stats, aes(x = log(qvalue.median))) + 
           geom_density() +
           geom_vline(xintercept = log(0.01))
plots.list = lappend(plots.list, p)

#Qvalue is based on median per batch
#peptides.peak_sums = merge(peptides.peak_sums, qvalues.stats.batch, by = c("R.Label", "EG.StrippedSequence", "batch.exp.n") ) 
peptides.peak_sums = merge(peptides.peak_sums, qvalues.stats, by = c("EG.StrippedSequence"))
peptides.peak_sums = tbl_df(peptides.peak_sums)
peptides.peak_sums = filter(peptides.peak_sums, qvalue.median <= 0.01)

## ---- transforming data to normal distribution ----
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

set.seed(123)
peptides.peak_sums.subset = peptides.peak_sums.trimmed[peptides.peak_sums.trimmed$R.Label %in% sample(unique(peptides.peak_sums.trimmed$R.Label),50),]
peptides.peak_sums.subset.norm_stats = ddply(peptides.peak_sums.subset, .(R.Label, batch.exp.n, sample), 
                                             .fun = function(z) {
                                               v<<-z
                                               tmp.data = z$T_signal
                                               #                                         tmp.data = z$T_signal[z$T_signal < quantile(z$T_signal,0.975) & 
                                               #                                                                         z$T_signal > quantile(z$T_signal,0.025)]
                                               tmp.dist = rnorm(n=length(tmp.data), mean=mean(tmp.data), sd=sd(tmp.data))
                                               tmp = ks.test(tmp.data, tmp.dist)
                                               y <- quantile(na.omit(tmp.data), c(0.25, 0.75))
                                               x <- qnorm(c(0.25, 0.75))
                                               
                                               slope =  diff(y)/diff(x)
                                               int  = y[1L] - slope * x[1L]
                                               ret = data.frame(slope =  slope, 
                                                                intercept  = int,
                                                                pval = tmp$p.value,
                                                                x_min = -3,
                                                                y_max = max(tmp.data)
                                               )
                                               return(ret)
                                             })


1 - sum(peptides.peak_sums.subset.norm_stats$pval < (0.05/length(peptides.peak_sums.subset.norm_stats$pval)))/length(peptides.peak_sums.subset.norm_stats$pval)
1 - sum(peptides.peak_sums.subset.norm_stats$pval < 0.05)/length(peptides.peak_sums.subset.norm_stats$pval)

message("Plotting QQ plots")
p = ggplot(peptides.peak_sums.subset) +  
          geom_point(aes(sample=T_signal), alpha=0.1, stat="qq") + 
          geom_abline(data=peptides.peak_sums.subset.norm_stats, aes(intercept=intercept, slope=slope)) +
          geom_text(data=peptides.peak_sums.subset.norm_stats, aes(x=x_min, y=y_max, label=round(pval,2))) +
          facet_wrap(~R.Label, scales="free") +
          theme(axis.title=element_text(size=20))

file_name = "qqplots.T_all_samples.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width=11.7)
plots.list = lappend(plots.list, p)


p = ggplot(peptides.peak_sums.subset) +  
      geom_density(aes(x=T_signal)) + 
      #geom_abline(data=peptides.peak_sums.subset.norm_stats, aes(intercept=intercept, slope=slope)) +
      #geom_text(data=peptides.peak_sums.subset.norm_stats, aes(x=x_min, y=y_max, label=round(pval,2))) +
      facet_wrap(~R.Label, scales="free") +
      theme(axis.title=element_text(size=20), aspect.ratio = 1)

file_name = "density_plots.T_all_samples.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width=11.7)
plots.list = lappend(plots.list, p)





## QC statistics CV
message("Performing QC stats")
peptides.peak_stats = group_by(peptides.peak_sums, batch.exp.n, sample, EG.StrippedSequence) %>%
                              dplyr::summarise(count = n(),
                                               mean.signal = mean(signal),
                                               CV = sd(signal)/mean.signal,
                                               mean.T = mean(T_signal),
                                               CV.T = sd(T_signal)/mean.T)

to_plot = filter(peptides.peak_stats, count>=2)
p = ggplot(to_plot, aes(x=as.numeric(EG.StrippedSequence), y=CV.T)) + 
           geom_point(alpha=0.1) +
           facet_wrap(~sample, scales="free") +
           ggtitle(paste("Grouped by", paste(attr(to_plot, which="vars"), collapse="."))) + 
           stat_smooth()

plots.list = lappend(plots.list, p)
file_name = "CV_all_samples.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width=11.7)


p = ggplot(to_plot, aes(x=as.numeric(EG.StrippedSequence), y=CV)) + 
          geom_point(alpha=0.1) +
          ggtitle(paste("Grouped by", paste(attr(to_plot, which="vars"), collapse="."))) + 
          stat_smooth() +
          theme(text = element_text(size=20))

plots.list = lappend(plots.list, p)
file_name = "CV_all_samples.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width=11.7)


p = ggplot(to_plot, aes(x=as.numeric(EG.StrippedSequence), y=CV.T)) + 
          geom_point(alpha=0.1) +
          ggtitle(paste("Grouped by", paste(attr(to_plot, which="vars"), collapse="."))) + 
          stat_smooth() +
          theme(text = element_text(size=20))
plots.list = lappend(plots.list, p)

file_name = "CV.T_all_samples.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width=11.7)


# ---- using ComBat to correct for Batch effects ----
message("correcting for batch effects")

peptides.df = dcast(data=peptides.peak_sums.trimmed, formula=EG.StrippedSequence~sample+replicate+batch_date+batch.exp.n+batch, value.var="T_signal")
peptides.matrix = as.matrix(peptides.df[,-1])

rownames(peptides.matrix) = peptides.df$EG.StrippedSequence

sp_batches = sub(x=colnames(peptides.matrix), pattern=".*?_([A-Za-z0-9]+)$", perl=T, replacement="\\1")
exp_batches = sub(x=colnames(peptides.matrix), pattern=".*?_([A-Za-z0-9]+)_[A-Za-z0-9]+$", perl=T, replacement="\\1")
batch_date = sub(x=colnames(peptides.matrix), pattern=".*?_([0-9]+_[0-9]+_[0-9]+)_[A-Za-z0-9]+_[A-Za-z0-9]+$", perl=T, replacement="\\1")
sample_name = sub(x=colnames(peptides.matrix), pattern="(.*?)_[0-9]+_[0-9]+_[0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+$", perl=T, replacement="\\1")

pheno = data.frame(exp_batches = exp_batches,
                   sample_name = sample_name,
                   batch_date = batch_date)
rownames(pheno) = colnames(peptides.matrix)
pheno$ORF = droplevels(sample_exp.map$ORF[match(pheno$sample_name, sample_exp.map$SampleName)])

pheno$ORF[pheno$sample_name == "KL_Try_027_c"] = "WT"
pheno$exp_batches[pheno$sample_name == "KL_Try_027_c"] = 5
pheno = droplevels(pheno[which(!is.na(pheno$ORF)),])
tmp.factor_size = ddply(pheno, .(batch_date, ORF), summarise, factor_size = length(ORF))
tmp.factor_size$batch_date.ORF = paste(tmp.factor_size$batch_date, tmp.factor_size$ORF, sep=".")

pheno$batch_date.ORF = paste(pheno$batch_date, pheno$ORF, sep=".")
pheno = droplevels(pheno[pheno$batch_date.ORF %in% tmp.factor_size$batch_date.ORF[tmp.factor_size$factor_size >=2],])
pheno = droplevels(pheno[pheno$batch_date %in% levels(droplevels(pheno$batch_date[grep(pattern="mix", pheno$sample_name, ignore.case=T)])),])

peptides.matrix.f = peptides.matrix[,match(rownames(pheno), colnames(peptides.matrix))]



# saving raw data ##########################
message("Tidying up batch raw data", appendLF=F)
peptides.long = melt(peptides.matrix.f, id.vars=rownames)
names(peptides.long) = c("EG.StrippedSequence", "variable", "value")

peptides.long = peptides.long %>% extract(variable, 
                                          into=c("R.Label", "batch_date", "batch.exp.n", "batch"), 
                                          regex="(.*?)_([0-9]+_[0-9]+_[0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$")
file_name = "peptides.long.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.long, file=file_path)  
message("...Done")
#############################################
########HERE£££££££££££££
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "batch_effects"

load("./R/objects/peptides.peak_sums.trimmed.RData")
load("./R/objects/protein_annotations.RData")
load("./R/objects/peptides.cor.stats.top.RData")

peptides.long = peptides.peak_sums.trimmed
peptides.selected = tbl_df(droplevels(filter(peptides.cor.stats.top, top == "3")))

peptides.long.selected = tbl_df(peptides.long[peptides.long$EG.StrippedSequence %in% peptides.selected$EG.StrippedSequence,])
peptides.long.selected$ORF = peptides.selected$ORF[match(peptides.long.selected$EG.StrippedSequence, peptides.selected$EG.StrippedSequence)]

proteins.long = peptides.long.selected  %>% group_by(ORF, R.Label, batch_date, batch.exp.n, batch) %>% summarise(mean_signal = mean(T_signal, na.rm=T),
                                                                                                       sum_signal = sum(T_signal, na.rm=T))



proteins.df = dcast(data=proteins.long, formula=ORF~R.Label+batch_date+batch.exp.n+batch, value.var="mean_signal")

proteins.matrix = proteins.df[,-1]
rownames(proteins.matrix) = proteins.df$ORF


pattern.p = "(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$"
matches = stringr::str_match_all(pattern=pattern.p, colnames(proteins.matrix))

stopifnot(sum(lapply(matches,length)!=0) == ncol(proteins.matrix))

pheno = data.frame(matrix(unlist(matches), ncol=length(matches[[1]]), byrow=T))
colnames(pheno) = c("name", "R.Label", "batch_date", "batch.exp.n", "batch" )
rownames(pheno) = colnames(proteins.matrix)

pheno$ORF = droplevels(sample_exp.map$ORF[match(pheno$R.Label, sample_exp.map$SampleName)])

pheno$ORF[pheno$R.Label == "KL_Try_027_c"] = "WT"
pheno$batch.exp.n[pheno$R.Label == "KL_Try_027_c"] = 5

pheno = droplevels(pheno[which(!is.na(pheno$ORF)),])



pheno$group = pheno$batch.exp.n #grouping variable to estimate batch effects

tmp.factor_size = ddply(pheno, .(group, ORF), summarise, factor_size = length(ORF))
tmp.factor_size$group.ORF = paste(tmp.factor_size$group, tmp.factor_size$ORF, sep=".")

pheno$group.ORF = paste(pheno$group, pheno$ORF, sep=".")
pheno = droplevels(pheno[pheno$group.ORF %in% tmp.factor_size$group.ORF[tmp.factor_size$factor_size >=2],])
pheno = droplevels(pheno[pheno$group %in% levels(droplevels(pheno$group[grep(pattern="mix", pheno$R.Label, ignore.case=T)])),])

proteins.matrix.f = proteins.matrix[,match(rownames(pheno), colnames(proteins.matrix))]
mod = model.matrix(~as.factor(ORF), data=pheno)

tmp.size_factors = DESeq::estimateSizeFactorsForMatrix(exp(proteins.matrix.f))
proteins.matrix.f.deseq = log(exp(proteins.matrix.f)/tmp.size_factors, base=2)
#proteins.matrix.f.deseq = log(exp(proteins.matrix.f), base=2)

proteins.matrix.f.combat = ComBat(na.omit(proteins.matrix.f), batch=pheno$group, mod=mod)
proteins.matrix.f.deseq.combat = ComBat(na.omit(proteins.matrix.f.deseq), batch=pheno$group, mod=mod)

file_name = "proteins.matrix.f.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.matrix.f.combat,file=file_path)  

file_name = "proteins.matrix.f.deseq.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.matrix.f.deseq.combat,file=file_path)  



# pca
message("plotting PCA results")
# pca with quantile normalization

file_name = "PCA_batch_effects.png"
file_path = paste(figures_dir, file_name, sep="/")
png(file_path, width=297, height=210, units="mm", res=150)
par(pty="s", mfrow=c(1,3), cex=0.75)

pca = prcomp(t(proteins.matrix.f[complete.cases(proteins.matrix.f),]), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], col=pheno$batch_date, pch=16, main="Before adjustments for batch effects",
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=pheno$group, cex=0.66)

pca = prcomp(t(proteins.matrix.f.combat), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], col=pheno$batch_date, pch=16, main="After adjustments for batch effects",
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=pheno$group, cex=0.66)


pca = prcomp(t(proteins.matrix.f.deseq.combat), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], col=pheno$batch_date, pch=16, main="After normalization and adjustments for batch effects",
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=pheno$group, cex=0.66)

p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()

# dendrograms 
file_name = "Clustering_batch_effects.png"
file_path = paste(figures_dir, file_name, sep="/")
png(file_path, width=297, height=210, units="mm", res=150)
resetPar()
par(mfrow=c(2,1))
h = hclust(dist(scale(t(proteins.matrix.f))))
plot(h, labels=pheno$group, cex=0.5, main="Before batch correction")

h = hclust(dist(scale(t(proteins.matrix.f.deseq.combat))))
plot(h, labels=pheno$group, cex=0.5, main="After batch correction")

p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()



#tidying batch corrected data
message("Tidying up batch corrected data", appendLF=F)

tmp.wide = proteins.matrix.f.combat
tmp.wide$ORF = rownames(proteins.matrix.f.combat)
proteins.deseq.combat.long = melt(tmp.wide, id.vars="ORF")
names(proteins.deseq.combat.long) = c("ORF", "variable", "value")

proteins.deseq.combat.long = proteins.deseq.combat.long %>% extract(variable, 
                                                                   into=c("R.Label", "batch_date", "batch.exp.n", "batch"), 
                                                                   regex="(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$")

col_names <- names(proteins.deseq.combat.long)[names(proteins.deseq.combat.long) != "value"]
proteins.deseq.combat.long[,col_names] <- lapply(proteins.deseq.combat.long[,col_names] , factor)
proteins.deseq.combat.long = tbl_df(proteins.deseq.combat.long)

file_name = "proteins.deseq.combat.long.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.deseq.combat.long, file=file_path)  
message("...Done")



tmp.wide = proteins.matrix.f.combat
tmp.wide$ORF = rownames(proteins.matrix.f.combat)
proteins.combat.long = melt(tmp.wide, id.vars="ORF")
names(proteins.combat.long) = c("ORF", "variable", "value")

proteins.combat.long = proteins.combat.long %>% extract(variable, 
                                                          into=c("R.Label", "batch_date", "batch.exp.n", "batch"), 
                                                          regex="(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$")

col_names <- names(proteins.combat.long)[names(proteins.combat.long) != "value"]
proteins.combat.long[,col_names] <- lapply(proteins.combat.long[,col_names] , factor)
proteins.combat.long = tbl_df(proteins.combat.long)

file_name = "proteins.combat.long.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.combat.long, file=file_path)  
message("...Done")


tmp.wide = proteins.matrix.f
tmp.wide$ORF = rownames(proteins.matrix.f)
proteins.long = melt(tmp.wide, id.vars="ORF")
names(proteins.long) = c("ORF", "variable", "value")

proteins.long = proteins.long %>% extract(variable, 
                                          into=c("R.Label", "batch_date", "batch.exp.n", "batch"), 
                                          regex="(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$")

col_names <- names(proteins.long)[names(proteins.long) != "value"]
proteins.long[,col_names] <- lapply(proteins.long[,col_names] , factor)
proteins.long = tbl_df(proteins.long)

file_name = "proteins.long.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.long, file=file_path)  
message("...Done")



set.seed(123)
toSelect = droplevels(sample(size=20 , x=unique(proteins.deseq.combat.long$R.Label)))
toPlot = proteins.deseq.combat.long[proteins.deseq.combat.long$R.Label %in% toSelect,]

p = ggplot(toPlot, aes(x=value)) +
       geom_histogram(aes(y=..density..), 
                      breaks=seq(1,10, 0.2),
                      colour="black", 
                      fill="white") +
       stat_function(fun=dnorm, args=list(mean=mean(toPlot$value), sd=sd(toPlot$value)))+
       facet_wrap(~R.Label, scales="free")


p1 = ggplot(proteins.long, aes(x=batch.exp.n, y=value, fill=batch.exp.n)) +
            geom_boxplot() +
            theme(legend.position="none")

p3 = ggplot(proteins.deseq.combat.long, aes(x=batch.exp.n, y=value, fill=batch.exp.n)) +
            geom_boxplot() + 
            theme(legend.position="top")

g = arrangeGrob(p1, p3)
plots.list = lappend(plots.list, g)

