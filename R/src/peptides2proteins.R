#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "explore"

## ---- data_load ----
load("./R/objects/peptides.data.RData") # clean.R::prepare_peptides()
load("./R/objects/experiment.map.RData")# load.R::load_batch_map()
load("./R/objects/dates_map.RData") # load.R::load_dates_map()
load("./R/objects/sample_exp.map.RData") # load.R::load_sample_map()
load("./R/objects/peptides.deseq.combat.long.RData")
load("./R/objects/protein_annotations.RData")


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


# ---- tyding data, preparing for batch effect corrections ----
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

file_name = "pheno.RData"
file_path = paste(output_dir, file_name, sep="/")
save(pheno, file=file_path)  

message("...Done")




## ---- checking for correlations of peptides per protein ---- 
peptides.dataset = tbl_df(peptides.long)
protein_annot.selected = droplevels(filter(protein_annotations, specificity == "unique"))

peptide2orfs = select(protein_annot.selected, strippedSequence, SystName)
names(peptide2orfs) = c("EG.StrippedSequence", "ORF")

peptides.dataset.f = droplevels(peptides.dataset[peptides.dataset$EG.StrippedSequence %in% protein_annot.selected$strippedSequence,])

peak_sums = peptides.dataset.f
peak_sums.f = tbl_df(droplevels(merge(peak_sums, peptide2orfs, by="EG.StrippedSequence")))
peak_sums.f = peak_sums.f %>% distinct(R.Label, EG.StrippedSequence) %>% arrange(R.Label, ORF, EG.StrippedSequence)

peptides.wide = dcast(data=peak_sums.f, formula=ORF+EG.StrippedSequence~R.Label+batch.exp.n+batch, value.var="value")
peptides.wide.selected = dcast(data=peak_sums.f, formula=ORF+EG.StrippedSequence+batch.exp.n~R.Label+batch, value.var="value")


#peptides.long %>% group_by(ORF, variable, value, R.Label, batch.exp.n, batch_date, batch) %>% dplyr::summarise(sum(value))

batch.lengths = peak_sums.f %>% group_by(batch.exp.n) %>% summarize(counts = length(unique(R.Label))) %>% arrange(-counts)

#selecting top by size of batch
peptides.wide.selected = droplevels(peptides.wide.selected[peptides.wide.selected$batch.exp.n %in% batch.lengths$batch.exp.n[1:5],]) 


peptides.sel.cor = ddply(peptides.wide.selected, .(ORF, batch.exp.n), 
                        .fun = function(z) {
                          v <<- z
                          rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                          if (nrow(z) == 1) {
                            return(data.frame(id = z$EG.StrippedSequence, variable=row.names(z), value=NA))
                          }
                          
                          rownames(z) = paste(z$EG.StrippedSequence, sep=".")  
                          tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                          tmp.matrix[!upper.tri(tmp.matrix)] = NA
                          
                          tmp.df = data.frame(tmp.matrix)
                          tmp.df$id = rownames(tmp.df)
                          tmp.long = melt(tmp.df, id.vars=c("id"))
                          tmp.long = tmp.long[complete.cases(tmp.long),]
                          return(tmp.long)
                        })

peptides.sel.cor[sapply(peptides.sel.cor, is.character)] <- lapply(peptides.sel.cor[sapply(peptides.sel.cor, is.character)],  as.factor)

p = ggplot(peptides.sel.cor, aes(x=value)) +
           geom_density() +
           facet_wrap(~batch.exp.n)

plots.list = lappend(plots.list, p)


peptides.correlations = ddply(peptides.wide, .(ORF), 
                              .fun = function(z) {
                                v <<- z                      
                                rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                if (nrow(z) == 1) {
                                  return(data.frame(id = z$EG.StrippedSequence, variable=row.names(z), value=NA))
                                }
                                rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                tmp.matrix = cor(t(z[,-c(1:2)]), method="spearman", use="pairwise.complete.obs")
                                tmp.matrix[!upper.tri(tmp.matrix)] = NA
                                
                                tmp.df = data.frame(tmp.matrix)
                                tmp.df$id = rownames(tmp.df)
                                tmp.long = melt(tmp.df, id.vars=c("id"))
                                tmp.long = tmp.long[complete.cases(tmp.long),]
                                return(tmp.long)
                              })


file_name = "correlation_histogram.png"
file_path = paste(figures_dir, file_name, sep="/")
png(file_path, width=210, height=297, units="mm", res=150)

par(pty="s")
hist(peptides.correlations$value, breaks=50, xlab="Spearman's correaltion, rho", main="", cex.axis=1.2, cex.lab=2)
legend("topleft", bg=NULL, bty="n", 
       legend=paste("Mean correlation ", round(mean(peptides.correlations$value, na.rm=T),2)))
dev.off()


peptides.correlations.subset = droplevels(peptides.correlations[(peptides.correlations$value > 0.5 | peptides.correlations$value < 0.5) & !is.na(peptides.correlations$value), ])
tmp1 = tbl_df(melt(peptides.wide[match(peptides.correlations.subset$id, peptides.wide$EG.StrippedSequence),], id.vars=c("ORF", "EG.StrippedSequence")))
tmp2 = tbl_df(melt(peptides.wide[match(peptides.correlations.subset$variable, peptides.wide$EG.StrippedSequence),], id.vars=c("ORF", "EG.StrippedSequence")))

peptides.correlations.long = tbl_df(cbind(tmp1,tmp2[,c("EG.StrippedSequence","value")]))
names(peptides.correlations.long)[2] = "EG.StrippedSequence1"
names(peptides.correlations.long)[5] = "EG.StrippedSequence2"
names(peptides.correlations.long)[4] = "value1"
names(peptides.correlations.long)[6] = "value2"

peptides.correlations.long$peptides = factor(paste(peptides.correlations.long$EG.StrippedSequence1, peptides.correlations.long$EG.StrippedSequence2, sep="."))
peptides.correlations.long = peptides.correlations.long %>% 
    group_by(ORF, EG.StrippedSequence1, EG.StrippedSequence2) %>% 
    arrange()

peptides.correlations.long = peptides.correlations.long %>% 
  extract(variable, into=c("R.Label", "batch.exp.n", "batch"), regex="(.*?)_([0-9]+)_([0-9]+)$")


# selecting top lowest correlated peptides
toSelect = droplevels(unique(peptides.correlations.subset[with(peptides.correlations.subset, order(value))[1:20],c("id","variable")])) 
toSelect$peptides = with(toSelect, paste(id, variable, sep="."))

toPlot = droplevels(peptides.correlations.long[peptides.correlations.long$peptides %in% toSelect$peptides,])
p = ggplot(toPlot, aes(x=value1, y=value2, color=batch.exp.n)) + 
  geom_point() +
  ggtitle("top 20 low correlation peptides") + 
  facet_wrap(~peptides, scales="free") + theme(aspect.ratio=1) +
  geom_smooth(method="lm", fill=NA)

file_name = "top20_lowCorrelation.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, width=11.69+0.1*11.69, height=8.27+0.1*8.27)
plots.list = lappend(plots.list, p)


# selecting top highest correlated peptides
toSelect = droplevels(unique(peptides.correlations.subset[with(peptides.correlations.subset, order(-value))[1:20],c("id","variable")])) 
toSelect$peptides = with(toSelect, paste(id, variable, sep="."))

toPlot = droplevels(peptides.correlations.long[peptides.correlations.long$peptides %in% toSelect$peptides,])
p = ggplot(toPlot, aes(x=value1, y=value2, color=batch.exp.n)) + 
  geom_point() +
  ggtitle("top 20 high correlation peptides") + 
  facet_wrap(~peptides, scales="free") +
  geom_smooth(method="lm", fill=NA)

file_name = "top20_highCorrelation.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, width=11.69+0.1*11.69, height=8.27+0.1*8.27)
plots.list = lappend(plots.list, p)


#ranking peptides and additional stats

peak_sums.f[sapply(peak_sums.f, is.character)] <- lapply(peak_sums.f[sapply(peak_sums.f, is.character)],  as.factor)

peptides.ranks = peak_sums.f %>% group_by(R.Label, ORF) %>% 
                                    mutate(rev.rank = rank(-value)) %>% # peptide rank per ORF
                                 group_by(R.Label) %>%
                                    mutate(chunk = cut(value, 5, labels=c(1:5))) %>% # devide sample signals into chunks 
                                 group_by(R.Label, ORF) %>%
                                    mutate(peptide.count = n()) # how many peptides per protein

peptides.ranks.summary =  peptides.ranks %>% group_by(ORF, EG.StrippedSequence) %>%  
                                                summarize(mean.rank = mean(rev.rank),
                                                          mean.chunk= mean(chunk),
                                                          peptide.count = unique(peptide.count)) 

#checking if there is a strong effect if batches per peptide pair 
#apparently there is a strong effect
batch.cor.var = peptides.sel.cor %>% group_by(id,variable) %>% summarize(CV = sd(abs(value), na.rm=T)/mean(abs(value), na.rm=T))
#hist(batch.cor.var$CV, xlim=c(0,2))


peptides.cor.stats = ddply(peptides.wide, .(ORF), 
                            .fun = function(z) {
                              v <<- z                      
                              z = v
                              rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                              if (nrow(z) == 1) {
                                return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                              }
                              rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                              tmp.matrix = cor(t(z[,-c(1:2)]), method="spearman", use="pairwise.complete.obs")
                              
                              cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                              cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                              cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                              
                              ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                              return(ret)
                            })


p = ggplot(melt(peptides.cor.stats, id.vars = c("ORF", "EG.StrippedSequence")), aes(x=value, color=variable)) +
           geom_density() +
           xlab("Expected (average,median,shorth) Spearman's correlation coefficient between peptides") +
           theme(aspect.ratio = 1)
plots.list = lappend(plots.list, p)

peptides.ranks.summary = merge(peptides.ranks.summary, peptides.cor.stats[,c("ORF", "EG.StrippedSequence", "cor.means")], by = c("ORF", "EG.StrippedSequence"))
peptides.ranks.summary = peptides.ranks.summary %>% group_by(ORF) %>% mutate(topN = rank(mean.rank))


## ---- correlations of peptides considering their ranks ---- 
peptides.wide.topN = merge(peptides.wide, peptides.ranks.summary[,c("ORF", "EG.StrippedSequence", "topN")], by=c("ORF", "EG.StrippedSequence"))
peptides.wide.topN = peptides.wide.topN[,c(colnames(peptides.wide.topN)[ncol(peptides.wide.topN)], colnames(peptides.wide.topN)[-ncol(peptides.wide.topN)])]

peptides.cor.stats.top5 = ddply(filter(peptides.wide.topN, topN <=5), .(ORF), 
                           .fun = function(z) {
                             v <<- z                      
                             z = v
                             rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                             if (nrow(z) == 1) {
                               return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                             }
                             rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                             tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                             
                             cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                             cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                             cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                             
                             ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                             return(ret)
                           })


peptides.cor.stats.gt_top5 = ddply(filter(peptides.wide.topN, topN > 5), .(ORF), 
                                .fun = function(z) {
                                  v <<- z                      
                                  z = v
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  if (nrow(z) == 1) {
                                    return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                                  }
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                                  
                                  cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                                  cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                                  cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                                  
                                  ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                                  return(ret)
                                })


peptides.cor.stats.top3 = ddply(filter(peptides.wide.topN, topN <=3), .(ORF), 
                                .fun = function(z) {
                                  v <<- z                      
                                  z = v
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  if (nrow(z) == 1) {
                                    return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                                  }
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                                  
                                  cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                                  cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                                  cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                                  
                                  ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                                  return(ret)
                                })

hist(peptides.cor.stats.top3$cor.means)
t.test(peptides.cor.stats.top3$cor.means, peptides.cor.stats.top5$cor.means)

peptides.ranks.summary = merge(peptides.ranks.summary, peptides.cor.stats[,c("ORF", "EG.StrippedSequence", "cor.means")], by = c("ORF", "EG.StrippedSequence"))

peptides.cor.stats.top3$top = 3
peptides.cor.stats.top5$top = 5
peptides.cor.stats.gt_top5$top = "gt_5"
tmp = peptides.cor.stats
tmp$top = "all"

peptides.cor.stats.top = rbind(peptides.cor.stats.top3, peptides.cor.stats.top5, peptides.cor.stats.gt_top5, tmp)
peptides.cor.stats.top.long = melt(peptides.cor.stats.top, id.vars=c("ORF", "EG.StrippedSequence", "top"))

p = ggplot(peptides.cor.stats.top.long, aes(x=value, fill=top)) +
           geom_density(alpha=0.5) +
           facet_wrap(~variable) +
           xlim("Spearman's rank correlation, rho")

plots.list = lappend(plots.list, p)


