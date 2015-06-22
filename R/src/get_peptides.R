#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "get_peptides"

## ---- data_load ----
load("./R/objects/peptides.data.RData") # clean.R::prepare_peptides()
load("./R/objects/experiment.map._load_.RData")# load.R::load_batch_map()
load("./R/objects/dates_map._load_.RData") # load.R::load_dates_map()
load("./R/objects/sample.map._load_.RData") # load.R::load_sample_map()
load("./R/objects/protein_annotations._load_.RData")


## ---- selecting peptides based on spectronaut Q-value
peptides.data$EG.StrippedSequence = factor(peptides.data$EG.StrippedSequence)
peptides.data$R.Label = factor(peptides.data$R.Label)

peptides.peak_sums <- group_by(peptides.data, batch_date, batch, batch.exp.n, R.Label, EG.StrippedSequence) %>%
  dplyr::summarise(count = n(),
                   signal = FG.TotalPeakArea[1],
                   EG.Qvalue = EG.Qvalue[1]  )


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

peptides.peak_sums.f = filter(peptides.peak_sums, qvalue.median <= 0.01)
peptides.peak_sums.stats = peptides.peak_sums %>% group_by(R.Label, batch.exp.n) %>% summarize(sum = sum(signal),
                                                                                               median = median(signal),
                                                                                               shorth = genefilter::shorth(signal))
peptides.peak_sums.stats$category = "raw"
  
peptides.peak_sums.f.stats = peptides.peak_sums.f %>% group_by(R.Label, batch.exp.n) %>% summarize(sum = sum(signal),
                                                                                                   median = median(signal),
                                                                                                   shorth = genefilter::shorth(signal))
peptides.peak_sums.f.stats$category = "filtered"


peptides.peak_sums.f.stats.mix = peptides.peak_sums.f.stats[grep(x=peptides.peak_sums.f.stats$R.Label, pattern="mix", ignore.case=T),]
peptides.peak_sums.f.stats.wt  = peptides.peak_sums.f.stats[peptides.peak_sums.f.stats$R.Label %in% sample_map$SampleName[sample_map$ORF == "WT"],]

toPlot = rbind(peptides.peak_sums.stats,peptides.peak_sums.f.stats)
toPlot$category = factor(toPlot$category)
p = ggplot(toPlot, aes(x=batch.exp.n, y=shorth)) +
      geom_boxplot(aes(fill =category)) +
      geom_point(data=peptides.peak_sums.f.stats.mix,aes(x=jitter(as.numeric(batch.exp.n)), y=shorth),col="red") +
      geom_point(data=peptides.peak_sums.f.stats.wt,aes(x=jitter(as.numeric(batch.exp.n)), y=shorth),col="blue")
plots.list = lappend(plots.list, p)




## ---- transforming data to normal distribution ----
peptides.peak_sums.f$T_signal = with(peptides.peak_sums, log(signal))

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
                   CV.T = sqrt(exp(sd(T_signal)^2)-1))

to_plot = filter(peptides.peak_stats, count>=2)
p = ggplot(to_plot, aes(x=as.numeric(EG.StrippedSequence), y=CV.T)) + 
  geom_point(alpha=0.1) +
  facet_wrap(~sample, scales="free") +
  ggtitle(paste("Grouped by", paste(attr(to_plot, which="vars"), collapse="."))) + 
  stat_smooth()

plots.list = lappend(plots.list, p)
file_name = "CV_all_samples_log.png"
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
  ylim(0,2)+
  theme(text = element_text(size=20))

plots.list = lappend(plots.list, p)

file_name = "CV.T_all_samples.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p, height=8.27, width=11.7)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 

if (FALSE) {
  
  ## ---- peptides to proteins Mapping statistics---- 
  file_name = "mapping_stats_histogram.png"
  file_path = paste(figures_dir, file_name, sep="/")
  png(file_path, width=297, height=210, units="mm", res=150)
  par(pty="s", mfrow=c(1,2))
  
  toSelect = names(table(peptide2orfs$EG.StrippedSequence)[!table(peptide2orfs$EG.StrippedSequence)>1])
  nr_peptides_all = length(unique(peptide2orfs$EG.StrippedSequence))
  nr_proteins_all = length(unique(peptide2orfs$ORF))
  nr_peptides1 = length(unique(peptide2orfs$EG.StrippedSequence)) - sum(table(peptide2orfs$EG.StrippedSequence) > 1)
  nr_proteins1 = length(levels(droplevels(peptide2orfs[peptide2orfs$EG.StrippedSequence %in% toSelect, "ORF"])))
  
  
  hist(table(peptide2orfs$EG.StrippedSequence), xlab="Number of peptides per ORF", main="")
  legend("topright", 
         legend = c(
           paste("Number of peptides total:", nr_peptides_all ),
           paste("Number of proteins total:", nr_proteins_all ),
           paste("Number of peptides unique mapping:", nr_peptides1 ),
           paste("Number of proteins with unique mapping:", nr_proteins1)))
  
  toSelect = names(table(peptide2orfs$EG.StrippedSequence)[!table(peptide2orfs$EG.StrippedSequence)>1])
  peptide2orfs.lt1 = droplevels(peptide2orfs[peptide2orfs$EG.StrippedSequence %in% toSelect, ])
  
  nr_ORF_all = length(unique(peptide2orfs$ORF))
  fr_ORF_gt1 = round(sum(table(peptide2orfs$ORF) > 1) /length(unique(peptide2orfs$ORF)),2)
  fr_ORF_1 = round(sum(table(peptide2orfs.lt1$ORF) > 1)/length(unique(peptide2orfs$ORF)),2)
  
  hist(table(peptide2orfs$ORF), breaks=50, xlab="Number of peptides mapped per ORF", main="")
  legend("topright", 
         legend = c(
           paste("Number of proteins total:", nr_ORF_all ),
           paste("Fraction of proteins >1 peptide:", fr_ORF_gt1 ),
           paste("Fraction of proteins >1 peptide\nafter removal non-unique peptides", fr_ORF_1)))
  
  p = recordPlot()
  plots.list = lappend(plots.list, p)
  dev.off()
  
  ## ---- dataset_Stats ----
  data.stats = data.frame(Kinases =  length(levels(droplevels(sample_map[sample_map$Type == "Kinase", "ORF"]))),
                          Metabolic = length(levels(droplevels(sample_map[sample_map$Type == "Metabolic", "ORF"]))),
                          Standard_mix = length(levels(droplevels(sample_map[sample_map$Type == "Standard Mix", "ORF"]))),
                          ZWF1_interaction = length(levels(droplevels(sample_map[sample_map$Type == "ZWF1 interaction", "ORF"]))),
                          n_samples = length(levels(sample_map$SampleName)))
  
  
  library("gridExtra")
  p_table = tableGrob(data.stats)
  
  file_name = paste("mutant_Stats", "table.pdf", sep=".")
  file_path = paste(figures_dir, file_name, sep="/")
  pdf(file_path, paper="a4")
  grid.arrange(p_table)
  
  dev.off()
}



