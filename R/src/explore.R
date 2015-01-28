#!/usr/bin/env Rscript
rm(list = ls())

source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "explore"

## ---- data_load ----
load("./R/objects/peptides.data.RData") # clean.R::prepare_peptides()
load("./R/objects/experiment.map.RData")# load.R::load_batch_map()
load("./R/objects/dates_map.RData") # load.R::load_dates_map()
load("./R/objects/sample.map.RData") # load.R::load_sample_map()
#View(peptides.data)

#peptides.data.mix = droplevels(peptides.data[grep(pattern="mix", ignore.case=T, peptides.data$R.Label),])


##batch/date signal statistics
grouped.data <- group_by(peptides.data, batch, batch_date, R.Label)
grouped_stats <- dplyr::summarise(grouped.data,
                                  sum_FG.TotalPeakArea = sum(FG.TotalPeakArea, na.rm=T),
                                  mean_FG.TotalPeakArea = mean(FG.TotalPeakArea, na.rm=T))

grouped_stats <- dplyr::summarise(grouped.data,
                                  sum_FG.TotalPeakArea = sum(FG.TotalPeakArea, na.rm=T),
                                  mean_FG.TotalPeakArea = mean(FG.TotalPeakArea, na.rm=T))

p = ggplot(grouped_stats, aes(x=batch, y = mean_FG.TotalPeakArea, colour=batch_date)) + 
  geom_boxplot() + 
  theme(aspect.ratio = 1) +
  ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))

plots.list = lappend(plots.list, p)

grouped_stats$sample_number = 1:nrow(grouped_stats)
p = ggplot(grouped_stats, aes(x=sample_number, y = mean_FG.TotalPeakArea, col=batch_date)) + 
  geom_point() + 
  theme(aspect.ratio = 1) +
  ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))

plots.list = lappend(plots.list, p)

p = ggplot(grouped_stats, aes(x=sample_number, y = mean_FG.TotalPeakArea, col=batch)) + 
  geom_boxplot() +
  geom_point() + 
  theme(aspect.ratio = 1) +
  ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))

plots.list = lappend(plots.list, p)

#PCA of peptides and batch effects

grouped.data <- group_by(peptides.data, batch, batch_date, sample, EG.StrippedSequence)
grouped_stats <- dplyr::summarise(grouped.data,
                                  count = n(),
                                  sum_FG.PeakArea = sum(F.PeakArea, na.rm=T))

peptides.df = dcast(data=grouped_stats, formula=EG.StrippedSequence~sample, value.var="sum_FG.PeakArea", fun.aggregate=mean)


peptides.matrix = as.matrix(peptides.df[,-1])
rownames(peptides.matrix) = peptides.df$EG.StrippedSequence
experiment_map$sample = factor(sub(x=experiment_map$SampleName, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
experiment_map$replicate = factor(sub(x=experiment_map$SampleName, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
batches = as.numeric(factor(experiment_map[match(colnames(peptides.matrix), experiment_map$sample), "batch"]))

#plotting
# ---- dendrogram ----
par(pty="s")
h = hclust(dist(scale(t(peptides.matrix))))
plot(h, labels=batches, cex=0.66 )
p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()

# ---- pca ----
par(pty="s")
pca = prcomp(t(peptides.matrix), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], col=batches, pch=16,
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=batches, cex=0.66)
p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()


# peptides.matrix = peptides.matrix[,which(!is.na(match(colnames(peptides.matrix), sample_map$SampleName)))]
# peptides.matrix = peptides.matrix[, which(experiment_map[match(colnames(peptides.matrix), experiment_map$SampleName), "batch"] != "")]


#   grouped_stats$sample_number = as.numeric(grouped_stats$R.Label)
#   
# length(unique(grouped_stats$sample_number))
# p = ggplot(grouped_stats, aes(x=log(sum_FG.PeakArea,base=10), colour=sample_number, group=sample_number)) + 
#   geom_density(fill=NA) +
#   theme(aspect.ratio = 1) +
#   ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
# p
# plots.list = lappend(plots.list, p)
# 
# p = ggplot(grouped_stats, aes(x=log(sum_FG.PeakArea,base=10), colour=date, group=sample_number)) + 
#   geom_density(fill=NA) +
#   theme(aspect.ratio = 1) +
#   ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
# plots.list = lappend(plots.list, p)

####################################
#plotting standart mix distributions
####################################
# grouped.data <- group_by(peptides.data.mix, batch, date, R.Label, EG.StrippedSequence)
# grouped_stats <- dplyr::summarise(grouped.data,
#                                   count = n(),
#                                   sum_FG.PeakArea = sum(F.PeakArea, na.rm=T))
# grouped_stats$sample_number = as.numeric(grouped_stats$R.Label)
# 
# p = ggplot(grouped_stats, aes(y=log(sum_FG.PeakArea,base=10), colour=date, x=sample_number)) + 
#   geom_boxplot(fill=NA) +
#   theme(aspect.ratio = 1) +
#   ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
# 
# plots.list = lappend(plots.list, p)

#   tmp = dcast(grouped_stats, formula=batch+R.Label~EG.StrippedSequence, value.var="sum_FG.PeakArea")
#   tmp = tmp[!duplicated(tmp$R.Label),]
#   
#   tmp.matrix = as.matrix(tmp[,-c(1,2)])
#   rownames(tmp.matrix) = with(tmp, paste(batch,R.Label,sep="."))  
#   tmp.matrix = t(tmp.matrix)
#   norm.matrix = tmp.matrix/estimateSizeFactorsForMatrix(tmp.matrix)
#   tmp.matrix = norm.matrix
#   cor.tmp.matrix = cor(tmp.matrix)

#   plot(cor.tmp.matrix[upper.tri(cor.tmp.matrix, diag=F)])


file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path) 