#!/usr/bin/env Rscript

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

peptides.data$batch.exp.n = factor(as.numeric(peptides.data$batch.exp))

##batch/date signal statistics
grouped.data <- group_by(peptides.data, batch.exp.n, batch_date, R.Label)
grouped_stats <- dplyr::summarise(grouped.data,
                                  sum_FG.TotalPeakArea = sum(FG.TotalPeakArea, na.rm=T),
                                  mean_FG.TotalPeakArea = mean(FG.TotalPeakArea, na.rm=T))

grouped_stats <- dplyr::summarise(grouped.data,
                                  sum_FG.TotalPeakArea = sum(FG.TotalPeakArea, na.rm=T),
                                  mean_FG.TotalPeakArea = mean(FG.TotalPeakArea, na.rm=T))

p = ggplot(grouped_stats, aes(x=batch.exp.n, y = mean_FG.TotalPeakArea, colour=batch_date)) + 
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

p = ggplot(grouped_stats, aes(x=sample_number, y = mean_FG.TotalPeakArea, col=batch.exp.n)) + 
  geom_boxplot() +
  geom_point() + 
  theme(aspect.ratio = 1) +
  ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))

plots.list = lappend(plots.list, p)

#PCA of peptides and batch effects

grouped.data <- group_by(peptides.data, batch.exp.n, batch, batch_date, sample, replicate, EG.StrippedSequence)
grouped_stats <- dplyr::summarise(grouped.data,
                                  count = n(),
                                  sum_FG.PeakArea = sum(F.PeakArea, na.rm=T))

peptides.df = dcast(data=grouped_stats, formula=EG.StrippedSequence~sample+replicate+batch.exp.n+batch, value.var="sum_FG.PeakArea")


peptides.matrix = as.matrix(peptides.df[,-1])
rownames(peptides.matrix) = peptides.df$EG.StrippedSequence

# experiment_map$sample = factor(sub(x=experiment_map$SampleName, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
# experiment_map$replicate = factor(sub(x=experiment_map$SampleName, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))

batches = sub(x=colnames(peptides.matrix), pattern=".*?_([A-Za-z0-9]+)$", perl=T, replacement="\\1")
sp_batches = sub(x=colnames(peptides.matrix), pattern=".*?_([A-Za-z0-9]+)$", perl=T, replacement="\\1")
exp_batches = sub(x=colnames(peptides.matrix), pattern=".*?_([A-Za-z0-9]+)_[A-Za-z0-9]+$", perl=T, replacement="\\1")

pheno = data.frame(batches = batches)
rownames(pheno) = colnames(peptide.matrix)

# samples = sub(x=colnames(peptides.matrix), pattern="(.*?)_[A-Za-z0-9]+_[A-Za-z0-9]+$", perl=T, replacement="\\1")
# samples = as.numeric(factor(samples))

#plotting
# ---- dendrogram ----
par(pty="s")
h = hclust(dist(scale(t(peptides.matrix))))
plot(h, labels=batches, cex=0.66, main= "uncorrected")
p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()

# ---- pca ----
par(pty="s")
pca = prcomp(t(peptides.matrix), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)


plot(pca$x[,1], pca$x[,2], col=sp_batches, pch=16,
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=paste(sp_batches,exp_batches,sep="."), cex=0.66)
p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()


pca.batches = data.frame(pc1 = pca$x[,1],
                         pc2 = pca$x[,2],
                         sp_batches = sp_batches,
                         exp_batches = exp_batches)


# ---- pca MIX ----
peptides.matrix.mix = peptides.matrix[,grep(pattern="mix", x=colnames(peptides.matrix), ignore.case=T)]
batches.mix = sub(x=colnames(peptides.matrix.mix), pattern=".*?_([A-Za-z0-9]+)$", perl=T, replacement="\\1")

pca = prcomp(t(peptides.matrix.mix), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], col=batches, pch=16,
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=batches.mix, cex=0.66)




# ###################################
# #plotting standart mix distributions
# ###################################
# 
# 
# peptides.data.mix = droplevels(peptides.data[grep(pattern="mix", x=peptides.data$sample, ignore.case=T),])
# peptides.data.mix$R.Label = factor(peptides.data.mix$R.Label)
# 
# grouped.data <- group_by(peptides.data.mix, batch_date, sample, replicate,batch.exp.n, R.Label, EG.StrippedSequence)
# grouped_stats <- dplyr::summarise(grouped.data,
#                                   count = n(),
#                                   sum_FG.PeakArea = sum(F.PeakArea, na.rm=T))
# 
# grouped_stats$sample_number = as.numeric(grouped_stats$R.Label)
# 
# p = ggplot(grouped_stats, aes(y=log(sum_FG.PeakArea,base=10), colour=batch_date, x=batch.exp.n)) + 
#           geom_boxplot() +
#           theme(aspect.ratio = 1) +
#           ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
# 
# plots.list = lappend(plots.list, p)
# 
# ## ---- linear batch fits ---- 
# str(grouped_stats$batch.exp.n)
# peptides.mix.grouped_stats = grouped_stats
# batch.fits = ddply(peptides.mix.grouped_stats, .(EG.StrippedSequence),
#                   .fun = function(z) {
#                     v <<- z
#                     fit = lm(sum_FG.PeakArea~batch.exp.n, data=z)
#                     ret = data.frame(batch.exp.n = fit$xlevels,
#                                      batch_factor = names(fit$coefficients),
#                                      coefficients = fit$coefficients,
#                                      EG.StrippedSequence = unique(z$EG.StrippedSequence))
#                     return(ret)
#                     
#                   })
# 
# 
# 
# peptides.grouped.mix.batches = merge(peptides.mix.grouped_stats, batch.fits, by=c("batch.exp.n", "EG.StrippedSequence"), all=T)
# peptides.grouped.mix.batches$corrected = with(peptides.grouped.mix.batches, sum_FG.PeakArea - coefficients)
# 
# ## ---- correcting for batch effects ----
# 
# grouped.data <- group_by(peptides.data, batch_date, batch.exp.n, R.Label,sample, replicate, EG.StrippedSequence)
# peptides.grouped_stats <- dplyr::summarise(grouped.data,
#                                   count = n(),
#                                   sum_FG.PeakArea = sum(F.PeakArea, na.rm=T))
# 
# peptides.grouped.batches = merge(peptides.grouped_stats, batch.fits, by=c("batch.exp.n", "EG.StrippedSequence"), all=T)
# peptides.grouped.batches$corrected = with(peptides.grouped.batches, sum_FG.PeakArea - coefficients)
# 
# p = ggplot(peptides.grouped.batches, aes(y=log(corrected,base=10), colour=batch_date, x=batch.exp.n)) + 
#           geom_boxplot() +
#           theme(aspect.ratio = 1) +
#           ggtitle(paste("Grouped by", paste(attr(grouped.data, which="vars"), collapse=".")))
# 
# plots.list = lappend(plots.list, p)
# 
# peptides.grouped.batches = tbl_df(peptides.grouped.batches)
# peptides.grouped.batches = peptides.grouped.batches[which(!is.na(peptides.grouped.batches$corrected)),]
# peptides.grouped.batches[peptides.grouped.batches$batch.exp.n == 1, "corrected"] = with(peptides.grouped.batches[peptides.grouped.batches$batch.exp.n == 1,], coefficients)
# 
# peptides.df.corrected = dcast(data=peptides.grouped.batches, formula=EG.StrippedSequence~sample+replicate+batch.exp.n, value.var="corrected")
# 
# peptides.matrix.corrected = as.matrix(peptides.df.corrected[,-1])
# rownames(peptides.matrix.corrected) = peptides.df.corrected$EG.StrippedSequence
# 
# batches.after = sub(x=colnames(peptides.matrix.corrected), pattern=".*?_([A-Za-z0-9]+)$", perl=T, replacement="\\1")
# 
# 
# 
# 
# # ---- dendrogram with linear batch correction----
# par(pty="s")
# h = hclust(dist(scale(t(peptides.matrix.corrected))))
# plot(h, labels=batches.after, cex=0.66 )
# p = recordPlot()
# plots.list = lappend(plots.list, p)
# dev.off()
# 
# 
# # ---- pca with bg correction----
# 
# par(pty="s", mfrow=c(1,2))
# pca = prcomp(t(peptides.matrix), scale.=T)
# x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
# y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
# 
# plot(pca$x[,1], pca$x[,2], col=batches, pch=16, main="Before adjustments for batch effects",
#      xlab=paste("PC1,", x_var), 
#      ylab=paste("PC2,", y_var))
# text(pca$x[,1], pca$x[,2], labels=batches, cex=0.66)
# 
# pca = prcomp(t(peptides.matrix.corrected), scale.=T)
# x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
# y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
# 
# plot(pca$x[,1], pca$x[,2], col=batches.after, pch=16, main="After adjustments for batch effects",
#      xlab=paste("PC1,", x_var), 
#      ylab=paste("PC2,", y_var),)
# text(pca$x[,1], pca$x[,2], labels=batches.after, cex=0.66)
# p = recordPlot()
# plots.list = lappend(plots.list, p)
# dev.off()
# 
# ## ---- batch adjustments using sva ----
# # sub(x=colnames(peptides.matrix), pattern="(.*?)_[A-Za-z0-9]+$", replacement="\\1")
# # mutants = as.character(sample_map[match(colnames(protein.matrix), sample_map$SampleName), "ORF" ])
# # mutants[mutants == "none"] = "Mix"
# # 
# # 
# # pheno = data.frame(row.names=colnames(protein.matrix), 
# #                    batch  = factor(batches), 
# #                    mutant = factor(mutants))
# # 
# # 
# # mod = model.matrix(~as.factor(batches), data=pheno)
# # mod0 = model.matrix(~1,data=pheno)
# # 
# # n.sv = num.sv(peptides.matrix, mod, method="leek")
# # 
# # View(peptides.matrix)
# # svobj = sva(peptides.matrix, mod, mod0, n.sv=n.sv)
# # sum(is.na(peptides.matrix)
# 
# # ---- Spectronaut batch effects ----
# 
# load(file="./R/objects/sample_exp.map.RData")
# sample_exp.map[with(sample_exp.map, order(SampleName)),]
# 
# peptides.data.s = droplevels(unique(select(peptides.data, R.Label, batch)))
# peptides.data.s = peptides.data.s[with(peptides.data.s, order(R.Label)),]
# 
# set.seed(123)
# to_select = sample(peptides.data.s[duplicated(peptides.data.s$R.Label),c("R.Label")], 20)
# 
# peptides.data.f = droplevels(peptides.data[peptides.data$R.Label %in% to_select,])
# 
# #peptides.data.f = droplevels(filter(peptides.data, R.Label == "KL_St_Mix_19" | R.Label == "KL_St_Mix_17"))
# 
# peptides.data.f$R.Label = factor(peptides.data.f$R.Label)
# p = ggplot(peptides.data.f, aes(x=log(FG.TotalPeakArea), colour=batch)) + 
#            geom_density() +
#            facet_wrap(~R.Label, scales="free") + 
#            theme(aspect.ratio = 1)
# 
# file_name = "spectronaut.effects.pdf"
# file_path = paste(figures_dir, file_name, sep="/")
# ggsave(filename=file_path, width=11.69+0.1*11.69, height=8.27+0.1*8.27)
# plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 
