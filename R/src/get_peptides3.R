#!/usr/bin/env Rscript
# batch correction based on fragments
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "get_peptides3"

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

p = ggplot(fragments.data, aes(x = log(qvalue.median))) + 
  geom_density() +
  geom_vline(xintercept = log(0.01))
plots.list = lappend(plots.list, p)


fragments.data.f = filter(fragments.data, qvalue.median <= 0.01, F.IsQuantified == T)

exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
exp_metadata$batch_kmeans = kmeans(exp_metadata$aquisition_date.str, 7)$cl


fragments.df = dcast(fragments.data.f, formula=FG.Id+fragment~R.Label, value.var="F.PeakArea")

fragments.matrix = as.matrix(fragments.df[,-c(1,2)])
rownames(fragments.matrix) = paste(fragments.df$FG.Id, fragments.df$fragment, sep=":")

pheno = exp_metadata[match(colnames(fragments.matrix), exp_metadata$sample_name),]
mod = model.matrix(~as.factor(ORF), data=pheno)

fragments.matrix.combat = ComBat(log(fragments.matrix), batch=pheno$batch_kmeans, mod=mod, par.prior=T)



before = log(fragments.matrix)
after = fragments.matrix.combat

# pca
message("plotting PCA results")

file_name = paste("PCA_batch_effects", fun_name, "pdf", sep=".")
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


plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 




