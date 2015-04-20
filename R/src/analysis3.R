#!/usr/bin/env Rscript
# getting peaks ratios based on all vs all graphs

rm(list=ls())
source("./R/functions.R")

plots.list = list()
fun_name = "analysis3"

## ---- data_load ----
load("./R/objects/peptides.peak_sums.trimmed.RData")
load("./R/objects/protein_annotations.RData")
load("./R/objects/sample_exp.map.RData")
load("./R/objects/experiment.map.RData")
load("./R/objects/peak_sums.ratios.RData")

peak_sums.ratios$ORF.pair =  paste(peak_sums.ratios$ORF, peak_sums.ratios$pair, sep=":")
peak_sums.ratios.wide = dcast(peak_sums.ratios, formula=ORF.pair~R.Label, value.var="ratio")
ratios.matrix = as.matrix(peak_sums.ratios.wide[,-1])
rownames(ratios.matrix) = peak_sums.ratios.wide$ORF.pair
colnames(ratios.matrix)
pheno = data.frame(R.Label = colnames(ratios.matrix))
pheno$KO = droplevels(sample_exp.map$ORF[match(pheno$R.Label, sample_exp.map$SampleName)])
pheno$KO[pheno$R.Label == "KL_Try_027_c"] = "WT"
pheno$KO[is.na(pheno$KO)] = "none"

# pheno$batch.exp = factor(experiment_map$batch[match(pheno$R.Label, experiment_map$SampleName)])
# pheno$batch_date = factor(experiment_map$date[match(pheno$R.Label, experiment_map$SampleName)])
# pheno$batch.exp.n = factor(as.numeric(pheno$batch.exp))



X = model.matrix(~pheno$KO + 0)
colnames(X) = levels(pheno$KO)

reference = "WT"
matrix = log(ratios.matrix)

lm.fit_model = lmFit(matrix, X)
ph = unique(as.character(pheno$KO))

contrasts = paste0( ph[ph !=reference] ,"-", reference)  

mc = makeContrasts(contrasts=contrasts, levels=X)    
c.fit = contrasts.fit(lm.fit_model, mc)
eb = eBayes(c.fit)

folds = rowFolds(data=exp(matrix), groups=pheno$KO, reference=reference)
folds = log(folds, 2)


#folds_tmp = melt(as.matrix(folds), id.vars="row.names")

#merging results
folds_tmp = melt(eb$coefficients, id.vars="row.names")
#folds_tmp$contrasts = factor(paste(folds_tmp$contrasts, "-", reference, sep=""))
pvals_tmp = melt(eb$p.value, id.vars="row.names")


names(folds_tmp) = c("ORF.pair", "contrasts", "logFC")
names(pvals_tmp) = c("ORF.pair", "contrasts", "p.value")

folds_tmp$contrasts = factor(folds_tmp$contrasts)
pvals_tmp$contrasts = factor(pvals_tmp$contrasts)

ratios.FC = merge(folds_tmp, pvals_tmp, all=T,
                  by=c("ORF.pair", "contrasts"))

##multiple testing correction
ratios.FC$p.value_BH = p.adjust(ratios.FC$p.value, method="BH")
ratios.FC$p.value_bonferroni = p.adjust(ratios.FC$p.value, method="bonferroni")

ratios.FC$KO = sub(x = ratios.FC$contrasts, pattern=paste("(.*?)-", reference, sep=""), replacement="\\1")
ratios.FC$reference = reference


ratios.FC.f = ratios.FC[ratios.FC$KO %in% droplevels(unique(sample_exp.map$ORF[grep("Kinase", sample_exp.map$Type)])),]

p_thr = 0.05
fc_thr = 1
top = 1000000
tmp.f = table(droplevels(filter(ratios.FC.f, p.value_BH < p_thr, (logFC < -fc_thr | logFC > fc_thr) )$ORF.pair))
selected = c()
if (top >= length(tmp.f)) {
  selected = names(sort(-tmp.f))
} else {
  selected = names(sort(-tmp.f)[1:top])
}

ratios.FC.f = tbl_df(ratios.FC.f) %>% arrange(contrasts, ORF.pair)

ratios.FC.f.selected = droplevels(ratios.FC.f[ratios.FC.f$ORF.pair %in% selected,])
ratios.FC.wide = droplevels(dcast(ratios.FC.f.selected, formula=ORF.pair~KO, value.var="logFC"))

annotation = data.frame(batch_date = sample_exp.map$date[match(colnames(ratios.FC.wide)[-1], sample_exp.map$ORF)])
rownames(annotation) = colnames(ratios.FC.wide)[-1]

col_breaks = c(-5,-4,-3,-2,-1,1,2,3,4,5)

aheatmap(ratios.FC.wide[,-1],
         annCol=annotation,
         #annRow=factor(cl_res2[[cl_H]]$consensusClass),
         #Colv=cl_res[[cl_V]]$consensusTree,
         #Rowv=cl_res2[[cl_H]]$consensusTree,
         color = paste("-RdBu", length(col_breaks)-1, sep=":"),
         breaks = col_breaks)


grep(x=filter(ratios.FC.f, contrasts=="YAL017W-WT")$ORF.pair, pattern="YAL003W:")


library("igraph")

pattern.p = "([A-Za-z0-9-]+):(\\w+).(\\w+)"

tmp = strsplit2(ratios.FC.f$ORF.pair, split="[:]")
tmp2 = strsplit2(tmp[,2], split="[.]")
tmp.matrix = data.frame(cbind(tmp[,1], tmp2))
names(tmp.matrix) = c("ORF", "peptide1", "peptide2")

ratios.data.forw = cbind(tmp.matrix, ratios.FC.f)
ratios.data.rev = ratios.data.forw
ratios.data.forw$type = factor("direct")


ratios.data.rev[,"peptide1"] = ratios.data.forw[,"peptide2"]
ratios.data.rev[,"peptide2"] = ratios.data.forw[,"peptide1"]
ratios.data.rev$logFC = log(1/2^(ratios.data.rev$logFC),2)
ratios.data.rev$type = factor("reverse")
ratios.data = rbind(ratios.data.forw, ratios.data.rev)
ratios.data = ratios.data %>% arrange(ORF, contrasts) %>% 
                              group_by(ORF, contrasts) %>% 
                              distinct(ORF, contrasts, peptide1) %>% 
                              mutate(counts = n())


ratios.stats = ratios.data %>% group_by(ORF,contrasts, peptide1 ) %>% summarize(change = mean(logFC),
                                                                 meta_p.value = Stouffer.test(p.value))

