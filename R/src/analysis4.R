#!/usr/bin/env Rscript
# getting peaks ratios based on all vs all graphs

rm(list=ls())
source("./R/functions.R")

plots.list = list()
fun_name = "analysis4"

## ---- data_load ----
load("./R/objects/peptides.peak_sums.trimmed.RData")
load("./R/objects/protein_annotations.RData")
load("./R/objects/sample_exp.map.RData")
load("./R/objects/experiment.map.RData")
load("./R/objects/peak_sums.ratios.RData")



peptides.dataset = tbl_df(peptides.peak_sums.trimmed)
protein_annot.selected = droplevels(filter(protein_annotations, specificity == "unique"))

peptide2orfs = dplyr::select(protein_annot.selected, strippedSequence, SystName)
names(peptide2orfs) = c("EG.StrippedSequence", "ORF")

peptides.dataset.f = droplevels(peptides.dataset[peptides.dataset$EG.StrippedSequence %in% protein_annot.selected$strippedSequence,])

peak_sums = peptides.dataset.f
peak_sums.f = tbl_df(droplevels(merge(peak_sums, peptide2orfs, by="EG.StrippedSequence")))
peak_sums.f = peak_sums.f %>% distinct(R.Label, batch.exp.n, EG.StrippedSequence) %>% arrange(R.Label, ORF, EG.StrippedSequence)

peptides.wide = dcast(data=peak_sums.f, formula=ORF+EG.StrippedSequence~R.Label+batch.exp.n, value.var="T_signal")
peptides.matrix = as.matrix(peptides.wide[,-c(1,2)])

rownames(peptides.matrix) = paste(peptides.wide$ORF, peptides.wide$EG.StrippedSequence, sep=".")

pattern.p = "(.*?)_([A-Za-z0-9]+)$"
matches = stringr::str_match_all(pattern=pattern.p, colnames(peptides.matrix))

stopifnot(sum(lapply(matches,length)!=0) == ncol(peptides.matrix))
pheno = data.frame(matrix(unlist(matches), ncol=length(matches[[1]]), byrow=T))

colnames(pheno) = c("name", "R.Label", "batch.exp.n" )
rownames(pheno) = colnames(peptides.matrix)
pheno$KO = droplevels(sample_exp.map$ORF[match(pheno$R.Label, sample_exp.map$SampleName)])
# pheno$KO[pheno$R.Label == "KL_Try_027_c"] = "WT"
# pheno$batch.exp.n[pheno$R.Label == "KL_Try_027_c"] = 5
pheno[is.na(pheno[,"KO"]),"KO"] = "none"
pheno = droplevels(pheno[pheno$KO != "none",]) #removing none-phenotype samples
peptides.matrix.f = peptides.matrix[,match(pheno$name, colnames(peptides.matrix))]

X = model.matrix(~pheno$KO + 0)
colnames(X) = levels(pheno$KO)

reference = "WT"
matrix = peptides.matrix.f

lm.fit_model = lmFit(matrix, X)
ph = unique(as.character(pheno$KO))
contrasts = paste0( ph[ph !=reference] ,"-", reference)  

mc = makeContrasts(contrasts=contrasts, levels=X)    
c.fit = contrasts.fit(lm.fit_model, mc)
eb = eBayes(c.fit)


folds = rowFolds(data=exp(matrix), groups=pheno$KO, reference=reference)
folds = log(folds, 2)

folds_tmp = melt(as.matrix(folds), id.vars="row.names")

#merging results
#folds_tmp = melt(eb$coefficients, id.vars="row.names")
#folds_tmp$contrasts = factor(paste(folds_tmp$contrasts, "-", reference, sep=""))
pvals_tmp = melt(eb$p.value, id.vars="row.names")


names(folds_tmp) = c("ORF.peptide", "contrasts", "logFC")
names(pvals_tmp) = c("ORF.peptide", "contrasts", "p.value")

folds_tmp$contrasts = factor(folds_tmp$contrasts)
pvals_tmp$contrasts = factor(pvals_tmp$contrasts)


peptides.FC = merge(folds_tmp, pvals_tmp, all=T,
                    by=c("ORF.peptide", "contrasts"))

peptides.FC = peptides.FC %>% extract(ORF.peptide, 
                                     into=c("ORF", "peptide"), 
                                     regex="([A-Za-z0-9-]+).([A-Z]+)$")

peptides.FC$KO = sub(x = peptides.FC$contrasts, pattern=paste("(.*?)-", reference, sep=""), replacement="\\1")
peptides.FC$reference = reference
peptides.FC = tbl_df(peptides.FC)



### making relative to reference

peptides.long = tbl_df(melt(log(exp(peptides.matrix.f),2), id.vars="row.names"))



names(peptides.long) = c("ORF.peptide", "variable", "value")

peptides.long = peptides.long %>% extract(variable, 
                                          into=c("R.Label", "batch.exp.n"), 
                                          regex="(.*?)_([A-Za-z0-9]+)$")


peptides.long = peptides.long %>% extract(ORF.peptide, 
                                          into=c("ORF", "peptide"), 
                                          regex="([A-Za-z0-9-]+).([A-Z]+)$")


peptides.long$KO = sample_exp.map$ORF[match(peptides.long$R.Label, sample_exp.map$SampleName)]
# peptides.long.tmp = peptides.long
# 
# peptides.long = peptides.long.tmp


ref_constant = filter(peptides.long, KO == "WT")
ref_constant.summary = ref_constant %>% group_by(ORF,peptide) %>% summarise(mean.value.ref = log(mean(2^(value), na.rm=T),2) )


peptides.long = tbl_df(merge(peptides.long, ref_constant.summary, by = c("ORF","peptide")))


peptides.FC.min = peptides.FC %>% group_by(KO, ORF) %>% summarise(min.peptide = peptide[which(abs(logFC) == min(abs(logFC), na.rm=T))],
                                                                  min.logFC   = logFC[which(abs(logFC) == min(abs(logFC), na.rm=T))])

peptides.long = peptides.long %>% mutate(peptide.change = value - mean.value.ref)

peptides.merged = tbl_df(merge(peptides.long, peptides.FC.min, by = c("KO", "ORF")))
peptides.merged = peptides.merged %>% group_by(KO, R.Label, ORF) %>% mutate(ratio = 2^(peptide.change - min.logFC[1]))


# u.test = filter(peptides.merged, (KO == "YAL017W" ) & ORF == "YAL003W" & peptide == "AFQSAYPEFSR")
# u.wt = filter(peptides.merged, (KO == "WT" ) & ORF == "YAL003W" & peptide == "AFQSAYPEFSR")

# log(mean(exp(matrix[rownames(matrix) == "YAL003W.AFQSAYPEFSR", colnames(matrix) %in% as.character(pheno[pheno$KO == "YAL017W", "name"])])),2) - 
# log(mean(exp(matrix[rownames(matrix) == "YAL003W.AFQSAYPEFSR", colnames(matrix) %in% as.character(pheno[pheno$KO == "WT", "name"])])),2)
# 
# log(mean(2^(u.test$value)),2) - log(mean(2^u.wt$value),2)


peptides.ratios = ddply(peptides.merged, .(KO, ORF), 
      .fun=function(z) {
        v <<- z
        z = v
        
        peptides = unique(z$peptide)
        if (length(peptides) == 1) {
          return(data.frame(peptide=peptides, ratio.FC = 1, p.value=NA, reference=peptides))
        }
        
        peptides = peptides[peptides != z$min.peptide[1]]
        res = matrix(ncol=2, nrow=length(peptides))
        rownames(res) = peptides
        
        reference.ratio = z$ratio[z$peptide == z$min.peptide[1]]
        
        for (i in 1:length(peptides) ) {
          compare2 = z$ratio[z$peptide == peptides[i]]
          test = t.test(reference.ratio, compare2)  
          
          res[i,1] = mean(compare2, na.rm=T)/mean(reference.ratio, na.rm=T)
          res[i,2] = test$p.value
        }
        
        ret = data.frame(res)
        ret$peptide = rownames(ret)
        ret$reference = z$min.peptide[1]
        names(ret) = c("ratio.FC","p.value", "peptide","reference")  
        return(ret[,c(3,1,2,4)])
      })


Ph = c("S","T","Y")
Ac_Ub_SM = c("K")
Me = c("K","R","C")
Ca = c("E")
NG = c("N")
any.ptm = unique(c(Ph, Ac_Ub_SM, Me, Ca, NG))

peptides.ratios$peptide = as.character(peptides.ratios$peptide)
peptide.ratios.PTM = peptides.ratios %>% group_by(KO, ORF, peptide) %>% mutate(Ph = sum(unlist(strsplit(peptide, "")) %in% Ph),
                                                                               Ac_Ub_SM = sum(unlist(strsplit(peptide, "")) %in% Ac_Ub_SM),
                                                                               Me = sum(unlist(strsplit(peptide, "")) %in% Me),
                                                                               Ca = sum(unlist(strsplit(peptide, "")) %in% Ca),
                                                                               NG = sum(unlist(strsplit(peptide, "")) %in% NG),
                                                                               any.ptm = sum(unlist(strsplit(peptide, "")) %in% any.ptm))


Ph = c("S","T","Y")
Ac_Ub_SM = c("K")
Me = c("K","R","C")
Ca = c("E")
NG = c("N")
all = unique(c(Ph, Ac_Ub_SM, Me, Ca, NG))

# peptides.ratios$peptide = factor(peptides.ratios$peptide)
# peptides.ratios$ORF = factor(peptides.ratios$ORF)
# peptides.ratios$KO = factor(peptides.ratios$KO)
# peptides.ratios.PTM = ddply(peptides.ratios, .(KO, ORF, peptide),
#                            .fun = function(z) {
#                              v <<- z
#                              z = v
#                              
#                              peptide = unlist(strsplit(as.character(z$peptide), ""))
#                              
#                              ret = data.frame(Ph = sum(peptide %in% Ph),
#                                         Ac_Ub_SM = sum(peptide %in% Ac_Ub_SM),
#                                         Me = sum(peptide %in% Me),
#                                         Ca = sum(peptide %in% Ca),
#                                         NG = sum(peptide %in% NG),
#                                         any.ptm = sum(peptide %in% all))
#                              return(ret)
#                            })



peptide.ratios.PTM$padj = p.adjust(peptide.ratios.PTM$p.value, method="BH")


## -- PTM enrichment ----
padj_thr = 0.05
toPlot = filter(peptide.ratios.PTM, is.na(p.value) == F)
# toPlot$changed = factor(ifelse(toPlot$padj <= padj_thr & abs(log(toPlot$ratio.FC,2)) >= 1 , 1, 0))
# toPlot$hasPTM = factor(ifelse(toPlot$Ph <=1, 1, 0))


a = sum(ifelse(toPlot$padj <= padj_thr & toPlot$Ph > 0, 1, 0)) #changed and has PTM
b = sum(ifelse(toPlot$padj <= padj_thr & toPlot$Ph == 0, 1, 0)) #changed and no PTM
c = sum(ifelse(toPlot$padj > padj_thr & toPlot$Ph > 0, 1, 0)) #not changed and has PTM
d = sum(ifelse(toPlot$padj > padj_thr & toPlot$Ph == 0, 1, 0)) #not changed and no PTM

PTM.check = matrix(c(a,b,c,d), ncol=2, byrow=T)
fisher.test(PTM.check, alternative="greater")


toPlot.long = melt(data.frame(toPlot), id.vars = c("KO","ORF", "peptide","ratio.FC","p.value", "reference", "padj") )
toPlot.long$value = factor(toPlot.long$value)

p = ggplot(filter(toPlot.long, padj <= padj_thr), aes(x=value, y=log(ratio.FC,2))) +
       geom_boxplot() +
       geom_hline(yintercept = 0, linetype = 3) +
       facet_wrap(~variable, scale="free") +
       theme(aspect.ratio = 1) +
       xlab("Number of potentially modified residues per peptide") +
       ylab("Ratio of peptide change (in comparison to WT) to minimally changed peptide per ORF, log2")    

plots.list = lappend(plots.list, p)



tmp = dcast(data=peptides.ratios, formula=ORF+peptide~KO, value.var="ratio.FC")
tmp[is.na(tmp)] = 1

annotation = data.frame(batch_date = sample_exp.map$date[match(colnames(tmp)[-c(1,2)], sample_exp.map$ORF)])
rownames(annotation) = colnames(tmp)[-c(1,2)]

col_breaks = c(-5,-4,-3,-2,-1,1,2,3,4,5)
aheatmap(log(tmp[,-c(1,2)]),
         annCol=annotation,
         color = paste("-RdBu", length(col_breaks)-1, sep=":"),
         breaks = col_breaks)

plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")



