#!/usr/bin/env Rscript
# peptide peak changes based on minimally chagned peptide 

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


distinct(peptides.dataset, R.Label, EG.StrippedSequence, batch_date)
peptide2orfs = droplevels(dplyr::select(protein_annot.selected, strippedSequence, SystName))
names(peptide2orfs) = c("EG.StrippedSequence", "ORF")

peptides.dataset.f = droplevels(peptides.dataset[peptides.dataset$EG.StrippedSequence %in% peptide2orfs$EG.StrippedSequence,])

peak_sums.f = peptides.dataset.f
peak_sums.f$ORF = peptide2orfs$ORF[match(peptides.dataset.f$EG.StrippedSequence, peptide2orfs$EG.StrippedSequence)]
peak_sums.f = droplevels(peak_sums.f)
peak_sums.f = peak_sums.f %>% distinct(R.Label, EG.StrippedSequence) %>% arrange(R.Label, ORF, EG.StrippedSequence)

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
names(folds_tmp) = c("ORF.peptide", "contrasts", "logFC")

folds_tmp$contrasts = factor(paste(folds_tmp$contrasts, "-", reference, sep=""))

pvals_tmp = melt(eb$p.value, id.vars="row.names")

names(pvals_tmp) = c("ORF.peptide", "contrasts", "p.value")

folds_tmp$contrasts = factor(folds_tmp$contrasts)
pvals_tmp$contrasts = factor(pvals_tmp$contrasts)

#merging results

peptides.FC = merge(folds_tmp, pvals_tmp, ,
                    by=c("ORF.peptide", "contrasts"))

peptides.FC = peptides.FC %>% extract(ORF.peptide, 
                                     into=c("ORF", "peptide"), 
                                     regex="([A-Za-z0-9-]+).([A-Z]+)$")

peptides.FC$KO = sub(x = peptides.FC$contrasts, pattern=paste("(.*?)-", reference, sep=""), replacement="\\1")
peptides.FC$reference = reference
peptides.FC = tbl_df(peptides.FC)



### making relative to reference
peptides.long = tbl_df(melt(log(exp(matrix),2), id.vars="row.names"))
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

tmp = dcast(peptides.merged, formula=ORF+peptide~R.Label, value.var="ratio")
stopifnot(sum(is.na(tmp)) == 0)

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
        #z = filter(peptides.merged, KO == "YAR018C", ORF=="YAL003W" )
        
        peptides = unique(z$peptide)
        if (length(peptides) == 1) {
          return(data.frame(peptide=peptides, ratio.FC = 1, p.value=NA, reference=peptides))
        }
        
        #peptides = peptides[peptides != z$min.peptide[1]]
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

stopifnot(sum(is.na(dcast(peptides.ratios, formula=peptide~KO, value.var="ratio.FC"))) == 0)
peptides.ratios$padj = p.adjust(peptides.ratios$p.value, method="BH")

file_name = "peptides.ratios.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.ratios, file=file_path)


## -- robust peptides ----
padj_thr = 0.1
par(pty="s")
peptides.ratios.stats = peptides.ratios %>% group_by(ORF, peptide) %>% summarise(count = sum(padj > padj_thr | is.na(p.value)))

selected.peptides = peptides.ratios.stats$peptide[peptides.ratios.stats$count >= median(peptides.ratios.stats$count)]

protein.stats.before = peak_sums.f %>% group_by(ORF) %>% summarise(pep.count = length(unique(EG.StrippedSequence)))
protein.stats.before$category = factor(ifelse(protein.stats.before$pep.count >=5,">5", protein.stats.before$pep.count))

protein.stats.after = peptides.ratios.stats %>% filter(count >= median(peptides.ratios.stats$count)) %>% group_by(ORF) %>% summarize(pep.count = n())                                                                                             
protein.stats.after$category = factor(with(protein.stats.after, ifelse(pep.count >=5 ,">5", pep.count)))

protein.stats.before$type = "before"
protein.stats.after$type = "after"
protein.stats = rbind(protein.stats.before, protein.stats.after)
protein.stats$type = factor(protein.stats$type, levels =c("before", "after"))

      
p1 = ggplot(data=protein.stats, aes(x=type, fill=factor(category))) +
    geom_bar() +
    ggtitle("Robust peptides")

p2 = ggplot(peptides.ratios.stats, aes(x=count)) +
       geom_histogram() +
       geom_vline(xintercept=median(peptides.ratios.stats$count), linetype=2, colour="red")

g = arrangeGrob(p2, p1)
plots.list = lappend(plots.list, g)

file_name = "peptides.ratios.stats.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.ratios.stats, file=file_path)



## -- PTM enrichments ---- 
Ph = c("S","T","Y")
Ac_Ub_SM = c("K")
Me = c("K","R","C")
Ca = c("E")
NG = c("N")
any.ptm = unique(c(Ph, Ac_Ub_SM, Me, Ca, NG))

peptides.ratios$peptide = as.character(peptides.ratios$peptide)
peptides.ratios.PTM = peptides.ratios %>% group_by(KO, ORF, peptide) %>% mutate(Ph = sum(unlist(strsplit(peptide, "")) %in% Ph),
                                                                               Ac_Ub_SM = sum(unlist(strsplit(peptide, "")) %in% Ac_Ub_SM),
                                                                               Me = sum(unlist(strsplit(peptide, "")) %in% Me),
                                                                               Ca = sum(unlist(strsplit(peptide, "")) %in% Ca),
                                                                               NG = sum(unlist(strsplit(peptide, "")) %in% NG),
                                                                               any.ptm = sum(unlist(strsplit(peptide, "")) %in% any.ptm))



peptides.ratios.PTM.f = filter(peptides.ratios.PTM, ratio.FC != 1 , p.value != 1)
peptides.ratios.PTM.f$padj = p.adjust(peptides.ratios.PTM.f$p.value, method="BH")


## -- PTM enrichment ----
padj_thr = 0.05
toPlot = filter(peptides.ratios.PTM.f, is.na(p.value) == F)

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

col_breaks = c(-4,-3,-2,-1,-0.5,0.5,1,2,3,4)
pheatmap(log(tmp[,-c(1,2)]), show_rownames=F, 
         breaks=col_breaks, 
         color=brewer.pal(n=length(col_breaks)-1, "PuOr"),
         annotation=annotation) 

p = recordPlot()
plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l")
