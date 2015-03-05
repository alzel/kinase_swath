rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "batch_effects"

load("./R/objects/proteins.deseq.combat.long.RData")
load("./R/objects/proteins.matrix.f.deseq.combat.RData")
load("./R/objects/sample_exp.map.RData")
load("./R/objects/peptides.matrix.f.RData")


proteins.matrix = proteins.matrix.f.deseq.combat
proteins.long   = proteins.deseq.combat.long

pattern.p = "(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$"
matches = stringr::str_match_all(pattern=pattern.p, colnames(proteins.matrix))

stopifnot(sum(lapply(matches,length)!=0) == ncol(proteins.matrix))
pheno = data.frame(matrix(unlist(matches), ncol=length(matches[[1]]), byrow=T))
colnames(pheno) = c("name", "R.Label", "batch_date", "batch.exp.n", "batch" )
rownames(pheno) = colnames(proteins.matrix)
pheno$ORF = droplevels(sample_exp.map$ORF[match(pheno$R.Label, sample_exp.map$SampleName)])
pheno$ORF[pheno$R.Label == "KL_Try_027_c"] = "WT"
pheno$batch.exp.n[pheno$R.Label == "KL_Try_027_c"] = 5


pheno = droplevels(pheno[pheno$ORF != "none",]) #removing none-phenotype samples
proteins.matrix.f = proteins.matrix[,match(pheno$name, colnames(proteins.matrix))]



X = model.matrix(~pheno$ORF + 0)
colnames(X) = levels(pheno$ORF)

reference = "WT"
matrix = proteins.matrix.f

lm.fit_model = lmFit(matrix, X)
ph = unique(as.character(pheno$ORF))
contrasts = paste0( ph[ph !=reference] ,"-", reference)  

mc = makeContrasts(contrasts=contrasts, levels=X)    
c.fit = contrasts.fit(lm.fit_model, mc)
eb = eBayes(c.fit)
  
eb$coefficients
folds = rowFolds(data=2^matrix, groups=pheno$ORF, reference=reference)
folds = log(folds, 2)

#folds_tmp = melt(as.matrix(folds), id.vars="row.names")

#merging results
folds_tmp = melt(eb$coefficients, id.vars="row.names")
#folds_tmp$contrasts = factor(paste(folds_tmp$contrasts, "-", reference, sep=""))
pvals_tmp = melt(eb$p.value, id.vars="row.names")


names(folds_tmp) = c("ORF", "contrasts", "logFC")
names(pvals_tmp) = c("ORF", "contrasts", "p.value")

folds_tmp$contrasts = factor(folds_tmp$contrasts)
pvals_tmp$contrasts = factor(pvals_tmp$contrasts)


proteins.FC = merge(folds_tmp, pvals_tmp, all=T,
                    by=c("ORF", "contrasts"))

##multiple testing correction
proteins.FC$p.value_BH = p.adjust(proteins.FC$p.value, method="BH")
proteins.FC$p.value_bonferroni = p.adjust(proteins.FC$p.value, method="bonferroni")




lb = -4
ub = 4
toPlot = proteins.FC
toPlot[toPlot$logFC < 0,]$logFC = ifelse(toPlot[toPlot$logFC < 0,]$logFC < lb, lb, toPlot[toPlot$logFC < 0,]$logFC)
toPlot[toPlot$logFC > 0,]$logFC = ifelse(toPlot[toPlot$logFC > 0,]$logFC > ub, ub, toPlot[toPlot$logFC > 0,]$logFC)

FC_thr = 1
pval_thr = 0.05
proteins.FC.stats = data.frame( ratio_sign = round(sum(proteins.FC$logFC < -FC_thr & proteins.FC$p.value_BH < pval_thr)/sum(proteins.FC$logFC > FC_thr & proteins.FC$p.value_BH < pval_thr),2),
                                ratio =  round(sum(proteins.FC$logFC < 0)/sum(proteins.FC$logFC > 0),2),
                                n_prot = length(unique(proteins.FC$ORF)),
                                n_sign = round(sum(proteins.FC$p.value_BH<0.05)/length(unique(proteins.FC$contrasts)),2),
                                x_min = -4,
                                y_max =max(-log10(proteins.FC$p.value_BH)))

toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)


p1 = ggplot(toPlot, aes(y=-log10(p.value_BH), x=logFC)) +
         geom_point(aes(color = sign), alpha=0.5) + 
         xlim(c(lb,ub)) +
         geom_hline(y=-log(pval_thr,10),linetype=3) +
         geom_vline(x=c(FC_thr,-FC_thr),linetype=3) +
         geom_text(data=proteins.FC.stats, aes(x=x_min, y=y_max, hjust=0, vjust=1,
                                             label=paste(paste0("#_mean_sign = ", n_sign),
                                                         paste0("up/down_ratio = ",ratio),
                                                         paste0("ratio_sign = ",   ratio_sign),
                                                         paste0("#_prot = ", n_prot), sep="\n"))) +
        xlab(paste("Log2(fold-change)")) +
        theme(aspect.ratio = 1, legend.position = "none", text=element_text(size=18))

p2 = ggplot(toPlot, aes(x=logFC)) +
            geom_histogram(colour = "white", fill = "black", binwidth = 0.15) +
            xlab(paste("Log2(fold-change)")) +
            theme(aspect.ratio = 1, text=element_text(size=18))

g = arrangeGrob(p1,p2 , ncol=2, main=paste("Reference:", reference))
plots.list = lappend(plots.list, g)

library("RColorBrewer")

proteins.FC.wide = dcast(proteins.FC, formula=ORF~contrasts, value.var="logFC")

col_breaks = c(-4,-2,-1,1,2,4)

pheatmap(proteins.FC.wide[,-1],
         color = rev(brewer.pal(name="RdBu", n=(length(col_breaks))-1)),
         breaks = col_breaks)



s = princomp(proteins.FC.wide[,-1])

plot(s$scores[,])
pam()



display.brewer.all()
display.brewer.all()
dev.off()
#colorRampPalette(c("navy", "white", "firebrick3"))(length(seq(-4,4,0.5)))


?pheatmap
