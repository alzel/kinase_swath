rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "batch_effects"

load("./R/objects/peptides.peak_sums.trimmed.RData")
load("./R/objects/protein_annotations._load_.RData")
load("./R/objects/peptides.cor.stats.top.RData")
load("./R/objects/sample_map._load_.RData")
#load("./R/objects/sample_exp.map.RData")

peptides.long = peptides.peak_sums.trimmed
peptides.selected = tbl_df(droplevels(filter(peptides.cor.stats.top, top == "3")))

peptides.long.selected = tbl_df(peptides.long[peptides.long$EG.StrippedSequence %in% peptides.selected$EG.StrippedSequence,])
peptides.long.selected$ORF = peptides.selected$ORF[match(peptides.long.selected$EG.StrippedSequence, peptides.selected$EG.StrippedSequence)]

proteins.long = peptides.long.selected  %>% group_by(ORF, R.Label, batch_date, batch.exp.n, batch) %>% summarise(mean_signal = mean(T_signal, na.rm=T),
                                                                                                                 sum_signal = log(sum(exp(T_signal), na.rm=T)))

peptides.df = dcast(data=peptides.peak_sums.trimmed, formula=EG.StrippedSequence~R.Label+batch_date+batch.exp.n+batch, value.var="T_signal")
proteins.df = dcast(data=proteins.long, formula=ORF~R.Label+batch_date+batch.exp.n+batch, value.var="sum_signal")

proteins.matrix = as.matrix(proteins.df[,-1])
rownames(proteins.matrix) = proteins.df$ORF


pattern.p = "(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$"
matches = stringr::str_match_all(pattern=pattern.p, colnames(proteins.matrix))

stopifnot(sum(lapply(matches,length)!=0) == ncol(proteins.matrix))

pheno = data.frame(matrix(unlist(matches), ncol=length(matches[[1]]), byrow=T))
colnames(pheno) = c("name", "R.Label", "batch_date", "batch.exp.n", "batch" )
rownames(pheno) = colnames(proteins.matrix)
pheno$ORF = droplevels(sample_map$ORF[match(pheno$R.Label, sample_map$SampleName)])
pheno$ORF[pheno$R.Label == "KL_Try_027_c"] = "WT"
pheno$batch.exp.n[pheno$R.Label == "KL_Try_027_c"] = 5
pheno = droplevels(pheno[which(!is.na(pheno$ORF)),])



pheno$group = pheno$batch_date #grouping variable to estimate batch effects


tmp.factor_size = ddply(pheno, .(group, ORF), summarise, factor_size = length(ORF))
tmp.factor_size$group.ORF = paste(tmp.factor_size$group, tmp.factor_size$ORF, sep=".")

pheno$group.ORF = paste(pheno$group, pheno$ORF, sep=".")
pheno = droplevels(pheno[pheno$group.ORF %in% tmp.factor_size$group.ORF[tmp.factor_size$factor_size >=2],])
pheno = droplevels(pheno[pheno$group %in% levels(droplevels(pheno$group[grep(pattern="mix", pheno$R.Label, ignore.case=T)])),])


proteins.matrix.f = proteins.matrix[,match(rownames(pheno), colnames(proteins.matrix))]

mod = model.matrix(~as.factor(ORF), data=pheno)

tmp.size_factors = DESeq::estimateSizeFactorsForMatrix(exp(proteins.matrix.f))
proteins.matrix.f.deseq = log(exp(proteins.matrix.f)/tmp.size_factors)


proteins.matrix.f.combat = ComBat(na.omit(proteins.matrix.f), batch=pheno$group, mod=mod, par.prior=T)
proteins.matrix.f.deseq.combat = ComBat(na.omit(proteins.matrix.f.deseq), batch=pheno$group, mod=mod, par.prior=T)
proteins.matrix.f.quantiles.combat = ComBat(normalizeQuantiles(na.omit(proteins.matrix.f.deseq)), batch=pheno$group, mod=mod, par.prior=T)




file_name = "proteins.matrix.f.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.matrix.f.combat,file=file_path)  

file_name = "proteins.matrix.f.deseq.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.matrix.f.deseq.combat,file=file_path)  

file_name = "proteins.matrix.f.quantiles.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.matrix.f.quantiles.combat,file=file_path)  

## ---- using SVA to adjust for batch effects----
if (FALSE) {
  
  
  proteins.matrix = as.matrix(proteins.df[,-1])
  
  peptides.matrix = as.matrix(peptides.df[,-1])
  rownames(proteins.matrix) = proteins.df$ORF
  rownames(peptides.matrix) = peptides.df$EG.StrippedSequence
  
  pattern.p = "(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$"
  matches = stringr::str_match_all(pattern=pattern.p, colnames(proteins.matrix))
  matches = stringr::str_match_all(pattern=pattern.p, colnames(peptides.matrix))
  
  
  stopifnot(sum(lapply(matches,length)!=0) == ncol(proteins.matrix))
  
  pheno.sva = data.frame(matrix(unlist(matches), ncol=length(matches[[1]]), byrow=T))
  colnames(pheno.sva) = c("name", "R.Label", "batch_date", "batch.exp.n", "batch" )
  rownames(pheno.sva) = colnames(proteins.matrix)
  pheno.sva$ORF = droplevels(sample_map$ORF[match(pheno.sva$R.Label, sample_map$SampleName)])
  pheno.sva$ORF[pheno$R.Label == "KL_Try_027_c"] = "WT"
  pheno.sva$batch.exp.n[pheno.sva$R.Label == "KL_Try_027_c"] = 5
  pheno.sva = droplevels(pheno.sva[which(!is.na(pheno.sva$ORF)),])
    
  #pheno.sva = droplevels(pheno.sva[pheno.sva$ORF != "none",])
  proteins.matrix.sva = proteins.matrix[,match(rownames(pheno.sva), colnames(proteins.matrix))]
  
  peptides.matrix.sva = peptides.matrix[,match(rownames(pheno.sva), colnames(peptides.matrix))]
  
#   tmp.size_factors = DESeq::estimateSizeFactorsForMatrix(exp(proteins.matrix.sva))
#   proteins.matrix.sva.deseq = log(exp(proteins.matrix.sva)/tmp.size_factors)
#   
  
  mod.sva = model.matrix(~ORF, data=pheno.sva)
  mod0.sva = model.matrix(~1, data=pheno.sva)
  
  n.sv = num.sv(proteins.matrix.sva, mod=mod.sva)
  
  n.sv = num.sv(peptides.matrix.sva, mod=mod.sva)
  
  svaobj <- sva::sva(dat=proteins.matrix.sva, mod=mod.sva, n.sv=n.sv)
  svaobj <- sva::sva(dat=peptides.matrix.sva, mod=mod.sva, n.sv=n.sv)
  
  svaobj <- sva.patched(dat=proteins.matrix.sva, mod=mod.sva, mod0=mod0.sva, n.sv=n.sv)
  
  svaobj <- sva.patched(dat=peptides.matrix.sva, mod=mod.sva, mod0=mod0.sva, n.sv=n.sv, method="two-step")
  
  
  svaX = model.matrix(~ORF + svaobj$sv, data=pheno.sva)
    
  #fit = lmFit(proteins.matrix.sva, svaX)
  fit = lmFit(peptides.matrix.sva, svaX)
  
  Batch <- fit$coef[,(nlevels(pheno.sva$ORF)+1):(nlevels(pheno.sva$ORF)+svaobj$n.sv)] %*% t(svaX[,(nlevels(pheno.sva$ORF)+1):((nlevels(pheno.sva$ORF))+svaobj$n.sv)])
  
  
  # and implement it like this:
    
  proteins.matrix.sva.adjusted = proteins.matrix.sva-Batch
  peptides.matrix.sva.adjusted = peptides.matrix.sva-Batch
  sum(proteins.matrix.sva.adjusted <0)
  
  file_name = "proteins.matrix.sva.adjusted.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(proteins.matrix.sva.adjusted,file=file_path) 
  
}



before = proteins.matrix.f
after  = proteins.matrix.f.deseq.combat

# pca
message("plotting PCA results")

file_name = "PCA_batch_effects.pdf"
file_path = paste(figures_dir, file_name, sep="/")
pdf(file_path, width=11.7+0.1*11.7, height=8.27+0.1*8.27)
par(pty="s", mfrow=c(1,2))

pca = prcomp(t(before), scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], cex=1.5, cex.lab=1.5, col=pheno$batch_date, pch=16, main="Before adjustments for batch effects", 
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=pheno$batch.exp.n, cex=0.5)

pca = prcomp(t(after), scale.=T)
#pca = prcomp(t(proteins.matrix.f., scale.=T)
x_var = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca$x[,1], pca$x[,2], cex=1.5, cex.lab=1.5, col=pheno$batch_date, pch=16, main="After adjustments for batch effects",
     xlab=paste("PC1,", x_var), 
     ylab=paste("PC2,", y_var))
text(pca$x[,1], pca$x[,2], labels=pheno$batch.exp.n, cex=0.5)

p = recordPlot()
plots.list = lappend(plots.list, p)

dev.off()

# dendrograms 
# file_name = "Clustering_batch_effects.png"
# file_path = paste(figures_dir, file_name, sep="/")
# png(file_path, width=297, height=210, units="mm", res=150)
par(mfrow=c(2,1))
h = hclust(dist(scale(t(proteins.matrix.f))))
plot(h, labels=pheno$group, cex=0.5, main="Before batch correction")

h = hclust(dist(scale(t(proteins.matrix.sva.adjusted))))
plot(h, labels=pheno$group, cex=0.5, main="After batch correction")

p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()


peptides.matrix.sva.adjusted
#tidying batch corrected data
message("Tidying up batch corrected data", appendLF=F)

tmp.wide = data.frame(proteins.matrix.f.deseq.combat)
tmp.wide$ORF = rownames(proteins.matrix.f.deseq.combat)
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



tmp.wide = data.frame(proteins.matrix.f.combat)
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


tmp.wide = data.frame(proteins.matrix.f)
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
  stat_function(fun=dnorm, args=list(mean=mean(toPlot$value), sd=sd(toPlot$value)))
plots.list = lappend(plots.list, p)

p = p + facet_wrap(~R.Label, scales="free")
plots.list = lappend(plots.list, p)

p1 = ggplot(proteins.long, aes(x=batch.exp.n, y=value, fill=batch.exp.n)) +
  geom_boxplot() +
  theme(legend.position="none")

p3 = ggplot(proteins.deseq.combat.long, aes(x=batch.exp.n, y=value, fill=batch.exp.n)) +
  geom_boxplot() + 
  theme(legend.position="top")

g = arrangeGrob(p1, p3)
plots.list = lappend(plots.list, g)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 
