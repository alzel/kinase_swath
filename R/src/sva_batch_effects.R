rm(list = ls())
set.seed(123)

source("./R/boot.R")
source("./R/functions.R")


fun_name = "sva_batch_effects"

load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/peptides.matrix.RData")

pheno = droplevels(exp_metadata[match(colnames(peptides.matrix), exp_metadata$sample_name),])
plots.list = list()


mod = model.matrix(~1, data=pheno)
mod0 = mod
fraction = 0.5
n_genes = floor(nrow(peptides.matrix)*fraction)
# choose ONE surrogate variable
n.sv  <- 1
# estimate surrogate variables

# Here, we use the lowly variable genes as control genes,
# i.e. we exclude the top 500 most variable genes from the
# estimation of the correction factors.
# This might prevent most of the biological signal from being regressed out.

edata = as.matrix(peptides.matrix)

varsEdata <- rowSds(edata)/rowMeans(edata)

selectEdata <- order(varsEdata, decreasing = TRUE)[seq(n_genes, length(varsEdata))]
controls <- as.numeric(rownames(edata) %in% rownames(edata)[selectEdata])
# as we do not include experimental conditions,  we use the supervised method
svobj <- sva(edata, mod = mod, mod0 = mod0, n.sv = n.sv,
             method = "supervised", controls = controls)
newV = NULL ## neccessary due to bug in the sva pacakge
fsvaobj <- fsva(edata, mod, svobj, newdat = NULL)

# get corrected data
edata_Adj <- fsvaobj$db

object_name = paste("peptides.matrix.sva", fraction,n.sv, sep = ".")
assign(object_name, edata_Adj)
file_name = paste(object_name, fun_name, "RData", sep = ".")
file_path = paste(output_dir, file_name, sep="/")
save(list = eval(object_name), file = file_path)

before = peptides.matrix
after  = get(eval(object_name))


pca = prcomp(t(before), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = data.frame(x_var, y_var, type="before")

scores = as.data.frame(pca$x[,1:5])
scores$type = "before"
scores$sample.id = rownames(scores)

pca = prcomp(t(after), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = rbind(annot,data.frame(x_var, y_var, type="after"))

scores = rbind(scores, data.frame(sample.id = rownames(scores), pca$x[,1:5], type = "after"))
scores$batch = factor(pheno$batch.exp.n[match(scores$sample.id, pheno$sample_name)])
scores.mix = scores[grepl(pattern="mix", ignore.case=T, x=rownames(scores)),]
scores$type = factor(scores$type, levels=c("before", "after"))

annot$text = paste(annot$x_var, annot$y_var)

library(cowplot)

p = ggplot(scores, aes(x=PC1, y=PC2)) + 
  geom_point(size=3, aes(col=batch) )+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(data=scores.mix, aes(x=PC1, y=PC2),size=3,col="black", shape=17) +
  geom_text(data = annot, aes(x=-50, y=-50, label=text)) +
  facet_wrap(~type, scales="fixed") + 
  theme_bw() +
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = rel(1.5)))

plots.list <- lappend(plots.list, p)

before.long = tbl_df(melt(before, id.vars="rownames"))
names(before.long) = c("EG.StrippedSequence", "R.Label", "signal")
before.long$aquisition_date =exp_metadata$aquisition_date[match(before.long$R.Label, exp_metadata$sample_name)]
before.long$batch = factor(exp_metadata$batch.exp.n[match(before.long$R.Label, exp_metadata$sample_name)])
before.long$category = "before"

after.long = tbl_df(melt(after, id.vars="rownames"))
names(after.long) = c("EG.StrippedSequence", "R.Label", "signal")
after.long$aquisition_date =exp_metadata$aquisition_date[match(after.long$R.Label, exp_metadata$sample_name)]
after.long$batch = factor(exp_metadata$batch.exp.n[match(after.long$R.Label, exp_metadata$sample_name)])
after.long$category = "after"


toPlot = rbind(before.long, after.long)
toPlot$category = factor(toPlot$category, levels=c("before", "after"))

set.seed(123)
tmp.peptides = unique(toPlot$EG.StrippedSequence)
tmp.selected = sample(tmp.peptides,5)
toPlot = toPlot[toPlot$EG.StrippedSequence %in% tmp.selected,]

toPlot$aquisition_date.str = as.Date(strptime(toPlot$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
toPlot.mix = toPlot[grep(x=toPlot$R.Label, pattern="mix", ignore.case=T),]
toPlot.mix.stats = toPlot.mix %>% group_by(EG.StrippedSequence, batch, category) %>% summarize(mean = mean(exp(signal)))

toPlot.wt  = toPlot[toPlot$R.Label %in% exp_metadata$sample_name[exp_metadata$ORF == "WT"],]

library(scales)
p = ggplot(toPlot, aes(x=aquisition_date.str, y=exp(signal), col=batch)) + 
  geom_point(size=3) +
  geom_point(data=toPlot.mix,aes(x=aquisition_date.str, y=exp(signal)),col="red") + #MIX
  geom_point(data=toPlot.wt, aes(x=aquisition_date.str, y=exp(signal)),col="blue") + #WT   
  scale_x_date(breaks = date_breaks("1 week"), minor_breaks = date_breaks("1 day"), labels=date_format("%m-%d"))+
  facet_grid(EG.StrippedSequence~category, scales="free")
plots.list <- lappend(plots.list, p)



file_name = paste("supplementary", fun_name, sep = ".")
file_path = paste(figures_dir, file_name, sep="/")

lapply(seq_along(plots.list) , 
       function(x) {
         
         tryCatch({
           p <- plots.list[[x]]
           scale = 1
           if (length(p$toScale) != 0 && p$toScale == T  ){
             scale = 2
           }
           ggplot2::ggsave(filename = paste(file_path, x , "pdf", sep = "."), device = NULL,
                           plot = p, width = 210 , height = 297, units = "mm", scale = scale)
           
           ggplot2::ggsave(filename = paste(file_path, x , "png", sep = "."), device = NULL,
                           plot = p, width = 210 , height = 297, dpi = 150, units = "mm", scale = scale)
           
         }, error = function(e) {
           message(paste("Plot", "x", "sucks!" ))
           return(NULL)
         }, finally = {
           message(paste("processed plot", x))
         })
         
         
       })
