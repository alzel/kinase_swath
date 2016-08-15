#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "sentinels_analysis"

load("./R/objects/sentinelsSRM._clean_.RData")
load("./R/objects/sentinels.table._clean_.RData")
load("./R/objects/sentinels.raw._load_.RData")
load("./R/objects/gene.annotations._load_.RData")
load("./R/objects/exp_metadata._clean_.RData")
orf2name = unique(data.frame(ORF = gene.annotations$V4,
                             sgd = gene.annotations$V5,
                             gene_name = gene.annotations$V6))
orf2name$ORF = as.character(orf2name$ORF)
orf2name$gene_name = as.character(orf2name$gene_name)
orf2name$gene_name[orf2name$gene_name ==""] = orf2name$ORF[orf2name$gene_name ==""]

sentinels.data = dplyr::select(sentinels.raw, R.Label, EG.StrippedSequence, EG.Label, FG.Id, EG.Qvalue, FG.TotalPeakArea, F.PeakArea )

stopifnot(length(unique(sentinels.data$EG.Label)) ==  length(unique(sentinels.data$FG.Id)))
  
fragments.data <- sentinels.data %>% group_by(R.Label, EG.Label) %>%
  mutate(FG.TotalPeakArea_new = sum(F.PeakArea)) %>% 
  group_by(R.Label, EG.Label, FG.Id ) %>% arrange()

#making precursor-fragment summaries
fragments.data <- fragments.data %>% group_by(FG.Id) %>%  mutate(qvalue.median = median(unique(EG.Qvalue)))

fragments.data.f <- fragments.data %>% filter(qvalue.median < 0.01)
peptides.peak_sums <- fragments.data.f %>% 
  group_by(R.Label, EG.Label, FG.Id) %>%
  dplyr::summarise(count = n(),  signal = FG.TotalPeakArea_new[1]) %>%
  group_by(R.Label, EG.Label) %>%
  dplyr::summarise(count = sum(count),  
                   signal = sum(signal))



load("./R/objects/exp_metadata._clean_.RData")
exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
exp_metadata$batch_kmeans = pam(exp_metadata$aquisition_date.str, 7)$clustering


peptides.df = dcast(peptides.peak_sums, formula=EG.Label~R.Label, value.var="signal")
sentinels.peptides.matrix = as.matrix(peptides.df[,-1])
rownames(sentinels.peptides.matrix) = peptides.df$EG.Label

pheno = exp_metadata[match(colnames(sentinels.peptides.matrix), exp_metadata$sample_name),]

mod = model.matrix(~as.factor(ORF), data=pheno)

sentinels.peptides.matrix.combat = ComBat(log(sentinels.peptides.matrix), batch=pheno$batch_kmeans, mod=mod, par.prior=T)
sentinels.peptides.matrix.quant.combat = ComBat(normalizeQuantiles(log(sentinels.peptides.matrix)), batch=pheno$batch_kmeans, mod=mod, par.prior=T)

pheno = exp_metadata[match(colnames(sentinels.peptides.matrix), exp_metadata$sample_name),]
mod = model.matrix(~1, data=pheno)
mod0 = mod

fraction = 0.5
n_genes = floor(nrow(sentinels.peptides.matrix)*fraction)
# choose ONE surrogate variable
n.sv  <- 1
# estimate surrogate variables

# Here, we use the lowly variable genes as control genes,
# i.e. we exclude the top 500 most variable genes from the
# estimation of the correction factors.
# This might prevent most of the biological signal from being regressed out.

edata = as.matrix(log(sentinels.peptides.matrix))
varsEdata <- rowSds(edata)/rowMeans(edata)

selectEdata <- order(varsEdata, decreasing = TRUE)[seq(n_genes, length(varsEdata))]
controls <- as.numeric(rownames(edata) %in% rownames(edata)[selectEdata])
# as we do not include experimental conditions,  we use the supervised method
svobj <- sva(sentinels.peptides.matrix, mod = mod, mod0 = mod0, n.sv = n.sv,
             method = "supervised", controls = controls)
newV = NULL ## neccessary due to bug in the sva pacakge
fsvaobj <- fsva(edata, mod, svobj, newdat = NULL)
# get corrected data
edata_Adj <- fsvaobj$db

object_name = paste("sentinels.peptides.matrix.sva", fraction,n.sv, sep = ".")
assign(object_name, edata_Adj)
file_name = paste(object_name, fun_name, "RData", sep = ".")
file_path = paste(output_dir, file_name, sep="/")
save(list = eval(object_name), file = file_path)




# file_name = "sentinels.peptides.matrix.combat.RData"
# file_path = paste(output_dir, file_name, sep="/")
# save(sentinels.peptides.matrix.combat,file=file_path) 
# 
# file_name = "sentinels.peptides.matrix.combat.RData"
# file_path = paste(output_dir, file_name, sep="/")
# save(sentinels.peptides.matrix.combat,file=file_path) 

file_name = "sentinels.peptides.matrix.quant.combat.RData"
file_path = paste(output_dir, file_name, sep="/")
save(sentinels.peptides.matrix.quant.combat,file=file_path) 



before = edata
after  = sentinels.peptides.matrix.sva.0.25.1


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
scores$batch_kmeans = factor(pheno$batch_kmeans[match(scores$sample.id, pheno$sample_name)])
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
  geom_text(data = annot, aes(x=-5, y=-5, label=text)) +
  facet_wrap(~type, scales="fixed") + 
  theme_bw() +
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = rel(1.5)))

#making proteins 

make_pseudo <- function(feature.matrix = sentinels.peptides.matrix.quant.combat) {
  
  feature.matrix = sentinels.peptides.matrix.quant.combat
  feature.matrix.long <- tbl_df(melt(feature.matrix))
  feature.matrix.long$ORF <- exp_metadata$ORF[match(feature.matrix.long$X2, exp_metadata$sample_name)]
  feature.matrix.long$type <- exp_metadata$type[match(feature.matrix.long$X2, exp_metadata$sample_name)]
  feature.matrix.long.f <- feature.matrix.long[grepl(feature.matrix.long$type, pattern = "Kinase|Wild Type"),]
  
  set.seed(123)
  tmp.features <- bind_rows(feature.matrix.long.f %>% # sampling
                              arrange(X1, ORF) %>% 
                              group_by(X1, ORF) %>%
                              sample_n(size=1) %>%
                              mutate(sample = 1),
                            feature.matrix.long.f %>% 
                              arrange(X1, ORF) %>% 
                              group_by(X1, ORF) %>%
                              sample_n(size=1) %>%
                              mutate(sample = 2),
                            feature.matrix.long.f %>% 
                              arrange(X1, ORF) %>% 
                              group_by(X1, ORF) %>%
                              sample_n(size=1) %>%
                              mutate(sample = 3))
  
  pseudo.strain <- tmp.features %>% group_by(X1, sample) %>% arrange(sample, X1) %>%
    summarise(X2 = paste("pseudo", unique(sample), sep = "_"),
              value = median(value, na.rm = T),
              ORF = "PSEUDO",
              type = "PSEUDO") %>%
    dcast(formula = "X1~X2", value.var = "value")
  
  
  pseudo_metadata <- data.frame(sample_name = c("pseudo_1", "pseudo_2", "pseudo_3"),
                                file_name = NA, 
                                ORF = "PSEUDO",
                                gene = "PSEUDO",
                                type  = "PSEUDO",
                                aquisition_date = NA,
                                batch.exp = NA,
                                batch_date = NA, 
                                batch.exp.n = NA)
  
  feature.matrix.df = data.frame(feature.matrix)
  feature.matrix.df$ORF <- rownames(feature.matrix)
  
  names(pseudo.strain)[1] = "ORF"
  
  tmp.df <- left_join(pseudo.strain, feature.matrix.df)
  tmp.matrix <- as.matrix(tmp.df[,-1])
  rownames(tmp.matrix) <- tmp.df$ORF
  return(list(tmp.matrix, pseudo_metadata))
  
}

tmp.res <- make_pseudo(feature.matrix = sentinels.peptides.matrix.quant.combat)
sentinels.peptides.matrix.quant.combat_pseudo <- tmp.res[[1]]

exp_metadata <- bind_rows(exp_metadata, tmp.res[[2]])


makeProteins = function(..., reference = "WT") {
  
  #... = peptides.matrix.combat.quant
  input_list = list(...)
  #input_list = list()
  #input_list[[1]] = peptides.matrix.combat.quant
  
  if (class(input_list[[1]]) != "matrix" & length(input_list) != 1) {
    stop("Matrix single argument has to be provided")
  }
  
  peptides.matrix_name = deparse(substitute(...))
  #peptides.matrix_name = deparse(substitute(peptides.matrix.combat.quant))
  proteins.matrix_name = paste("proteins", peptides.matrix_name, sep=".")
  
  if(grep(pattern="peptide", x=peptides.matrix_name)) {
    proteins.matrix_name = sub(x=peptides.matrix_name, pattern="peptide[s]?", "proteins")    
  }
  
  peptides.matrix = input_list[[1]]
  
  peptides.long = tbl_df(melt(peptides.matrix, id.vars="rownames"))
  names(peptides.long) = c("EG.Label", "R.Label", "signal")
  
  
  peptides.long.selected <- peptides.long
  peptides.long.selected <- peptides.long.selected %>% separate(col = "EG.Label", into = c("ORF", "modifiedSequence", "charge"), sep = "\\.")
  
  proteins.long = peptides.long.selected  %>% group_by(ORF, R.Label) %>% summarise(mean_signal = mean(exp(signal), na.rm=T),
                                                                                   sum_signal = log(sum(exp(signal), na.rm=T)))
  
  proteins.df = dcast(data=proteins.long, formula=ORF~R.Label, value.var="sum_signal")
  proteins.matrix = as.matrix(proteins.df[,-1])
  rownames(proteins.matrix) = proteins.df$ORF
  
  message(paste("number of proteins", nrow(proteins.matrix)))
  file_name = paste(proteins.matrix_name, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  
  assign(eval(proteins.matrix_name), proteins.matrix)
  save(list=eval(proteins.matrix_name), file=file_path)
  
  ## -- protein fold-changes ----
  
  pheno = exp_metadata[match(colnames(proteins.matrix), exp_metadata$sample_name),]
  pheno = droplevels(filter(pheno, type != "Standard Mix"))
  
  proteins.matrix.f = proteins.matrix[,match(pheno$sample_name, colnames(proteins.matrix))]
  
  X = model.matrix(~pheno$ORF + 0)
  colnames(X) = unique(pheno$ORF)
  
  #reference = "WT"
  matrix = proteins.matrix.f
  
  lm.fit_model = lmFit(matrix, X)
  ph = unique(as.character(pheno$ORF))
  contrasts = paste0( ph[ph !=reference] ,"-", reference)  
  
  mc = makeContrasts(contrasts=contrasts, levels=X)    
  c.fit = contrasts.fit(lm.fit_model, mc)
  eb = eBayes(c.fit)
  
  
  folds = rowFolds(data=exp(matrix), groups=pheno$ORF, reference=reference)
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
  
  proteins.FC$KO = sub(x = proteins.FC$contrasts, pattern=paste("(.*?)-", reference, sep=""), replacement="\\1")
  proteins.FC$reference = reference
  
  proteins.FC_name = paste(proteins.matrix_name, "FC", sep=".")
  
  file_name = paste(proteins.FC_name, reference,  "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  
  assign(eval(proteins.FC_name), proteins.FC)
  save(list=eval(proteins.FC_name), file=file_path)
  
}

makeProteins(sentinels.peptides.matrix.quant.combat, reference =  "WT")
makeProteins(sentinels.peptides.matrix.quant.combat_pseudo, reference = "PSEUDO")
makeProteins(sentinels.peptides.matrix.sva.0.5.1, reference =  "WT")




### ---- plotting sentinels ------
load("./R/objects/sentinels.proteins.matrix.quant.combat.FC.RData")
load("./R/objects/sentinels.proteins.matrix.quant.combat.RData")

load("./R/objects/proteins.matrix.combat.quant.FC.RData")
load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/iMM904._load_.RData")

protein.matrix = proteins.matrix.combat.quant
proteins.FC = sentinels.proteins.matrix.quant.combat.FC

proteins.FC$gene_name <- orf2name$gene_name[match(proteins.FC$ORF, orf2name$ORF)]
reference = unique(as.character(proteins.FC$reference))
proteins.FC.f = proteins.FC[proteins.FC$KO %in% unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"])),]

pval_thr = 0.01
set.seed(123)
FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)


sentinels.table$short_info <- sub(pattern = "Marker for ", replacement = "", x = sentinels.table$Sentinel)
sentinels.table$short_info <- trimws(sentinels.table$short_info)
sentinels.table$ORF.ID <- trim(sentinels.table$ORF.ID)

proteins.FC.f$KO.gene <- exp_metadata$gene[match(proteins.FC.f$KO, exp_metadata$ORF)]

toPlot <- proteins.FC.f %>% 
  filter(ORF != "iRT", abs(logFC) > FC_thr, p.value_BH < pval_thr) %>%
  group_by(KO) %>%
  mutate(sign = n()) %>%
  group_by(ORF) %>%
  mutate(popularity = n()) %>%
  filter(popularity >=10, sign >= 15)

sentinels_desc_collapsed = sentinels.table %>% dplyr::select(ORF.ID, short_info) %>% distinct() %>%
  filter(ORF.ID %in% toPlot$ORF) %>%
  group_by(short_info) %>%
  mutate(info_n = n()) %>%
  group_by(ORF.ID) %>%
  summarise(description_full = paste(sort(short_info),  collapse = "|"),
            description = ifelse(length(short_info) > 2, "multiple", paste(sort(short_info),  collapse = "|")),
            popular_info = paste(max(info_n), short_info[which.max(info_n)], collapse = ":" )) %>% 
  ungroup %>% 
  arrange(popular_info)

sentinels_desc_single = sentinels.table %>% dplyr::select(ORF.ID, short_info) %>%  filter(ORF.ID %in% toPlot$ORF) %>% distinct()


lb = -2.5
ub = 2.5

my_breaks <- seq(lb, ub, 0.5)
my_breaks[5] <- FC_thr
my_breaks[7] <- -FC_thr
my_levels = levels(cut(toPlot$logFC, breaks = my_breaks))

my_colours1 <- brewer.pal(name = "Reds", n = 4)
my_colours2 <- rev(brewer.pal(name = "Blues", n = 4))
my_colours = c(my_colours2, my_colours1)

toPlot$my_fill <- as.character(cut(toPlot$logFC, breaks = my_breaks))
toPlot$my_fill[which(is.na(toPlot$my_fill))] <- ifelse(toPlot$logFC[which(is.na(toPlot$my_fill))] <= lb, lb, ub)
toPlot$my_fill <- factor(toPlot$my_fill, levels = my_levels)

toPlot$info <- sentinels_desc_collapsed$description[match(toPlot$ORF, sentinels_desc_collapsed$ORF.ID)]
toPlot$sentinel <- paste(toPlot$info, toPlot$gene_name, sep = "::")

toPlot$popular_info <- sentinels_desc_collapsed$popular_info[match(toPlot$ORF, sentinels_desc_collapsed$ORF.ID)]
toPlot$sentinel_popular <- paste(toPlot$popular_info, toPlot$gene_name, sep = "::")

# sorted = sort(unique(toPlot$sentinel))
# sorted <- c(sorted[grep("multiple", x = sorted, invert = T)], sorted[grep("multiple", x = sorted)])
#toPlot$sentinel <- factor(toPlot$sentinel, levels = rev(sorted))
sorted <- sort(unique(toPlot$sentinel_popular))
sorted <- sorted[c(grep(12, sorted,invert =T), grep(12, sorted))]
toPlot$sentinel_popular <- factor(toPlot$sentinel_popular, levels = sorted)

p.heatmap_sentinels <- ggplot(toPlot) +  
  geom_tile(aes(x = KO.gene, y = sentinel_popular, fill = my_fill ), colour="grey") +
  #   scale_fill_gradient2(low="#1F78B4",high="#E31A1C",mid ="white",
  #                        breaks = seq(-0.75, 0.75, 0.25),
  #                        midpoint=0)  +
  scale_fill_manual(values = my_colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = c(0.2, 0.8) ) +
  labs(x="", y = "")


### ---- sentinel coverage ---- ####

sum(unique(sentinelsSRM$ORF) %in% rownames(sentinels.proteins.matrix.quant.combat))/length(unique(sentinelsSRM$ORF))

#sum((sentinels.table %>% filter(Sentinel.Grade == "A") %>% select(ORF.ID) %>% distinct())$ORF.ID %in% rownames(proteins.matrix.combat.quant))
  



