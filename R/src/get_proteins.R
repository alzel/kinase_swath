rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "get_proteins"

load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/peptides.matrix.combat.RData")
load("./R/objects/peptides.matrix.combat.quant.RData")
load("./R/objects/peptides.cor.stats.top.RData")
load("./R/objects/protein_annotations._load_.RData")

#write.table(x=rownames(peptides.matrix.combat), file="peptide.txt", quote=F, row.names=F, col.names=F)

makeProteins = function(...) {
  
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
    
  
  unique_peptides = tbl_df(droplevels(filter(protein_annotations, specificity == "unique") %>% dplyr::select(strippedSequence, SystName)))
    
  peptides.matrix = input_list[[1]]
  peptides.long = tbl_df(melt(peptides.matrix, id.vars="rownames"))
  names(peptides.long) = c("EG.StrippedSequence", "R.Label", "signal")
  
  peptides.selected = tbl_df(droplevels(filter(peptides.cor.stats.top, top == "3")))
  peptides.selected = peptides.selected[peptides.selected$EG.StrippedSequence %in% unique_peptides$strippedSequence,]
  
  peptides.long.selected     = peptides.long[peptides.long$EG.StrippedSequence %in% peptides.selected$EG.StrippedSequence,]
  peptides.long.selected$ORF = unique_peptides$SystName[match(peptides.long.selected$EG.StrippedSequence, unique_peptides$strippedSequence)]
  
  proteins.long = peptides.long.selected  %>% group_by(ORF, R.Label) %>% summarise(mean_signal = mean(exp(signal), na.rm=T),
                                                                                   sum_signal = log(sum(exp(signal), na.rm=T)))
  
  proteins.df = dcast(data=proteins.long, formula=ORF~R.Label, value.var="sum_signal")
  proteins.matrix = as.matrix(proteins.df[,-1])
  rownames(proteins.matrix) = proteins.df$ORF
  
  
  file_name = paste(proteins.matrix_name, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  
  assign(eval(proteins.matrix_name), proteins.matrix)
  save(list=eval(proteins.matrix_name), file=file_path)
  
  ## -- protein fold-changes ----
  
  pheno = exp_metadata[match(colnames(proteins.matrix), exp_metadata$sample_name),]
  pheno = droplevels(filter(pheno, type != "Standard Mix"))
  
  proteins.matrix.f = proteins.matrix[,match(pheno$sample_name, colnames(proteins.matrix))]
  
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
  
  file_name = paste(proteins.FC_name, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  
  assign(eval(proteins.FC_name), proteins.FC)
  save(list=eval(proteins.FC_name), file=file_path)
  
}



#peptides.long = tbl_df(melt(peptides.matrix.combat, id.vars="rownames"))
makeProteins(peptides.matrix.combat)
makeProteins(peptides.matrix.combat.quant)



# ## peptides/proteins differential expression
# 
# ## -- protein fold-changes ----
# proteins.matrix = proteins.matrix.combat.quant
# 
# pheno = exp_metadata[match(colnames(proteins.matrix), exp_metadata$sample_name),]
# pheno = droplevels(filter(pheno, type != "Standard Mix"))
# 
# proteins.matrix.f = proteins.matrix[,match(pheno$sample_name, colnames(proteins.matrix))]
# 
# X = model.matrix(~pheno$ORF + 0)
# colnames(X) = levels(pheno$ORF)
# 
# reference = "WT"
# matrix = proteins.matrix.f
# 
# lm.fit_model = lmFit(matrix, X)
# ph = unique(as.character(pheno$ORF))
# contrasts = paste0( ph[ph !=reference] ,"-", reference)  
# 
# mc = makeContrasts(contrasts=contrasts, levels=X)    
# c.fit = contrasts.fit(lm.fit_model, mc)
# eb = eBayes(c.fit)
# 
# 
# folds = rowFolds(data=exp(matrix), groups=pheno$ORF, reference=reference)
# folds = log(folds, 2)
# 
# #folds_tmp = melt(as.matrix(folds), id.vars="row.names")
# 
# #merging results
# folds_tmp = melt(eb$coefficients, id.vars="row.names")
# #folds_tmp$contrasts = factor(paste(folds_tmp$contrasts, "-", reference, sep=""))
# pvals_tmp = melt(eb$p.value, id.vars="row.names")
# 
# 
# names(folds_tmp) = c("ORF", "contrasts", "logFC")
# names(pvals_tmp) = c("ORF", "contrasts", "p.value")
# 
# folds_tmp$contrasts = factor(folds_tmp$contrasts)
# pvals_tmp$contrasts = factor(pvals_tmp$contrasts)
# 
# 
# proteins.FC = merge(folds_tmp, pvals_tmp, all=T,
#                     by=c("ORF", "contrasts"))
# 
# ##multiple testing correction
# proteins.FC$p.value_BH = p.adjust(proteins.FC$p.value, method="BH")
# proteins.FC$p.value_bonferroni = p.adjust(proteins.FC$p.value, method="bonferroni")
# 
# proteins.FC$KO = sub(x = proteins.FC$contrasts, pattern=paste("(.*?)-", reference, sep=""), replacement="\\1")
# proteins.FC$reference = reference
# 
# proteins.FC.combat.quant = proteins.FC
# file_name = "proteins.FC.combat.quant.RData"
# file_path = paste(output_dir, file_name, sep="/")
# save(proteins.FC.combat.quant, file=file_path)
# 





