#!/usr/bin/env Rscript
# linear models of metabolism for TCA metabolites from selected strains

rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "create_datasets"


load("./R/objects/metabolite_metadata._clean_.RData")
load("./R/objects/exp_metadata._clean_.RData")
load("./R/objects/metabolites.data._clean_.RData")
load("./R/objects/protein_annotations._load_.RData")
load("./R/objects/proteins.matrix.combat.quant.RData")
load("./R/objects/proteins.matrix.combat.RData")
load("./R/objects/exp_metadata._clean_.RData")

load("./R/objects/aa_metadata._clean_.RData")
load("./R/objects/aa.data._clean_.RData")

#load("./R/objects/yeast.model._load_.RData")
library("Amelia")


load("./R/objects/iMM904._load_.RData")
load("./R/objects/metabolite2iMM904._load_.RData")


orf2name = unique(data.frame(ORF = protein_annotations$SystName, 
                             gene_name = protein_annotations$sgdName))



getResPred = function(response.matrix, predictors.matrix, i,  B, order=NULL, include.metabolites = F) {
  #B - bipartite graph with type == 1 of metabolites 
  
  stopifnot(!any(is.null(predictors.matrix)| is.null(response.matrix) | !is.bipartite(B)))
  
  if(is.null(order)) {
    order = 1  
  }
  
  current_nodes = as.character(metabolite2iMM904$model_name[metabolite2iMM904$id == i])
  stopifnot(i %in% colnames(response.matrix))
  
  if (length(current_nodes) == 0) {
    message(paste("No ", i, "found in network"))
    return(-1)
  }
  
  
  SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=order, nodes=current_nodes))))  
  #SUB = induced.subgraph(B, unique(unlist(neighborhood(B, order=2, nodes=current_nodes))))  
  
  genes.pred = unique(V(SUB)$name[V(SUB)$type == 0])
  genes.pred = colnames(predictors.matrix)[colnames(predictors.matrix) %in% genes.pred]
  
  if (length(genes.pred) == 0) {
    return(-1)
  }
  additional.present = integer(0)
  
  if (order > 1 & include.metabolites) {
    metabolites.pred = V(SUB)$name[V(SUB)$type == 1]
    metabolites.pred = metabolites.pred[!(metabolites.pred %in% current_nodes)] #removing response
    additional.pred = as.vector(metabolite2iMM904$id[metabolite2iMM904$model_name %in% metabolites.pred])
    additional.present = colnames(response.matrix)[which(colnames(response.matrix) %in% additional.pred)]
  }
  
  if (length(additional.present) == 0) {
    tmp.df = as.data.frame(cbind(response.matrix[,i], predictors.matrix[,genes.pred]))
    colnames(tmp.df) = c(i, genes.pred)
    return(tmp.df)
  } else {
    tmp.df = as.data.frame(cbind(response.matrix[,i], response.matrix[,additional.present], predictors.matrix[,genes.pred]))
    colnames(tmp.df) = c(i,additional.present, genes.pred)
    return(tmp.df)
  }
}


## -- selected kinases TCA ----
clean_data_TCA = function(imputed = F) {
    
  load("./R/objects/metabolitesTCA_metadata._clean_.RData")
  load("./R/objects/metabolitesTCA.data._clean_.RData")
  
  ## -- TCA metabolites batch effects ----
  metabolitesTCA.data$date = metabolitesTCA_metadata$measure_date[match(metabolitesTCA.data$sample_id, metabolitesTCA_metadata$sample_id)]
  metabolitesTCA.data$batch = metabolitesTCA_metadata$measure_batch[match(metabolitesTCA.data$sample_id, metabolitesTCA_metadata$sample_id)]
  
  metabolitesTCA.data = droplevels(metabolitesTCA.data[grep(pattern="_[123]+$", x=metabolitesTCA.data$sample_id, perl=T),])
  metabolitesTCA.data = droplevels(filter(metabolitesTCA.data, variable != "Glu"))
  
  wt_points = metabolitesTCA.data[grep(x=metabolitesTCA.data$sample_id, pattern="WT", ignore.case=T),]
  
#   library(scales)
#   p1 = ggplot(metabolitesTCA.data, aes(x=batch, y=value)) +
#     geom_point() +
#     geom_point(data=wt_points, aes(x=batch, y=value), col="red") +
#     #geom_text(data=wt_points, hjust=1, vjust = 0, aes(x=batch, y=value, label=sample_id), col="red") +
#     ggtitle("Before correction") +
#     facet_wrap(~variable, scales="free")  
  
  metabolitesTCA.df = dcast(metabolitesTCA.data, formula=variable~sample_id)
  metabolitesTCA.matrix = as.matrix(metabolitesTCA.df[,-1])
  rownames(metabolitesTCA.matrix) = metabolitesTCA.df$variable
  
  
  metabolitesTCA.matrix.t = t(metabolitesTCA.matrix)
  
  set.seed(123)
  metabolitesTCA.imputed = amelia(metabolitesTCA.matrix.t, logs=colnames(metabolitesTCA.matrix.t), m=5)
  metabolitesTCA.imputed.matrix = Reduce("+",metabolitesTCA.imputed$imputations)/length(metabolitesTCA.imputed$imputations)
  
  
  pheno = metabolitesTCA_metadata[metabolitesTCA_metadata$sample_id %in% rownames(metabolitesTCA.imputed.matrix),]
  pheno = as.data.frame(droplevels(pheno[match(rownames(metabolitesTCA.imputed.matrix), pheno$sample_id),]))
  pheno$measure_batch = factor(pheno$measure_batch)
  
  mod = model.matrix(~as.factor(ORF), data=pheno)
  metabolitesTCA.matrix.combat = ComBat(log(t(metabolitesTCA.imputed.matrix)), batch=pheno$measure_batch, mod=mod, par.prior=T)
  metabolitesTCA.matrix.combat.long = melt(t(metabolitesTCA.matrix.combat), id.vars="rownames")
  names(metabolitesTCA.matrix.combat.long) = c("sample_id","variable", "value")
  
  metabolitesTCA.matrix.combat.long$batch = metabolitesTCA_metadata$measure_batch[match(metabolitesTCA.matrix.combat.long$sample_id, metabolitesTCA_metadata$sample_id)]
  
  #wt_points.combat = metabolitesTCA.matrix.combat.long[grep(x=metabolitesTCA.matrix.combat.long$sample_id, pattern="WT", ignore.case=T),]
  
#   p2 = ggplot(metabolitesTCA.matrix.combat.long, aes(x=batch, y=value)) +
#     geom_point() +
#     geom_point(data=wt_points.combat, aes(x=batch, y=value), col="red") +
#     #geom_text(data=wt_points, hjust=1, vjust = 0, aes(x=batch, y=value, label=sample_id), col="red") +
#     ggtitle("After correction") +
#     facet_wrap(~variable, scales="free")
#   g = arrangeGrob(p1,p2, ncol=2)
#   plots.list = lappend(plots.list, g)
  
  ## -- cleaning TCA data for linear models ----
  metabolitesTCA.long = merge(metabolitesTCA.matrix.combat.long, metabolitesTCA.data, by=c("sample_id", "variable", "batch"), suffixes=c(".combat", ".raw"))
  if (!imputed) {
    metabolitesTCA.long$value.combat[is.na(metabolitesTCA.long$value.raw)] = NA
  }
  
  proteins.matrix = proteins.matrix.combat
  
  proteins.long = melt(proteins.matrix, id.vars="rownames")
  names(proteins.long) = c("ORF", "R.Label", "signal")
  proteins.long$KO = exp_metadata$ORF[match(proteins.long$R.Label, exp_metadata$sample_name)]
  
  set.seed(123)
  metabolitesTCA_metadata.wt = metabolitesTCA_metadata %>% filter(ORF== "WT") %>% group_by(ORF) %>% sample_n(3)
  metabolitesTCA_metadata.NOwt = metabolitesTCA_metadata %>% filter(ORF != "WT") %>% group_by(sample, replicate) %>% distinct()
  
  TCAmetadata = rbind(metabolitesTCA_metadata.wt, metabolitesTCA_metadata.NOwt) %>% 
    group_by(measure_batch, ORF) %>% 
    mutate(prot.sample_id = exp_metadata$sample_name[base::sample(x=which(exp_metadata$ORF %in% ORF),replace=F)[1:length(ORF)]],
           aa.sample_id = pheno$sample_id[base::sample(x=which(pheno$ORF %in% ORF),replace=F)[1:length(ORF)]])
  TCAmetadata = droplevels(TCAmetadata)
  
  metabolitesTCA.final.df = dcast(metabolitesTCA.long, formula=sample_id~variable, value.var="value.combat")
  metabolitesTCA.final.matrix = as.matrix(metabolitesTCA.final.df[,-1])
  rownames(metabolitesTCA.final.matrix) = metabolitesTCA.final.df$sample_id
  metabolitesTCA.final.matrix = metabolitesTCA.final.matrix[rownames(metabolitesTCA.final.matrix) %in% TCAmetadata$sample_id,]
  rownames(metabolitesTCA.final.matrix) = TCAmetadata$prot.sample_id[match(rownames(metabolitesTCA.final.matrix), as.vector(TCAmetadata$sample_id))]
  
  #TCA proteins
  proteins.matrix.combat.t = t(proteins.matrix)
  proteins.matrix.combat.t.f = proteins.matrix.combat.t[rownames(proteins.matrix.combat.t) %in% as.vector(TCAmetadata$prot.sample_id),]
  
  both.present = intersect(rownames(proteins.matrix.combat.t.f), rownames(metabolitesTCA.final.matrix))
  
  proteins.matrix.combat.t.f = proteins.matrix.combat.t.f[rownames(proteins.matrix.combat.t.f) %in% both.present,]
  metabolitesTCA.final.matrix = metabolitesTCA.final.matrix[rownames(metabolitesTCA.final.matrix) %in% both.present,]
  
  proteinsTCA.present = proteins.matrix.combat.t.f[match(both.present, rownames(proteins.matrix.combat.t.f)),]
  metabolitesTCA.present = metabolitesTCA.final.matrix[match(both.present, rownames(metabolitesTCA.final.matrix)),]

  return(list(proteins = exp(proteinsTCA.present),
              metabolites = exp(metabolitesTCA.present),
              proteins.log = proteinsTCA.present,
              proteins.log.quant = normalizeQuantiles(proteinsTCA.present)))
}


clean_data_AA = function() {
      
  aa.data$date = aa_metadata$date[match(aa.data$sample_id, aa_metadata$sample_id)]
  aa.data$batch = aa_metadata$Batch[match(aa.data$sample_id, aa_metadata$sample_id)]
  
  aa_metadata.f = droplevels(aa_metadata[grep(pattern="^A_|cryoStar|cryo\\+|C_WT_[ABC]+.*?", x=aa_metadata$sample_id, invert=T),])
  
  aa.data = aa.data[aa.data$sample_id %in% aa_metadata.f$sample_id,]
#   qc_points = aa.data[grep(x=aa.data$sample_id, pattern="WT", ignore.case=T),]
  
  
#   p1 = ggplot(aa.data, aes(x=batch, y=value)) +
#     geom_point() +
#     geom_point(data=qc_points, aes(x=batch, y=value), col="red") +
#     geom_text(data=qc_points, aes(x=batch, y=value, label=sample_id), col="red") +
#     ggtitle("Before correction")+
#     #scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels=date_format("%m-%d")) +
#     facet_wrap(~variable, scales="free")
  
  pheno = aa_metadata.f
  
  aa.data.df = dcast(aa.data, sample_id~variable, value.var="value")
  aa.data.matrix = as.matrix(aa.data.df[,-1])
  rownames(aa.data.matrix) = aa.data.df$sample_id
  
  
  #correcting for batch effects AA data
  aa.data.matrix = aa.data.matrix[complete.cases(aa.data.matrix),]
  pheno = droplevels(pheno[match(rownames(aa.data.matrix), pheno$sample_id),])
  aa.data.matrix.t = t(aa.data.matrix)
  
  mod = model.matrix(~ORF, data=pheno)
  stopifnot(length(pheno$Batch) == ncol(aa.data.matrix.t))
  
  aa.data.matrix.combat = ComBat(log(aa.data.matrix.t), batch=pheno$Batch, mod=mod, par.prior=T)
  aa.data.combat.long = melt(t(aa.data.matrix.combat), id.vars="rownames")
  names(aa.data.combat.long) = c("sample_id","variable","value" )
  
  
  aa.data.combat.long$batch = aa_metadata.f$Batch[match(aa.data.combat.long$sample_id, aa_metadata.f$sample_id)]
  
  #qc_points.combat = aa.data.combat.long[grep(x=aa.data.combat.long$sample_id, pattern="WT", ignore.case=T),]
  
#   p2 = ggplot(aa.data.combat.long, aes(x=batch, y=value)) +
#     geom_point() +
#     geom_point(data=qc_points.combat, aes(x=batch, y=value), col="red") +
#     ggtitle("After correction") +
#     facet_wrap(~variable, scales="free")
#   
#   g = arrangeGrob(p1,p2, ncol=1)
#   plots.list = lappend(plots.list, g)
  
  set.seed(123)
  aa_metadata.wt = aa_metadata.f %>% filter(ORF== "WT") %>% group_by(ORF) %>% sample_n(3)
  aa_metadata.NOwt = aa_metadata.f %>% filter(ORF != "WT") %>% group_by(ORF, Replicate) %>% distinct()
  
  AAmetadata = rbind(aa_metadata.wt, aa_metadata.NOwt) %>% 
    group_by(Batch, ORF) %>% 
    mutate(prot.sample_id = exp_metadata$sample_name[base::sample(x=which(exp_metadata$ORF %in% ORF),replace=F)[1:length(ORF)]])
  
  #aa.sample_id = pheno$sample_id[base::sample(x=which(pheno$ORF %in% ORF),replace=F)[1:length(ORF)]])
  
  aa.data.combat.long$ORF = aa_metadata$ORF[match(aa.data.combat.long$sample_id, aa_metadata$sample_id)]
  aa.data.combat.final.df = dcast(aa.data.combat.long, formula = sample_id~variable)
  aa.data.combat.final.df = aa.data.combat.final.df[aa.data.combat.final.df$sample_id %in% AAmetadata$sample_id,]
  aa.data.combat.final.matrix = as.matrix(aa.data.combat.final.df[,-1])
  
  rownames(aa.data.combat.final.matrix) = AAmetadata$prot.sample_id[match(aa.data.combat.final.df[,1], AAmetadata$sample_id)]
  aa.data.combat.final.matrix = aa.data.combat.final.matrix[!is.na(rownames(aa.data.combat.final.matrix)),]
  
  ##AA proteins
  proteins.matrix = proteins.matrix.combat
  proteins.matrix.combat.t = t(proteins.matrix)
  proteins.matrix.combat.t.f = proteins.matrix.combat.t[rownames(proteins.matrix.combat.t) %in% as.vector(AAmetadata$prot.sample_id),]
  
  both.present = intersect(rownames(proteins.matrix.combat.t.f), rownames(aa.data.combat.final.matrix))
  
  proteins.matrix.combat.t.f = proteins.matrix.combat.t.f[rownames(proteins.matrix.combat.t.f) %in% both.present,]
  aa.data.combat.final.matrix = aa.data.combat.final.matrix[rownames(aa.data.combat.final.matrix) %in% both.present,]
  
  proteinsAA.present = proteins.matrix.combat.t.f[match(both.present, rownames(proteins.matrix.combat.t.f)),]
  
  metabolitesAA.present = aa.data.combat.final.matrix[match(both.present, rownames(aa.data.combat.final.matrix)),]  
  metabolitesAA.present = metabolitesAA.present[,-which(colnames(metabolitesAA.present) == "isoleucine")]
  
  return(list(proteins = exp(proteinsAA.present),
              metabolites = exp(metabolitesAA.present),
              proteins.log = proteinsAA.present,
              proteins.log.quant = normalizeQuantiles(proteinsAA.present)))
  
}

clean_data_PPP_AA = function(imputed = F) {
  
  metabolites.data$batch = factor(droplevels(metabolite_metadata$batch[match(metabolites.data$sample_id, metabolite_metadata$sample_id)]))
  
   WTs = as.character(droplevels(metabolites.data$sample_id[grep(pattern="wt", x=metabolite_metadata$ORF, ignore.case=T)]))
#   WTs.subset = metabolites.data[metabolites.data$sample_id %in% WTs,]
#   
#   p = ggplot(metabolites.data, aes(x=batch, y=value)) +
#     geom_boxplot() +
#     geom_point(data=WTs.subset, aes(x=batch, y=value), col="red") +
#     facet_wrap(~variable, scales="free")
#   plots.list = lappend(plots.list, p)
  
  
  #imputing missing values
  #TODO: check validity of imputation
  
  metabolites.df = dcast(metabolites.data, formula=variable~sample_id)
  metabolites.matrix = as.matrix(metabolites.df[,-1])
  rownames(metabolites.matrix) = metabolites.df$variable
  
  

  metabolites.matrix.t = t(metabolites.matrix)
  set.seed(123)
  metabolites.imputed = amelia(metabolites.matrix.t, logs=colnames(metabolites.matrix.t), m=5)
  metabolites.imputed.matrix = Reduce("+",metabolites.imputed$imputations)/length(metabolites.imputed$imputations)
  
  
  metabolites.imputed.long = melt(metabolites.imputed.matrix, id.vars="rownames")
  names(metabolites.imputed.long) = c("sample_id","variable","value" )
  
  metabolites.imputed.long$batch = factor(droplevels(metabolite_metadata$batch[match(metabolites.imputed.long$sample_id, metabolite_metadata$sample_id)]))
  
  
  # -- batch correction for metabolite data ----
  
   WTs.batch = unique(as.character(metabolite_metadata$batch[grep(pattern="wt", x=metabolite_metadata$ORF, ignore.case=T)]))
   WTs.imputed.subset = metabolites.imputed.long[metabolites.imputed.long$sample_id %in% WTs,]
   WTs.batch = WTs.batch[WTs.batch != 8]
#   
#   
#   p1 = ggplot(metabolites.imputed.long, aes(x=batch, y=value)) +
#     geom_boxplot() +
#     geom_point(data=WTs.imputed.subset, aes(x=batch, y=value), col="red") +
#     ggtitle("Before correction") +
#     facet_wrap(~variable, scales="free")
#   
  
  pheno = droplevels(metabolite_metadata[metabolite_metadata$batch %in% WTs.batch,])
  metabolites.imputed.matrix.t = t(metabolites.imputed.matrix)
  metabolites.imputed.matrix.t.f = metabolites.imputed.matrix.t[,match(pheno$sample_id, colnames(metabolites.imputed.matrix.t))]
  
  mod = model.matrix(~as.factor(ORF), data=pheno)
  
  metabolites.matrix.combat = ComBat(log(metabolites.imputed.matrix.t.f), batch=pheno$batch, mod=mod, par.prior=F)
  
  metabolites.matrix.combat.long = melt(t(metabolites.matrix.combat), id.vars="rownames")
  names(metabolites.matrix.combat.long) = c("sample_id","variable", "value")
  
  metabolites.matrix.combat.long$batch = droplevels(metabolite_metadata$batch[match(metabolites.matrix.combat.long$sample_id, metabolite_metadata$sample_id)])
#   WTs.combat.subset = droplevels(metabolites.matrix.combat.long[metabolites.matrix.combat.long$sample_id %in% WTs,])
#   
#   p2 = ggplot(metabolites.matrix.combat.long, aes(x=batch, y=value)) +
#     geom_boxplot() +
#     geom_point(data=WTs.combat.subset, aes(x=batch, y=value), col="red") +
#     ggtitle("After correction") +
#     facet_wrap(~variable, scales="free")
#   
#   g = arrangeGrob(p1,p2, ncol=1)
#   plots.list = lappend(plots.list, g)
#   
  metabolites.matrix.combat.long$ORF = pheno$ORF[match(metabolites.matrix.combat.long$sample_id, pheno$sample_id)]
  metabolites.long.mean = tbl_df(metabolites.matrix.combat.long) %>% group_by(variable, ORF) %>% summarize(mean = mean(value))
  
  metabolites.long.merged = merge(metabolites.matrix.combat.long, metabolites.data, by=c("sample_id", "variable", "batch"), suffixes=c(".combat", ".raw"))
  if(!imputed) {
    metabolites.long.merged$value.combat[is.na(metabolites.long.merged$value.raw)] = NA
  }
    
  metabolites.long.mean.models = tbl_df(metabolites.long.merged) %>% group_by(variable, ORF) %>% summarize(mean = mean(value.combat))
  
  
  proteins.matrix = proteins.matrix.combat
  
  proteins.long = melt(proteins.matrix, id.vars="rownames")
  names(proteins.long) = c("EG.StrippedSequence", "R.Label", "signal")
  proteins.long$ORF = exp_metadata$ORF[match(proteins.long$R.Label, exp_metadata$sample_name)]
  
  proteins.long.mean = tbl_df(proteins.long) %>% group_by(EG.StrippedSequence, ORF) %>% summarize(mean = mean(signal))
  
  
  proteins.mean.df = dcast(proteins.long.mean, formula=EG.StrippedSequence~ORF, value.var="mean")
  proteins.mean.matrix = as.matrix(proteins.mean.df[,-1])
  rownames(proteins.mean.matrix) = as.matrix(proteins.mean.df$EG.StrippedSequence)
  
  
  metabolites.mean.models.df = dcast(metabolites.long.mean.models, formula=variable~ORF, value.var="mean")
  metabolites.mean.matrix = as.matrix(metabolites.mean.models.df[,-1])
  rownames(metabolites.mean.matrix) = metabolites.mean.models.df$variable

  metabolites.binded = metabolites.long.mean.models
  metabolites.all.mean.df = dcast(metabolites.binded, formula=variable~ORF, value.var="mean")
  
  metabolites.all.mean.matrix = as.matrix(metabolites.all.mean.df[,-1])
  rownames(metabolites.all.mean.matrix) = metabolites.all.mean.df$variable
  
  both.present = intersect(colnames(proteins.mean.matrix), colnames(metabolites.all.mean.matrix))
  
  proteins.mean.matrix.present = t(proteins.mean.matrix[,match(both.present, colnames(proteins.mean.matrix))])
  metabolites.all.mean.matrix.present = t(metabolites.all.mean.matrix[,match(both.present, colnames(metabolites.all.mean.matrix))])
  
  return(list(proteins = exp(proteins.mean.matrix.present),
              metabolites = exp(metabolites.all.mean.matrix.present),
              proteins.log = proteins.mean.matrix.present,
              proteins.log.quant = normalizeQuantiles(proteins.mean.matrix.present)))
}



createDataset = function(response.matrix, predictors.matrix, order, include.metabolites, output_dir, preffix="dataset") {
  
#   response.matrix = dataAA$metabolites
#   predictors.matrix = dataAA$proteins
#   order = 1
#   include.metabolites = F
#   preffix = "test"  
#   output_dir = "./results/2015-08-03/data.AA/"
  
  stopifnot(!any(is.null(predictors.matrix)| is.null(response.matrix) | is.null(output_dir) | is.null(preffix)))
  
   
  yeast.model = iMM904
  yeast.model = yeast.model[grep("t", yeast.model$reaction, invert=T),] #removing all tranporters
  edgelist = unique(droplevels(na.omit(subset(yeast.model, metabolite != "h"  & metabolite !="h2o" , select = c("metabolite", "gene")))))
  
  B <- graph.data.frame(edgelist)
  V(B)$type <- V(B)$name %in% edgelist$metabolite
  stopifnot(is.bipartite(B))
  
  for(ji in 1:ncol(response.matrix)) {
    
    i <<- colnames(response.matrix)[ji]
    
        
    tmp.df = getResPred(response.matrix=response.matrix, 
                        predictors.matrix=predictors.matrix, 
                        B=B, i=i, order=order, include.metabolites=include.metabolites)
    
    if (length(tmp.df) == 1 && tmp.df  == -1) {
      next
    }
    
    if (ncol(tmp.df) == 2 ) {
      next
    }
    
    file_name = paste(preffix, i, order, as.numeric(include.metabolites),  "RData", sep=".")
    file_path = paste(output_dir, file_name, sep="/")
    save(tmp.df, file=file_path)
  }
}


## -- TCA ----
dataTCA = clean_data_TCA(imputed=F)

output_dir = "./results/2015-08-03/data.TCA"
unlink(output_dir, recursive = T, force = FALSE)
dir.create(output_dir, recursive=T)

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.TCA")

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.TCA")

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.TCA")

## -- TCA proteins log ----
createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins.log, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.TCA.log")

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins.log, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.TCA.log")

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins.log, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.TCA.log")

## -- TCA proteins log quant ----
createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins.log.quant, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.TCA.log.quant")

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins.log.quant, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.TCA.log.quant")

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins.log.quant, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.TCA.log.quant")



## -- TCA imputed ----
dataTCA.imputed = clean_data_TCA(imputed=T)

output_dir = "./results/2015-08-03/data.TCA.imputed"
unlink(output_dir, recursive = T, force = FALSE)
dir.create(output_dir, recursive=T)

createDataset(response.matrix=dataTCA.imputed$metabolites, 
              predictors.matrix=dataTCA$proteins, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.TCA.imputed")

createDataset(response.matrix=dataTCA.imputed$metabolites, 
              predictors.matrix=dataTCA$proteins, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.TCA.imputed")

createDataset(response.matrix=dataTCA$metabolites, 
              predictors.matrix=dataTCA$proteins, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.TCA.imputed")


intersect(rownames(dataTCA$metabolites), rownames(dataAA$metabolites))

## -- AA ----
dataAA = clean_data_AA()

output_dir = "./results/2015-08-03/data.AA"
unlink(output_dir, recursive = T, force = FALSE)
dir.create(output_dir, recursive=T)

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.AA")

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.AA")


createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.AA")

## -- AA proteins log ----

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins.log, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.AA.log")

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins.log, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.AA.log")

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins.log, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.AA.log")

## -- AA proteins log quant ----

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins.log.quant, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.AA.log.quant")

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins.log.quant, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.AA.log.quant")

createDataset(response.matrix=dataAA$metabolites, 
              predictors.matrix=dataAA$proteins.log.quant, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.AA.log.quant")





## -- PPP ----

dataPPP_AA = clean_data_PPP_AA(imputed=F)
output_dir = "./results/2015-08-03/data.PPP_AA"
unlink(output_dir, recursive = T, force = FALSE)
dir.create(output_dir, recursive=T)

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA")

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.PPP_AA")

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA")

## -- proteins PPP log ----

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins.log, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA.log")

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins.log, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.PPP_AA.log")

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins.log, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA.log")


## -- proteins PPP log quant ----

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins.log.quant, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA.log.quant")

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins.log.quant, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.PPP_AA.log.quant")

createDataset(response.matrix=dataPPP_AA$metabolites, 
              predictors.matrix=dataPPP_AA$proteins.log.quant, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA.log.quant")






## -- PPP imputed ----
dataPPP_AA.imputed = clean_data_PPP_AA(imputed=T)
output_dir = "./results/2015-08-03/data.PPP_AA.imputed"
unlink(output_dir, recursive = T, force = FALSE)
dir.create(output_dir, recursive=T)

createDataset(response.matrix=dataPPP_AA.imputed$metabolites, 
              predictors.matrix=dataPPP_AA.imputed$proteins, order=1, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA.imputed")


createDataset(response.matrix=dataPPP_AA.imputed$metabolites, 
              predictors.matrix=dataPPP_AA.imputed$proteins, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.PPP_AA.imputed")


createDataset(response.matrix=dataPPP_AA.imputed$metabolites, 
              predictors.matrix=dataPPP_AA.imputed$proteins, order=3, include.metabolites=F, output_dir=output_dir, preffix="data.PPP_AA.imputed")



## -- TCA & AA ---- 
output_dir = "./results/2015-08-03/data.TCA_AA"
common = intersect(rownames(dataAA$proteins), rownames(dataTCA$proteins))
unlink(output_dir, recursive = T, force = FALSE)
dir.create(output_dir, recursive=T)


common_metabolites = cbind(dataTCA$metabolites[match(common, rownames(dataTCA$metabolites)),], dataAA$metabolites[match(common, rownames(dataAA$metabolites)),])
common_proteins = dataTCA$proteins[match(rownames(common_metabolites), rownames(dataTCA$proteins)),]

createDataset(response.matrix=common_metabolites,
              predictors.matrix=common_proteins, order=2, include.metabolites=T, output_dir=output_dir, preffix="data.TCA_AA")







  


