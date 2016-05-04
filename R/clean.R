
rm(list=ls())
source("./R/boot.R")
output_dir = "./R/objects"
dir.create(output_dir)
suffix = "_clean_"



create.peptides = function() {
  load("./R/objects/data.frame.peptides.raw.RData") # load.R::load_proteome()
  load("./R/objects/dates_map.RData")               # load.R::
  load("./R/objects/experiment.map.RData")
  
  peptides.raw <- tbl_df(dataset.peptides.raw)  
  peptides.filtered = peptides.raw
  stopifnot(identical(peptides.filtered$R.FileName, peptides.filtered$R.Label))
  
  
  
  peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  
  peptides.filtered$fragment = factor(with(peptides.filtered,  paste(F.FrgType, F.FrgNum, F.FrgLossType, sep=".")))
  
#   peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, 
#                          EG.Label, EG.StrippedSequence, EG.Qvalue,  EG.Id, fragment, FG.PrecursorMz, FG.TotalPeakArea, F.PeakArea )

  peptides.data = dplyr::select(peptides.filtered, R.Label, R.FileName, sample, replicate, batch,
                         EG.StrippedSequence, FG.Id, EG.Qvalue, FG.TotalPeakArea,
                         fragment,
                         F.InterferenceScore, 
                         F.PossibleInterference, 
                         F.PeakArea )
  
  peptides.data$EG.StrippedSequence = factor(peptides.data$EG.StrippedSequence)
  peptides.data$R.Label = factor(peptides.data$R.Label)
  


  #peptides.data$date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
#   peptides.data$batch.exp = factor(experiment_map$batch[match(peptides.data$R.Label, experiment_map$SampleName)])
#   peptides.data$batch.exp.n = factor(as.numeric(peptides.data$batch.exp))
#   peptides.data$batch_date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
#   
#   peptides.data$acquisition_date = factor(dates_map$Date.of.acquisition[match(peptides.data$R.Label, dates_map$Data_file_name)])
  

#   tmp.tbl = peptides.data %>% group_by(R.Label, batch) %>% distinct(R.Label, batch) %>% dplyr::select(R.Label, batch)
#   file_name = paste("R.Label_batch", suffix, "tsv", sep=".")
#   file_path = paste(output_dir, file_name, sep="/")
#   write.table(x=tmp.tbl, file=file_path, quote=F, sep="\t", row.names=F, col.names=T)

  
  file_name = paste("peptides.data", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(peptides.data,file=file_path)  
  
  
#   tmp.peptides.subset = droplevels(subset(peptides.data, R.Label == peptides.data$R.Label[1]))
#   
#   protein.list = list()
#   u.proteins = unique(tmp.peptides.subset$EG.Label)
#   u.proteins.peptides = unique(select(tmp.peptides.subset, EG.Label, EG.StrippedSequence))
#     
#   for (i in 1:length(u.proteins) ) {
#     
#     proteins = strsplit(x=u.proteins[i], split="[.]")[[1]][1]
#     orfs = strsplit(x=proteins, split="[;]")
#     tmp_split = grep("^[A-Za-z]+$", unlist(strsplit(orfs[[1]], split="[|]")), invert=T, value=T, perl=T)
#     tmp = as.data.frame(matrix(tmp_split, ncol=2, byrow=T))
#     names(tmp) = c("sgd.id", "ORF")
#     tmp$initial.id = u.proteins[i]
#     protein.list[[i]] = tmp
#   }
#   
#   peptide2orfs = do.call(rbind.data.frame, protein.list)
#   peptide2orfs$EG.StrippedSequence = u.proteins.peptides$EG.StrippedSequence[match(peptide2orfs$initial.id, u.proteins.peptides$EG.Label )]
#   
#   
#   #length(unique(grep(".*", perl=T, peptide2orfs$sgd.id, value=T)))
#   file_name = "peptides2orfs.RData"
#   file_path = paste(output_dir, file_name, sep="/")
#   save(peptide2orfs,file=file_path)  
    
}

create.peptides2 = function() {
  
  load("./R/objects/data.frame.peptides.raw_2016_01_16._load_.RData") # load.R::load_proteome()
  load("./R/objects/dates_map.RData")               # load.R::
  load("./R/objects/experiment.map.RData")
  
  peptides.raw <- tbl_df(dataset.peptides.raw)  
  peptides.filtered = peptides.raw
  #stopifnot(identical(peptides.filtered$R.FileName, peptides.filtered$R.Label))
  
  
  
  #peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  #peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  
  peptides.filtered$fragment = factor(with(peptides.filtered,  paste(F.FrgType, F.FrgNum, F.FrgLossType, sep=".")))
  
  #   peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, 
  #                          EG.Label, EG.StrippedSequence, EG.Qvalue,  EG.Id, fragment, FG.PrecursorMz, FG.TotalPeakArea, F.PeakArea )
  
  peptides.data = dplyr::select(peptides.filtered, R.Label, file,
                                EG.StrippedSequence, FG.Id, EG.Qvalue, FG.TotalPeakArea,
                                fragment,
                                F.PeakArea )
  
  peptides.data$EG.StrippedSequence = factor(peptides.data$EG.StrippedSequence)
  peptides.data$R.Label = factor(peptides.data$R.Label)
  

  file_name = paste("peptides.data_2016_01_16", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(peptides.data,file=file_path)  
}

create.exp_annotations = function() {
  #creating metadate
  load("./R/objects/experiment.map._load_.RData")
  load("./R/objects/sample.map._load_.RData")
  load("./R/objects/acqusition_times._load_.RData")
  
  
  experiment_map.f = experiment_map[which(experiment_map$SampleName %in% acqusition_times$sample_name),]
  sample_exp.map = merge(unique(sample_map), subset(experiment_map.f, select=c("SampleName", "batch", "date"), by="SampleName", all=T))

  file_name = paste("sample_exp.map", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(sample_exp.map,file=file_path)  
  
}


createMetadata = function() {
  
  #creating metadate
  load("./R/objects/experiment.map._load_.RData")
  load("./R/objects/sample.map._load_.RData")
  load("./R/objects/acqusition_times._load_.RData")
  load("./R/objects/dates_map._load_.RData")
  load("./R/objects/gene.annotations._load_.RData")
  
  
  
  orf2gene_name = unique(data.frame(ORF = gene.annotations$V4, gene_name = gene.annotations$V6))
  orf2gene_name$ORF = as.character(orf2gene_name$ORF)
  orf2gene_name$gene_name = as.character(orf2gene_name$gene_name)
  
  stopifnot(!duplicated(orf2gene_name$ORF))
  
  
  
  exp_annotation = data.frame(sample_name = acqusition_times$sample_name,
                              filename = acqusition_times$filename,
                              ORF = sample_map$ORF[match(acqusition_times$sample_name, sample_map$SampleName)],
                              gene = sample_map$gene[match(acqusition_times$sample_name, sample_map$SampleName)],
                              type = sample_map$Type[match(acqusition_times$sample_name, sample_map$SampleName)],
                              aquisition_date = acqusition_times$AcquisitionDate)
  
  exp_annotation$ORF  = as.character(exp_annotation$ORF)
  exp_annotation$gene = as.character(exp_annotation$gene)
  

  exp_annotation$ORF[exp_annotation$sample_name == "KL_Try_016_c"] = "WT"
  exp_annotation$ORF[exp_annotation$sample_name == "KL_Try_027_c"] = "WT"
  exp_annotation$ORF[exp_annotation$sample_name == "KL_Try_013_c"] = "YPL031C"
  exp_annotation$ORF[exp_annotation$sample_name == "KL_Try_051_b"] = "YLR354C"
  exp_annotation$ORF[exp_annotation$sample_name == "KL_Try_116_c"] = "YKL168C"
  exp_annotation$ORF[exp_annotation$sample_name == "KL_Try_125_c"] = "YKL161C"
  
  exp_annotation = exp_annotation[grep(x=exp_annotation$sample_name, pattern="wash", ignore.case=T, invert=T),]
  
  exp_annotation$ORF[grep(x=exp_annotation$sample_name, pattern="mix", ignore.case=T)] = "none"
  exp_annotation$type[grep(x=exp_annotation$sample_name, pattern="mix", ignore.case=T)] = "Standard Mix"
  
  exp_annotation$gene = orf2gene_name$gene_name[match(exp_annotation$ORF, orf2gene_name$ORF)]
  exp_annotation$gene[grep(x=exp_annotation$sample_name, pattern="mix", ignore.case=T)] = "none"
  exp_annotation$gene[exp_annotation$ORF =="WT"] = "WT"
  
  exp_annotation$type[exp_annotation$ORF =="WT"] = "Wild Type"
  exp_annotation = exp_annotation[exp_annotation$filename !="KL_Try_028_a.wiff.1.~idx2",]
  exp_annotation = exp_annotation[exp_annotation$filename !="KL_batches_JV_v01.xlsx",]
  
  exp_annotation$batch.exp  = factor(experiment_map$batch[match(exp_annotation$sample_name, experiment_map$SampleName)])
  exp_annotation$batch_date = factor(experiment_map$date[match(exp_annotation$sample_name, experiment_map$SampleName)])
  exp_annotation$batch.exp.n = factor(as.numeric(exp_annotation$batch.exp))
  
  exp_annotation = exp_annotation %>% group_by(ORF) %>% mutate (type = type[!is.na(type)][1])
  exp_annotation = exp_annotation %>% group_by(ORF) %>% mutate (batch.exp = batch.exp[!is.na(batch.exp)][1])
  
  exp_annotation = droplevels(exp_annotation)
  exp_annotation$gene[exp_annotation$gene == ""] = exp_annotation$ORF[exp_annotation$gene == ""]
  exp_annotation$ORF = factor(exp_annotation$ORF)
  exp_annotation$gene = factor(exp_annotation$gene)
  
  exp_annotation[exp_annotation$sample_name == "KL_Try_027_c",]$batch_date = "2014_04_14"
  exp_annotation[exp_annotation$sample_name == "KL_Try_027_c",]$batch.exp.n = 5
  
  
  #stopifnot(!is.na(exp_annotation))
  exp_metadata = tbl_dt(exp_annotation)
    
  file_name = paste("exp_metadata", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(exp_metadata,file=file_path)  
}



# createMetabolites = function() {
#   load("./R/objects/metabolites.raw.RData")
#   load("./R/objects/sample_exp.map.RData")
#   load("./R/objects/dates_map._load_.RData")
#     
#     
#   pattern.p = "(\\d+)_KL_([A-Za-z0-9]+)_?(\\d?)"
#   matches.tmp = stringr::str_match_all(string=metabolites.raw[,1], pattern=pattern.p)
#   
#   pheno.met = data.frame(matrix(unlist(matches.tmp), byrow=T, ncol=length(matches.tmp[[1]])))
#   names(pheno.met) = c("SampleName", "batch", "sample.id", "replicate")
#   
#   pheno.met$sample.id = paste("KL",as.character(pheno.met$sample.id), sep="")
#   pheno.met$sample.id[grep(pattern="WT", ignore.case=T, x=pheno.met$sample.id)] = "WT"
#   pheno.met$sample.id = factor(pheno.met$sample.id)
#   
#   
#   pheno.met$ORF = as.character(dates_map$ORF[match(pheno.met$sample.id, dates_map$Nr)])
#   pheno.met$ORF[grep(pattern="WT", ignore.case=T, x=pheno.met$sample.id)] = "WT"
#   pheno.met$ORF = factor(pheno.met$ORF)
#   phenotypes.f = droplevels(pheno.met[!is.na(pheno.met$ORF),])
#   metabolite.matrix = as.matrix(t(metabolites.raw[,-1]))
#   colnames(metabolite.matrix) = metabolites.raw$Sample
#   
#   
#   metabolite.matrix.long = melt(metabolite.matrix, id.vars=row.names)
#   names(metabolite.matrix.long) = c("id", "variable", "value")
#   
#   metabolite.matrix.long = tbl_df(metabolite.matrix.long %>% extract(variable, 
#                                             into=c("batch", "sample"),  
#                                             regex="(\\d+)_(\\w+)"))
#     
#   
#   
#   
#   #metabolites.matrix.f1 = metabolite.matrix[,colnames(metabolite.matrix) %in% phenotypes.f$SampleName]
#   metabolites.matrix.f = metabolite.matrix[,match(phenotypes.f$SampleName, colnames(metabolite.matrix))]
#   
#   metabolites.folds = rowFolds(data=metabolites.matrix.f, groups=phenotypes.f$ORF, reference="WT")
#   metabolites.folds[is.na(metabolites.folds)] = NA
#   
#   metabolite.folds.long = reshape2::melt(as.matrix(metabolites.folds), id.vars=row.names)
#   names(metabolite.folds.long) = c("id", "variable", "value")
#   
# 
#   metabolites.data = metabolite.matrix.long
#   file_name = paste("metabolites.data", suffix, "RData", sep=".")
#   file_path = paste(output_dir, file_name, sep="/")
#   save(metabolites.data,file=file_path)  
#   
#   metabolites.folds = metabolite.folds.long
#   file_name = paste("metabolites.folds", suffix, "RData", sep=".")
#   file_path = paste(output_dir, file_name, sep="/")
#   save(metabolites.folds,file=file_path)  
# 
# }


createMetabolites2 = function() {
  load("./R/objects/dataset.metabolites.raw._load_.RData")
  load("./R/objects/dates_map._load_.RData")
    
  dataset.metabolites.raw$sample_id = factor(paste(dataset.metabolites.raw$Experiment, dataset.metabolites.raw$Sample, sep="_"))
    
  pattern.p = "(\\d+)_KL_([A-Za-z0-9]+)_?(\\d?)"
  matches.tmp = stringr::str_match_all(string=dataset.metabolites.raw$sample_id, pattern=pattern.p)
  
  pheno.met = data.frame(matrix(unlist(matches.tmp), byrow=T, ncol=length(matches.tmp[[1]])))
  
  names(pheno.met) = c("sample_id", "batch", "sample", "replicate")
  
  pheno.met$sample = paste("KL",as.character(pheno.met$sample), sep="")
  pheno.met$sample[grep(pattern="WT", ignore.case=T, x=pheno.met$sample)] = "WT"
  pheno.met$sample = factor(pheno.met$sample)
  
  
  pheno.met$ORF = as.character(dates_map$ORF[match(pheno.met$sample, dates_map$Nr)])
  pheno.met$ORF[grep(pattern="WT", ignore.case=T, x=pheno.met$sample)] = "WT"
  pheno.met$ORF = factor(pheno.met$ORF)
  
  metabolite_metadata = pheno.met
  
  file_name = paste("metabolite_metadata", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolite_metadata,file=file_path)  
  
  
  metabolites.data = melt(dplyr::select(dataset.metabolites.raw, -file, -Sample, -Experiment,-batch), id.vars=c("sample_id"))
  
  file_name = paste("metabolites.data", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolites.data,file=file_path)  
  
}


createMetabolitesTCA = function() {
  
  load("./R/objects/metabolitesTCA.raw._load_.RData")
  load("./R/objects/dates_map._load_.RData")
  metabolitesTCA.raw$Name = trimWhiteSpace(metabolitesTCA.raw$Name)
  
  
  pattern.p = "Sample_(.*?)_(\\w+)(.*)"
  
  matches.tmp = stringr::str_match_all(string=metabolitesTCA.raw$Name, pattern=pattern.p)
  pheno = data.frame(matrix(unlist(matches.tmp), byrow=T, ncol=length(matches.tmp[[1]])))
    
  names(pheno) =  c("sample_name", "sample", "replicate", "measure_date")
  pheno$measure_date = trimWhiteSpace(pheno$measure_date)
  pheno$measure_date = gsub(pattern="[\\(\\)]+", replacement="", perl=T, x=pheno$measure_date)
  pheno$measure_batch = as.numeric(factor(pheno$measure_date))
  
  pheno = pheno %>% group_by(measure_batch, sample) %>% mutate(number = 1:length(sample))
  
  
  pheno$sample_id = with(pheno, paste(measure_batch, sample, replicate, number, sep="_" ))
  pheno$ORF = as.character(dates_map$ORF[match(pheno$sample, dates_map$Nr)])
  pheno$ORF[grep(pattern="WT",ignore.case=T, x=pheno$sample_name)] = "WT"
  
  
  metabolitesTCA.data = cbind(pheno$sample_id, metabolitesTCA.raw[,-c(1,2)]) #pheno is sorted in the same way as metabolitesTCA.raw
  metabolitesTCA.data[metabolitesTCA.data == 0] = NA
  names(metabolitesTCA.data)[1] = "sample_id"
  
  metabolitesTCA_metadata = pheno
      
  file_name = paste("metabolitesTCA_metadata", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolitesTCA_metadata,file=file_path)  
    
  metabolitesTCA.data = melt(metabolitesTCA.data, id.vars=c("sample_id"))
  
  #removing weird values
  metabolitesTCA.data$value[which(with(metabolitesTCA.data, variable == "S7P" & (sample_id == "3_WT_A_1" | sample_id == "3_WT_B_2" | sample_id == "3_WT_C_3" |
                                                 sample_id == "5_WT_A_1" | sample_id == "5_WT_B_2" | sample_id == "5_WT_C_3" )))] = NA
  
  file_name = paste("metabolitesTCA.data", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolitesTCA.data,file=file_path)
  
}

createAA = function() {
  load("./R/objects//aa.raw._load_.RData")
  #load("./R/objects/exp_metadata._clean_.RData")
  load("./R/objects/metabolite_metadata._clean_.RData")
  load("./R/objects/dates_map._load_.RData")
  
#   file_name = paste("kinases", suffix, "txt", sep=".")
#   file_path = paste(output_dir, file_name, sep="/") 
#   write.table(x=unique(exp_metadata$ORF[exp_metadata$type == "Kinase"]),col.names=F, row.names=F, file=file_path, quote=F)
#   
  aa.raw$Strain = sub(pattern="^0(.*?)", replacement="\\1", x=as.character(aa.raw$Strain))
  aa.raw$Strain = factor(sub(pattern="^(\\d+)", replacement="KL\\1", x=as.character(aa.raw$Strain), perl=T))
  
       
  aa.raw$sample_id = factor(with(aa.raw, paste(Batch, Strain, Replicate,Date, sep="_")))
  aa_metadata = aa.raw %>% dplyr::select(sample_id, Batch, Strain, Replicate, Date)
  
  aa_metadata$ORF = as.character(dates_map$ORF[match(aa_metadata$Strain, dates_map$Nr)])
  aa_metadata$ORF[is.na(aa_metadata$ORF)] = as.character(aa_metadata$Strain[which(is.na(aa_metadata$ORF))])
  aa_metadata$ORF[grep(pattern="Wt", aa_metadata$ORF)] = "WT"
  
  file_name = paste("aa_metadata", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(aa_metadata,file=file_path)  
  
  aa.data = melt(dplyr::select(aa.raw, -Batch, -Strain, -Replicate, -Date), id.vars=c("sample_id"))
  file_name = paste("aa.data", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(aa.data,file=file_path)  

}

createAA_michael = function() {
  load("./R/objects//aa_michael.raw._load_.RData")
  
  aa_michael.metadata = aa_michael.raw %>% dplyr::select(Name, gene, ORF, Batch, Time)
  names(aa_michael.metadata) = c("sample_id", "gene", "ORF", "batch", "date")
  aa_michael.data = melt(aa_michael.raw %>% dplyr::select(-gene, -ORF, -Batch, -Time), id.vars="Name")
  names(aa_michael.data) = c("sample_id", "variable", "value")
  
  aa_michael.metadata$ORF = as.character(aa_michael.metadata$ORF)
  aa_michael.metadata$ORF[grep(pattern="QC", x=aa_michael.metadata$sample_id, ignore.case=T)] = "QC"
  aa_michael.metadata$ORF = factor(aa_michael.metadata$ORF)
  
  aa_michael.metadata$batch = factor(aa_michael.metadata$batch)
  aa_michael.metadata$date = as.Date(aa_michael.metadata$date)

  file_name = paste("aa_michael.metadata", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(aa_michael.metadata,file=file_path)  
  
  file_name = paste("aa_michael.data", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(aa_michael.data,file=file_path)  
  
} 


createPicotti_Data_1 = function () {
  picotti.raw = read.xlsx2("./data/2015-04-16/nature11835-s2.xlsx", 
                           sheetIndex=5, startRow=3, header=T) # discovery data
  picotti.raw[picotti.raw == ""] = NA
  ind_cols = 5:ncol(picotti.raw)
  picotti.raw[,ind_cols] = apply(picotti.raw[,ind_cols], 2, FUN=function(x) {as.numeric(as.character(x))})
  
  peptides.LIT = picotti.raw[,c(-3,-4)]
  names(peptides.LIT)[c(1,2)] = c("ORF", "EG.StrippedSequence")
    
  
  file_name = paste("peptides.LIT", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(peptides.LIT,file=file_path)
    
}

createKinaseClasses = function() {
  load("./R/objects/gene.annotations._load_.RData")
    
  
  kinase_classes.raw <- read.delim("./data/2015-04-28/kinase.classes.txt", header=F)
  list.splited = strsplit(as.character(kinase_classes.raw$V2), split=" ",fixed=T)  
  tmp.list = list()
  for (i in 1:length(list.splited)) {
    tmp = data.frame(kinase = list.splited[[i]])
    tmp$class = kinase_classes.raw$V1[i]
    tmp.list[[i]] = tmp
  }
  kinase_classes = do.call(rbind.data.frame, tmp.list)
  
  kinase_classes$ORF = as.character(gene.annotations$V4[match(kinase_classes$kinase, gene.annotations$V6)])
  kinase_classes$ORF[kinase_classes$kinase == "ABC1" ] = "YGL119W"
  
  
  
  kinase_classes$ORF[is.na(kinase_classes$ORF)] = as.character(kinase_classes$kinase[is.na(kinase_classes$ORF)])
    
  file_name = paste("kinase_classes", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(kinase_classes,file=file_path)
}
  
create.BRENDA = function () {
  ec.brenda.raw <- read.delim("./results/2015-10-01/ec.brenda.tsv", header=F)
  ec.brenda.raw = ec.brenda.raw[, -length(ec.brenda.raw)]
  names_cols = seq(1,16, 2)
  data_cols  = seq(2,16, 2)
  col_names = ec.brenda.raw[1,names_cols]
  data = data.frame(ec.brenda.raw[,data_cols])
  colnames(data) = as.character(t(col_names))
  ec.brenda = data
  file_name = paste("ec.brenda", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(ec.brenda,file=file_path)
}

create.kinase_transcrtiptome = function() {
  
  data.raw <- read.delim("./data/2015-09-27/DataS1.txt", header=F)
  
  load("./R/objects/gene.annotations._load_.RData")
  
  orf2name = droplevels(unique(gene.annotations[,c("V4", "V6")]))
  orf2name$V4 = as.character(orf2name$V4)
  orf2name$V6 = as.character(orf2name$V6)
  orf2name$V6[orf2name$V6 == ""] = orf2name$V4[orf2name$V6 == ""]
  names(orf2name) = c("ORF", "gene_name")
  
  
  colNames = gsub(pattern="-del", replacement="", t(data.raw[1,]))
  colNames = gsub(pattern=" |vs.|wt", replacement="", colNames)
  colNames = gsub(pattern="[+-]", replacement="_", colNames)
  
  colnames(data.raw) = toupper(colNames)
  
  select_pvalues = grep(x=as.character(t(data.raw[2,])), pattern="p_value", perl=T)
  select_folds = which(as.character(t(data.raw[2,])) == "M")
  
  folds.matrix = data.raw[-c(1,2),c(1,select_folds)] 
  pvalue.matrix = data.raw[-c(1,2),c(1,select_pvalues)] 
  
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  
  data = data.raw
  folds.matrix[, -1] <- sapply(folds.matrix[, -1], as.numeric.factor)
  pvalue.matrix[, -1] <- sapply(pvalue.matrix[, -1], as.numeric.factor)
  
  folds.long  = melt(folds.matrix, id.vars="SYSTEMATICNAME")
  pvalue.long = melt(pvalue.matrix, id.vars="SYSTEMATICNAME")
  
  names(folds.long)  = c("ORF", "KO", "logFC")
  names(pvalue.long) = c("ORF", "KO", "p.value")
  
  transcriptome.FC = merge(folds.long, pvalue.long, by = c("ORF", "KO"))
  
  transcriptome.FC$KO.ORF = orf2name$ORF[match(transcriptome.FC$KO, orf2name$gene_name)]
  transcriptome.FC$KO.ORF[which(is.na(transcriptome.FC$KO.ORF))] = orf2name$ORF[match(transcriptome.FC$KO[which(is.na(transcriptome.FC$KO.ORF))], orf2name$ORF)]
  transcriptome.FC = droplevels(transcriptome.FC[grep(pattern="_", x=transcriptome.FC$KO, invert=T),])
  transcriptome.FC$KO.ORF[transcriptome.FC$KO == "ABC1"] = "YGL119W"
  transcriptome.FC$KO.ORF[transcriptome.FC$KO == "CDK8"] = "YPL042C"
  transcriptome.FC$KO.ORF[transcriptome.FC$KO == "SHA3"] = "YPL026C"
  transcriptome.FC = droplevels(transcriptome.FC)
  transcriptome.FC$KO.ORF = factor(transcriptome.FC$KO.ORF)
  transcriptome.FC$p.value_BH= p.adjust(transcriptome.FC$p.value, method="BH")
  
  file_name = paste("transcriptome.FC", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(transcriptome.FC,file=file_path)
}

create.sentinel_list = function() {
  sentinels.table = read.xlsx2("./data/2015-10-20//nmeth.3101-S2.xlsx", 
                           sheetIndex=1, startRow=2, header=T) # discovery data
  
  file_name = paste("sentinels.table", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(sentinels.table,file=file_path)
}


create.sentinelSRM_list = function() {
  sentinelsSRM = read.xlsx2("./data/2015-10-20/nmeth.3101-S4.xlsx", 
                            sheetIndex=1, startRow=12, header=T) # discovery data
  
  sentinelsSRM$ID.ORF = sub(pattern="^([A-Za-z0-9]+)..*", replacement="\\1", x=sentinelsSRM$ID, perl=T)
    
  file_name = paste("sentinelsSRM", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(sentinelsSRM,file=file_path)
}

create_irt_dataset <- function() {
  
  
  db_path = "./results/2016-01-17/results/OpenSwathWorkflow_with_dscore_unagrigated_filtered.db"
  my_db <- src_sqlite(db_path, create = F)
  dataset_sqlite <- tbl(my_db, "my_table")
  
  main_table = "my_table"
  
  my_query <- paste0("SELECT * FROM my_table
                        WHERE decoy = 0 AND 
                        peak_group_rank = 1 AND 
                        Peak_Area > 0 AND
                        ProteinName LIKE '%Biognosys%'")
  
  dataset_irt_openswath <- tbl(my_db,  sql(my_query)) %>% collect()
  
  file_name = paste("dataset_irt_openswath", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(dataset_irt_openswath,file=file_path)  
  
}


create.Sentinels_dataset <- function() {
  load("./R/objects/sentinels.raw._load_.RData") # load.R::load_sentinels_dataset()
  
  sentinels.raw <- tbl_df(sentinels.raw)  
  #stopifnot(identical(peptides.filtered$R.FileName, peptides.filtered$R.Label))
  
  #peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  #peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  
  sentinels.raw$fragment = with(sentinels.raw,  paste(F.FrgType, F.FrgNum, F.FrgLossType, sep="."))
  
  sentinels.data = dplyr::select(sentinels.raw, R.Label, EG.StrippedSequence, FG.Id, EG.Qvalue, FG.TotalPeakArea, fragment,  F.PeakArea )
  
  file_name = paste("sentinels.data", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(sentinels.data,file=file_path)  
}



main = function() {
  create.peptides()
  create.exp_annotations()
  createMetabolites()
  #createPicotti_Data_1()
  createMetadata()
  createAA() 
  createAA_michael()
  createMetabolites2()
  
  createMetabolitesTCA()
  create.BRENDA()
  create.kinase_transcrtiptome()
  create.sentinel_list()
  create.sentinelSRM_list()
  create.Sentinels_dataset()
}


#main()

