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
  str(peptides.filtered)
  
  peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  #peptides.filtered$fragment = with(peptides.filtered,  paste(F.FrgType, F.FrgNum,  sep="."))
  
#   peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, 
#                          EG.Label, EG.StrippedSequence, EG.Qvalue,  EG.Id, fragment, FG.PrecursorMz, FG.TotalPeakArea, F.PeakArea )

  peptides.data = dplyr::select(peptides.filtered, R.Label, R.FileName, sample, replicate, batch,
                         EG.StrippedSequence, EG.Qvalue, FG.TotalPeakArea,
                         F.InterferenceScore, 
                         F.PossibleInterference, 
                         F.PeakArea )
  
    
  #peptides.data$date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  peptides.data$batch.exp = factor(experiment_map$batch[match(peptides.data$R.Label, experiment_map$SampleName)])
  peptides.data$batch.exp.n = factor(as.numeric(peptides.data$batch.exp))
  peptides.data$batch_date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  
  peptides.data$acquisition_date = factor(dates_map$Date.of.acquisition[match(peptides.data$R.Label, dates_map$Data_file_name)])
  

  
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
  load("./R/objects/sample_exp.map._load_.RData")
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
  
  exp_annotation$batch.exp  = factor(sample_exp.map$batch[match(exp_annotation$sample_name, sample_exp.map$SampleName)])
  exp_annotation$batch_date = factor(sample_exp.map$date[match(exp_annotation$sample_name, sample_exp.map$SampleName)])
  exp_annotation$batch_date = factor(as.numeric(exp_annotation$batch.exp))
  
  exp_annotation = exp_annotation %>% group_by(ORF) %>% mutate (type = type[!is.na(type)][1])
  exp_annotation = exp_annotation %>% group_by(ORF) %>% mutate (batch.exp = batch.exp[!is.na(batch.exp)][1])
  
  exp_annotation = droplevels(exp_annotation)
  
  
  #stopifnot(!is.na(exp_annotation))
  exp_metadata = tbl_dt(exp_annotation)
    
  file_name = paste("exp_metadata", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(exp_metadata,file=file_path)  
}



createMetabolites = function() {
  load("./R/objects/metabolites.raw.RData")
  load("./R/objects/sample_exp.map.RData")
  load("./R/objects/dates_map._load_.RData")
    
  
  
  pattern.p = "(\\d+)_KL_([A-Za-z0-9]+)_?(\\d?)"
  matches.tmp = stringr::str_match_all(string=metabolites.raw[,1], pattern=pattern.p)
  
  pheno.met = data.frame(matrix(unlist(matches.tmp), byrow=T, ncol=length(matches.tmp[[1]])))
  names(pheno.met) = c("SampleName", "batch", "sample.id", "replicate")
  
  pheno.met$sample.id = paste("KL",as.character(pheno.met$sample.id), sep="")
  pheno.met$sample.id[grep(pattern="WT", ignore.case=T, x=pheno.met$sample.id)] = "WT"
  pheno.met$sample.id = factor(pheno.met$sample.id)
  
  
  pheno.met$ORF = as.character(dates_map$ORF[match(pheno.met$sample.id, dates_map$Nr)])
  pheno.met$ORF[grep(pattern="WT", ignore.case=T, x=pheno.met$sample.id)] = "WT"
  pheno.met$ORF = factor(pheno.met$ORF)
  phenotypes.f = droplevels(pheno.met[!is.na(pheno.met$ORF),])
  metabolite.matrix = as.matrix(t(metabolites.raw[,-1]))
  colnames(metabolite.matrix) = metabolites.raw$Sample
  
  
  metabolite.matrix.long = melt(metabolite.matrix, id.vars=row.names)
  names(metabolite.matrix.long) = c("id", "variable", "value")
  
  metabolite.matrix.long = tbl_df(metabolite.matrix.long %>% extract(variable, 
                                            into=c("batch", "sample"),  
                                            regex="(\\d+)_(\\w+)"))
    
  
  
  
  #metabolites.matrix.f1 = metabolite.matrix[,colnames(metabolite.matrix) %in% phenotypes.f$SampleName]
  metabolites.matrix.f = metabolite.matrix[,match(phenotypes.f$SampleName, colnames(metabolite.matrix))]
  
  metabolites.folds = rowFolds(data=metabolites.matrix.f, groups=phenotypes.f$ORF, reference="WT")
  metabolites.folds[is.na(metabolites.folds)] = NA
  
  metabolite.folds.long = reshape2::melt(as.matrix(metabolites.folds), id.vars=row.names)
  names(metabolite.folds.long) = c("id", "variable", "value")
  
  

  metabolites.data = metabolite.matrix.long
  file_name = paste("metabolites.data", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolites.data,file=file_path)  
  
  metabolites.folds = metabolite.folds.long
  file_name = paste("metabolites.folds", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolites.folds,file=file_path)  

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
  



main = function() {
  #create.peptides()
  create.exp_annotations()
  createMetabolites()
  createPicotti_Data_1()
  createMetadata()
}


#main()

