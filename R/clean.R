source("./R/boot.R")
output_dir = "./R/objects"
dir.create(output_dir)


create.peptides = function() {
  load("./R/objects/data.frame.peptides.raw.RData") # load.R::load_proteome()
  load("./R/objects/dates_map.RData")               # load.R::
  load("./R/objects/experiment.map.RData")
  peptides.raw <- tbl_df(dataset.peptides.raw)  
  peptides.filtered = peptides.raw
  
  
  peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  peptides.filtered$fragment = with(peptides.filtered,  paste(F.FrgType, F.FrgNum,  sep="."))
  
#   peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, 
#                          EG.Label, EG.StrippedSequence, EG.Qvalue,  EG.Id, fragment, FG.PrecursorMz, FG.TotalPeakArea, F.PeakArea )
  peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, R.FileName,
                         EG.StrippedSequence, EG.Qvalue,  fragment, FG.TotalPeakArea,
                         F.InterferenceScore, 
                         F.PossibleInterference, 
                         F.PeakArea )
  
  #peptides.data$date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  peptides.data$batch.exp = factor(experiment_map$batch[match(peptides.data$R.Label, experiment_map$SampleName)])
  peptides.data$batch_date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  peptides.data$sampling_date = factor(dates_map$Sampling.date[match(peptides.data$sample, dates_map$Data_file_name)])
  peptides.data$processing_date = factor(dates_map$Processing.date[match(peptides.data$sample, dates_map$Data_file_name)])
  peptides.data$acquisition_date = factor(dates_map$Date.of.acquisition[match(peptides.data$sample, dates_map$Data_file_name)])
  peptides.data$batch.exp.n = factor(as.numeric(peptides.data$batch.exp))
  
  file_name = "peptides.data.RData"
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
  load("./R/objects/experiment.map.RData")
  load("./R/objects/sample.map.RData")
  
  experiment_map.f = experiment_map[which(experiment_map$SampleName %in% sample_map$SampleName),]
  sample_exp.map = merge(unique(sample_map), subset(experiment_map.f, select=c("SampleName", "batch", "date"), by="SampleName", all=T))
  
  
  sample_exp.map
  file_name = "sample_exp.map.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(sample_exp.map,file=file_path)  
  
}


createMetabolites = function() {
  load("./R/objects/metabolites.raw.RData")
  load("./R/objects/sample_exp.map.RData")
    
#   pattern.p = "(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$"
#   matches = stringr::str_match_all(pattern=pattern.p, colnames(proteins.matrix))
#   
#   stopifnot(sum(lapply(matches,length)!=0) == ncol(proteins.matrix))
#   pheno = data.frame(matrix(unlist(matches), ncol=length(matches[[1]]), byrow=T))
#   colnames(pheno) = c("name", "R.Label", "batch_date", "batch.exp.n", "batch" )
#   rownames(pheno) = colnames(proteins.matrix)
#   pheno$ORF = droplevels(sample_exp.map$ORF[match(pheno$R.Label, sample_exp.map$SampleName)])
#   pheno$ORF[pheno$R.Label == "KL_Try_027_c"] = "WT"
#   pheno$batch.exp.n[pheno$R.Label == "KL_Try_027_c"] = 5
#   
#   tmp = sub(pattern="(KL_Try_\\d+).*", replacement="\\1", x=sample_exp.map$SampleName)
  
  pattern.p = "(\\d+)_KL_([A-Za-z0-9]+)_?(\\d?)"
  matches.tmp = stringr::str_match_all(string=metabolites.raw[,1], pattern=pattern.p)
  
  pheno.met = data.frame(matrix(unlist(matches.tmp), byrow=T, ncol=length(matches.tmp[[1]])))
  names(pheno.met) = c("SampleName", "batch", "sample.id", "replicate")
    
  pheno.met$sample.id = as.character(pheno.met$sample.id)
  pheno.met$sample.id[grep(pattern="WT", ignore.case=T, x=pheno.met$sample.id)] = "WT"
  pheno.met$sample.id = factor(pheno.met$sample.id)
  
  
  pheno.met$ORF = sample_exp.map$ORF[match(pheno.met$sample.id, sample_exp.map$ProtSampleNo)]
  pheno.met$sample.id[grep(pattern="WT", ignore.case=T, x=pheno.met$sample.id)] = "WT"
  pheno.met$ORF[grep(pattern="WT", ignore.case=T, x=pheno.met$sample.id)] = "WT"
  
  phenotypes.f = droplevels(pheno.met[!is.na(pheno.met$ORF),])
  metabolite.matrix = as.matrix(t(metabolites.raw[,-1]))
  colnames(metabolite.matrix) = metabolites.raw$Sample
  
  
  metabolite.matrix.long = melt(metabolite.matrix, id.vars=row.names)
  names(metabolite.matrix.long) = c("id", "variable", "value")

  metabolite.matrix.long = tbl_df(metabolite.matrix.long %>% extract(variable, 
                                            into=c("batch", "sample"),  
                                            regex="(\\d+)_(\\w+)"))
    
  
  metabolites.matrix.f = metabolite.matrix[,colnames(metabolite.matrix) %in% phenotypes.f$SampleName]
  
  metabolites.folds = rowFolds(data=metabolites.matrix.f, groups=phenotypes.f$ORF, reference="WT")
  metabolites.folds[is.na(metabolites.folds)] = NA
  
  metabolite.folds.long = reshape2::melt(as.matrix(metabolites.folds), id.vars=row.names)
  names(metabolite.folds.long) = c("id", "variable", "value")
  
  

  metabolites.data = metabolite.matrix.long
  file_name = "metabolites.data.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolites.data,file=file_path)  
  
  metabolites.folds = metabolite.folds.long
  file_name = "metabolites.folds.RData"
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
    
  
  file_name = "peptides.LIT.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(peptides.LIT,file=file_path)
  
  
}



main = function() {
  create.peptides()
  create.exp_annotations()
  createMetabolites()
  createPicotti_Data_1()
}


#main()

