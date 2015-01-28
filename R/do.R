rm(list=ls())
source("./R/functions.R")
figures_dir = "./figures"
output_dir = "./R/objects"
dir.create(figures_dir)
dir.create(output_dir)


explore = function() {
  
  plots.list = list()
  fun_name = "explore"
  
  load("R/objects/data.frame.peptides.raw.RData")
  load("R/objects/experiment.map.RData")
  
  peptides.raw <- tbl_df(dataset.peptides.raw)  
  peptides.filtered = peptides.raw
  
  peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_(\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  
  peptides.filtered$fragment = with(peptides.filtered,  paste(F.FrgType, F.FrgNum,  sep="."))
  
  peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, EG.Label, EG.StrippedSequence, EG.Id, fragment, FG.PrecursorMz, FG.TotalPeakArea, F.PeakArea )
  peptides.data$date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  peptides.data.mix = droplevels(peptides.data[grep(pattern="mix", ignore.case=T, peptides.data$R.Label),])
  
  
  ##batch/date signal statistics
  grouped.data <- group_by(peptides.data, batch, date, R.Label)
  grouped_stats <- dplyr::summarise(grouped.data,
                                    sum_FG.TotalPeakArea = sum(FG.TotalPeakArea, na.rm=T),
                                    mean_FG.TotalPeakArea = mean(FG.TotalPeakArea, na.rm=T))
  
  grouped_stats <- dplyr::summarise(grouped.data,
                                    sum_FG.TotalPeakArea = sum(FG.TotalPeakArea, na.rm=T),
                                    mean_FG.TotalPeakArea = mean(FG.TotalPeakArea, na.rm=T))
  
  p = ggplot(grouped_stats, aes(x=batch, y = mean_FG.TotalPeakArea, colour=date)) + 
    geom_boxplot() + 
    theme(aspect.ratio = 1) +
    ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
  
  plots.list = lappend(plots.list, p)
  
  grouped_stats$sample_number = 1:nrow(grouped_stats)
  p = ggplot(grouped_stats, aes(x=sample_number, y = mean_FG.TotalPeakArea, col=date)) + 
    geom_point() + 
    theme(aspect.ratio = 1) +
    ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
  
  plots.list = lappend(plots.list, p)
  
  p = ggplot(grouped_stats, aes(x=sample_number, y = mean_FG.TotalPeakArea, col=batch)) + 
    geom_boxplot() +
    geom_point() + 
    theme(aspect.ratio = 1) +
    ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
  
  plots.list = lappend(plots.list, p)
  
  #peptide distributions per samples
  
  grouped.data <- group_by(peptides.data, batch, date, R.Label, EG.StrippedSequence)
  grouped_stats <- dplyr::summarise(grouped.data,
                                    count = n(),
                                    sum_FG.PeakArea = sum(F.PeakArea, na.rm=T))
  grouped_stats$sample_number = as.numeric(grouped_stats$R.Label)
  #peptide.matrix = dcast(data=grouped_stats, formula=EG.StrippedSequence~batch+R.Label, value.var="sum_FG.PeakArea")
#   peptide.matrix$R.Label[(table(peptide.matrix$R.Label) > 1)]
#   grouped_stats$sample_number = as.numeric(grouped_stats$R.Label)
#   
  length(unique(grouped_stats$sample_number))
  p = ggplot(grouped_stats, aes(x=log(sum_FG.PeakArea,base=10), colour=sample_number, group=sample_number)) + 
    geom_density(fill=NA) +
    theme(aspect.ratio = 1) +
    ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
  
  plots.list = lappend(plots.list, p)
  
  p = ggplot(grouped_stats, aes(x=log(sum_FG.PeakArea,base=10), colour=date, group=sample_number)) + 
    geom_density(fill=NA) +
    theme(aspect.ratio = 1) +
    ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
  plots.list = lappend(plots.list, p)
  
  ####################################
  #plotting standart mix distributions
  ####################################
  grouped.data <- group_by(peptides.data.mix, batch, date, R.Label, EG.StrippedSequence)
  grouped_stats <- dplyr::summarise(grouped.data,
                                    count = n(),
                                    sum_FG.PeakArea = sum(F.PeakArea, na.rm=T))
  grouped_stats$sample_number = as.numeric(grouped_stats$R.Label)
  
  p = ggplot(grouped_stats, aes(y=log(sum_FG.PeakArea,base=10), colour=date, x=sample_number)) + 
    geom_boxplot(fill=NA) +
    theme(aspect.ratio = 1) +
    ggtitle(paste("Grouped by", paste(attr(grouped.data,which="vars"), collapse=".")))
  
  plots.list = lappend(plots.list, p)
  
#   tmp = dcast(grouped_stats, formula=batch+R.Label~EG.StrippedSequence, value.var="sum_FG.PeakArea")
#   tmp = tmp[!duplicated(tmp$R.Label),]
#   
#   tmp.matrix = as.matrix(tmp[,-c(1,2)])
#   rownames(tmp.matrix) = with(tmp, paste(batch,R.Label,sep="."))  
#   tmp.matrix = t(tmp.matrix)
#   norm.matrix = tmp.matrix/estimateSizeFactorsForMatrix(tmp.matrix)
#   tmp.matrix = norm.matrix
#   cor.tmp.matrix = cor(tmp.matrix)
  
#   plot(cor.tmp.matrix[upper.tri(cor.tmp.matrix, diag=F)])
  

  file_name = paste(fun_name, "report.pdf", sep=".")
  file_path = paste(figures_dir, file_name, sep="/")
  save_plots(plots.list, filename=file_path) 
  
    
}

#compare batches using standarts 
compare_mix_normalizations = function() {
  
  load("R/objects/sample.map.RData")
  
  standards = unique(as.character(with(sample_map, droplevels(SampleName[grep(pattern="mix", ignore.case=T, x=SampleName)]))))
  
  load("R/objects/data.frame.peptides.raw.RData")
  load("R/objects/experiment.map.RData")
  load("R/objects/sample.map.RData")
  
  peptides.raw <- tbl_df(dataset.peptides.raw)  
  peptides.filtered = peptides.raw
  
  peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_(\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  
  peptides.filtered$fragment = with(peptides.filtered,  paste(F.FrgType, F.FrgNum,  sep="."))
  
  peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, EG.Label, EG.StrippedSequence, EG.Id, fragment, FG.PrecursorMz, FG.TotalPeakArea, F.PeakArea )
  peptides.data$date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  
  sample_map
  View(experiment_map)
  peptides.data
  peptides.data.mix = droplevels(peptides.data[grep(pattern="mix", ignore.case=T, peptides.data$R.Label),])
  
  
  ##batch/date signal statistics, finding confounding
  
  grouped.data <- group_by(peptides.data, R.Label, EG.Label, EG.StrippedSequence)
  grouped_stats <- dplyr::summarise(grouped.data,
                                    sum_FG.TotalPeakArea = sum(F.PeakArea, na.rm=T))
  
  
  #adjusting for batches
  proteins.df = dcast(grouped_stats, EG.Label~R.Label, value.var="sum_FG.TotalPeakArea")
  

  View(sample_map)
  protein.matrix = as.matrix(proteins.df[,-1])
  rownames(protein.matrix) = proteins.df$EG.Label
  
  protein.matrix = protein.matrix[,which(!is.na(match(colnames(protein.matrix), sample_map$SampleName)))]
  protein.matrix = protein.matrix[, which(experiment_map[match(colnames(protein.matrix), experiment_map$SampleName), "batch"] != "")]
    
  batches = as.numeric(factor(experiment_map[match(colnames(protein.matrix), experiment_map$SampleName), "batch"]))
  mutants = as.character(sample_map[match(colnames(protein.matrix), sample_map$SampleName), "ORF" ])
  mutants[mutants == "none"] = "Mix"
  
  
  pheno = data.frame(row.names=colnames(protein.matrix), 
                     batch  = factor(batches), 
                     mutant = factor(mutants))
  
  set.seed(1)
  ind = sample(1:100)
  
  ind
  par(mfrow=c(1,2))
  plot(pca$x[,1], pca$x[,2], col=batches)
  plot(pca$sdev[1:10]^2/sum(pca$sdev^2))
  
  pheno
  mod = model.matrix(~droplevels(pheno$mutant[ind]))
  mod
  
  protein.matrix.clean = ComBat(dat=normalize.quantiles(protein.matrix[,ind]), batch=droplevels(pheno$batch[ind]), mod=mod,)
  
}



main = function() {
  explore()
  
  
}







