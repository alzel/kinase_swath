source("./R/boot.R")
output_dir = "./R/objects"
dir.create(output_dir)


prepare_peptides = function() {
  load("./R/objects/data.frame.peptides.raw.RData") # load.R::load_proteome()
  load("./R/objects/dates_map.RData")               # load.R::
  peptides.raw <- tbl_df(dataset.peptides.raw)  
  peptides.filtered = peptides.raw
  
  peptides.filtered$sample = factor(sub(x=peptides.filtered$R.Label, pattern="^(KL_\\w+_\\w+)_\\w+", replacement="\\1", perl=T))
  peptides.filtered$replicate = factor(sub(x=peptides.filtered$R.Label, pattern="^KL_\\w+_\\w+_(\\w+)", replacement="\\1", perl=T))
  peptides.filtered$fragment = with(peptides.filtered,  paste(F.FrgType, F.FrgNum,  sep="."))
  
  peptides.data = select(peptides.filtered, batch, sample, replicate, R.Label, EG.Label, EG.StrippedSequence, EG.Id, fragment, FG.PrecursorMz, FG.TotalPeakArea, F.PeakArea )
  
  #peptides.data$date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  
  peptides.data$batch_date = factor(experiment_map$date[match(peptides.data$R.Label, experiment_map$SampleName)])
  peptides.data$sampling_date = factor(dates_map$Sampling.date[match(peptides.data$sample, dates_map$Data_file_name)])
  peptides.data$processing_date = factor(dates_map$Processing.date[match(peptides.data$sample, dates_map$Data_file_name)])
  peptides.data$acquisition_date = factor(dates_map$Date.of.acquisition[match(peptides.data$sample, dates_map$Data_file_name)])
  
  
  
  
  file_name = "peptides.data.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(peptides.data,file=file_path)  
}



main = function() {
  prepare_peptides()
}


main()

