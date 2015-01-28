source("./R/boot.R")

##################################################
###loading proteome data  
##################################################
load_proteome = function() {
  
  input_path = "./data/2015-01-24/"
  message("Loading proteome")
  
  filesToProcess = dir(path=input_path, pattern = "KL_Try_Spectronaut_Batch_*", recursive=F)
  filesToProcess = filesToProcess[grep(pattern="*.xls", x=filesToProcess)]
  pattern.p = "KL_Try_Spectronaut_Batch_([0-9]+)_(\\d+_[A-Za-z0-9 ]+)_([0-9]+)_([0-9]+).xls"
  matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

  read_peptides = function(x) {
    file_name = paste(input_path,x[[1]], sep="/") 
    table = read.table(file_name, sep="\t", header=T)
    table$batch = factor(rep(as.integer(x[[2]]), nrow(table)))
    table$type  = factor(rep(x[[3]], nrow(table)))
    table$date =  factor(rep(x[[4]], nrow(table)))
    table$time =  factor(rep(x[[5]], nrow(table)))
    table$file =  factor(rep(x[[1]], nrow(table)))
    return(table)
  }
  
    
  file.list = lapply(matches, FUN=read_peptides)
  dataset.peptides.raw = do.call(rbind.data.frame, file.list)
  
  file_name = "data.frame.peptides.raw.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(dataset.peptides.raw,file=file_path)
  
}

load_batch_map = function(x) {
  experiment_map.raw <- read.delim("./data/2014-11-03/1411_Batches/KL_batches_JV_v01.txt")
  experiment_map.raw$date = sub(x=experiment_map.raw$batch, pattern="^(2014_[0-9_]+)_KL_Try.*?$", replacement="\\1", perl=T)
  experiment_map = experiment_map.raw
  
  file_name = "experiment.map.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(experiment_map,file=file_path)

}


load_sample_map = function(x) {
  sample_map.raw <- read.delim("./data/2014-12-23/KL_Sample_Code_v01.csv", sep=",", header=T)
  sample_map = sample_map.raw
  
  file_name = "sample.map.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(sample_map,file=file_path)
}

load_dates_map = function(x) {
  tmp = read.xlsx(file="./data/2015-01-24/Kinase List_FC_sortedProtNo_JV.xls", header=T, sheetName="Sheet2")
  dates_map = select(tmp, Data_file_name, Sampling.date, Processing.date, Date.of.acquisition)
  dates_map
  file_name = "dates_map.RData"
  file_path = paste(output_dir, file_name, sep="/")
  save(dates_map,file=file_path)
}

main = function() {
  load_proteome()
  load_batch_map()
  load_sample_map()
  load_dates_map()
}

main()




