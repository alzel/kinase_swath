source("./R/boot.R")
rm(list=ls())
output_dir = "./R/objects"
dir.create(output_dir)
suffix = "_load_"

##################################################
###loading proteome data  
##################################################
load_proteome = function() {
  
  input_path = "./data/2015-03-05"
  message("Loading proteome")
  
  filesToProcess = dir(path=input_path, pattern = "KL_Try_Spectronaut_Batch*", recursive=F)
  filesToProcess = filesToProcess[grep(pattern="*.xls", x=filesToProcess)]
  pattern.p = "KL_Try_Spectronaut_Batch02_([0-9]+)_([A-Za-z0-9 ]+)_([0-9]+)_([0-9]+).xls"
  
  matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

  read_peptides = function(x) {
    file_name = paste(input_path,x[[1]], sep="/") 
    table = read.table(file_name, sep="\t", header=T)
    table$batch = factor(rep(as.integer(x[[2]]), nrow(table))) #Spectronaut batch
    table$type  = factor(rep(x[[3]], nrow(table)))
    #table$Sp_date =  factor(rep(x[[4]], nrow(table))) #Spectranaut date
    #table$Sp_time =  factor(rep(x[[5]], nrow(table))) 
    table$file =  factor(rep(x[[1]], nrow(table)))
    return(table)
  }
  
    
  file.list = lapply(matches, FUN=read_peptides)
  dataset.peptides.raw = do.call(rbind.data.frame, file.list)
  
  file_name = paste("data.frame.peptides.raw", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(dataset.peptides.raw,file=file_path)
  
}

load_batch_map = function() {
  experiment_map.raw <- read.delim("./data/2014-11-03/1411_Batches/KL_batches_JV_v01.txt")
  
  experiment_map.raw$date = sub(x=experiment_map.raw$batch, pattern="^(2014_[0-9_]+)_KL_Try.*?$", replacement="\\1", perl=T)
  experiment_map = experiment_map.raw
  
  file_name = paste("experiment.map", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(experiment_map,file=file_path)

}


load_sample_map = function() {
  sample_map.raw <- read.delim("./data/2014-12-23/KL_Sample_Code_v01.csv", sep=",", header=T)
  sample_map = sample_map.raw
  
  file_name = paste("sample_map", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(sample_map,file=file_path)
}

load_dates_map = function() {

  tmp = read.xlsx(file="./data/2015-01-24/Kinase List_FC_sortedProtNo_JV.xls", header=T, sheetName="Sheet2")
  tmp = tbl_df(tmp)
  View(tmp)
  dates_map = dplyr::select(tmp, ORF, gene, Type, Nr, Data_file_name, Sampling.date, Processing.date, Date.of.acquisition)
  
  file_name = paste("dates_map", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(dates_map,file=file_path)
}


load_acqusition_times = function() {
  
  KL_Acqusition_dates <- read.csv("./data/2015-02-24/KL_Acqusition_dates.csv")
  KL_Acqusition_dates$sample_name = sub(pattern="(.*).wiff", replacement="\\1", x=KL_Acqusition_dates$filename)
  KL_Acqusition_dates = KL_Acqusition_dates %>% arrange(AcquisitionDate)
  #KL_Acqusition_dates$acquisition_time.str = strptime(KL_Acqusition_dates$AcquisitionDate, "%Y-%m-%d %H:%M:%S")
  acqusition_times = KL_Acqusition_dates
  
  file_name = paste("acqusition_times", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(acqusition_times,file=file_path)  
  
}


load_protein_annotations = function() {
  tmp = read.csv(file="./data/2015-02-24/protein_annotation_expanded_uniprot_v01.csv", header=T)
  protein_annotations = tmp
  
  file_name = paste("protein_annotations", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(protein_annotations,file=file_path)
}

loadKEGG = function() {
  library("RCurl")
  
  my_splitter = function(response){
    return(matrix(unlist(strsplit(strsplit(x=response[1], "\n")[[1]], split="\\t")), ncol = 2, byrow=T) )
  }
  
  url_string = "http://rest.kegg.jp/list/pathway/sce"
  response = getURL(url=url_string)
  
  pathways = data.frame(my_splitter(response))
  names(pathways) = c("path", "description")
  
  
#   url_string = "http://rest.kegg.jp/conv/uniprot/sce"
#   response = getURL(url=url_string)
#   org2uniprot = data.frame(my_splitter(response))
#   names(org2uniprot) = c("org_id", "uniprot")
#   
#   url_string = "http://rest.kegg.jp/conv/ncbi-geneid/hsa"
#   response = getURL(url=url_string)
#   hsa2ncbi_geneid = data.frame(my_splitter(response))
#   names(hsa2ncbi_geneid) = c("hsa", "ncbi_geneid")
  
  url_string = "http://rest.kegg.jp/link/sce/pathway"
  response = getURL(url=url_string)
  
  pathway2orf = data.frame(my_splitter(response))
  names(pathway2orf) = c("pathway", "kegg_ORF")
    
  pathway2orf$ORF = factor(sub(x=pathway2orf$kegg_ORF, pattern="sce:(\\w+)", replacement="\\1"))
  pathway2orf$description = pathways$description[match(pathway2orf$pathway, pathways$path)]

  ##modules
  url_string = "http://rest.kegg.jp/list/module/sce"
  response = getURL(url=url_string)
  
  modules = data.frame(my_splitter(response))
  names(modules) = c("md", "description")
  

  
  url_string = "http://rest.kegg.jp/link/sce/module"
  response = getURL(url=url_string)
  module2orf = data.frame(my_splitter(response))
  
  names(module2orf) = c("md", "kegg_ORF")
  module2orf$ORF = factor(sub(x=module2orf$kegg_ORF, pattern="sce:(\\w+)", replacement="\\1"))
  
  module2orf$description = modules$description[match(module2orf$md, modules$md)]

  file_name = paste("pathway2orf", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(pathway2orf, file=file_path)
  
  file_name = paste("pathways", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(pathways, file=file_path)

  
  file_name = paste("modules", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(modules, file=file_path)
  
  file_name = paste("module2orf", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(module2orf, file=file_path)
  
}


loadGO_slim = function() {
  GO_slim.raw = tbl_df(read.delim("./data/2015-03-19/go_slim_mapping.tab", header=F))
  GO_slim.selected =  GO_slim.raw %>% dplyr::select(V6,V1,V5,V4)
  
  names(GO_slim.selected) = c("pathway", "ORF", "description", "type")
  
  GO_slim.compartment =  GO_slim.selected %>% filter(type=="C")
  GO_slim.process     =  GO_slim.selected %>% filter(type=="P")
  GO_slim.function    =  GO_slim.selected %>% filter(type=="F")
  
  
  
  file_name = paste("GO_slim.compartment", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(GO_slim.compartment, file=file_path)
  
  file_name = paste("GO_slim.process", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(GO_slim.process, file=file_path)
  
  file_name = paste("GO_slim.function", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(GO_slim.function, file=file_path)
  
}

loadGO = function() {
  GO.raw = tbl_df(read.delim("./data/2015-03-19/gene_association.sgd", header=F, comment.char="!"))
  GO.raw.f = GO.raw %>% filter(V15 == "SGD")
  
  GO.raw.f %>% dplyr::select(V5, V2, V10, V9)
  
  GO_slim.selected =  GO_slim.raw %>%select(V6,V1,V5,V4)
  
  names(GO_slim.selected) = c("pathway", "ORF", "description", "type")
  
  GO_slim.compartment =  GO_slim.selected %>% filter(type=="C")
  GO_slim.process     =  GO_slim.selected %>% filter(type=="P")
  GO_slim.function    =  GO_slim.selected %>% filter(type=="F")
  
  file_name = paste("GO_slim.compartment", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(GO_slim.compartment, file=file_path)
  
  file_name = paste("GO_slim.process", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(GO_slim.process, file=file_path)
  
  file_name = paste("GO_slim.function", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(GO_slim.function, file=file_path)
  
}

load_metabolites = function() {
  metabolites.raw = read.table("./data/2014-03-05/KL_screen_normalized_PPP_results.csv", header=T, sep=",")
  file_name = paste("metabolites.raw", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(metabolites.raw, file=file_path)
}

load_gene_annotations = function() {
  gene.annotations <- read.delim("./data/2015-04-26/dbxref.tab", header=F)
  
  file_name = paste("gene.annotations", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(gene.annotations, file=file_path)
  
}

main = function() {
  #load_proteome()
  
  load_batch_map()
  load_sample_map()
  load_dates_map()
  load_acqusition_times()
  load_protein_annotations()
  load_gene_annotations()
  load_metabolites()
  loadKEGG()
}

main()




