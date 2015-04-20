rm(list=ls())
source("./R/functions.R")
figures_dir = "./figures"
output_dir = "./R/objects"
dir.create(figures_dir)
dir.create(output_dir)



main = function() {
  
  # ---- data_explore  ----
  source("./R/src/explore.R")
    
  source("./R/src/get_peptides.R")
  source("./R/src/peptides_correlations.R")
  source("./R/src/batch_effects.R")
  
  
  
}
