# manual https://mshadbolt.github.io/Ae_aegypti-toolset/AccessingVectorbaseinBiomaRt.html
library(biomaRt)

source('./R/Query/functions/config_ds.R')
source('./R/Query/functions/find_ds.R')
source('./R/Query/functions/getGenesFromMart.R')

getOrthologs <- function(species){
  
  true_sp_names <- unlist(lapply(species, find_ds))
  
  orths <- getGenesFromMart(dsName=true_sp_names[1], attributes = c('ensembl_gene_id', true_sp_names[2:length(true_sp_names)]))
  
  names(orths) <- true_sp_names
  
  return(orths[!apply(orths, 1, function(x){any(x=='')}),])
}