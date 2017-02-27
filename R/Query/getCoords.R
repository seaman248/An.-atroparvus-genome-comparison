source('./R/Query/functions/config_ds.R')
source('./R/Query/functions/getGenesFromMart.R')

if(!file.exists('./R/Query/output_data/gambiae_genes.csv')){
  lapply(sp_of_interest, function(sp){
    geneTable <- getGenesFromMart(sp)
    pathToFile <- paste0('./R/Query/output_data/', sp, '_genes.csv')
    write.csv(geneTable, pathToFile, row.names = FALSE)
    return(pathToFile)
  })
  rm(list=ls())
}