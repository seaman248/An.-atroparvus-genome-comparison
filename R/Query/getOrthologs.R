# manual https://mshadbolt.github.io/Ae_aegypti-toolset/AccessingVectorbaseinBiomaRt.html
if(!file.exists('./R/Query/output_data/orthologs.csv')){
  source('./R/Query/functions/config_ds.R')
  source('./R/Query/functions/getOrthologs.R')
  
  orthologs <- getOrthologs(sp_of_interest)
  
  write.csv2(orthologs, './R/Query/output_data/orthologs.csv', row.names = F)
  
  rm(list=ls())
}