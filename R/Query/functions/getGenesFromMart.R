library(biomaRt)

source('./R/Query/functions/config_ds.R')
source('./R/Query/functions/find_ds.R')

getGenesFromMart <- function(
  dsName, # name of dataset _ example : aatroparvus_eg_gene
  attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position','strand')
){
  trueDSName <- find_ds(dsName)
  
  thisMart <- useMart(
    gbase,
    host=vb_host,
    dataset=trueDSName
  )
  
  return(getBM(attributes=attributes, mart=thisMart))
  
}