# manual https://mshadbolt.github.io/Ae_aegypti-toolset/AccessingVectorbaseinBiomaRt.html
library(biomaRt)
#vb_host <- 'biomart.vectorbase.org'

find_ds <- function(sp_string){
  vb_host <- 'biomart.vectorbase.org'
  gbase <- 'vb_gene_mart_1612'
  sp_ds <- as.character(
    listDatasets(
      useMart(
        gbase, 
        host=vb_host
      )
    )$dataset
  )
  return(
    sp_ds[grep(sp_string, sp_ds)]
  )
}

getGenesFromMart <- function(
  geneIDs=c(''),
  filter=F,
  dsName, # name of dataset _ example : aatroparvus_eg_gene
  attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position','strand')
){
  vb_host <- 'biomart.vectorbase.org'
  gbase <- 'vb_gene_mart_1612'
  trueDSName <- find_ds(dsName)
  thisMart <- useMart(
    gbase,
    host=vb_host,
    dataset=trueDSName
  )
  if(filter){
    return(
      getBM(attributes=attributes, mart=thisMart, filters = 'ensembl_gene_id', values = geneIDs)
    )
  } else {
    getBM(attributes=attributes, mart=thisMart)
  }
  
}


getOrthologs <- function(species){
  vb_host <- 'biomart.vectorbase.org'
  gbase <- 'vb_gene_mart_1612'
  true_sp_names <- unlist(lapply(species, find_ds))
  orths <- getGenesFromMart(dsName=true_sp_names[1], attributes = c('ensembl_gene_id', true_sp_names[2:length(true_sp_names)]))
  names(orths) <- true_sp_names
  return(orths[!apply(orths, 1, function(x){any(x=='')}),])
}