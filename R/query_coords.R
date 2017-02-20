source('./R/martCFG.R')

orths <- read.csv2('./data/processed/orthsList.csv', stringsAsFactors = FALSE)
sps <- colnames(orths)

### Query gene coordinates

gene_coords <- lapply(sps, function(sp){
  
  attributes <- c('ensembl_gene_id','chromosome_name','start_position','end_position','strand')
  
  MART <- useMart(gbase, vb_host, dataset=sp)
  
  genes <- getBM(
    attributes = attributes,
    mart = MART,
    filters = 'ensembl_gene_id',
    values = orths[,sp]
  )
  
  return(as.list(genes))
  
})

names(gene_coords) <- sps

### Find genes in orths table

lapply(sps, function(sp){
  geneSet <- do.call('cbind', gene_coords[[sp]])[1:10,]
  
  
})