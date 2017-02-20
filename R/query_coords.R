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

### Convert to list of data.frames

gene_coords.frames <- lapply(sps, function(sp){
  geneSet <- data.frame(do.call('cbind', gene_coords[[sp]]), stringsAsFactors = FALSE)
  geneSet.new <- data.frame(matrix(NA, ncol=5, nrow = nrow(orths)))
  
  geneSet.new[match(geneSet[,1], orths[,sp]),] <- geneSet
  geneSet.new
})

### Convert to one data frame to remove NA's

gene_coords.data_frame <- na.omit(do.call('cbind', gene_coords.frames))

### Save result data.frame

write.csv(gene_coords.data_frame, './data/processed/alb_atr_gam_sin(with coords).csv', row.names = FALSE)
