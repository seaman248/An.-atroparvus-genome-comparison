# clean orthologs table
# we need only genes from KI421924, KI421922, KI421894, KI421900

# define species of interest
sp_of_interest <- c('albimanus', 'atroparvus', 'gambiae', 'steph', 'funestus')
# load the data
orthologs <- read.csv2('./R/Query/output_data/orthologs.csv')
sps_genes <- lapply(sp_of_interest, function(sp){
  path <- paste0('./R/Query/output_data/', sp, '_genes.csv')
  read.csv(path, stringsAsFactors = F)
})
names(sps_genes) <- sp_of_interest

# Chose the atroparvus ids from KI421924, KI421922, KI421894, KI421900 atroparvus contig

atr_peric_genes_coords <- which(!is.na(match(sps_genes$atroparvus$chromosome_name, c('KI421924', 'KI421922', 'KI421894', 'KI421900'))))
atr_peric_ids <- sps_genes$atroparvus$ensembl_gene_id[atr_peric_genes_coords]

# Chose orthologs only for this ids
orth_coords <- which(!is.na(match(
  orthologs$aatroparvus_eg_gene, atr_peric_ids
)))

orthologs <- orthologs[orth_coords,]
colnames(orthologs) <- sp_of_interest
write.csv2(orthologs, './R/Secondary_stat/pericentromeric/output_data/peric_orthologs.csv', row.names=F)

lapply(sp_of_interest, function(sp){
  coords <- which(!is.na(match(sps_genes[[sp]]$ensembl_gene_id, orthologs[, sp])))
  table <- sps_genes[[sp]][coords,]
  path <- paste0('./R/Secondary_stat/pericentromeric/output_data/', sp, '_peric_dirty.csv')
  write.csv2(table, path, row.names = F)
})


rm(list=ls())
