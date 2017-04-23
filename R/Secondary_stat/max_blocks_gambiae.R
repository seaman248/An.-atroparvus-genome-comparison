library(biomaRt) 
library(dplyr)

UCBs <- read.table('./R/Secondary_stat/max_blocks.txt', header = T)

vb_host <- 'biomart.vectorbase.org' 
genes_mart <- useMart(listMarts(host=vb_host)$biomart[1], host = vb_host, dataset = 'aatroparvus_eg_gene')


gam_ids <- apply(UCBs, 1, function(row){
  
  coord <- paste0(c(row[2], row[3], row[4]), collapse = ':')
  getBM(
    attributes = c('agambiae_eg_gene'),
    filters = c('chromosomal_region'),
    values = c(coord),
    mart = genes_mart
  )
})

genes_mart <- useMart(listMarts(host=vb_host)$biomart[1], host = vb_host, dataset = 'agambiae_eg_gene')

gam_bands <- lapply(gam_ids, function(gam_id_s){
  getBM(
    attributes = c('ensembl_gene_id', 'band'),
    filters = c('ensembl_gene_id'),
    values = c(gam_id_s),
    mart = genes_mart
  )
})

names(gam_bands) <- UCBs[, 1]

lapply(names(gam_bands), function(n){
  df <- gam_bands[[n]]
  bands <- df[, 2]
  bt <- table(bands)
  bt['block'] <- n
  bt['angenes'] <- UCBs[UCBs$block == n, 5]
  bt
})