library(biomaRt)
library(dplyr)

vb_host <- 'biomart.vectorbase.org' 
genes_mart <- useMart(listMarts(host=vb_host)$biomart[1], host = vb_host, dataset = 'aatroparvus_eg_gene')

IRP <- c('PF00756', 'PF13417', 'PF00067', 'PF06512')
UCB <- read.table('./R/Secondary_stat/max_blocks.txt', header = T)
atr_genes <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'), mart = genes_mart)


atr_blocks <- apply(UCB, 1, function(ucb_row){
  
  coord_string <- paste0(c(ucb_row[-1], 1), collapse = ':')
  print(ucb_row[1])
  
  data.frame(
    block = ucb_row[1],
    getBM(
    attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'agambiae_eg_gene'),
    filters = c('chromosomal_region'),
    values = coord_string,
    mart = genes_mart
  ))
  
})

atr_blocks <- bind_rows(atr_blocks)

exp_mart <- useMart('expression', host=vb_host, dataset = 'vb')
dev_exp <- getBM(
  attributes = c('p-value', 'ensembl_gene_id'),
  filters = c('hyb_annotation_type', 'experiment_id'),
  values = list(
    c('DevelopmentalStage'),
    c(29)
  ),
  mart = exp_mart
)


