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

# dev_exp <- getBM(
#   attributes = c('p-value', 'ensembl_gene_id'),
#   filters = c('hyb_annotation_type', 'experiment_id'),
#   values = list(
#     c('DevelopmentalStage'),
#     c(29)
#   ),
#   mart = exp_mart
# )


find_coords <- function(apx, scf){
  atr_scfs <- read.csv2('./R/Clean/input_data/atr_order.csv')
  length <- atr_scfs[atr_scfs[,2] == scf, 4]
  start <- UCB[UCB$scf == scf, 3]
  end <- UCB[UCB$scf == scf, 4]
  data.frame(
    start = apx * start / length,
    end = apx * end / length
  )
}

agam_scf_loc <- read.csv('./data/agam_scfs_location.csv', stringsAsFactors = F)
agam_scf_loc <- agam_scf_loc[agam_scf_loc[, 2] != '' & agam_scf_loc[, 2] != 'Unknown', ]
agam_scf_loc[, 1] <- paste0('supercont1.', agam_scf_loc[, 1])

apply(UCB, 1, function(row){
  df <- getBM(
    attributes = c('aaegypti_eg_chromosome'),
    filters = c('chromosomal_region'),
    values = paste0(c(row[-1], 1), collapse = ':'),
    mart = genes_mart
  )
  res <- table(agam_scf_loc[which(!is.na(match(agam_scf_loc[, 1], unique(df[, 1])))), 2])
  res['block'] <- row[1]
  res
})

