library(dplyr)


###


atr_blocks <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')[,c(6:9)]
atr_grimm <- read.table('./R/Clean/output_data/GRIMM.txt', stringsAsFactors = F)[, c(1, 6:9)]
orthologs <- read.csv2('./R/Query/output_data/orthologs.csv', stringsAsFactors = F)

alb_genes <- read.csv('./R/Query/output_data/albimanus_genes.csv')
atr_genes <- read.csv('./R/Query/output_data/atroparvus_genes.csv')
gam_genes <- read.csv('./R/Query/output_data/gambiae_genes.csv')

topmost_blocks <- atr_blocks[atr_blocks[,3] > 1000000, ]

atr_grimm <- transform(atr_grimm, V7 = as.numeric(V7), V8 = as.numeric(V8), V9 = as.numeric(V9))



###


# Находим какие ортологи генов альбимануса задействованы в крупных блоках
alb_ids_set <- apply(topmost_blocks, 1, function(row){
  
  start = as.numeric(row[2])
  end = start + as.numeric(row[3])
  
  alb_ids <- filter(atr_grimm, V7 >= start, V7 <= end)[, 1]
})

# Находим по таблице ортологов соответствующие гены у An. atroparvus
atr_genes_blocks <- lapply(alb_ids_set, function(alb_ids){
  ms_orth <- match(alb_ids, orthologs[, 1])
  atr_ids <- orthologs[ms_orth, 2]
  ms_genes <- match(atr_ids, atr_genes$ensembl_gene_id)
  atr_genes[ms_genes, ]
})

# Находим точные координаты консервативных блоков у An. atroparvus

block_coords <- bind_rows(lapply(names(atr_genes_blocks), function(n){
  block_genes <- atr_genes_blocks[[n]]
  scfs <- table(as.character(block_genes[, 2]))
  scfs <- scfs[scfs>2]
  bind_rows(lapply(names(scfs), function(scf){
    data.frame(
      block = n,
      scf = scf,
      start = min(block_genes[block_genes[, 2] == scf, 3]),
      end = max(block_genes[block_genes[, 2] == scf, 4])
    )
  }))
}))

write.table(block_coords, file='./R/Secondary_stat/max_blocks.txt', row.names = F, quote = F)
