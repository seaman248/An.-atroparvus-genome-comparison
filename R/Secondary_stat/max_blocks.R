atr_blocks <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')[,c(6:9)]
atr_genes <- read.csv2('./R/Clean/output_data/atr_clean.csv')
atr_genes_mart <- read.csv('./R/Query/output_data/atroparvus_genes.csv')

topmost_blocks <- atr_blocks[atr_blocks[,3] > 1000000, ]

conservative_genes <- apply(topmost_blocks, 1, function(row){
  end <- as.numeric(row[2]) + as.numeric(row[3])
  genes <- atr_genes[atr_genes$chromosome_name==row[1] & atr_genes$start_position > row[2] & atr_genes$end_position < end,]
})

top_scfs <- lapply(conservative_genes, function(genset){
  n <- which(!is.na(match(atr_genes_mart$ensembl_gene_id, genset$ensembl_gene_id)))
  genes <- atr_genes_mart[n, ]
  paste0(
    '(', as.character(unique(genset[,2])), ')',
    as.character(unique(genes[,2])), ':',
    min(genes[,3]),'-',
    max(genes[,4])
  )
})

unname(unlist(top_scfs))
