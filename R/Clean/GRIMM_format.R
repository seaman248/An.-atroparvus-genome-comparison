library(dplyr)

# load clean table for every sp
active_sps <- c('alb', 'atr', 'gam')

genes.list <- lapply(active_sps, function(sp){
  data.folder <- './R/Clean/output_data/'
  output_data.files <- list.files(data.folder)
  genes.csv <- output_data.files[grep(sp, output_data.files)]
  read.csv2(paste0(data.folder, '/', genes.csv))
})

# load table of orthologs
orths <- read.csv2('./R/Query/output_data/orthologs.csv')
colnames(orths) <- c('alb', 'atr', 'gam')
names(genes.list) <- active_sps

# combine tables into one
genes.trueorder <- bind_cols(lapply(active_sps, function(sp_name){
  order <- orths[,sp_name]
  genes <- genes.list[[sp_name]]
  
  true_order <- match(order, genes$ensembl_gene_id)
  genes <- genes[true_order, ]
  
  # change end column to length
  genes[,4] <- genes[,4] - genes[,3]
  
  return(genes)
}))

genes.trueorder <- na.omit(genes.trueorder)

columns <- unlist(lapply(1:length(active_sps), function(n){
  end <- n * 5
  start <- end - 3
  start:end
}))

GRIMM.table <- genes.trueorder[,c(
  1,
  columns
)]

write.table(GRIMM.table, './R/Clean/output_data/GRIMM.txt', row.names = F, quote = F)

# rm(list=ls())
