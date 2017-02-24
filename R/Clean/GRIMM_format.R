library(dplyr)

# load clean table for every sp
genes.list <- lapply(c('alb', 'atr', 'gam'), function(sp){
  data.folder <- './R/Clean/output_data/'
  output_data.files <- list.files(data.folder)
  genes.csv <- output_data.files[grep(sp, output_data.files)]
  read.csv2(paste0(data.folder, '/', genes.csv))
})

# load table of orthologs
orths <- read.csv2('./R/Query/output_data/orthologs.csv')

names(genes.list) <- names(orths)

# combine tables into one
genes.trueorder <- bind_cols(lapply(names(orths), function(sp_name){
  order <- orths[,sp_name]
  genes <- genes.list[[sp_name]]
  
  true_order <- match(order, genes$ensembl_gene_id)
  genes[true_order, ]
}))

# rm excess
rm(orths, genes.list)

genes.trueorder <- na.omit(genes.trueorder)

GRIMM.table <- genes.trueorder[,c(1:5, 7:10, 12:15)]

write.table(GRIMM.table, './R/Clean/output_data/GRIMM.txt', row.names = F, quote = F)

rm(list=ls())
