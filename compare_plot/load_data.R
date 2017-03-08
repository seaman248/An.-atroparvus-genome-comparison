library(plyr)
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

# deduplicate orths by hands
orths <- as.data.frame(apply(orths, 2, function(col){
  col[duplicated(col) & duplicated(col, fromLast = T)] <- NA
  col
}))

orths <- na.omit(orths)


# combine tables into one
genes.trueorder <- bind_cols(lapply(names(orths), function(sp_name){
  order <- orths[,sp_name]
  genes <- genes.list[[sp_name]]
  
  true_order <- match(order, genes$ensembl_gene_id)
  genes <- genes[true_order, ]
  
  # change end column to length
  genes[,4] <- genes[,4] - genes[,3]
  
  genes[,2] <- mapvalues(genes[,2], from=paste0('e', 1:5), to=1:5)
  
  return(genes)
}))

genes.trueorder <- na.omit(genes.trueorder)

rm(orths, genes.list)
