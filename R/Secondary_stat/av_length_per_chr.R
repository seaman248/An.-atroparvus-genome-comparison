library(dplyr)
library(ggplot2)

source('./R/Visualize/functions/parse_blocks.R')

blocks <- parse_blocks(
  read.table('./R/GRIMM/output_data/blocks/blocks.txt')
)

blocks_length <- lapply(blocks, function(blocks_sp){
  t <- as.list(tapply(blocks_sp$end, blocks_sp$chr, sd))
})

blocks_length <- as.data.frame(bind_rows(blocks_length))/1000

blocks_length[4,] <- round(colMeans(blocks_length), digits=2)
