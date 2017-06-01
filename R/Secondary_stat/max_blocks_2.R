library(dplyr)
library(ggplot2)
library(biomaRt)

report <- readLines('./R/GRIMM/output_data/blocks/report.txt')
results <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')
anchors <- read.table('./R/GRIMM/output_data/anchors/unique_coords.txt')

atr_scfs_coords <- read.csv('./R/Query/output_data/atroparvus_genes.csv', stringsAsFactors = F)

orthologs <- read.csv2('./R/Query/output_data/orthologs.csv')

genes <- lapply(c('alb', 'atr', 'gam'), function(sp){
  file <- paste0(c('./R/Clean/output_data/', sp, '_clean.csv'), collapse = '')
  read.csv2(file)
})



names(genes) <- c('alb', 'atr', 'gam')

sp1 <- 'gam'
sp2 <- 'atr'

# find "block 1: 5 anchors"

length <- report[grep('block \\d+: \\d+ anchors', report)]

# extract number of anchors

length <- sub('block \\d+: ', '', length)
length <- as.numeric(sub(' anchors', '', length))

results$length <- length

quant <- quantile(results$length, 0.99)

big_blocks <- results %>%
  filter(length >= quant) %>%
  mutate(V8 = V7 + V8, V4 = V3+V4)

# find atr_genes

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

apply(big_blocks, 1, function(row){
  block_genes <- genes[[sp2]] %>%
    filter(start_position > row[7], end_position < row[8]) %>%
    filter(chromosome_name == getmode(chromosome_name)) %>%
    arrange(start_position)
  
  scfs <- atr_scfs_coords[match(block_genes[, 1], atr_scfs_coords[, 1]), c(2, 3, 4)]
  
  scfs
})


# run
# source('./R/Clean/GRIMM_format.R')
# source('R/GRIMM/make_blocks.R')
