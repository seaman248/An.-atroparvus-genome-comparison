library(dplyr)
library(ggplot2)
library(biomaRt)

report <- readLines('./R/GRIMM/output_data/blocks/report.txt')
results <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')
anchors <- read.table('./R/GRIMM/output_data/anchors/unique_coords.txt')


orthologs <- read.csv2('./R/Query/output_data/orthologs.csv')
sps <- c('alb', 'atr', 'gam')

scfs_coords <- lapply(sps, function(sp){
  query_dir <- './R/Query/output_data'
  sps_in_dir <- dir(query_dir)
  file <- paste0(query_dir, '/', sps_in_dir[grep(sp, sps_in_dir)])
  read.csv(file)
})

genes <- lapply(sps, function(sp){
  file <- paste0(c('./R/Clean/output_data/', sp, '_clean.csv'), collapse = '')
  read.csv2(file)
})

names(scfs_coords) <- names(genes) <- sps

sp1 <- 'alb'
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

atr_coords <- bind_rows(apply(big_blocks, 1, function(row){
  block_genes <- genes[[sp2]] %>%
    filter(start_position >= row[7], end_position <= row[8]) %>%
    filter(chromosome_name == getmode(chromosome_name)) %>%
    arrange(start_position)
  
  scfs <- scfs_coords[[sp2]][match(block_genes[, 1], scfs_coords[[sp2]][, 1]), c(2, 3, 4)]
  
  coords <- scfs %>%
    group_by(chromosome_name) %>%
    summarise(start = min(start_position), end = max(end_position), block = row[1]) %>%
    mutate(coords = paste0(chromosome_name, ':', start, '-',end)) %>%
    group_by(block) %>%
    summarise(coords = paste0(coords, collapse=', '))
}))

sp1_coords <- bind_rows(apply(big_blocks, 1, function(row){
  block_genes <- genes[[sp1]] %>%
    filter(start_position >= row[3], end_position <= row[4]) %>%
    filter(chromosome_name == getmode(chromosome_name)) %>%
    arrange(start_position) %>%
    summarise(chr = unique(chromosome_name), start = min(start_position), end = max(end_position), block = row[1]) %>%
    mutate(coords = paste0(chr, ':', start, '-', end))
  block_genes['coords']
}))

result_table <- bind_cols(atr_coords, sp1_coords, big_blocks)[, c(1, 2, 3, 13)]
colnames(result_table) <- c('block_id', 'atr_coords', 'gam_coords', 'anchors')

write.csv2(result_table, paste0('./R/Secondary_stat/output_data/', sp2, '-', sp1, '-big_blocks.csv'), quote = F)

# run
# source('./R/Clean/GRIMM_format.R')
# source('R/GRIMM/make_blocks.R')
