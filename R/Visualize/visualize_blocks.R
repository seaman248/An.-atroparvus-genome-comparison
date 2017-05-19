lapply(c('genoPlotR', 'dplyr'), require, character.only=T)

source('./R/Visualize/functions/parse_blocks.R')

MIN_BLOCK_SIZE = 1000000

blocks <- parse_blocks(read.table('./R/GRIMM/output_data/blocks/blocks.txt'))

names(blocks) <- c('alb', 'atr', 'gam')

c_left_arms <- list(
  alb = c(1, -1, 1, 1, -1),
  atr = c(-1, -1, 1, -1, 1),
  gam = c(1, -1, 1, -1, 1)
)

c_left_arms <- lapply(c_left_arms, function(orientation_vector){
  names(orientation_vector) <- paste0('e', 1:5)
  orientation_vector
})

seqs <- lapply(names(blocks), function(sp){
  sp_blocks <- blocks[[sp]]
  
  sp_blocks <- sp_blocks %>%
    mutate(end = start + (end * strand)) # length to end coordinates
  
  sp_blocks <- data.frame(sp_blocks)
  
  sp_blocks$strand <- sp_blocks$strand * blocks$atr$strand # atroparvus strand as basic


  dna_seg(data.frame(
    name = row.names(sp_blocks),
    start = sp_blocks$start,
    end = sp_blocks$end,
    strand = sp_blocks$strand,
    col = 'black',
    gene_type = 'side_blocks',
    chr = sp_blocks$chr
  ))
})

names(seqs) <- names(blocks)

xlims <- lapply(names(blocks), function(sp){
  el_to_reverse <- which(c_left_arms[[sp]] == -1)

  lims_table <- seqs[[sp]] %>%
    group_by(chr) %>%
    summarise(start = min(start), end = max(end)) %>%
    select(chr, start, end)
  
  lims_table[el_to_reverse, c(2, 3)] <- lims_table[el_to_reverse, c(3, 2)]
  # 
  as.vector(as.matrix(t(lims_table[, c(2, 3)])))
})

comparisons <- lapply(1:(length(seqs)-1), function(n){
  comp <- data.frame(
    start1 = seqs[[n]]$start, end1 = seqs[[n]]$end,
    start2 = seqs[[n+1]]$start, end2 = seqs[[n+1]]$end
  )
  comp$col <- 'grey'
  comparison(comp)
})

# names(blocks) <- c('alb', 'atr', 'gam')
# 
# # create seq
# seqs <- lapply(blocks, create_seqs)
# 
# # create xlims
# 
# xlims <- lapply(names(blocks), function(sp){
#   coords <- blocks[[sp]] %>%
#     group_by(chr) %>%
#     summarise(start = min(start), end = as.numeric(max(end)))
#   
#   coords <- as.data.frame(coords)
#   
#   rows_to_change <- arm_strand[[sp]]
# 
#   coords[rows_to_change, c(2, 3)] <- coords[rows_to_change, c(3, 2)]
# 
#   as.vector(t(as.matrix(coords[, c(2, 3)])))[c(9, 10)]
# })
# 
# names(xlims) <- names(blocks)
# 
# # comparisons
# comparisons <- make_comparisons(seqs)
# 
# comparisons <- lapply(1:(length(blocks)-1), function(n){
#   coords <- data.frame(
#     start1 = blocks[[n]]$start, end1 = blocks[[n]]$end,
#     start2 = blocks[[n+1]]$start, end2 = blocks[[n+1]]$end
#   )
#   coords$col <- 'grey'
#   coords$col[blocks[[n]]$chr != blocks[[n+1]]$chr] <- 'red'
#   coords$direction <- 1
#   comparison(coords)
# })
# 
# # visualize

plot_gene_map(
  dna_segs=seqs,
  xlims = xlims,
  comparisons = comparisons
)

# rm(list=ls())
