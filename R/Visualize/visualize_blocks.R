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
  bold <- 200
  
  start1 <- middle(seqs[[n]]) - bold/2
  end1 <- middle(seqs[[n]]) + bold/2
  start2 <- middle(seqs[[n+1]]) - bold/2
  end2 <- middle(seqs[[n+1]]) + bold/2
  
  comp <- data.frame(
    start1 = start1, end1 = end1,
    start2 = start2, end2 = end2
  )
  comp$col <- 'blue'
  comp$col[seqs[[n]]$strand != seqs[[n+1]]$strand] <- 'red'
  comparison(comp)
})


# # visualize
# tiff('./output/full.tiff', width = 5000, height = 1000, units = 'px', pointsize = 40)

plot_gene_map(
  dna_segs=seqs,
  xlims = xlims,
  comparisons = comparisons
)

# dev.off()

# rm(list=ls())
