lapply(c('genoPlotR', 'dplyr'), require, character.only=T)

source('./R/Visualize/functions/parse_blocks.R')

MIN_BLOCK_SIZE = 1000000


# rows_to_bind <- list(
#   c(334, 335, 336, 337), # e2 e4  e4
#   c(768, 773, 774), # e5 e3  e3
#   c(769, 770, 771, 772, 775) # e5 e4  e3
# )

blocks_file <- read.table('./R/GRIMM/output_data/blocks/blocks.txt', row.names = 1)

V3_to_remove <- unlist(c(
V3_alb_overlaps <- blocks_file %>%
  mutate(alb_end = V3 + V4) %>%
  arrange(V3) %>%
  mutate(er_alb = c(alb_end[1:length(alb_end)-1] - V3[2:length(V3)], 0)) %>%
  filter(er_alb > 0) %>%
  dplyr::select(V3)
,
V3_atr_overlaps <- blocks_file %>%
  mutate(atr_end = V7 + V8)%>%
  arrange(V7) %>%
  mutate(er_atr = c(atr_end[1:length(atr_end)-1] - V7[2:length(V7)], 0)) %>%
  filter(er_atr > 0) %>%
  dplyr::select(V3)
,
V3_gam_overlaps <- blocks_file %>%
  mutate(gam_end = V11 + V12) %>%
  arrange(V11) %>%
  mutate(er_gam = c(gam_end[1:length(gam_end)-1] - V11[2:length(V11)], 0)) %>%
  filter(er_gam > 0) %>%
  dplyr::select(V3)
))

rows_to_remove <- c(match(V3_to_remove, blocks_file$V3), 410, 335, 336, 334, 773, 774, 770, 771, 772, 775)

# blocks_file2 <- read.table('./data/blocks.txt', row.names = 1)
# rows_to_remove <- as.numeric(row.names(blocks_file[(blocks_file[, 1] != blocks_file[, 5]) | (blocks_file[, 5] != blocks_file[, 9]), ]))
# blocks_file2[(blocks_file2[, 1] != blocks_file2[, 5]) | (blocks_file2[, 5] != blocks_file2[, 9]), ]


blocks <- parse_blocks(read.table('./R/GRIMM/output_data/blocks/blocks.txt')[-rows_to_remove, ])

names(blocks) <- c('alb', 'atr', 'gam')

c_left_arms <- list(
  alb = c(1, -1, 1, 1, -1),
  atr = c(1, -1, 1, -1, 1),
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
    gene_type = 'blocks',
    chr = sp_blocks$chr,
    cex = 0.0001,
    lwd = 0.0001
    
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
  bold <- 10000
  
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
tiff('./output/full.tiff', width = 5000, height = 1000, units = 'px', pointsize = 40)

plot_gene_map(
  dna_segs=seqs,
  xlims = xlims,
  comparisons = comparisons
)

dev.off()

# rm(list=ls())
