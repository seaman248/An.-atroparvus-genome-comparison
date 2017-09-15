lapply(c('genoPlotR', 'dplyr'), require, character.only=T)

source('./R/Visualize/functions/parse_blocks.R')

MIN_BLOCK_SIZE = 1000000

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

V3_intra_chromosome_translocations <- blocks_file %>%
  filter(V2 != V6 | V6 != V10 | V2 != V10) %>%
  dplyr::select(V3)

rows_to_remove <- c(match(V3_to_remove, blocks_file$V3), 410)

# blocks_file2 <- read.table('./data/blocks.txt', row.names = 1)
# rows_to_remove <- as.numeric(row.names(blocks_file[(blocks_file[, 1] != blocks_file[, 5]) | (blocks_file[, 5] != blocks_file[, 9]), ]))
# blocks_file2[(blocks_file2[, 1] != blocks_file2[, 5]) | (blocks_file2[, 5] != blocks_file2[, 9]), ]

blocks_file <- blocks_file[-rows_to_remove, ]
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
    mutate(end = start + end) # length to end coordinates
  
  sp_blocks <- data.frame(sp_blocks)
  
  sp_blocks$strand <- sp_blocks$strand * blocks$atr$strand # atroparvus strand as basic


  seg <- dna_seg(data.frame(
    name = row.names(sp_blocks),
    start = sp_blocks$start,
    end = sp_blocks$end,
    strand = sp_blocks$strand,
    col = 'black',
    gene_type = 'blocks',
    chr = sp_blocks$chr,
    lwd = 0.05
    
  ))
  seg
})


names(seqs) <- names(blocks)

xlims <- lapply(names(blocks), function(sp){
  el_to_reverse <- which(c_left_arms[[sp]] == -1)

  lims_table <- seqs[[sp]] %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(start = min(start), end = max(end)) %>%
    dplyr::select(chr, start, end)

  lims_table[el_to_reverse, c(2, 3)] <- lims_table[el_to_reverse, c(3, 2)]

  as.vector(as.matrix(t(lims_table[, c(2, 3)])))
})

comparisons <- lapply(1:(length(seqs)-1), function(n){
  bold <- 10000

  start1 <- middle(seqs[[n]]) - bold/2
  end1 <- middle(seqs[[n]]) + bold/2
  start2 <- middle(seqs[[n+1]]) - bold/2
  end2 <- middle(seqs[[n+1]]) + bold/2

  # start1 <- seqs[[n]]$start
  # end1 <- seqs[[n]]$start + bold
  # start2 <- seqs[[n+1]]$start
  # end2 <- seqs[[n+1]]$start + bold

  comp <- data.frame(
    start1 = start1, end1 = end1,
    start2 = start2, end2 = end2
  )
  comp$col <- 'blue'
  comp$col[seqs[[n]]$strand != seqs[[n+1]]$strand] <- 'red'
  comparison(comp)
})


# # visualize
tiff('./output/full.tiff', width = 5000, height = 1000, units = 'px', res = 600, pointsize = 4)

plot_gene_map(
  dna_segs=seqs,
  comparisons = comparisons,
  xlims = xlims
)

dev.off()

# rm(list=ls())


# Зависимость длинны синтенного блока от удаленности от теломеры
names(xlims) <- names(blocks)

distance_from_telomere_table <- lapply(names(xlims), function(sp){
    lims <- xlims[[sp]]
    
    coords <- bind_rows(lapply(seq(1, length(lims), 2), function(n){
      data.frame(el = paste0('e', (n+1)/2), centromere = lims[n], telomere = lims[n+1])
    }))
    
    bind_rows(apply(coords, 1, function(row){
      blocks_from_el <- blocks[[sp]] %>% filter(chr == row[1])
      
      data.frame(
        sp = sp,
        el = row[1],
        distance_from_telomere = abs(as.numeric(row[3]) - blocks_from_el$start),
        length = blocks_from_el$end
      )
      
    }))
})

names(distance_from_telomere_table) <- names(blocks)

tiff('./output/dist-length_cor.tiff', width = 3500, height = 1000, units = 'px', res = 600, pointsize = 4)

par(mfrow=c(1, 3))
lapply(names(distance_from_telomere_table), function(sp){
  
  dist_table <- distance_from_telomere_table[[sp]]
  quant_99 <- quantile(dist_table$length, probs = c(0.99)) # 99 percentile
    
  tab <- dist_table %>%
    mutate(col = ifelse(length > quant_99, 'blue', 'black'), pch = ifelse(length > quant_99, 19, 1))
  tab_99 <- tab %>% filter(length > quant_99)
  
  coeff <- coef(lm(length~distance_from_telomere, data = tab))
  coeff_99 <- coef(lm(length~distance_from_telomere, data = tab%>%filter(length > quant_99)))
  
  cor <- round(cor(tab$length, tab$distance_from_telomere), 2)
  cor_99 <- round(cor(tab_99$length, tab_99$distance_from_telomere), 2)
  
  
  plot(tab$distance_from_telomere, tab$length, xlab = 'Distance from telomere, bp', ylab = 'Block length, bp', main = paste0(sp), lwd = 0.5, col = tab$col, pch = tab$pch)
  abline(coef = coeff, col = 'red', lwd = 1)
  abline(coef = coeff_99, col = 'blue')
  text(x = max(tab$distance_from_telomere), y = max(tab$length), pos = 2, labels = paste0('r = ', cor), col = 'red')
  text(x = max(tab$distance_from_telomere), y = max(tab$length)-120000, pos = 2, labels = paste0('r = ', cor_99), col = 'blue')
  
})

dev.off()


# Comparison with random model
random_breaks <- floor(runif(nrow(blocks$atr)*2, min=0, max = max(blocks$atr$start)))
random_breaks <- random_breaks[order(random_breaks)]

random_blocks<- bind_rows(lapply(1:(length(random_breaks)-1), function(m){
  data.frame(
    # start = random_breaks[n],
    # end = random_breaks[n+1],
    length = random_breaks[m+1] - random_breaks[m]
  )
}))$length

random_rows_to_remove <- sample(1:nrow(blocks_file), replace = F)
random_blocks <- random_blocks[random_rows_to_remove]
wilcox_p <- wilcox.test(random_blocks, blocks$atr$end)

random_99 <- quantile(random_blocks, probs = .99)
atr_99 <- quantile(blocks$atr$end, probs = .99)

# tiff('./output/compare_with_random_model.tiff', width = 1000, height = 1000, units = 'px', res = 200, pointsize = 4)
ggplot() +
  geom_histogram(aes(random_blocks, y=..count../sum(..count..)), fill = 'white', col = 'black', alpha = .5) +
  geom_vline(aes(xintercept = random_99)) +

  geom_histogram(aes(blocks$atr$end, y=..count../sum(..count..)), fill = 'red', col = 'red', alpha = .5, binwidth = 30000) +
  geom_vline(aes(xintercept = atr_99), col = 'red')+

  ylab('Density') +
  xlab('Blocks length') +
# dev.off()




### Create table of translocations
translocations <- blocks_file %>% filter(V2!=V6 | V6!=V10 | V2!=V10)
colnames(translocations) <- c('alb_el', 'alb_start', 'alb_length', 'alb_strand', 'atr_el', 'atr_start', 'atr_length', 'atr_strand', 'gam_el', 'gam_start', 'gam_length', 'gam_strand')

translocation_rows <- match(translocations$alb_start, blocks_file$V3)

translocations <- bind_cols(lapply(names(xlims), function(sp){
  lims <- xlims[[sp]]
  
  coords <- bind_rows(lapply(seq(1, length(lims), 2), function(n){
    data.frame(el = paste0('e', (n+1)/2), centromere = lims[n])
  }))
  translocation_blocks <- blocks[[sp]] %>%
    slice(translocation_rows) %>%
    mutate(centromere_coord = coords$centromere[match(chr, coords$el)]) %>%
    mutate(from_cent = abs(centromere_coord - start))%>%
    dplyr::select(chr, end, strand, from_cent) %>%
   setNames(c(chr = paste0(sp, '_chr'), end = paste0(sp, '_length'), strand = paste0(sp, '_strand'), from_cent = paste0(sp, '_bp-from-cent')))
  
}))

write.csv2(translocations, './output/intrachr.csv')
