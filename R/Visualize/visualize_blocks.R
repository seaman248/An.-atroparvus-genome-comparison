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

tiff('./output/dist-length_cor-2.tiff', width = 5000, height = 2500, units = 'px', res = 650, pointsize = 6)

par(mfrow=c(3, 5))
lapply(names(distance_from_telomere_table), function(sp){
  
  quant_99 <- quantile(distance_from_telomere_table[[sp]]$length, probs = c(.99)) # 99 percentile

  lapply(paste0('e', c(1:5)), function(elt){
    
    
    tab <- distance_from_telomere_table[[sp]] %>%
      filter(el == elt) %>%
      mutate(col = ifelse(length > quant_99, 'red', 'grey'), pch = ifelse(length > quant_99, 17, 20))

    tab_99 <- tab %>% filter(length > quant_99)

    coeff <- coef(lm(length~distance_from_telomere, data = tab))

    cor <- round(cor(tab$length, tab$distance_from_telomere), 2)
    
    cent <- max(tab$distance_from_telomere)
    xaxis_labels <- c(0, cent*0.25, cent*0.5, cent*0.75 ,cent)
    

    n_big_blocks_each_quantile <- tab %>%
      group_by(gr = cut(distance_from_telomere, breaks = quantile(distance_from_telomere))) %>%
      summarise(all_blocks = n(), big_blocks = length(col[col == 'red'])) %>%
      mutate(big_blocks_percents = round(big_blocks/all_blocks*100, 1)) %>%
      filter(!is.na(gr))
    
    num_of_big_blocks <- data.frame(
      x_pos = max(tab$distance_from_telomere) * c(.125, .375, .625, .875),
      labels = n_big_blocks_each_quantile$big_blocks
    )
    
    num_of_big_blocks <- num_of_big_blocks %>%
      mutate(col = ifelse(labels > 0, 'grey', 'white'))
    
    par(
      las = 1, #labels always horizontal
      cex.axis = 1.6, 
      cex.lab = 1.8,
      cex.main = 2
    )
    
    plot(
      tab$distance_from_telomere, tab$length,
      xlab = 'Relative chromosome length, %',
      ylab = 'Length, Mb',
      main = paste0(sp, '/', unique(tab$el)),
      lwd = 0.5,
      col = tab$col,
      pch = tab$pch,
      xaxt="n",
      yaxt = 'n',
      ylim = c(200, 2100000),
      cex = 2
    )
    axis(1, at = xaxis_labels, labels = c('T', '25%', '50%', '75%', 'C'))
    axis(2, at = c(500000, 1000000, 1500000, 2000000), labels = c(0.5, 1, 1.5, 2))
    abline(v=xaxis_labels, col = 'grey', lty = 2)
    abline(h=quant_99, col = 'red', lty = 2)
    # text(x = max(tab$distance_from_telomere), y = quant_99-100000, col = 'red', labels = '99p level', pos = 2)
    # text(x = num_of_big_blocks$x_pos, y = 1300000, labels = num_of_big_blocks$labels)
    

  })
})




dev.off()


# Comparison with random model

generate_random_blocks <- function(n, genome_length){
  random_breaks <- floor(runif(n*2, min=0, max = genome_length))
  random_breaks <- random_breaks[order(random_breaks)]
  
  random_blocks <- lapply(1:(length(random_breaks)-1), function(m){
    random_breaks[m+1] - random_breaks[m]
  })
  
  random_rows_to_remove <- sample(1:n, replace = F)
  unlist(random_blocks[random_rows_to_remove])
}

## Generate 10 samples of random blocks

random_lengths <- lapply(1:1, function(n){
  generate_random_blocks(nrow(blocks_file), sum(blocks$atr$end))
})

names(random_lengths) <- paste0('random_', 1:length(random_lengths))

## Extract observed lengths from 3 species

observed_lengths <- lapply(names(blocks), function(sp){
  blocks[[sp]]$end
})

names(observed_lengths) <- names(blocks)
# Find 99 percentile for random and observed blocks length

find_percentiles <- function(list_of_lengths){
  unlist(lapply(list_of_lengths, function(lengths){
    quantile(lengths, probs = c(0.99))
  }))
}

random_99p_mean <- mean(find_percentiles(random_lengths))
observed_99p_mean <- mean(find_percentiles(observed_lengths))

std <- function(x) sd(x)/sqrt(length(x)) # se function

random_stat <- stack(bind_rows(random_lengths)) %>%
  group_by(ind) %>%
  summarise(p99 = quantile(values, probs=.99), mean_all = mean(values), mean_over99 = mean(values[values>p99]), se_over99 = std(values[values>p99]))

observed_stat <- stack(bind_rows(observed_lengths)) %>%
  group_by(ind) %>%
  summarise(p99 = quantile(values, probs = .99), mean_all = mean(values), mean_over99 = mean(values[values>p99]), se_over99 = std(values[values>p99]))

observed_lims_x_start <- seq(0.5, nrow(observed_stat)-0.5, 1)
observed_lims_x_end <- observed_lims_x_start + 1
random_lims_x_start <- seq(observed_lims_x_end[length(observed_lims_x_end)], nrow(random_stat)-0.5 + length(observed_lims_x_start), 1)
random_lims_x_end <- random_lims_x_start + 1

tiff('./output/compare_with_random_model2.tiff', width = 2500, height = 1500, units = 'px', res = 470, pointsize = 4)
ggplot() +
  geom_boxplot(data = stack(bind_cols(random_lengths)), aes(x = ind, y = values), fill = 'grey', colour = 'black') +
  geom_boxplot(data = stack(bind_cols(observed_lengths)), aes(x = ind, y = values), fill = 'orange', colour = 'black') +
  xlab('') +
  ylab('Length, Mb') +
  
  geom_segment(aes(x = random_lims_x_start, xend = random_lims_x_end, y = random_stat$p99, yend = random_stat$p99), col = 'grey', linetype = 'dashed') +
  geom_rect(aes(xmin = random_lims_x_start, xmax = random_lims_x_end, ymin = random_stat$mean_over99 - random_stat$se_over99, ymax = random_stat$mean_over99 + random_stat$se_over99), alpha = .5, fill = 'grey') +
  geom_segment(aes(x = random_lims_x_start, xend = random_lims_x_end, y = random_stat$mean_over99, yend = random_stat$mean_over99), col = 'grey') +

  geom_segment(aes(x = observed_lims_x_start, xend = observed_lims_x_end, y = observed_stat$p99, yend = observed_stat$p99), col = 'orange', linetype = 'dashed') +
  geom_rect(aes(xmin = observed_lims_x_start, xmax = observed_lims_x_end, ymin = observed_stat$mean_over99 - observed_stat$se_over99, ymax = observed_stat$mean_over99 + observed_stat$se_over99), alpha = .4, fill = 'orange') +
  geom_segment(aes(x = observed_lims_x_start, xend = observed_lims_x_end, y = observed_stat$mean_over99, yend = observed_stat$mean_over99), col = 'orange') +
  
  scale_x_discrete(labels = c('An. albimanus', 'An. atroparvus', 'An. gambiae', 'Random model')) +
  scale_y_continuous(breaks = c(500000, 1000000, 1500000, 2000000), labels = c(0.5, 1, 1.5, 2)) +
  annotate('text', x = c(observed_lims_x_end, random_lims_x_end), y = c(observed_stat$p99, random_stat$p99) - 20000, label = '99p level', colour = c(rep('orange', 3), 'grey'), size = 3, hjust = 1, vjust = 1)

dev.off()

wilcox_heatmap <- as.data.frame(bind_rows(lapply(random_lengths, function(random_length){
  bind_cols(lapply(observed_lengths, function(observed_length){
    wilcox.test(random_length, observed_length)$p.value < 0.01
  }))
})))

rownames(wilcox_heatmap) <- paste0('random_', 1:length(random_lengths))

# tiff('./output/compare_with_random_model.tiff', width = 1000, height = 1000, units = 'px', res = 200, pointsize = 4)
ggplot() +
  geom_histogram(aes(random_blocks, y=..count../sum(..count..)), fill = 'white', col = 'black', alpha = .5) +
  geom_vline(aes(xintercept = random_99)) +

  geom_histogram(aes(blocks$atr$end, y=..count../sum(..count..)), fill = 'red', col = 'red', alpha = .5, binwidth = 30000) +
  geom_vline(aes(xintercept = atr_99), col = 'red')+

  ylab('Density') +
  xlab('Blocks length')
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
