library(ggplot2)
testGRIMM_g <- function(m='6000', g){
  grimm_path <- '~/Documents/GRIMM_SYNTENY-2.02/grimm_synt'
  input_path <- '~/rproj/full_gene/data/processed/tableForGRIMM.txt'
  output_anchors <- '~/rproj/full_gene/data/processed/test_anchors'
  output_blocks <- '~/rproj/full_gene/data/processed/test_blocks'
  
  anchor_command <- paste(grimm_path, '-A -f', input_path, '-d', output_anchors)
  system(anchor_command, wait = TRUE)
  
  min_block = m
  
  results <- list()
  pb <- txtProgressBar(min = 1, max = length(g), style = 3)
  for(i in 1:length(g)){
    Sys.sleep(0.1)
    gap=g[i]
    
    blocks_command <- paste0(
      grimm_path, ' -f ', 
      output_anchors, '/unique_coords.txt', ' -d ', 
      output_blocks, 
      ' -c -m', min_block, ' -g ', gap, ' -Q'
    )
    system(blocks_command, wait=TRUE)
    
    results[[i]] <- data.frame(read.table('./data/processed/test_blocks/blocks.txt'))
    setTxtProgressBar(pb, i)
  }
  return(results)
}

# experiment

gap_range <- seq(10, 250000, 5000)

grimm_g <- testGRIMM_g(g=gap_range)

# analyze

## Number of blocks

col_blocks <- unlist(lapply(grimm_g, nrow))

plot(x=gap_range, y=col_blocks, type='l', xlab = 'Gap size bp', ylab='Number of blocks')

## Mean length of blocks

mean_block_length <- lapply(lapply(grimm_g, '[', c(4, 8, 12)), function(lengths){
  return(c(
    mean(lengths[[1]]),
    mean(lengths[[2]]),
    mean(lengths[[3]])
  ))
})

mean_block_length <- do.call(rbind.data.frame, mean_block_length)
names(mean_block_length) <- c('alb', 'atr', 'gam')

sp_opt_gap <- lapply(c(1, 2, 3), function(n){
  gap_range[which(mean_block_length[,n]>142000)[1]]
})

plot(gap_range, mean_block_length$alb, type = 'l', col='green', xlab = 'Gap size, bp', ylab = 'Mean block length, bp')
lines(gap_range, mean_block_length$atr, col='orange')
lines(gap_range, mean_block_length$gam, col='red')
lines(c(0, 250000), c(142000, 142000), col='black', pch=22, lty=2)
lapply(sp_opt_gap, function(gap){
  lines(c(gap, gap), c(0, 142000), pch=22, lty=2)
})
legend(140000, 80000, legend=c('alb', 'atr', 'gam'), col=c('green', 'orange', 'red'), lty=1, box.lty = 0, cex=0.8)
axis(side = 1, at= seq(0, 250000, 5000),  labels = seq(0, 250000, 5000), padj = 1)

