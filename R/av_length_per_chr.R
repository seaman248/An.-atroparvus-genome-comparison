library(dplyr)
library(ggplot2)
blocks <- read.table('./data/processed/blocks.txt')

names(blocks) <- c('id', rep(c('chr', 'start', 'end', 'strand'),3))
blocks[,1]<- as.character(blocks[,1])

ablocks <- lapply(seq(5, 13, 4), function(n){
  cols <- c(1, (n-3):(n))
  blocks[,cols]
})

mean_by_chr <- lapply(ablocks, function(blocks.table){
  as.list(with(blocks.table, tapply(blocks.table$end, blocks.table$chr, mean)))
})

names(mean_by_chr) <- c('alb', 'atr', 'gam')

chr_order <- list(
  alb = c(5, 2, 3, 1, 4),
  atr = c(5, 4, 1, 2, 3),
  gam = c(5, 2, 1, 4, 3)
)
  
mean_by_chr.order <- lapply(names(chr_order), function(sp){
  means <- rbind(mean_by_chr[[sp]])
  means <- means[,c(chr_order[[sp]])]
  means <- rbind(means)
  #means[1,] <- paste0(means[1,], c('e1', 'e2', 'e3', 'e4', 'e5'))
  colnames(means) <- c('e1', 'e2', 'e3', 'e4', 'e5')
  rownames(means) <- sp
  return(combine(means))
})

names(mean_by_chr.order) <- names(chr_order)

mean_by_chr.dataframe <- data.frame(do.call('rbind', mean_by_chr.order))

colnames(mean_by_chr.dataframe) <- c('e1', 'e2', 'e3', 'e4', 'e5')


mean_by_chr.means <- as.data.frame(colMeans(mean_by_chr.dataframe))
colnames(mean_by_chr.means) <- c('average block length')
