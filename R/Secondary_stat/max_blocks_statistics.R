library(ggplot2)
library(dplyr)
atr_blocks <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')[,c(6:9)]
atr_grimm <- read.table('./R/Clean/output_data/GRIMM.txt', stringsAsFactors = F)[, c(1, 6:9)]
orthologs <- read.csv2('./R/Query/output_data/orthologs.csv', stringsAsFactors = F)

o <- function(size){
  size[size >= quantile(size, probs = c(0.98))[1]]
}

ggplot(atr_blocks, aes(x = V6, y = V8))+
  geom_boxplot(outlier.size  = 0, outlier.alpha = 0)+
  stat_summary(fun.y = 'o', geom = 'point', col = 'red', size = 1.5)+
  theme_bw() +
  xlab('Arm') +
  ylab('Size of SB, bp')

topmost_blocks <- atr_blocks %>%
  group_by(V6) %>%
  filter(V8 >= o(V8))

