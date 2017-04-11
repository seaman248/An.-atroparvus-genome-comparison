library(ggplot2)
library(dplyr)
atr_blocks <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')[,c(6:9)]

#mean(atr_blocks$V8) + sd(atr_blocks$V8) * 3

lengths <- as.numeric(cut(atr_blocks$V8, breaks = seq(0, 1700000, 100000))) / 10

atr_blocks$lengths <- lengths

plot(lengths)
axis(1, at = seq(0, max(lengths), 0.1))

ggplot(atr_blocks, aes(x = V6, y = lengths))+
  geom_boxplot()


atr_blocks %>%
  group_by(V6) %>%
  filter(V8 >= mean(V8) + 3 * sd(V8))



