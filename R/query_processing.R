source('./R/functions/rowsWithBadChromosomes.R')
source('./R/functions/arrange_chr.R')

orths.fulltable <- read.csv('./data/processed/alb_atr_gam_sin(with coords).csv')

# Arrange atroparvus
atr_order <- read.csv2('./data/atr_order_with_x.csv')
atr_chr <- arrange_chr(orths.fulltable[,c('X2.1')], atr_order[,c(1, 2)])


# Remove all rows with strange chromosomes in albimanus and gambiae
goodChrs <- c('X', '2R', '2L', '3R', '3L')
orths.goodChrs <- orths.fulltable[rowsWithBadChromosomes(orths.fulltable$X2, goodChrs),]
orths.goodChrs <- orths.goodChrs[rowsWithBadChromosomes(orths.goodChrs[,'X2.2'], goodChrs),]


