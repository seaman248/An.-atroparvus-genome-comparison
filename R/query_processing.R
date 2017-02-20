source('./R/functions/arrange_chr.R')
source('./R/functions/through.R')

orths.fulltable <- read.csv('./data/processed/alb_atr_gam_sin(with coords).csv')

# Arrange atroparvus
atr_order <- read.csv2('./data/atr_order_with_x.csv')
atr_chr <- arrange_chr(orths.fulltable[,c('X2.1')], atr_order[,c(1, 2)])
orths.fulltable[,c(6:10)] <- through_num(orths.fulltable[,c(6:10)], atr_order)

# Remove all rows with strange chromosomes in albimanus and gambiae

