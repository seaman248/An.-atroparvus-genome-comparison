source('./R/functions/arrange_chr.R')
source('./R/functions/through.R')
source('./R/functions/through_chr.R')
orths.fulltable <- read.csv('./data/processed/alb_atr_gam_sin(with coords).csv')

# Arrange atroparvus
atr_order <- read.csv2('./data/atr_order_with_x.csv')
atr_chr <- arrange_chr(orths.fulltable[,c('X2.1')], atr_order[,c(1, 2)])
orths.fulltable[,c(6:10)] <- through_num(orths.fulltable[,c(6:10)], atr_order)

# Arrange albimanus and gambiae
  
  goodChrs <- c("X",  "2R", "2L", "3R", "3L")
  
  # albimanus
  orths.fulltable[,c(3, 4)] <- through_chr(orths.fulltable[,c(2:4)], chrs=goodChrs)

  # gambiae
  orths.fulltable[,c(13, 14)] <- through_chr(orths.fulltable[,c(12:14)], chrs=goodChrs)
  
  # make NA rows instead bad chrs in albimanus and gambiae
  orths.fulltable[which(is.na(match(orths.fulltable[,c(2)], goodChrs))),] <- NA
  orths.fulltable[which(is.na(match(orths.fulltable[,c(12)], goodChrs))),] <- NA

  # Work with Sinensis
  SinC.scaffolds <- read.csv2('https://raw.githubusercontent.com/seaman248/sinensis/master/data/SinC_Scaffolds.csv')[,c(4, 1, 5, 2, 3)]
  SinC.scaffolds <- SinC.scaffolds[c(1:49, 51),c(1, 2, 3, 5)]
  SinC.SinS <- read.csv2('https://github.com/seaman248/sinensis/raw/master/data/S2_to_C2_bridge.csv')[,c(3, 2)]
  SinC.fulltable <- read.csv('https://raw.githubusercontent.com/seaman248/sinensis/master/data/sinC_genes.csv')
  # Replace SinS2 gene_id to SinC2
  # Add SinC2 coordinates (Scf, start, end, strand)
  # Make through num

# Remove all rows with strange chromosomes in albimanus and gambiae

