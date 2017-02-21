source('./R/functions/arrange_chr.R')
source('./R/functions/through.R')
source('./R/functions/through_chr.R')
source('./R/functions/changeSinCoords.R')
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
  
  #
  
  orths.fulltable[,c(16:20)] <- changeSin(orths.fulltable[,c(16:20)], SinC.fulltable, SinC.SinS)
  # Make through num
  orths.fulltable[,20] <- unlist(lapply(orths.fulltable[,20], function(n){
    if(is.na(n)){
      return(NA)
    }else{
      n <- paste0(n, '1')
      as.numeric(n)
    }
  }))
  orths.fulltable[,c(16:20)] <- through_num(orths.fulltable[,c(16:20)], SinC.scaffolds)

# Remove all NAs
  
  orths.naomit <- na.omit(orths.fulltable)
  
# Make names
  names <- c('id', 'chr', 'start', 'end', 'strand')

  names(orths.naomit) <- unlist(lapply(c('alb', 'atr', 'gam', 'sin'), function(sp){
    paste0(sp, '_',names)
  }))
  
# write csv
  
  write.csv2(orths.naomit, './data/alb_atr_gam_sin.csv', row.names = FALSE)
  