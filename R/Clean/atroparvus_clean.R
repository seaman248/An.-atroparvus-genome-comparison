library(dplyr)
library(plyr)
source('./R/Clean/functions/scf_to_chr.R')
source('./R/Clean/functions/rmBadChr.R')

atr_genes <- read.csv('./R/Query/output_data/atroparvus_genes.csv')

# translate coordinates of scaffold to coordinate of chromosomes
atr_order <- read.csv2('./R/Clean/input_data/atr_order.csv')

for(chr in c('X', '2L', '3L')){
  m <- which(atr_order$chr == chr)
  atr_order[m, ] <- atr_order[rev(m), ]
}

atr_genes <- na.omit(scf_to_chr(atr_genes, atr_order))

rm(atr_order)
# remove rows with bad chromosome names

atr_genes <- atr_genes[rmBadChr(atr_genes$chromosome_name),]

# rename X, 2R... to e1, e2...

atr_genes$chromosome_name <- revalue(atr_genes$chromosome_name, c(
  'X' = 'e1',
  '3R' = 'e2',
  '2L' = 'e3',
  '2R' = 'e4',
  '3L' = 'e5'
))

# save to file

write.csv2(atr_genes, './R/Clean/output_data/atr_clean.csv', row.names = F)

rm(list=ls())
