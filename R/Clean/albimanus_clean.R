library(plyr)
source('./R/Clean/functions/rmBadChr.R')
source('./R/Clean/functions/seq_num.R')

alb_genes <- read.csv('./R/Query/output_data/albimanus_genes.csv')

# remove rows with bad chromosome names

alb_genes <- alb_genes[rmBadChr(alb_genes$chromosome_name),]

# rename chromosome arms from X, 2R... to e1, e2...

alb_genes$chromosome_name <- revalue(alb_genes$chromosome_name, c(
  'X' = 'e1',
  '2R' = 'e2',
  '3L' = 'e3',
  '2L' = 'e4',
  '3R' = 'e5'
))

# make sequential numbering of gene coordinates from e1 to e5 (e1:1-100, e2:100-200 ...)

alb_genes[,c(3, 4)] <- seq_num(alb_genes[,c(2, 3, 4)], chrs = c(paste0('e', c(1:5))))

write.csv2(alb_genes, './R/Clean/output_data/alb_clean.csv', row.names = F)
