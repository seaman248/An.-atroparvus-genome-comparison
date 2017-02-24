library(plyr)
source('./R/Clean/functions/rmBadChr.R')
source('./R/Clean/functions/seq_num.R')

gam_genes <- read.csv('./R/Query/output_data/gambiae_genes.csv')


# remove rows with bad chromosome names

gam_genes <- gam_genes[rmBadChr(gam_genes$chromosome_name),]

# rename chromosome arms from X, 2R... to e1, e2...

gam_genes$chromosome_name <- revalue(gam_genes$chromosome_name, c(
  'X' = 'e1',
  '2R' = 'e2',
  '2L' = 'e3',
  '3R' = 'e4',
  '3L' = 'e5'
))

# make sequential numbering of gene coordinates from e1 to e5 (e1:1-100, e2:100-200 ...)

gam_genes[,c(3, 4)] <- seq_num(gam_genes[,c(2, 3, 4)], chrs = c(paste0('e', c(1:5))))

write.csv2(gam_genes, './R/Clean/output_data/gam_clean.csv', row.names = F)

rm(list=ls())

