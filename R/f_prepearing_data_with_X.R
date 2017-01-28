library('plyr')
library('dplyr')

genes <- na.omit(
  read.csv('./data/genes.txt', na.strings = c('', NA))
)
Xgenes <- na.omit(
  read.csv('./data/Xgenes.txt', na.strings = c('', NA))
)
genes <- rbind(genes, Xgenes)
rm(Xgenes)
names(genes) <- c(
  'alb_ID',
  'alb_chr',
  'alb_start',
  'alb_end',
  'alb_strand',
  
  'atr_ID',
  'atr_scf',
  'atr_start',
  'atr_end',
  
  'gam_ID',
  'gam_chr',
  'gam_start',
  'gam_end'
)

genes <- genes %>%
  filter(gam_chr!='UNKN')

genes <- genes[!(duplicated(genes$alb_ID, fromLast = TRUE) | duplicated(genes$alb_ID)), ]
genes <- genes[!(duplicated(genes$atr_ID, fromLast = TRUE) | duplicated(genes$atr_ID)), ]
genes <- genes[!(duplicated(genes$gam_ID, fromLast = TRUE) | duplicated(genes$gam_ID)), ]

write.csv(genes$atr_ID, './data/processed/atrID.csv', row.names = FALSE, quote = FALSE)
write.csv(genes$gam_ID, './data/processed/gamID.csv', row.names = FALSE, quote = FALSE)

atr_strands <- read.csv('./data/atr_strands+X.txt')
gam_strands <- read.csv('./data/gam_strands+X.txt')

genes$atr_strand <- atr_strands$Strand[match(genes$atr_ID, atr_strands$Gene.stable.ID)]
genes$gam_strand <- gam_strands$Strand[match(genes$gam_ID, gam_strands$Gene.stable.ID)]
rm(atr_strands, gam_strands)

genes <- genes[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 10, 11, 12, 13, 15)]
