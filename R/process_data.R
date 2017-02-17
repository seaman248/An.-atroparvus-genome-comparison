#source('./R/f_prepearing_data.R')
source('./R/f_prepearing_data_with_X.R')
source('./R/functions/arrange_chr.R')
source('./R/functions/through.R')
source('./R/functions/through_chr.R')
source('./R/functions/endToLength.R')
atr_order <- read.csv2('./data/atr_order_with_x.csv')

genes$atr_chr <- arrange_chr(genes[,c(6, 7)], atr_order[,c(1, 2)])
genes <- na.omit(genes)

genes[, c(6:10)] <- through_num(genes[,c(6:10)], atr_order)
genes <- na.omit(genes)
order <- as.character(unique(atr_order$chr))
rm(atr_order)

# make through order for gam and alb

genes[,c(13, 14)] <- through_chr(genes[,c(12, 13, 14)], order)
genes[,c(3, 4)] <- through_chr(genes[,c(2, 3, 4)], order)
rm(order)
genes$atr_chr <- NULL

# end coords to length

genes[,4] <- endToLength(genes[,c(3, 4)])
genes[,9] <- endToLength(genes[,c(8, 9)])
genes[,14] <- endToLength(genes[c(13, 14)])
genes <- plyr::rename(genes, c('alb_end'='alb_length', 'atr_end'='atr_length', 'gam_end'='gam_length'))

tableForGRIMM <- genes[, c(1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 15)]
write.table(tableForGRIMM,
            file='./data/processed/tableForGRIMM.txt',
            quote = FALSE,
            row.names = FALSE
            )


# data into grim

grimm_path <- '~/Documents/GRIMM_SYNTENY-2.02/grimm_synt'
input_path <- '~/rproj/full_gene/data/processed/tableForGRIMM.txt'
output_anchors <- '~/rproj/full_gene/data/processed/anchors'
output_blocks <- '~/rproj/full_gene/data/processed/'

anchor_command <- paste(grimm_path, '-A -f', input_path, '-d', output_anchors)

system(anchor_command, wait = TRUE)

min_block = '2'
min_gap= '1000'

blocks_command <- paste0(
  grimm_path, ' -f ', 
  output_anchors, '/unique_coords.txt', ' -d ', 
  output_blocks, 
  ' -c -n', min_block, ' -g ', min_gap, ' -Q'
  )

system(blocks_command, wait=TRUE)
