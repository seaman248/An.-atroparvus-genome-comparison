source('./R/functions/endToLength.R')

orths.table <- read.csv2('./data/alb_atr_gam_sin.csv', stringsAsFactors = FALSE)

orths.table[,seq(4, 19, 5)] <- lapply(seq(4, 19, 5), function(endcol){
  startcol <- endcol-1
  endToLength(orths.table[,c(startcol, endcol)])
})

write.table(orths.table[,c(1:5, 7:10, 12:15, 17:20)], './data/processed/tableForGRIMMquery.txt', quote=FALSE, row.names = FALSE)

grimm_path <- '~/Documents/GRIMM_SYNTENY-2.02/grimm_synt'
input_path <- '~/rproj/full_genome_comparison/data/processed/tableForGRIMMquery.txt'
output_anchors <- '~/rproj/full_genome_comparison/data/query_processed/anchors'
output_blocks <- '~/rproj/full_genome_comparison/data/query_processed/GRIMM'

anchor_command <- paste(grimm_path, '-A -f', input_path, '-d', output_anchors)

system(anchor_command, wait = TRUE)

min_block = '2'
min_gap= '115000'

blocks_command <- paste0(
  grimm_path, ' -f ', 
  output_anchors, '/unique_coords.txt', ' -d ', 
  output_blocks, 
  ' -c -n', min_block, ' -g ', min_gap, ' -Q'
)

system(blocks_command, wait=TRUE)

