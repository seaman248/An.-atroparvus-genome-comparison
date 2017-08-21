library(ggplot2)
library(dplyr)
atr_blocks <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')[,c(6:9)]
atr_grimm <- read.table('./R/Clean/output_data/GRIMM.txt', stringsAsFactors = F)[, c(1, 6:9)]
orthologs <- read.csv2('./R/Query/output_data/orthologs.csv', stringsAsFactors = F)

alb_genes <- read.csv('./R/Query/output_data/albimanus_genes.csv')
atr_genes <- read.csv('./R/Query/output_data/atroparvus_genes.csv')
gam_genes <- read.csv('./R/Query/output_data/gambiae_genes.csv')
atr_grimm <- transform(atr_grimm, V7 = as.numeric(V7), V8 = as.numeric(V8), V9 = as.numeric(V9))

#mean(atr_blocks$V8) + sd(atr_blocks$V8) * 3

ggplot(g_blocks, aes(x = V6, y = V8))+
  geom_boxplot(outlier.size  = 0, outlier.alpha = 0)+
  stat_summary(fun.y = 'o', geom = 'point', col = 'red', size = 1.5)+
  theme_bw() +
  xlab('Arm') +
  ylab('Size of SB, bp')

atr_topmost_blocks <- g_blocks %>%
  group_by(V6) %>%
  filter(V8 >= o(V8))

##### Вычисляем ids An. albimanus

alb_ids <- apply(atr_topmost_blocks, 1, function(row){
  start_pos <- as.numeric(row[3])
  end_pos <- start_pos + as.numeric(row[4])
  
  g_table %>%
    mutate(end = start_position + end_position) %>%
    filter(start_position >= start_pos, end <= end_pos) %>%
    select(ensembl_gene_id)
  
})

names(alb_ids) <- atr_topmost_blocks$V1

UCBs <- atr_blocks %>%
  group_by(V6) %>%
  filter(V8 >= mean(V8) + 3 * sd(V8))

gam_ids <- lapply(alb_ids, function(ids){
  ids <- ids[, 1]
  q_orthologs[which(!is.na(match(q_orthologs$aalbimanus_eg_gene, ids))), 3]
})
