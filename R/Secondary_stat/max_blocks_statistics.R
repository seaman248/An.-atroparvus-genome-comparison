library(ggplot2)
library(dplyr)
g_table <- read.table('./R/Clean/output_data/GRIMM.txt', header = T)
g_blocks <- read.table('./R/GRIMM/output_data/blocks/blocks.txt')
q_orthologs <- read.csv2('./R/Query/output_data/orthologs.csv', stringsAsFactors = F)

o <- function(size){
  size[size >= quantile(size, probs = c(0.98))[1]]
}

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

##### Вычисляем ids An. gambiae

gam_ids <- lapply(alb_ids, function(ids){
  ids <- ids[, 1]
  q_orthologs[which(!is.na(match(q_orthologs$aalbimanus_eg_gene, ids))), 3]
})

gam_ids <- bind_rows(lapply(names(gam_ids), function(n){
  data.frame(
    block = n,
    ids = gam_ids[[n]]
  )
}))

block_stat <- gam_ids %>%
  group_by(block) %>%
  summarize(n())

# agam_expression <- read.table('./R/Secondary_stat/Anopheles-gambiae_MSM-STATS_VB-2017-04.txt', sep = '\t', header = T)

exp_blocks <- bind_cols(apply(agam_expression[, 2:5], 2, function(col){
  exp_data <- data.frame(agam_expression[, 1], col) %>%
    filter(!is.na(col), col<= 0.05)
  exp_ids <- as.character(exp_data[, 1])
  
  matches <- which(!is.na(match(
    exp_ids, gam_ids$ids
  )))
  
  exp_stat <- gam_ids %>%
    slice(matches) %>%
    group_by(block) %>%
    summarize(n())
  
  matches <- match(
    block_stat$block, exp_stat$block
  )
  
  stat <- data.frame(block_stat, exp_stat[matches, 2])
  
  colnames(stat) <- c('block', 'exp', 'obs')
  
  n_exp <- stat %>%
    mutate(n_exp = obs / exp * 100) %>%
    select(n_exp)
  
  if(any(na.omit(n_exp) > 90)){
    which(!is.na(n_exp[n_exp > 90]))
  }
}))


