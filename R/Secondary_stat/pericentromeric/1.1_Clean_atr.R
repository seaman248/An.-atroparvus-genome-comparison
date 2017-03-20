library(plyr)
library(dplyr)
atr_genes <- read.csv2('R/Secondary_stat/pericentromeric/output_data/atroparvus_peric_dirty.csv')

atr_order <- data.frame(
  scf = c('KI421900', 'KI421924', 'KI421922', 'KI421894'),
  chr = c('e4', 'e4', 'e4', 'e3'),
  strand = c(-1, 1, -1, 1),
  length = c(3534686, 300483, 450785, 5662806),
  orientation = c('R', 'R', 'R', 'L')
)

atr_genes <- lapply(1:nrow(atr_order), function(order_n){
  sc_genes <- atr_genes[atr_genes$chromosome_name == atr_order[order_n, 1],]
  if(atr_order[order_n, 3] == -1){
    sc_genes$end_position <- atr_order[order_n, 4] - sc_genes$start_position
    sc_genes$start_position <- atr_order[order_n, 4] - sc_genes$end_position
    sc_genes$strand <- sc_genes$strand * -1
  }
  if(order_n > 1){
    sc_genes$start_position <- sc_genes$start_position + atr_order[order_n-1, 4]
    sc_genes$end_position <- sc_genes$end_position + atr_order[order_n-1, 4]
  }
  sc_genes$orientation <- atr_order[order_n,5]
  return(sc_genes)
})

write.csv2(bind_rows(atr_genes), file='./R/Secondary_stat/pericentromeric/output_data/atroparvus_peric_clean.csv')
