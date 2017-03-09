source('./compare_plot/functions/index_f.R')
source('./compare_plot/load_data.R')

genes_to_compare <- genes.trueorder[genes.trueorder[, 2]==4 & genes.trueorder[, 7]==4 & genes.trueorder[, 12]==4, ]
genes_to_compare <- genes_to_compare[order(genes_to_compare[,8]*-1),]
genes_to_compare <- genes_to_compare[1:500,]

e4_orders <- list(
  alb = 1,
  atr = -1,
  gam = -1
)

tlist <- devide_tables(genes_to_compare)

names(tlist) <- c('alb', 'atr', 'gam')

matches <- orth_matches(tlist, main_sp = 2, order = e4_orders)

compare_plot(tlist[c(2, 3)], matches[2], main_sp = 1)
write.csv2(bind_cols(tlist[c(2, 3)]), './compare_plot/output_data/test_output/atr_gam.csv')
