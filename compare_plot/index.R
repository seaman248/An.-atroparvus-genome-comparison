source('./compare_plot/functions/index_f.R')
source('./compare_plot/load_data.R')

genes_to_compare <- tail(genes.trueorder[genes.trueorder[,7]==4, ], 5000)
genes_to_compare <- genes_to_compare[genes_to_compare[, 2] != 4 & genes_to_compare[, 12] != 4, ]
compare_plot(genes_to_compare, 2, sps=c('alb', 'atr', 'gam'))


