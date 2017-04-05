source('./compare_plot/functions/index_f.R')
#source('./compare_plot/load_data.R')

# genes_to_compare <- genes.trueorder[genes.trueorder[, 2]==4 & genes.trueorder[, 7]==4 & genes.trueorder[, 12]==4, ]
# genes_to_compare <- genes_to_compare[order(genes_to_compare[,8]*-1),]
# genes_to_compare <- genes_to_compare[1:500,]
# 
# e4_orders <- list(
#   alb = 1,
#   atr = -1,
#   gam = -1
# )

genes_to_compare <- read.csv2('https://raw.githubusercontent.com/seaman248/scf_0043_location/master/input_data/un_table.csv')[, c(-1)]

par(mfrow=c(1, 5))
lapply(unique(genes_to_compare$chromosome_name.1), function(gam_chr){
  if(gam_chr != 'UNKN' & gam_chr != 'X'){
    local_gen_to_comp <- genes_to_compare[genes_to_compare$chromosome_name.1 == gam_chr, ]
    local_gen_to_comp <- local_gen_to_comp[order(local_gen_to_comp$start_position.1), ]
    
    local_tlist <- devide_tables(local_gen_to_comp)
    
    names(local_tlist) <- c('steph', 'gam')
    
    local_matches <- orth_matches(local_tlist, main_sp = 1)
    
    compare_plot(local_tlist, local_matches, main_sp = 1)
    gam_chr
  }
})

# tlist <- devide_tables(genes_to_compare)
# 
# names(tlist) <- c('steph', 'gam')
# 
# matches <- orth_matches(tlist, main_sp = 1)
# 
# compare_plot(tlist, matches, main_sp = 1)


