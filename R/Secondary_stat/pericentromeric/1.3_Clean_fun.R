# load orths, atr, fun
full_orth <- read.csv2('./R/Secondary_stat/pericentromeric/output_data/peric_orthologs.csv')
atr_genes <- read.csv2('./R/Secondary_stat/pericentromeric/output_data/atroparvus_peric_clean.csv')
fun_genes <- read.csv2('./R/Secondary_stat/pericentromeric/output_data/funestus_peric_dirty.csv')

fun_order <- read.table('./R/Secondary_stat/pericentromeric/input_data/fun_order.csv')[, 2:4]
