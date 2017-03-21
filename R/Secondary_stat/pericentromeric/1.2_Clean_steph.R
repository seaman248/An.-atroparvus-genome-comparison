library(rjson)

# convert coords to indian coords
steph_convert <- fromJSON(file = 'R/Secondary_stat/pericentromeric/input_data/full_orth.json')
steph_convert <- bind_rows(steph_convert)

steph_genes <- read.csv2('R/Secondary_stat/pericentromeric/output_data/steph_peric_dirty.csv')

steph_convert_peric <- steph_convert[match(steph_genes$ensembl_gene_id, steph_convert$orth), c(4, 1, 5, 6)]

steph_genes[, c(2, 3, 4, 5)] <- steph_convert_peric

table(steph_genes$chromosome_name)

# make true order
steph_order <- read.csv2('R/Secondary_stat/pericentromeric/input_data/steph_order.csv')

steph_order <- na.omit(steph_order[steph_order$X != 'UNKN',])

quant_table <- table(steph_genes$chromosome_name)
quant_table <- as.data.frame(quant_table)

steph_scaffolds <- quant_table[match(steph_order$scf....unique.steph_genes.chromosome_name.,quant_table$Var1), ][, 1]

# atr_order
full_orth <- read.csv2('R/Secondary_stat/pericentromeric/output_data/peric_orthologs.csv')
atr_genes <- read.csv2('./R/Secondary_stat/pericentromeric/output_data/atroparvus_peric_clean.csv')

lapply(steph_scaffolds, function(steph_scf){
  steph_scf_genes <- na.omit(steph_genes[steph_genes$chromosome_name == steph_scf,])
  atr_scf_genes_coords <- match(steph_scf_genes$ensembl_gene_id, full_orth$steph)
  atr_scf_genes <- full_orth[atr_scf_genes_coords,]$atroparvus
  atr_scfs <- atr_genes[match(atr_scf_genes, atr_genes$ensembl_gene_id), ]$chromosome_name
  print(steph_scf)
  print(unique(atr_scfs))
})
