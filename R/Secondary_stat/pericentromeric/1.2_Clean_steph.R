library(rjson)

# convert coords to indian coords
steph_convert <- fromJSON(file = 'R/Secondary_stat/pericentromeric/input_data/full_orth.json')
steph_convert <- bind_rows(steph_convert)

steph_genes <- read.csv2('R/Secondary_stat/pericentromeric/output_data/steph_peric_dirty.csv')

steph_convert_peric <- steph_convert[match(steph_genes$ensembl_gene_id, steph_convert$orth), c(4, 1, 5, 6)]

steph_genes[, c(2, 3, 4, 5)] <- steph_convert_peric


# make true order
steph_order <- read.csv2('R/Secondary_stat/pericentromeric/input_data/steph_order.csv')

