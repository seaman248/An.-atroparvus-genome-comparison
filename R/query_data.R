source('./R/functions/getGenesFromMart.R')
source('./R/functions/genesAsaOrth.R')
source('./R/functions/through.R')


# define species of interest
sp_of_interest <- c(
  'albimanus', 'atroparvus', 'gambiae', 'sinensis'
)

# get table of orthologs from biomart
orths <- as.list(getOrthologs(sp_of_interest))

# get all coordinates for genes and regroup that
genes <- mapply(getGenesFromMart, geneIDs = orths, as.list(TRUE), dsName=names(orths))
genes <- as.list(as.data.frame(genes))
genes <-lapply(genes, as.data.frame)

# reorder genes accoring to ortologs order
wgenes <- as.list(as.data.frame(mapply(genesAsOrth, genes=genes, orths=orths)))
wgenes <- lapply(wgenes, as.data.frame)
# scf -> chr in atroparvus and stephensi
atr_order <- read.csv2('./data/atr_order_with_x.csv')

wgenes$aatroparvus_eg_gene <- through_num(wgenes$aatroparvus_eg_gene, atr_order)
