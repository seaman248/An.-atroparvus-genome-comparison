library(biomaRt)
library(dplyr)
vb_host <- 'biomart.vectorbase.org'
gbase <- listMarts(host=vb_host)$biomart[1]

thisMart <- useMart(
  gbase,
  host=vb_host,
  dataset='aatroparvus_eg_gene'
)

table_886_938_scfs <- getBM(
  filters = c('chromosome_name'),
  values = c('KI421886', 'KI421938'),
  attributes = c('chromosome_name',
                 'ensembl_gene_id',
                 'chromosome_name',
                 'start_position',
                 'end_position',
                 'aalbimanus_eg_homolog_ensembl_gene',
                 'aalbimanus_eg_homolog_chromosome',
                 'aalbimanus_eg_homolog_chrom_start',
                 'aalbimanus_eg_homolog_chrom_end',
                 'agambiae_eg_homolog_ensembl_gene',
                 'agambiae_eg_homolog_chromosome',
                 'agambiae_eg_homolog_chrom_start',
                 'agambiae_eg_homolog_chrom_end'
  ),
  mart = thisMart
)

c_table_886_938_scfs <- na.omit(table_886_938_scfs)

c_table_886_938_scfs %>%
  arrange(chromosome_name, agambiae_eg_homolog_chrom_start) %>%
  slice((n()-5):n())
