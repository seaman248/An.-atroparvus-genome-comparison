### Query data
  ## Get table of orthologs sp1_gene_id / sp2_gene_id ... / ... spX_gene_id
  source('./R/Query/getOrthologs.R')
  
  ## Get coorinates for every gene of every species sp1: gene_id / chr / start / stop / strand
  source("./R/Query/getCoords.R")
  
  ## All data saved in ./R/Query/output_data