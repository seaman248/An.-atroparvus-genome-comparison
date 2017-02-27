### Query data
  ## Get table of orthologs sp1_gene_id / sp2_gene_id ... / ... spX_gene_id
  source('./R/Query/getOrthologs.R')
  
  ## Get coorinates for every gene of every species sp1: gene_id / chr / start / stop / strand
  source("./R/Query/getCoords.R")
  
  ## All data saved in ./R/Query/output_data
  
### Clean data
  ## Separate file per species
  lapply(c('albimanus', 'atroparvus', 'gambiae'), function(sp){
    files <- list.files('./R/Clean')
    file <- paste0('./R/Clean/' ,files[grep(sp, files)])
    source(file, echo=F)
    return(NULL)
  })
  
  ## Conmbine all sps in one table for GRIMM_synteny
  source('./R/Clean/GRIMM_format.R')

### Make GRIMM
  source('./R/GRIMM/make_blocks.R')

### Visualize Results