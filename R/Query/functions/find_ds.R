library(biomaRt)

source('./R/Query/functions/config_ds.R')

# Find true dataset name in vectorbase
find_ds <- function(sp_string){
  sp_ds <- as.character(
    listDatasets(
      useMart(
        gbase, 
        host=vb_host
      )
    )$dataset
  )
  return(
    sp_ds[grep(sp_string, sp_ds)]
  )
}