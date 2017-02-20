source('./R/functions/getGenesFromMart.R')
source('./R/functions/genesAsaOrth.R')
source('./R/functions/through.R')


# define species of interest
sp_of_interest <- c(
  'albimanus', 'atroparvus', 'gambiae', 'sinensis'
)

# get table of orthologs from biomart
orths <- as.list(getOrthologs(sp_of_interest))

# make 
orths <- data.frame(mapply(function(or){
  or
}, or=orths))

write.csv2(orths, './data/processed/orthsList.csv')
