library(biomaRt)
vb_host <- 'biomart.vectorbase.org'
gbase <- listMarts(host=vb_host)$biomart[1]

sp_of_interest <- c(
  'albimanus','atroparvus', 'gambiae'
)