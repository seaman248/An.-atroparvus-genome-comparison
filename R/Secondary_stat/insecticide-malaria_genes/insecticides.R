library(biomaRt)
library(dplyr)
source('R/Visualize/functions/parse_blocks.R')

vb_host <- 'biomart.vectorbase.org' 
mart <- useMart(listMarts(host=vb_host)$biomart[1], host = vb_host, dataset = 'agambiae_eg_gene') 

genes <- read.table('./R/Clean/output_data/GRIMM.txt', stringsAsFactors = F)
blocks <- parse_blocks(read.table('./R/GRIMM/output_data/blocks/blocks.txt'))
orthologs <- read.csv2('./R/Query/output_data/orthologs.csv', stringsAsFactors = F)


pfam_ids <- getBM(
  attributes = c('ensembl_gene_id', 'pfam_pf'),
  mart = mart
)

gam_regions <- getBM(
  attributes = c('ensembl_gene_id', 'band'),
  mart = mart
)

MIN_BIG_BLOCK <- 500000

atr_blocks <- blocks[[2]]
atr_big_blocks <- atr_blocks[atr_blocks$end > MIN_BIG_BLOCK, ]


gam_ids <- apply(atr_big_blocks, 1, function(row){
  end <- as.numeric(row[2]) + as.numeric(row[3])
  alb_ids <- genes[genes$V7 == anchor_coords$V7, 1]

  res <- data.frame(
    ids = orthologs[which(!is.na(match(orthologs[, 1], alb_ids))), 3]
  )
  res$band <- gam_regions[match(res$ids, gam_regions[, 1]), 2]
  res$pfam <- pfam_ids[match(res$ids, pfam_ids[, 1]), 2]
  res
})

lapply(gam_ids, function(df){
  df <- df %>% slice(which(!is.na(match(df$pfam, IRP)))) %>%
    select(pfam, band) %>%
    group_by(pfam, band)
})

# insecticide-resistent protein ids:

IRP <- c('PF00756', 'PF13417', 'PF00067', 'PF06512')

expected_freq <- length(which(!is.na(match(pfam_ids[, 2], IRP)))) / length(pfam_ids[, 2])
IRP_blocks <- unlist(lapply(gam_ids, function(df){
  actual_freq <- length(which(!is.na(match(df$pfam, IRP)))) / nrow(df)
  actual_freq / expected_freq
}))


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

unlist(lapply(gam_ids[IRP_blocks > 5], function(IRP_ids){
  getmode(gam_regions[match(IRP_ids, gam_regions[, 1]), 2])
}))


