library(genoPlotR)

source('./R/functions/create_seg.R')
source('./R/functions/create_xlims.R')
source('./R/functions/make_chr_annots.R')
source('./R/functions/block_annots.R')
source('./R/functions/make_comparison.R')

blocks <- read.table('./data/query_processed/GRIMM/blocks.txt')
names(blocks) <- c('id', rep(c('chr', 'start', 'end', 'strand'),4))
blocks[,1]<- as.character(blocks[,1])
blocks[,seq(4, 16, 4)] <-blocks[,seq(3, 15, 4)]+blocks[,seq(4, 16, 4)]

ablocks <- lapply(seq(5, 17, 4), function(n){
  cols <- c(1, (n-3):(n))
  blocks[,cols]
})

ablocks <- ablocks[c(1, 2, 4, 3)]

names(ablocks) <- c('alb', 'atr', 'sin', 'gam')

### Create dna_seqs
dna_segs <- lapply(ablocks, create_seqs)

### Create xlims
chr_order <- list(
  alb=c('X','2R', '3L', '2L', '3R'),
  atr=c('X', '3R', '2L', '2R', '3L'),
  sin=c('X','3R', '2L', '2R', '3L'),
  gam=c('X','2R', '2L', '3R', '3L')
)

cent_right <- list(
  alb=c(F, T, F, F, T),
  atr=c(T, T, T, T, T),
  sin=c(F, T, F, T, F),
  gam=c(F, T, F, T, F)
)

xlims <- as.list(
  data.frame(
    mapply(create_xlims, genes=ablocks, strand=cent_right, chrs=chr_order)
  )
)

### Create comparisons
comparisons <- make_comparisons(dna_segs)

### Create annotations
annotations <- as.list(data.frame(mapply(make_chr_annot, sp_xlims=xlims, chrs=chr_order)))


plot_gene_map(
  dna_segs = dna_segs,
  xlims=xlims,
  comparisons = comparisons
)
