library('genoPlotR')
source('./R/process_data.R')
source('./R/functions/create_seg.R')
source('./R/functions/create_xlims.R')
source('./R/functions/make_chr_annots.R')
source('./R/functions/block_annots.R')
source('./R/functions/make_comparison.R')
blocks <- read.table('./data/processed/blocks.txt')
names(blocks) <- c('id', 
                   'chr', 'start', 'end', 'strand',
                   'chr', 'start', 'end', 'strand',
                   'chr', 'start', 'end', 'strand'
                   )


blocks[,1]<- as.character(blocks[,1])
# length to end coordinates
blocks[,c(4, 8, 12)] <- blocks[,c(3, 7, 11)]+blocks[,c(4, 8, 12)]

alb_blocks <- blocks[,c(1, 2, 3, 4, 5)]
atr_blocks <- blocks[,c(1, 6, 7, 8, 9)]
gam_blocks <- blocks[,c(1, 10, 11, 12, 13)]

ablocks <- list(alb_blocks, atr_blocks, gam_blocks)
names(ablocks) <- c('alb', 'atr', 'gam')
rm(blocks, alb_blocks, atr_blocks, gam_blocks)

# make atroparvus main
ablocks <- lapply(ablocks, function(blocks){
  blocks[,5] <- blocks[,5]*ablocks[[2]][5]
  return(blocks)
})

# create segs
dna_segs <- lapply(ablocks, create_seqs)

# create xlims

chr_order <- list(
  alb=c('X','2R', '3L', '2L', '3R'),
  atr=c('X', '3R', '2L', '2R', '3L'),
  gam=c('X','2R', '2L', '3R', '3L')
)

cent_right <- list(
  alb=c(F, T, F, F, T),
  atr=c(T, T, T, T, T),
  gam=c(F, T, F, T, F)
)


xlims <- as.list(
  data.frame(
    mapply(create_xlims, genes=ablocks, strand=cent_right, chrs=chr_order)
  )
)

# create comparisons
comparisons <- make_comparisons(dna_segs)

#create annotations

alb_order <- c('X','2L(e4)', '2R(e2)', '3L(e3)','3R(e5)')
atr_order <- c('X','2L(e3)', '2R(e4)', '3L(e5)', '3R(e2)')
gam_order <- c('X','c2L(e3)', '2R(e2)', '3L(e5)', '3R(e4)')

alb_annot <- make_chr_annot(xlims[[1]], alb_order)
atr_annot <- make_chr_annot(xlims[[2]], atr_order)
gam_annot <- make_chr_annot(xlims[[3]], gam_order)
# alb_annot <- blockAnnots(alb_seg)
# atr_annot <- blockAnnots(atr_seg)
# gam_annot <- blockAnnots(gam_seg)

annots <- list(alb_annot, atr_annot, gam_annot)
rm(alb_order, atr_order, gam_order, alb_annot, atr_annot, gam_annot)
# create phyl tree
tree <- newick2phylog("((An_atroparvus,An_gambiae), An_albimanus);")

#tiff('./output/plot_c_p_m5_g100(main atr).tiff', width=15000, height=2500, units='px', compression = 'none', res=600)
plot_gene_map(
  dna_segs = dna_segs,
  xlims = xlims,
  comparisons = comparisons,
  dna_seg_labels = NULL,
  dna_seg_label_cex = 0
)
#dev.off()