library(genoPlotR)

source('./R/Visualize/functions/parse_blocks.R')
source('./R/Visualize/functions/create_seg.R')
source('./R/Visualize/functions/create_xlims.R')
source('R/Visualize/functions/make_comparison.R')

blocks <- parse_blocks(read.table('./R/GRIMM/output_data/blocks/blocks.txt'))


# length to end
blocks <- lapply(blocks, function(set){
  set$col[set$end > 400000] <- 'black'
  set$end <- set$start+set$end
  set
})

# create seq
seqs <- lapply(blocks, create_seqs)


# create xlims
arm_strand <- list(
  alb=c(F, T, F, F, T),
  atr=c(T, T, F, T, F),
  gam=c(F, T, F, T, F)
)
xlims <- as.list(data.frame(mapply(create_xlims, genes=blocks, strand=arm_strand)))
print(xlims)

# comparisons
comparisons <- make_comparisons(seqs)


# visualize
plot_gene_map(
  dna_segs=seqs,
  comparisons = comparisons,
  xlims = xlims
)

rm(list=ls())
