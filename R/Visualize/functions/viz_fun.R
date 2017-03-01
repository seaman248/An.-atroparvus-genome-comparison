# Function for visualize blocks

# Blocks_list: (chr, start, length, strand)*sp

viz_fun <- function(path_to_blocks_table){
  library(genoPlotR)
  source('./R/Visualize/functions/create_seg.R')
  source('./R/Visualize/functions/create_xlims.R')
  source('./R/Visualize/functions/make_comparison.R')
  source('./R/Visualize/functions/parse_blocks.R')
  
  # Read from file
  blocks_list <- parse_blocks(read.table(path_to_blocks_table))
  
  # Length to strand
  blocks_list <- lapply(blocks_list, function(set){
    set$end <- set$start+set$end
    set
  })
  
  # Create segs
  seqs <- lapply(blocks_list, create_seqs)
  
  # Create xlims
  arm_strand <- list(
    alb=c(F, T, F, F, T),
    atr=c(T, T, T, T, T),
    gam=c(F, T, F, T, F)
  )
  xlims <- as.list(
    data.frame(
      mapply(create_xlims, genes=blocks_list, strand=arm_strand)
    )
  )
  
  # Create comparisons
  comparisons <- make_comparisons(seqs)
  
  plot <- plot_gene_map(
    dna_segs = seqs,
    xlims = xlims,
    comparisons = comparisons
  )
  return(plot)
}