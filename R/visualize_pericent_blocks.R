#atr_2R_cgenes <- atr_2R_cgenes[,c(1:16, 18)]
peric_genes <- lapply(c('alb', 'atr', 'gam'), function(sp){
  spframe <- atr_2R_cgenes[,c(grep(sp, names(atr_2R_cgenes)), 17)]
  spframe <- rbind(spframe, atr_2l_cgenes[,c(grep(sp, names(atr_2l_cgenes)), 17)])
  names(spframe) <- c('name', 'chr', 'start', 'end', 'strand', 'block')
  return(spframe)
})

block_palete <- data.frame(
  block_id=unique(peric_genes[[2]][,6]),
  color=rainbow_hcl(length(unique(peric_genes[[2]][,6])))
)

peric_pallete <- lapply(peric_genes, function(genes){
  block_palete[match(genes[,6], block_palete[,1]),2]
})

peric_segs <- lapply(peric_genes, create_seqs)
peric_segs <- as.list(data.frame(mapply(function(seg, col){
  seg$col <- col
  return(dna_seg(seg))
}, seg=peric_segs, col=peric_pallete)))
peric_segs <- lapply(peric_segs, dna_seg)
names(peric_segs) <- c('albimanus', 'atroparvus', 'gambiae')

rm(block_palete, peric_pallete)

chr_order_peric_comparison <- list(
  alb=c('X','2R', '3L', '2L', '3R'),
  atr=c('2L','2R'),
  gam=c('2R', '2L', '3R', '3L')
)

cent_right_peric_comparison <- list(
  alb=c(F, T, F, F, T),
  atr=c(T, T),
  gam=c(T, F, T, F)
)

peric_xlims <- mapply(create_xlims, genes=peric_genes, strand=cent_right_peric_comparison, chrs=chr_order_peric_comparison)


peric_comparisons <- make_comparisons(peric_segs)

tiff('./output/peric_reg.tiff', width=8000, height=5000, units='px', compression = 'none', res=600)
plot_gene_map(
  dna_segs = peric_segs,
  xlims = peric_xlims,
  comparisons = peric_comparisons
)
dev.off()

pdf('./output/peric_reg.pdf', width=50, height=20)
plot_gene_map(
  dna_segs = peric_segs,
  xlims=peric_xlims,
  comparisons = peric_comparisons
)
dev.off()