source('./R/Visualize/functions/genCol.R')
#genes(data.frame): id, start, end, strand
create_seqs <- function(genes){
  return(dna_seg(data.frame(
    name=row.names(genes),
    start=genes[,2],
    end=genes[,3],
    strand=genes[,4],
    gene_type='side_blocks'
    #, col=genes$col
  )))
}