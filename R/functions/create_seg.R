source('./R/functions/genCol.R')
#genes(data.frame): id, start, end, strand
create_seqs <- function(genes){
  return(dna_seg(data.frame(
    name=genes[,1],
    start=genes[,3],
    end=genes[,4],
    strand=genes[,5],
    gene_type='blocks',
    col=unlist(lapply(genes[,5], genCol))
  )))
}