# genes: chr, start, end
seq_num <- function(genes, chrs=as.character(unique(genes[,1]))){
  prev_max_coord <- max(genes[genes[,1]==chrs[1],3])
  
  for(i in 2:length(chrs)){
    genes[genes[,1]==chrs[i],c(2, 3)] <-
      genes[genes[,1]==chrs[i],c(2, 3)] + prev_max_coord
    
    prev_max_coord <- max(genes[genes[,1]==chrs[i],3])
  }
  return(genes[,c(2,3)])
}