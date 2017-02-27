#genes(data.frame): chr, start, end
#chrs: examp: c('2L', '2R')
create_xlims <- function(genes, chrs=paste0('e', 1:5), strand){
  xlims <- c()
  for(i in 1:length(chrs)){
    min_start <- min(genes[genes[,1]==chrs[i], 2])
    max_end <- max(genes[genes[,1]==chrs[i], 3])
    if(strand[i]){
      xlims <- append(xlims, c(max_end, min_start))
    } else {
      xlims <- append(xlims, c(min_start, max_end))
    }
    
  }
  return(xlims)
}
