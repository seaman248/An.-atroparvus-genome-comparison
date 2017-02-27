source('R/Visualize/functions/genCol.R')

make_comparisons <- function(segs){
  comparisons <- list()
  for(i in 1:(length(segs)-1)){
    comparison_i <- comparison(data.frame(
      start1=segs[[i]][,2], end1=segs[[i]][,3],
      start2=segs[[(i+1)]][,2], end2=segs[[(i+1)]][,3]
    ))
    comparison_i$direction=segs[[i]][,4]*segs[[(i+1)]][,4]
    comparison_i$col <- unlist(lapply(comparison_i$direction, genCol))
    comparisons[[i]] <- comparison_i
  }
  return(comparisons)
}