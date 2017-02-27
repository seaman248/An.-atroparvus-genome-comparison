changeSin <- function(sinS2, sinC2, bridgeCS){
  sinC2_geneIDs <- apply(sinS2, 1, function(generow){
    bridgeCS[match(generow[1], bridgeCS[,2]), 1]
  })
  
  sinC2_rows <- unlist(lapply(sinC2_geneIDs, function(sinC2_geneID){
    match(sinC2_geneID, sinC2[, 3])
  }))
  
  return(sinC2[sinC2_rows, c(3, 2, 4, 5, 6)])
}