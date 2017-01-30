genesAsOrth <- function(genes, orths){
  matchOrth <- match(orths, genes[,1])
  newGenes <- genes[matchOrth,]
  return(newGenes)
}
