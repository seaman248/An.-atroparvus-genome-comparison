atr2lblocks <- subset(ablocks$atr, end<=103070046 & start>=61701894)
atr2lblocks <- atr2lblocks[order(-atr2lblocks$start),]
atr2lblocks <- atr2lblocks[1:5,]

atr_2l_cgenes <- apply(atr2lblocks, 1, function(block){
  cgenes <- subset(genes, atr_start>=block[3]& atr_end<=block[4]&atr_scf=='2L')
  cgenes$block<- rep(block[1], nrow(cgenes))
  return(cgenes)
})

atr_2l_cgenes <- ldply(atr_2l_cgenes, data.frame)

write.csv2(atr_2l_cgenes, './data/processed/pericentromeric2Lblocks.csv')


atr2Rblocks <- subset(ablocks$atr, end<=61200434&start>=16121022)
atr2Rblocks <- atr2Rblocks[order(-atr2Rblocks$start),]
atr2Rblocks <- atr2Rblocks[1:10,]

atr_2R_cgenes <- apply(atr2Rblocks, 1, function(block){
  cgenes <- subset(genes, atr_start>=block[3]&atr_end<=block[4]&atr_scf=='2R')
  cgenes$block <- rep(block[1], nrow(cgenes))
  return(cgenes)
})

atr_2R_cgenes <- ldply(atr_2R_cgenes, data.frame)

write.csv2(atr_2R_cgenes, './data/processed/pericentromeric2Rblocks.csv')
