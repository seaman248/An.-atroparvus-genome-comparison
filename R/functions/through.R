# genes: data.frame - ID, scf, start, end, strand
# order: data.frame - chr, scf, strand, size

through_num <- function (genes, order){
  # define incrementor
  IR <- 0
  
  order <- order[which(order[,2] %in% genes[,2]),]
  tGenes <- data.frame(matrix(NA, ncol = 5, nrow = nrow(genes)))
  names(tGenes) <- c('tID', 'tChr', 'tStart', 'tEnd', 'tStrand')
  tGenes$tID <- genes[,1]
  for(i in 1:nrow(order)){
    # find all cells with those scf
    scf_coords <- which(
      !is.na(
        match(genes[,2], order[i, 2])
      )
    )
    
    # if scf exist
    if(is.finite(max(genes[scf_coords,4]))){
      # if direct order
      if(order[i,3]==1){
        tGenes$tStart[scf_coords] <- genes[scf_coords,3] + IR
        tGenes$tEnd[scf_coords] <- genes[scf_coords,4] + IR
        tGenes$tChr[scf_coords] <- as.character(order[i, 1])
        tGenes$tStrand[scf_coords] <- genes[scf_coords,5]
      } # end if direct
      # if reverse order
      if(order[i, 3]==-1){
        tGenes$tStart[scf_coords] <- order[i, 4] - genes[scf_coords, 4] + IR
        tGenes$tEnd[scf_coords] <- order[i, 4] - genes[scf_coords, 3] + IR
        tGenes$tChr[scf_coords] <- as.character(order[i, 1])
        tGenes$tStrand[scf_coords] <- !tGenes$tStrand[scf_coords]
        
        tGenes$tStrand[scf_coords] <- genes[scf_coords,5]*(-1)
        
      } # end if reverse
      IR<- IR + order[i, 4]
    } # end if exist
  } # end for
  names(tGenes) <- names(genes)
  return(tGenes)
}
