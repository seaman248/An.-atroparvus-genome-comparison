rowsWithBadChromosomes <- function(chrs, goodCrhs){
  which(!is.na(match(chrs, goodCrhs)))
}