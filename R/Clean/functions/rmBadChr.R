rmBadChr <- function(chrs, goodChrs=c('X', '2R', '2L', '3R', '3L')){
  which(!is.na(match(chrs, goodChrs)))
}