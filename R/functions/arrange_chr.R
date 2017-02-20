#gtable ID, scf
#ctable chr, scf
arrange_chr <- function (gtable, ctable){
  chr_string <- rep(NA, length(gtable))
  
  for(i in 1:nrow(ctable)){
    chr_string[
      which(
        !is.na(
          match(gtable, ctable[i, 2])
        )
      )
      ] <- as.character(ctable[i,1])
  }
  
  return(chr_string)
}
