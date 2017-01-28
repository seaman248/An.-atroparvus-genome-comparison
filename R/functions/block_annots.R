
blockAnnots <- function(dna_seg){
  coords <- middle(dna_seg)
  text <- dna_seg$name
  return(annotation(
    x1=coords, text = text, rot=90
  ))
}