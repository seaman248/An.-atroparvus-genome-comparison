make_chr_annot <- function(sp_xlims, chrs){
  return(annotation(
    x1=sp_xlims[seq(1, length(sp_xlims), 2)],
    x2=sp_xlims[seq(2, length(sp_xlims), 2)],
    text=as.character(chrs)
  ))
}