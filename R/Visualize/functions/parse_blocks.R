parse_blocks <- function(blocks){
  sp_count <- (ncol(blocks)-1) / 4
  blocks <- lapply(seq(2, sp_count*4, 4), function(start_col){
    end_col <- start_col+3
    blocks_set <- blocks[,start_col:end_col]
    names(blocks_set) <- c('chr', 'start', 'end', 'strand')
    return(blocks_set)
  })
}