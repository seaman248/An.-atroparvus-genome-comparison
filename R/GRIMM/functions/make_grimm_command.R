make_grimm_command <- function(input='/R/Clean/output_data/GRIMM.txt',
                               output='/R/GRIMM/output_data',
                               A=T, min_block = 2, gap = '115000',
                                grimm_path = '~/Documents/GRIMM_SYNTENY-2.02/grimm_synt'){
  input <- paste0(getwd(), input)
  output <- paste0(getwd(), output)
  
  base_command <- paste0(grimm_path, ' -f ', input)
  
  if(A){
    anchors_dir <- paste0(output, '/anchors')
    if(!dir.exists(anchors_dir)) dir.create(anchors_dir)
    
    anchor_command <- paste0(base_command, ' -d ', anchors_dir, ' -A')
    return(anchor_command)
  } else {
    blocks_dir <- paste0(output, '/blocks')
    if(!dir.exists(blocks_dir)) dir.create(blocks_dir)
    
    blocks_command <- paste0(base_command, ' -d ', blocks_dir, ' -c ', '-n ', min_block, ' -g ', gap, ' -Q')
    return(blocks_command)
  }
}

