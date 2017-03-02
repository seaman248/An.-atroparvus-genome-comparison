library(gridExtra) 
library(grid) 

source('./R/GRIMM/functions/make_grimm_command.R')

GRIMMtable.full <- read.table('./R/Clean/output_data/GRIMM.txt', header = TRUE) 

elements <- paste0('e', c(1:5)) 
# prepare table for GRIMM_synteny by chromosome 

GRIMM_tables.elX <- lapply(elements, function(elX){ 
  GRIMMtable.chr <- GRIMMtable.full[GRIMMtable.full[,6]==elX,] 
})

names(GRIMM_tables.elX) <- elements

GRIMM_reports <- lapply(elements, function(el){
  # set dirs
    work_dir <- '/R/Secondary_stat/_/'
    input_tab <- paste0(work_dir, el, '.txt')
    
  # write table
    write.table(GRIMM_tables.elX[[el]], paste0('.', input_tab), row.names = F)
    
  # make anchor
    anch_command <- make_grimm_command(input=input_tab, output = paste0(work_dir))
    system(anch_command, wait = T)
  # make grimm
    blocks_command <- make_grimm_command(input=paste0(work_dir, 'anchors/unique_coords.txt'), output=paste0(work_dir), A=F)
    system(blocks_command, wait = T)
  # read reports
    list(readLines(paste0('.', work_dir, 'blocks/report.txt')),
         read.table(paste0('.', work_dir, 'blocks/blocks.txt'))
    )
})
names(GRIMM_reports) <- elements

# Extract distance tables 
distance_tables <- lapply(GRIMM_reports, function(report){
  report <- report[[1]]
  distance_row <- grep('Distance Matrix', report) 
  rows_with_table <- (distance_row+1):(distance_row+3) 
  distance_table <- report[rows_with_table] 
  con <- textConnection(distance_table) 
  table <- read.table(con) 
  close(con) 
  names(table) <- c('alb', 'atr', 'gam') 
  rownames(table) <- c('alb', 'atr', 'gam') 
  return(table) 
}) 
names(distance_tables) <- elements

# Normalize distance tables

distance_tables_normalize <- lapply(elements, function(el){
  nblocks <- nrow(GRIMM_reports[[el]][[2]])
  dist_table <- distance_tables[[el]]
  as.data.frame(t(data.frame(
    atr_gam = round(dist_table[2, 3] / nblocks, digits = 2),
    atr_alb = round(dist_table[1, 2] / nblocks, digits = 2),
    alb_gam = round(dist_table[1, 3] / nblocks, digits = 2)
  )))
})

# Calculate mean reverse distance 
mean_distance <- lapply(distance_tables, function(table){ 
  mean <- (table[1, 2]+table[1, 3]+table[2, 3])/3 
  round(mean, digits=1) 
}) 

# Normalize mean
mean_dist_normalize <- lapply(elements, function(el){
  m <- mean_distance[[el]]/nrow(GRIMM_reports[[el]][[2]])
  round(m, digits = 2)
})

# Visualize discance_tables 
distance_plot <- lapply(distance_tables, tableGrob) 
text_grobs <- lapply(paste0(elements, ', (mean=', mean_distance,')'), textGrob) 

grid.arrange(grobs=c(distance_plot, text_grobs), nrow=2, ncol=5, heights=unit(c(1,10), c("in", "mm"))) 
