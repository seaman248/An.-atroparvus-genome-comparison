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
    readLines(paste0('.', work_dir, 'blocks/report.txt'))
})

# Extract distance tables 
distance_tables <- lapply(GRIMM_reports, function(report){ 
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

# Calculate mean reverse distance 
mean_distance <- lapply(distance_tables, function(table){ 
  mean <- (table[1, 2]+table[1, 3]+table[2, 3])/3 
  round(mean, digits=1) 
}) 

# Visualize discance_tables 
distance_plot <- lapply(distance_tables, tableGrob) 
text_grobs <- lapply(paste0(elements, ', (mean=', mean_distance,')'), textGrob) 

grid.arrange(grobs=c(distance_plot, text_grobs), nrow=2, ncol=5, heights=unit(c(1,10), c("in", "mm"))) 
