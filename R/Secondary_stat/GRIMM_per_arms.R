library(gridExtra) 
library(grid)
library(dplyr)

source('./R/GRIMM/functions/make_grimm_command.R')

GRIMMtable.full <- read.table('./R/Clean/output_data/GRIMM.txt', header = TRUE) 
genome_size <- read.csv('./R/Secondary_stat/input_data/genomes_size.csv')[,-1]

elements <- paste0('e', c(1:5))
row.names(genome_size) <- elements
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

lengths <- lapply(elements, function(el){
  length <- list(
    alb = genome_size[el, 1],
    atr = genome_size[el, 2]
  )
})

names(lengths) <- elements

distance_tables_normalize <- lapply(elements, function(el){
  dist_table <- distance_tables[[el]]
  data.frame(
    atr_gam = round(dist_table[2, 3] / lengths[[el]]$atr, digits = 2),
    atr_alb = round(dist_table[1, 2] / lengths[[el]]$atr, digits = 2),
    alb_gam = round(dist_table[1, 3] / lengths[[el]]$alb, digits = 2)
  )
})

# Calculate mean normalize reverse distance 
RDperMB <- bind_cols(lapply(distance_tables_normalize, function(tab){
  tab$mean <- mean(unlist(tab[1, 1:3]))
  tab <- round(tab, digits = 2)
  as.data.frame(t(tab))
}))

RDperMBperMY <- bind_cols(lapply(distance_tables_normalize, function(tab){
  tab[1,1] <- tab[1,1] / 58
  tab[1,2] <- tab[1,2] / 100
  tab[1,3] <- tab[1,3] / 100
  tab$mean <- mean(unlist(tab[1, 1:3]))
  tab <- round(tab, digits = 3)
  as.data.frame(t(tab))
}))

# Visualize discance_tables 
distance_plot <- lapply(distance_tables, tableGrob) 
text_grobs <- lapply(paste0(elements, ', (mean=', mean_distance,')'), textGrob) 

grid.arrange(grobs=c(distance_plot, text_grobs), nrow=2, ncol=5, heights=unit(c(1,10), c("in", "mm"))) 
