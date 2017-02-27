library(gridExtra)
library(grid)
GRIMMtable.full <- read.table('./data/processed/tableForGRIMM.txt', header = TRUE)

# Rename chromosomes

GRIMMtable.full$alb_chr <- revalue(GRIMMtable.full$alb_chr, c(
  'X'= 'e1',
  '2R'= 'e2',
  '2L'='e4',
  '3R'='e5',
  '3L'='e3'
))

GRIMMtable.full$atr_scf <- revalue(GRIMMtable.full$atr_scf, c(
  'X'= 'e1',
  '2R'= 'e4',
  '2L'='e3',
  '3R'='e2',
  '3L'='e5'
))

GRIMMtable.full$gam_chr <- revalue(GRIMMtable.full$gam_chr, c(
  'X'= 'e1',
  '2R'= 'e2',
  '2L'='e3',
  '3R'='e4',
  '3L'='e5'
))

elements <- paste0('e', c(1:5))
# prepare table for GRIMM_synteny by chromosome
GRIMM_tables.elX <- lapply(elements, function(elX){
  GRIMMtable.chr <- GRIMMtable.full[GRIMMtable.full$atr_scf==elX,]
})

# Make and set dirs
setDir <- function(){
  dirs <- list(
    GRIMM_by_chr.dir = './data/processed/GRIMM_by_chr',
    GRIMM_by_chr.dir.anchors = './data/processed/GRIMM_by_chr/anchors',
    GRIMM_by_chr.dir.blocks = './data/processed/GRIMM_by_chr/blocks'
  )
  
  dirs <- lapply(dirs, function(dir){
    if(!dir.exists(dir)){
      system(paste0('mkdir ', dir))
      dir
    } else {
      dir
    }
  })
}

function(){
  report <- readLines(paste0(dirs$GRIMM_by_chr.dir.blocks, '/report.txt'))
  return(report[1])
}()

# Make grimm
make_grimm <- function(
  GRIMM_tables.elX,
  grimm_path = '~/Documents/GRIMM_SYNTENY-2.02/grimm_synt',
  min_block = '2', min_gap='115000',
  dirs = setDir()
){
  
  reports <- lapply(GRIMM_tables.elX, function(GRIMM_table){
    
    # Write table with only one element
    grimm_table.path <- paste0(dirs$GRIMM_by_chr.dir, '/grimm_table.txt')
    write.table(GRIMM_table, file=grimm_table.path, quote=FALSE, row.names= FALSE)
    
    # Make GRIMM_anchors
    anchor_command <- paste(grimm_path, '-A -f', grimm_table.path, '-d', dirs$GRIMM_by_chr.dir.anchors)
    system(anchor_command, wait=TRUE)

    # Make GRIMM_synteny
    blocks_command <- paste0(
      grimm_path, ' -f ',
      dirs$GRIMM_by_chr.dir.anchors, '/unique_coords.txt', ' -d ',
      dirs$GRIMM_by_chr.dir.blocks,
      ' -c -n', min_block, ' -g ', min_gap, ' -Q'
    )
    system(blocks_command, wait=TRUE)
    # 
    # # Read blocks.txt
    report <- readLines(paste0(dirs$GRIMM_by_chr.dir.blocks, '/report.txt'))

    return(report)
  })
  return(reports)
}

GRIMMreports <- make_grimm(GRIMM_tables.elX)

# Extract distance tables
distance_tables <- lapply(GRIMMreports, function(report){
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
