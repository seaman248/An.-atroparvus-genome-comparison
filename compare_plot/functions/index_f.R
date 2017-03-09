library(lattice)
library(randomcoloR)
library(fields)
# ############################################################
devide_tables <- function(table){
  sp_tables <- lapply(seq(1, ncol(table), 5), function(range){
    range <- seq(range, range+4)
    table[, range]
  })
  names(sp_tables) <- paste0('Sp_', 1:length(sp_tables))
  return(sp_tables)
}

# ############################################################
orth_matches <- function(tlist, main_sp){
  # Declare x axis of matrix
  x_genes <- tlist[[main_sp]]
  
  # Sort x by start coordinate
  x_sorted <- x_genes[order(x_genes[,3]), ]
  
  # Generate compare matrix for every sp for compare with main sp
  orth_matches <- lapply(tlist[-main_sp], function(y_genes){
    # Sort y by start coordinate
    y_sorted <- y_genes[order(y_genes[,3]), ]
    
    # Compute matches
    matches <- match(row.names(y_sorted), row.names(x_sorted))
    
  })
  return(orth_matches)
}

# ############################################################
compare_matrixes <- function (tlist, main_sp = 1){
  # Declare x axis of matrix
  x_genes <- tlist[[main_sp]]
  
  # Sort x by start coordinate
  x_sorted <- x_genes[order(x_genes[,3]), ]
  
  # Generate compare matrix for every sp for compare with main sp
  matrixes_list <- lapply(tlist[-main_sp], function(y_genes){
    # Sort y by start coordinate
    y_sorted <- y_genes[order(y_genes[,3]), ]
    
    # Compute matches
    matches <- match(row.names(y_sorted), row.names(x_sorted))
    
    # # Generate list of row that will be combine into matrix
    # matrix_dim <- nrow(x_sorted)
    # list_of_row <- lapply(matches, function(match){
    #   row <- rep(0, matrix_dim)
    #   row[match] <- y_sorted[match, 2]
    #   return(row)
    # })
    # 
    # do.call('rbind', list_of_row)
  })
  
  return(matrixes_list)
}

# ############################################################
compare_plot <- function (table, main_sp=1, sps){
  # Make list of gene tables for each sp
  tlist <- devide_tables(table)
  names(tlist) <- sps
  # Make matrixes for generate plot
  matrixes_list <- compare_matrixes(tlist, main_sp)
  
  # Define parameters for plot
  titles <- paste0(names(tlist[main_sp]), '/', names(tlist[-main_sp]))
  colours <- c('white', randomColor(5))
  x_labels <- tlist[[main_sp]][,1]
  x_label_coords <- seq(0, 1, 1/(length(x_labels)-1))
  y_labels <- tlist[-main_sp]
  
  # Generate plot for each matrix
  par(mfrow=c(1, length(matrixes_list)))
  
  mapply(function(matrix, title, y_label){
    image(matrix, col = colours, axes = F, asp=1)
    axis(1, labels=x_labels, at=x_label_coords, las=2, cex=0.3, line=0.3)
    axis(2, labels=y_label[,1], at=x_label_coords, las=1, cex=0.3, line=0.3)
    title(title)
    image.plot(matrix, legend.only = TRUE, col = colours)
  }, matrixes_list, titles, y_labels)
}