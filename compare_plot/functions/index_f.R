library(lattice)
library(randomcoloR)
library(fields)
library(SDMTools)
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
orth_matches <- function(tlist, main_sp, order){
  # Declare x axis of matrix
  x_genes <- tlist[[main_sp]]
  x_order <- order[[main_sp]]
  # Sort x by start coordinate
  x_sorted <- x_genes[order(x_genes[,3]* x_order), ]
  
  # Generate compare matrix for every sp for compare with main sp
  orth_matches <- mapply(function(y_genes, y_order){
    # Sort y by start coordinate
    y_sorted <- y_genes[order(y_genes[,3]*y_order), ]
    
    # Compute matches
    matches <- match(row.names(y_sorted), row.names(x_sorted))
    
  }, tlist[-main_sp], order[-main_sp], SIMPLIFY = FALSE)
  
  return(orth_matches)
}

# ############################################################
compare_matrix <- function (matches, chrs){
  # Define size of matrix
  matrix_size <- length(matches)
  
  # Generate matrix as list of rows
  result_matrix <- sapply(matches, function(match){
    row <- rep(0, matrix_size)
    row[match] <- 1
    return(row)
  })
  
  return(as.matrix(result_matrix))
}


# ############################################################
compare_plot <- function (tlist, matches, main_sp=1){
  
  # Make matrixes for generate plot
  matrixes_list <- lapply(names(tlist[-main_sp]), function(sp_name){
    m <- compare_matrix(matches[[sp_name]], tlist[[sp_name]][,2])
    #m <- ConnCompLabel(m)
    m
  })
  
  # Define parameters for plot
  titles <- paste0(names(tlist[main_sp]), '/', names(tlist[-main_sp]))
  par(mfrow=c(1, length(matrixes_list)))

  mapply(function(matrix, title){
    image(matrix, col=c('white', 'black'), axes = F, asp=1)
    # axis(1, labels=x_labels$labels, at=x_labels$coords, las=2, cex=0.3, line=0.3)
    # axis(2, labels=y_labels$labels, at=y_labels$coords, las=1, cex=0.3, line=0.3)
    title(title)
    #image.plot(matrix, legend.only = TRUE, col = colours)
  }, matrixes_list, titles)
}