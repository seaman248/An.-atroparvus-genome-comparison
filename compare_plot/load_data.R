sps <- read.csv2('./compare_plot/input_data/coords.csv')

colnames(sps) <- 1:ncol(sps)

sps[, c(3, 8, 13, 18)] <- lapply(c(3, 8, 13, 18), function(col_n){
  gsub('.+-', '', sps[, col_n])
})

sps[, c(6, 11, 16, 21)] <- lapply(c(6, 11, 16, 21), function(col_n){
  as.numeric(paste0(sps[, col_n], '1'))
})


sps_e4 <- na.omit(sps[sps[, c(3, 8, 13, 18)] == 'e4', ])

sps_e4 <- sps_e4[order(sps_e4[, 9]*-1), ]

devide_tables(sps_e4)
