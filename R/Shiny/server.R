library(shiny)
library(plyr)
library(dplyr)
setwd('/home/seaman/rproj/full_genome_comparison/')

source('./R/Visualize/functions/viz_fun.R')
source('./R/Visualize/functions/parse_blocks.R')

sps <- c('alb', 'atr', 'gam')
blocks_list <- parse_blocks(read.table('./R/GRIMM/output_data/blocks/blocks.txt'))
names(blocks_list) <- sps

shinyServer(function(input, output, session){
  values <- reactiveValues(blocks_list = blocks_list, els = 1:5)
  
  observe({
    values$blocks_list <- lapply(blocks_list[input$sps_checkbox], function(blocks){
      blocks[!is.na(match(blocks_list$atr$chr, paste0('e', input$el_checkbox))),]
    })
  })
  
  
  output$genPlot <- renderPlot({
    withProgress(message = 'Rendering plot...', {
      viz_fun(values$blocks_list, as.numeric(values$els))
    })
  })
  
  # output$genTable <- renderTable({
  #   withProgress({
  #     bind_cols(values$blocks_list)
  #   })
  # })
  
})