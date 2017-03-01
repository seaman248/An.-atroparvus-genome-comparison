library(shiny)

shinyUI(fluidPage(
  column(4, wellPanel(
    checkboxGroupInput(
      'sps_checkbox', 'Species',
      c('An. albimanus'='alb', 'An. atroparvus'='atr', 'An. gambiae'='gam'),
      selected = c('alb', 'atr', 'gam')
      ),
    checkboxGroupInput(
      'el_checkbox', 'Chromosome arms',
      c('e1'=1, 'e2'=2, 'e3'=3, 'e4'=4, 'e5'=5),
      selected = c(1:5)
    )
  )),
  column(8,
         plotOutput('genPlot')
         #tableOutput('genTable')
         )
))