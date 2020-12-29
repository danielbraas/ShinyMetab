library(tidyverse)
library(ShinyMetab)
library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(DT)

conf <- jsonlite::fromJSON(read_lines('~/config.json'))

ui <- navbarPage(
  title = 'Modules',
  theme = shinythemes::shinytheme('cerulean'),

  tabPanel(
    title = 'Medium',
    medium_analysisUI(id = 'exp')
  )
)

server <- function(input, output, session){

  volumes <- jsonlite::fromJSON(read_lines('~/config.json'))

  callModule(
    module = medium_analysis,
    id = 'exp',
    volumes = volumes$Shiny
  )
}

shinyApp(ui, server)
