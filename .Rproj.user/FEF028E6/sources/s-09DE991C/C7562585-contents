library(tidyverse)
library(ShinyMetab)
library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(DT)

ui <- navbarPage(
  title = 'Metabolomics Workflow',
  theme = shinythemes::shinytheme('cerulean'),

# Clear memory ------------------------------------------------------------

  # tabPanel(
  #   title = 'New',
  #   actionButton(inputId = 'new',
  #                label = 'Clear memory',
  #                icon = icon('file'))
  # ),

# Data Converter ----------------------------------------------------------

  tabPanel(
    title = 'Data Converter',
    icon = icon('atom')
  ),

  # Data creator ------------------------------------------------------------

  navbarMenu(
    title = 'Data Creator',
    icon = icon('hammer'),

    tabPanel(
      title = 'Extracts',
      extract_analysisUI(id = 'extract')
    ),
    tabPanel(
      title = 'Medium',
      medium_analysisUI(id = 'medium')
    )
  ),

  # Plot Viewer -------------------------------------------------------------

  navbarMenu(
    title = 'Plot Viewer',
    icon = icon('eye'),

    tabPanel(
      title = 'Relative Amounts',
      column(
        width = 2,
        data_ViewerUI(id = 'RelA',
                      label = 'Choose relative amount file',
                      Title = 'Relative Amounts')
      ),
      column(
        width = 10,
        plotOutput('RelA_plot')
      )
    ),
    tabPanel(
      title = 'MID',
      column(
        width = 2,
        data_ViewerUI(id = 'MID',
                      label = 'Choose MID file',
                      Title = 'MIDs')
      ),
      column(
        width = 10,
        plotOutput('MID_plot')
      )
    ),
    tabPanel(
      title = 'Fractional contribution',
      column(
        width = 2,
        data_ViewerUI(id = 'FC',
                      label = 'Choose FC file',
                      Title = 'FC')
      ),
      column(
        width = 10,
        plotOutput('FC_plot')
      )
    ),
    tabPanel(
      title = 'Percent labeled',
      column(
        width = 2,
        data_ViewerUI(id = 'labeled',
                      label = 'Choose percent labeled file',
                      Title = 'labeled')
      ),
      column(
        width = 10,
        plotOutput('labeled_plot')
      )
    )
  )
)


# server ------------------------------------------------------------------

server <- function(input, output, session){

  volumes <- jsonlite::fromJSON(read_lines('~/config.json'))$Shiny

# clear memory ------------------------------------------------------------

  observeEvent(input$new,{
    print('this works')
    rm(list = ls())
  })

# extract_analysis --------------------------------------------------------

  callModule(
    module = extract_analysis,
    id = 'extract',
    volumes = volumes
  )

# medium_analysis --------------------------------------------------------

  callModule(
    module = medium_analysis,
    id = 'medium',
    volumes = volumes
  )


# data_viewer -------------------------------------------------------------


  output$RelA_plot <- callModule(module = data_Viewer,
                                 id = 'RelA',
                                 volumes = volumes)

  output$MID_plot <- callModule(module = data_Viewer,
                                id = 'MID',
                                volumes = volumes)

  output$FC_plot <- callModule(module = data_Viewer,
                                id = 'FC',
                                volumes = volumes)

  output$labeled_plot <- callModule(module = data_Viewer,
                               id = 'labeled',
                               volumes = volumes)

}

shinyApp(ui, server)
