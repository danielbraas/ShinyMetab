
# data_viewer -------------------------------------------------------------
#' This UI function contains allows both viewing as well as production of graphs.
#' The finished UI will contain multiple possibilities starting with bar
#' graphs of relative or absolute data, as well as
#' @author Daniel Braas
#' @param id The module ID.
#' @param label The label to be displayed depending on the type of data.
#' @param Title The text field to be displayed when opening the file.
#' @return A bar graph depicting metabolomics data based on the selected pathway.
#' @export

data_ViewerUI <- function(id, label, Title){

  ns <- NS(id)

  tagList(

    shinyFilesButton(id = ns('file'),
                   label = label,
                   title = Title,
                   multiple = F,
                   icon = icon('folder-open')),
    hr(),
    selectizeInput(inputId = ns('pathway'),
                   label = 'Select pathway',
                   choices = c('Choose'='','Glycolysis','Amino Acids','TCA','Adenine',
                               'CoAs','Currency Metabolites','Cysteine','Cytosine',
                               'Fatty Acid Metabolism','Fructose','Guanine','Hexose',
                               'Neurotransmitter','Pentose Phosphate Pathway','Thymine',
                               'Uracil'),
                   multiple = F),

    )
}

#' The corresponding server function to data_ViewerUI().
#' @author Daniel Braas
#' @export
data_Viewer <- function(input, output, session, volumes = volumes){

  shinyFileChoose(input = input,
                  id = 'file',
                  filetypes = c('rdata','rds'),
                  roots = volumes)

  dat <- reactive({

    if (!is.integer(input$file[1])){

      load(parseFilePaths(roots = volumes, input$file)$datapath)
      name <- ls()
      get(name) %>%
        ungroup() %>%
        return()
    }
  })

  output$plot <- renderPlot({
    req(dat())
    switch(input$pathway,
           Glycolysis = bar('glycolysis', dat()),
           'Amino Acids' = bar('AA', dat()),
           'TCA' = bar('TCA', dat()),
           'Adenine' = bar('Adenine', dat()),
           'CoAs' = bar('CoAs', dat()),
           'Currency Metabolites' = bar('Curr', dat()),
           'Cysteine' = bar('Cys', dat()),
           'Cytosine' = bar('Cytosine', dat()),
           'Fatty Acid Metabolism' = bar('FA', dat()),
           'Fructose' = bar('Fru', dat()),
           'Guanine' = bar('Guanine', dat()),
           'Hexose' = bar('Hex', dat()),
           'Neurotransmitter' = bar('Neurotrans', dat()),
           'Pentose Phosphate Pathway' = bar('PPP', dat()),
           'Thymine' = bar('Thymine', dat()),
           'Uracil' = bar('Uracil', dat()))
  }, height = function() {
    session$clientData$output_plot_width / 12 * 10
  }
  )
}

