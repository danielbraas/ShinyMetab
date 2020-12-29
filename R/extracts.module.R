# extract_analysisUI ------------------------------------------------------
#' This function produces the UI to process data analyzed with TraceFinder.
#' @author Daniel Braas
#' @param id The module id used in the app.
#' @return The processed data in the selected folder.
#' @export

extract_analysisUI <- function(id){

  ns <- NS(id)

  tagList(

    tabPanel(
      title = 'Extract',
      column(
        width = 2,

# Experiment setup --------------------------------------------------------

        p(style = "font-weight:bold",
          'Select experiment folder'),
        shinyDirButton(id = ns('Dir'),
                       label = 'Select folder',
                       title = 'Experiment folder',
                       icon = icon('folder-open'),
                       style = "width:100%"),

        shinyFilesButton(id = ns('info'),
                         label = 'Choose experiment info file',
                         title = 'Experiment info',
                         multiple = F,
                         icon = icon('info-circle'),
                         style = "width:100%"),
        actionButton(inputId = ns('show_info'),
                     label = 'Show/hide info',
                     icon = icon('eye'),
                     width = '100%'),
        checkboxInput(inputId = ns('save_info'),
                      label = 'Save info',
                      value = TRUE),
        checkboxInput(inputId = ns('save_dropbox'),
                      label = 'Save raw on Dropbox',
                      value = TRUE),
        checkboxInput(inputId = ns('save_raw'),
                      label = 'Save raw data',
                      value = TRUE),
        checkboxInput(inputId = ns('load_raw'), # In case this step has been performed before
                      label = 'Load gathered data',
                      value = FALSE),
        textInput(inputId = ns('title'),
                  label = 'Experiment title',
                  value = '',
                  width = "100%"),
        actionButton(inputId = ns('use_title'),
                     label = 'Use',
                     icon = icon('thumbs-up'),
                     width = '100%'),
        hr(),

# Prepare MS data ------------------------------------------------------------

        p(style = "font-weight:bold",
          'Load MS data'),
        checkboxInput(inputId = ns('create_Std'),
                      label = 'Produce Standards plots',
                      value = TRUE),
        selectInput(inputId = ns('normalizer'),
                    label = 'Select normalizer',
                    choices = 'none'),
        actionButton(inputId = ns('prep_data'),
                     label = 'Prepare MS data',
                     icon = icon('flask'),
                     width = '100%'),
        hr(),

# Process data ------------------------------------------------------------

        p(style = "font-weight:bold",
          'Process metabolomics data'),
        actionButton(inputId = ns('make_RelAmounts'),
                     label = 'Calculate relative amounts',
                     icon = icon('balance-scale-left'),
                     width = '100%'),
        actionButton(inputId = ns('make_MID'),
                     label = 'Calculate Mass Isotopologue Distribution',
                     icon = icon('code-branch'),
                     width = '100%'),
        actionButton(inputId = ns('make_FC'),
                     label = 'Calculate Fractional Contribution',
                     icon = icon('gulp'),
                     width = '100%'),
        actionButton(inputId = ns('make_PercentLabeled'),
                     label = 'Calculate percent labeled',
                     icon = icon('glass-whiskey'),
                     width = '100%'),
        hr(),

# Graphical output --------------------------------------------------------

        p(style = "font-weight:bold",
          'Graphical output'),
        actionButton(inputId = ns('make_Plots'),
                     label = 'Bar graphs',
                     icon = icon('chart-bar'),
                     width = '100%'),
        # actionButton(inputId = ns('make_volcano'),
        #              label = 'Volcano plot',
        #              icon = icon('brain'),
        #              width = '100%'),
        actionButton(inputId = ns('make_Heatmap'),
                     label = 'Heatmap',
                     icon = icon('braille'),
                     width = '100%'),
        selectInput(inputId = ns('PC1'),
                    label = '1st PC:',
                    choices = 1:5,
                    selected = 1),
        selectInput(inputId = ns('PC2'),
                    label = '2nd PC:',
                    choices = 1:5,
                    selected = 2),
        actionButton(inputId = ns('make_PCA'),
                     label = 'PCA',
                     icon = icon('hubspot'),
                     width = '100%'),
        actionButton(inputId = ns('save_all'),
                     label = 'Save all data',
                     icon = icon('save'),
                     width = '100%')
      ),

# Console -----------------------------------------------------------------

      column(
        width = 10,
        h5('Experiment folder:'),
        verbatimTextOutput(ns('DIR')),
        h5('Experiment title:'),
        verbatimTextOutput(ns('TITLE')),
        DTOutput(ns('info')),
        hr(),
        h5('Raw data files in experiment folder:'),
        verbatimTextOutput(ns('raw_files')),
        plotOutput(ns('NORM')) %>%
          withSpinner(type = 1),
        verbatimTextOutput(ns('preprocessed')) %>%
          withSpinner(type = 1),
        hr(),
        verbatimTextOutput(ns('dat.RelA')) %>%
          withSpinner(type = 1),
        verbatimTextOutput(ns('dat.MID')) %>%
          withSpinner(type = 5),
        verbatimTextOutput(ns('dat.FC')) %>%
          withSpinner(type = 1),
        verbatimTextOutput(ns('dat.lab')) %>%
          withSpinner(type = 1),
        hr(),
        verbatimTextOutput(ns('done.plots')) %>%
          withSpinner(type = 1),
        verbatimTextOutput(ns('done.heatmaps')) %>%
          withSpinner(type = 1),
        verbatimTextOutput(ns('done.PCA')) %>%
          withSpinner(type = 1),
        verbatimTextOutput(ns('saved.all')) %>%
          withSpinner(type = 1)
      )

    )
  )
}

# extract_analysis server -------------------------------------------------
#' This is the corresponding server function to the extract_analysisUI function.
#' @author Daniel Braas
#' @return The processed data in the selected folder.
#' @export
extract_analysis <- function(input, output, session,
                             volumes = volumes){

# Experiment setup --------------------------------------------------------

  # Navigate to project folder
  shinyDirChoose(input = input,
                 id = 'Dir',
                 roots = volumes,
                 defaultRoot = 'R_Projects')

  # Select and change into directory
  observe({
    req(input$Dir)
    dir <- parseDirPath(roots = volumes, input$Dir)
    if(length(dir > 0)) setwd(dir)
    output$DIR <- renderText(dir)
  })

  # Pick sample information file
  shinyFileChoose(input = input,
                  id = 'info',
                  filetypes = c('csv', 'xls', 'xlsx', 'rds'),
                  roots = c('Exp' = '.')) # Have to find a way to pass folder name (or just top level?) into quotes

  # Load sample info and prepare Sample.Name column
  info_dat <- reactive({
    req(input$Dir)
    if (!is.integer(input$info[1])){
      name <- parseFilePaths(volumes, input$info)$name

      if(grepl('csv|CSV', name)){
        info <- read_csv(parseFilePaths(volumes, input$info)$name,
                         #col_names = input$heading,
                         #na = input$na.string)
                         na = c('','NA'))
      }

      if(str_detect(name, 'xls|XLS')){
        info <- readxl::read_excel(parseFilePaths(volumes, input$info)$name)
      }

      if(str_detect(name, 'rds')){
        info <- read_rds(parseFilePaths(volumes, input$info)$name)
      }

      # Process info data table
      info <- mutate(info, Condition = str_replace_all(Condition, '^\\s+|\\s+$', ''),  # remove leading or trailing white space
                     Condition = str_replace_all(Condition, '/|_', '-'),
                     Condition = factor(Condition, levels = unique(Condition))) %>%
        group_by(Condition) %>%
        mutate(Sample.Name = str_c(Condition, row_number(), sep = '_')) %>%
        ungroup()
      return(info)
    }
  })

  # Display sample info in console
  output$info <- renderDT({
    req(input$info)

    # Decide if to show info or not
    if(input$show_info %% 2 == 0){
      datatable(info_dat(),
                editable = 'cell')
    }
  })

  # Experiment title
  observe({
    req(input$Dir, input$info)

    val <- parseFilePaths(volumes, input$info)$name %>%
      str_remove(., '.xls(x)?') %>% # Remove extensions
      str_remove(., '.csv') %>% # Remove extensions
      str_remove(., '.rds') %>%
      str_remove(., '_sample_info') # Remove sample info appendix if sample info file had been created

    updateTextInput(session = session,
                    inputId = 'title',
                    label = 'Title',
                    value = val)
  })

  # Print title
  output$TITLE <- renderText(input$title)

# Load MS data ------------------------------------------------------------
  # Does raw data file exist and update selection if it does
  observe({
    req(input$Dir)
    if(file.exists('all_data_raw.csv') | file.exists('all_data_raw.rds')){
      updateCheckboxInput(session = session,
                          inputId = 'load_raw',
                          value = TRUE)
    }

  })

  # Show raw data files if available
  output$raw_files <- renderPrint({
    req(input$Dir)
    list.files(pattern='.csv', recursive=T) %>%
      .[str_detect(., '(.)*/.(.)*.csv')] # Selects only csv files in subfolders; might be suffiecient to select the right files
  })

  # Gather data
  data <- eventReactive(input$use_title, {
    req(input$Dir, input$title, info_dat())

    if(input$load_raw){
      if(file.exists('all_data_raw.rds')){
        read_rds('all_data_raw.rds') %>%
          shape_raw()
      } else {
        suppressWarnings(read_csv('all_data_raw.csv')) %>%
          shape_raw()
      }
    } else {
      print(volumes$Raw_data$Extracts)
      gather_data(Title = input$title,
                  info = info_dat(),
                  save_raw = input$save_raw,
                  save_dropbox = input$save_dropbox,
                  drop_vol = volumes$Raw_data$Extracts) %>%
        shape_raw()
    }

  })

# Prepare data for analysis -----------------------------------------------

  # Choose normalizer
  observe({
    req(data())

    # Find external normalizer
    std <- c(data()$Compound[str_detect(data()$Compound, '[Nn]orv')],
             data()$Compound[str_detect(data()$Compound, '[Ss]td')],
             'none') %>%
      unique()

    updateSelectInput(session = session,
                      inputId = 'normalizer',
                      choices = std)
  })

  # Make a line plot of the response of the normalizer
  output$NORM <- renderPlot({
    req(data(), input$normalizer)

    data() %>%
      filter(Compound == input$normalizer) %>%
      select(matches('-[0-9]+$')) %>%
      gather(Sample, Value) %>%
      mutate(Sample=factor(Sample, levels=unique(Sample))) %>%
      select(Sample, Value) %>%
      ggplot(., aes(Sample, Value)) +
      geom_point(size=3) +
      geom_line(aes(as.integer(Sample), Value), color='blue') +
      theme_bw()+
      theme(text = element_text(size = 16, face = 'bold'))+
      labs(x='Sample Name', y='Response', title=paste0(input$normalizer, ' Response')) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  })

  # Pull normalizer values for heatmap
  norm <- reactive({
    req(data())

    if(input$normalizer == 'none'){
      return(1)
    } else {
      filter(data(), Compound == input$normalizer) %>%
        select(matches('-[0-9]+$')) %>%
        gather(Sample, Value) %>%
        pull(Value)
    }
  })

  # Preprocess data
  data2 <- eventReactive(input$prep_data, {
    req(data(), input$normalizer)

    prepare_data(data = data(),
                 Title = input$title,
                 Sample.Name = info_dat()$Sample.Name,
                 numbers = info_dat()$Cell.Number,
                 std = FALSE,
                 bkg = FALSE,
                 norm = input$normalizer)
  })

  output$preprocessed <- renderText(
    paste0('Data pre-processed: ',
           dim(data2())[1], ' observations and ',
           dim(data2())[2], ' variables.')
    )

# Process data ------------------------------------------------------------

  dat.RelA <- eventReactive(input$make_RelAmounts,{
    req(info_dat(), data2(), input$title)
    make_RelAmounts(data2(), Title = input$title)
  })

  output$dat.RelA <- renderText({
    paste0('Relative Amounts processed: ',
           dim(dat.RelA())[1], ' observations and ',
           dim(dat.RelA())[2], ' variables.')
    })

  dat.MID <- eventReactive(input$make_MID,{
    req(info_dat(), data2(), input$title)
    make_MID(data2(), Title = input$title)
  })

  output$dat.MID <- renderText({
    req(dat.MID())
    paste0('Mass Isotopologue Distributions processed: ',
           dim(dat.MID())[1], ' observations and ',
           dim(dat.MID())[2], ' variables.')
  })

  dat.FC <- eventReactive(input$make_FC,{
    req(info_dat(), dat.MID(), input$title)
    # Should add feedback if dat.MID doesn't exist
    make_FC(dat.MID(), Title = input$title)
  })

  output$dat.FC <- renderText({ # Insert spinner logic
    paste0('Fractional contributions processed: ',
           dim(dat.FC())[1], ' observations and ',
           dim(dat.FC())[2], ' variables.')
  })

  dat.lab <- eventReactive(input$make_PercentLabeled,{
    req(info_dat(), dat.MID(), input$title)
    make_labeled(dat.MID(), Title = input$title)
  })

  output$dat.lab <- renderText({ # Insert spinner logic
    paste0('Percent labeled processed: ',
           dim(dat.lab())[1], ' observations and ',
           dim(dat.lab())[2], ' variables.')
  })

# Make heatmaps -----------------------------------------------------------

  observeEvent(input$make_Heatmap, {
    req(input$title)

    # Relative amounts
    if(!is.null(dat.RelA())){

      mat <- make_matrix(dat.RelA(), Title = input$title)

      make_pheatmap(matrix = mat,
                    samples = info_dat(),
                    heat.color = c('navy','grey90','firebrick3'),
                    Norv = norm(),
                    Title = input$title,
                    folder = getwd())
    }

    # MIDs
    if(!is.null(dat.MID())){

      mat <- make_matrix(dat.MID(), Title = input$title)

      make_pheatmap(matrix = mat,
                    samples = info_dat(),
                    heat.color = c('navy','grey90','firebrick3'),
                    Norv = norm(),
                    Title = input$title,
                    folder = getwd())
    }

    # FC
    if(!is.null(dat.FC())){

      mat <- make_matrix(dat.FC(), Title = input$title)

      make_pheatmap(matrix = mat,
                    samples = info_dat(),
                    heat.color = c('navy','grey90','firebrick3'),
                    Norv = norm(),
                    Title = input$title,
                    folder = getwd())
    }

    # Percent Labeled
    if(!is.null(dat.lab())){

      mat <- make_matrix(dat.lab(), Title = input$title)

      make_pheatmap(matrix = mat,
                    samples = info_dat(),
                    heat.color = c('navy','grey90','firebrick3'),
                    Norv = norm(),
                    Title = input$title,
                    folder = getwd())
    }


  })

  output$done.heatmaps <- renderText({
    input$make_Heatmap
    if(sum(str_detect(list.files(), '_Heatmap_')) > 0) {
      return('Heatmaps have been created.')
    }
  })


# Make PCA plots ----------------------------------------------------------

  observeEvent(input$make_PCA, {
    req(input$title)

    # Relative amounts
    if(!is.null(dat.RelA())){

      mat <- make_matrix(dat.RelA(), Title = input$title)

      make_PCA(matrix = mat,
               info = info_dat(),
               a = as.numeric(input$PC1),
               b = as.numeric(input$PC2),
               cutoff = 0.5,
               Title = input$title,
               folder = getwd())
    }

    # MIDs
    if(!is.null(dat.MID())){

      mat <- make_matrix(dat.MID(), Title = input$title)

      make_PCA(matrix = mat,
               info = info_dat(),
               a = as.numeric(input$PC1),
               b = as.numeric(input$PC2),
               cutoff = 0.5,
               Title = input$title,
               folder = getwd())
    }

    # FC
    if(!is.null(dat.FC())){

      mat <- make_matrix(dat.FC(), Title = input$title)

      make_PCA(matrix = mat,
               info = info_dat(),
               a = as.numeric(input$PC1),
               b = as.numeric(input$PC2),
               cutoff = 0.5,
               Title = input$title,
               folder = getwd())
    }

    # Percent Labeled
    if(!is.null(dat.lab())){

      mat <- make_matrix(dat.lab(), Title = input$title)

      make_PCA(matrix = mat,
               info = info_dat(),
               a = as.numeric(input$PC1),
               b = as.numeric(input$PC2),
               cutoff = 0.5,
               Title = input$title,
               folder = getwd())
    }
  })

  output$done.PCA <- renderText({
    input$make_PCA
    if(sum(str_detect(list.files(), '_PC[1-5] vs. PC[1-5]_')) > 0) {
      return('PCA plots have been created.')
    }
  })


# Making bar graphs -------------------------------------------------------

  observeEvent(input$make_Plots, {
    req(input$title)

    plotname <- paste0(input$title, '_plots.pdf')

    pdf(file = plotname, width=14, height=10, pointsize=12)
    try(print(bar('glycolysis', dat.RelA())))
    try(print(bar('glycolysis', dat.MID())))
    try(print(bar('glycolysis', dat.lab())))
    try(print(bar('glycolysis', dat.FC())))
    try(print(bar("TCA", dat.RelA())))
    try(print(bar("TCA", dat.MID())))
    try(print(bar("TCA", dat.lab())))
    try(print(bar("TCA", dat.FC())))
    try(print(bar("PPP", dat.RelA())))
    try(print(bar("PPP", dat.MID())))
    try(print(bar("PPP", dat.lab())))
    try(print(bar("PPP", dat.FC())))
    try(print(bar("Curr", dat.RelA())))
    try(print(bar("Curr", dat.MID())))
    try(print(bar("Curr", dat.lab())))
    try(print(bar("Curr", dat.FC())))
    try(print(bar("Cys", dat.RelA())))
    try(print(bar("Cys", dat.MID())))
    try(print(bar("Cys", dat.lab())))
    try(print(bar("Cys", dat.FC())))
    try(print(bar("AA", dat.RelA())))
    try(print(bar("AA", dat.MID())))
    try(print(bar("AA", dat.lab())))
    try(print(bar("AA", dat.FC())))
    try(print(bar("FA", dat.RelA())))
    try(print(bar("FA", dat.MID())))
    try(print(bar("FA", dat.lab())))
    try(print(bar("FA", dat.FC())))
    try(print(bar("Hex", dat.RelA())))
    try(print(bar("Hex", dat.MID())))
    try(print(bar("Hex", dat.lab())))
    try(print(bar("Hex", dat.FC())))
    try(print(bar("Adenine", dat.RelA())))
    try(print(bar("Adenine", dat.MID())))
    try(print(bar("Adenine", dat.lab())))
    try(print(bar("Adenine", dat.FC())))
    try(print(bar("Cytosine", dat.RelA())))
    try(print(bar("Cytosine", dat.MID())))
    try(print(bar("Cytosine", dat.lab())))
    try(print(bar("Cytosine", dat.FC())))
    try(print(bar("Guanine", dat.RelA())))
    try(print(bar("Guanine", dat.MID())))
    try(print(bar("Guanine", dat.lab())))
    try(print(bar("Guanine", dat.FC())))
    try(print(bar("Thymine", dat.RelA())))
    try(print(bar("Thymine", dat.MID())))
    try(print(bar("Thymine", dat.lab())))
    try(print(bar("Thymine", dat.FC())))
    try(print(bar("Uracil", dat.RelA())))
    try(print(bar("Uracil", dat.MID())))
    try(print(bar("Uracil", dat.lab())))
    try(print(bar("Uracil", dat.FC())))
    try(print(bar('Fru', dat.RelA())))
    try(print(bar('Fru', dat.MID())))
    try(print(bar("Fru", dat.lab())))
    try(print(bar('Fru', dat.FC())))
    try(print(bar('CoAs', dat.RelA())))
    try(print(bar('CoAs', dat.MID())))
    try(print(bar("CoAs", dat.lab())))
    try(print(bar('CoAs', dat.FC())))

    dev.off()

  })

  output$done.plots <- renderText({
    input$make_Plots
    if(sum(str_detect(list.files(), '_plots.pdf')) > 0) {
      return('Bar graphs have been created.')
    }
  })

# Saving results ----------------------------------------------------------

  observeEvent(input$save_all, {

    # Move all data to global environment
    try(info <<- info_dat())
    try(dat2 <<- data2())
    try(data4 <<- dat.RelA())
    try(data3 <<- dat.MID())
    try(FC <<- dat.FC())
    try(data_labeled <<- dat.lab())
    try(numbers <<- pull(info_dat(), 'Cell.Number'))
    try(Sample.Name <<- info_dat()$Sample.Name)

    # Replace saved image if it exists; not sure why I have to explicitly state this
    if(file.exists(paste0(input$title, '.rdata'))) file.remove(paste0(input$title, '.rdata'))
    save.image(paste0(input$title, '.rdata'))

    output$saved.all <- renderText({
      if(sum(str_detect(list.files(), str_c(input$title, '.rdata'))) > 0) {
            return('All data have been saved.')
          }
    })
  })

  # Make sure spinner of save all is not constantly on
  observe(output$saved.all <- renderText(''))

}
