#' This function produces the UI to process medium data analyzed with TraceFinder.
#' @author Daniel Braas
#' @param id The module id used in the app.
#' @return The processed data in the selected folder.
#' @export

medium_analysisUI <- function(id){

  ns <- NS(id)

  tagList(

    tabPanel(
      title = 'Medium',
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

# Process data -----------------------------------------------------------

        actionButton(inputId = ns('prep_all_RelA'),
                     label = 'Produce all relative amounts',
                     icon = icon('chart-bar'),
                     width = '100%'),
        actionButton(inputId = ns('prep_RelA'),
                     label = 'Produce relative amounts',
                     icon = icon('gitter'),
                     width = '100%'),
        actionButton(inputId = ns('make_MID'),
                     label = 'Calculate Mass Isotopologue Distribution',
                     icon = icon('code-branch'),
                     width = '100%'),
        hr(),
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
        hr(),
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
        verbatimTextOutput(ns('all_RelA')) %>%
          withSpinner(type = 1),
        verbatimTextOutput(ns('RelA')) %>%
          withSpinner(type = 1),
        hr(),
        verbatimTextOutput(ns('dat.MID')) %>%
          withSpinner(type = 1),
        hr(),
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


#' This is the corresponding server function to medium_analysisUI.
#' @author Daniel Braas
#' @return The processed data in the selected folder.
#' @export
medium_analysis <- function(input, output, session,
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
      print(volumes$Raw_data$Medium)
      gather_data(Title = input$title,
                  info = info_dat(),
                  save_raw = input$save_raw,
                  save_dropbox = input$save_dropbox,
                  drop_vol = volumes$Raw_data$Medium) %>%
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

# Produce relative amounts of all samples ---------------------------------

  amounts <- eventReactive(input$prep_all_RelA,{
    req(data2(), input$title)

    medium_RelA(dat = data2(),
                Title = input$title)
  })

  output$all_RelA <- renderText(
    paste0('All relative amounts processed: ',
           dim(amounts())[1], ' observations and ',
           dim(amounts())[2], ' variables.')
  )

# Produce relative amounts relative to fresh medium -----------------------

  relamounts <- eventReactive(input$prep_RelA,{
    req(amounts(), input$title)

    medium_relative(dat = amounts(),
                Title = input$title,
                info = info_dat())
  })

  output$RelA <- renderText(
    paste0('All relative amounts processed: ',
           dim(relamounts())[1], ' observations and ',
           dim(relamounts())[2], ' variables.')
  )

# Produce mass isotopologue distributions ---------------------------------

  dat.MID <- eventReactive(input$make_MID,{
    req(data2(), input$title)
    MID <- make_MID(data2(), Title = input$title)

    # Select only metabolites with isotopologue data
    metab <- filter(MID, Iso == 'M1') %>%
        .$Name %>%
        as.character()

    pdf(paste0(input$title, '_MID plots.pdf'), width = 14, height = 10)
    plot <- bar(unique(MID$Name), filter(MID, Name %in% metab))
    print(plot)
    dev.off()

    return(MID)
  })

  output$dat.MID <- renderText({
    req(dat.MID())
    paste0('Mass Isotopologue Distributions processed: ',
           dim(dat.MID())[1], ' observations and ',
           dim(dat.MID())[2], ' variables.')
  })
# Make heatmaps -----------------------------------------------------------

  observeEvent(input$make_Heatmap, {
    req(input$title)

    # Relative amounts
    if(!is.null(relamounts())){

      mat <- make_matrix(relamounts(), Title = input$title)

      # Only pass in cell samples not fresh medium
      info <- filter(info_dat(),
                     !str_detect(Condition,
                                 'unspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia'))

      norm <- norm()[info$Sample]

      make_pheatmap(matrix = mat,
                    samples = info,
                    heat.color = c('navy','grey90','firebrick3'),
                    Norv = norm,
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

    # # FC
    # if(!is.null(dat.FC())){
    #
    #   mat <- make_matrix(dat.FC(), Title = input$title)
    #
    #   make_pheatmap(matrix = mat,
    #                 samples = info_dat(),
    #                 heat.color = c('navy','grey90','firebrick3'),
    #                 Norv = norm(),
    #                 Title = input$title,
    #                 folder = getwd())
    # }

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
    if(!is.null(relamounts())){

      mat <- make_matrix(relamounts(), Title = input$title)

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

    # # FC
    # if(!is.null(dat.FC())){
    #
    #   mat <- make_matrix(dat.FC(), Title = input$title)
    #
    #   make_PCA(matrix = mat,
    #            info = info_dat(),
    #            a = as.numeric(input$PC1),
    #            b = as.numeric(input$PC2),
    #            cutoff = 0.5,
    #            Title = input$title,
    #            folder = getwd())
    # }

    output$done.PCA <- renderText({
      input$make_PCA
      if(sum(str_detect(list.files(), '_PC[1-5] vs. PC[1-5]_')) > 0) {
        return('PCA plots have been created.')
      }
    })
  })
  observe(output$done.PCA <- renderText(''))

# Saving results ----------------------------------------------------------

  observeEvent(input$save_all, {

    # Move all data to global environment
    try(info <<- info_dat())
    try(dat2 <<- data2())
    try(amount <<- amounts())
    try(amounts3 <<- relamounts())
    try(data3 <<- dat.MID())
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
