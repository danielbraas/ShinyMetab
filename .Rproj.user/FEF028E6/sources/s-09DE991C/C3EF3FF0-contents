
#' This function gathers the individual integration files produced by TraceFinder.
#' @author Daniel Braas
#' @param Title Project title
#' @param info A data frame with information about the metabolomics experiment perfromed.
#' @param save_raw (boolean)Should the raw data be saved.
#' @param save_dropbox (boolean) Should the raw data be saved on your Dropbox.
#' This parameter expects that the path to the Dropbox folder is known.
#' @param drop_vol The path to the Dropbox project folder to be used.
#' @return A data frame with integrated metabolomics data.
#' @export

gather_data <- function(Title,
                        info = info,
                        save_raw = TRUE,
                        save_dropbox = TRUE,
                        drop_vol = drop_vol){

  # Gather raw MS data ----------------------------------------------------------
  files <- list.files(pattern='.csv', recursive=T) %>%
    .[str_detect(., '(.)*/.(.)*.csv')] # Selects only csv files in subfolders; might be suffiecient to select the right files

  # Collect data
  dat <- map(files, read_csv) %>%
    map_df(~bind_rows(.)) %>%
    mutate(Experiment = Title)

  # Save raw data in experiment folder
  if(save_raw){
    write_rds(dat, 'all_data_raw.rds')
    write_rds(info, paste(Title, 'sample_info.rds', sep = '_'))
  }

  # Save raw data to Dropbox folder
  if(save_dropbox){
    # Make sure folder for all raw data exists
    if(!dir.exists(drop_vol)){
      dir.create(drop_vol)
    }
    write_rds(dat, paste(drop_vol,
                         paste(Title, 'all_data_raw.rds', sep = '_'),
                         sep = '/'))

    write_rds(info, paste(drop_vol,
                          paste(Title, 'sample_info.rds', sep = '_'),
                          sep = '/'))
  }

  return(dat)
}

#' This function separates the compound names into compound and isotopologue names.
#' @author Daniel Braas
#' @param dat A data frame with integrated metabolomics data.
#' @return A data frame where the compound names and isotopologue data have been separated and formated.
#' @export

shape_raw <- function(dat){

  # Reshape and format data
  select(dat, Compound, Filename, Area) %>%
    spread(Filename, Area) %>%
    mutate(Iso = if_else(str_detect(Compound, ' M[0-9]+$'),
                         str_extract(Compound, 'M[0-9]+$'),
                         'M0'),
           Iso = factor(Iso, levels=paste0('M',0:50)),
           Used_ID = str_remove(Compound, ' M[0-9]+$| Std')) %>%
    select(Used_ID, Iso, everything()) %>%
    arrange(Used_ID, Iso) %>%
    inner_join(., Abbrev, by = 'Used_ID')

}

#' prepare_data is a function gathers the individual integration files produced by TraceFinder.
#' @author Daniel Braas
#' @param data A formated data frame that contains the non-normalized MS data.
#' @param Title Project title
#' @param Sample.Name A vector containing the names for each sample. Usually,
#' this will be a concatenation of the condition and the sample repeat.
#' @param numbers A vector of cell numbers or other specimen normalization
#' quantity.
#' @param std (boolean) Should a plot of all internal standards be created.
#' @param bkg (boolean) Should the background be subtracted from all samples.
#' @param norm The standard (internal or external) to normalize to.
#' @return A data frame with integrated metabolomics data.
#' @export

prepare_data <- function(data,
                         Title,
                         Sample.Name = '',
                         numbers = 1,
                         std = FALSE,
                         bkg = FALSE,
                         norm = 'none'){

  # analyzing the metabolite standards --------------------------------------
  # This should be expanded significantly if normalizations to standards should happen
  if(std){
    Std <- filter(data, grepl('Std', data$Compound)) %>%
      select(-Iso) %>%
      gather(Sample, Value, -Used_ID) %>%
      filter(!(grepl('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]', Sample))) %>%
      group_by(Used_ID) %>%
      mutate(Av = mean(Value, na.rm=T),
             Std = sd(Value, na.rm=T),
             CV = Std / Av,
             ID = paste0(gsub(' Std','',Used_ID), ' ',round(CV*100,0),'%')) %>%
      ungroup()

    pdf('Relative response of ISTDs with CV.pdf', width=14, height=10)
    std <- ggplot(Std, aes(Sample, log(Value,10)))+
      geom_point(size=3)+
      facet_wrap(~ID, scales='free')+
      theme_bw()+
      theme(text = element_text(face='bold',size=14),
            axis.text.x = element_text(angle=90, vjust=.3, hjust=1))+
      labs(x='',y='Relative response of ISTD (A.U., log10)')
    print(std)
    dev.off()
  }

  # Background subtraction; needs work
  if(bkg){

    if (sum(grepl('blank|water|buffer', names(data))) > 0) {
      blank <- data[,grep('blank|water|buffer', names(data))]
      data <- data[,-grep('blank|water|buffer', names(data))]
      blank[is.na(blank)] <- 0
    }

    if (sum(grepl('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]|200[Kk]', names(data))) > 0) {                          ##collect all standard samples and remove from regular samples
      Standard <- data[,c(1:3, grep('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]|200[Kk]', names(data)))]
      data <- data[,-grep('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]|200[Kk]', names(data))]
    }

  }

  # Remove samples used for background
  data <- select(data, -which(str_detect(names(data), 'blank|water|buffer')),
                 -which(str_detect(names(data), 'k_cells')))

  # Sample normalization to external standard
  if(norm != 'none'){
    #Norm <- filter(data, Used_ID == norm) %>%
    Norm <- filter(data, Compound == norm) %>%
      select(matches('-[0-9]+$')) %>%
      gather(Sample, Norm) %>%
      .$Norm

    # Make normalizer plot
    Norm.title = paste0(norm, " Response-", Title, ".pdf")
    pdf(file = Norm.title, width=12, height=9, pointsize=12)
    plot <- data %>%
      filter(Compound == norm) %>%
      select(matches('-[0-9]+$')) %>%
      gather(Sample, Value) %>%
      mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
      select(Sample, Value) %>%
      ggplot(., aes(Sample, Value)) +
      geom_point(size = 3) +
      geom_line(aes(as.integer(Sample), Value), color = 'blue') +
      theme_bw()+
      theme(text = element_text(size = 16, face = 'bold'))+
      labs(x = 'Sample Name',
           y = 'Response',
           title = paste0(norm, ' Response')) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(plot)
    dev.off()
  } else {
    # If data is not to be normalized to a standard then assign 1 as normalizer
    Norm <- rep(1, length(numbers))
  }

  # Normalize data to normalizer and cell number
  mat <- select(data, matches('-[0-9]+$')) %>%
    as.matrix() %>%
    t(.) / Norm / numbers

  # Make data frame and replace sample columns in data
  data <- t(mat) %>%
    data.frame() %>%
    bind_cols(data, .) %>%
    select(-matches('-[0-9]+$')) %>%
    filter(!Used_ID == norm,
           !str_detect(Used_ID, 'Std'),
           !str_detect(Used_ID, '[Nn]orvaline')) # remove normalizer

  # reshape data for downstream analysis ------------------------------------
  data2 <- rename(data, Name = Abb) %>%
    select(Name, Iso, Used_ID, KEGG.ID, Nr.C, matches('.[0-9]+$')) %>%
    gather(Exp, Value, -Name, -Iso, -Used_ID, -KEGG.ID, -Nr.C) %>%
    mutate(Exp = plyr::mapvalues(Exp, # Label samples with sample name
                                 from = names(data)[str_detect(names(data), '.[0-9]+$')],
                                 to = Sample.Name)) %>%
    separate(Exp, c('Condition', 'Exp'), sep='_') %>%
    select(Name, Iso, Condition, Exp, Value, Used_ID, KEGG.ID, Nr.C) %>%
    mutate(Condition = factor(Condition, levels = unique(Condition)),
           Exp = paste0('Exp', Exp)) %>%
    group_by(Name) %>%
    mutate(NA_list = sum(is.na(Value)),
           NA_potential = n(),
           Name = as.character(Name)) %>%
    ungroup() %>%
    filter(NA_list != NA_potential) %>%
    mutate(NA_list = NULL,
           NA_potential = NULL,
           Name = if_else(is.na(Name), Used_ID, Name)) # Makes sure that names that didn't match Abbrev

  write.csv(data2, file = paste0(Title,"_all data normalized.csv"),
            row.names = FALSE)
  write_rds(data2, file = 'data2.rds')
}

# Prepare relative amounts including fresh medium------------------------------------------------

#' medium_RelA calculates statistical metrics from all replicates (including fresh medium blanks).
#' @author Daniel Braas
#' @param dat A formated data frame produced by prepare_data().
#' @param Title Project title.
#' @return A data frame with integrated metabolomics data as well as
#' bar graphs of all measured metabolites.
#' @export

medium_RelA <- function(dat,
                        Title = ''){
  amounts <- dat %>%
    select(Name, KEGG.ID, Condition, Iso, Exp, Value) %>%
    group_by(Name, Condition, Exp) %>%
    mutate(Amount = sum(Value, na.rm = T)) %>%
    ungroup() %>%
    filter(Iso == 'M0')

  # Need to convert 0 to NA for Av calculation
  amounts$Amount[amounts$Amount == 0] <- NA

  amounts <- amounts %>%
    group_by(Name, Condition) %>%
    mutate(Av = mean(Amount, na.rm = T),
      Std = sd(Amount, na.rm = T),
      CV = sd(Amount, na.rm = T) / mean(Amount, na.rm = T)) %>%
    select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
    ungroup() %>%
    arrange(Condition, Name)

  test1 <- split(amounts, amounts[c('Name', 'Condition')])

  NA_function <- function(x) if (sum(is.na(x)) < length(x))return(x = x)
  else return(x = rep(0, length(x)))

  new.Value <- as.vector(sapply(test1, function(x)NA_function(x$Amount)))
  if (class(new.Value) == 'list')new.Value <- unlist(new.Value)

  amounts <- amounts %>%
    mutate(Amount = new.Value)

  data8 = split(amounts, amounts[, 1])
  ANOVA = suppressWarnings(sapply(data8, function(x)
    anova(aov(x$Amount ~ x$Condition))$Pr[1]))
  ANOVA = rep(ANOVA, each = length(unique(amounts$Condition)))

  amounts <- amounts %>%
    arrange(Name) %>%
    spread(Exp, Amount) %>%
    cbind('ANOVA' = ANOVA, 'Sig' = NA)

  for (i in 1:nrow(amounts)){
    if (amounts$ANOVA[i] == "NaN") amounts$Sig[i] = ""
    else if (amounts$ANOVA[i] <= 0.001) amounts$Sig[i] = "***"
    else if (amounts$ANOVA[i] <= 0.01) amounts$Sig[i] = "**"
    else if (amounts$ANOVA[i] <= 0.05) amounts$Sig[i] = "*"
    else amounts$Sig[i] = ""
  }

  plotname = paste(Title,"_plots_with_fresh.pdf", sep='')
  title <- paste0(Title, '_all relative amounts')
  plot <- amounts %>%
    mutate(Name = paste(Name, Sig, sep = ' ')) %>%
    ggplot(., aes(Condition, Av, group = Condition, fill = Condition)) +
    geom_bar(position = "dodge", stat = "identity", colour = "black", width = 0.9) +
    facet_wrap(~ Name, scales = "free", nrow = floor(sqrt(length(unique(amounts$Name))))) +
    theme_bw() +
    labs(x = "", y = "Relative amounts", title = title, fill = element_blank()) +
    theme(
      plot.title = element_text(size = 20, face = "bold", vjust = 2),         #sets title properties
      axis.title = element_text(size = 16, lineheight = 20, face = "bold"),   #sets theme for axis font
      axis.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_blank(),
      legend.text = element_text(face = "bold",size = 12),                  #sets legend text
      strip.text = element_text(face = "bold", size = 15),           #sets theme for title in facets
      panel.grid.major=element_blank()) +
    scale_fill_manual("Conditions", values = colors)  +
    geom_errorbar(aes(ymin=Av, ymax=Av+Std), position=position_dodge(0.9), width=.2)

  pdf(file = plotname, width=14, height=10, pointsize=12)
  print(plot)
  dev.off()

  # Save data
  write.csv(amounts, paste0(Title, "_Amounts", '.csv'), row.names = FALSE)
  write_rds(amounts, file = paste0(Title, "_all_Amounts", '.rds'))

  # Return data
  return(amounts)
}

# Preparing relative amounts relative to fresh medium ---------------------

#' medium_relative subtracts the blank/fresh medium samples from all samples which contain
#' the corresponding medium. The output is a data frame and bar graphs of
#' all metabolites where negative numbers represent consumed metabolites
#' and positive numbers produced metabolites.
#' @author Daniel Braas
#' @param dat A formated data frame produced by prepare_data().
#' @param Title Project title.
#' @param info A data frame with information about the experiment.
#' @return A data frame with integrated metabolomics data as well as
#' bar graphs of all measured metabolites.
#' @export

medium_relative <- function(dat,
                            Title = '',
                            info){

  amounts2 <- dat %>%
    select(Name, KEGG.ID, Condition, contains('Exp')) %>%
    gather(Exp, Amount, -Name, -KEGG.ID, -Condition) %>%
    arrange(Name, Condition) %>%
    unite(Condition_Exp, c(Condition,Exp), sep='_') %>%
    mutate(Condition_Exp = factor(Condition_Exp, levels = unique(Condition_Exp))) %>%
    spread(Condition_Exp, Amount)

  # Find blank medium samples
  fresh <- info[grep('unspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia', info$Condition),]$Sample

  # Remove columns with only NAs; not sure when this became necessary...maybe if there were uneven numbers of samples
  if (length(which(sapply(amounts2[,3:length(amounts2)], sum, na.rm=T) == 0)) > 0){
    amounts2 <- amounts2[-(which(sapply(amounts2[,3:length(amounts2)], sum, na.rm=T) == 0)+2)]
  }

  # Create variable with average amount in blank medium
  for (i in seq_len(length(unique(info$Medium)))) {
    assign(LETTERS[i],
           apply(amounts2[,fresh[c(i*3-2, i*3-1, i*3)]+2], 1, mean, na.rm = T))
  }

  # change NAs in 'fresh' variables to 0
  samples <- info[-grep('unspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia', info$Condition),]
  amounts3 <- amounts2

  # Subtract fresh medium from samples
  for (i in seq_len(length(unique(samples$Medium)))){
    spent <- samples[grep(LETTERS[i], samples$Medium),]$Sample + 2
    amounts3[,spent] <- amounts3[,spent] - get(LETTERS[i])
  }

  # Remove columns with fresh media
  amounts3 <- amounts3[,-(fresh+2)]

  condition <- as.character(samples$Sample.Name)

  amounts3 <- amounts3 %>%
    gather(Condition_Exp, Amount, -Name, -KEGG.ID) %>%
    separate(Condition_Exp, c('Condition', 'Exp'), sep='_') %>%
    mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
    group_by(Name, Condition) %>%
    mutate(Av = mean(Amount, na.rm = T),
           Std = sd(Amount, na.rm = T),
           CV = Std / Av) %>%
    ungroup()

  data8 = split(amounts3, amounts3[,1])
  ANOVA = suppressWarnings(sapply(data8, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
  ANOVA = rep(ANOVA,1,each=length(unique(amounts3$Condition)))

  amounts3 <- amounts3 %>%
    spread(Exp, Amount) %>%
    arrange(Name) %>%
    mutate(ANOVA = ANOVA,
           Sig = '')

  for (i in 1:nrow(amounts3)){
    if (amounts3$ANOVA[i] == "NaN") amounts3$Sig[i] = ""
    else if (amounts3$ANOVA[i] <= 0.001) amounts3$Sig[i] = "***"
    else if (amounts3$ANOVA[i] <= 0.01) amounts3$Sig[i] = "**"
    else if (amounts3$ANOVA[i] <= 0.05) amounts3$Sig[i] = "*"
    else amounts3$Sig[i] = ""
  }
  write.csv(amounts3, paste0(Title,'_Amounts normalized to unspent.csv'), row.names=FALSE)
  write_rds(amounts3, file = paste0(Title, '_Amounts normalized to unspent.rds'))

# Make plot --------------------------------------------------------------

  plotname = paste(Title, "_plots.pdf", sep = '')
  pdf(file = plotname, width=14, height=10, pointsize=12)

  plot <- amounts3 %>%
    mutate(Name = paste(Name, Sig, sep=' '),
           Std = if_else(Av < 0, -1 * Std, Std)) %>% # Make standard deviations negative for negative averages
    ggplot(., aes(Condition, Av, group = Condition, fill = Condition)) +
    geom_bar(position = "dodge", stat = "identity", colour = "black", width = 0.9) +
    facet_wrap(~Name, scales = "free", nrow = floor(sqrt(length(unique(amounts3$Name))))) +
    theme_bw() +
    labs(list(x = "", y = "Relative amounts", title = Title, fill = element_blank())) +
    theme(
      plot.title = element_text(size = 20, face = "bold", vjust = 2),         #sets title properties
      axis.title=element_text(size = 16, lineheight = 20, face = "bold"),   #sets theme for axis font
      axis.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_blank(),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(face = "bold", size = 12),                  #sets legend text
      strip.text = element_text(face = "bold", size = 15),           #sets theme for title in facets
      panel.grid.major = element_blank()) +
    scale_fill_manual("Conditions", values = colors)  +
    geom_errorbar(aes(ymin = Av, ymax = Av + Std), position = position_dodge(0.9), width = .2)
  print(plot)
  dev.off()

  return(amounts3)
}
