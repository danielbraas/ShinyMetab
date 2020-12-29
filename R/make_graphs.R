#' This function produces a bar plot using data from function bar.
#' @author Daniel Braas
#' @param a is the ggplot variable generated with the bar function
#' @param met is the data frame that contains what metabolites should be plotted
#' @param Title is the title for the plot
#' @param x is the title for the x-axis
#' @param y is the title for the y-axis
#' @param axis.text.x tells ggplot whether or not to label x-axis ticks
#' @param scales fixed or free scales
#' @return a data frame of all potential isomers of that particular chemical formula
#' @export

bar_plot = function(a, met, Title, x, y, axis.text.x, scales){
  a + geom_bar(position="dodge", stat="identity", width=0.9) +
    geom_bar(position="dodge", stat="identity", colour="black", width=0.9) +
    facet_wrap( ~ Name, scales=scales) +
    theme_bw() +
    labs(list(x=x, y=y, title=Title, fill=element_blank())) +
    theme(
      plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
      axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
      axis.text=element_text(size=11, face="bold"),
      #axis.text.x=element_blank(),
      axis.text.x=axis.text.x,
      legend.title=element_text(face="bold", size=12),
      legend.text=element_text(face="bold",size=12),                  #sets legend text
      strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
      panel.grid.major=element_blank()) +
    scale_fill_manual("Conditions", values = colors)  +
    geom_errorbar(aes(ymin=Norm_Av, ymax=Norm_Av+Norm_Std), position=position_dodge(0.9), width=.2)+
    geom_text(vjust=-0.1065, color='darkblue', fontface='bold')
}


#' This function produces a bar plot using data from function bar.
#' @author Daniel Braas
#' @param metabolites Metabolic pathway and what type of data to be plotted, for instance: 'glycolysis' for glycolytic metabolites. Individual metabolites can also be entered as vectors.
#' @param bar.type The type of data to be plotted. Can be 'tot' (relative amounts), 'iso' (isotopologue data), 'lab' (percent labeled) or 'FC' (fractioanl contribution)
#' @param repeats (Only needed in previous versions of ggplot2) Number of conditions to be plotted
#' @return a plot initialized with ggplot
#' @import tidyverse
#' @export

bar <- function(metabolites, bar.type, repeats){
  if (length(metabolites) > 1) {
    ending <- 'select metabolites'
  }
  else if (metabolites == 'glycolysis') {
    metabolites = glycolysis
    ending='glycolytic metabolites'
  }
  else if (metabolites=='TCA') {
    metabolites = TCA
    ending <- 'TCA metabolites'
  }
  else if (metabolites=='PPP') {
    metabolites = PPP
    ending <- 'PPP metabolites'
  }
  else if (metabolites=='Curr') {
    metabolites <- Curr
    ending <- 'Currency metabolites'
  }
  else if (metabolites=='Cys'){
    metabolites <- Cys
    ending <- 'Cysteine metabolites'
  }
  else if (metabolites=='Adenine'){
    metabolites <- Adenine
    ending <- 'Adenosine derivatives'
  }
  else if (metabolites=='Cytosine'){
    metabolites <- Cytosine
    ending <- 'Cytidine derivatives'
  }
  else if (metabolites=='Guanine'){
    metabolites <- Guanine
    ending <- 'Guanine derivatives'
  }
  else if (metabolites=='Thymine'){
    metabolites <- Thymine
    ending <- 'Thymine derivatives'
  }
  else if (metabolites=='Uracil'){
    metabolites <- Uracil
    ending <- 'Uracil derivatives'
  }
  else if (metabolites=='AA'){
    metabolites <- AA
    ending <- 'Amino Acids'
  }
  else if (metabolites=='Hex'){
    metabolites <- Hex
    ending <- 'Hexosamine metabolites'
  }
  else if (metabolites=='FA'){
    metabolites <- FA
    ending <- 'Fatty Acids intermediates'
  }
  else if (metabolites=='Fru'){
    metabolites <- Fru
    ending <- 'Fructose Metabolism'
  }
  else if (metabolites=='CoAs'){
    metabolites <- CoAs
    ending <- 'CoA metabolism'
  }
  else if (metabolites=='Neurotrans'){
    metabolites <- Neurotrans
    ending <- 'Neurotransmitter levels'
  }
  else ending = ''

  if (sum(grepl('MID', names(bar.type))) >= 1){
    met = subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Iso=paste(Iso, Sig, sep='\n'),
                   Sig='') %>%
            mutate(Iso = factor(Iso, levels = paste(rep(paste('M', 0:50, sep=''), each=4),
                                              c('','*','**','***'), sep='\n')))

    Title = paste0("Isotopologue distribution of ",ending)
    x <- 'Isotopologue'
    y <- '% Labeled'
    a <-ggplot(met, aes(Iso, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_text(size=11, face="bold")
    bar_plot(a, met, Title, x, y, axis.text.x, scales='free')
  }
  else if (sum(grepl('Exp', names(bar.type))) >= 1){
    met = subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title = paste0("Relative amounts of ",ending)
    x=''
    y='Relative Amounts'
    a <-ggplot(met, aes(Condition, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot(a, met, Title, x, y, axis.text.x, scales='free')
  }

  else if (sum(grepl('Labeled', names(bar.type))) >= 1){
    met <- subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title <- paste0('Percent labeled in ', ending)
    x <- ''
    y <- '% Labeled'
    a <- ggplot(met, aes(Condition, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot(a, met, Title, x, y, axis.text.x, scales='fixed')
  }
  else if (sum(grepl('FC', names(bar.type))) >= 1){
    met <- subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title <- paste0('Fractional Contribution to ', ending)
    x <- ''
    y <- '% Fractional Contribution'
    a <- ggplot(met, aes(Condition, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot(a, met, Title, x, y, axis.text.x, scales='fixed')
  }
}

#' This function takes in a data matrix and produces an annotated heatmap using the pheatmap package.
#' @author Daniel Braas
#' @param matrix The data matrix to be used
#' @param samples A data frame with information about the Conditions and Cell.Number
#' @param Norv The internal standard
#' @param Title The title to be used. This title will also be part of the file name.
#' @param heat.color The color scheme to be used for the heatmap.
#' @param folder The folder where the data is supposed to be saved. Default is the current directory.
#' @return An annotated heatmap saved as a pdf file.
#' @export

make_pheatmap <- function(matrix, samples = info, heat.color = c('navy','grey90','firebrick3'),
                          Norv = 1, Title, folder = getwd()){

  # Make normalizing metabolite = 1 if none is presented
  if(length(Norv) == 1) Norv = rep(1, ncol(matrix))

  # Find out what kind of input matrix is, ie relative amounts, MIDs etc
  type <- str_remove(colnames(matrix), '^(.)*_') %>%
    str_remove('( )?[0-9]+$') %>%
    unique()

  switch(type,
         Exp = ext <- 'Relative Amounts',
         MID = ext <- 'MIDs',
         FC = ext <- 'Fractional Contribution',
         Labeled = ext <- 'Percent Labeled'
         )

  # Make sure the folder has a "/" at the end
  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  #if (exists('samples')==F) samples <- info   # Important for heatmap annotation
  ann <- select(samples, Condition, Cell.Number) %>%
    as.data.frame()
  rownames(ann) <- colnames(matrix)

  ann_colors = list(    # This names list defines the colors for sample groups
    Condition = colors[1:length(unique(gsub('_(.)*','', colnames(matrix))))],
    Norvaline = c("white", "blue"),
    Cell.Number = c("white", "green")
  )

  if (exists('Norv')==T) {    # Show normalization value in each sample if set
    ann$Norvaline <- Norv     # This could become an anker for Shiny
  } else {
    ann_colors$Norvaline <- NULL  # Don't show norvaline trace
  }

  names(ann_colors[['Condition']]) <- unique(gsub('_(.)*','',colnames(matrix)))
  heatmap.title = paste0(folder, Title, '_Heatmap_',ext,'.pdf', sep='')

  if(nrow(matrix) < 2){
    pheatmap::pheatmap(matrix, cluster_row = F, cluster_col = T,
                       color = colorRampPalette(heat.color)(100),
                       border_color = "black", scale = "row",
                       cellwidth = 20, cellheight = 10,
                       annotation = ann, annotation_colors = ann_colors,
                       show_colnames = F, main = paste(Title,ext,sep = '-'),
                       filename = heatmap.title)

  } else{
    pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = T,
                       clustering_distance_rows = 'correlation',
                       clustering_distance_cols = 'correlation',
                       color = colorRampPalette(heat.color)(100),
                       border_color = "black", scale = "row",
                       cellwidth = 20, cellheight = 10,
                       annotation = ann, annotation_colors = ann_colors,
                       show_colnames = F, main = paste(Title,ext,sep = '-'),
                       filename = heatmap.title)
  }

}

#' This function takes in a data matrix and produces a PCA plot and a correlation circle plot.
#' @author Daniel Braas
#' @param matrix The data matrix to be used.
#' @param info A data frame with information about the Conditions and Cell.Number
#' @param a The PC used for the x-axis.
#' @param b The PC used for the y-axis.
#' @param cutoff The cutoff for correlation.
#' @param Title The title to be used. This title will also be part of the file name.
#' @param folder The folder where the data is supposed to be saved. Default is the current directory.
#' @return A pdf file with a scree plot, a pair plot showing the first five PCs as well as a PCA plot with PCs a and b as specified in the function call and the corresponding top 30 loadings as bar plot.
#' @export

make_PCA <- function(matrix,
                     info,
                     a = 1,
                     b = 2,
                     cutoff = 0.5,
                     Title = '',
                     folder = getwd()) {

  if (exists('Title')==F) stop('Title not specified')
  if (!exists(colors)) colors <- c(RColorBrewer::brewer.pal(9, 'Set1'), 'white','black')

  type <- str_remove(colnames(matrix), '^(.)*_') %>%
    str_remove('( )?[0-9]+$') %>%
    unique()

  switch(type,
         Exp = ext <- 'Relative Amounts',
         MID = ext <- 'MIDs',
         FC = ext <- 'Fractional Contribution',
         Labeled = ext <- 'Percent Labeled'
  )

  # Make sure the folder has a "/" at the end
  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  pca <- prcomp(t(matrix), center=T, scale=T)
  var_PCs <- round(summary(pca)$imp[2,]*100,1)
  CCP <- cor(scale(t(matrix), center=T), pca$x, use='pairwise') %>%
    data.frame() %>%
    mutate(Metabolite = rownames(.))

  CCP$Corr <- sqrt((CCP[,a])^2 + (CCP[,b])^2)

  PC <- pca$x %>%
    as.data.frame() %>%
    mutate(Sample.Name = rownames(.),
           Sample.Name = gsub('Exp( )?|MID( )?|FC( )?|Labeled( )?','', Sample.Name)) %>%
    left_join(., info, by='Sample.Name')

   loadings <- data.frame(pca$rotation)
  loadings$Name <- rownames(loadings)
  loadings <- select(loadings, Name, everything())
  write_csv(loadings, paste0(folder, Title,'_Loadings_',ext,'.csv'))

  scores <- data.frame(pca$x)
  write_csv(scores, paste0(folder, Title,'_Scores_', ext, '.csv'))

  PC.title <- paste0(folder, Title,'_PC',a, ' vs. PC', b, '_PCA Plots_', ext, '.pdf')
  pdf(file = PC.title, width=16, height=10)

  # Scree plot
  plot(var_PCs[1:min(10, length(pca$sdev))], type='b', pch=20, col='blue', ylab='Variance explained (%)', xlab='Principal Component', main='Screeplot')

  # Pair plot
  print(lattice::splom(data=PC, ~PC[,1:5],
                       groups = PC$Condition,
                       par.settings=list(superpose.symbol=list(col=colors, pch=19)),
                       auto.key=list(columns=4), pch=19))

  x.lab <- paste0('PC', a, ' (',var_PCs[a],'%)')
  y.lab <- paste0('PC', b, ' (',var_PCs[b],'%)')

  # PCA and Correlation circle plot
  plot <- ggplot(PC, aes(PC[,a], PC[,b], fill=Condition, label=Sample.Name))+
    geom_text(color='black', fontface='bold', size=4, vjust=-0.2)+
    geom_point(size=5, shape=21, color='black')+
    labs(x = x.lab,
         y = y.lab,
         title=paste('PC',a,' vs. PC',b,': All samples', sep=''))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          text = element_text(face='bold'))+
    scale_fill_manual('Condition', values = colors)

  CCP_plot <- filter(CCP, Corr >= cutoff) %>%
    ggplot(., aes(.[,a], .[,b], label=Metabolite))+
    geom_point(size=2, color='grey90')+
    geom_text(vjust=-1, color='navy', size=3)+
    labs(x = x.lab,
         y = y.lab,
         title = 'Correlation circle plot')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          text = element_text(face='bold'))+
    xlim(-1,1)+
    ylim(-1,1)

  gridExtra::grid.arrange(plot, CCP_plot, nrow=1)
  dev.off()

  CCP1 <- suppressWarnings(CCP %>%
                             right_join(., Abbrev, by=c('Metabolite'='Abb')) %>%
                             select(KEGG.ID, paste0('PC',a)))
  names(CCP1)[2] <- 'Norm_Av'
  CCP1$Norm_Av[is.na(CCP1$Norm_Av)] <- 0
  write_csv(CCP1, paste0(folder, 'CCP_PC', a, '_', ext,'.csv'))

  CCP2 <- suppressWarnings(CCP %>%
                             right_join(., Abbrev, by=c('Metabolite'='Abb')) %>%
                             select(KEGG.ID, paste0('PC',b)))
  names(CCP2)[2] <- 'Norm_Av'
  CCP2$Norm_Av[is.na(CCP2$Norm_Av)] <- 0
  write_csv(CCP2, paste0(folder, 'CCP_PC', b, '_', ext,'.csv'))
}
