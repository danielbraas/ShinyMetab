#' This function produces a data frame with relative metabolite amounts.
#' @author Daniel Braas
#' @param DF The input data. This should be a data frame with Name, Iso,
#' Condition, Exp, Value, Used_ID, KEGG.ID, Nr.C, Rt and Formula columns.
#' @param Title The title to be used.
#' @param folder The folder for the data to be saved. Default is the current directory.
#' @return A data frame with relative amount information both averaged as well as averaged and normalized to the first condition.
#' @export

make_RelAmounts <- function(DF, Title = '', folder = getwd()){

  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  if (exists('Title')==F) stop('Title not specified')

  data4 <- DF %>%
    select(Name, Condition, Iso, KEGG.ID, Exp, Value) %>%
    mutate(Exp = factor(Exp, levels = unique(Exp))) %>%
    group_by(Name, Condition, Exp) %>%
    mutate(Amount=sum(Value, na.rm=T)) %>%
    ungroup() %>%
    filter(Iso=='M0') %>%
    select(Name, Condition, KEGG.ID, Exp, Amount) %>%
    spread(Exp, Amount)

  ATP_ADP=try(cbind(Name="ADP/ATP",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                    data4[data4$Name=="ADP",4:length(data4)]/data4[data4$Name=="ATP",4:length(data4)]),
              silent=T)
  if (exists('ATP_ADP')==T & class(ATP_ADP) != 'try-error') data4 <- rbind(data4, ATP_ADP)
  ATP_AMP=try(cbind(Name="AMP/ATP",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                    data4[data4$Name=="AMP",4:length(data4)]/data4[data4$Name=="ATP",4:length(data4)]),
              silent = T)
  if (exists('ATP_AMP')==T & class(ATP_AMP) != 'try-error') data4 <- rbind(data4, ATP_AMP)
  GSH_GSSG=try(cbind(Name="GSH/GSSG",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                     data4[data4$Name=="GSH",4:length(data4)]/data4[data4$Name=="GSSG",4:length(data4)]),
               silent = T)
  if (exists('GSH_GSSG')==T & class(GSH_GSSG) != 'try-error') data4 <- rbind(data4, GSH_GSSG)
  Creatine_PCreatine=try(cbind(Name="Creatine/P-Creatine",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                               data4[data4$Name=="Creatine",4:length(data4)]/data4[data4$Name=="P-Creatine",4:length(data4)]),
                         silent = T)
  if (exists('Creatine_PCreatine')==T & class(Creatine_PCreatine) != 'try-error') data4 <- rbind(data4, Creatine_PCreatine)

  data4 <- data4 %>%
    gather(Exp, Amount, -Name, -Condition,-KEGG.ID)
  data4$Amount <- suppressMessages(mapvalues(data4$Amount, c('Inf','NaN'), c(NA,NA)))
  data4$Amount[data4$Amount==0] <- NA
  test1 <- split(data4, data4[c('Name','Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Amount)))

  data4 <- suppressWarnings(data4 %>%
    arrange(Condition, Name) %>%
    mutate(Amount=new.Value) %>%
    group_by(Name, Condition) %>%
    mutate(Av=mean(Amount, na.rm=T),
           Std=sd(Amount, na.rm=T),
           CV=sd(Amount, na.rm=T)/mean(Amount, na.rm=T)) %>%
    ungroup())

  data8=split(data4, data4[,1])
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
  ANOVA=rep(ANOVA,1,each=length(levels(data4$Condition)))

  data4 <- data4 %>%
    select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
    spread(Exp, Amount) %>%
    arrange(Name, Condition) %>%
    cbind('ANOVA'=ANOVA, 'Sig'=NA) %>%
    group_by(Name) %>%
    mutate(Norm_Av=Av/Av[1],
           Norm_Std=Std/Av[1]) %>%
    ungroup()

  for (i in 1:nrow(data4)){
    if (is.na(data4$ANOVA[i])==T) data4$Sig[i]=""
    else if (data4$ANOVA[i] <= 0.001) data4$Sig[i]="***"
    else if (data4$ANOVA[i] <= 0.01) data4$Sig[i]="**"
    else if (data4$ANOVA[i] <= 0.05) data4$Sig[i]="*"
    else data4$Sig[i]=""
  }
  write_csv(data4, paste0(folder, Title,"_Amounts.csv"))
  write_rds(data4, paste0(folder, 'RelAmounts.rds'))
  return(data4)
}

#' This function produces a data frame with mass isotopologue distribution (MID) data that is corrected for naturally occurring 13C.
#' @author Daniel Braas
#' @param DF The input data. This should be a data frame with Name, Iso, Condition, Exp,
#' Value, Used_ID, KEGG.ID, Nr.C, Rt and Formula columns.
#' @param Title The title to be used.
#' @param folder The folder where the data is supposed to be saved. Default is thecurrent directory.
#' @return A data frame with MID data.
#' @export

make_MID <- function(DF, Title = '', folder = getwd()){

  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  if (exists('Title')==F) stop('Title not specified')

  percent_label <- DF %>%
    as.data.frame() %>%
    select(Name, KEGG.ID, Condition, Iso, Nr.C, Exp, Value) %>%
    mutate(Exp = str_replace(Exp, 'Exp', 'Samp')) %>%
    spread(Exp, Value)
  # INPUT: Max carbons to calculate and #iterations (2 is min number)
  na =.01109
  MaxCarbonCalculate = 50
  MaxIter = 5000

  # Make metabolite name column factors (was not before in test set)
  #percent_label[,1] = as.factor(percent_label[,1])
  percent_label$Name <- factor(percent_label$Name)

  #Duplicate data frame to input corrected values into
  percent_label_corrected=percent_label

  #Make empty data frame for abundance correction parameters (calculated abundance, and %error)
  correctiondata <- data.matrix(matrix(0,nrow(percent_label),16))
  colnames(correctiondata)=c("Metabolite",
                             "rep1 corrected","rep2 corrected","rep3 corrected",
                             "rep1 na added","rep2 na added","rep3 na added",
                             "rep 1 %error","rep 2 %error","rep 3 %error",
                             "rep 1 numiter","rep 2 numiter","rep 3 numiter",
                             "rep 1 errorminfromiter","rep 2 errorminfromiter","rep 3 errorminfromiter")

  correctiondata[,1] <- paste(as.character(percent_label$Name),as.character(percent_label[,3]),as.character(percent_label[,4]))

  #iterate through conditions, assumes conditions is in column 3
  for(c in 1:length(levels(percent_label$Condition))){

    #iterate through metabolites, assumes metabolites are in column 1
    for(n in 1:length(levels(percent_label$Name))){

      #see what metabolite you are currently calculating
      #print(levels(percent_label$Name)[n])

      #find the rows that satisfy condition c and metabolite n
      index <- which((percent_label$Condition == levels(percent_label$Condition)[c]) &
                       (percent_label$Name == levels(percent_label$Name)[n]))

      #Only continue if there are datapoints that satisfy condition c and metabolite n
      if(length(index) != 0){

        #read in the condition + metabolite into matrix temp
        temp <- percent_label[index,]
        #find Cmax for the metabolite
        Cmax <- temp[1,5]

        #limit the max number of carbons to calculate.
        if(Cmax>MaxCarbonCalculate){
          Cmax <- MaxCarbonCalculate
        }

        #Iterate correction across all replicates starting from column 6 where replicate 1 is supposed to be.
        for(z in 6:ncol(temp)){

          #Read natural abundance isotopomers of replicate z from temp and set NA to 0
          Iob=rep(0,Cmax+1)
          Iob[1:nrow(temp)]=temp[1:nrow(temp),z]
          Iob[is.na(Iob)]=0

          #Set observed = Ina(Ina will be Iob updated with Icalc values where Iob=0 later)
          Ina=Iob

          #variables to initalize
          iter=0
          error=1e11
          errornext=1e10

          I=rep(0,Cmax+1)
          Iprev=rep(0,Cmax+1)
          Icalcprev=rep(0,Cmax+1)
          Icalc=rep(0,Cmax+1)
          erroritermin=0

          #Re-update Ina with Iob+Icalc(where Iob=0) as long as error is decreasing (Iob-Icalc)
          while(errornext < error && iter <= MaxIter){

            if(iter >= 2){

              erroritermin = error - errornext + erroritermin

            }

            #save previous rounds values
            error <- errornext
            Iprev <- I
            Icalcprev <- Icalc

            #calculate error minimized from iterations

            #Subtract natural abundance Ina-->I
            for(i in 0:Cmax){

              #initialize toright(carbons that are naturally labeled from the true isotopomer and decrease the measured amount) and fromleft (added
              #from natural abundance labeling of lower isotope, increases the measured amount)

              toright = rep(0,Cmax+1)
              fromleft= rep(0,Cmax+1)

              if((i+1)<=Cmax){
                for(k in (i+1):Cmax){
                  toright[k+1]=(choose(Cmax-i,k-i))*(1-na)^(Cmax-k)*na^(k-i)
                }
              } else {
                toright[k+1]=0
              }
              x=0
              while(x<i){
                fromleft[x+1]=I[x+1]*(choose(Cmax-x,i-x))*(1-na)^(Cmax-i)*na^(i-x)
                x=x+1
              }
              I[i+1]=0
              I[i+1]=(Ina[i+1]-sum(fromleft))/(1-sum(toright))

            }

            #Flatten  + renormalize I to Iob
            I[I<0]=0;
            I=I/(sum(I))*sum(Iob[1:(Cmax+1)]);

            #Add back natural abundance I-->Icalc
            Icalc=rep(0,Cmax+1)
            for(i in 0:Cmax){

              toright = rep(0,Cmax+1)
              fromleft= rep(0,Cmax+1)

              if((i+1)<=Cmax){
                for(k in (i+1):Cmax){
                  toright[k+1]=(choose(Cmax-i,k-i))*(1-na)^(Cmax-k)*na^(k-i)
                }
              } else {
                toright[k+1]=0
              }
              x=0
              while(x<i){
                fromleft[x+1]=I[x+1]*(choose(Cmax-x,i-x))*(1-na)^(Cmax-i)*na^(i-x)
                x=x+1
              }

              Icalc[i+1]=I[i+1]*(1-sum(toright))+sum(fromleft)
            }

            #Calculate error
            errornext=sum(abs(Iob[1:(Cmax+1)]-Icalc))

            #update Ina with Iob and Icalc for 0 values in Iob
            Itemp=Iob
            replacement=which(Iob[1:(Cmax+1)]==0)
            Itemp[replacement]=Icalc[replacement]

            if(is.na(sum(Itemp))){
              Itemp=0
            }

            if(sum(Itemp)!=0){
              Itemp=Itemp/(sum(Itemp))*sum(Iob[1:(Cmax+1)])
            }

            Ina=Itemp
            #print(iter)
            iter=iter+1;

            #if error =NA b/c zero values, then stop iterating
            if(is.na(errornext)){
              errornext=2*error
            }
          }

          #Update values depending on #measurements > or < Cmax
          if(length(index)<Cmax){
            start=index[1]
            end=index[length(index)]
          } else {
            start=index[1]
            end=index[1]+Cmax
          }

          #input error and re-calculated natural abundance MID into correction data
          #correctiondata[start:end,z+3]=as.numeric(error/sum(Iob[1:(Cmax+1)]))*100
          #correctiondata[start:end,z]=Icalcprev[1:(end-start+1)]
          #correctiondata[start:end,z+6]=iter
          #correctiondata[start:end,z+9]=as.numeric(erroritermin/sum(Iob[1:(Cmax+1)]))*100

          #update new values into isotopomer matrix
          percent_label_corrected[start:end,z]=Iprev[1:(end-start+1)]
        }
      }
    }
  }


  #change 0 values to NA
  percent_label_corrected[percent_label_corrected == 0] = NA
  label = percent_label_corrected

  MID <- label %>%
    gather(Exp, Value,-Name,-Condition,-Iso,-KEGG.ID,-Nr.C) %>%
    group_by(Name, Condition, Exp) %>%
    arrange(Name, Condition, Exp) %>%
    mutate(Sum = sum(Value, na.rm = T),
           Fraction = Value * 100 / Sum) %>%
    ungroup() %>%
    group_by(Name, Condition, Iso) %>%
    mutate(
      Norm_Av = mean(Fraction, na.rm = T),
      Norm_Std = sd(Fraction, na.rm = T),
      CV = sd(Fraction, na.rm = T) / mean(Fraction, na.rm = T),
      Av = Norm_Av
    ) %>%
    ungroup() %>%
    select(Name, Condition, Iso, Exp, Fraction, Norm_Av, Norm_Std, CV, Av, Nr.C)
  MID$Exp <- gsub('Samp', 'MID', MID$Exp)

  #This part needs to be inserted at some point to account for NAs
  #test1 <- split(MID, MID[c('Iso','Condition', 'Name')])
  #NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  #new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Fraction)))

  MID$Fraction[is.na(MID$Fraction)] <- 0
  # Anova analysis instead of Ttest
  data8 = split(MID, MID[,c(3,1)], drop=TRUE)
  #ANOVA=sapply(data8, function(x) {if(sum(x$Fraction, na.rm=T)==0) return(NA) else {anova(aov(x$Fraction~x$Condition))$Pr[1]}})
  ANOVA = suppressWarnings(sapply(data8, function(x) anova(aov(x$Fraction~x$Condition))$Pr[1]))

  data3 <- MID %>%
    spread(Exp, Fraction) %>%
    inner_join(label, ., by = c("Name", "Condition", "Iso","Nr.C")) %>%
    arrange(Name, Iso)

  #add indicator of significance
  data3$ANOVA = rep(ANOVA, each = length(unique(data3$Condition)))
  for (i in 1:nrow(data3)){
    if (data3$ANOVA[i] == "NaN") data3$Sig[i] = ""
    else if (data3$ANOVA[i] <= 0.001) data3$Sig[i] = "***"
    else if (data3$ANOVA[i] <= 0.01) data3$Sig[i] = "**"
    else if (data3$ANOVA[i] <= 0.05) data3$Sig[i] = "*"
    else data3$Sig[i] = ""
  }

  write_csv(data3, paste0(folder, Title,"_Mass Isotopologue Distribution.csv"))
  write_rds(data3, paste0(folder, 'MID.rds'))
  return(data3)
}

#' This function produces a data frame with fractional contribution data.
#' @author Daniel Braas
#' @param DF The input data. This should be a data frame that contains MID data generated
#' with the 'make_MID' function.
#' @param Title The title to be used.
#' @param folder The folder where the data is supposed to be saved. Default is the current directory.
#' @return A data frame with fractional contribution data.
#' @export

make_FC <- function(DF, Title = '', folder = getwd()){

  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  if (exists('Title')==F) stop('Title not specified')

  FC <- suppressWarnings(
    DF %>%
      as.data.frame() %>%
      select(Name, KEGG.ID, Condition, Iso, Nr.C, starts_with('MID')) %>%
      gather(Exp, MID, -Name, -KEGG.ID, -Condition, -Iso,-Nr.C) %>%
      mutate(i = as.numeric(gsub('M','',.$Iso)),
             iSi = i * MID) %>%
      group_by(Name, Condition, Exp) %>%
      mutate(FC = sum(iSi, na.rm = T)/Nr.C) %>%
      filter(Iso=='M0') %>%
      ungroup() %>%
      mutate(FC=mapvalues(FC, 0, NA)) %>%                #this is important when a sample is missing or all MID values were 0
      select(Name, KEGG.ID, Condition, Exp, FC) %>%
      group_by(Name, Condition) %>%
      mutate(Norm_Av=mean(FC, na.rm=T),
             Norm_Std=sd(FC, na.rm=T),
             CV=Norm_Std/Norm_Av,
             Av = Norm_Av) %>%
      ungroup()
  )
  FC$Exp <- gsub('MID','FC',FC$Exp)
  test1 <- split(FC, FC[c('Name', 'Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$FC)))

  FC <- FC %>%
    arrange(Condition, Name) %>%
    mutate(FC=new.Value)

  data8=split(FC, FC[,1])
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$FC~x$Condition))$Pr[1]))
  ANOVA=rep(ANOVA,1,each=length(unique(FC$Condition)))

  FC <- spread(FC, Exp, FC) %>%
    arrange(Name) %>%
    mutate(Sig='NA')
  FC$ANOVA <- ANOVA

  for (i in 1:nrow(FC)){
    if (FC$ANOVA[i] == "NaN") FC$Sig[i]=""
    else if (FC$ANOVA[i] <= 0.001) FC$Sig[i]="***"
    else if (FC$ANOVA[i] <= 0.01) FC$Sig[i]="**"
    else if (FC$ANOVA[i] <= 0.05) FC$Sig[i]="*"
    else FC$Sig[i]=""
  }

  write_csv(FC, paste0(folder, Title, '_fractional contribution.csv'))
  write_rds(FC, paste0(folder, 'FC.rds'))
  return(FC)
}

#' This function produces a data frame with percent labeled data.
#' @author Daniel Braas
#' @param DF The input data. This should be a data frame that contains MID data generated
#' with the 'make_MID' function.
#' @param Title The title to be used.
#' @param folder The folder where the data is supposed to be saved. Default is the current directory.
#' @return A data frame with percent labeled data.
#' @export

make_labeled <- function(DF, Title = '', folder = getwd()){

  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  if (exists('Title')==F) stop('Title not specified')

  data_labeled <- DF %>%
    as.data.frame() %>%
    select(Name, KEGG.ID, Condition, Iso, starts_with('MID')) %>%
    filter(Iso=='M0') %>%
    gather(Exp, Value, -Name, -KEGG.ID, -Condition, -Iso) %>%
    mutate(Value=mapvalues(Value, 0, NA),
           Labeled=(1-Value/100)*100) %>%
    group_by(Name, Condition) %>%
    mutate(Norm_Av=mean(Labeled, na.rm=T),
           Norm_Std=sd(Labeled, na.rm=T),
           CV=Norm_Std/Norm_Av,
           Av = Norm_Av) %>%
    select(-Iso, -Value) %>%
    ungroup()
  data_labeled$Exp <- gsub('MID','Labeled', data_labeled$Exp)
  test1 <- split(data_labeled, data_labeled[c('Name', 'Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Labeled)))

  data_labeled <- data_labeled %>%
    arrange(Condition, Name) %>%
    mutate(Labeled=new.Value)

  data8=split(data_labeled, data_labeled[,1])
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Labeled~x$Condition))$Pr[1]))
  ANOVA=rep(ANOVA,1,each=length(unique(data_labeled$Condition)))

  data_labeled <- spread(data_labeled, Exp, Labeled) %>%
    arrange(Name)

  data_labeled <- cbind(data_labeled, ANOVA, 'Sig'=NA)
  for (i in 1:nrow(data_labeled)){
    if (data_labeled$ANOVA[i] == "NaN") data_labeled$Sig[i]=""
    else if (data_labeled$ANOVA[i] <= 0.001) data_labeled$Sig[i]="***"
    else if (data_labeled$ANOVA[i] <= 0.01) data_labeled$Sig[i]="**"
    else if (data_labeled$ANOVA[i] <= 0.05) data_labeled$Sig[i]="*"
    else data_labeled$Sig[i]=""
  }

  write_csv(data_labeled, paste0(folder, Title,'_labeled data.csv'))
  write_rds(data_labeled, paste0(folder, 'percent_labeled.rds'))
  return(data_labeled)
}

#' This function produces a matrix for hierarchical clustering or PCA
#' @author Daniel Braas
#' @param MS_data The data frame to be used for the matrix. Can be relative amounts, fractional
#' contribution or percent labeled.
#' @param Title The title to be used.
#' @param Type Which type of data is used. Right now this is not a used variable, but will become
#' important if I add something for isotopologues.
#' @param anova A cutoff value: usually a p-value of some sort.
#' @param folder The folder where the data is supposed to be saved. Default is the current directory.
#' @return A matrix that can be used as input for hierarchical clustering or PCA
#' @export

make_matrix <- function(MS_data, Title = '', Type, anova = 0.05, folder = getwd()){

  # Check if any ANOVA values are equal or below the cutoff; set cutoff to 1 if all values are below
  if(sum(MS_data$ANOVA <= anova, na.rm = T) == 0) {
    message('All ANOVA values are above the cutoff. Setting ANOVA cutoff to 1!')
    anova <- 1
    Title <- paste0(Title, '_ANOVA_cutoff_1')
  }

  # Find out what kind of input matrix is, ie relative amounts, MIDs etc
  type <- names(MS_data)[str_detect(names(MS_data), 'Exp|MID|FC|Labeled')] %>%
    str_remove('( )?[0-9]+$') %>%
    unique()

  switch(type,
         Exp = ext <- 'Relative Amounts',
         MID = ext <- 'MIDs',
         FC = ext <- 'Fractional Contribution',
         Labeled = ext <- 'Percent Labeled'
  )

  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  if (sum(grepl('MID', names(MS_data))) >= 1) {
    data9 <- MS_data %>%
      filter(ANOVA <= anova) %>%
      select(Name, Condition, Iso, contains('MID')) %>%
      gather(Exp, Value,-Name,-Condition,-Iso) %>%
      arrange(Name, Condition) %>%
      unite(Condition_Exp, c(Condition, Exp), sep = '_') %>%
      mutate(Condition_Exp = factor(Condition_Exp, levels = unique(Condition_Exp))) %>%
      unite(Name_Iso, c(Name, Iso), sep = '_') %>%
      spread(Condition_Exp, Value)

    data9 <- data9[!(rowSums(data9[2:length(data9)], na.rm = T) == 0), ]

  } else {

    data9 <- MS_data %>%
      filter(ANOVA <= anova) %>%
      select(Name, Condition, grep('Exp|FC|Labeled', names(MS_data))) %>%
      gather(Exp, Value,-Name,-Condition) %>%
      arrange(Name, Condition) %>%
      group_by(Name) %>%
      mutate(Nr.NA = sum(is.na(Value)),
             Nr.Samples = n()) %>%
      filter(Nr.NA < Nr.Samples - 1) %>%
      ungroup() %>%
      select(-Nr.NA,-Nr.Samples) %>%
      unite(Condition_Exp, c(Condition, Exp), sep = '_') %>%
      mutate(Condition_Exp = factor(Condition_Exp, levels = unique(Condition_Exp))) %>%
      spread(Condition_Exp, Value)
  }

  data9[is.na(data9)] <- 0

  # Filter out rows and columns without any data
  data9 <- filter(data9, apply(data9[,2:length(data9)], 1, sum) != 0)
  data9 <- select(data9,
                  -(which(sapply(data9[,2:length(data9)], sum) == 0) + 1))

  data5 <- as.matrix(data9[2:length(data9)])
  rownames(data5) <- data9$Name

  # Filter out rows and columns without any data
  #data5 <- data5[!(rowSums(data5)) == 0, ]
  #data5 <- data5[,which(apply(data5, 2, sum) != 0)]

  # Save matrix data
  write.csv(data5, file = paste0(Title, '_Heatmap_', ext, '.csv'))

  return(data5)
}

#' This function produces a data frame with fractional contribution data.
#' @author Daniel Braas
#' @param DF The input data. This should be a data frame that contains the MID envelope.
#' The data frame should have the following columns: Used_ID, Iso and each sample separated
#' into individual columns.
#' @param Title The title to be used.
#' @param max.iso The maximal number of isotopologues in the data (default is 50).
#' @param folder The folder where the data is supposed to be saved. Default is the current directory.
#' @return A data frame with fractional contribution data.
#' @export

make_unlabeled <- function(DF, Title = '', max.iso = 50, folder = getwd()){

  if(!grepl('/$', folder)) folder <- paste0(folder, '/')

  if (exists('Title')==F) stop('Title not specified')

  flat <- function(dat){
    if (nrow(dat)==1) return(dat)
    else {
      for (i in 1:nrow(dat)){
        if (is.na(dat$Value[i])) {
          dat$Value[i:nrow(dat)] <- NA
          return(dat)
        }
      }
    }
  }

  DF$Iso <- factor(DF$Iso, levels=paste0('M', 0:max.iso))
  DF$Value[DF$Value==0] <- NA

  if (sum(grepl('Name', names(DF))) > 0) {
    Order <- levels(DF$Condition)
    DF <- rename(DF, Used_ID=Name) %>%
      as.data.frame() %>%
      unite(Sample, Condition, Exp, sep='_') %>%
      split(.[c('Used_ID','Sample')]) %>%
      map(~ flat(.)) %>%
      do.call(rbind,.) %>%
      separate(Sample, c('Condition','Exp'), sep='_') %>%
      mutate(Condition = factor(Condition, levels = Order)) %>%
      rename(Name=Used_ID)
    rownames(DF) <- NULL

  } else {
    DF <- DF %>%
      as.data.frame() %>%
      gather(Sample, Value,-Used_ID, -Iso) %>%
      arrange(Used_ID, Sample, Iso) %>%
      split(.[c('Used_ID','Sample')]) %>%
      map(~ flat(.)) %>%
      do.call(rbind,.) %>%
      spread(Sample, Value)
  }

  write_csv(DF, paste0(folder, Title, '-Uncorrected raw data.csv'))
  return(DF)
}

