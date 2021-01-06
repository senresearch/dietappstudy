# File: analysisFunc.R
# Author: Courtney Gale, Gregory Farage
# Date: October 2020
# ----------




#-------------
#************
# Functions *
#************
#************
weightLossPred <- function(logWeightChange) {
  if (is.na(logWeightChange)) {
    return(NA_character_)
  }
  else if (logWeightChange <= log(.95)) {
    return("Success")
  }
  else if (logWeightChange <= log(.977) & 
           logWeightChange > log(.95)) {
    return("Some Loss")
  }
  else {
    return("No Loss")
  }
}


#************
filt4PCA <- function(myDF, weeks) {
  # Convert categorical variables into numeric.
  # Create a dummy data frame
  if (!is.null(myDF$Treatment)){
    newDF <- dummy_cols(myDF, select_columns = c("Treatment", "Gender"), remove_first_dummy = TRUE)
  }
  else {
    newDF <- myDF
  }
  # Delete unnecessary columns
  newDF$Treatment <- NULL
  newDF$Gender <- NULL
  newDF$WeightLoss12Month <- NULL
  newDF$WeightLoss4Month <- NULL
  newDF$NumWeigh <- NULL
  #newDF <- data.frame(newDF)
  rownames(newDF) <- newDF$id
  newDF$id <- NULL
  newDF$Class4Month <- NULL
  newDF$Class12Month <- NULL
  newDF$X <- NULL
  return(newDF)
}


#-------------
#****************
# MAIN FUNCTION *
#****************

GetLosePCA4 <- function(df.data, numWeeks, numPC, indsup = NULL, graph = FALSE) {
  
  # FILTER DATA FOR PCA
  # Get filtered data
  filtDF <- filt4PCA(df.data, numWeeks)
  # Use the sqrt or the log for results
  filtDF.active <- filtDF
  # START PCA
  # Get the PCA results
  rslt.pca <- PCA(filtDF.active, ncp = numPC, graph = FALSE, ind.sup = indsup) # ncp is number of dimensions kept in the final results 
  #record individual supplements if we have them
  if (is.null(indsup)) {
    individualSup = NULL
  }
  else {
    individualSup = rslt.pca$ind.sup$coord
  }
  # Get the eigenvalues of the PCA
  eig.val <- get_eigenvalue(rslt.pca)
  # Generate the scree plot which plots the eigenvalues ordered from largest to the smallest
  # The number of component could be determined by the point beyond which the remaining eigenvalues are 
  # all relatively small and of comparable size
  scree.plot <- fviz_eig(rslt.pca, addlabels = TRUE, ylim = c(0, 65))
  # Get matrices containing all the results for the active variables, as the correlation of variables
  varPCA <- get_pca_var(rslt.pca)
  if (isTRUE(graph)) { 
    # Display the quality of representation of the variables on factor map
    cos2Corr.plot <- corrplot(varPCA$cos2, is.corr=FALSE, title = "Quality of Representation of Variables")
    # Highlight the most contributing variables for each dimension
    # The larger the value of the contribution, the more the variable contributes to the component
    contribCorr.plot <- corrplot(varPCA$contrib, is.corr=FALSE, title = "Contribution of Var")
  }
  else {
    cos2Corr.plot <- NULL
    contribCorr.plot <- NULL
  }
  # Get matrices containing all the results for the individuals 
  indPCA <- get_pca_ind(rslt.pca)
  indPCA <- data.frame(indPCA$coord)
  #Add Weight Loss and Class Columns into df with individual PC scores
  if (!is.null(indsup)) {
    df.data <- df.data[-indsup,]
  }
  indPCA$WeightLoss4Month <- df.data$WeightLoss4Month
  indPCA$Class4Month <- df.data$Class4Month
  #Filter out NA for the weight loss at 4 or 12 Months
  No.NA <- indPCA
  No.NA <- na.omit(indPCA)
  rslt.list <- list("pca.individual" =  indPCA, "pca.individual.noNA" <- No.NA,
                    "varLoadings" = rslt.pca$var$coord,"indSup" = individualSup,
                    "eigenval" = data.frame(eig.val),
                    "scree" = scree.plot, "QualInterpret" = cos2Corr.plot,
                    "Contrib" = contribCorr.plot)
  return(rslt.list)
}

GetLosePCA12 <- function(df.data, numWeeks, numPC, indsup = NULL, graph = FALSE) {
  
  # FILTER DATA FOR PCA
  # Get filtered data
  filtDF <- filt4PCA(df.data, numWeeks)
  # Use the sqrt or the log for results
  filtDF.active <- filtDF
  # START PCA
  # Get the PCA results
  rslt.pca <- PCA(filtDF.active, ncp = numPC, graph = FALSE,ind.sup = indsup) # ncp is number of dimensions kept in the final results 
  #record individual supplements if we have them
  if (is.null(indsup)) {
    individualSup = NULL
  }
  else {
    individualSup = rslt.pca$ind.sup$coord
  }
  # Get the eigenvalues of the PCA
  eig.val <- get_eigenvalue(rslt.pca)
  # Generate the scree plot which plots the eigenvalues ordered from largest to the smallest
  # The number of component could be determined by the point beyond which the remaining eigenvalues are 
  # all relatively small and of comparable size
  scree.plot <- fviz_eig(rslt.pca, addlabels = TRUE, ylim = c(0, 65))
  # Get matrices containing all the results for the active variables, as the correlation of variables
  varPCA <- get_pca_var(rslt.pca)
  if (isTRUE(graph)) { 
    # Display the quality of representation of the variables on factor map
    cos2Corr.plot <- corrplot(varPCA$cos2, is.corr=FALSE, title = "Quality of Representation of Variables")
    # Highlight the most contributing variables for each dimension
    # The larger the value of the contribution, the more the variable contributes to the component
    contribCorr.plot <- corrplot(varPCA$contrib, is.corr=FALSE, title = "Contribution of Var")
  }
  else {
    cos2Corr.plot <- NULL
    contribCorr.plot <- NULL
  }
  # Get matrices containing all the results for the individuals 
  indPCA <- get_pca_ind(rslt.pca)
  indPCA <- data.frame(indPCA$coord)
  #Add Weight Loss and Class Columns into df with individual PC scores
  if (!is.null(indsup)) {
    df.data <- df.data[-indsup,]
  }
  indPCA$WeightLoss12Month <- df.data$WeightLoss12Month
  indPCA$Class12Month <- df.data$Class12Month
  #Filter out NA for the weight loss at 4 or 12 Months
  No.NA <- indPCA
  No.NA <- na.omit(indPCA)
  rslt.list <- list("pca.individual" =  indPCA, "pca.individual.noNA" <- No.NA,
                    "varLoadings" = rslt.pca$var$coord, "indSup" = individualSup,
                    "eigenval" = data.frame(eig.val),
                    "scree" = scree.plot, "QualInterpret" = cos2Corr.plot,
                    "Contrib" = contribCorr.plot)
}


#**************************************************
crossValNewPCA4M <- function(data, numWeeks, numPC, numOfFolds, model, title,  df.extra = NULL) {
  #find NA response vars
  NAs <- subset(data, is.na(WeightLoss4Month))
  #get k Folds without NAs
  RNGkind(sample.kind = "Rounding")
  set.seed(578)
  data.noNA <- subset(data, !is.na(WeightLoss4Month))
  folds <- createFolds(data.noNA[, "WeightLoss4Month"], k = numOfFolds, returnTrain = TRUE)
  predVec <- vector()
  seedOpts <- c(85,672,827,333,275,875,22,399,222,801)
  for(i in 1:numOfFolds) { 
      #get the index numbers of fold i
      index <- as.numeric(unlist(folds[i]))
      #Conduct a PCA on the fold and data with NA
      test <- data.noNA[-index,]
      testIndex <- which(data$id %in% test$id)
      newPCA <- GetLosePCA4(data, numWeeks, numPC, testIndex)
      newdf <- newPCA$pca.individual
      # Merge the PCA data frame with the traditional variables for training values
      if (!is.null(df.extra)) {
        newDfExtra <- df.extra[-testIndex,]
        newdf$id <- newDfExtra$id
        newdf$Treatment <- newDfExtra$Treatment
        newdf$Gender <- newDfExtra$Gender
        newdf$Age <- newDfExtra$Age
        newdf$logInitialWeight <- newDfExtra$logInitialWeight
        newdf$NumWeigh <- newDfExtra$NumWeigh
        newdf <- subset(newdf, !is.na(WeightLoss4Month))
      }
      #create the first linear model on the fold i
      cv <- lm(model, data = newdf) 
      
      # Merge the PCA data frame with the traditional variables for test values
      if (!is.null(df.extra)) {
        newDfExtra <- df.extra[testIndex,]
        testDF <- as.data.frame(newPCA$indSup)
        testDF$id <- newDfExtra$id
        testDF$Treatment <- newDfExtra$Treatment
        testDF$Gender <- newDfExtra$Gender
        testDF$Age <- newDfExtra$Age
        testDF$logInitialWeight <- newDfExtra$logInitialWeight
        testDF$NumWeigh <- newDfExtra$NumWeigh
        testDF$WeightLoss4Month <- test$WeightLoss4Month
      }
      
      if (!is.null(df.extra)) {
        #predict on test values not in fold i
        preds <- data.frame(predict(cv, testDF))
      }
      else {
        #predict on test values not in fold i
        preds <- data.frame(predict(cv, data.frame(newPCA$indSup)))
      }
      
      #combine predicted values
      predVec <- rbind(predVec, preds) 
  }
  names(predVec) = "Prediction"
  rownames(data.noNA) <- data.noNA$id
  #merge the predicted values into the data set by id number
  newdf <- merge(data.noNA, predVec, by = "row.names") 
  #calculate r-squared 
  rsquared <- round(calcRSqr(newdf$Prediction, newdf[, "WeightLoss4Month"]), digits = 4) 
  spearman <- cor.test(newdf$WeightLoss4Month, newdf$Prediction, method = "spearman", exact = FALSE)
  rho <- round(spearman$estimate, digits = 4)
  #Calculate the AUC for the multi-class prediction
  mauc = GetCstats(newdf$WeightLoss4Month, newdf$Prediction)
  mauc = round(mauc, digits = 4)
  #return a plot of predicted vs actual values
  return(list("rsquared" = rsquared,
              "rho" = unname(rho),
              "mAUC" = mauc))
  # return("plot" = plotPreds(newdf, newdf[, "WeightLoss4Month"], rsquared, rho, mauc, title))
}

#**************************************************
crossValNewPCA12M <- function(data, numWeeks, numPC, numOfFolds, model, title, df.extra = NULL) {
  #find NA response vars
  NAs <- subset(data, is.na(WeightLoss12Month))
  #get k Folds without NAs
  RNGkind(sample.kind = "Rounding")
  set.seed(578)
  data.noNA <- subset(data, !is.na(WeightLoss12Month))
  folds <- createFolds(data.noNA[, "WeightLoss12Month"], k = numOfFolds, returnTrain = TRUE)
  predVec <- vector()
  seedOpts <- c(85,672,827,333,275,875,22,399,222,801)
  for(i in 1:numOfFolds) { 
    #get the index numbers of fold i
    index <- as.numeric(unlist(folds[i]))
    #Conduct a PCA on the fold and data with NA
    test <- data.noNA[-index,]
    testIndex <- which(data$id %in% test$id)
    newPCA <- GetLosePCA12(data, numWeeks, numPC, testIndex)
    newdf <- newPCA$pca.individual
    # Merge the PCA data frame with the traditional variables for training values
    if (!is.null(df.extra)) {
      newDfExtra <- df.extra[-testIndex,]
      newdf$id <- newDfExtra$id
      newdf$Treatment <- newDfExtra$Treatment
      newdf$Gender <- newDfExtra$Gender
      newdf$Age <- newDfExtra$Age
      newdf$logInitialWeight <- newDfExtra$logInitialWeight
      newdf$NumWeigh <- newDfExtra$NumWeigh
      newdf <- subset(newdf, !is.na(WeightLoss12Month))
    }
    #create the first linear model on the fold i
    cv <- lm(model, data = newdf) 
    
    # Merge the PCA data frame with the traditional variables for test values
    if (!is.null(df.extra)) {
      newDfExtra <- df.extra[testIndex,]
      testDF <- as.data.frame(newPCA$indSup)
      testDF$id <- newDfExtra$id
      testDF$Treatment <- newDfExtra$Treatment
      testDF$Gender <- newDfExtra$Gender
      testDF$Age <- newDfExtra$Age
      testDF$logInitialWeight <- newDfExtra$logInitialWeight
      testDF$NumWeigh <- newDfExtra$NumWeigh
      testDF$WeightLoss12Month <- test$WeightLoss12Month
    }
    
    if (!is.null(df.extra)) {
    #predict on test values not in fold i
    preds <- data.frame(predict(cv, testDF))
    }
    else {
      #predict on test values not in fold i
      preds <- data.frame(predict(cv, data.frame(newPCA$indSup)))
    }
    #combine predicted values
    predVec <- rbind(predVec, preds) 
  }
  names(predVec) = "Prediction"
  rownames(data.noNA) <- data.noNA$id
  #merge the predicted values into the data set by id number
  newdf <- merge(data.noNA, predVec, by = "row.names") 
  #calculate r-squared 
  rsquared <- round(calcRSqr(newdf$Prediction, newdf[, "WeightLoss12Month"]), digits = 4)  
  spearman <- cor.test(newdf$WeightLoss12Month, newdf$Prediction, method = "spearman", exact = FALSE)
  rho = round(spearman$estimate, digits = 4)
  #Calculate the AUC for the multi-class prediction
  mauc = GetCstats(newdf$WeightLoss12Month, newdf$Prediction)
  mauc = round(mauc, digits = 4)
  #return a plot of predicted vs actual values
  return(list("rsquared" = rsquared,
              "rho" = unname(rho),
              "mAUC" = mauc))
  # return("plot" = plotPreds(newdf, newdf[, "WeightLoss12Month"], rsquared, rho, mauc, title))
}


#************************************************** 
plotPreds <- function (df, actualWeightLossVar, rsqr, Spearman, mAUC, title) {
  ggplot(data= df, aes(x = Prediction, y = actualWeightLossVar)) + geom_point() +
    scale_y_continuous(breaks = seq(-.26, .18, .04), limits = c(-.26, .18)) +
    scale_x_continuous(breaks = seq(-.08, .06, .02), limits = c(-.08, .065)) +
    labs(y = "Actual", x = "Predicted", title = paste(title, ": " ,"R-squared = ", rsqr, "%    ", "Rho = ", Spearman, "mAUC = ", mAUC)) +
        
      annotate("rect", xmin = -Inf, xmax = log(.95), ymin = -Inf, ymax = log(.95), 
             alpha = .3, fill = "green") +
    annotate("rect", xmin = log(.95), xmax = log(.977), ymin = log(.95), ymax = log(.977), 
             alpha = .3, fill = "blue") +
    annotate("rect", xmin = log(.977), xmax = Inf, ymin = log(.977), ymax = Inf,
             alpha = .3, fill = "red") +
    geom_hline(yintercept = log(.95)) + geom_hline(yintercept = log(.977)) + 
    geom_vline(xintercept = log(.95)) + geom_vline(xintercept = log(.977)) +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
}
#Save as .tiff
#ggsave(paste(title,".tiff", sep = ""), units="in", width = 6, height = 4.5, dpi=600, compression = 'lzw')

#Save as .png
#ggsave(paste(title,".png", sep = ""), units="in", width = 6, height = 4.5, dpi=600)


#************************************************** 
GetCstats <- function (actualVal, predVal) {
   
ActualClass = as.integer(as.logical(actualVal <=log(0.95))) +
                  as.integer(as.logical(actualVal >log(0.95) & 
                                          actualVal <=log(0.977)))*2 +
                  as.integer(as.logical(actualVal > log(0.977)))*3
  
df.class = data.frame(ActualClass)
df.class$PredClass = as.integer(as.logical(predVal <=log(0.95))) +
                        as.integer(as.logical(predVal >log(0.95) & 
                                                predVal <=log(0.977)))*2 +
                        as.integer(as.logical(predVal > log(0.977)))*3

Cstats =  multiclass.roc(response = df.class$Actual, predictor = df.class$PredClass, direction = "<")

return(Cstats$auc[1])

}

#**************************************************
calcRSqr <- function(pred, actual) {
  SSE <- sum((actual - pred)^2) #Calculate error sum of squares
  SSTO <- sum((actual - mean(actual))^2) #Calculate total sum of squares
  rsqr <- 1 - (SSE/SSTO) #r-squared
  return(rsqr)
}

#**************************************************
linearCV <- function(data, varName, numOfFolds, model, title) {
  NAs <- subset(data, is.na(data[, varName]))
  RNGkind(sample.kind = "Rounding")
  set.seed(578)
  data.noNA <- subset(data, !is.na(data[, varName]))
  folds <- createFolds(data.noNA[, varName], k = numOfFolds, returnTrain = TRUE)
  predVec <- vector()
  seedOpts <- c(85,672,827,333,275,875,22,399,222,801)
  for(i in 1:numOfFolds) { 
    #get the index numbers of fold i
    index <- as.numeric(unlist(folds[i]))
    training <- rbind(data.noNA[index,], NAs)
    #create a new linear model based on fold i
    cv <- lm(model, data = training) 
    #predict on test values not in fold i
    preds <- data.frame(predict(cv, data.noNA[-index,])) 
    #combine predicted values
    predVec <- rbind(predVec, preds) 
  }
  
  names(predVec) = "Prediction"
  #merge the predicted values into the data set by id number
  newdf <- merge(data.noNA, predVec, by = 0) 
  #calculate r-squared 
  rsquared <- round(calcRSqr(newdf$Prediction, newdf[, varName]), digits = 4) 
  spearman <- cor.test(newdf[, varName], newdf$Prediction, method = "spearman")
  rho = round(spearman$estimate, digits = 4)
  #Calculate the AUC for the multi-class prediction
  mauc = GetCstats(newdf[, varName], newdf$Prediction)
  mauc = round(mauc, digits = 4)
  #return a plot of predicted vs actual values
  return(list("rsquared" = rsquared,
              "rho" = unname(rho),
              "mAUC" = mauc)) 
  # return(list("rsquared" = rsquared, "plot" = plotPreds(newdf, newdf[, varName], rsquared, rho, mauc,
  #                                                       title)))
}


#**************************************************
allTypes4Month <- function(data, numWeeks, numPC, numFolds, varName) {
  #Create formulas for each group of variables
  tradModel <- as.formula("WeightLoss4Month ~ Age + Treatment")
  
  
  appPCAModel <- as.formula("WeightLoss4Month ~ Dim.1 + Dim.2 + NumWeigh")
  tradAppPCAModel <- as.formula("WeightLoss4Month ~ Dim.1 + Dim.2+
                                  Age + Treatment + NumWeigh")
  
  #Get reduced data frame for calculating PCA with traditional and NumWeigh
  myData <- subset(data, select = c(id, WeightLoss4Month, Treatment, Age, NumWeigh))

  # Get reduced data frame for calculating PCA with only app variables
  # exclude traditional variables and numweigh
  myvars <- names(data) %in% c("WeightLoss12Month", "Treatment", "Age", "Gender", "logInitialWeight", "NumWeigh") 
  appData <- data[!myvars]
  
  NumWeighData <- subset(data, select = c(id, NumWeigh))
  
  #Find actual vs. predicted values 

  # For only baseline variables
  tradNoPCA <- linearCV(data, "WeightLoss4Month", numFolds, tradModel, "No_PCA_Baseline")
  # For only app variables and NumWeigh (i.e. all variables except baseline variables)
  appVarPCA <- crossValNewPCA4M(appData, numWeeks, numPC, numFolds, appPCAModel, "PCAApp_NumWeigh", NumWeighData)
  # For only app variables for PCA then baseline and numWeigh variables are added to the model
  tradAppVarPCA <- crossValNewPCA4M(appData, numWeeks, numPC, numFolds, tradAppPCAModel,
                                    "Baseline_PCAApp_NumWeigh", myData)

  return(list("appPCA" = appVarPCA, "tradAppPCA" = tradAppVarPCA,"tradNoPCA" = tradNoPCA))
  
}

#**************************************************
allTypes12Month <- function(data, numWeeks, numPC, numFolds, varName) {
  
  #Create formulas for each group of variables
  tradModel <- as.formula("WeightLoss12Month ~ Age + Treatment")
  appPCAModel <- as.formula("WeightLoss12Month ~ Dim.1 + Dim.2 + NumWeigh")
  tradAppPCAModel <- as.formula("WeightLoss12Month ~ Dim.1 + Dim.2+
                                  Age + Treatment + NumWeigh")
  
  #Get reduced data frame for calculating PCA with traditional and NumWeigh
  myData <- subset(data, select = c(id, WeightLoss12Month, Treatment, Age, NumWeigh))
  
  # Get reduced data frame for calculating PCA with only app variables
  # exclude traditional variables and NumWeigh
  myvars <- names(data) %in% c("WeightLoss4Month", "Treatment", "Age", "Gender", "logInitialWeight", "NumWeigh") 
  appData <- data[!myvars]
  
  NumWeighData <- subset(data, select = c(id, NumWeigh))
  
  #Find actual vs. predicted values

  # For only baseline variables
  tradNoPCA <- linearCV(data, "WeightLoss12Month", numFolds, tradModel, "No_PCA_Baseline")
  # For only app variables and NumWeigh (i.e. all variables except baseline variables)
  appVarPCA <- crossValNewPCA12M(appData, numWeeks, numPC, numFolds, appPCAModel, "PCAApp_NumWeigh", NumWeighData)
  # For only app variables for PCA then baseline variables and NumWeigh are added to the model
  tradAppVarPCA <- crossValNewPCA12M(appData, numWeeks, numPC, numFolds, tradAppPCAModel, "Baseline_PCAApp_NumWeigh", myData)


  return(list("appPCA" = appVarPCA, "tradAppPCA" = tradAppVarPCA,"tradNoPCA" = tradNoPCA))
}

