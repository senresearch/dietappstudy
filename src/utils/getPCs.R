# File: getPCs.R
# Author: Gregory Farage
#----------


# Libraries
# library(FactoMineR)
# library(ggplot2)
# library(factoextra)
# library(corrplot)

# Name: getPCs
# Synopsis:
# getPCs # getPCs returns results from the principal component analysis.
#
# Args:
# - dat                  Data frame            Contains input data to be analysed
# - idxApp               Vector of Integer     Contains selected variable indexes of dat for PCA
# - numPCs               Integer               Contains number of PCs to display in contribution figures 
#
# Returns:
# - dfData             Data frame     Contains data frame of lose-it data ready for analysis
# 
getPCs <- function(dat, idxApp, numPCs) {    
    #**************
    #* Scree plot *
    #**************
    # Get the PCA results
    # ncp is number of dimensions kept in the final results
    # numPCSrslts is the number of PC showing up in the results
    numPCSrslts <-  numPCs
    rsltPCA <- PCA(dat[,idxApp], ncp = numPCSrslts, graph = FALSE, ind.sup = NULL)
    # Set default names for the 5 first PCs
    pcNames <- c("PC 1", "PC 2", "PC 3", "PC 4", "PC 5")
    # Get the eigenvalues of the PCA
    eigVal <- get_eigenvalue(rsltPCA)
    # Generate the scree plot which plots the eigenvalues ordered from largest to the smallest 
    screePlot <- fviz_eig(rsltPCA, geom="line", addlabels = TRUE, ylim = c(0, 65))

    #*********************************
    #* Cumulative Explained Variance *
    #*********************************
    # Generate Cumulative Explained Variances
    cum.var.per <- c(0:10)
    cum.var.per[2:11] <- eigVal[1:10,"cumulative.variance.percent"]
    df.Var <- data.frame(dim = factor(0:10), cum.var.per)

    # Plot setting
    linecolor = "black" 
    main <- "Cumulative Explained Variances"
    xlab <- "Dimensions"
    ylab <- "Cumulative Percentage of explained variances"
    linetype <- "solid"
    eig <- df.Var$cum.var.per
    text_labels <- paste0(round(eig,1), "%")

    # Generate the cumulative variance explained
    p <- ggplot(df.Var, aes(dim, cum.var.per, group =1))+ geom_line(color = linecolor, linetype = linetype)+
    geom_point(shape=19, color=linecolor)+geom_text(label = text_labels, vjust=+2.75, hjust = 0.2)+ 
    labs(title = main, x = xlab, y = ylab)
    
    #****************
    #* Contribution *
    #****************
    # Get correlation of variables and contribution information
    varPCA <- get_pca_var(rsltPCA)
    
    # Display the quality of representation of the variables on factor map
    options(repr.plot.width = 10, repr.plot.height = 5)
    colnames(varPCA$cos2) <- pcNames[1:numPCSrslts]
    
    # Save display plot in within the closure 
    cos2Corr.plot <- function(label.size){
        corrplot(t(varPCA$cos2), method = "circle", is.corr=FALSE, title = "Quality of Representation of Variables",
                              tl.cex = label.size, mar=c(0,0,2,0))
    }
    
    # Highlight the most contributing variables for each dimension
    # The larger the value of the contribution, the more the variable contributes to the component
    colnames(varPCA$contrib) <- pcNames[1:numPCSrslts]
    contribCorr.plot <- function(label.size){
        corrplot(t(varPCA$contrib), is.corr=FALSE, title = "Contribution of Variables",tl.cex = label.size, mar=c(0,0,2,0))
    }   
    # Save display with rate
    cos2CorrNum.plot <- function(){
        corrplot(t(varPCA$cos2), method = "number", is.corr=FALSE, title = "Quality of Representation of Variables",
                              mar=c(0,0,2,0))
    }
    # Highlight the most contributing variables for each dimension
    # The larger the value of the contribution, the more the variable contributes to the component
    contribCorrNum.plot <- function(){
        corrplot(t(varPCA$contrib), method = "number", tl.cex = 1.5,  is.corr=FALSE, title = "Contribution of Variables",mar=c(0,0,2,0))
    }
    
    # Add PC 1 individuals scores to the data
    indPCA <- get_pca_ind(rsltPCA)
    dat$App1 <- indPCA$coord[,1]
    dat$App2 <- indPCA$coord[,2]
    dat$App3 <- indPCA$coord[,3]
    dat$App4 <- indPCA$coord[,3]
    
    list(screePlot = screePlot, CV.plot = p, 
         cos2Corr.plot= cos2Corr.plot, contribCorr.plot = contribCorr.plot,
         cos2CorrNum.plot= cos2CorrNum.plot, contribCorrNum.plot = contribCorrNum.plot,
         data = dat)
}

