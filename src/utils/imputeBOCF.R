# File: imputeBOCF.R
# Author: Gregory Farage
# Date: October 2020
# ----------


# Name: imputeBOCF
# Synopsis:
# imputeBOCF imputes NAs with a Baseline Observation Carried Forward.
#
# Args:
# - data             Data frame    Contains input data to be analysed
#
# Returns:
# - data             Data frame    Contains data frame without NAs in weightloss at 4 and 12 months
# 
imputeBOCF <-function(data, response = "both"){ 
    
    # Impute NA in WeightLoss4Month with 0
    # which means no change between 4 month and baseline.
    if(response == "both" || response == "4-month"){
        idxNA4mo <- which(is.na(data$WeightLoss4Month))
        data$WeightLoss4Month[idxNA4mo] <- 0
    }
    # Impute NA in WeightLoss12Month with same value in WeightLoss4Month
    if(response == "both" || response == "12-month"){
        idxNA12mo <- which(is.na(data$WeightLoss12Month))
        data$WeightLoss12Month[idxNA12mo] <- 0
    }
    
    return(data)

}



