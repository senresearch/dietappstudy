# File: getData.R
# Author: Gregory Farage
# Date: March 2019
#----------


# Name: getData
# Synopsis:
# getData reads and returns data frame containing lose-it data summary from a csv file.
#
# Args:
# - dirInput             String         Contains full path of data directory
# - i                    Integer        Contains element number of file name list
#
# Returns:
# - dfData             Data frame     Contains data frame of lose-it data ready for analysis
# 
getData <- function(dirInput, i) {
    # File names
    myFileName<-c("SummaryMonths1to4.csv", "SummaryMonths1to12.csv", "SummaryMonths5to12Weight.csv",
                 "SummaryMonths5to12All.csv", "SummaryWeeks1to4_2019.csv", "SummaryWeeks1to8.csv")  
    inputFile <- paste(dirInput, myFileName[i], sep = "")
    
    # Read input file ignoring all single quotes
    dfData <- read.csv(inputFile)
    
    # Remove the first column
    dfData <- dfData[,-1]
    
    # Assign 0 for Self-Paced participants (controlled)
    # Assign 1 for Counselor Initiated participants (treated)
    dfData$Treatment[dfData$Treatment == 2] <- 0
    
    # Change sign for weight loss so that: positive is a loss and negative is a gain.
    dfData$WeightLoss4Month <- -dfData$WeightLoss4Month
    dfData$WeightLoss12Month <- -dfData$WeightLoss12Month
    
    # Change sign for intake calorie so that: intake calories are negative and burnt calories are positive 
    dfData$AvgTotalCal <- -dfData$AvgTotalCal
    dfData$AvgBreakfastCal <- -dfData$AvgBreakfastCal
    dfData$AvgLunchCal <- -dfData$AvgLunchCal
    dfData$AvgDinnerCal <- -dfData$AvgDinnerCal
    dfData$AvgSnackCal <- -dfData$AvgSnackCal
        
    return(dfData)
}
