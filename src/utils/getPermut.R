# File: getPermut.R
# Author: Gregory Farage
# Reference: https://statslectures.com/r-scripts-datasets 
#----------


# Name: getPermut
# Synopsis:
# getPermut apply a permutation test based on T-test stats and display the results .
#
# Args:
# - feat             Vector         Contains variable's values to be tested
# - pop              Vector         Contains groupable population variable
# - feat             String         Contains full path of data directory
# - P                Integer        Contains number of permutations
#
# Returns:
# - dfData             Data frame     Contains data frame of lose-it data ready for analysis

getPermut <- function(feat, pop, P = 100000, IVname = "Y"){
    # Set seed for reproducability of results
    set.seed(1234) 
    
    # Set colors for violin plot
    myColors <- c("dodgerblue1","orange2")
    
    # Display violin plot
    vioplot(split(feat ,pop),
            main= IVname, #"Weight Loss",
            xlab="Treatment",ylab= IVname, drawRect=FALSE,
            col = myColors)
    
    # Get reference T-test statistic 
    testStat <- testTstat(feat, pop, T, IVname)
    
    # Number of permutation samples
    cat("Number of permutation sample:", P, "\n")
    
    # Number of observations to sample
    n <- length(pop)  
    
    # Initialize the permutation data matrix
    PermSamples <- matrix(0, nrow=n, ncol=P)
    
    # Generate permutation samples
    for(i in 1:P){
      PermSamples[,i] <- sample(feat, size= n, replace=FALSE)
    }

    # Initialize the Test-stats vector
    PermTestStat <- rep(0, P)
    
    # Estimate T-test statistic for each permutation
    for (i in 1:P){
      PermTestStat[i] <- testTstat(PermSamples[,i], pop)
    }
    
    # Calculate the p-value, for all permuations
    cat("P-value with t-test for", P, "samples is",
        mean(PermTestStat >= testStat), ".\n")
    
    # Generate a density plot of all the permutation 
    # test-stats, and add in our observed test stat
    plot(density(PermTestStat), 
         xlab=expression( group("|", Tc - Tm, "|") ) , 
         main="Permutation Test Stats", las=1)
    abline(v=testStat, col="blue", lty="dotted")
    text(60,0.0005, "p-value", col="blue", cex=0.7)
    
}




# Name: testTstat
# Synopsis:
# testTstat apply a permutation test based on T-test stats and display the results .
#
# Args:
# - feat             Vector         Contains variable's values to be tested
# - pop              Vector         Contains groupable population variable
# - feat             String         Contains full path of data directory
# - rslt             Logical        True returns population's mean and absolute T-test stat
#
# Returns:
# - testRslt          Numeric       Contains T-test stat


testTstat <- function(feat, pop, rslt = F, varName = "var Y"){
    # calculate the difference in sample MEANS
    meanPop1 <- mean(feat[pop==1])    
    meanPop0 <- mean(feat[pop==0])    

    # calculate the absolute diff in T-test stat
    testRslt <- abs((t.test(feat ~ pop))$statistic)
        
    if (rslt) {
        cat("The mean", varName ,"for intervention is:", meanPop1, "\n")
        cat("The mean", varName ,"for no intervention is:", meanPop0, "\n")
        cat("The absolute t-test statistic is: ", testRslt, "\n")
    }
    
    return(testRslt)
    }