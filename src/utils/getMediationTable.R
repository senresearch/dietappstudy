# File: getMediationTable.R
# Author: Gregory Farage
# Date: March 2019
# ----------


# Libraries


# External Functions
# source("imputeBOCF.R")
# source("imputeLOCF.R")

# Name: getMediationTable
# Synopsis:
# getMediationTable returns results from different mediation analysis.
#
# Args:
# - pathFile              String            Contains data directory
# - numFile               Integer           Contains selected data file
#
# Returns:
# - med.tbl               Data Frame        Contains mediation results as a table
# - rgrs.rslts            List              Contains models regression results
# - med.rslts             List              Contains mediation results


# options: 1: "SummaryMonths1to4_2019_02_27.csv"
#          2: "SummaryMonths1to12_2019_02_27.csv"
#          3: "SummaryMonths5to12Weight.csv
#          4: "SummaryMonths5to12All_2019_03_13.csv"

getMediationTable <- function(dat, numFile, myBoot = FALSE, 
                                myRobustSE = FALSE,  numSims = 10000, confLvl = 0.95,
                                imputation = "") {
    
    oldw <- getOption("warn")
    options(warn = -1)


    ##################
    # Data imputation#
    ##################
    
    
    if (imputation == "BOCF"){
        dat <- imputeBOCF(dat)
    } else {
        if (imputation == "LOCF"){
            dat <- imputeLOCF(dat)
        }
    }   


    ##########
    # Models #
    ##########

    # Model with confounders
    if (numFile == 1){

        # Regress weight loss at 4 months on treatment 
        out.fit <- lm(WeightLoss4Month ~ Treatment + Age , data = dat)


        # Regress weight loss at 4 months on mediators: app usage and number of weigh-in
        outApp.fit <- lm(WeightLoss4Month ~ Treatment + App1 + Age , data = dat)
        outNumW.fit <- lm(WeightLoss4Month ~ Treatment + NumWeigh + Age, data = dat)


    } else {

        # Regress weight loss at 12 months on treatment
        out.fit <- lm(WeightLoss12Month ~ Treatment + Age, data = dat)

        # Regress Weight loss at 4 months on Treatment with confounders
        med4Month.fit <- lm(WeightLoss4Month ~ Treatment + Age, data = dat)

        # Regress weight loss at 12 months on mediators: app usage, number of weigh-in, and weight loss at 4 months
        outApp.fit <- lm(WeightLoss12Month ~ Treatment + App1 + Age, data = dat)
        outNumW.fit <- lm(WeightLoss12Month ~ Treatment + NumWeigh + Age, data = dat)
        out4Month.fit <- lm(WeightLoss12Month ~ Treatment + WeightLoss4Month + Age, data = dat)

    }

    # Regress App usage mediator on Treatment with confounders
    medApp.fit <- lm(App1 ~ Treatment + Age, data = dat)

    # Regress Number of weigh-in mediator on Treatment with confounders
    medNumW.fit <- lm(NumWeigh ~ Treatment + Age, data = dat)


    # Save regressions
    if (numFile == 1){
            rgrs.rslts <- list(out.fit = out.fit,                       # 1 Y on X
                               medApp.fit = medApp.fit,                 # 2 M1 on X
                               medNumW.fit = medNumW.fit,               # 3 M2 on X 
                               outApp.fit = outApp.fit,                 # 4 Y on X + M1
                               outNumW.fit = outNumW.fit)               # 5 Y on X + M2
    } else {
            rgrs.rslts <- list(out.fit = out.fit,                       # 1 Y on X
                               medApp.fit = medApp.fit,                 # 2 M1 on X
                               medNumW.fit = medNumW.fit,               # 3 M2 on X
                               med4Month.fit = med4Month.fit,           # 4 M3 on X
                               outApp.fit = outApp.fit,                 # 5  Y on X + M1
                               outNumW.fit = outNumW.fit,               # 6  Y on X + M2
                               out4Month.fit = out4Month.fit)           # 7  Y on X + M3
    }

    ##############
    # Mediations #
    ##############

    if (numFile == 1){

        # Get mediation results with App1 as mediator
        medApp.out <- mediate(rgrs.rslts[[2]], rgrs.rslts[[4]], treat = "Treatment", mediator = "App1",
                              boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

        # Get mediation results with NumW as mediator
        medNumW.out <- mediate(rgrs.rslts[[3]], rgrs.rslts[[5]], treat = "Treatment", mediator = "NumWeigh",
                              boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

    } else {

        # Get mediation results with App1 as mediator
        medApp.out <- mediate(rgrs.rslts[[2]], rgrs.rslts[[5]], treat = "Treatment", mediator = "App1",
                              boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

        # Get mediation results with NumW as mediator
        medNumW.out <- mediate(rgrs.rslts[[3]], rgrs.rslts[[6]], treat = "Treatment", mediator = "NumWeigh",
                               boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

        # Get mediation results with NumW as mediator
        med4Month.out <- mediate(rgrs.rslts[[4]], rgrs.rslts[[7]], treat = "Treatment", mediator = "WeightLoss4Month",
                               boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

    }
    

    
    ###################
    # Save mediations #
    ###################
    
    if (numFile == 1){
             med.rslts <- list(medApp.out = medApp.out,                 # 1 M1 
                               medNumW.out = medNumW.out)               # 2 M2  
            
    } else {
             med.rslts <- list(medApp.out = medApp.out,                 # 1 M1 
                               medNumW.out = medNumW.out,               # 2 M2
                               med4Month.out =med4Month.out)            # 3 M3
    }
    
    
    ########################
    # Sensitivity Analysis #
    ########################

    flag = TRUE
    if (flag){
    if (numFile == 1){
    
        # Get sensitivity analysis results with App1 as mediator
        sensApp.out <- medsens(medApp.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with NumW as mediator
        sensNumW.out <- medsens(medNumW.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

    } else {

        # Get sensitivity analysis results with App1 as mediator
        sensApp.out <- medsens(medApp.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with NumW as mediator
        sensNumW.out <- medsens(medNumW.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with WeightLoss4Month as mediator
        sens4Month.out <- medsens(med4Month.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

    }
    }
    ####################
    # Save Sensitivity #
    ####################
    if (flag) { 
    if (numFile == 1){
         sens.rslts <- list(sensApp.out = sensApp.out,                 # 1 M1 
                            sensNumW.out = sensNumW.out)               # 2 M2  
    } else {
             sens.rslts <- list(sensApp.out = sensApp.out,                 # 1 M1 
                                sensNumW.out = sensNumW.out,               # 2 M2
                                sens4Month.out = sens4Month.out)           # 3 M3
    }
    }
    # !!!!!!!!!!!!Need to be added to results!!!!!!!!!!!!!!!!!!!!!!!
    
    ###########################
    # Mediation Results Table #
    ###########################    

    if (numFile == 1){
        # For 4 months weight loss
        k = 2
        # Initialize results table
        med.tbl <- matrix(0, k, 10)
        for (i in 1:(k)){
            med.tbl[i,] <- c(med.rslts[[i]]$tau.coef, as.double(med.rslts[[i]]$tau.ci),
                             med.rslts[[i]]$z.avg, as.double(med.rslts[[i]]$z.avg.ci),
                             med.rslts[[i]]$d.avg, as.double(med.rslts[[i]]$d.avg.ci),
                             med.rslts[[i]]$n.avg) %>% round(6)
        }     


    } else { 
        # For 12 months weight loss
        k = 3
        # Initialize results table
        med.tbl <- matrix(0, k, 10)
        for (i in 1:(k)){
            med.tbl[i,] <- c(med.rslts[[i]]$tau.coef, as.double(med.rslts[[i]]$tau.ci),
                             med.rslts[[i]]$z.avg, as.double(med.rslts[[i]]$z.avg.ci),
                             med.rslts[[i]]$d.avg, as.double(med.rslts[[i]]$d.avg.ci),
                             med.rslts[[i]]$n.avg) %>% round(6)
        }     

    }


    # Give names to the results
    med.tbl <- as.data.frame(med.tbl)
    names(med.tbl)  <- c("Total Effect", "Lower CI", "Upper CI", 
                         "ADE", "Lower CI", "Upper CI",
                         "ACME",  "Lower CI", "Upper CI",
                         "Prop. Mediated")
    if (numFile == 1){
    rownames(med.tbl) <- c("App-usage", "Self-weight")

    } else {

        rownames(med.tbl) <- c("App-usage", "Self-weight", "4 Months Weight Loss")

    } 
    
    
       
    #############################
    # Sensitivity Results Table #
    #############################

    if (flag){
    
    if (numFile == 1){
        # For 4 months weight loss
        k = 2
        # Initialize results table
        sens.tbl <- matrix(0, k, 10)
        for (i in 1:(k)){
            # Select the rho at which ACME vanishes
            rho_ACME <- sens.rslts[[i]]$rho[which(abs(sens.rslts[[i]]$d0)==min(abs(sens.rslts[[i]]$d0)))]   
            sens.tbl[i,] <- c( rho_ACME, sens.rslts[[i]]$R2tilde.d.thresh) %>% round(6)
        }     
        

    } else { 
        # For 12 months weight loss
        sens.tbl <- matrix(0, k, 10)
        for (i in 1:(k)){
            # Select the rho at which ACME vanishes
            rho_ACME <- sens.rslts[[i]]$rho[which(abs(sens.rslts[[i]]$d0)==min(abs(sens.rslts[[i]]$d0)))]   
            sens.tbl[i,] <- c( rho_ACME, sens.rslts[[i]]$R2tilde.d.thresh) %>% round(6)
        }     
    }


    # Give names to the results
    sens.tbl <- as.data.frame(sens.tbl)
    names(sens.tbl)  <- c("Rho", "VarianceProduct")
    if (numFile == 1){
    rownames(sens.tbl) <- c("App-usage ", "Self-weight")

    } else {

        rownames(sens.tbl) <- c("App-usage ", "Self-weight", "4 Months Weight Loss")

    } 
    
    }  
                                
     options(warn = oldw)
                                
     list(med.tbl = med.tbl, rgrs.rslts = rgrs.rslts, med.rslts= med.rslts, sens.tbl = sens.tbl )

}    
