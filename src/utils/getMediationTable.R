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


        # Regress weight loss at 4 months on mediators with interactions: app usage and number of weigh-in
        outAppInter.fit <- lm(WeightLoss4Month ~ App1 * Treatment + Age, data = dat)
        outNumWInter.fit <- lm(WeightLoss4Month ~  NumWeigh  * Treatment + Age, data = dat)
        

    } else {

        # Regress weight loss at 12 months on treatment
        out.fit <- lm(WeightLoss12Month ~ Treatment + Age, data = dat)

        # Regress Weight loss at 4 months on Treatment with confounders
        med4Month.fit <- lm(WeightLoss4Month ~ Treatment + Age, data = dat)

        # Regress weight loss at 12 months on mediators: app usage, number of weigh-in, and weight loss at 4 months
        outApp.fit <- lm(WeightLoss12Month ~ Treatment + App1 + Age, data = dat)
        outNumW.fit <- lm(WeightLoss12Month ~ Treatment + NumWeigh + Age, data = dat)
        out4Month.fit <- lm(WeightLoss12Month ~ Treatment + WeightLoss4Month + Age, data = dat)

        # Regress weight loss at 12 months on mediators with interactions: app usage, number of weigh-in, and weight loss at         # 4 months
        outAppInter.fit <- lm(WeightLoss12Month ~ App1 * Treatment + Age, data = dat)
        outNumWInter.fit <- lm(WeightLoss12Month ~ NumWeigh  * Treatment + Age, data = dat)
        out4MonthInter.fit <- lm(WeightLoss12Month ~ WeightLoss4Month * Treatment + Age, 
                                 data = dat)

    }

    # Regress App usage mediator on Treatment with confounders
    medApp.fit <- lm(App1 ~ Treatment + Age, data = dat)

    # Regress Number of weigh-in mediator on Treatment with confounders
    medNumW.fit <- lm(NumWeigh ~ Treatment + Age, data = dat)

    ##################################################
    # Models with Treatment and mediator Interaction #
    ##################################################

    # Outcome model with confounders
    if (numFile == 1){

        # Regress weight loss at 4 months on mediators: app usage and number of weigh-in
        outAppInter.fit <- lm(WeightLoss4Month ~ App1 * Treatment + Age, data = dat)
        outNumWInter.fit <- lm(WeightLoss4Month ~  NumWeigh  * Treatment + Age, data = dat)

    } else {

        # Regress weight loss at 12 months on mediators: app usage, number of weigh-in, and weight loss at 4 months
        outAppInter.fit <- lm(WeightLoss12Month ~ App1 * Treatment + Age, data = dat)
        outNumWInter.fit <- lm(WeightLoss12Month ~ NumWeigh  * Treatment + Age, data = dat)
        out4MonthInter.fit <- lm(WeightLoss12Month ~ WeightLoss4Month * Treatment + Age,
                                 data = dat)

    }

    # Save regressions
    if (numFile == 1){
            rgrs.rslts <- list(out.fit = out.fit,                       # 1 Y on X
                               medApp.fit = medApp.fit,                 # 2 M1 on X
                               medNumW.fit = medNumW.fit,               # 3 M2 on X 
                               outApp.fit = outApp.fit,                 # 4 Y on X + M1
                               outNumW.fit = outNumW.fit,               # 5 Y on X + M2
                               outAppInter.fit = outAppInter.fit,       # 6 Y on M1 * X
                               outNumWInter.fit = outNumWInter.fit)     # 7 Y on M2 * X
    } else {
            rgrs.rslts <- list(out.fit = out.fit,                       # 1 Y on X
                               medApp.fit = medApp.fit,                 # 2 M1 on X
                               medNumW.fit = medNumW.fit,               # 3 M2 on X
                               med4Month.fit = med4Month.fit,           # 4 M3 on X
                               outApp.fit = outApp.fit,                 # 5  Y on X + M1
                               outNumW.fit = outNumW.fit,               # 6  Y on X + M2
                               out4Month.fit = out4Month.fit,           # 7  Y on X + M3
                               outAppInter.fit = outAppInter.fit,       # 8  Y on M1 * X
                               outNumWInter.fit = outNumWInter.fit,     # 9  Y on M2 * X
                               out4MonthInter.fit = out4MonthInter.fit) # 10 Y on M3 * X
    }

    ##############
    # Mediations #
    ##############

    if (numFile == 1){

        # Get mediation results with App1 as mediator
        medApp.out <- mediate(rgrs.rslts[[2]], rgrs.rslts[[4]], treat = "Treatment", mediator = "App1",
                              boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)
        # Get mediation results with App1 as mediator assuming interaction with Treatment
        medAppInter.out <- mediate(rgrs.rslts[[2]], rgrs.rslts[[6]], treat = "Treatment", mediator = "App1",
                      boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

        # Get mediation results with NumW as mediator
        medNumW.out <- mediate(rgrs.rslts[[3]], rgrs.rslts[[5]], treat = "Treatment", mediator = "NumWeigh",
                              boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)
        # Get mediation results with NumW as mediator assuming interaction with Treatment
        medNumWInter.out <- mediate(rgrs.rslts[[3]], rgrs.rslts[[7]], treat = "Treatment", mediator = "NumWeigh",
                      boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

    } else {

        # Get mediation results with App1 as mediator
        medApp.out <- mediate(rgrs.rslts[[2]], rgrs.rslts[[5]], treat = "Treatment", mediator = "App1",
                              boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)
        # Get mediation results with App1 as mediator assuming interaction with Treatment
        medAppInter.out <- mediate(rgrs.rslts[[2]], rgrs.rslts[[8]], treat = "Treatment", mediator = "App1",
                                   boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

        # Get mediation results with NumW as mediator
        medNumW.out <- mediate(rgrs.rslts[[3]], rgrs.rslts[[6]], treat = "Treatment", mediator = "NumWeigh",
                               boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)
        # Get mediation results with NumW as mediator assuming interaction with Treatment
        medNumWInter.out <- mediate(rgrs.rslts[[3]], rgrs.rslts[[9]], treat = "Treatment", mediator = "NumWeigh",
                                    boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

        # Get mediation results with NumW as mediator
        med4Month.out <- mediate(rgrs.rslts[[4]], rgrs.rslts[[7]], treat = "Treatment", mediator = "WeightLoss4Month",
                               boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)
        # Get mediation results with NumW as mediator assuming interaction with Treatment
        med4MonthInter.out <- mediate(rgrs.rslts[[4]], rgrs.rslts[[10]], treat = "Treatment", mediator = "WeightLoss4Month",
                                    boot = myBoot, robustSE = myRobustSE, sims = numSims, conf.level = confLvl)

    }
    

    ####################
    # Multi Mediation  #
    ####################

    # Covariates
    myCovar <- c("Age")#, "Gender", "logInitialWeight")
    myMed <- c("App1", "NumWeigh", "WeightLoss4Month")

    # Mediation
    if (numFile == 1){

        # Get multi mediation results with App1 and NumWeigh
        m.medAppNumW <- multimed(outcome = "WeightLoss4Month", med.main = myMed[1], med.alt = myMed[2],
                                 treat = "Treatment", covariates = myCovar, data = dat, sims = numSims, conf.level = confLvl)

    } else {

        # Get multi mediation results with WeightLoss4Month and App1
        m.med4MoApp1 <- multimed(outcome = "WeightLoss12Month", med.main = myMed[3], med.alt = myMed[1],
                                 treat = "Treatment", covariates = myCovar, data = dat, sims = numSims, conf.level = confLvl)

        # Get multi mediation results with App1 and NumWeigh
        m.medAppNumW <- multimed(outcome = "WeightLoss12Month", med.main = myMed[1], med.alt = myMed[2],
                                 treat = "Treatment", covariates = myCovar, data = dat, sims = numSims, conf.level = confLvl)

        # Get multi mediation results with WeightLoss4Month and NumWeigh
        m.med4MoNumW <- multimed(outcome = "WeightLoss12Month", med.main = myMed[3], med.alt = myMed[2],
                                 treat = "Treatment", covariates = myCovar, data = dat, sims = numSims, conf.level = confLvl)

        # Get multi mediation results with WeightLoss4Month, App1 and NumWeigh
        m.med4MoAppNumW <- multimed(outcome = "WeightLoss12Month", med.main = myMed[3], med.alt = myMed[1:2],
                                 treat = "Treatment", covariates = myCovar, data = dat, sims = numSims, conf.level = confLvl)

    }
    
    ###################
    # Save mediations #
    ###################
    
    if (numFile == 1){
             med.rslts <- list(medApp.out = medApp.out,                 # 1 M1 
                               medNumW.out = medNumW.out,               # 2 M2  
                               medAppInter.out = medAppInter.out,       # 3 M1 * X
                               medNumWInter.out = medNumWInter.out,     # 4 M2 * X
                               m.medAppNumW = m.medAppNumW)             # 5 M1 & M2
    } else {
             med.rslts <- list(medApp.out = medApp.out,                 # 1 M1 
                               medNumW.out = medNumW.out,               # 2 M2
                               med4Month.out =med4Month.out,            # 3 M3
                               medAppInter.out = medAppInter.out,       # 4 M1 * X
                               medNumWInter.out = medNumWInter.out,     # 5 M2 * X
                               med4MonthInter.out =med4MonthInter.out,  # 6 M3 * X 
                               m.medAppNumW = m.medAppNumW,             # 7 M1 & M2
                               m.med4MoApp1 = m.med4MoApp1,             # 8 M3 & M1
                               m.med4MoNumW = m.med4MoNumW,             # 9 M3 & M2
                               m.med4MoAppNumW = m.med4MoAppNumW)       # 10 M3 & M1 & M2
    }
    
    
    ########################
    # Sensitivity Analysis #
    ########################

    flag = FALSE
    if (flag){
    if (numFile == 1){
    
        # Get sensitivity analysis results with App1 as mediator
        sensApp.out <- medsens(medApp.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with App1 as mediator assuming interaction with Treatment
        sensAppInter.out <- medsens(medAppInter.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with NumW as mediator
        sensNumW.out <- medsens(medNumW.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with NumW as mediator assuming interaction with Treatment
        sensNumWInter.out <- medsens(medNumWInter.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

    } else {

        # Get sensitivity analysis results with App1 as mediator
        sensApp.out <- medsens(medApp.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with App1 as mediator assuming interaction with Treatment
        sensAppInter.out <- medsens(medAppInter.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with NumW as mediator
        sensNumW.out <- medsens(medNumW.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with NumW as mediator assuming interaction with Treatment
        sensNumWInter.out <- medsens(medNumWInter.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with WeightLoss4Month as mediator
        sens4Month.out <- medsens(med4Month.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

        # Get sensitivity analysis results with WeightLoss4Month as mediator assuming interaction with Treatment
        sens4MonthInter.out <- medsens(med4MonthInter.out, rho.by = 0.1, effect.type = "indirect", sims = numSims)

    }
    }
    ####################
    # Save Sensitivity #
    ####################
    if (flag) { 
    if (numFile == 1){
         sens.rslts <- list(sensApp.out = sensApp.out,                 # 1 M1 
                            sensNumW.out = sensNumW.out,               # 2 M2  
                            sensAppInter.out = sensAppInter.out,       # 3 M1 * X
                            sensNumWInter.out = sensNumWInter.out)     # 4 M2 * X
    } else {
             sens.rslts <- list(sensApp.out = sensApp.out,                 # 1 M1 
                                sensNumW.out = sensNumW.out,               # 2 M2
                                sens4Month.out = sens4Month.out,           # 3 M3
                                sensAppInter.out = sensAppInter.out,       # 4 M1 * X
                                sensNumWInter.out = sensNumWInter.out,     # 5 M2 * X
                                sens4MonthInter.out = sens4MonthInter.out) # 6 M3 * X
    }
    }
    # !!!!!!!!!!!!Need to be added to results!!!!!!!!!!!!!!!!!!!!!!!
    
    ###########################
    # Mediation Results Table #
    ###########################    

    if (numFile == 1){
        # For 4 months weight loss
        k = 5
        # Initialize results table
        med.tbl <- matrix(0, 5, 10)
        for (i in 1:(k-1)){
            med.tbl[i,] <- c(med.rslts[[i]]$tau.coef, as.double(med.rslts[[i]]$tau.ci),
                             med.rslts[[i]]$z.avg, as.double(med.rslts[[i]]$z.avg.ci),
                             med.rslts[[i]]$d.avg, as.double(med.rslts[[i]]$d.avg.ci),
                             med.rslts[[i]]$n.avg) %>% round(6)
        }     

        med.tbl[k,] <- c(med.rslts[[k]]$tau,med.rslts[[k]]$tau.ci,
                         0.5*(med.rslts[[k]]$z0.lb[1]+med.rslts[[k]]$z1.lb[1]), c(as.double(med.rslts[[k]]$z.ave.ci[,1])),
                         0.5*(med.rslts[[k]]$d0.lb[1]+med.rslts[[k]]$d1.lb[1]), c(as.double(med.rslts[[k]]$d.ave.ci[,1])),
                         0.5*(med.rslts[[k]]$d0.lb[1]+med.rslts[[k]]$d1.lb[1])/med.rslts[[k]]$tau) %>% round(6)

    } else { 
        # For 12 months weight loss
        k = 10
        # Initialize results table
        med.tbl <- matrix(0, 10, 10)
        for (i in 1:(k-4)){
            med.tbl[i,] <- c(med.rslts[[i]]$tau.coef, as.double(med.rslts[[i]]$tau.ci),
                             med.rslts[[i]]$z.avg, as.double(med.rslts[[i]]$z.avg.ci),
                             med.rslts[[i]]$d.avg, as.double(med.rslts[[i]]$d.avg.ci),
                             med.rslts[[i]]$n.avg) %>% round(6)
        }     

        for (i in 7:k){
            med.tbl[i,] <- c(med.rslts[[i]]$tau,med.rslts[[i]]$tau.ci,
                            0.5*(med.rslts[[i]]$z0.lb[1]+med.rslts[[i]]$z1.lb[1]), c(as.double(med.rslts[[i]]$z.ave.ci[,1])),
                            0.5*(med.rslts[[i]]$d0.lb[1]+med.rslts[[i]]$d1.lb[1]), c(as.double(med.rslts[[i]]$d.ave.ci[,1])),
                            0.5*(med.rslts[[i]]$d0.lb[1]+med.rslts[[i]]$d1.lb[1])/med.rslts[[i]]$tau) %>% round(6)
        }
    }


    # Give names to the results
    med.tbl <- as.data.frame(med.tbl)
    names(med.tbl)  <- c("Total Effect", "Lower CI", "Upper CI", 
                         "ADE", "Lower CI", "Upper CI",
                         "ACME",  "Lower CI", "Upper CI",
                         "Prop. Mediated")
    if (numFile == 1){
    rownames(med.tbl) <- c("App usage only", "Number Weigh-in only", "App usage with interaction",
                          "Number Weigh-in with interaction", "App usage & Number Weigh-in")

    } else {

        rownames(med.tbl) <- c("App usage only", "Number Weigh-in only", "4 Months Weight Loss",
                              "App usage with interaction", "Number Weigh-in with interaction",
                              "4 Months Weight Loss with interaction", "App usage & Number Weigh-in",
                              "4 Months Weight Loss & App usage", "4 Months Weight Loss & Number Weigh-in",
                              "4 Months Weight Loss & App usage & Number Weigh-in")

    } 
    
    
       
    #############################
    # Sensitivity Results Table #
    #############################

    if (flag){
    
    if (numFile == 1){
        # For 4 months weight loss
        k = 5
        # Initialize results table
        sens.tbl <- matrix(0, 5, 10)
        for (i in 1:(k-1)){
            # Select the rho at which ACME vanishes
            rho_ACME <- cbind(sens.rslts[[i]]$d0,sens.rslts[[i]]$rho) %>% as_tibble() %>% filter(V1>0) %>% filter(V1 == min(V1))
            sens.tbl[i,] <- c(sens.rslts[[i]]$R2tilde.d.thresh, rho_ACME[[2]]) %>% round(6)
        }     
        # Select the R2 at which ACME vanishes
        R2_sigma_ACME <- cbind(med.rslts[[k]]$d.ave.lb, med.rslts[[k]]$R2tilde, med.rslts[[k]]$sigma) %>% as_tibble()%>% filter(V1<0) %>% 
                         abs() %>% filter(V1 == min(V1))
        sens.tbl[k,] <- c(R2_sigma_ACME[[2]], R2_sigma_ACME[[3]]) %>% round(6)

    } else { 
        # For 12 months weight loss
        k = 10
        # Initialize results table
        med.tbl <- matrix(0, 10, 10)
        for (i in 1:(k-4)){
            med.tbl[i,] <- c(med.rslts[[i]]$tau.coef, as.double(med.rslts[[i]]$tau.ci),
                             med.rslts[[i]]$z.avg, as.double(med.rslts[[i]]$z.avg.ci),
                             med.rslts[[i]]$d.avg, as.double(med.rslts[[i]]$d.avg.ci),
                             med.rslts[[i]]$n.avg) %>% round(6)
        }     

        for (i in 7:k){
            med.tbl[i,] <- c(med.rslts[[i]]$tau,med.rslts[[k]]$tau.ci,
                             med.rslts[[i]]$z.ave.lb[1], c(mean(med.rslts[[i]]$z.ave.ci[1,]),mean(med.rslts[[i]]$z.ave.ci[2,])),
                             med.rslts[[i]]$d.ave.lb[1], c(mean(med.rslts[[i]]$d.ave.ci[1,]),mean(med.rslts[[i]]$d.ave.ci[2,])),
                             med.rslts[[i]]$d.ave.lb[1]/med.rslts[[i]]$tau) %>% round(6)
        }
    }


    # Give names to the results
    med.tbl <- as.data.frame(med.tbl)
    names(med.tbl)  <- c("Total Effect", "Lower CI", "Upper CI", 
                         "ADE", "Lower CI", "Upper CI",
                         "ACME",  "Lower CI", "Upper CI",
                         "Prop. Mediated")
    if (numFile == 1){
    rownames(med.tbl) <- c("App usage only", "Number Weigh-in only", "App usage with interaction",
                          "Number Weigh-in with interaction", "App usage & Number Weigh-in")

    } else {

        rownames(med.tbl) <- c("App usage only", "Number Weigh-in only", "4 Months Weight Loss",
                              "App usage with interaction", "Number Weigh-in with interaction",
                              "4 Months Weight Loss with interaction", "App usage & Number Weigh-in",
                              "4 Months Weight Loss & App usage", "4 Months Weight Loss & Number Weigh-in",
                              "4 Months Weight Loss & App usage & Number Weigh-in")

    } 
    
    }    
     list(med.tbl = med.tbl, rgrs.rslts = rgrs.rslts, med.rslts= med.rslts)

}    
