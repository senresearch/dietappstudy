# File: demofunc.R
# Author: Courtney Gale
# ----------

demo <- function(demographic, characteristic) {
  num <- sum(demographic == characteristic)
}

percentdemo <- function(demographic, characteristic) {
  demo <- vector(length = 2)
  num <- sum(demographic == characteristic)
  percent <- (num/length(demographic))* 100
  percent <- sprintf("%.1f %%", percent)
  demo <- paste(num,"(", percent, ")", sep = "")
  t(demo)
}

getBMIGroup <- function(initialBMI) {
  if (initialBMI < 18.5) {
    return("Underweight")
  }
  else if (initialBMI > 18.5 & initialBMI < 24.9) {
    return("Normal Weight")
  }
  else if (initialBMI > 25 & initialBMI < 29.9) {
    return("Overweight") 
  }
  else {
    return("Obese")
  }
}

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

demoTest <- function(dat, dfAge){
  # subset number of participants
  dem <- dat[1:9, c(1, 2)]
  
  # convert to numeric type
  dem[, 1] <- as.numeric(dem[, 1])
  dem[, 2] <- as.numeric(dem[, 2])
  
  # apply chi-squared test
  XtestGender <- chisq.test(x = (dem[1:2,c(1,2)]))
  XtestRace <- chisq.test(x = dem[3:5,c(1,2)])
  XtestEducation <- chisq.test(x = dem[6:7,c(1,2)])
  XtestBMI <- chisq.test(x = dem[8:9,c(1,2)])
  
  # apply t-test
  TtestAge <- t.test(Age~Treatment, data = dfAge)
  
  # generate data.frame to return results
  rslts <- data.frame(Xtest = c("XtestGender", "XtestRace", "XtestEducation", "XtestBMI", "TtestAge"),
                      Pvalue = round(c(XtestGender$p.value, XtestRace$p.value, XtestEducation$p.value, 
                                       XtestBMI$p.value, TtestAge$p.value),2))
  
  
  retun(rslts)
  
}
  
  
