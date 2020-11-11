# File: demofunc.R
# Author: Courtney Gale
# Date: October 2020
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
