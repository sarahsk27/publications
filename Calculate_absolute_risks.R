##################################
#### Calculate absolute risks ####
##################################

#### Definitions ####

# logRR_G.ERS.PRS = linear predictor calculated using parameter estimates from predictive model
# RR_G.ERS.PRS = exp(logRR_G.ERS.PRS)
# controls = dataset of all controls (non-cases)
# seer_inc_mort = dataset of annual age-specific breast cancer incidence from SEER data in 2008-2012 and age-specific all-cause mortality data from CDC WONDER database in 2008-2012
# Status = case-control status (0: control, 1: case)
# famhxBC = family history in first-degree relative (0: no, 1: yes)
# age = age at reference date
# max age = 80 years in this analysis
# inc = marginal age-specific incidence of breast cancer 
# mort = age-specific all-cause mortality rate


#### Formula to calculate relative risks ####

# Each Beta is a parameter estimate corresponding to a variable in the predictive model

# Example for carriers of a pathogenic variant in a gene other than BRCA1 or BRCA2 (ATM carrier as an example):
# logRR_G.ERS.PRS = Beta_ATM + Beta_logERS * logERS + Beta_PRS * PRS + Beta_PRS.AgeDecade.interaction * PRS * AgeDecade + Beta_PRS.logERS.interaction * PRS * logERS + Beta_famhxBC * famhxBC

# Example for carriers of a pathogenic variant in BRCA1 or BRCA2 (BRCA2 carrier as an example):
# logRR_G.ERS.PRS = Beta_BRCA2 + Beta_logERS * logERS + Beta_PRS * PRS + Beta_BRCA2.AgeDecade.interaction *  AgeDecade + Beta_PRS.BRCA2.interaction * PRS + Beta_PRS.AgeDecade.interaction * PRS * AgeDecade + Beta_PRS.logERS.interaction * PRS * logERS + Beta_famhxBC * famhxBC


#### Formula to calculate absolute risks ####

# LifetimeAbsoluteRisk(G,ERS,PRS) = (annual age-specific incidence/(avgRR_G.ERS.PRS.controls)) * RR_G.ERS.PRS * (1- probability diagnosed at age 0) * (1-probability diagnosed at age 1) * ... * (1-probability diagnosed at max age)


#### Cumulative absolute risks among controls ####

# Calculate average RR among controls 

sumRR_G.ERS.PRS.controls <- sum(controls$RR_G.ERS.PRS)

avgRR_G.ERS.PRS.controls <- sumRR_G.ERS.PRS.controls / length(controls$Status)


# Function to calculate cumulative probability of NOT being diagnosed with breast cancer by age x and surviving to age x

fn_surv.brca <- function(age, inc, mort){
  surv.brca_age <- rep(1, length(age))
  surv.allcause_age <- rep(1, length(age))
  cum_surv <- rep(1, length(age))
  
  for (x in (2:length(age))) {
    surv.brca_age[x] <- 1-inc[x-1]
    surv.allcause_age[x] <- 1-mort[x-1]
    cum_surv[x] <- cum_surv[x-1] * surv.brca_age[x] * surv.allcause_age[x]
  }
  return(cum_surv)
}

seer_inc_mort$cum_surv <- NA
seer_inc_mort$cum_surv <- fn_surv.brca(seer_inc_mort$age, seer_inc_mort$inc_annual, seer_inc_mort$mort_annual)


# Function to calculate lifetime absolute risk (AR)

fn_LR <- function(age, inc, cum_surv, avgRR, RR){
  LR <- rep(0, length(RR))
  AR_age <- rep(0, length(age))
  cum_AR <- rep(0, length(age))
  
  for (i in (1:length(RR))) {
    AR_age[1] <- (inc[1]/avgRR) * RR[i] * cum_surv[1]
    cum_AR[1] <- AR_age[1]
    
    for (x in (2:length(age))) {
      AR_age[x] <- (inc[x]/avgRR) * RR[i] * cum_surv[x]
      cum_AR[x] <- cum_AR[x-1] + AR_age[x]
    }
    
    LR[i] <- cum_AR[length(age)]
  }
  return(LR)
}

controls$LR_G.ERS.PRS <- NA
controls$LR_G.ERS.PRS <- fn_LR(seer_inc_mort$age, seer_inc_mort$inc, seer_inc_mort$cum_surv, avgRR_G.ERS.PRS.controls, controls$RR_G.ERS.PRS)
summary(controls$LR_G.ERS.PRS)
