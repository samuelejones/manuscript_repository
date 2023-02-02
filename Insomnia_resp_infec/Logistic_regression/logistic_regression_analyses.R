# ANALYSES:
#
# 1) J10_UPPERINFEC ~ F5_INSOMNIA
# 2) J10_ACUTEUPPERINFEC ~ F5_INSOMNIA
# 3) J10_INFLUENZA ~ F5_INSOMNIA
#
# FOLLOWS LEFT-CENSORING TO DATE WHERE ALL REGISTRIES ARE PRESENT
# RIGHT-CENSORING NOW AT 2020-12-31
#
# EXPOSURE->OUTCOME SURVIVAL CONSIDERS ONLY INCIDENT CASES SO ANY PEOPLE
# WITH PREVALENT OUTCOME SHOULD BE EXCLUDED
#
# SAMPLES WITH EXPOSURE BEFORE START DATE ARE ASSUMED TO HAVE THEIR EXPOSURE
# EVENT AT THE START DATE
#
# NEED TO CHECK NUMBER OF CASES AT START OF RECORD
#

# load libraries
library(tidyverse)
library(lubridate)

# set end of records (maximum date in cohort diagnosis data?)
REC_END <- "2020-12-31"

###
#
# INPUT data frame: df_logreg
#
# COLUMNS:
#  - SAMPLEID:                       IDs of cohort samples
#  - SEX_F:                          binary variable; =1 when sample is female, =0 when sample is male
#  - AGE_END_FOLLOWUP:               age in years at REC_END or age at death if earlier than REC_END
#  - BMI:                            BMI of sample
#  - EXPOSURE_BINARY_INSOMNIA:       binary variable; =1 if sample has insomnia dx, =0 if sample has no insomnia dx
#  - OUTCOME_BINARY_ACUTEUPPERINFEC: binary variable; =1 if sample has URI dx, =0 if sample has no URI dx
#  - OUTCOME_BINARY_INFLUENZA:       binary variable; =1 if sample has influenza dx, =0 if sample has no influenza dx
#
###

#
# LOGISTIC REGRESSION MODEL
#
# DESIGN: DEPENDENT   - ACUTEUPPERINFEC, INFLUENZA
#         INDEPENDENT - INSOMNIA
#	  COVARIATES  - SEX, AGE AT FOLLOW-UP END
#

#### RUN MODELS
# 1) ACUTEUPPERINFEC ~ INSOMNIA
# 2) INFLUENZA ~ INSOMNIA

### WITHOUT BMI
# loop over each exposure
for(exposure in c("INSOMNIA")) {

   # sink to log file
   sink(paste0("logistic_regression_results_",exposure,"_exposure.log"))

   # run models
   for(outcome in c("ACUTEUPPERINFEC","INFLUENZA")) {

      cat(paste0("\n\n### Logistic regressions results for ",exposure," indep. var. vs. ",outcome," dep. var. ###\n\n"))

      eval(parse(text=paste0("res.logreg.",exposure,".vs.",outcome,"<- glm(OUTCOME_BINARY_",outcome," ~ EXPOSURE_BINARY_",exposure," + SEX_F + AGE_END_FOLLOWUP, data = df_logreg, family=\"binomial\")")))

      eval(parse(text=paste0("print(summary(res.logreg.",exposure,".vs.",outcome,"))")))

      eval(parse(text=paste0("cis <- confint(res.logreg.",exposure,".vs.",outcome,")")))
      print(cis)

      cat("\n\nP-values:\n")
      eval(parse(text=paste0("pvals <- coef(summary(res.logreg.",exposure,".vs.",outcome,"))[,4]")))
      print(pvals)

   }

   # close sink
   sink()

}

### WITH BMI
# loop over each exposure
for(exposure in c("INSOMNIA")) {

   # sink to log file
   sink(paste0("logistic_regression_results_",exposure,"_exposure_inc.BMI.log"))

   # run models
   for(outcome in c("ACUTEUPPERINFEC","INFLUENZA")) {

      cat(paste0("\n\n### Logistic regressions results for ",exposure," indep. var. vs. ",outcome," dep. var. ###\n\n"))

      eval(parse(text=paste0("res.logreg.",exposure,".vs.",outcome,"<- glm(OUTCOME_BINARY_",outcome," ~ EXPOSURE_BINARY_",exposure," + SEX_F + AGE_END_FOLLOWUP + BMI, data = df_logreg, family=\"binomial\")")))

      eval(parse(text=paste0("print(summary(res.logreg.",exposure,".vs.",outcome,"))")))
      
      eval(parse(text=paste0("cis <- confint(res.logreg.",exposure,".vs.",outcome,")")))
      print(cis)

      cat("\n\nP-values:\n")
      eval(parse(text=paste0("pvals <- coef(summary(res.logreg.",exposure,".vs.",outcome,"))[,4]")))
      print(pvals)

   }

   # close sink
   sink()

}


