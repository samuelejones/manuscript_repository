#!/bin/bash
#

# activate the PYSURV environment - see provided .yml file
conda activate PYSURV

# NEED TO LOOP OVER EACH EXPOSURE-OUTCOME PAIR
# AND CREATE SPECIFIC INPUT_DENSE_FEVENTS FOR EACH
# WITH INDIVIDUALS REMOVED WHO HAVE MISSING DX DATE
# FOR EITHER (OR BOTH) EXPOSURE AND OUTCOME

# set exposure(s)
exposures="insomnia"

# set outcome(s)
outcomes="acuteupperinfec influenza"

# remove combined results file if it exists
rm -f survival_results_all_EXP_all_OUT_incBMI.csv
rm -f survival_results_all_EXP_all_OUT_incBMI.timings

# loop over exposures
for exposure in $exposures
do

   # loop over outcomes
   for outcome in $outcomes
   do

      # extract the relevant input pairs lines
      {
         head -n 1 X_input_pairs.csv
         grep $exposure X_input_pairs.csv | grep $outcome
      } > X_input_pairs_current.csv
   
      # remove individuals with missing age dx for either exposure or outcome
      awk -v expo=$exposure -v outc=$outcome 'BEGIN{FS=","}NR==FNR&&$2=="NA"&&($3==expo||$3==outc){r[$1]++;next}NR==FNR{next}FNR==1||(!($1 in r)){print}' X_input_dense_fevents.csv X_input_dense_fevents.csv > X_input_dense_fevents_current.csv

      # export specific variables for this exposure-outcome pair
      export INPUT_PAIRS="X_input_pairs_current.csv"
      export INPUT_DENSE_FEVENTS="X_input_dense_fevents_current.csv"
      export OUTPUT="survival_results_${exposure}_EXP_${outcome}_OUT_incBMI.csv"
      export TIMINGS="survival_results_${exposure}_EXP_${outcome}_OUT_incBMI.timings"
      # exposure variables which are same for all exposure-outcome pairs
      export INPUT_DEFINITIONS="X_input_definitions.csv"
      export INPUT_INFO="X_input_info_incBMI.csv"

      # set start and end date of survival study for left- and right-censoring
      export STUDY_STARTS="2002.162" # study start date is 2002-03-01 = 2002+(60/365)
      export STUDY_ENDS="2019.631" # study end date is 2019-08-18 = 2019+(230/365)

      # optional - set subcohort size to reduce computing time
      # if set, will only run analysis in a subset of the cohort
      #export N_SUBCOHORT="10_000"

      # remove output files if they exist
      [ -f $OUTPUT ] && rm -v $OUTPUT
      [ -f $TIMINGS ] && rm -v $TIMINGS

      # run python script
      python surv_analysis_incBMI.py

      # combine results (create new files if combined results files don't exist)
      if [ ! -s survival_results_all_EXP_all_OUT_incBMI.csv ]
      then
         cp survival_results_${exposure}_EXP_${outcome}_OUT_incBMI.csv survival_results_all_EXP_all_OUT_incBMI.csv
         cp survival_results_${exposure}_EXP_${outcome}_OUT_incBMI.timings survival_results_all_EXP_all_OUT_incBMI.timings
      else
         tail -n +2 survival_results_${exposure}_EXP_${outcome}_OUT_incBMI.csv >> survival_results_all_EXP_all_OUT_incBMI.csv
         tail -n +2 survival_results_${exposure}_EXP_${outcome}_OUT_incBMI.timings >> survival_results_all_EXP_all_OUT_incBMI.timings
      fi
 
   done

done

# remove exposure-outcome pair specific input files
rm -f X_input_pairs_current.csv X_input_dense_fevents_current.csv

# deactivate the PYSURV environment
conda deactivate
