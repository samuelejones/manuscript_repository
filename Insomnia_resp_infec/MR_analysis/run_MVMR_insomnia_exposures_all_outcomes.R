# load TwoSampleMR library
library(MendelianRandomization)
library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# source file containing winner's curse functions
source("mr_winners_curse_correction.R")


##################################
# load in exposure data manually #
##################################

# read in list of instruments
instr.list <- read.table("mv_mr_varlist_insomniaWatanabe_insomniaLane_smokingWootton_bmiNeale.txt",header=T)
instr.list$CPRA_b37 <- tolower(as.character(instr.list$CPRA_b37))

# emulate the mv_extract_exposures function output
exp_dat_raw <- read.table("mv_mr_summstats_insomniaWatanabe_insomniaLane_smokingWootton_bmiNeale.txt",header=T)

colnames(exp_dat_raw) <- c("exposure","SNP","RSID","CHR","POS","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure")

exp_dat_raw$id.exposure[grepl("Insomnia_",exp_dat_raw$exposure)] <- "Primary"
exp_dat_raw$id.exposure[exp_dat_raw$exposure=="BMI_UKB_Neale_EUR_unrelated"] <- "Secondary1"
exp_dat_raw$id.exposure[exp_dat_raw$exposure=="LifetimeSmokingIndex_Wootton2019"] <- "Secondary2"

exp_dat_raw$exposure <- as.character(exp_dat_raw$exposure)
exp_dat_raw$SNP <- tolower(as.character(exp_dat_raw$SNP))
exp_dat_raw$exposure.instrument <- FALSE
for(exposure in unique(exp_dat_raw$exposure)) {
   exp_dat_raw$exposure.instrument[exp_dat_raw$exposure==exposure&(exp_dat_raw$SNP %in% instr.list$CPRA_b37[instr.list$Phenotype==exposure])] <- TRUE
}
exp_dat_raw$mr_keep <- TRUE

# add ncases, ncontrols and nsamples
exp_dat_raw$ncases.exposure[exp_dat_raw$exposure=="Insomnia_Lane2019"] <- 129270
exp_dat_raw$ncontrols.exposure[exp_dat_raw$exposure=="Insomnia_Lane2019"] <- 108357
exp_dat_raw$nsamples.exposure[exp_dat_raw$exposure=="Insomnia_Lane2019"] <- 129270+108357

exp_dat_raw$ncases.exposure[exp_dat_raw$exposure=="Insomnia_Watanabe2021"] <- 593724
exp_dat_raw$ncontrols.exposure[exp_dat_raw$exposure=="Insomnia_Watanabe2021"] <- 1771286
exp_dat_raw$nsamples.exposure[exp_dat_raw$exposure=="Insomnia_Watanabe2021"] <- 593724+1771286

exp_dat_raw$ncases.exposure[exp_dat_raw$exposure=="LifetimeSmokingIndex_Wootton2019"] <- NA
exp_dat_raw$ncontrols.exposure[exp_dat_raw$exposure=="LifetimeSmokingIndex_Wootton2019"] <- NA
exp_dat_raw$nsamples.exposure[exp_dat_raw$exposure=="LifetimeSmokingIndex_Wootton2019"] <- 462690

exp_dat_raw$ncases.exposure[exp_dat_raw$exposure=="BMI_UKB_Neale_EUR_unrelated"] <- NA
exp_dat_raw$ncontrols.exposure[exp_dat_raw$exposure=="BMI_UKB_Neale_EUR_unrelated"] <- NA
exp_dat_raw$nsamples.exposure[exp_dat_raw$exposure=="BMI_UKB_Neale_EUR_unrelated"] <- 359983

exp_dat_mv <- exp_dat_raw %>% select(-one_of(c("CHR","POS","RSID")))

# make list of variants to extract
expSNP <- unique(exp_dat_mv$SNP)

########################
# read in outcome data #
########################

# COVID A2
out_COVIDA2_dat <- read_outcome_data(
   filename = "20210415_results_20210607_COVID19_HGI_A2_ALL_leave_23andme_20210607.b37_CPRA.txt.gz",
   snps = expSNP,
   snp_col = "rsid",
   sep = "\t",
   beta_col = "all_inv_var_meta_beta",
   se_col = "all_inv_var_meta_sebeta",
   effect_allele_col = "ALT",
   other_allele_col = "REF",
   pval_col = "all_inv_var_meta_p",
   eaf_col = "all_meta_AF"
)
out_COVIDA2_dat$outcome <- "COVIDA2_COVIDHGI"
out_COVIDA2_dat$units.outcome <- "logOR"
out_COVIDA2_dat$type.outcome <- "binary"

# COVID B2
out_COVIDB2_dat <- read_outcome_data(
   filename = "20210415_results_20210607_COVID19_HGI_B2_ALL_leave_ukbb_23andme_20210607.b37_CPRA.txt.gz",
   snps = expSNP,
   snp_col = "rsid",
   sep = "\t",
   beta_col = "all_inv_var_meta_beta",
   se_col = "all_inv_var_meta_sebeta",
   effect_allele_col = "ALT",
   other_allele_col = "REF",
   pval_col = "all_inv_var_meta_p",
   eaf_col = "all_meta_AF"
)
out_COVIDB2_dat$outcome <- "COVIDB2_COVIDHGI"
out_COVIDB2_dat$units.outcome <- "logOR"
out_COVIDB2_dat$type.outcome <- "binary"

# COVID C2
out_COVIDC2_dat <- read_outcome_data(
   filename = "20210415_results_20210607_COVID19_HGI_C2_ALL_leave_ukbb_23andme_20210607.b37_CPRA.txt.gz",
   snps = expSNP,
   snp_col = "rsid",
   sep = "\t",
   beta_col = "all_inv_var_meta_beta",
   se_col = "all_inv_var_meta_sebeta",
   effect_allele_col = "ALT",
   other_allele_col = "REF",
   pval_col = "all_inv_var_meta_p",
   eaf_col = "all_meta_AF"
)
out_COVIDC2_dat$outcome <- "COVIDC2_COVIDHGI"
out_COVIDC2_dat$units.outcome <- "logOR"
out_COVIDC2_dat$type.outcome <- "binary"

# FinnGen R9 URI
out_URI_R9_ALL_dat <- read_outcome_data(
   filename = "FinnGen_R9_summstats/CX_ACUTEUPPERINFEC/finngen_R9_CX_ACUTEUPPERINFEC_b37_CPRA.regenie.gz",
   snps = expSNP,
   snp_col = "ID",
   sep = "\t",
   beta_col = "BETA",
   se_col = "SE",
   effect_allele_col = "ALLELE1",
   other_allele_col = "ALLELE0",
   pval_col = "LOG10P",
   eaf_col = "A1FREQ",
   log_pval=TRUE
)
out_URI_R9_ALL_dat$outcome <- "URI_R9"
out_URI_R9_ALL_dat$units.outcome <- "logOR"
out_URI_R9_ALL_dat$type.outcome <- "binary"

# FinnGen R9 Influenza
out_Influenza_R9_ALL_dat <- read_outcome_data(
   filename = "FinnGen_R9_summstats/CX_INFLUENZA/finngen_R9_CX_INFLUENZA_b37_CPRA.regenie.gz",
   snps = expSNP,
   snp_col = "ID",
   sep = "\t",
   beta_col = "BETA",
   se_col = "SE",
   effect_allele_col = "ALLELE1",
   other_allele_col = "ALLELE0",
   pval_col = "LOG10P",
   eaf_col = "A1FREQ",
   log_pval=TRUE
)
out_Influenza_R9_ALL_dat$outcome <- "Influenza_R9"
out_Influenza_R9_ALL_dat$units.outcome <- "logOR"
out_Influenza_R9_ALL_dat$type.outcome <- "binary"

#################################
# Create directories for output #
#################################

system("mkdir -p MR_plots")
system("mkdir -p MR_results")
system("mkdir -p MR_SNPlists")


#######################################
# Analyse all exposure-outcomes pairs #
#######################################

# create empty data frame for flat results file
mr_results_flat_mv <- data.frame("exposure"=character(),"outcome"=character(),"Nsnp"=numeric(),"Insomnia_IVW_beta"=numeric(),"Insomnia_IVW_SE"=numeric(),"Insomnia_IVW_P"=numeric(),"Insomnia_WM_beta"=numeric(),"Insomnia_WM_SE"=numeric(),"Insomnia_WM_P"=numeric(),"Insomnia_MREgger_beta"=numeric(),"Insomnia_MREgger_SE"=numeric(),"Insomnia_MREgger_P"=numeric(),"Insomnia_R2Sum"=numeric(),"Insomnia_MeanF"=numeric(),"BMI_IVW_beta"=numeric(),"BMI_IVW_SE"=numeric(),"BMI_IVW_P"=numeric(),"BMI_WM_beta"=numeric(),"BMI_WM_SE"=numeric(),"BMI_WM_P"=numeric(),"BMI_MREgger_beta"=numeric(),"BMI_MREgger_SE"=numeric(),"BMI_MREgger_P"=numeric(),"BMI_R2Sum"=numeric(),"BMI_MeanF"=numeric(),"LifetimeSmoking_IVW_beta"=numeric(),"LifetimeSmoking_IVW_SE"=numeric(),"LifetimeSmoking_IVW_P"=numeric(),"LifetimeSmoking_WM_beta"=numeric(),"LifetimeSmoking_WM_SE"=numeric(),"LifetimeSmoking_WM_P"=numeric(),"LifetimeSmoking_MREgger_beta"=numeric(),"LifetimeSmoking_MREgger_SE"=numeric(),"LifetimeSmoking_MREgger_P"=numeric(),"LifetimeSmoking_R2Sum"=numeric(),"LifetimeSmoking_MeanF"=numeric(),"MREggerInt_beta"=numeric(),"MREggerInt_SE"=numeric(),"MREggerInt_P"=numeric())

# loop over primary insomnia exposures
for(ins.exp in c("Insomnia_Lane2019","Insomnia_Watanabe2021")) {
#for(ins.exp in c("Insomnia_Watanabe2021")) {

   ### Insomnia GWAMA (Watanabe et al) ###
   exp_cur_mv <- subset(exp_dat_mv,exposure %in% c(ins.exp,"BMI_UKB_Neale_EUR_unrelated","LifetimeSmokingIndex_Wootton2019"))

   ## create SNP list and add summary stats for each exposure, aligning to lowest alphabetical allele as effect
   snplist <- data.frame(SNP=unique(exp_cur_mv$SNP))
   snplist <- snplist %>% separate(SNP,c("CHR","POS","EFFECT","OTHER"),"_",remove=F)

   for(trait in c(ins.exp,"BMI_UKB_Neale_EUR_unrelated","LifetimeSmokingIndex_Wootton2019"))
   {
      # add BETA, SE, P, PRESENT and EXP_INST columns for each trait
      snplist[[paste0(trait,"__BETA")]] <- NA
      snplist[[paste0(trait,"__SE")]] <- NA
      snplist[[paste0(trait,"__P")]] <- NA
      snplist[[paste0(trait,"__EAF")]] <- NA
      snplist[[paste0(trait,"__R")]] <- NA
      snplist[[paste0(trait,"__R2")]] <- NA
      snplist[[paste0(trait,"__PRESENT")]] <- 0
      snplist[[paste0(trait,"__EXP_INST")]] <- NA
      # take subset of exp_dat_raw for this trait to simplify commands
      exp_cur_mv_trait <- subset(exp_cur_mv,exposure==trait)
      # loop over each variant and
      for(i in 1:length(snplist$SNP)) {
         mi <- match(snplist$SNP[i],exp_cur_mv_trait$SNP)
         # if variant is found for this trait, match won't be NA
         if(!is.na(mi)&!is.na(exp_cur_mv_trait$effect_allele.exposure[mi])) {
            # check if effect allele is correct way round
            if(snplist$EFFECT[i]==tolower(exp_cur_mv_trait$effect_allele.exposure[mi])) {
               # alleles are aligned
               snplist[[paste0(trait,"__BETA")]][i] <- exp_cur_mv_trait$beta.exposure[mi]
               snplist[[paste0(trait,"__SE")]][i] <- exp_cur_mv_trait$se.exposure[mi]
               snplist[[paste0(trait,"__P")]][i] <- exp_cur_mv_trait$pval.exposure[mi]
	       snplist[[paste0(trait,"__EAF")]][i] <- exp_cur_mv_trait$eaf.exposure[mi]
               snplist[[paste0(trait,"__PRESENT")]][i] <- 1
               snplist[[paste0(trait,"__EXP_INST")]][i] <- exp_cur_mv_trait$exposure.instrument[mi]
            } else if(!is.na(exp_cur_mv_trait$other_allele.exposure[mi])&snplist$EFFECT[i]==tolower(exp_cur_mv_trait$other_allele.exposure[mi])) {
               # alleles are not aligned
               snplist[[paste0(trait,"__BETA")]][i] <- -exp_cur_mv_trait$beta.exposure[mi]
               snplist[[paste0(trait,"__SE")]][i] <- exp_cur_mv_trait$se.exposure[mi]
               snplist[[paste0(trait,"__P")]][i] <- exp_cur_mv_trait$pval.exposure[mi]
	       snplist[[paste0(trait,"__EAF")]][i] <- 1-exp_cur_mv_trait$eaf.exposure[mi]
               snplist[[paste0(trait,"__PRESENT")]][i] <- 1
               snplist[[paste0(trait,"__EXP_INST")]][i] <- exp_cur_mv_trait$exposure.instrument[mi]
            } else if (snplist$OTHER[i]==tolower(exp_cur_mv_trait$effect_allele.exposure[mi])) {
               snplist[[paste0(trait,"__BETA")]][i] <- -exp_cur_mv_trait$beta.exposure[mi]
               snplist[[paste0(trait,"__SE")]][i] <- exp_cur_mv_trait$se.exposure[mi]
               snplist[[paste0(trait,"__P")]][i] <- exp_cur_mv_trait$pval.exposure[mi]
	       snplist[[paste0(trait,"__EAF")]][i] <- 1-exp_cur_mv_trait$eaf.exposure[mi]
               snplist[[paste0(trait,"__PRESENT")]][i] <- 1
               snplist[[paste0(trait,"__EXP_INST")]][i] <- exp_cur_mv_trait$exposure.instrument[mi]
            }
         }
      }
      ## calculate r and r2 for trait (then drop r and eaf columns)
      # use "get_r_from_lor" if binary or "get_r_from_bsen" if continuous
      if(sum(is.na(exp_cur_mv_trait$ncases.exposure))==0) {
	 snplist[[paste0(trait,"__R")]][!is.na(snplist[[paste0(trait,"__BETA")]])] <- as.numeric(get_r_from_lor(lor=snplist[[paste0(trait,"__BETA")]][!is.na(snplist[[paste0(trait,"__BETA")]])],af=snplist[[paste0(trait,"__EAF")]][!is.na(snplist[[paste0(trait,"__BETA")]])],ncase=rep(unique(exp_cur_mv_trait$ncases.exposure),length(snplist[[paste0(trait,"__EAF")]][!is.na(snplist[[paste0(trait,"__BETA")]])])),ncontrol=rep(unique(exp_cur_mv_trait$ncontrols.exposure),length(snplist[[paste0(trait,"__EAF")]][!is.na(snplist[[paste0(trait,"__BETA")]])])),prevalence=rep(unique(exp_cur_mv_trait$ncases.exposure)/(unique(exp_cur_mv_trait$ncases.exposure)+unique(exp_cur_mv_trait$ncontrols.exposure)),length(snplist[[paste0(trait,"__EAF")]][!is.na(snplist[[paste0(trait,"__BETA")]])])),model="logit",correction=F))
         #snplist[[paste0(trait,"__R")]] <- as.numeric(get_r_from_lor(lor=snplist[[paste0(trait,"__BETA")]],af=snplist[[paste0(trait,"__EAF")]],ncase=rep(unique(exp_cur_mv_trait$ncases.exposure),length(snplist[[paste0(trait,"__EAF")]])),ncontrol=rep(unique(exp_cur_mv_trait$ncontrols.exposure),length(snplist[[paste0(trait,"__EAF")]])),prevalence=rep(unique(exp_cur_mv_trait$ncases.exposure)/(unique(exp_cur_mv_trait$ncases.exposure)+unique(exp_cur_mv_trait$ncontrols.exposure)),length(snplist[[paste0(trait,"__EAF")]])),model="logit",correction=F))
      } else {
         snplist[[paste0(trait,"__R")]] <- as.numeric(get_r_from_bsen(b=snplist[[paste0(trait,"__BETA")]],se=snplist[[paste0(trait,"__SE")]],n=rep(unique(exp_cur_mv_trait$nsamples.exposure),length(snplist[[paste0(trait,"__BETA")]]))))
      }
      # calculate r2 column
      snplist[[paste0(trait,"__R2")]] <- snplist[[paste0(trait,"__R")]]^2
      # drop r and eaf columns
      snplist <- snplist %>% select(-one_of(c(paste0(trait,"__R"),paste0(trait,"__EAF"))))
   }


   # make list of variants not present in all exposures
   uniq_vars_all_exp <- unlist(dimnames(table(exp_cur_mv$SNP)[table(exp_cur_mv$SNP)==length(unique(exp_cur_mv$exposure))]))
   attr(uniq_vars_all_exp, "names") <- NULL
   # drop those that have NA for any exposure
   uniq_vars_all_exp_no_NA <- uniq_vars_all_exp[! uniq_vars_all_exp %in% unique(subset(exp_cur_mv,is.na(eaf.exposure)|is.na(beta.exposure)|is.na(se.exposure)|is.na(pval.exposure))$SNP)]

   # subset current exposure stats for those variants not missing data
   exp_cur_mv <- subset(exp_cur_mv,SNP %in% uniq_vars_all_exp_no_NA)

   # flip betas, effect alleles and eaf if beta<0
   exp_cur_mv$ean  <- as.character(exp_cur_mv$effect_allele.exposure)
   exp_cur_mv$oan  <- as.character(exp_cur_mv$other_allele.exposure)
   exp_cur_mv$eafn <- exp_cur_mv$eaf.exposure
   exp_cur_mv$bn   <- exp_cur_mv$beta.exposure

   exp_cur_mv$ean[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8]  <- as.character(exp_cur_mv$other_allele.exposure[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8])
   exp_cur_mv$oan[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8]  <- as.character(exp_cur_mv$effect_allele.exposure[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8])
   exp_cur_mv$eafn[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8] <- 1-exp_cur_mv$eaf.exposure[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8]
   exp_cur_mv$bn[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8]   <- -exp_cur_mv$beta.exposure[exp_cur_mv$beta.exposure<0&exp_cur_mv$pval.exposure<=5e-8]

   exp_cur_mv$effect_allele.exposure <- exp_cur_mv$ean
   exp_cur_mv$other_allele.exposure  <- exp_cur_mv$oan
   exp_cur_mv$eaf.exposure           <- exp_cur_mv$eafn
   exp_cur_mv$beta.exposure          <- exp_cur_mv$bn

   exp_cur_mv[,c("ean","oan","eafn","bn")] <- list(NULL)

   # if exposure is binary, calculate r^2 from logOR, otherwise use beta, SE and sample size
   exp_cur_mv$r.exposure <- NA
   exp_cur_mv$r.exposure[!is.na(exp_cur_mv$ncases.exposure)] <- as.numeric(get_r_from_lor(exp_cur_mv$beta.exposure[!is.na(exp_cur_mv$ncases.exposure)],exp_cur_mv$eaf.exposure[!is.na(exp_cur_mv$ncases.exposure)],exp_cur_mv$ncases.exposure[!is.na(exp_cur_mv$ncases.exposure)],exp_cur_mv$ncontrols.exposure[!is.na(exp_cur_mv$ncases.exposure)],exp_cur_mv$ncases.exposure[!is.na(exp_cur_mv$ncases.exposure)]/(exp_cur_mv$ncases.exposure[!is.na(exp_cur_mv$ncases.exposure)]+exp_cur_mv$ncontrols.exposure[!is.na(exp_cur_mv$ncases.exposure)]),model="logit",correction=F))
   exp_cur_mv$r.exposure[is.na(exp_cur_mv$ncases.exposure)] <- as.numeric(get_r_from_bsen(b=exp_cur_mv$beta.exposure[is.na(exp_cur_mv$ncases.exposure)],se=exp_cur_mv$se.exposure[is.na(exp_cur_mv$ncases.exposure)],n=exp_cur_mv$nsamples.exposure[is.na(exp_cur_mv$ncases.exposure)]))

   exp_cur_mv$r2.exposure <- exp_cur_mv$r.exposure^2

   ### loop over all outcomes
   for (outcome in c("COVIDA2","COVIDB2","COVIDC2","URI_R9","Influenza_R9")) {

      # harmonise data (use eval as data frame name changes)
      harmonised_data <- eval(parse(text=paste0("mv_harmonise_data(exp_cur_mv,out_",outcome,"_dat)")))

      # add empty variables to store r2
      harmonised_data$exposure_r2 <- harmonised_data$exposure_beta
      harmonised_data$exposure_r2[,] <- NA

      # ensure that expname variable in harmonised_data data frame is in the same order
      # as the other variables
      if(sum(harmonised_data$expname$id.exposure==unlist(dimnames(harmonised_data$exposure_beta)[2]))!=length(harmonised_data$expname$id.exposure)) {
         exp_reorder <- match(harmonised_data$expname$id.exposure,unlist(dimnames(harmonised_data$exposure_beta)[2]))
	 harmonised_data$expname$id.exposure <- harmonised_data$expname$id.exposure[exp_reorder]
	 harmonised_data$expname$exposure <- harmonised_data$expname$exposure[exp_reorder]
      }

      # perform winner's curse correction for each exposure individually
      for(id.exp in harmonised_data$expname$id.exposure) {

         # get exp name for this id.exp
         expname <- harmonised_data$expname$exposure[harmonised_data$expname$id.exposure==id.exp]
      
         # pull out exposure and outcome betas only for variants that are instruments
         # for this exposure
         bxg <- as.numeric(harmonised_data$exposure_beta[(dimnames(harmonised_data$exposure_beta)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]),c(id.exp)])
         sebxg <- as.numeric(harmonised_data$exposure_se[(dimnames(harmonised_data$exposure_se)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]),c(id.exp)])
         byg <- harmonised_data$outcome_beta[dimnames(harmonised_data$exposure_beta)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]]
         sebyg <- harmonised_data$outcome_se[dimnames(harmonised_data$exposure_beta)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]]

         wcc_stats <- wcc2(BetaXG=bxg,seBetaXG=sebxg,BetaYG=byg,seBetaYG=sebyg)

         # store wcc stats back in separate variables
         bxg <- wcc_stats$BetaXGc
         sebxg <- wcc_stats$seBetaXGc
         byg <- wcc_stats$BetaYGc
         sebyg <- wcc_stats$seBetaYGc

         # map these wcc corrected stats back into correct location of harmonised_data
         # data frame
         harmonised_data$exposure_beta[(dimnames(harmonised_data$exposure_beta)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]),c(id.exp)] <- bxg
         harmonised_data$exposure_se[(dimnames(harmonised_data$exposure_se)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]),c(id.exp)] <- sebxg
         harmonised_data$outcome_beta[dimnames(harmonised_data$exposure_beta)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]] <- byg
         harmonised_data$outcome_se[dimnames(harmonised_data$exposure_beta)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]] <- sebyg

         # add r2 for variants that are instruments for each exposure
         r2.tmp <- subset(exp_cur_mv[c("exposure","SNP","r2.exposure")],exposure==expname)
         instr.list.exp.cur <- data.frame("SNP"=dimnames(harmonised_data$exposure_r2)[[1]][dimnames(harmonised_data$exposure_r2)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]])
         r2.instr.list.cur <- merge(r2.tmp,instr.list.exp.cur,by="SNP",all.y=TRUE)
         harmonised_data$exposure_r2[(dimnames(harmonised_data$exposure_r2)[[1]] %in% instr.list$CPRA_b37[instr.list$Phenotype==expname]),c(id.exp)] <- r2.instr.list.cur$r2.exposure

      }

      # add column to snplist to show whether variant is used for this outcome
      snplist[[paste0(outcome,"__USED")]] <- 0
      snplist[[paste0(outcome,"__USED")]][snplist$SNP %in% unlist(dimnames(harmonised_data$exposure_beta)[1])] <- 1 

      # run MV IVW in TwoSampleMR
      #mr_results_2smr_ivw_cur <- mv_multiple(
      #   harmonised_data,
      #   intercept = FALSE,
      #   instrument_specific = TRUE,
      #   pval_threshold = 5e-08,
      #   plots = FALSE
      #)

      # create MRMVInputObject for Mendelian Randomization package functions
      MRMVInputObject_cur <- mr_mvinput(bx=cbind(harmonised_data$exposure_beta[,c("Primary")],harmonised_data$exposure_beta[,c("Secondary1")],harmonised_data$exposure_beta[,c("Secondary2")]),bxse=cbind(harmonised_data$exposure_se[,c("Primary")],harmonised_data$exposure_se[,c("Secondary1")],harmonised_data$exposure_se[,c("Secondary2")]),by=harmonised_data$outcome_beta,byse=harmonised_data$outcome_se)

      mr_results_MR_ivw_cur <- mr_mvivw(
         MRMVInputObject_cur,
         model = "random",
         robust = FALSE,
         correl = FALSE,
         distribution = "normal",
         alpha = 0.05
      )

      mr_results_MR_median_cur <- mr_mvmedian(
         MRMVInputObject_cur,
         distribution = "normal",
         alpha = 0.05,
         iterations = 10000,
         seed = 314159265
      )

      mr_results_MR_egger_cur <- mr_mvegger(
         MRMVInputObject_cur,
         orientate = 1,
         correl = FALSE,
         distribution = "normal",
         alpha = 0.05
      )

      # calculate sum r2 and mean F for each exposure
      r2.sum.insomnia <- sum(harmonised_data$exposure_r2[,1],na.rm=T)
      r2.sum.BMI <- sum(harmonised_data$exposure_r2[,2],na.rm=T)
      r2.sum.lifetimesmoking <- sum(harmonised_data$exposure_r2[,3],na.rm=T)

      ## calculate mean (median?) F statistic
      # insomnia
      ncases.insomnia <- subset(exp_cur_mv,exposure==ins.exp)$ncases.exposure[1]
      ncontrols.insomnia <- subset(exp_cur_mv,exposure==ins.exp)$ncontrols.exposure[1]
      nvars.insomnia <- sum(!is.na(harmonised_data$exposure_r2[,1]))
      meanF.insomnia <- (r2.sum.insomnia*((ncases.insomnia+ncontrols.insomnia)-1-nvars.insomnia))/((1-r2.sum.insomnia)*nvars.insomnia)
      # BMI
      nsamples.BMI <- subset(exp_cur_mv,exposure=="BMI_UKB_Neale_EUR_unrelated")$nsamples.exposure[1]
      nvars.BMI <- sum(!is.na(harmonised_data$exposure_r2[,2]))
      meanF.BMI <- (r2.sum.BMI*(nsamples.BMI-1-nvars.BMI))/((1-r2.sum.BMI)*nvars.BMI)
      # lifetime smoking
      nsamples.lifetimesmoking <- subset(exp_cur_mv,exposure=="LifetimeSmokingIndex_Wootton2019")$nsamples.exposure[1]
      nvars.lifetimesmoking <- sum(!is.na(harmonised_data$exposure_r2[,3]))
      meanF.lifetimesmoking <- (r2.sum.lifetimesmoking*(nsamples.lifetimesmoking-1-nvars.lifetimesmoking))/((1-r2.sum.lifetimesmoking)*nvars.lifetimesmoking)

      # store flattened results in the data frame
      mr_results_flat_mv <- rbind(mr_results_flat_mv,data.frame("exposure"=harmonised_data$expname$exposure[harmonised_data$expname$id.exposure=="Primary"],"outcome"=harmonised_data$outname$outcome,"Nsnp"=mr_results_MR_ivw_cur$SNPs,"Insomnia_IVW_beta"=as.numeric(mr_results_MR_ivw_cur$Estimate[1]),"Insomnia_IVW_SE"=as.numeric(mr_results_MR_ivw_cur$StdError[1]),"Insomnia_IVW_P"=as.numeric(mr_results_MR_ivw_cur$Pvalue[1]),"Insomnia_WM_beta"=mr_results_MR_median_cur$Estimate[1],"Insomnia_WM_SE"=mr_results_MR_median_cur$StdError[1],"Insomnia_WM_P"=mr_results_MR_median_cur$Pvalue[1],"Insomnia_MREgger_beta"=mr_results_MR_egger_cur$Estimate[1],"Insomnia_MREgger_SE"=mr_results_MR_egger_cur$StdError.Est[1],"Insomnia_MREgger_P"=mr_results_MR_egger_cur$Pvalue.Est[1],"Insomnia_R2Sum"=r2.sum.insomnia,"Insomnia_MeanF"=meanF.insomnia,"BMI_IVW_beta"=as.numeric(mr_results_MR_ivw_cur$Estimate[2]),"BMI_IVW_SE"=as.numeric(mr_results_MR_ivw_cur$StdError[2]),"BMI_IVW_P"=as.numeric(mr_results_MR_ivw_cur$Pvalue[2]),"BMI_WM_beta"=mr_results_MR_median_cur$Estimate[2],"BMI_WM_SE"=mr_results_MR_median_cur$StdError[2],"BMI_WM_P"=mr_results_MR_median_cur$Pvalue[2],"BMI_MREgger_beta"=mr_results_MR_egger_cur$Estimate[2],"BMI_MREgger_SE"=mr_results_MR_egger_cur$StdError.Est[2],"BMI_MREgger_P"=mr_results_MR_egger_cur$Pvalue.Est[2],"BMI_R2Sum"=r2.sum.BMI,"BMI_MeanF"=meanF.BMI,"LifetimeSmoking_IVW_beta"=as.numeric(mr_results_MR_ivw_cur$Estimate[3]),"LifetimeSmoking_IVW_SE"=as.numeric(mr_results_MR_ivw_cur$StdError[3]),"LifetimeSmoking_IVW_P"=as.numeric(mr_results_MR_ivw_cur$Pvalue[3]),"LifetimeSmoking_WM_beta"=mr_results_MR_median_cur$Estimate[3],"LifetimeSmoking_WM_SE"=mr_results_MR_median_cur$StdError[2],"LifetimeSmoking_WM_P"=mr_results_MR_median_cur$Pvalue[3],"LifetimeSmoking_MREgger_beta"=mr_results_MR_egger_cur$Estimate[3],"LifetimeSmoking_MREgger_SE"=mr_results_MR_egger_cur$StdError.Est[3],"LifetimeSmoking_MREgger_P"=mr_results_MR_egger_cur$Pvalue.Est[3],"LifetimeSmoking_R2Sum"=r2.sum.lifetimesmoking,"LifetimeSmoking_MeanF"=meanF.lifetimesmoking,"MREggerInt_beta"=mr_results_MR_egger_cur$Intercept,"MREggerInt_SE"=mr_results_MR_egger_cur$StdError.Int,"MREggerInt_P"=mr_results_MR_egger_cur$Pvalue.Int))

   }

   # write snplist to file for this ins.exp
   write.csv(snplist,file=paste0("MR_SNPlists/snplist_",ins.exp,"_PRIM_EXP.csv"),row.names=F,quote=F)

}

# write full set of flat results
write.csv(mr_results_flat_mv,file=paste0("MR_results/mr_results_mv_compact_all.csv"),row.names=F,quote=F)

