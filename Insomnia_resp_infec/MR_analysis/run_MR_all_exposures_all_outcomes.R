# load TwoSampleMR library
library(TwoSampleMR)
library(ggplot2)
library(viridis)

# source file containing winner's curse functions
source("mr_winners_curse_correction.R")


#########################
# load in exposure data #
#########################

# insomnia (Watanabe)
exp_insomnia_dat <- read_exposure_data(
   filename = "watanabe_insomnia_correctEAF_CPRA.txt",
   sep = "\t",
   snp_col = "RSID",
   beta_col = "BETA",
   se_col = "SE",
   effect_allele_col = "EA",
   other_allele_col = "OA",
   eaf_col = "EAF",
   pval_col = "P",
)
exp_insomnia_dat$exposure <- "Insomnia_Watanabe2021"
exp_insomnia_dat$units.exposure <- "logOR"
exp_insomnia_dat$type.exposure <- "binary"
exp_insomnia_dat$ncases.exposure <- 593724
exp_insomnia_dat$ncontrols.exposure <- 1771286

# insomnia (Lane)
exp_insomniafreq_dat <- read_exposure_data(
   filename = "lane_insomnia_CPRA.txt",
   sep = "\t",
   snp_col = "SNP",
   beta_col = "BETA",
   se_col = "SE",
   effect_allele_col = "EA",
   other_allele_col = "OA",
   eaf_col = "EAF",
   pval_col = "P",
)
exp_insomniafreq_dat$exposure <- "Insomnia_Lane2019"
exp_insomniafreq_dat$units.exposure <- "logOR"
exp_insomniafreq_dat$type.exposure <- "binary"
exp_insomniafreq_dat$ncases.exposure <- 129270
exp_insomniafreq_dat$ncontrols.exposure <- 108357

#short sleep
exp_shortsleep_dat <- read_exposure_data(
   filename = "dashti_short_sleep_CPRA.txt",
   sep = "\t",
   snp_col = "SNP",
   beta_col = "BETA",
   se_col = "SE",
   effect_allele_col = "EA",
   other_allele_col = "OA",
   eaf_col = "EAF",
   pval_col = "P",
)
exp_shortsleep_dat$exposure <- "ShortSleep_Dashti2019"
exp_shortsleep_dat$units.exposure <- "logOR"
exp_shortsleep_dat$type.exposure <- "binary"
exp_shortsleep_dat$ncases.exposure <- 106192
exp_shortsleep_dat$ncontrols.exposure <- 305742

# number sleep episodes
exp_nosleepepisodes_dat <- read_exposure_data(
   filename = "jones_number_sleep_episodes_CPRA.txt",
   sep = "\t",
   snp_col = "SNP",
   beta_col = "BETA",
   se_col = "SE",
   effect_allele_col = "EA",
   other_allele_col = "OA",
   eaf_col = "EAF",
   pval_col = "P",
)
exp_nosleepepisodes_dat$exposure <- "NoSleepEpisodes_Jones2019"
exp_nosleepepisodes_dat$units.exposure <- "n_episodes"
exp_nosleepepisodes_dat$type.exposure="continuous"
exp_nosleepepisodes_dat$nsamples.exposure <- 84810

# build full list of SNPs to extract from outcome data
expSNP <- unique(c(exp_insomnia_dat$SNP,exp_insomniafreq_dat$SNP,exp_sleepdur_dat$SNP,exp_shortsleep_dat$SNP,exp_nosleepepisodes_dat$SNP))



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

## FinnGen R9 phenotypes
r9phenos=c("CX_ACUTEUPPERINFEC","CX_INFLUENZA")

for(pheno in r9phenos) {
   eval(parse(text=paste0("out_",pheno,"_dat <- read_outcome_data(filename=\"FinnGen_R9_summstats/",pheno,"/finngen_R9_",pheno,"_b37_CPRA.regenie.gz\",snps=expSNP,snp_col=\"ID\",sep=\"\\t\",beta_col=\"BETA\",se_col=\"SE\",effect_allele_col=\"ALLELE1\",other_allele_col=\"ALLELE0\",pval_col=\"LOG10P\",eaf_col=\"A1FREQ\",log_pval=TRUE)")))
   eval(parse(text=paste0("out_",pheno,"_dat$outcome <- \"",pheno,"\"")))
   eval(parse(text=paste0("out_",pheno,"_dat$units.outcome <- \"logOR\"")))
   eval(parse(text=paste0("out_",pheno,"_dat$type.outcome <- \"binary\"")))
}

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
mr_results_flat <- data.frame("exposure"=character(),"outcome"=character(),"nsnp"=numeric(),"IVW_beta"=numeric(),"IVW_SE"=numeric(),"IVW_P"=numeric(),"WM_beta"=numeric(),"WM_SE"=numeric(),"WM_P"=numeric(),"MREgger_beta"=numeric(),"MREgger_SE"=numeric(),"MREgger_P"=numeric(),"MREggerInt_beta"=numeric(),"MREggerInt_SE"=numeric(),"MREggerInt_P"=numeric(),"R2Sum"=numeric(),"MeanF"=numeric())

# loop over all exposures
for(exposure in c("insomnia","insomniafreq","shortsleep","nosleepepisodes")) {

   # flip betas, effect alleles and eaf if beta<0
   eval(parse(text=paste0("exp_cur <- exp_",exposure,"_dat")))

   exp_cur$ean  <- exp_cur$effect_allele.exposure
   exp_cur$oan  <- exp_cur$other_allele.exposure
   exp_cur$eafn <- exp_cur$eaf.exposure
   exp_cur$bn   <- exp_cur$beta.exposure

   exp_cur$ean[exp_cur$beta.exposure<0]  <- exp_cur$other_allele.exposure[exp_cur$beta.exposure<0]
   exp_cur$oan[exp_cur$beta.exposure<0]  <- exp_cur$effect_allele.exposure[exp_cur$beta.exposure<0]
   exp_cur$eafn[exp_cur$beta.exposure<0] <- 1-exp_cur$eaf.exposure[exp_cur$beta.exposure<0]
   exp_cur$bn[exp_cur$beta.exposure<0]   <- -exp_cur$beta.exposure[exp_cur$beta.exposure<0]

   exp_cur$effect_allele.exposure <- exp_cur$ean
   exp_cur$other_allele.exposure  <- exp_cur$oan
   exp_cur$eaf.exposure           <- exp_cur$eafn
   exp_cur$beta.exposure          <- exp_cur$bn

   exp_cur[,c("ean","oan","eafn","bn")] <- list(NULL)

   eval(parse(text=paste0("exp_",exposure,"_dat <- exp_cur")))

   # if exposure is binary, calculate r^2 from logOR, otherwise use beta, SE and sample size
   if(exp_cur$type.exposure[1] == "binary") {
      exp_cur$r.exposure <- as.numeric(get_r_from_lor(exp_cur$beta.exposure,exp_cur$eaf.exposure,exp_cur$ncases.exposure,exp_cur$ncontrols.exposure,exp_cur$ncases.exposure/(exp_cur$ncases.exposure+exp_cur$ncontrols.exposure),model="logit",correction=F))
   } else if(exp_cur$type.exposure[1] == "continuous") {
      exp_cur$r.exposure <- as.numeric(get_r_from_bsen(b=exp_cur$beta.exposure,se=exp_cur$se.exposure,n=exp_cur$nsamples.exposure))
   }

   exp_cur$r2.exposure <- exp_cur$r.exposure^2

   # create SNP list using input exposure variants
   exp_snplist <- data.frame("SNP"=exp_cur$SNP,"R2"=exp_cur$r2.exposure)

   ### loop over all outcomes
   for (outcome in c("COVIDA2","COVIDB2","COVIDC2",r9phenos)) {

      # harmonise data (use eval as data frame name changes)
      harmonised_data <- eval(parse(text=paste0("harmonise_data(exp_",exposure,"_dat,out_",outcome,"_dat)")))

      # perform winner's curse correction
      wcc_stats <- wcc2(BetaXG=harmonised_data$beta.exposure,seBetaXG=harmonised_data$se.exposure,BetaYG=harmonised_data$beta.outcome,seBetaYG=harmonised_data$se.outcome)

      # NEED TO EITHER MAKE WCC2 HANDLE NEGATIVE EXP BETAS OR FLIP ALL BEFORE RUNNING WCC2

      # merge wcc stats back in to harmonised data
      harmonised_data$beta.exposure <- wcc_stats$BetaXGc
      harmonised_data$se.exposure <- wcc_stats$seBetaXGc
      harmonised_data$beta.outcome <- wcc_stats$BetaYGc
      harmonised_data$se.outcome <- wcc_stats$seBetaYGc

      # flag variants with missing exposure stats for removal
      # (as these variants had flipped or massively attenuated corrected effect size)
      harmonised_data$mr_keep[is.na(harmonised_data$beta.exposure)] <- FALSE

      # run analyses
      mr_results_all <- mr(harmonised_data,parameters=list(nboot=1000),method_list=c("mr_ivw_mre","mr_weighted_median","mr_egger_regression_bootstrap"))

      # extract only variants used in analysis for plotting
      harmonised_data_keep <- subset(harmonised_data,mr_keep)

      # get egger intercept stats and add to main results data frame
      mr_results_egger <- mr_egger_regression_bootstrap(harmonised_data_keep$beta.exposure,harmonised_data_keep$beta.outcome,harmonised_data_keep$se.exposure,harmonised_data_keep$se.outcome,parameters=list(nboot=1000))
      mr_results_egger_int <- mr_results_all[mr_results_all$method=="MR Egger (bootstrap)",]
      mr_results_egger_int$method <- "MR Egger (bootstrap) intercept"
      mr_results_egger_int$b      <- mr_results_egger$b_i
      mr_results_egger_int$se     <- mr_results_egger$se_i
      mr_results_egger_int$pval   <- mr_results_egger$pval_i
      rownames(mr_results_egger_int) <- max(as.numeric(rownames(mr_results_all))) + 1
      
      mr_results_all <- rbind(mr_results_all,mr_results_egger_int)

      # calculate total r2 for this exposure-outcome pair
      r2_sum <- sum(exp_cur$r2.exposure[exp_cur$SNP %in% harmonised_data_keep$SNP])

      # calculate mean (median?) F statistic
      if(exp_cur$type.exposure[1] == "binary") {
         meanF <- (r2_sum*((exp_cur$ncases.exposure[1]+exp_cur$ncontrols.exposure[1])-1-length(harmonised_data_keep$SNP)))/((1-r2_sum)*length(harmonised_data_keep$SNP))
      } else if(exp_cur$type.exposure[1] == "continuous") {
         meanF <- (r2_sum*((exp_cur$nsamples.exposure[1])-1-length(harmonised_data_keep$SNP)))/((1-r2_sum)*length(harmonised_data_keep$SNP))
      } else {
         meanF=NA
      }

      # add total r2 and F stat to the results data frame
      mr_results_all$R2Sum <- r2_sum
      mr_results_all$meanF <- meanF

      # get units and phenotype names for plot
      eval(parse(text=paste0("exp_units <- unique(exp_",exposure,"_dat$units.exposure)")))
      eval(parse(text=paste0("out_units <- unique(out_",outcome,"_dat$units.outcome)")))

      eval(parse(text=paste0("exp_name <- unique(exp_",exposure,"_dat$exposure)")))
      eval(parse(text=paste0("out_name <- unique(out_",outcome,"_dat$outcome)")))

      # create scatter plot
      #mr_scatter_plot(mr_results,harmonised_data)
      g <- ggplot(data=harmonised_data_keep,aes(x=beta.exposure,y=beta.outcome)) + geom_hline(yintercept=0,colour="black") + geom_point() + geom_errorbar(aes(ymin=beta.outcome-1.96*se.outcome,ymax=beta.outcome+1.96*se.outcome),width=0.02*diff(range(harmonised_data_keep$beta.exposure))) + geom_errorbarh(aes(xmin=beta.exposure-1.96*se.exposure,xmax=beta.exposure+1.96*se.exposure),height=0.02*diff(range(harmonised_data_keep$beta.outcome))) + geom_abline(aes(intercept=0,slope=mr_results_all$b[mr_results_all$method=="Inverse variance weighted (multiplicative random effects)"],colour="ivwcol")) + geom_abline(aes(intercept=mr_results_egger$b_i,slope=mr_results_all$b[mr_results_all$method=="MR Egger (bootstrap)"],colour="eggercol")) + geom_abline(aes(intercept=0,slope=mr_results_all$b[mr_results_all$method=="Weighted median"],colour="wmcol")) + ggtitle(paste0("Exposure: ",exp_name," || Outcome: ",out_name)) + xlab(paste0(exposure," beta (",exp_units,")")) + ylab(paste0(outcome," beta (",out_units,")")) + theme_bw() + scale_colour_manual(name="Method",values=c("ivwcol"=viridis(3)[1],"eggercol"=viridis(3)[2],"wmcol"=viridis(3)[3]),labels=c("IVW","MR Egger","WM")) + theme(legend.background=element_rect(fill="white",colour="black"),legend.position=c(1,1),legend.justification=c(1,1))
      ggsave(paste0("MR_plots/mr_plot_",exposure,"_EXP_",outcome,"_OUT.png"),plot=g,device="png")

      # write results to full file
      write.csv(mr_results_all,paste0("MR_results/mr_results_full_",exposure,"_EXP_",outcome,"_OUT.csv"),row.names=F,quote=F)

      # add column to SNPlist table for this outcome
      eval(parse(text=paste0("exp_snplist$",out_name," <- 0")))
      eval(parse(text=paste0("exp_snplist$",out_name,"[exp_snplist$SNP %in% harmonised_data_keep$SNP] <- 1")))

      # write flattened file (single line per exp-out pair)
      mr_results_flat_cur <- data.frame("exposure"=exp_name,"outcome"=out_name,"nsnp"=mr_results_all$nsnp[1],"IVW_beta"=mr_results_all$b[mr_results_all$method=="Inverse variance weighted (multiplicative random effects)"],"IVW_SE"=mr_results_all$se[mr_results_all$method=="Inverse variance weighted (multiplicative random effects)"],"IVW_P"=mr_results_all$pval[mr_results_all$method=="Inverse variance weighted (multiplicative random effects)"],"WM_beta"=mr_results_all$b[mr_results_all$method=="Weighted median"],"WM_SE"=mr_results_all$se[mr_results_all$method=="Weighted median"],"WM_P"=mr_results_all$pval[mr_results_all$method=="Weighted median"],"MREgger_beta"=mr_results_all$b[mr_results_all$method=="MR Egger (bootstrap)"],"MREgger_SE"=mr_results_all$se[mr_results_all$method=="MR Egger (bootstrap)"],"MREgger_P"=mr_results_all$pval[mr_results_all$method=="MR Egger (bootstrap)"],"MREggerInt_beta"=mr_results_all$b[mr_results_all$method=="MR Egger (bootstrap) intercept"],"MREggerInt_SE"=mr_results_all$se[mr_results_all$method=="MR Egger (bootstrap) intercept"],"MREggerInt_P"=mr_results_all$pval[mr_results_all$method=="MR Egger (bootstrap) intercept"],"R2Sum"=r2_sum,"MeanF"=meanF)

      write.csv(mr_results_flat_cur,paste0("MR_results/mr_results_compact_",exp_name,"_EXP_",out_name,"_OUT.csv"),row.names=F,quote=F)

      # concatenate results to existing flat data frame
      mr_results_flat <- rbind(mr_results_flat,mr_results_flat_cur)

   }

   # write the SNP list file for this exposure
   write.csv(exp_snplist,paste0("MR_SNPlists/snplist_",exp_name,"_EXP.csv"),row.names=F,quote=F)

}

# write full set of flat results
write.csv(mr_results_flat,paste0("MR_results/mr_results_compact_all.csv"),row.names=F,quote=F)

