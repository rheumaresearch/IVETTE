#################################################################################################################
#################################################################################################################
### Project: jcalvet_202110_marato
### Description: Proteomics study in plasma of COVID-19 patients in collaboration with IRSICaixa ad Sant Pau -
###              La Marato TV3 - Design computations in jcalvet_202003_proteomics
### Author: Toni B.
### Date: 2021/10/06
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.wd <- "/home/arodrigo/Documents/ivette_marato_publicRepo/"
.dat <- paste0(.wd, "data/")
.res0 <- paste0(.wd, "reports/res/")

setwd(.wd)
dir.create(.res0, F)

library(affy)
library(cluster)
library(hwriter)
library(multcomp)
library(amap)
library(Minirand)
library(gtools)
library(parallel)
#library(devtools)
#install_github("marchtaylor/sinkr")
#library(githubinstall)
#githubinstall("sinkr")
library(sinkr)
library(quantreg)
library(limma)
library(lme4)
library(lmerTest)
library(psych)
# library(sva)
library(iq)
library(ggplot2)
library(ggrepel)
library(eulerr)
library(VennDiagram)
library(colorspace)
library(ppcor)
library(randomForest)
library(pROC)
library(glmnet)
library(minpack.lm)
library(drc)
library(roastgsa)
library(combinat)

library(ggpubr)
library(rstatix)
library(ggtext)
library(openxlsx)
library(Hmisc)
library(gt)
library(ggplot2)
library(emmeans)
library(gtsummary)
library(patchwork)



#################################################################################################################
# G. Read and analyze DIA data - FragPipe + MAXLFQ

### G2. Differential expression - MaxLFQ
### G4. Differential expression with ordinal models
### G5. Enrichment analysis with ROAST-GSA
### G5.2. Enrichment analysis with ROAST-GSA - lineal

source("code/G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression.R")
source("code/G.Discovery_DIA_FP_MaxLFQ/G4.Diff_expression_ordinal.R")
source("code/G.Discovery_DIA_FP_MaxLFQ/G5.ROAST_GSA.R")
source("code/G.Discovery_DIA_FP_MaxLFQ/G5.2.ROAST_GSA_linear.R")
source("code/G.Discovery_DIA_FP_MaxLFQ/G5.3.ROAS_GSA_summary.R")



#################################################################################################################
# M. Validation experiments

### M5.3. Differential expression for PRM data - 1st round - no DIA samples - adjusted by sex, age and vaccine
###      (no controls)
### M12.3. Differential expression for PRM data - 2nd round - no DIA samples - adjusted by sex, age and vaccine (no controls)
### M20.3. Differential expression for PRM data - 3rd round - no DIA samples - adjusted by sex, age and vaccine (no controls)
### M25.2. Differential expression for PRM data - 4th round - no DIA samples - adjusted by sex, age and vaccine (no controls)
### M26. Differential expression for ordinal model - all rounds - adjusted by sex, age and vaccine (no controls)

source("code/M.Validation_exp/M5.3.Diff_expression_PRM_1st.round_MaxLQF_no.DIA_adj.sex.age.vaccine.R")
source("code/M.Validation_exp/M12.3.Diff_expression_PRM_2nd.round_MaxLQF_no.DIA_adj.sex.age.vaccine.R")
source("code/M.Validation_exp/M20.3.Diff_expression_PRM_3rd.round_MaxLQF_no.DIA_adj.sex.age.vaccine.R")
source("code/M.Validation_exp/M25.2.Diff_expression_PRM_4th.round_MaxLQF_no.DIA_adj.sex.age.vaccine.R")
source("code/M.Validation_exp/M26.3.Diff_expression_PRM_all_ordinal_adjSexAgeVacc.R")
source("code/M.Validation_exp/R5.2.PRM_figs_tables_adjSex.age.vacc")


#################################################################################################################
  # O. Formal analyses of concordance between technologies - PRM all rounds

### O2. Differential expression Hospitalized vs Not hospitalized - PRM
### O3. Differential expression Hospitalized vs Not hospitalized - DIA, proteins in PRM
### O4. Concordance analysis - DIA - PRM

source("code/O.Concordance_DIA_PRM_complete/O2.Diff_exprs_hosp_PRM.R")
source("code/O.Concordance_DIA_PRM_complete/O3.Diff_exprs_hosp_DIA_sel.R")
source("code/O.Concordance_DIA_PRM_complete/O4.Concordance_DIA_PRM.R")

#################################################################################################################
  # P. Multivariate predictors - DIA -> PRM

### P1. Multivariate predictor
### P2. Predictor for Not-hospitalized vs Hospitalized
### P3. Predictor for Severe vs Mild
### P5. Prediction ability of validated proteins - univariate

source("code/P.Multiv_predictor/P1.Data_setup.R")
source("code/P.Multiv_predictor/P2.Multiv_pred_hosp_2.R")
source("code/P.Multiv_predictor/P3.2.Multiv_pred_sev.mild_2.R")
source("code/P.Multiv_predictor/P5.Predict_univ.R")

#################################################################################################################
  # S. ELISA analyses

### S2.7.Diff_expression_validated_adjSexAgeVacc_draft
#     Comparisons between severity groups for proteins analyzed by ELISA
### S2.8.Diff_expression_validated_hospVsNotHosp_adjSexAgeVacc_draft
#     Comparisons between hospitalized and non.hospitalized  for proteins analyzed by ELISA

source("code/S.ELISA_validation/S2.7.Diff_expression_validated_adjSexAgeVacc_draft.R")
source("code/S.ELISA_validation/S2.8.Diff_expression_validated_hospVsNotHosp_adjSexAgeVacc_draft.R")

#################################################################################################################
  # A. Tables and Figures for drafts

### A1. Descriptives tables for Discovery series
### A2. Descriptives tables for the PRM validation series
### A3. Descriptives tables for the ELISA validation series
### A4. Descriptives tables for all patients in the biomarkers study
### R5. Figures and tables for PRM study
### R6. Overall view of the DIA data

source("code/A.Read_data/A1.Descriptives_disc.R")
source("code/A.Read_data/A2.Descriptives_prm.R")
source("code/A.Read_data/A3.Descriptives_elisa.R")
source("code/A.Read_data/A4.Descriptives_all.R")
source("code/M.Validation_exp/R5.2.PRM_figs_tables_adjSex.age.vacc.R")
source("code/G.Discovery_DIA_FP_MaxLFQ/R6.DIA_overall_view.R")



#################################################################################################################
#################################################################################################################
#################################################################################################################







