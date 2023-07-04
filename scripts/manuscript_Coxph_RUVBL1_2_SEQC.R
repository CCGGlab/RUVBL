#### Cox Proportional Hazards Regression Analysis of NB-SEQC dataset

# load packages
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)


# load data
exprs <- read_rds("data/exprs_SEQC_RNA_log2RPM.rds")
clin <- read_rds("data/clin_SEQC.rds")

### find the median expression value for RUVBL1 and RUVBL2 and categorise RUVBL expression as high or low in each sample 
genes <- c("RUVBL1", "RUVBL2")

## make clin data for each gene
clin_list <- list()
for(gene in genes){
      # median expression value and isHigh
      exprs_median<- median(as.numeric(exprs[gene,]), na.rm=T) 
      isHigh<- as.numeric(exprs[gene,])>exprs_median
      names(isHigh)<- colnames(exprs)
      
      rownames(clin)<- clin$title
      clin$isHigh<- isHigh[rownames(clin)]
      clin<- clin[!is.na(clin$isHigh),]
      
      ## format/rename variables
      clin  <- clin %>% dplyr::filter(MYCN_status!="unknown")  %>% # remove unwanted MYCN variables
            dplyr::mutate(MYCN=case_when(MYCN_status=="Amplified" ~ 1, TRUE~ 0)) %>% #set MYCN Amplified to 1 and not amplified to 0
            dplyr::mutate(isHigh=case_when(isHigh=="TRUE" ~ 1, TRUE~ 0)) %>% # set isHigh=="TRUE" to 1 and FALSE to 0
            dplyr::mutate(Stage=case_when(INSS_stage == 4 ~ 1,  TRUE~0))  %>% # set INSS stage 4 to 1 and all others to 0
            dplyr::mutate(age_months=as.numeric(`Age at Diagnosis in Days`)/31) %>% # convert age from days to months
            dplyr::mutate(Age=case_when(age_months>=18~1, TRUE~0)) # set age >=18 to 1 and <18months to 0
      # save clin to list
      clin_list[[gene]] <- clin
}
names(clin_list) <- genes


# create survival object

surv_obj_list <-clin_list %>%   map(function(clin)
      Surv(time = as.numeric(clin$overall_survival), event = clin$`Vital Status`=="Dead"))

### run Coxph  and plot result as forest plot
# covariates to consider
# isHigh refers to where RUVBL1 or RUVBL2 is highly expressed
vars <- c("isHigh","MYCN","Age","Stage")

# coxph
source("/home/joachim/projects/RUVBL2/scripts/functions/plot_CPH_forest_adapt.R")# function to make forest plot of Coxph result
cph_fplot_list <- purrr::map2(surv_obj_list,clin_list, function(surv_obj,clin)
      plot_CPH_forest(clin_merged = clin,surv_object = surv_obj,vars = vars))

## rename plot variables
# rename isHigh to gene name (RUVBL1 or RUVBL2) in the plot data and relevel the covariates
cph_fplot_list <- purrr::map2(cph_fplot_list, genes, function(plot, gene){
      # rename
      plot$data <- plot$data %>% dplyr::mutate(cond=case_when(cond=="isHigh"~ gene,TRUE ~ as.character(cond))) 
      # relevel covariates
      plot$data$cond  <-  factor(as.factor(plot$data$cond), levels = c("Stage","Age","MYCN",paste0(gene)))
      plot
})

# patch plots
(cph_fplot_list[[1]]|cph_fplot_list[[2]])+plot_annotation("Multivariate Coxph (NB-SEQC)")
         