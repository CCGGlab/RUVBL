### This script plots overal and Event-free survival curves for RUVBL1 and RUVBL2 expression 
# in SEQC neuroblastoma patient cohort (GEO; accession nr GSE49710)

# load libraries
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)


## load data
exprs <- read_rds("data/exprs_SEQC_RNA_log2RPM.rds")
clin <- read_rds("data/clin_SEQC.rds")

### find the median expression value for RUVBL1 and RUVBL2 
# and categorise RUVBL expression as high or low in each sample 
genes <- c("RUVBL1", "RUVBL2")


for(gene in genes){
      exprs_median<- median(as.numeric(exprs[gene,]), na.rm=T) # median expression value
      isHigh<- as.numeric(exprs[gene,])>exprs_median
      names(isHigh)<- colnames(exprs)
      
      rownames(clin)<- clin$title
      clin$isHigh<- isHigh[rownames(clin)]
      clin<- clin[!is.na(clin$isHigh),] %>% 
            dplyr::mutate(isHigh=case_when(isHigh==TRUE ~ 1, TRUE~ 0)) 
      clin <- clin %>%  dplyr::rename("{gene}":=isHigh)
}

## Kaplan Meier overall and event-free survival plots (method-median)
# genes to plot
genes <- c("RUVBL1","RUVBL2")

# OS 
OS_list <-  purrr::map(genes,function(gene){
      gene <- as.symbol(gene)
      # fit
      fit <- rlang::inject(survfit(Surv(as.numeric(clin$overall_survival),
                          event = clin$`Vital Status`=="Dead")~!!gene, data=clin))
      # surv plot
      surv_plot <- ggsurvplot(fit, data = clin, pval = TRUE, pval.method=T, 
                              legend.title = as.character(gene), 
                              legend.labs = c(paste0("Low (n=",sum(clin[[gene]]),")"),
                                              paste0("High (n=",sum(!clin[[gene]]),")")), 
                              risk.table = TRUE, xlab = "Time (years)",
                              ylab="Overall survival probability",
                              legend = c(0.7,0.3),
                              palette  = c("#332FD0","#FF1818" ))
      
      surv_plot
})

# name list
OS_list <- setNames(OS_list,genes)
#plot
OS_list[[1]] # RUVBL1
OS_list[[2]] # RUVBL2

## EF survival
EF_list <-  purrr::map(genes,function(gene){
      gene <- as.symbol(gene)
      # fit
      fit <- rlang::inject(survfit(Surv(time = as.numeric(clin$EF_survival), 
                          event = clin$`First Event`=="Event")~!!gene, data=clin))
      # surv plot
      surv_plot <- ggsurvplot(fit, data = clin, pval = TRUE, pval.method=T, 
                              legend.title = as.character(gene), 
                              legend.labs = c(paste0("Low (n=",sum(clin[[gene]]),")"),
                                              paste0("High (n=",sum(!clin[[gene]]),")")), 
                              risk.table = TRUE, xlab = "Time (years)",
                              ylab="Event-free survival probability",
                              legend = c(0.7,0.2),
                              palette  = c("#332FD0","#FF1818" ))
      surv_plot
})
      

## Multivariate Kaplan-Meier analysis

# remove unwanted variables and regroup covariates
clin <- clin %>% filter(MYCN_status!="unknown") %>% 
      dplyr::mutate(MYCN=case_when(MYCN_status=="Amplified" ~ 1,  TRUE~ 0)) %>% 
      dplyr::mutate(INSS_group=case_when(INSS_stage == 4 ~ 1, TRUE~0)) %>%# set INSS stage 4 to 1 and all others to 0
      dplyr::mutate(age_months=as.numeric(`Age at Diagnosis in Days`)/31) %>% 
      dplyr::mutate(Age=case_when(age_months>=18~1, TRUE~0)) # set age >=18 to 1 and < 18months to 0

# surv plot
# RUVBL1+MYCN

fit <- survfit(Surv(as.numeric(clin$overall_survival), event = clin$`Vital Status`=="Dead")~ RUVBL1+MYCN, data=clin)

surv_plot2 <- ggsurvplot(fit, data = clin, pval = TRUE, pval.method=T,risk.table = T,
                        xlab="Time (years)", ylab="Overall survival probability",legend=c(0.8,0.2),
                        palette = c("#D4ADFC","#5C469C","#1D267D","#1B9C85"))
surv_plot2
