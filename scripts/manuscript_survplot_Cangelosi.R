### Analysis of neuroblastoma tumor survival data from Cangelosis intergrated platform dataset
#data source: Davide Cangelosi et al., 2020 (PMID: 32825087): Data from integrated NB tumour platforms.

library(tidyverse)
library(survival)
library(survminer)
library(patchwork)

# load data
clin <- read_csv("data/clinical_data_Cangelosi_formated.csv")
exprs <- read_rds("data/exprs_data_Cangelosi_formated.rds")

### filtering dataset
## Cangelosi dataset is an integrated dataset of multiple publicly available datasets which include NB-SEQC
## NB-SEQC dataset is being used in this study as an independent dataset hence needs to filtered out of Cangelosi's

# load clincal data for SEQC datatset
clin_seqc <- read_rds("data/clin_SEQC.rds")

# gene IDs of samples common in both Cangelosi and SEQC
common_with_SEQC <- intersect(toupper(clin$`Patient ID`), clin_seqc$GSE)

# filter out SEQC samples from Cangelosi clinical and exprs data
clin <- clin %>%  dplyr::filter(!`Patient ID` %in% tolower(common_with_SEQC))
exprs <- exprs %>%  dplyr::select(!tolower(common_with_SEQC))

#### surv analysis

## find the median expression value for RUVBL1 and RUVBL2 and categorise RUVBL expression as high or low in each sample 

# genes of interest
genes <- c("RUVBL1", "RUVBL2")

for(gene in genes){
      exprs_median<- median(as.numeric(exprs[gene,]), na.rm=T) # median expression value
      isHigh<- as.numeric(exprs[gene,])>exprs_median
      names(isHigh)<- colnames(exprs)
      rownames(clin)<- clin$`Patient ID`
      clin$isHigh<- isHigh[rownames(clin)]
      clin<- clin[!is.na(clin$isHigh),] %>% dplyr::mutate(isHigh=case_when(isHigh==TRUE ~ 1, TRUE~ 0)) 
      clin <- clin %>%  dplyr::rename("{gene}":=isHigh)
}


## KM overall and event-free survival plots (method-median)
OS_list <-  purrr::map(genes,function(gene){
      gene <- as.symbol(gene)
      # fit
      fit <- rlang::inject(survfit(Surv(as.numeric(clin$overall_survival),event = clin$Event_overall=="yes")~ !!gene, data=clin))
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

# plot
OS_list[[1]] # RUVBL1
OS_list[[2]] # RUVBL2


## EF survival
EF_list <-  purrr::map(genes,function(gene){
      gene <- as.symbol(gene)
      # fit
      fit <- rlang::inject(survfit(Surv(time = as.numeric(clin$EF_survival), 
                                        event = clin$Event_free=="yes")~ !!gene, data=clin))
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
EF_list[[1]]
EF_list[[2]]
