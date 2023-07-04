
##### Multivariate Coxph analysis on all genes  
## Covariates used (MYCN-amplification(Amp vs nonAmp),Age at diagnosis(>18months vs <18months) and INSS stage (Stage 4 vs stages 1,2,3))
## data : SEQC neuroblastoma patient cohort (GEO; accession nr GSE49710)

library(rlang)
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
library(survivalAnalysis) #use the  cox_as_data_frame function to export coxph results
library(openxlsx)


## load data
#SEQC
exprs  <- read_rds("data/exprs_SEQC_RNA_log2RPM.rds")
clin  <- read_rds("data/clin_SEQC.rds")


## format clin and or exprs data
# SEQC
clin   <- clin  %>% filter(MYCN_status!="unknown") %>%
      dplyr::mutate(MYCN_status=case_when(MYCN_status=="Amplified" ~ 1, #set MYCN Amplified to 1 and not amplified to 0
                                          TRUE~ 0)) %>%
      dplyr::mutate(INSS_group=case_when(INSS_stage == 4 ~ 1, # set INSS stage 4 to 1 and all others to 0
                                         TRUE~0))  %>%
      dplyr::mutate(age_months=as.numeric(`Age at Diagnosis in Days`)/31) %>%
      dplyr::mutate(Age_group=case_when(age_months>=18~1, TRUE~0)) %>%  # set age >=18 to 1 and < 18months to 0
      dplyr::mutate(Age=Age_group)

exprs  <- exprs  %>%  dplyr::select(all_of(clin$title))



## classify gene expression as high or low (isHigh)  based on MEDIAN expression of each gene.
## do this for all genes in the dataset

# genes of interest
all_genes<- rownames(exprs)

## for each gene determine isHigh and append result to clin data
#SEQC
future::plan(future::multisession(workers = 100))
isHigh  <-  furrr::future_map(all_genes, function(gene){
      exprs_median<- median(as.numeric(exprs[gene,]), na.rm=T) # median expression value
      isHigh<- as.numeric(exprs[gene,])>exprs_median
      tibble({{gene}} := isHigh)
})
# append to clin data
clin  <- bind_cols(clin , isHigh )

#### Coxph
# surv object
surv_object <-Surv(time = as.numeric(clin$overall_survival), event = clin$`Vital Status`=="Dead")

## multivariate Coxph
# for paralleling
future::plan(future::multisession(workers = 100))
# run coxph
res_coxph_df_seqc  <-  furrr::future_map(all_genes , safely(function(gene){
      genesym <- as.symbol(gene)
      res <-  rlang::inject(coxph(surv_object ~ !!genesym+MYCN_status+Age+INSS_group, data = clin ))
      # convert coxph result to dataframe
      res <- cox_as_data_frame(res)
      # format  dataframe
      res <- res  %>% dplyr::mutate(across(c("HR","Lower_CI","Upper_CI"), round,2)) %>%
            dplyr::mutate("95% CI"=paste(Lower_CI,Upper_CI, sep = " - ")) %>%
            dplyr::mutate("95% CI"=paste0("(",`95% CI`,")")) %>%
            dplyr::rename("P value" = p) %>% 
            # select columns of interest and filter out only row with gene info
            dplyr::select(factor.name, HR, "95% CI","P value" ) %>%
            dplyr::filter(factor.name==gene)
      # save result as tibble
      tibble({{gene}} := res)})) %>%
      # filter for genes that have been sucessfully analyzed
      map("result") %>% # extract result layer from the tibble
      # remove all empty objects from tibble
      compact() %>%
      #bind  coxph results of all genes by rows into single dataframe
      map_dfr(flatten)


## save data
# rds
readr::write_rds(res_coxph_df_seqc,"data/res_coxph_all_genes_seqc.rds")
# # xlsx
# openxlsx::write.xlsx(res_coxph_df_seqc,"results/tables/res_coxph_all_genes_seqc.xlsx" )


