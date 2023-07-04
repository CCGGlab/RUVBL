##### Multivariate Coxph analysis on all genes (13709) 
##### Covariates used (MYCN-amplification(Amp vs nonAmp),Age at diagnosis(>18months vs <18months) and INSS stage (Stage 4 vs stages 1,2,3))
##### data source: Davide Cangelosi et al., 2020 (PMID: 32825087): Data from integrated NB tumour platforms.


library(rlang)
library(tidyverse)
library(survival)
library(survminer)
library(survivalAnalysis) #use the  cox_as_data_frame function to export coxph results



# load data
clin <- read_csv("data/clinical_data_Cangelosi_formated.csv")
exprs <- read_rds("data/exprs_data_Cangelosi_formated.rds")

########## filtering dataset
## Cangelosi dataset is an integrated dataset of multiple publicly available datasets which include SEQC-RNAseq
## SEQC-RNAseq dataset is being used in this study as an independent dataset hence needs to filtered out of Canglosis
# load clincal data of SEQC datatset
clin_seqc <- read_rds("data/clin_SEQC.rds")

# gene IDs of samples common in both Cangelosi and SEQC
common_with_SEQC <- intersect(toupper(clin$`Patient ID`), clin_seqc$GSE)

## filter out SEQC samples from Cangelosi clinical and exprs data
clin <- clin %>%  dplyr::filter(!`Patient ID` %in% tolower(common_with_SEQC)) 
# order the rownames (sample names)
clin <- clin %>%  dplyr::arrange(as.character(clin$`Patient ID`)) 
exprs <- exprs %>%  dplyr::select(!tolower(common_with_SEQC))   
# order columns (sample names)
exprs <- exprs %>% select(order(as.character(colnames(exprs)))) 

# reformat data for coxph analysis
clin <- clin %>%
      # set Age_group >=18 to 1 and <18 to 0
      dplyr::mutate(Age=case_when(Age_group==">=18"~ 1,TRUE~0)) %>%
      # set MYCN Amplified to 1 and not amplified to 0
      dplyr::mutate(MYCN_status=case_when(MYCN_status=="Amplified" ~ 1,TRUE~ 0)) %>%
      dplyr::mutate(INSS_group=case_when(INSS_stage==4 ~ 1,TRUE~0))

#### find the median expression value gene and categorise gene expression as high or low in each sample
# genes of interest
all_genes <- rownames(exprs)
# median expr and isHigh
future::plan(future::multisession(workers = 20))
isHigh  <-  furrr::future_map(all_genes, function(gene){
      exprs_median<- median(as.numeric(exprs[gene,]), na.rm=T) # median expression value
      isHigh<- as.numeric(exprs[gene,])>exprs_median
      tibble({{gene}} := isHigh)
})

# append to clin data
clin  <- bind_cols(clin , isHigh )

### Cox Proportional Hazards Regression
# surv object
surv_object <-Surv(time = as.numeric(clin$overall_survival), event = clin$Event_overall=="yes")

# for paralleling
future::plan(future::multisession(workers = 20))
# run coxph
res_coxph_df_cangelosi  <-  furrr::future_map(all_genes , safely(function(gene){
      genesym <- as.symbol(gene)
      res <-  rlang::inject(coxph(surv_object ~ !!genesym+MYCN_status+Age+INSS_group, data = clin ))
      # convert coxph result to dataframegene
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
readr::write_rds(res_coxph_df_cangelosi,"data/res_coxph_all_genes_cangelosi.rds")

# # xlsx
# openxlsx::write.xlsx(res_coxph_df_cangelosi,"/results/tables/res_coxph_all_genes_cangelosi.xlsx" )
