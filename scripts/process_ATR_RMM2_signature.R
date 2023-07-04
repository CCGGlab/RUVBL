### Create ATR and RMM2 signature genesets using data ATR and RMM2 data from Szydzik et al,2021 and Nunes et al,2022 respectively.

# Data from NB cell lines either treated with ATRi(Elimusertib (BAY1895344)) or RMM2i(3AP) for 48hrs
# ATRi data from CLB-BA cells
# RRM2 data from CLB-GA and IMR-32 cells

library(enrichplot)
library(ggnewscale)
library(clusterProfiler)
library(msigdbr)
library(Hmisc)
library(scales)
library(openxlsx)
library(tidyverse)


# load prcoess DESeq2 data for ATRi and RMM2i treated neuroblastoma cell lines 
ATR_DE <- readr::read_rds("data/ATR_DE.rds")
RRM2_DE <- readr::read_rds("data/RRM2_Nunes_2022.rds")

# ATR upregulated regulated DEGs  
ATR_up <-as.data.frame(ATR_DE$BAR$`48`) %>%  
      dplyr::filter(padj<= 0.01 ) %>% 
      dplyr::filter(log2FoldChange>=1) %>% 
      dplyr::pull(HGNC) 

# ATR downregulated regulated DEGs      
ATR_down <-as.data.frame(ATR_DE$BAR$`48`) %>%  
      dplyr::filter(padj<= 0.01 ) %>% 
      dplyr::filter(log2FoldChange<=-1) %>% 
      dplyr::pull(HGNC)   

# all ATR DEGs
ATR_all <- as.data.frame(ATR_DE$BAR$`48`) %>%  
      dplyr::filter(padj<= 0.01 ) %>% 
      dplyr::filter(abs(log2FoldChange)>=1) %>% 
      dplyr::pull(HGNC)

## Nune RRM2 data
Nunes_rmm2 <- list(clbGa=RRM2_DE$RNA$CLBGA$`3AP_IC50`,imr=RRM2_DE$RNA$IMR32$`3AP_IC50`)

## reduce list levels in data
Nunes_rmm2 <- Nunes_rmm2 %>% map(~flatten(.) %>%  as.data.frame(.)) 

# get DEGs for RRM2i at 48 hrs
RRM2_sig <- map2(Nunes_rmm2,names(Nunes_rmm2), function(res,name){
      # RRM2 upregulated DEGs
      RRM2_up <-res %>%  
            dplyr::filter(padj<= 0.01 ) %>% 
            dplyr::filter(log2FoldChange>=1) %>% 
            dplyr::pull(HGNC)
      #RRM2 downregulated DEGs
      RRM2_down <-res %>%  
            dplyr::filter(padj<= 0.01 ) %>% 
            dplyr::filter(log2FoldChange<=-1) %>% 
            dplyr::pull(HGNC)
      # all RRM2 DEGs
      RRM2_all <- res %>%  
            dplyr::filter(padj<= 0.01 ) %>% 
            dplyr::filter(abs(log2FoldChange)>=1) %>% 
            dplyr::pull(HGNC)
      # list of DEGs
      res2 <- list(  RRM2_down,  RRM2_up, RRM2_all) %>% 
            setNames(., c(paste0("RRM2_down_",name),paste0("RRM2_up_",name),
                         paste("RRM2_all_",name)))
} ) 
      

# find common genes for between clbGa and imr32 for each geneset category

RRM2_sig_ClbGa <- RRM2_sig$clbGa 
RRM2_sig_imr <- RRM2_sig$imr

RRM2_sig <- map2(RRM2_sig_ClbGa,RRM2_sig_imr, function(ga, imr)
      intersect(ga,imr)) %>% 
      setNames(., c("RRM2_down","RRM2_up","RRM2_all"))


# make dataframe of all gene sets
max_length <- max(length(ATR_all),length(ATR_down),length(ATR_up), 
                  length(RRM2_sig$RRM2_down),length(RRM2_sig$RRM2_all),
                  length(RRM2_sig$RRM2_up))

ATR_RRM2_sig <- data.frame(ATR_down=c(ATR_down, rep(NA, max_length-length(ATR_down))),
                      ATR_up=c(ATR_up, rep(NA, max_length-length(ATR_up))),
                      ATR_all=c(ATR_all, rep(NA, max_length-length(ATR_all))),
                      RRM2_down=c(RRM2_sig$RRM2_down,rep(NA, max_length-length(RRM2_sig$RRM2_down))),
                      RRM2_up=c(RRM2_sig$RRM2_up,rep(NA, max_length-length(RRM2_sig$RRM2_up))),
                      RRM2_all=c(RRM2_sig$RRM2_all,rep(NA, max_length-length(RRM2_sig$RRM2_all))))
# covert to long data format
ATR_RRM2_sig <- ATR_RRM2_sig %>%  pivot_longer(cols  = "ATR_down":"RRM2_all",
                              names_to = "gs_name",values_to = "gene_symbol",
                              values_drop_na = TRUE)

# save
# write_rds(ATR_RRM2_sig, "data/ATR_RRM2_sig_geneset_clusterPformat.rds")
