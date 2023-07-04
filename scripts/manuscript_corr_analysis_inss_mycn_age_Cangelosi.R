
### Correlation of RUVBL1 and RUVBL2  with INSS stages, MYCN amplification and Age
#data source: Davide Cangelosi et al., 2020 (PMID: 32825087): Data from integrated NB tumour platforms.

library(tidyverse)
library(patchwork)

# load data

clin <- read_csv("data/clinical_data_Cangelosi_formated.csv")

exprs <- read_rds("data/exprs_data_Cangelosi_formated.rds")

########## filtering dataset
## Cangelosi dataset is an integrated dataset of multiple publicly available datasets which include SEQC-RNAseq
## SEQC-RNAseq dataset is being used in this study as an independent dataset hence needs to filtered out of Cangelosi's

# load clincal data for SEQC datatset
clin_seqc <- read_rds("data/clin_SEQC.rds")

# gene IDs of samples common in both Cangelosi and SEQC
common_with_SEQC <- intersect(toupper(clin$`Patient ID`), clin_seqc$GSE)
# #save common samples
# write_rds(common_with_SEQC, "~/data/sample_ids_common_to_SEQC_Cangelosi.rds")

# filter out SEQC samples from Cangelosi clinical and exprs data
clin <- clin %>%  dplyr::filter(!`Patient ID` %in% tolower(common_with_SEQC))
clin_ori <- clin

exprs <- exprs %>%  dplyr::select(!tolower(common_with_SEQC))
exprs_ori <- exprs

######Correlation with INSS stages, MYCN amplification and Age

# transpose expression data
exprs <- as.data.frame(t(exprs))

exprs <-rownames_to_column(exprs,var = "Patient ID")


exprs <- exprs %>%  dplyr::select("Patient ID", RUVBL1, RUVBL2) %>% 
      dplyr::mutate_at(c("RUVBL1", "RUVBL2"), as.numeric)

# join clin data and exprs data
clin <- left_join(clin, exprs, by="Patient ID")

# remove na from INSS_stage  in data
clin <- clin %>%  dplyr::filter(!INSS_stage=="na")

# genes of interest
genes <- c("RUVBL1","RUVBL2")

# correlation with INSS_stage

inss_cor_list_cangelosi <- list()

for(i in 1:length(genes)){
      p1 <- ggplot(data=clin, mapping=aes(x=INSS_stage, y=.data[[genes[i]]]))+
            geom_violin(fill="#AF0171", color="#AF0171")+theme_bw()+
            #geom_jitter(width = 0.1)+
            geom_boxplot(width=0.1, outlier.shape = NA)+
            
            labs(x="INSS stage", y=paste0("log2",genes[i]), title = "NB-Cangelosi")+
            theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
                  plot.title = element_text(size = 8, face = "bold",hjust = 0.5))#+
      scale_y_continuous(breaks = seq(-3,2,by=0.5))
      
      my_comparisons <- list(c("1","2"), c("1","3"), c("1","4"), c("4", "4s"))
      
      p2 <- p1 + stat_compare_means(comparisons = my_comparisons, label.y = c(5.9,6.1,6.8,6.55), label = "p.format")
      
      inss_cor_list_cangelosi[[i]] <- p2
}


(inss_cor_list_cangelosi[[1]]|inss_cor_list_cangelosi[[2]])


### Correlation with MYCN status

# reset clinical data
clin <- clin_ori

# join clin data and exprs data

clin <- left_join(clin, exprs, by="Patient ID")

clin <- clin %>%  dplyr::filter(!MYCN_status=="NA")

MYCN_status_cor_list_cangelosi <- list() # initialize list

for( i in 1:length(genes)){
      
      p3 <- ggplot(data = clin, mapping = aes(x=MYCN_status, y=clin[[genes[i]]]))+
            theme_bw()+
            geom_violin(fill="#790252", color="#790252")+
            #geom_jitter(width = 0.1)+
            geom_boxplot(width=0.1, outlier.shape = NA)+
            labs(x="MYCN status", y=paste0("log2",genes[i]), title = "NB-Cangelosi")+
            theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
                  plot.title = element_text(size = 8, face = "bold",hjust = 0.5))#+
      #scale_y_continuous(breaks = seq(-3,2,by=0.5))
      
      p4 <- p3 + stat_compare_means(comparisons = list(c("Amplified" , "Not Amplified")), 
                                    label.y = 6.7, label = "p.format") 
      
      MYCN_status_cor_list_cangelosi[[i]] <- p4
      
}

(MYCN_status_cor_list_cangelosi[[1]]|MYCN_status_cor_list_cangelosi[[2]])


##### correlation with Age

# plot

# join clin data and exprs data
clin <- clin_ori
clin <- left_join(clin, exprs, by="Patient ID")

Age_corr_plot_list_Cangelosi <- list()

for( i in 1:length(genes)){
      
      p3 <- ggplot(data = clin, mapping = aes(x=as.factor(Age_group), y=.data[[genes[i]]]))+
            theme_bw()+
            geom_violin(fill="steelblue", color="steelblue")+
            #geom_jitter(width = 0.1)+
            geom_boxplot(width=0.1, outlier.shape = NA)+
            labs(x="Age at diagnosis", y=paste0("log2",genes[i]), title = "Cangelosi")+
            theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
                  plot.title = element_text(size = 8, hjust = 0.5))
      
      #label_y <- c(10.3,14) # specify the position of p_value labels on plot
      
      # statistical analysis
      p4 <- p3 + stat_compare_means(comparisons = list(c("<18" , ">=18")), 
                                    label = "p.format") 
      
      Age_corr_plot_list_Cangelosi[[i]] <- p4
      
}

(Age_corr_plot_list_Cangelosi[[1]]|Age_corr_plot_list_Cangelosi[[2]])


##### correlation between RUVBL1 and RUVBL2

# reset exprs data
exprs <- exprs_ori
# extract RUVBL1 and RUVBL2 expr data from main data
df <- exprs %>% filter(row.names(exprs) %in% c("RUVBL1", "RUVBL2")) 
df <- as.data.frame(t(df)) %>% 
      remove_rownames() %>% 
      mutate_all(as.numeric)

# 
RUV1_vs_RUV2_corr_plt <- ggscatter(df, x="RUVBL1", y="RUVBL2", conf.int = TRUE, 
                                   add="reg.line", cor.coef = TRUE, 
                                   color = "steelblue", size = 0.5)+
      #cor.coef.coord=c(0.5,-2))+
      theme_bw()+ theme(axis.text= element_text(size=8),axis.title = element_text(size = 8, face="italic"))

RUV1_vs_RUV2_corr_plt


