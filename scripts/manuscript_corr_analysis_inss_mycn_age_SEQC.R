### Correlation of RUVBL1 and RUVBL2  with INSS stages, MYCN amplification and Age
### data: NB-SEQC 

library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
library(ggpubr)

## load data
exprs <- read_rds("data/exprs_SEQC.rds")
clin <- read_rds("data/clin_SEQC.rds")

# transpose expression data
exprs <- as.data.frame(t(exprs))
exprs <-rownames_to_column(exprs,var = "GSE")

# select RUVBL1 and RUVBL2 expression data
exprs <- exprs %>%  dplyr::select("GSE", RUVBL1, RUVBL2) %>% 
      dplyr::mutate_at(c("RUVBL1", "RUVBL2"), as.numeric)

# join clin data and exprs data
clin <- left_join(clin, exprs, by="GSE")

###### correlation with INSS_stage
# genes
genes <- c("RUVBL1", "RUVBL2")
# position of p_value labels on plot
label_y <- list(c(9.7,10,10.4,10.3), c(13.9,14.1,14.3,13.8)) 

# plot
inss_cor_plts <-  map2(genes,label_y, function(gene,lab){
      p1 <- ggplot(data=clin, mapping=aes(x=INSS_stage, y=.data[[gene]], group=INSS_stage))+
            geom_violin(fill="#AF0171", color="#AF0171")+theme_bw()+
            #geom_jitter(width = 0.1)+
            geom_boxplot(width=0.1, outlier.shape = NA)+
            
            labs(x="INSS stage", y=paste0("log2",gene), title = "NB-SEQC")+
            theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
                  plot.title = element_text(size = 8, face = "bold",hjust = 0.5))
      my_comparisons <- list(c("1","2"), c("1","3"), c("1","4"), c("4", "4S"))
      p2 <- p1 + stat_compare_means(comparisons = my_comparisons, label.y =lab , label = "p.format")
      p2 })
# name plots in list
names(inss_cor_plts) <- c("RUVBL1", "RUVBL2")

inss_cor_plts[[1]]|inss_cor_plts[[2]]

###### Correlation with MYCN status

# position of p_value labels on plot
label_y <- c(10.3,14) 
# plots
mycn_cor_plts <- map2(genes,label_y, function(gene, lab){
      # filter out "unknown" MYCN category
      clin1 <- clin %>%  dplyr::filter(!MYCN_status=="unknown")
      #plot
      p1 <- ggplot(data = clin1, mapping = aes(x=MYCN_status, y=.data[[gene]]))+
            theme_bw()+
            geom_violin(fill="#790252", color="#790252")+
            #geom_jitter(width = 0.1)+
            geom_boxplot(width=0.1, outlier.shape = NA)+
            labs(x="MYCN status", y=paste0("log2",gene), title = "NB-SEQC")+
            theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
                  plot.title = element_text(size = 8, face = "bold",hjust = 0.5))
      # statistical analysis
      p2 <- p1 + stat_compare_means(comparisons = list(c("Amplified" , "Not Amplified")), 
                                    label.y = lab, label = "p.format") 
      p2
})
names(mycn_cor_plts) <- c("RUVBL1", "RUVBL2")
mycn_cor_plts[[1]]|mycn_cor_plts[[2]]

##### correlation with Age
# plots
age_cor_plt <- map(genes, function(gene){
      #  categorize age
      clin1 <-  clin %>%  dplyr::mutate(age_months=as.numeric(`Age at Diagnosis in Days`)/31) %>% 
            dplyr::mutate(Age_group=case_when(age_months>=18~">18 months", TRUE~"<18 months")) 
      # plot
      p1 <- ggplot(data = clin1, mapping = aes(x=as.factor(Age_group), y=.data[[gene]]))+
            theme_bw()+
            geom_violin(fill="steelblue", color="steelblue")+
            #geom_jitter(width = 0.1)+
            geom_boxplot(width=0.1, outlier.shape = NA)+
            labs(x="Age at diagnosis", y=paste0("log2",gene), title = "NB-SEQC")+
            theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
                  plot.title = element_text(size = 8, hjust = 0.5))
      # statistical analysis
      p2 <- p1 + stat_compare_means(comparisons = list(c("<18 months" , ">18 months")), 
                                    label = "p.format") 
      p2 
})

names(age_cor_plt) <- c("RUVBL1", "RUVBL2")

(age_cor_plt[[1]]|age_cor_plt[[2]])

### correlation between RUVBL1 and RUVBL2
# correlation plot
RUV1_vs_RUV2_corr_plt <- ggscatter(exprs, x="RUVBL1", y="RUVBL2", conf.int = TRUE, 
                                   add="reg.line", cor.coef = TRUE, 
                                   color = "steelblue", size=0.5)+
      theme_bw()+ theme(axis.text= element_text(size=8),
                        axis.title = element_text(size = 8, face="italic"))

RUV1_vs_RUV2_corr_plt
