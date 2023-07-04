#### Analysis of Dependecy score of RUVBL1 and RUVBL2 compared to median score of all other genes in neuroblastoma

library(tidyverse)
library(ggpubr)

# load data
CRISPR_gene_effect <- read_rds("downloads/DepMap/CRISPR_gene_effect.rds")
sample_info <- read_csv("downloads/DepMap/sample_info.csv")

# add cell lines name to CRISPR_gene_effect data
CRISPR_gene_effect <-  CRISPR_gene_effect %>% rownames_to_column(var ="stripped_cell_line_name" )

# subset cell line and primary disease columns from sample info
sample_info_sub <- sample_info %>%  dplyr::select(stripped_cell_line_name, primary_disease)

# add primary_disease column to CRISPR_gene_effect
CRISPR_gene_effect <- left_join(CRISPR_gene_effect,sample_info_sub, by="stripped_cell_line_name")

# rearrange columns
CRISPR_gene_effect <- CRISPR_gene_effect %>% 
      dplyr::select(stripped_cell_line_name,primary_disease, everything())

# subset only neuroblastoma cell lines from CRISPR_gene_effect
CRISPR_gene_effect_NB <- CRISPR_gene_effect %>%  filter(primary_disease=="Neuroblastoma")
CRISPR_gene_effect_NB <- column_to_rownames(CRISPR_gene_effect_NB, "stripped_cell_line_name")
CRISPR_gene_effect_NB <-CRISPR_gene_effect_NB %>%  dplyr::select(!primary_disease)

# format data
CRISPR_gene_effect_NB <- as.data.frame(CRISPR_gene_effect_NB) %>% 
  rownames_to_column(var = "cell_line")

# covert to longer data format
CRISPR_gene_effect_long <- CRISPR_gene_effect_NB %>% 
  pivot_longer(!cell_line, names_to = "gene_name", 
               values_to = "dependency_score")

# get list of DepMap predefined essential genes
CRISPR_common_essentials <- readRDS("downloads/DepMap/CRISPR_common_essentials.rds")

# remove RUVBL1 and RUVBL2 from list of essential genes
CRISPR_common_essentials_filt <- CRISPR_common_essentials %>%
  dplyr::filter(!gene %in% c("RUVBL1", "RUVBL2"))

# categorise gene names into RUVBL1, RUVBL2, essen_genes and non_essen genes
CRISPR_gene_effect_long <-  CRISPR_gene_effect_long %>%  
  dplyr::mutate(gene_category=case_when(gene_name=="RUVBL1"~"RUVBL1",
                                        gene_name=="RUVBL2"~ "RUVBL2",
                                        gene_name %in% CRISPR_common_essentials_filt$gene~ "essen_gene",
                                        TRUE~"non_essen_gene"))
                                  
# make violin plot
plt<- CRISPR_gene_effect_long %>%  dplyr::group_by(gene_category) %>% 
      dplyr::mutate(gene_category=factor(gene_category, 
                                         levels=c("essen_gene","RUVBL1","RUVBL2","non_essen_gene"))) %>% 
       ggplot(mapping = aes(x=gene_category, y=dependency_score))+
       geom_violin(fill="#EF5B0C", color="#EF5B0C")+
       theme_bw()+
       #geom_jitter(width = 0.01)+
       geom_boxplot(width=0.13, outlier.shape = NA)+
       theme(axis.text=element_text(size = 10), axis.title = element_text(size=10), 
             plot.title = element_text(size = 10, face = "bold",hjust = 0.5))+
       scale_y_continuous(breaks = seq(-4,3,by=1))
 
  
# statistical analysis
my_comparisons <- list(c("essen_gene", "RUVBL1"), c("essen_gene", "RUVBL2"),
                       c("non_essen_gene", "RUVBL1"), c("non_essen_gene", "RUVBL2"),
                       c("non_essen_gene", "essen_gene"))

plt <- plt + stat_compare_means(comparisons = my_comparisons, method = "t.test",label.y = c(2.4,2.55,1.4,1.6,2.2), 
                                  label = "p.signif")
plt
