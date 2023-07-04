
######## Correlate RUVBL1 and RUVBL2 with all genes in independent cohorts of NB tumors 
######## rank genes based on their correlation coefficient
######## run GSEA using the Hallmark gene sets 

## data source: Davide Cangelosi et al., 2020 (PMID: 32825087)

library(tidyverse)
library(ggrepel)
library(msigdbr)
library(clusterProfiler)
library(patchwork)

### Load tumor exprs data
# Cangelosi (use the complete Cangelosi (integrated platforms) dataset which contains about 360 samples from SEQC)
exprs <- read_rds("data/exprs_data_Cangelosi_formated.rds")

### Correlate RUVBL1 and RUVBL2 with all genes

# transpose data
exprs_df <- as.data.frame(t(exprs)) 

# expr data of genes of interest
exprs_goi <- exprs_df %>% dplyr::select(RUVBL1,RUVBL2)

# generate correlation analysis function
cor_fn <- function(x) map_dfr(1:ncol(exprs_df), function(vec){
      # correlation analysis
      res <- cor.test(as.numeric(unlist(x)),as.numeric(exprs_df[[vec]]),method = "pearson",alternative="two.sided")
      #extract results of cor.test into dataframe
      df <- data.frame(corr_coefficient=res$estimate,corr_pvalue=res$p.value,t_statistic=res$statistic)
})

# correaltion analysis
res_cor <- map(exprs_goi, cor_fn) %>%
      # assign gene names
      map(function(res) {rownames(res) <- colnames(exprs_df)
      res <- rownames_to_column(res, var = "GeneSymbol")
      res })


############# Gene Set Analysis
## Running GSEA
# rank genes based on t_statistic of the correlation analaysis
genelist_list <- res_cor %>%  map(function(corr_df) corr_df %>% dplyr::select(GeneSymbol,corr_coefficient) %>%
                                             deframe(.)  %>% sort(., decreasing = TRUE))
## Hallmark genesets
# extract hallmark geneset from msigdb
H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
      dplyr::select(gs_name, gene_symbol)

# remove HALLMARK_ prefix from geneset name
H_t2g$gs_name <- gsub("HALLMARK_", "", H_t2g$gs_name)

# GSEA ---- result is dataframe
res_gsea_df_list <- genelist_list %>%  map(function(genelist){
      res <- GSEA(genelist, TERM2GENE = H_t2g, eps=0)
      res <- res@result # return result as dataframe
})


res_gsea_com <- intersect(res_gsea_df_list[[1]]$Description, res_gsea_df_list[[2]]$Description)


##### make volcano plot
# pathways to highlight in volcano plot
pathways_highlight <-  c("DNA_REPAIR","E2F_TARGETS","G2M_CHECKPOINT","GLYCOLYSIS","MTORC1_SIGNALING",
                         "MYC_TARGETS_V1","MYC_TARGETS_V2", "OXIDATIVE_PHOSPHORYLATION","REACTIVE_OXYGEN_SPECIES_PATHWAY",
                         "APOPTOSIS", "TNFA_SIGNALING_VIA_NFKB","INFLAMMATORY_RESPONSE", 
                         "IL6_JAK_STAT3_SIGNALING","	P53_PATHWAY" )
# labels
label_gs_list <- res_gsea_df_list %>% map(function(res_gsea) res_gsea %>% dplyr::filter(ID %in% pathways_highlight ))

genes <- c("RUVBL1", "RUVBL2")

# volcano
volcPlt_list <- map2(res_gsea_df_list,genes, function(res_gsea, gene)
      res_gsea %>% 
            ggplot(aes(x=NES, y=-log10(p.adjust)))+
            geom_point(size=3)+
            theme_bw()+ xlab("normalized enrichemnt score")+
            xlim(-3,4)+ ylim(0,60)+
            ggtitle(paste("Preranked GSEA (Hallmark geneset)",gene, sep = "-"), subtitle = "Cangelosi (integrated platforms)")+
            theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
                  plot.title = element_text(size = 8, face = "bold",hjust = 0.5),
                  plot.subtitle =element_text(size = 7,face = "bold",hjust = 0.5))
)

# add labels to  points on plot
volcPlt_list <- map2(volcPlt_list, label_gs_list, function(volcPlt, label_gs)
      volcPlt + geom_point(data = label_gs,shape = 21,size = 4,fill = "#712B75",colour = "black") +  # add labels
            geom_label_repel(data = label_gs, # Add labels last to appear as the top layer
                             aes(label = Description),color = 'black',size=2, force = 4,nudge_y = 5, 
                             nudge_x = 0.05,box.padding = 0.5,color="black", size=2.5, fill=NA))

####### running score and preranked list of GSEA result

# GSEA ---- result is gsea object (to be used to make score plpots)
# gene sets common to RUVBL1 and RUVBL2 GSEA
res_gsea_list <- genelist_list %>%  map(function(genelist)
      GSEA(genelist, TERM2GENE = H_t2g, eps=0)
)
# list of pathways to plot (running score plot)
selected_pathways <- c("MYC_TARGETS_V1", "E2F_TARGETS", "G2M_CHECKPOINT", "APOPTOSIS" )
mytheme <- theme(plot.title = element_text(size=7), axis.title = element_text(size = 7),
                 axis.text.y = element_text(size = 6),
                 axis.text.x = element_text(size = 5)) 

## running score plot 
# runing score plot funtion
run_scoreplot_fn <- function(res_gsea) map(selected_pathways, function(pathway){
      # res_gsea is a result of GSEA run (GSEA object), selected_pathways are pathways to highlight
      plt <- gseaplot(res_gsea, geneSetID = pathway, title = paste(pathway),by = "runningScore")+ mytheme
      plt
})

# plots
score_plts_list <- map(res_gsea_list,run_scoreplot_fn)

# subset plots for each gene
score_plt_RUVBL1 <- score_plts_list[[1]]
score_plt_RUVBL2 <- score_plts_list[[2]]

# patch plots
plot_RUVBL1 <- ((volcPlt_list[[1]]|((score_plt_RUVBL1[[1]]|score_plt_RUVBL1[[2]])/
                                         (score_plt_RUVBL1[[3]]|score_plt_RUVBL1[[4]])))+plot_layout(width=c(2,1)))
plot_RUVBL1

plot_RUVBL2 <- (volcPlt_list[[2]]|((score_plt_RUVBL2[[1]]|score_plt_RUVBL2[[2]])/
                                        (score_plt_RUVBL2[[3]]|score_plt_RUVBL2[[4]])))+plot_layout(width=c(2,1))
plot_RUVBL2



