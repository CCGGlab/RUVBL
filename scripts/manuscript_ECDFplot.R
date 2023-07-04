## create an ECDF plot (or Empirical Cumulative Density Function) of results of multivariate Coxph analyses of all genes

library(tidyverse)
library(patchwork)

##  load data
res_coxph_df_seqc <- readr::read_rds("data/res_coxph_all_genes_seqc.rds")
res_coxph_df_cangelosi <-  readr::read_rds("data/res_coxph_all_genes_cangelosi.rds")

## create an ECDF plot (or Empirical Cumulative Density Function) of results of Coxph analyses
## lists of results of Coxph analyses
res_coxph_multi_list <- list(NB_SEQC=res_coxph_df_seqc, Cangelosi=res_coxph_df_cangelosi)

# Empirical Cumulative Distribution Function (CDF) Plots
plts_multi <- map_depth(res_coxph_multi_list,1,function(res_coxph){
  RUVBL2_HR = res_coxph %>% dplyr::filter(factor.name=="RUVBL2") %>% pull(HR)
  RUVBL1_HR = res_coxph %>% dplyr::filter(factor.name=="RUVBL1") %>% pull(HR)
  # create ecdf function
  e = ecdf(res_coxph$HR)
  # plot
  plt <- ggplot(res_coxph, aes(HR)) + 
    stat_ecdf(geom = "point")+
    #label for RUVBL2
    annotate(geom = "label", x = RUVBL2_HR,y = e(RUVBL2_HR), label = "RUVBL2")+
    geom_point(aes(x = RUVBL2_HR ,y = e(RUVBL2_HR)),shape = 21, size = 2, 
               fill = "firebrick",colour = "firebrick", show.legend = F)+
    # label for RUVBL1
    annotate(geom = "label", x = RUVBL1_HR,y = e(RUVBL1_HR), label = "RUVBL1")+
    geom_point(aes(x = RUVBL1_HR ,y = e(RUVBL1_HR)),shape = 21, size = 2, 
               fill = "firebrick", colour = "firebrick", show.legend = F)+
    ggtitle(label = "ECDF plot of hazard ratios of multivariate Coxph analysis ")+
           
    ylab("ECDF")+
    theme_bw()+
    theme(axis.title = element_text(size = 8), axis.text = element_text(size=8), 
          plot.title = element_text(size = 9),plot.subtitle = element_text(size = 8))
  plt })
# patch plots
ecdf_plot_SEQC <- plts_multi[[1]]+plot_annotation(title = "NB-SEQC")
                      
ecdf_plot_Cangelosi <- plts_multi[[2]]+plot_annotation(title = "NB-Cangelosi")


## lists of results of Coxph analyses
res_coxph_multi_list <- list(NB_SEQC=res_coxph_df_seqc, Cangelosi=res_coxph_df_cangelosi)

# Empirical Cumulative Distribution Function (CDF) Plots
plts_multi <- map_depth(res_coxph_multi_list,1,function(res_coxph){
      RUVBL2_HR = res_coxph %>% dplyr::filter(factor.name=="RUVBL2") %>% pull(HR)
      RUVBL1_HR = res_coxph %>% dplyr::filter(factor.name=="RUVBL1") %>% pull(HR)
      # create ecdf function
      e = ecdf(res_coxph$HR)
      # plot
      plt <- ggplot(res_coxph, aes(HR)) + 
            stat_ecdf(geom = "point")+
            #label for RUVBL2
            annotate(geom = "label", x = RUVBL2_HR,y = e(RUVBL2_HR), label = "RUVBL2")+
            geom_point(aes(x = RUVBL2_HR ,y = e(RUVBL2_HR)),shape = 21, size = 2, 
                       fill = "firebrick",colour = "firebrick", show.legend = F)+
            # label for RUVBL1
            annotate(geom = "label", x = RUVBL1_HR,y = e(RUVBL1_HR), label = "RUVBL1")+
            geom_point(aes(x = RUVBL1_HR ,y = e(RUVBL1_HR)),shape = 21, size = 2, 
                       fill = "firebrick", colour = "firebrick", show.legend = F)+
            ggtitle(label = "ECDF plot of hazard ratios of multivariate Coxph analysis ")+
            
            ylab("ECDF")+
            theme_bw()+
            theme(axis.title = element_text(size = 8), axis.text = element_text(size=8), 
                  plot.title = element_text(size = 9),plot.subtitle = element_text(size = 8))
      plt })
## patch plots

plts_multi[[1]]+plot_annotation(title = "NB-SEQC")

plts_multi[[2]]+plot_annotation(title = "NB-Cangelosi")


