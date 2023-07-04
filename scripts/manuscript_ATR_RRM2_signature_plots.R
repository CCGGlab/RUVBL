#### ATR and RRM2 signature

library(enrichplot)
library(clusterProfiler)
library(tidyverse)

#load ATR and RRM2 signature/genesets created from Szydzik et al,2021 and Nunes et al,2022 respectively
ATR_RRM2 <- read_rds("data/ATR_RRM2_sig_geneset_clusterPformat.rds")

## load CB-6644 neuroblastoma DEseq2 data from clb-ba _72hrs
# data
res_DEseq_list_sknas_clb_ba <- read_rds("data/DGE.rds" )
# extract dataset for clb-ba at 72hrs only
res_DEseq_ba_72h <- res_DEseq_list_sknas_clb_ba[[2]][[3]] %>% as.data.frame(.)

##### Running GSEA
# make named genelist, vector
genelist <- res_DEseq_ba_72h %>%  dplyr::select(HGNC, stat) %>% 
      deframe(.) %>% 
      sort(., decreasing = TRUE) # rank genelist

# GSEA
res_gsea_ATR_RRM2 <- GSEA(genelist, TERM2GENE = ATR_RRM2, eps=0)



## make barplot

# get result of ATR and RRM2 separately as df
# ATR
res_df_ATR <- res_gsea_ATR_RRM2@result %>%  dplyr::filter(str_detect(ID,"ATR_up|ATR_down"))
res_df_ATR$Description <- factor(res_df_ATR$Description, levels = res_df_ATR$Description[rev(order(res_df_ATR$p.adjust))])
# RRM2
res_df_RRM2 <- res_gsea_ATR_RRM2@result %>%  dplyr::filter(str_detect(ID,"RRM2_up|RRM2_down"))
res_df_RRM2$Description <- factor(res_df_RRM2$Description, levels = res_df_RRM2$Description[rev(order(res_df_RRM2$p.adjust))])

# plot
res_df <- list(ATR=res_df_ATR,RRM2=res_df_RRM2)
sig_names <- c("ATR_signature","RRM2_signature")
authors <- c("_Szydzik","_Nunes")

mytheme <- theme(plot.title = element_text(size=7), axis.title = element_text(size = 7),
                 axis.text.y = element_text(size = 6),
                 axis.text.x = element_text(size = 5)) 

bxplts <- pmap(list(res_df, sig_names,authors), function(res,name,author){
      plt <- ggplot(res, mapping = aes(x=-log10(p.adjust), y=Description , color=Description, fill=Description))+
            geom_bar(stat = "identity",width=0.5)+
            theme_bw()+ theme(axis.text= element_text(size=10),axis.title.x = element_text(size =10), axis.title.y = element_blank())+
            ggtitle(paste0("GSEA (", name, author,")"))+ 
            theme(legend.position = "none")+
            xlab("-log10(Padj)")+mytheme
      
})

bxplts[[1]]|bxplts[[2]]

##  make running score plot
# get result of ATR and RRM2 separately as gsea objects
# ATR
res_gsea_ATR <- res_gsea_ATR_RRM2 %>%  dplyr::filter(str_detect(ID,"ATR_up|ATR_down"))
# RRM2
res_gsea_RRM2 <- res_gsea_ATR_RRM2 %>%  dplyr::filter(str_detect(ID,"RRM2_up|RRM2_down"))

# plot ATR
ATR_sigs <- c("ATR_up", "ATR_down")
run_score_plt_ATR <- map(ATR_sigs, function(sig)
      gseaplot(res_gsea_ATR, title = paste(sig),geneSetID=sig,by = "runningScore")+mytheme
)

# plot RRM2
RRM2_sigs <- c("RRM2_up", "RRM2_down")
run_score_plt_RRM2 <- map(RRM2_sigs, function(sig)
      gseaplot(res_gsea_RRM2, title = paste(sig),geneSetID=sig,by = "runningScore")+mytheme
)


### patch plots
# ATR
ATR_plt <-(bxplts[[1]]|(run_score_plt_ATR[[1]]|run_score_plt_ATR[[2]]))+plot_layout(widths = c(1,2))

# RRM2
RRM2_plt <- (bxplts[[2]]|(run_score_plt_RRM2[[1]]|run_score_plt_RRM2[[2]]))+plot_layout(widths = c(1,2))

(ATR_plt/RRM2_plt)


