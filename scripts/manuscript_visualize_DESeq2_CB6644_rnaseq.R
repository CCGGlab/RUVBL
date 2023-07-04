# To perform GSEA anlysis on bulk CB-6644 treated RNAseq samples

library(enrichplot)
library(ggnewscale)
library(clusterProfiler)
library(msigdbr)
library(Hmisc)
library(scales)
library(tidyverse)
library(openxlsx)
library(patchwork)
library(ggrepel)



#### Load DESeq2 processed data
res_DEseq_list_sknas_clb_ba <- read_rds("data/DGE_cb6644_rnaseq.rds" )

# extract dataset for clb-ba at 72hrs only
res_DEseq_ba_72h <- res_DEseq_list_sknas_clb_ba[[2]][[3]] %>% as.data.frame(.)

##### Running GSEA 
# make named genelist, vector
genelist <- res_DEseq_ba_72h %>%  dplyr::select(HGNC, stat) %>% 
      deframe(.) %>% 
      sort(., decreasing = TRUE) # rank genelist

### compute GSEA on hallmark genesets
# extract hallmark geneset from msigdb
H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
      dplyr::select(gs_name, gene_symbol) %>%  
      dplyr::mutate(gs_name=gsub("HALLMARK_", "", gs_name))

# run GSEA
res_gsea_H <- GSEA(genelist, TERM2GENE = H_t2g, eps=0)
# View(res_gsea_H@result)

##### save results to excel
# make a dataframe to save to excel

res_gsea_H_df <- res_gsea_H@result

# write df to excel workbook . 
# write.xlsx(res_gsea_H_df, file = "~/results/tables/res_gsea_H_BA_72h_jvde.xlsx", firstRow=TRUE,
#            sheetName="CTRL vs CB-6644_72h")



### make vector of gensets to highlight

# geneset to not highlight

unselect <- c("UNFOLDED_PROTEIN_RESPONSE", "COAGULATION","ESTROGEN_RESPONSE_EARLY",
              "CHOLESTEROL_HOMEOSTASIS", "ESTROGEN_RESPONSE_LATE", "ADIPOGENESIS", 
              "EPITHELIAL_MESENCHYMAL_TRANSITION", "PEROXISOME", "MYOGENESIS")

# filter gene sets to highlight in plot
res_gsea_H_select <- res_gsea_H %>% dplyr::filter(!ID %in% unselect) 


# pathways to highlight
selected <- res_gsea_H@result %>% dplyr::filter(!ID %in% unselect) %>% 
      dplyr::pull(ID) 

#### make cnetplot of gsea results

# cnetplot
cnet_plt <- cnetplot(res_gsea_H_select, foldChange=genelist, showCategory = selected, 
                     cex_label_category = 0.45,cex_label_gene = 0.25, color_category = "#C70A80",
                     shadowtext="none", node_label="category")+
      #scale_color_gradient(name="Wald stat")+
      ggtitle("Preranked GSEA (Hallmark geneset)")+
      theme(plot.title = element_text(size = 10, face = "bold",hjust = 0.5),
            plot.subtitle =element_text(size = 9,face = "bold",hjust = 0.5))

## change the size parameter for the geneset category in the cnetplot to reflect
# p-adjust value of each geneset

# # view plot data
# cnet_plt$layers[[2]]$data

# covert padjust to -log10 scale
padj <- -log10(res_gsea_H_select@result$p.adjust)

# re-assign p.adjust values to geneset categories under size parameter in plot data
for(i in 1:length(padj)){
      cnet_plt$layers[[2]]$data$size[i] <- padj[i] 
}

# adjust size of gene nodes in cnetplot
for(i in 1:length(cnet_plt$layers[[3]]$data$size)){
      cnet_plt$layers[[3]]$data$size[i] <- as.numeric(0.001)
}

# scale  node sizes in cnetplots
cnet_plt <- cnet_plt + scale_size(
      name = "-log10(padj)",
      breaks = waiver(),
      labels = waiver(),
      limits = NULL,
      range = c(0.5,10),
      trans = "identity",
      guide = "legend")

### change color of category (geneset) nodes to reflect normalized enrichment score (NES) values

# get NES as named vector
nes_scores<- res_gsea_H@result %>% select(ID, NES) %>% deframe

# change values of color in plot data to NES vlaues
cnet_plt$layers[[2]]$data  <- cnet_plt$layers[[2]]$data %>% mutate(color = nes_scores[name])

# change  aesthetic mapping of new color from discrete to continues(which is expressed as NES under color)
cnet_plt$layers[[2]]$mapping$colour_new <- expr(color)

# set color scales for geneset category nodes

colour_new_idx <- which(cnet_plt$scales$find("colour_new")) # find the index of colour_scale in the plot data
cnet_plt$scales$scales[[colour_new_idx]] <-
      scale_colour_distiller(
            name = "NES",
            palette = "PiYG",
            aesthetics = "colour_new",
            guide = guide_colourbar(available_aes = "colour_new")
      )

# lengend name for gene nodes
colour_new_idx <- which(cnet_plt$scales$find("colour"))
cnet_plt$scales$scales[[colour_new_idx]]$name = "Wald stat"

# plot
cnet_plt


######## running score and preranked list of GSEA result

# list of pathways to plot (running score plot)
selected2 <- c("MYC_TARGETS_V1", "E2F_TARGETS", "G2M_CHECKPOINT", 
               "MTORC1_SIGNALING","P53_PATHWAY","APOPTOSIS" )

mytheme <- theme(plot.title = element_text(size=7), axis.title = element_text(size = 7),
                 axis.text.y = element_text(size = 6),
                 axis.text.x = element_text(size = 5)) 
# running score plot 
score_plt_list <- list()
for(i in 1:length(selected2)){
      plt <- gseaplot(res_gsea_H, geneSetID = selected2[i],
                      title = paste(selected2[i]),
                      by = "runningScore",)+ mytheme
      score_plt_list[[i]] <-   plt
      
}


## patch plots
margin <- theme(plot.margin = unit(c(0,20,0,0), "pt"))
margin2 <- theme(plot.margin = unit(c(0,10,30,0), "pt"))
margin3 <- theme(plot.margin = unit(c(40,0,0,0), "pt"))


cnet_score_plts1 <- (cnet_plt|((score_plt_list[[1]]+
                                    margin2|score_plt_list[[2]])/
                                (score_plt_list[[3]]+margin2|score_plt_list[[4]])/
                                    (score_plt_list[[5]]+margin2|score_plt_list[[6]])))+plot_layout(widths = c(2.5,1))
cnet_score_plts1

cnet_score_plts2 <- cnet_plt  /((score_plt_list[[1]]|score_plt_list[[2]]|score_plt_list[[3]]|score_plt_list[[4]]|
            score_plt_list[[5]]|score_plt_list[[6]]))+margin3 +plot_layout(heights = c(3.3,1))
# plot
cnet_score_plts2

#########################################################

### compute GSEA on REACTOME genesets

# extract REACTOME geneset from msigdb
REACTOME_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME" ) %>% 
      dplyr::select(gs_name, gene_symbol)


# remove REACTOME from geneset description
REACTOME_t2g  <- REACTOME_t2g %>% mutate(gs_name= gsub("REACTOME_", "", gs_name))

res_gsea_REACTOME <- GSEA(genelist, TERM2GENE = REACTOME_t2g, eps=0)



# # save results as excel
# write.xlsx(res_gsea_REACTOME@result, file = "~/results/tables/res_gsea_REACTOME_BA_72h_jvde.xlsx", firstRow=TRUE,
#            sheetName="CTRL vs CB-6644_72h")


## make volcano plot of REACTOME result

# genesets to highlight
gs_select <- c("NONSENSE_MEDIATED_DECAY_NMD","DNA_REPLICATION","S_PHASE","MITOTIC_G1_PHASE_AND_G1_S_TRANSITION",
               "METABOLISM_OF_NUCLEOTIDES","CELL_CYCLE_CHECKPOINTS","DNA_REPAIR","EXTENSION_OF_TELOMERES",
               "DNA_DOUBLE_STRAND_BREAK_REPAIR","CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY",
               "TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES","ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
               "CASPASE_ACTIVATION_VIA_DEATH_RECEPTORS_IN_THE_PRESENCE_OF_LIGAND","TELOMERE_MAINTENANCE",
               "TRANSLATION","MITOTIC_G2_G2_M_PHASES","MITOTIC_G1_PHASE_AND_G1_S_TRANSITION", "MISMATCH_REPAIR",
               "BASE_EXCISION_REPAIR","NUCLEOTIDE_EXCISION_REPAIR")
# geneset labels
label <- res_gsea_REACTOME@result %>%  
      dplyr::filter(ID %in% gs_select ) 

# volcano
volplt <- res_gsea_REACTOME@result %>% 
      ggplot(aes(x=NES, y=-log10(p.adjust)))+
      geom_point(size=3, color="steelblue")+
      theme_bw()+
      xlab("normalized enrichemnt score")+
      #xlim(-3,4)+
      #ylim(0,60)+
      ggtitle("Preranked GSEA (Reactome geneset)")+
      theme(axis.text=element_text(size = 8), axis.title = element_text(size=8), 
            plot.title = element_text(size = 8, face = "bold",hjust = 0.5),
            plot.subtitle =element_text(size = 7,face = "bold",hjust = 0.5))
# label genesets
plt <- volplt + geom_point(data = label, # add labels
                         shape = 21,size = 3,fill = "#712B75",colour = "black")+
      geom_label_repel(data = label, # Add labels last to appear as the top layer
                       aes(label = ID),color = 'black',size=1.7, force = 4,nudge_y = 2, 
                       nudge_x = 0.5,box.padding = 0.5,color="black", size=2.5, fill=NA,
                       max.overlaps = 50)


plt

###############################################################

###### Find transcription factors potentially regulating  genes

### compute ranked GSEA on TFT (transcription factor targets (Msigdb))

# extract TFT geneset from msigdb
TFT_t2g <- msigdbr(species = "Homo sapiens",category ="C3" ) %>% 
      dplyr::select(gs_name, gene_symbol)

res_gsea_TFT <- GSEA(genelist, TERM2GENE = TFT_t2g, eps=0)

# save unfiltered results
write_rds(res_gsea_TFT, "results/tables/res_gsea_TFT_unfilt_jvde.rds")

# filter out unknown TFs and MicroRNAs
res_gsea_TFT <- res_gsea_TFT %>%  dplyr::filter(!grepl("UNKNOWN|MIR", ID)) 

# remove all "_" from TF descriptions. same TF appears multiple times due to different
# binding motifs and hence are annotated to reflect such differences. removing all such 
# annotations will give information on number of specific TFs.

res_gsea_TFT <- res_gsea_TFT %>%  dplyr::mutate(Description=gsub("_.*", "", Description)) %>% 
      dplyr::mutate(Description=gsub("DP[1-2].*", "", Description)) %>%  
      dplyr::mutate(Description=str_replace_all(Description,c("GCCATNTTG"="YY1", # this renames other TFs in Description column with appropriate gene names
                                                              "SGCGSSAAA"="E2F1",
                                                              "CATTGTYY"="SOX9"))) 


# group the TFs by Description, arrange based on padjust values and select the most 
# significant TFs for each group

res_gsea_TFT_filt <- res_gsea_TFT %>%  dplyr::group_by(Description) %>% 
      dplyr::arrange(p.adjust) %>%  
      dplyr::filter(row_number()==1) %>% 
      dplyr::filter(p.adjust<= 0.05) # padj cut-off of 0.05

# save TF result
# write_rds(res_gsea_TFT_filt,"results/tables/res_gsea_TFT_filt.rds" )

## save results to excel

# # write list to excel workbook . This list contains all conditions for both sknas and clb-ba
# write.xlsx(res_gsea_TFT_filt@result, file = "~/results/tables/res_gsea_TFT_filt_BA_72h_jvde.xlsx", firstRow=TRUE,
#             sheetName="CTRL vs CB-6644_72h")

# select 
sig_Tfs <- res_gsea_TFT_filt@result %>% dplyr::select(ID,Description, p.adjust, NES) %>% 
      as.data.frame(.) 

# re-level the pathways
sig_Tfs$Description <- factor(sig_Tfs$Description, levels = sig_Tfs$Description[rev(order(sig_Tfs$p.adjust))])


# create new column to group TFs into MYCs, E2Fs, RUVBL1 and others
sig_Tfs <- sig_Tfs %>%  dplyr::mutate(label=if_else(str_detect(Description,"MYC"),"MYCs",
                                                    if_else(str_detect(Description,"E2F"), "E2Fs",
                                                            if_else(str_detect(Description,"RUVBL1"), "RUVBL1",
                                                                    if_else(str_detect(Description,"P53"), "P53",
                                                                            "others")))))
# set levels 
sig_Tfs$label <- factor(sig_Tfs$label, levels = c("MYCs","E2Fs","RUVBL1","P53","others"))

# make bar plot of top pathway
barplot <- ggplot(sig_Tfs, mapping = aes(x=-log10(p.adjust), y=Description , color=label, fill=label))+
      geom_bar(stat = "identity",width=0.5)+
      theme_bw()+ theme(axis.text= element_text(size=8),axis.title.x = element_text(size = 8), axis.title.y = element_blank())+
      ggtitle("Transcription factor prediction")+ scale_fill_discrete(name = "TF category")

barplot


#####################################################
## Make volcano plot of DESeq2-computed differentially expressed genes

## volcanoe plot for clb-b, 72 hrs only
# extract dataset for clb-ba at 72hrs only

# filter out genes with no gene symbols or those with NA pvalues
res_DEseq_ba_72h <- res_DEseq_ba_72h %>%  dplyr::filter(HGNC!="") %>%  
      dplyr::filter(padj !="NA" )

#Add a new column, DE, to table 
res_DEseq_ba_72h <- res_DEseq_ba_72h %>%  dplyr::mutate(DE="NO")

# if log2Foldchange >= 1 and padj < 0.05, set as "UP" 
# if log2Foldchange <= -1 and padj < 0.05, set as "DN" 
# else "NO"
res_DEseq_ba_72h <- res_DEseq_ba_72h %>%  dplyr::mutate(DE=if_else(log2FoldChange >= 1 &
                                                                         padj < 0.01 & 
                                                                         baseMean >10, "UP" ,
                                                                   if_else(log2FoldChange <= -1 &
                                                                                 padj < 0.01 & 
                                                                                 baseMean >10, "DN", "NO")))
# get number of UP and DN regulated genes
num_DE_genes_ba_72h <- res_DEseq_ba_72h  %>% dplyr::count(DE)
num_DE_genes_ba_72h

# make a new variable for labelling of selected genes
selected_genes <- c("MCM7","MCM6","MCM5","NME1","PDK1","NPM1","CHEK2","TERT",  
                    "PHOX2B","MYCN","RRM2","AURKB","TK1","FAS","GPX3",
                    "NTRK1","LDHA","CDKN1A","GADD45A","HK2","HMGA1","STMN2",  
                    "PAICS")
res_DEseq_ba_72h <- res_DEseq_ba_72h %>%  dplyr::mutate(DElabel=if_else(HGNC %in% selected_genes, HGNC, NA_character_))

# colour theme
mycolors <- c("steelblue", "steelblue", "gray")
names(mycolors) <- c("DN", "UP", "NO")

# genes to highlight
genes_highlight <- res_DEseq_ba_72h %>%  filter(HGNC %in% selected_genes)

# volcano plot
vol_plt_ba_72h <- ggplot(data=res_DEseq_ba_72h, aes(x=log2FoldChange, y=-log10(padj)))+
      geom_point(aes(color=DE), show.legend = F)+
      scale_colour_manual(values = mycolors)+
      geom_point(data = genes_highlight,
                 shape = 21,
                 size = 2, 
                 fill = "firebrick", 
                 colour = "firebrick", show.legend = F)+
      geom_label_repel(data = genes_highlight, # Add labels last to appear as the top layer  
                       aes(label = DElabel),
                       force = 2,
                       nudge_y = 3,
                       box.padding = 0.5,
                       color="black", size=2.5, fill=NA,
                       max.overlaps = 25)+
      scale_x_continuous(breaks = c(-8,-6,-4, -2,0, 2,4,6,8))+
      #scale_y_continuous(limits = c(0,250))+
      scale_colour_manual(values = mycolors) +
      theme_bw()
vol_plt_ba_72h 


############################################

####### make volcano plot of all time points  and cell lines

res_DEseq_list_sknas_clb_ba_filt <-  res_DEseq_list_sknas_clb_ba %>% 
      map_depth(2,function(x) as.data.frame(x) %>%  dplyr::filter(HGNC!="") %>%  
                      dplyr::filter(padj !="NA" ))



res_DEseq_list_sknas_clb_ba_filt <- res_DEseq_list_sknas_clb_ba_filt %>% 
      map_depth(2,function(x) dplyr::mutate(x,DE="NO") %>% #Add a new column, DE, to table
                      # if log2Foldchange >= 1 and padj < 0.05, set as "UP" 
                      # if log2Foldchange <= -1 and padj < 0.05, set as "DN" 
                      # else "NO"
                      dplyr::mutate(DE=if_else(log2FoldChange >= 1 &
                                                     padj < 0.01 & 
                                                     baseMean >10, "UP" ,
                                               if_else(log2FoldChange <= -1 &
                                                             padj < 0.01 & 
                                                             baseMean >10, "DN", "NO"))))


# get number of UP and DN regulated genes in each dataset
num_DE_genes_sknas_ba_list <- res_DEseq_list_sknas_clb_ba_filt %>%  map_depth(2,function(x)as.data.frame(x) %>% 
                                                                                    dplyr::count(DE))
# save result of number fo DE genes in each dataset
##### save results to excel
# library(openxlsx)
# 
# names_dataset <-  c("SK-N-AS_24hrs_vs_0hrs","SK-N-AS_48hrs_vs_0hrs","SK-N-AS_72hrs_vs_0hrs",
#                                 "CLB-BA_24hrs_vs_0hrs","CLB-BA_48hrs_vs_0hrs","CLB-BA_72hrs_vs_0hrs")
# # write list to excel workbook . This list contains all conditions for both sknas and clb-ba
# write.xlsx(flatten(num_DE_genes_sknas_ba_list), file = "results/tables/num_DE_genes_sknas_ba_list_jvde.xlsx", firstRow=TRUE,
#            sheetName=names_dataset)


# list of genes to highlight on volcano plots

selected_genes2 <- c("MCM7","MCM6","MCM5","NME1","PFKFB4","CHEK2",  
                     "MYCN","RRM2","AURKB","TK1","GPX3",
                     "NTRK1","HEXIM1","CDKN1A","GADD45A","HK2","HMGA1","VGF",
                     "LMO3","BIRC3","RACK1","ENO1",
                     "NME1-NME2","CEP162","BRCA2","STMN2", "PAICS")


#make a new variable for labeling of selected genes
res_DEseq_list_sknas_clb_ba_filt<- res_DEseq_list_sknas_clb_ba_filt %>% map_depth(2,function(x) 
      dplyr::mutate(x, DElabel=if_else(HGNC %in% selected_genes2, HGNC, NA_character_)))


# genes to highlight
genes_highlight_list <- res_DEseq_list_sknas_clb_ba_filt %>% map_depth(2,function(x) 
      dplyr::filter(x,HGNC %in% selected_genes2)) 

# reduce list levels 
genes_highlight_list <- flatten(genes_highlight_list)
res_DEseq_list_sknas_clb_ba_filt_flt <- flatten(res_DEseq_list_sknas_clb_ba_filt)

# plot names
plt_names <- c("SK-N-AS:24hrs","SK-N-AS:48hrs","SK-N-AS:72hrs",
               "CLB-BA:24hrs","CLB-BA:48hrs","CLB-BA:72hrs")
# plot
vol_plts_sknas_clb_ba_list <-  purrr::pmap(list(res_DEseq_list_sknas_clb_ba_filt_flt,genes_highlight_list,plt_names),
                                           function(res,genes,plt_name) 
                                                 ggplot(data=as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj), col=DE))+
                                                 geom_point(aes(color=DE), show.legend = F)+
                                                 scale_colour_manual(values = mycolors)+
                                                 geom_point(data = as.data.frame(genes),
                                                            shape = 21,
                                                            size = 2, 
                                                            fill = "firebrick", 
                                                            colour = "firebrick", show.legend = F)+
                                                 geom_label_repel(data = as.data.frame(genes), # Add labels last to appear as the top layer  
                                                                  aes(label = DElabel),
                                                                  force = 2,
                                                                  nudge_y = 3,
                                                                  box.padding = 0.5,
                                                                  color="black", size=2,fill=NA,
                                                                  max.overlaps = 45,
                                                                  label.size=0)+
                                                 scale_x_continuous(breaks = c(-8,-6,-4, -2,0, 2,4,6,8))+
                                                 #scale_y_continuous(limits = c(0,200))+
                                                 scale_colour_manual(values = mycolors) +
                                                 ggtitle(plt_name)+
                                                 theme_bw())  

# patch plots
plts_all_CLs_all_cond <- (vol_plts_sknas_clb_ba_list[[4]]|vol_plts_sknas_clb_ba_list[[5]]|vol_plts_sknas_clb_ba_list[[6]])/
      (vol_plts_sknas_clb_ba_list[[1]]|vol_plts_sknas_clb_ba_list[[2]]|vol_plts_sknas_clb_ba_list[[3]])

plts_all_CLs_all_cond