
### Aim: single-sample gene signature score using Hallmark geneset 
# sample: TARGET NB (phs000467) tumour gene expression data, downloadeded from https://portal.gdc.cancer.gov/projects/TARGET-NBL

library(singscore)
library(GSEABase)
library(SummarizedExperiment)
library(edgeR)
library(rtracklayer)
library(plyr)
library(tidyverse)
library(msigdbr)
library(patchwork)


## load data
counts <- read_rds("/home/joachim/projects/RUVBL2/downloads/TARGET/mRNA/STAR_counts/STAR_counts_unstranded_NB.rds.gz")
clin <- read_rds("/home/joachim/projects/RUVBL2/downloads/TARGET/clin/clin_NBL(161).rds")


## subset clinical data into High-risk and low/intermediate-risk groups
# high-risk
clin_HR <- clin %>% dplyr::filter(`COG Risk Group`=="High Risk") 
# low-risk
clin_LR <- clin %>% dplyr::filter(`COG Risk Group`=="Low Risk" | `COG Risk Group`== "Intermediate Risk") 

# subset count data into  High-risk and low/intermediate-risk groups
counts_HR <- semi_join(counts, clin_HR, by="TARGET USI") %>% 
      group_by(`TARGET USI`) %>% 
      dplyr::filter(row_number()==1) # remove duplicates
counts_LR <- semi_join(counts, clin_LR, by="TARGET USI")%>% 
      group_by(`TARGET USI`) %>% 
      dplyr::filter(row_number()==1) # remove duplicates

# change TARGET USI column into rownames on both the clinical and count data.

clin_HR <- column_to_rownames(clin_HR, "TARGET USI") # metadata
counts_HR <- column_to_rownames(counts_HR, "TARGET USI")  #count data

clin_LR <- column_to_rownames(clin_LR, "TARGET USI") # metadata
counts_LR <- column_to_rownames(counts_LR, "TARGET USI")

# arrange rows of each data
clin_HR <- clin_HR %>%  dplyr::arrange(rownames(.))
counts_HR <- counts_HR %>% dplyr::arrange(rownames(.))
clin_LR <- clin_LR %>%  dplyr::arrange(rownames(.))
counts_LR <- counts_LR %>%  dplyr::arrange(rownames(.))


# make SummarizedExperiment object
counts_HR.se <- SummarizedExperiment(t(counts_HR), colData = clin_HR)
counts_LR.se <- SummarizedExperiment(t(counts_LR), colData = counts_LR)

###### Transformation to FPKM values and normalisation
# The GDC mRNA quantification analysis pipeline measures gene level expression
# with   STAR as raw read counts. These data are generated through this
# by first aligning reads to the GRCh38 reference genome.
# Annotation based on GDC.h38 GENCODE v36 GTF

#download v36 of the GENCODE annotation

gencode_file = 'gencode.v36.annotation.gtf.gz'
gencode_link = paste(
  'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36',
  gencode_file,
  sep = '/'
)
download.file(gencode_link, gencode_file, method = 'libcurl')



gtf = import.gff(gencode_file, format = 'gtf', genome = 'GRCh38.p13', feature.type = 'exon')
#split records by gene to group exons of the same gene
grl =GenomicRanges::reduce(split(gtf, elementMetadata(gtf)$gene_name))
gene_lengths = ldply(grl, function(x) {
#   #sum up the length of individual exons
  return(c('gene_length' = sum(width(x))))
}, .id = 'gene_name')


# #extract information on gene biotype
genetype = unique(elementMetadata(gtf)[, c('gene_name', 'gene_type')])

colnames(genetype)[1] = 'gene_name'
gene_lengths = merge(genetype, gene_lengths, by="gene_name")

# # save gene_length
# write_rds(gene_lengths, "/home/joachim/projects/RUVBL2/downloads/gencode/gene_lengths_gencode.38.p13.rds")

gene_lengths <- read_rds("/home/joachim/projects/RUVBL2/downloads/gencode/gene_lengths_gencode.38.p13.rds")


# allocate rownames for ease of indexing
rownames(gene_lengths) <- gene_lengths$gene_name

# Filter out genes with low counts with egdeR and transform to FPKM

se_list <- list(counts_LR.se, counts_HR.se)
names(se_list) <- c("LR","HR")

se.filt_list <- se_list %>%  map(function(se){
      # make a DGEList (need for edgeR)
      data_dge <-  DGEList(counts = assay(se, genes = rowData(se)))
      # filter on the CPMs rather than raw counts as the former 
      # accounts for variation in library sizes, therefore, is unbiased
      prop_expressed <- rowMeans(cpm(data_dge) > 1)
      keep <- prop_expressed > 0.5
      # op = par(no.readonly = TRUE)
      # par(mfrow = c(1, 2))
      # hist(cpm(data_dge, log = TRUE), main = 'Unfiltered', xlab = 'logCPM')
      # abline(v = log(1), lty = 2, col = 2)
      # hist(cpm(data_dge[keep, ], log = TRUE), main = 'Filtered', xlab = 'logCPM')
      # abline(v = log(1), lty = 2, col = 2)
      
      # subset the data
      data_dge <-  data_dge[keep, , keep.lib.sizes = FALSE]
      se.filt <-  se[keep, ]
      rowData(se.filt)$gene_length <- gene_lengths[rownames(se.filt),'gene_length']
      #annotate gene lengths for the DGE object
      data_dge$genes$length <-  gene_lengths[rownames(data_dge), 'gene_length']
      
      ## Calculate RPKM
      # normalization factor
      data_dge_tmm <- calcNormFactors(data_dge, method = 'TMM')
      
      #compute FPKM values and append to assays
      assay(se.filt, 'logFPKM_TMM') <- rpkm(data_dge_tmm, log = TRUE)
      se.filt
})


# # save summarizedExperiment with count and logFPKM assays
# write_rds(se.filt_list, "/home/joachim/projects/RUVBL2/downloads/TARGET/mRNA/counts_NB_filt.se.list.rds.gz", compress = "gz")
# counts_NB_TCGA.se.list <- read_rds( "/home/joachim/projects/RUVBL2/downloads/TARGET/mRNA/counts_NB_filt.se.list.rds.gz")

###### find hallmark genesets that are highly enriched in neuroblasto tumors
## Transcriptional signatures

# get Hallmark genesets names
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)
hallmark_names <- hallmark_genesets %>% dplyr::select(gs_name) %>% distinct(.)
hallmark_names <- hallmark_names$gs_name

#generate URLs
hallmark_link = paste0(
      'http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=',
      hallmark_names,
      '&fileType=xml')

setwd("/home/joachim/projects/RUVBL2/downloads/mSigDB/genesets_xml")

#download files
hallmark_files <- paste0(hallmark_names, '.xml')
mapply(download.file, hallmark_link, hallmark_files, method = 'libcurl') 

hallmark_sigs = getBroadSets(hallmark_files, membersId = 'MEMBERS_SYMBOLIZED')

# change directory to home
setwd("/home/joachim/projects/RUVBL2")

# check the size of the gene sets
for(i in 1:length(hallmark_sigs)){
  length <-  length(GSEABase::geneIds(hallmark_sigs[[i]]))
  print(paste(names(hallmark_sigs)[i],length,"genes",sep = "-"))
}

### Sample scoring
# extract se objects from list
counts_LR.filt.se <- se.filt_list[[1]] # low risk neuroblastoma
counts_HR.filt.se <- se.filt_list[[2]] # high risk neuroblastoma

# rank genes in expr set
rankData_LR <- rankGenes(counts_LR.filt.se)
rankData_HR <- rankGenes(counts_HR.filt.se)

# score samples for hallmark geneset signatures
rankData_list <- list(rankData_LR,rankData_HR) # list of ranked expression set
score_list <- list() # to store signature scores for each dataset(expression set)
score_list.list <- list() # to store scoredf_list of both datasets
for(i in 1:length(rankData_list)){
      temp_rank_data <- rankData_list[[i]]
      for(j in 1:length(hallmark_sigs)){
            score <- simpleScore(temp_rank_data,  # score samples for each geneset signature
                                 upSet = hallmark_sigs[[j]],knownDirection=FALSE)
            score_list[[j]] <- score
            # print(j)
      }
      score_list.list[[i]] <- score_list
      # print(i)
}

### Put all geneset signature scores into dataframe 
scoredf_list <- list()
# extract signature score for each genset from scoredf_list and bind together into dataframe
for(i in 1:length(score_list.list)){
      score_list <- score_list.list[[i]]
      # create empty dataframe to initialize
      scoredf = data.frame(matrix("", nrow =nrow(score_list[[i]]), ncol=0)) #length(score)=number of samples
      for(j in 1:length(score_list)){
            # extract score for each geneset signature
            score <- as.data.frame(score_list[[j]]$TotalScore)
            # add signature name
            names(score) <- hallmark_names[j]
            # bind all scores
            scoredf <-  bind_cols(scoredf, score) 
            # print(hallmark_names[j])
      }
      
      scoredf_list[[i]] <- scoredf
      # print(i)
}

# extract geneset scores for each dataset (LR and HR)
scoredf_LR <- scoredf_list[[1]]
scoredf_HR <- scoredf_list[[2]]

# convert dataframe to long data format
scoredf_LR <- pivot_longer(scoredf_LR,cols = 1:50, names_to = "geneset", values_to = "score")
scoredf_HR <- pivot_longer(scoredf_HR,cols = 1:50, names_to = "geneset", values_to = "score")


## plot geneset median scores in neuroblastoma

# low-risk tumours
plot_LR = ggplot(scoredf_LR, aes(x=reorder(geneset, score, FUN=median ), y=score, fill=geneset)) +
      geom_boxplot(position = 'dodge', alpha = 0.6) +
      # scale_fill_brewer(palette = 'Set2') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.text.y = element_text(size = 8))+
      ggtitle("Hallmark geneset signature score in low-risk neuroblastoma")+
      ylab("geneset signature score")

plot_LR

# high-risk tumours
plot_HR = ggplot(scoredf_HR, aes(x=reorder(geneset, score, FUN=median ), y=score, fill=geneset)) +
      geom_boxplot(position = 'dodge', alpha = 0.6) +
      # scale_fill_brewer(palette = 'Set2') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.text.y = element_text(size = 8))+
      ggtitle("Hallmark geneset signature score in high-risk neuroblastoma")+
      ylab("geneset signature score")
plot_HR












