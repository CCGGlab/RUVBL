### Aim: To  performe single sample's gene expression profile of Hallmark geneset signature  in NB cell lines
# data source: Depmap: CCLE 

library(singscore)
library(GSEABase)
library(SummarizedExperiment)
library(edgeR)
library(rtracklayer)
library(plyr)
library(tidyverse)
library(msigdbr)
library(patchwork)

### load in gene expression count and sample_info data
counts <- read_rds("data/counts_NB_Depmap_CCLE.rds")
pdata <- read_rds("data/sample_info_NB_Depmap_CCLE.rds")


# make DepMapID column into rownames on both the pdata and counts object.
pdata <- column_to_rownames(pdata, "DepMap_ID")
counts <- column_to_rownames(counts, "DepMap_ID") # DepMap_ID

# order rows of each data
pdata <- pdata[ order(row.names(pdata)), ]
counts <- counts[order(row.names(counts)),]

# make SummarizedExperiment object
CL_counts.se <- SummarizedExperiment(t(counts), colData = pdata)

###### Filter out genes with low counts with egdeR
# make a DGEList (need for edgeR)
data_dge = DGEList(counts = assay(CL_counts.se), genes = rowData(CL_counts.se))

# filter on the CPMs rather than raw counts as the former 
prop_expressed = rowMeans(cpm(data_dge) > 1)
keep = prop_expressed > 0.5

op = par(no.readonly = TRUE)
par(mfrow = c(1, 2))
hist(cpm(data_dge, log = TRUE), main = 'Unfiltered', xlab = 'logCPM')
abline(v = log(1), lty = 2, col = 2)
hist(cpm(data_dge[keep, ], log = TRUE), main = 'Filtered', xlab = 'logCPM')
abline(v = log(1), lty = 2, col = 2)

# subset the data
data_dge = data_dge[keep, , keep.lib.sizes = FALSE]
CL_counts.se = CL_counts.se[keep, ]

##### Transformation to FPKM values and normalisation
# Depmap CCLE_expression gene expression count data has been obtained after alignment
# using STAR and quantified using RNA-SeQC 2, which quantifies reads mapping to the exons of each gene,
# therefore, effective gene lengths can be calculated as the sum of all exons spanning the gene
# The GENCODE v34 annotation file was used

# download v22 of the GENCODE annotation
gencode_file = 'gencode.v34.annotation.gtf.gz'
gencode_link = paste(
      'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34',
      gencode_file,
      sep = '/'
)
download.file(gencode_link, gencode_file, method = 'libcurl')
gtf = import.gff(gencode_file, format = 'gtf', genome = 'GRCh38.106', feature.type = 'exon')

#split records by gene to group exons of the same gene
grl =GenomicRanges::reduce(split(gtf, elementMetadata(gtf)$gene_name))
gene_lengths = ldply(grl, function(x) {
  #sum up the length of individual exons
  return(c('gene_length' = sum(width(x))))
}, .id = 'gene_name')


#extract information on gene biotype
genetype = unique(elementMetadata(gtf)[, c('gene_name', 'gene_type')])

colnames(genetype)[1] = 'gene_name'
gene_lengths = merge(genetype, gene_lengths, by="gene_name")

# allocate rownames for ease of indexing
rownames(gene_lengths) = gene_lengths$gene_name
rowData(CL_counts.se)$gene_length = gene_lengths[rownames(CL_counts.se), 'gene_length']
rowData(CL_counts.se)$gene_biotype = gene_lengths[rownames(CL_counts.se), 'gene_type']

#annotate gene lengths for the DGE object
data_dge$genes$length = gene_lengths[rownames(data_dge), 'gene_length']

## Calculate RPKM
# normalization factor
data_dge_tmm = calcNormFactors(data_dge, method = 'TMM')

#compute FPKM values and append to assays
assay(CL_counts.se, 'logFPKM_TMM') = rpkm(data_dge_tmm, log = TRUE)

####### Transcriptional signatures
# download signature and read into GeneSet objects from the GSEABase R package

# get Hallmark genesets 
hallmark_gs<- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)
hallmark_names <- hallmark_gs %>% dplyr::select(gs_name) %>% distinct(.) %>% pull(.)

#generate URLs
hallmark_link = paste0(
      'http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=',
      hallmark_names,
      '&fileType=xml'
)

setwd("/home/joachim/projects/RUVBL2/temp/data/genesets")

#download files
hallmark_files <- paste0(hallmark_names, '.xml')
mapply(download.file, hallmark_link, hallmark_files, method = 'libcurl') 

hallmark_sigs = getBroadSets(hallmark_files, membersId = 'MEMBERS_SYMBOLIZED')

setwd("/home/joachim/projects/RUVBL2")

### Sample scoring
# rank genes in expr set
rankData <- rankGenes(CL_counts.se)

# score samples for hallmark geneset signatures
scoredf_list <- hallmark_sigs %>%  map(function(gs)simpleScore(rankData, upSet = gs,knownDirection=FALSE))

## create an empty dataframe to record scores 
scoredf = data.frame(matrix(ncol = 0, nrow =31)) 

# extract signature score for each geneset from scoredf_list and bind together into dataframe
scores_H <- map2(scoredf_list, hallmark_names, function(s_df,name) {
      score <- s_df %>%  dplyr::select(TotalScore)
      names(score) <- name
      scoredf <- bind_cols(scoredf, score) 
      scoredf
}) %>% flatten() %>%  
      as.data.frame(.) 

# covert to long data format
score_H_long <- pivot_longer(scores_H, cols = everything(),names_to = "geneset", values_to = "score")

# plot geneset median scores in neuroblastoma CLs 
plot_score = ggplot(score_H_long, aes(x=reorder(geneset, score, FUN=median ), y=score, fill=geneset)) +
      geom_boxplot(position = 'dodge', alpha = 0.6) +
      ylab("signature score")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
            legend.position = "none", axis.title.x  = element_blank())

plot_score


