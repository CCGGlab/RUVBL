### Download and Process transcriptomic dataset of neuroblastoma Th-MYCN mice tumours, at time stages of tumourigenesis

library(limma)
library(Biobase)
library(tidyverse)

# Source of microarray data:
# DOI: 10.3390/cancers13194783

# Read AE file structure
ae <- ArrayExpress::ae2bioc(readRDS("downloads/arrayexpress/E-MTAB-3247/files.rds"))

# Metadata (phenotypic)
pheno <- pData(ae)

meta <- pheno %>%
  select(
    genotype = Factor.Value.genotype.,
    age = Factor.Value.age.,
    sample = Assay.Name
  ) %>%
  magrittr::set_rownames(NULL) %>%
  column_to_rownames("sample")
saveRDS(meta, "data/th_mycn_pheno.rds")

# Process expression data
path <- str_glue_data(pData(ae), "downloads/arrayexpress/E-MTAB-3247/{Array.Data.File}")
expression_raw <- read.maimages(path, source = "agilent", green.only = TRUE)

# Background correction
# Choose minimum method to combine with between arrays normalization
expression_bgcorr <- backgroundCorrect(expression_raw, method="minimum")

# Normalization:
# Between arrays needed (single channel)
# Using quantile normalization
expression_normalized = normalizeBetweenArrays(expression_bgcorr, method = "quantile")

exp_by_probe <- expression_normalized %>%
  # Set rownames in expression data to the probe names.
  # Retrieve the names from genes$ProbeName
  # update_list: apply function to an element of a list
  update_list(values = ~magrittr::set_rownames(E, genes$ProbeName)) %>%
  # Extract values
  pluck("values") %>%
  # Change column names to file names (instead of full paths)
  magrittr::set_colnames(colnames(.) %>% basename)

# -> Probe names converted to gene names after DE analysis
# Platform data of same microarray from GEO
# Agilent-028005 SurePrint G3 Mouse GE 8x60K Microarray (Probe Name version)
platform <- readRDS("downloads/arrayexpress/E-MTAB-3247/platform.rds") %>%
  GEOquery::Table()

exp_by_gene <- exp_by_probe %>%
  as_tibble(rownames = "ID") %>%
  # Add gene names from platform data
  left_join(select(platform, ID, gene = GENE_SYMBOL), by="ID") %>%
  # Remove rows not referring to genes
  filter(gene != "") %>%
  select(-ID) %>%
  select(gene, everything())
  
list(exp_by_probe = exp_by_probe,
     exp_by_gene = exp_by_gene,
     platform = platform) %>%
  saveRDS("data/th_mycn.rds")
