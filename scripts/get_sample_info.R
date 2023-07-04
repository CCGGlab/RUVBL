# Sample info
###############

sample_info<- read.table("raw//sample_info.csv", sep=";",header=T, stringsAsFactors = F)

sampleNames<- sample_info$run_id
CL_tmp<- sample_info$names
treatment_tmp<- "CB66"
time_tmp<-  gsub(".*_","",sample_info$names_full)
time_tmp<- gsub("ctrl","0",time_tmp)
time_tmp<- gsub("hr","",time_tmp)
conc_tmp<- "250nM"
repl_tmp<- rep(1:3,8)

sample_info<- data.frame(
  row.names = sampleNames,
  sampleNames = sampleNames,
  CL = CL_tmp,
  treatment = treatment_tmp,
  concentration = conc_tmp,
  time = time_tmp,
  replicate = repl_tmp
)

# Save
saveRDS(sample_info,file="data/sample_info_cb6644_rnaseq.rds")

