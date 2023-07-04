# Rpacks
#########

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(pkgs=c(
  "ggbreak",
  "survivalAnalysis", # use the  cox_as_data_frame function to export coxph results
  "WriteXLS"
))


