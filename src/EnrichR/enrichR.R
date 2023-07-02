# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html

library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)
################################################
setwd("~/Scrivania/Bioinformatics/Exam")
source("EnrichR/getEnrichmentPlot.R")
################################################

dirRes = "Results/"
if(!dir.exists(dirRes)){
  dir.create(dirRes)
}else { 
  print(paste("the directory ", dirRes,"already exist"))
}

dataset = "RA" # a folder for my case of study
dirDataset = paste0(dirRes,dataset,"/")
if(!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else { 
  print(paste("the directory ", dirDataset,"already exist"))
}

dirEnrich <- paste0(dirDataset,"Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}
################################################
top_term <- 10 
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirDataset,"DEG.txt")

dbs <- listEnrichrDbs() #lista di tutti i db

dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "DisGeNET")

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")
#input_list <- input_list$genes
input_list$genes <- sub("///.*", "", input_list$genes)

list <- split(input_list$genes, input_list$direction)

df_UP <- enrichr(list$UP, dbs)
df_DOWN <- enrichr(list$DOWN, dbs)

BP_UP <- df_UP$GO_Biological_Process_2021
MF_UP <- df_UP$GO_Molecular_Function_2021
DisGeNET_UP <- df_UP$DisGeNET
KEGG_UP <- df_UP$KEGG_2021_Human

BP_DOWN <- df_DOWN$GO_Biological_Process_2021
MF_DOWN <- df_DOWN$GO_Molecular_Function_2021
DisGeNET_DOWN <- df_DOWN$DisGeNET
KEGG_DOWN <- df_DOWN$KEGG_2021_Human

getEnrichmentPlot(BP_UP,"GO_BP_UP",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(MF_UP,"GO_MF_UP",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(KEGG_UP,"KEGG_UP",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(DisGeNET_UP,"DisGeNET_UP",top_term,thr_pval,dirEnrich)

getEnrichmentPlot(BP_DOWN,"GO_BP_DOWN",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(MF_DOWN,"GO_MF_DOWN",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(KEGG_DOWN,"KEGG_DOWN",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(DisGeNET_DOWN,"DisGeNET_DOWN",top_term,thr_pval,dirEnrich)


