#Downloading required packages 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

install.packages("survminer")
install.packages("survival") 

if (!require("BiocManager", quitely = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

install.packages("tidyverse")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#loading the packages 

library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

#Getting the clincial data from TCGA-BCRA cohort 

brca_data<- GDCquery_clinic("TCGA-BRCA")
any(colnames(brca_data) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(brca_data) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
brca_data[,c(9,39,45)]

#variables associated with survival 

table(brca_data$vital_status)

#changing values 

brca_data$deceased <- ifelse(brca_data$vital_status == "Alive" , FALSE , TRUE)
table(brca_data$deceased)


brca_data$overall_survival <- ifelse(brca_data$vital_status == "Alive", 
                                     brca_data$days_to_last_follow_up,
                                     brca_data$days_to_death)


#query for getting gene expression data for cohort

query_brca_all = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling" ,
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts" ,
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

output_brca<- getResults(query_brca_all)
output_brca

#getting 20 primary tissue sample barcodes 
tumor_barcode<- output_brca$cases[1:20]

#query for downloading thw 20 sample cases 

query_brca_samples = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling" ,
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts" ,
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor_barcode)


GDCdownload(query_brca_samples)

#counts

tcga_brca_data<- GDCprepare(query_brca_samples, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, "unstranded")
brca_matrix[1:10,1:10]

#extracting gene & metadata from summarizedExperiemnt from object

gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))

#variance stabilization transformation for survival analysis - countdata object 
dds<- DESeqDataSetFromMatrix(countData = brca_matrix,
                             colData = coldata,
                             design = ~ 1)


#removing genes with sum total of 10 across all reads

select_genes <- rowSums(counts(dds)) >= 10
dds <- dds[select_genes,]

vsd<- vst(dds, blind = FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10,1:10]


#getting metadata info and data for tp53 gene

brca_TP53<-brca_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id' , value = 'counts', -gene_id) %>%
 left_join(., gene_metadata , by = 'gene_id') %>%
filter(gene_name == 'TP53')

#median value 
 median_value <- median(brca_TP53$counts)

 #which cases have higher or lower expression than median count 
 brca_TP53$strata <- ifelse(brca_TP53$counts >= median_value, "HIGH" , "LOW")

 #adding clinial information to brca_TP53
brca_TP53$case_id <- gsub('-01.*','',brca_TP53$case_id)
brca_TP53 <- merge(brca_TP53,brca_data , by.x = 'case_id' , by.y = 'submitter_id')

#survival curve

brca_surv_curve <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_TP53)

ggsurvplot(brca_surv_curve,
           data = brca_TP53,
           pval = T,
           risk.table = T)


brca_curve2 <-survdiff(Surv(overall_survival, deceased) ~ strata, data = brca_TP53)

