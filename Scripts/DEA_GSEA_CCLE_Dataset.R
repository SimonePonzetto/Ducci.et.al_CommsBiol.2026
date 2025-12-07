#Bioinformatician:: Simone Ponzetto [simone.ponzetto@ifom.eu]

##DEA and GSEA Analysis - CCLE cell lines (downloaded from https://sites.broadinstitute.org/ccle)
#The comparison of biological interest is Mutate (6 cell lines) vs WT (11 cell lines) (regarding FGFR3 gene)

#for analysis full reproducibility, all packages, dependencies and versions are
#tracked using renv, to re-run the analysis please restore the environment with
#the renv.lock file

##R-Version details:
# platform       x86_64-pc-linux-gnu         
# arch           x86_64                      
# os             linux-gnu                   
# system         x86_64, linux-gnu           
# status                                     
# major          4                           
# minor          5.2                         
# year           2025                        
# month          10                          
# day            31                          
# svn rev        88974                       
# language       R                           
# version.string R version 4.5.2 (2025-10-31)
# nickname       [Not] Part in a Rumble  

################################################################################
##Loading libraries-------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(CePa)
library(readxl)
library(EnsDb.Hsapiens.v86)

##Loading Raw count data and Annotations----------------------------------------
#Raw data
rawdata <- read.gct("Raw/CCLE/CCLE_RNAseq_genes_counts_20180929.gct")

#Annotation data
rawdata_annot <- read_xlsx("Raw/CCLE/annotation_tpm_bladder_CCLE ripulito.xlsx", sheet = 1)

#Filter out unwanted cell lines (see the paper for the rationale)
unwanted <- c("253J_URINARY_TRACT", "253JBV_URINARY_TRACT", "639V_URINARY_TRACT",
              "HS172T_FIBROBLAST", "SCABER_URINARY_TRACT", "UMUC1_URINARY_TRACT",
              "KMBC2_URINARY_TRACT", "SW1710_URINARY_TRACT", "VMCUB1_URINARY_TRACT")
rawdata_annot <- rawdata_annot[!(rawdata_annot$`Cell line` %in% unwanted),]

#Filter raw data (rescue manually 5637 and 647V)
rawdata.filt <- rawdata[,colnames(rawdata) %in% rawdata_annot$`Cell line`]
indx <- grep("5637|647V", colnames(rawdata))
rawdata.filt <- cbind(rawdata.filt, rawdata[,indx])

##Polishing column names to match annotation------------------------------------
names <- colnames(rawdata.filt)
names <- gsub("^X(?=\\d)", "", names, perl = TRUE)
colnames(rawdata.filt) <- names
rawdata.filt <- as.data.frame(rawdata.filt)

## Convert Gene IDs and match---------------------------------------------------
edb <- EnsDb.Hsapiens.v86
edb <- genes(edb, columns = c("gene_id", "gene_name", "gene_biotype"))
edb <- as.data.frame(edb)[, c("gene_id", "gene_name", "gene_biotype")]

rawdata.filt$gene_id <- gsub("\\..*$", "", rownames(rawdata.filt))

rawdata.filt <- rawdata.filt %>%
  left_join(edb, by = c("gene_id"= "gene_id"))
rawdata.filt$gene_name <- ifelse(is.na(rawdata.filt$gene_name),
                                 rawdata.filt$gene_id,
                                 rawdata.filt$gene_name)

## Group counts from same gene--------------------------------------------------
rawdata.filt <- rawdata.filt %>%
  dplyr::select(-c(gene_id, gene_biotype)) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum, na.rm = T))

## Prepare design matrix--------------------------------------------------------
