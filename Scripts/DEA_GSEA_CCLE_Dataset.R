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
library(reformulas)
library(variancePartition)
library(sva)

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
indx.match <- match(colnames(rawdata.filt[2:18]), rawdata_annot$`Cell line`)
rawdata_annot <- rawdata_annot[indx.match,]

designMatrix <- data.frame(
  row.names = colnames(rawdata.filt[2:18]),
  sample = colnames(rawdata.filt[2:18]),
  Invasiveness = rawdata_annot$Invasiveness,
  FGFR3 = rawdata_annot$FGFR3
)

## Evaluate variance explained by variable of interest--------------------------
exprSet <- ExpressionSet(assayData = as.matrix(rawdata.filt[2:18]))

# Model formula (random effects)
form <- ~ (1|FGFR3) + (1|Invasiveness)
varPart <- fitExtractVarPartModel(exprSet, form, designMatrix,
                                  BPPARAM=SerialParam())

plotVarPart(varPart)

############################## Commentary ######################################
## Since the VPA revealed a moderate variance explanatory power for FGFR3 as
## variable of interest, a Surrogate Variable Analysis shall be performed to
## mitigate the masking effect for other confounding factors so to extract 
## clean and meaningful results from DEA and GSEA.
## Invasiveness, since a co-variate of interest and since its explained variance
## surpass FGFR3 will be kept as random effect for modelling.

## I dds object creation--------------------------------------------------------
count.matrix <- as.matrix(rawdata.filt[2:18])
rownames(count.matrix) <- rawdata.filt$gene_name
dds <-  DESeqDataSetFromMatrix(countData = count.matrix,
                               colData = designMatrix, design = ~Invasiveness + FGFR3)

## Perform SVA------------------------------------------------------------------
# Extract the dds object
vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

# Create the required models
mod <- model.matrix(~ Invasiveness + FGFR3, data = colData(dds))    
mod0 <- model.matrix(~ 1, data = colData(dds))

# Estimate optimal SVs
n.sv <- num.sv(expr_mat, mod, method = "be")  
svobj <- sva(expr_mat, mod, mod0, n.sv = n.sv) 

# Add SVs to dds object
for (i in 1:n.sv) {
  colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
}

## Update design to include SVs-------------------------------------------------
sv_formula <- paste0("~", paste(paste0("SV", 1:n.sv), collapse = " + "), "+ Invasiveness + FGFR3")
design(dds) <- as.formula(sv_formula)
dds$FGFR3 <- relevel(dds$FGFR3, ref = "wt")

## Filter dds and run DESeq2----------------------------------------------------
keep <- rowSums(counts(dds) >= 20) >= 4
dds <- dds[keep,]
dds <- DESeq(dds, test = c("Wald"), fitType = c("parametric"))

resultsNames(dds)

res <- results(dds, name = "FGFR3_mut_vs_wt", pAdjustMethod = "BH", test = "Wald")
res2 <- lfcShrink(dds, res = res, coef = "FGFR3_mut_vs_wt", type = c("ashr"))
res2 <- as.data.frame(res2)

norm_counts <- counts(dds, normalized = T)
norm_counts <- as.data.frame(norm_counts)

res2$Gene <- rownames(res2)
norm_counts$Gene <- rownames(norm_counts)

Final <- merge(res2, norm_counts, by = "Gene", all = T)
Final$rnk <- Final$log2FoldChange*(-log10(Final$pvalue))

write.csv(Final, "../FGFR3_MutvsWT_CCLEDataset_SVs_Corrected_07122025.csv", row.names = T)

##Create Volcano Plot-----------------------------------------------------------
volcano <- res2
volcano[is.na(volcano$padj),"padj"] <- 1
volcano$Logpval <- -log10(volcano$padj)

volcano <- volcano %>%
  mutate(direction = case_when(log2FoldChange >= 1.0 & Logpval > 1 ~ "UpReg",
                               log2FoldChange <= -1.0 & Logpval > 1.0 ~ "DownReg",
                               TRUE ~ "NotSign"))

topgenes <- volcano[volcano$Logpval >= 3 & abs(volcano$log2FoldChange) >= 2,]

ggplot(volcano, aes(x = log2FoldChange, y = Logpval, colour = direction, label = Gene)) +
  geom_point() +
  geom_text_repel(data = topgenes, force = 4, colour = "black", segment.alpha = 0.4) +
  scale_color_manual(values = c("blue4", "grey", "red4")) +
  geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
  geom_hline(yintercept = 1, colour = "black", linetype = "dashed") +
  xlab("Log2 (FoldChange)") +
  ylab("Log10 (Adj. p-Value)") +
  labs(title = "Volcano Plot", subtitle = "FGFR3 Mut vs Wt - CCLE Dataset",
       colour = "Mode of Regulation") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16))

##Extract and plot PCA----------------------------------------------------------
dds.trsf <- vst(dds)
dds.pca <- plotPCA(dds.trsf, intgroup = "sample", returnData = T)

ggplot(dds.pca, aes(x = PC1, y = PC2, colour = FGFR3, label = sample)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#9100d2", "#d9a111")) +
  geom_text_repel(colour = "black", segment.alpha = 0.4, force = 4) +
  labs(title = "PCA Plot", subtitle = "FGFR3 Mut vs Wt - Cell Lines",
       colour = "Condition") +
  xlab("PC1: 33%") +
  ylab("PC2: 11%") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16))
