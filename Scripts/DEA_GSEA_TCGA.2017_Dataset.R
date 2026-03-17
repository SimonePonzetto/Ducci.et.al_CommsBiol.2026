#Bioinformatician:: Simone Ponzetto [simone.ponzetto@ifom.eu]

##DEA and GSEA Analysis - TCGA Dataset downloaded from https://www.cancer.gov [BLCA Samples]
#The comparison of biological interest is Mutate (6 cell lines) vs WT (11 cell lines) (regarding FGFR3 gene)

#for analysis full reproducibility, all packages, dependencies and versions are
#tracked using renv, to re-run the analysis please restore the environment with
#the renv.lock file

##R-Version details:
#platform       x86_64-pc-linux-gnu         
#arch           x86_64                      
#os             linux-gnu                   
#system         x86_64, linux-gnu           
#status                                     
#major          4                           
#minor          5.2                         
#year           2025                        
#month          10                          
#day            31                          
#svn rev        88974                       
#language       R                           
#version.string R version 4.5.2 (2025-10-31)
#nickname       [Not] Part in a Rumble  

################################################################################
##Loading libraries-------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(variancePartition)
library(reformulas)
library(sva)

##Loading Raw count data and Annotations----------------------------------------
#Raw data
rawdata <- read.csv("Raw/TCGA/bladder_dataset_TCGA.csv", sep = ";")
rawdata[2:409] <- sapply(rawdata[2:409], as.numeric)

#Annotation data
rawdata_annot <- readxl::read_xlsx("Raw/TCGA/annotation_bladder_TCGA.xlsx")
rawdata_annot$Sample <- gsub("-", ".", rawdata_annot$Sample)

##Group Counts from same gene---------------------------------------------------
rawdata_grouped <- rawdata %>%
  group_by(X) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

  ##Prepare Design Matrix-------------------------------------------------------
design <- data.frame(Sample = colnames(rawdata_grouped)[-1])
design <- left_join(design, rawdata_annot, by = c("Sample" = "Sample"))
rownames(design) <- design$Sample

## Evaluate variance explained by variable of interest--------------------------
exprSet <- ExpressionSet(assayData = as.matrix(rawdata_grouped[,-1]))

# Model formula (random effects)
form <- ~ (1|FGFR3)
varPart <- fitExtractVarPartModel(exprSet, form, design,
                          BPPARAM=SerialParam())

plotVarPart(varPart)

##Plot PCA to visualize batch effects------------------------------------------
pca <- prcomp(t(rawdata_grouped[,-1]))
pcaData <- as.data.frame(pca$x)
pcaData$FGFR3 <- design$FGFR3
percentVar <- round(100 * attr(pcaData, "percentVar"))  
ggplot(pcaData, aes(PC1, PC2, color=FGFR3)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

############################## Commentary ######################################
## Since the VPA and PCA revealed a moderate variance explanatory power for FGFR3 as
## variable of interest, a Surrogate Variable Analysis shall be performed to
## mitigate the masking effect for other confounding factors so to extract 
## clean and meaningful results from DEA and GSEA.

##Prepare Deseq2 object--------------------------------------------------------
count.matrix <- as.matrix(rawdata_grouped[,-1])
## Ensure counts fit 32-bit integer range (DESeq2 requires integers)
max_val <- max(count.matrix, na.rm = TRUE)
if(max_val > .Machine$integer.max){
  factor <- ceiling(max_val / (.Machine$integer.max * 0.9))
  count.matrix <- floor(count.matrix / factor)
}
# Coerce to integer and validate
count.matrix <- matrix(as.integer(count.matrix), nrow = nrow(count.matrix),
                       dimnames = dimnames(count.matrix))
if(any(is.na(count.matrix))){
  stop("NA values present in count.matrix after coercion; inspect raw data before proceeding")
}
rownames(count.matrix) <- rawdata_grouped$X

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                              colData = design,
                              design = ~FGFR3)

##Perform SVA to estimate hidden batch effects-------------------------------
#Extract dds object counts
vsd <- vst(dds, blind=TRUE)
expr_mat <- assay(vsd)

#Create the required models
mod <- model.matrix(~ FGFR3, data = colData(dds))
mod0 <- model.matrix(~ 1, data = colData(dds))

#Estimate number of surrogate variables
n.sv <- num.sv(expr_mat, mod, method="leek")
svobj <- sva(expr_mat, mod, mod0, n.sv=n.sv)

# Add SVs to dds object
for (i in 1:n.sv) {
  colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
}

## Update design to include SVs-------------------------------------------------
sv_formula <- paste0("~", paste(paste0("SV", 1:n.sv), collapse = " + "), "+ FGFR3")
design(dds) <- as.formula(sv_formula)
dds$FGFR3 <- relevel(dds$FGFR3, ref = "wt")

## Filter dds and run DESeq2----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 56
dds <- dds[keep,]
dds <- DESeq(dds, test = c("Wald"), fitType = c("parametric"))

resultsNames(dds)

res <- results(dds, name = "FGFR3_mut_vs_wt", pAdjustMethod = "BH", test = "Wald")
res2 <- lfcShrink(dds, res = res, coef = "FGFR3_mut_vs_wt", type = c("apeglm"))
res2 <- as.data.frame(res2)

norm_counts <- counts(dds, normalized = T)
norm_counts <- as.data.frame(norm_counts)

res2$Gene <- rownames(res2)
norm_counts$Gene <- rownames(norm_counts)

Final <- merge(res2, norm_counts, by = "Gene", all = T)

Final$padj <- ifelse(is.na(Final$padj), 1, Final$padj)
Final$rnk <- Final$log2FoldChange*(-log10(Final$padj))

write.csv(Final, "../FGFR3_MutvsWT_TCGADataset_SVs_Corrected_11012026.csv", row.names = T)

##Create Volcano Plot-----------------------------------------------------------
volcano <- res2
volcano[is.na(volcano$padj),"padj"] <- 1
volcano$Logpval <- -log10(volcano$padj)

volcano <- volcano %>%
  mutate(direction = case_when(log2FoldChange >= 1.0 & Logpval > 1 ~ "UpReg",
                               log2FoldChange <= -1.0 & Logpval > 1.0 ~ "DownReg",
                               TRUE ~ "NotSign"))

topgenes <- volcano[volcano$Logpval >= 10 & abs(volcano$log2FoldChange) >= 2,]

pdf(file = "../FGFR3_Mut_vs_Wt_TCGA.BLCA_Volcano.Plot.pdf", width = 10, height = 8)
ggplot(volcano, aes(x = log2FoldChange, y = Logpval, colour = direction, label = Gene)) +
  geom_point() +
  geom_text_repel(data = topgenes, force = 4, colour = "black", segment.alpha = 0.4) +
  scale_color_manual(values = c("blue4", "grey", "red4")) +
  geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
  geom_hline(yintercept = 1, colour = "black", linetype = "dashed") +
  xlab("Log2 (FoldChange)") +
  ylab("Log10 (Adj. p-Value)") +
  labs(title = "Volcano Plot", subtitle = "FGFR3 Mut vs Wt - TCGA-BLCA",
       colour = "Mode of Regulation") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16))
dev.off()

##Extract and plot PCA----------------------------------------------------------
dds.trsf <- vst(dds)
dds.pca <- plotPCA(dds.trsf, returnData = T, intgroup = "FGFR3")

pdf(file = "../FGFR3_Mut_vs_Wt_TCGA.BLCA_PCA.Plot.pdf", width = 10, height = 8)
ggplot(dds.pca, aes(x = PC1, y = PC2, colour = FGFR3)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#9100d2", "#d9a111")) +
  labs(title = "PCA Plot", subtitle = "FGFR3 Mut vs Wt - TCGA-BLCA",
       colour = "Condition") +
  xlab("PC1: 23%") +
  ylab("PC2: 11%") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16))
dev.off()
