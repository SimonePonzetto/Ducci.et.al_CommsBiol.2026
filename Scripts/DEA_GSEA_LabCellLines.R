#Bioinformatician:: Simone Ponzetto [simone.ponzetto@ifom.eu]

##DEA and GSEA Analysis - in Lab cell lines [RT112, 5637, UMUC3, HT1376, J82]
#The comparison of biological interest is Mutate vs WT (regarding FGFR3 gene)

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

##Loading Raw count data--------------------------------------------------------
rawdata <- read.csv("Raw/Lab_Cell_Lines/Raw_Matrix.csv")
rawdata <- rawdata[1:37]
rawdata <- rawdata[-c(grep("RT4", colnames(rawdata)))]
rawdata <- rawdata[-c(grep("3D", colnames(rawdata)))]

##Polishing column names--------------------------------------------------------
names <- colnames(rawdata)
names <- gsub("_RNASeq_hg38ReadsPerGene.out", "", names)
names <- gsub("X5637", "5637", names)
names <- gsub("rowname", "Gene", names)

colnames(rawdata) <- names

##Prepare dds design matrix-----------------------------------------------------
#Create the count matrix
rawdata.mtx <- as.matrix(rawdata[2:16])
row.names(rawdata.mtx) <- rawdata$Gene

#Create the design matrix
condition <- ifelse(grepl("UMUC3|5637|HT1376",
                          colnames(rawdata.mtx)), "Wt", "Mut")
Invasiveness <- ifelse(grepl("5637|HT1376|J82|UMUC", colnames(rawdata.mtx)), "MIBC", "NMIBC")

design_matrix <- data.frame(
  row.names = colnames(rawdata.mtx),
  sample = colnames(rawdata.mtx),
  Invasiveness = Invasiveness,
  condition = condition  
)

##Prepare dds object------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = rawdata.mtx,
                              colData = design_matrix,
                              design = ~Invasiveness + condition)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

##Run DEA-----------------------------------------------------------------------
dds$condition <- relevel(dds$condition, "Wt")
dds <- DESeq(dds, test = c("Wald"), fitType = c("parametric"))

resultsNames(dds)

##Extract results---------------------------------------------------------------
res <- results(dds, name = "condition_Mut_vs_Wt",
               pAdjustMethod = "BH", test = "Wald")
res2 <- lfcShrink(dds, res = res, coef = "condition_Mut_vs_Wt", type = c("apeglm"))
res2 <- as.data.frame(res2)

norm_counts <- counts(dds, normalized = T)
norm_counts <- as.data.frame(norm_counts)

res2$Gene <- rownames(res2)
norm_counts$Gene <- rownames(norm_counts)

FinalResultsTable <- merge(res2, norm_counts, by = "Gene", all = T)
FinalResultsTable$padj <- ifelse(is.na(FinalResultsTable$padj), 1, FinalResultsTable$padj)
FinalResultsTable$rnk <- FinalResultsTable$log2FoldChange*(-log10(FinalResultsTable$padj))

write.csv(FinalResultsTable, "../ResultsTable_FGFR3_MUTvsWT_Wald_Parametric_08122025.csv", row.names = F)

## Save pre-rnk file------------------------------------------------------------
rnk.file <- FinalResultsTable[c(1,22)]
i <- order(rnk.file$rnk, decreasing = T)
rnk.file <- rnk.file[i,]
write.table(rnk.file$rnk, row.names = rnk.file$Gene, quote = F, col.names = F,
            file = "../LabCL.rnk", sep = "\t")

## Save reference and universe--------------------------------------------------
write.table(FinalResultsTable$Gene, quote = F, row.names = F, col.names = F,
            file = "../LabCL.Universe.WebGestalt.txt")
write.table(volcano$Gene[!(volcano$direction %in% "NotSign")], quote = F,
            row.names = F, col.names = F,
            file = "../LabCL.Universe.WebGestalt.txt")

##Create Volcano Plot-----------------------------------------------------------
volcano <- res2
volcano[is.na(volcano$padj),"padj"] <- 1
volcano$Logpval <- -log10(volcano$padj)

volcano <- volcano %>%
  mutate(direction = case_when(log2FoldChange >= 1.0 & Logpval > 1 ~ "UpReg",
                               log2FoldChange <= -1.0 & Logpval > 1.0 ~ "DownReg",
                               TRUE ~ "NotSign"))

topgenes <- volcano[volcano$Logpval >= 10 & abs(volcano$log2FoldChange) >= 2,]

ggplot(volcano, aes(x = log2FoldChange, y = Logpval, colour = direction, label = Gene)) +
  geom_point() +
  geom_text_repel(data = topgenes, force = 4, colour = "black", segment.alpha = 0.4) +
  scale_color_manual(values = c("blue4", "grey", "red4")) +
  geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
  geom_hline(yintercept = 1, colour = "black", linetype = "dashed") +
  xlab("Log2 (FoldChange)") +
  ylab("Log10 (Adj. p-Value)") +
  labs(title = "Volcano Plot", subtitle = "FGFR3 Mut vs Wt - Cell Lines",
       colour = "Mode of Regulation") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16))

##Extract and plot PCA----------------------------------------------------------
dds.trsf <- vst(dds)
dds.pca <- plotPCA(dds.trsf, returnData = T)

ggplot(dds.pca, aes(x = PC1, y = PC2, colour = condition, label = sample)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#9100d2", "#d9a111")) +
  geom_text_repel(colour = "black", segment.alpha = 0.4, force = 4) +
  labs(title = "PCA Plot", subtitle = "FGFR3 Mut vs Wt - Cell Lines",
       colour = "Condition") +
  xlab("PC1: 44%") +
  ylab("PC2: 23%") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16))
