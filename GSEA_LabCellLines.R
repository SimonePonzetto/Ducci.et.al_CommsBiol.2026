#Bioinformatician:: Simone Ponzetto [IFOM]

##GSEA Analysis - in Lab cell lines [RT112, 5637, UMUC3, HT1376, J82]
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
##Loading libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
