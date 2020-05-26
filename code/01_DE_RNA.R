
library(DESeq2)
library(tximport)
library(tidyverse)
library(biomaRt)
source('code/sox9_themes.R')

# Import design data
design <- read_csv('data/rna/biliary_sox9_agracz_samples.csv')
design$path <- file.path('results/quant/', design$Sample, 'quant.sf')

# Setup mart 
mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
mart_res <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_id', 'external_gene_name'), mart = mart)
mart_res <- mart_res[,c(1,3)]

# Import quant files
txi <- tximport(design$path, type = 'salmon', tx2gene = mart_res, ignoreTxVersion = T)
head(txi$counts)
colnames(txi$counts) <- design$Sample
colnames(txi$counts)
design$Group <- relevel(as.factor(design$Group), 'Low')

# Function for reading out results
conde <- function(obj) { 
  obj <- obj %>% as.data.frame() %>% rownames_to_column()  
  return(obj)
}

test_dds <- DESeqDataSetFromMatrix(round(txi$counts) , colData = design, design = ~ Group) 
test_dds_expanded <- DESeq(test_dds, modelMatrixType = 'expanded', betaPrior = T)
resultsNames(test_dds_expanded)

res_low <- results(test_dds_expanded, contrast=list("GroupLow",
                           c('GroupHigh', 'GroupNeg', 'GroupSL') ), 
                          listValues=c(1, -1/3)) 
res_hi <- results(test_dds_expanded, contrast=list("GroupHigh",
                           c('GroupLow', 'GroupNeg', 'GroupSL') ), 
                          listValues=c(1, -1/3)) 

res_sl <- results(test_dds_expanded, contrast=list("GroupSL",
                           c('GroupHigh', 'GroupNeg', 'GroupLow') ), 
                          listValues=c(1, -1/3)) 

res_neg <- results(test_dds_expanded, contrast=list("GroupNeg",
                           c('GroupHigh', 'GroupLow', 'GroupSL') ), 
                          listValues=c(1, -1/3))

# Summary/plot information

#summary(res_hi)
#summary(res_low) 
#summary(res_sl) 
#summary(res_neg) 

#plotMA(res_hi)
#plotMA(res_low)
#plotMA(res_sl) 
#plotMA(res_neg) 

# Save the significant results
dir.create('results/DESeq_results/')
write_csv(x = conde(res_low), 'results/DESeq_results/low_all.csv')
write_csv(x = conde(res_hi), 'results/DESeq_results/hi_all.csv')
write_csv(x = conde(res_sl), 'results/DESeq_results/sl_all.csv')
write_csv(x = conde(res_neg), 'results/DESeq_results/neg_all.csv')

#Save the transformed counts
norm_counts <- counts(test_dds_expanded, norm = T) 
m <- t(scale(t(as.matrix(norm_counts)), center = T, scale = T) )
write_csv(conde(m), 'results/DESeq_results/zscored_counts_norm.csv')
write_csv(conde(norm_counts), 'results/DESeq_results/norm_counts.csv') 
save(test_dds, file = 'results/DESeq_results/test_dds.R')  

#PCA of samples
vst <- varianceStabilizingTransformation(test_dds_expanded)  
pc_plot <- plotPCA(vst, intgroup = 'Group', returnData = T)   #PC 1 = 55% PC2 = 9%
plt <- pc_plot %>%
  mutate(Group = ifelse(Group == 'SL', 'Sub', as.character(Group))) %>% 
  mutate(Group = factor(Group, levels  = c("High", "Low", "Sub", 'Neg') ) ) %>% 
  ggplot(aes( x= PC1, y = PC2, fill = Group)) + 
  geom_point(size = 5, pch = 21) + 
  theme_sox9() + 
  theme(axis.text = element_text(size= 18), 
        axis.title = element_text(size = 20) ) + 
  scale_fill_manual (values = rev(sox9_cols)) + 
  xlab("PC1 55%") + ylab('PC2 9%') 
ggsave(plt, 'results/DESeq_results/pca_plot.pdf')
sessionInfo() 
#plt
