# Plot DESeq results as Heatmap
library(tidyverse)
library(ComplexHeatmap)
library(UpSetR)
source('code/sox9_themes.R')


# Import DESeq results (01_DE_RNA.R)r
res_hi  <- read_csv('results/DESeq_results/hi_all.csv')
res_low <- read_csv('results/DESeq_results/low_all.csv')
res_sl  <- read_csv('results/DESeq_results/sl_all.csv')
res_neg <- read_csv('results/DESeq_results/neg_all.csv')

# Get list of all signfiicant 
hi_sig <- res_hi %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0) 
low_sig <- res_low %>% filter(padj < 0.05)  %>% filter(log2FoldChange > 0 ) 
sl_sig  <- res_sl %>% filter(padj < 0.05)  %>% filter(log2FoldChange > 0 )
neg_sig <- res_neg %>% filter(padj < 0.05)  %>% filter(log2FoldChange > 0 )

all_sig <- union(hi_sig$rowname, union(low_sig$rowname, union(sl_sig$rowname, neg_sig$rowname) ) ) 

#Filter to those uniquely up-regulated
hi_up_only <- hi_sig %>% filter(!rowname %in% union(low_sig$rowname, union(sl_sig$rowname, neg_sig$rowname) ) )  %>% filter(log2FoldChange > 0) 
low_up_only <- low_sig %>% filter(!rowname %in% union(hi_sig$rowname, union(sl_sig$rowname, neg_sig$rowname) ) ) %>% filter(log2FoldChange > 0) 
sl_up_only <- sl_sig %>% filter(!rowname %in% union(hi_sig$rowname, union(low_sig$rowname, neg_sig$rowname) ) )  %>% filter(log2FoldChange > 0 ) 
neg_up_only <- neg_sig %>% filter(!rowname %in% union(hi_sig$rowname, union(low_sig$rowname, sl_sig$rowname) ) ) %>% filter(log2FoldChange > 0 ) 

#Import z-scored centered count data from normalized count data (01_DE_RNA.R)
m <- read_csv('results/DESeq_results/zscored_counts_norm.csv')
mn <- m$rowname
m <- as.matrix(m[,2:ncol(m)])
rownames(m) <- mn
m_sub <- m[rownames(m) %in% all_sig, ]

genes <- rownames(m_sub)
samples <- data.frame(samples = colnames(m_sub)[1:ncol(m_sub)]) 
samples <- samples %>% separate(samples, into = c('pop', 'rep'), sep = '_', remove = F) 
samples <- samples %>% 
  mutate(pop = factor(pop, levels = c('High', 'Low', 'SL', 'Neg') ) )
samples <- samples   %>% arrange(pop)

# Now filter the scaled data by those in our ordering and use only those genes found significantly upregulated in only one population
up.single.gx <- m_sub[rownames(m_sub) %in% c(hi_up_only$rowname, low_up_only$rowname, sl_up_only$rowname, neg_up_only$rowname), ]  
up.single.gx2 <- up.single.gx[, match(samples$samples, colnames(up.single.gx)) ] # arrange ordering

# combine replicates by averaging
up.avg <- as.data.frame(up.single.gx2) %>% 
  rownames_to_column() %>% 
  gather(sample, val, -rowname) %>% 
  separate(sample, into = c('Pop', 'Rep'), sep = '_', remove = F)  %>%
  group_by(Pop, rowname) %>%
  summarise(u = mean(val, na.rm = T) ) %>%
  spread(rowname, u)
up.groups <- up.avg$Pop
up.avg <- up.avg[,2:ncol(up.avg)]
up.avg <- t(up.avg) 
up.avg <- up.avg[,c(1,2,4,3)]
up.avg <- as.matrix(up.avg) 
colnames(up.avg) <- c('Hi', 'Low', 'SubLow', 'Neg') 

up.avg <- up.avg %>% as.data.frame() %>% dplyr::select(Neg, SubLow, Low, Hi) 

#Generate a dataframe for splitting the heatmap. 
row.splits <- data.frame(rownames = rownames(up.single.gx2) )
hi.up.only.df <- data.frame(gene = hi_up_only$rowname, group = 4) 
low.up.only.df <- data.frame(gene = low_up_only$rowname, group = 3) 
sl.up.only.df <- data.frame(gene = sl_up_only$rowname, group = 2) 
neg.up.only.df <- data.frame(gene = neg_up_only$rowname, group = 1) 
comb.uponly <- rbind(hi.up.only.df, low.up.only.df, sl.up.only.df, neg.up.only.df) 
row.splits <- inner_join(row.splits, comb.uponly, by = c('rownames' = 'gene') ) 
#head(row.splits)

# Change group names from numbers to real population names
row.splits <- row.splits %>%
  mutate(newgroup = case_when(group == 3 ~ 'Low', group == 4 ~ 'Hi' , group == 1 ~ 'Neg', group ==2 ~'SL')) %>%
  mutate(newgroup = factor(newgroup, levels= c('Neg', 'SL', 'Low', 'Hi')))

# Addition of genes to label on heatmap
# Add column names from addam
label_genes <- read_tsv('data/NEW_Sox9BEC_HeatmapGenes_05142020.txt', col_names = c('Gene_name', 'pop') )  
label_genes <- label_genes[label_genes$Gene_name %in% rownames(up.avg),   ]
index <- unlist(lapply(label_genes$Gene_name, function(x) as.character(which(rownames(up.avg) == x ) )) ,  )
label_genes$index <- index
label_genes$index <- as.numeric(label_genes$index)

# Plot the heatmap
ha <- HeatmapAnnotation(pop = rev(c('Hi', 'Low', 'SL', 'Neg') ), 
                        which = 'column', 
                        col = list(pop =  c("Hi" = sox9_cols[4], "Low" = sox9_cols[3], "SL" = sox9_cols[2], "Neg" = sox9_cols[1]) )  ) 

ha_anno <- rowAnnotation(genes = anno_mark(at = label_genes$index, labels = label_genes$Gene_name) ) 
# Heatmap rows are not actually arranged here - they are just split by the groupings - this looks much better than arranging within group which never works quite right in my hands
mycols <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100) 
myBreaks <- seq(-1.5, 1.5, length.out = 100)
hm <- Heatmap(as.matrix(up.avg) , 
              show_row_names = F, 
              cluster_columns = F, 
              show_column_names = F, 
              col = circlize::colorRamp2(myBreaks, colors = mycols), 
              cluster_rows = F, 
              split = row.splits$newgroup, 
              show_row_dend = F, 
              top_annotation = ha, 
              right_annotation = ha_anno) 

#hm
pdf(file= 'results/plots/heatmap_diff_genes.pdf', height = 12, width = 4, onefile = T, family = 'Helvetica') 
hm
dev.off()

# Write out the data frames for these genes in this heatmap . 
write_csv(hi_up_only, 'results/DESeq_results/hi_up_only.csv')
write_csv(low_up_only, 'results/DESeq_results/low_up_only.csv')
write_csv(sl_up_only, 'results/DESeq_results/sl_up_only.csv')
write_csv(neg_up_only, 'results/DESeq_results/neg_up_only.csv') 


# Upset plot of upregulated genes in each population
# Import up-regulated genes (can be upregulated in more than one population)
hi_up <- read_csv('results/DESeq_results/hi_all.csv') %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0 ) %>% mutate(group = 'Hi')    # n = 832
low_up <- read_csv('results/DESeq_results/low_all.csv') %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0 ) %>% mutate(group = 'Low') # n = 803
sl_up  <- read_csv('results/DESeq_results/sl_all.csv') %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0 )  %>% mutate(group = 'SL')  # n = 852
neg_up <- read_csv('results/DESeq_results/neg_all.csv') %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0 ) %>% mutate(group = 'Neg') # n = 1097



#Upset plot of interactions 
all_list <- list(Hi = hi_up$rowname, Low = low_up$rowname, SL = sl_up$rowname, Neg  = neg_up$rowname)  
pdf(file = 'results/plots/upset_plot.pdf', height = 8, width = 12, onefile = T, family = 'Helvetica')
upset(fromList(all_list), keep.order = T, text.scale = 3)  
dev.off()

sessionInfo()
