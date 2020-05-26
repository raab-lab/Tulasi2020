# Compare expression using single cell data.
# Data from - https://doi.org/10.1016/j.cell.2018.02.001

library(tidyverse)
library(singscore) 
library(ggbeeswarm)
library(tximport)
library(biomaRt)
library(ComplexHeatmap)
library(broom)
source('code/sox9_themes.R')


# Use gene expression siganatures derived from singscore
# Use normalized abundance data from our RNA seq
# Import design data
design <- read_csv('data/rna/biliary_sox9_agracz_samples.csv')
design$path <- file.path('results/quant/', design$Sample, 'quant.sf')

# Setup mart 
mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
mart_res <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_id', 'external_gene_name'), mart = mart)
mart_res <- mart_res[,c(1,3)]

# Import quant files
txi <- tximport(design$path, type = 'salmon', tx2gene = mart_res, ignoreTxVersion = T, countsFromAbundance ='lengthScaledTPM') # Using length scaled TPM from counts
colnames(txi$counts) <- design$Sample
design$Group <- relevel(as.factor(design$Group), 'Low')
ab <- txi$counts
colnames(ab) <- design$Sample

# Do Siganture analysis
# Remove any genes without a rowSum of 10) # thought about removing these, but ultimately left in to match as many genes from old papers as I can.
#ab <- ab[rowSums(ab) > 10, ]#]
# Get rankscores for each cell type
ranks <- rankGenes(ab, tiesMethod = 'min')  # from singscore package

#look at some specific siagnautres form other papers
# notice not all genes in published papers can be matched in our data set, some due to older annotations and some due to excel changing values.
# affects ~20% of the published list

yap <- read_csv('data/DongJ_YAPsignature_Cell_2007.csv')
yap_score <- simpleScore(ranks, upSet = unique(yap$Gene.Symbol)[1:600])   #Using top 600 yap genes as done in other paper
# We did this with the full set as well and it made little difference
#yap_score <- simpleScore(ranks, upSet = unique(yap$Gene.Symbol)) 

# Get Notch ScoreGet Notch ScoreGet Notch Scores
notch <- read_csv('data/VillanuevaA_NotchTargets_Gastro_2012_S7.csv')
notch_up <- notch %>% filter(Direction == 'Up') 
notch_up <- unique(notch_up$`Gene ID`)
notch_down <- notch %>% filter(Direction == 'Down') 
notch_down <- unique(notch_down$`Gene ID`)
notch_score <- simpleScore(ranks, upSet = notch_up, downSet = notch_down)

#Derive Signatures from GSE108315 - using their values for differences between cholangiocytes and hepatocytes
# Schaub 2018
subpops <- read_csv('data/Schaub2018_S1.csv')
subpops <- subpops %>% dplyr::select(`Gene Symbol`, pC1, pC2, pC3, H1, H2, H3, HpC1, HpC2, HpC3, HpC4) 
subpops <- subpops %>%
  gather(sample, score, -`Gene Symbol`)  %>%
  mutate(group = case_when(sample %in% c('pC1', 'pC2', 'pC3') ~ 'PC', 
                           sample %in% c('H1', 'H2', 'H3') ~ 'Hep', 
                           sample %in% c('HpC1', 'HpC2', 'HpC3', 'HpC4') ~ 'HpC') ) %>%
  group_by(group, `Gene Symbol`) %>%
  summarise(u_score  = mean(score, na.rm = T) ) 
subpops <- subpops %>% pivot_wider(names_from = group, values_from = u_score) 

# Now define a cut off for up vs down
# Use +/- 2 for hep +-1
hep_up <- subpops %>% filter(Hep > 1)  #Hep down 
hep_down <- subpops %>% filter(Hep < -1)  #HEP Up
pc_up   <- subpops %>% filter(PC > 1)  #BEC UP
pc_down <- subpops %>% filter(PC < -1)  #BEC down

# Now score Sox9 sampes based on these signatures
hep_score <- simpleScore(ranks, upSet = hep_up$`Gene Symbol`, downSet = hep_down$`Gene Symbol`)
pc_score <- simpleScore(ranks, upSet = pc_up$`Gene Symbol`, downSet = pc_down$`Gene Symbol`) 
# Convert scores to useable df for plotting
convert_scores <- function(score_obj, signature = NULL) { 
  s <- as.data.frame(score_obj) %>%
    rownames_to_column() %>%
    mutate(signature = signature)  %>% 
    separate(rowname, into = c('group', 'rep'), sep = '_', remove = F)  %>%
    mutate(group = factor(group,  levels = rev(c('High', 'Low', 'SL', 'Neg') ) ) )  %>%
    dplyr::select(rowname, group, TotalScore, signature)
  return(s)  
}

# Clean up score DF for each signature
hep_score <- convert_scores(hep_score, 'Hep')
pc_score  <- convert_scores(pc_score, 'PC') 
notch_score <- convert_scores(notch_score, 'Notch') 
yap_score <- convert_scores(yap_score, 'Yap') 

# Combine into one DF
paper_sigs <- rbind(hep_score, pc_score, notch_score, yap_score) 
paper_sigs  <- paper_sigs %>% mutate(group = ifelse(group == 'SL', 'Sub', as.character(group) )) %>%
  mutate(group = factor(group, levels = rev(c('High', 'Low', 'Sub', 'Neg') ) ) ) 

#Plot as line
paper_sigs_mean <- paper_sigs %>%
  group_by(group, signature) %>%
  summarise(u_t = mean(TotalScore), 
            ci = 1.96 * sd(TotalScore)/(sqrt(n() ) ) ) %>%
  mutate(ymin = u_t-ci, ymax = u_t+ci) 

paper_sigs <- paper_sigs %>%
  mutate(signature = factor(signature, levels = c('Hep','PC', 'Notch', 'Yap')))

paper_sigs_plot <- paper_sigs %>%
  ggplot(aes(x = group, y = TotalScore, fill = group) ) +
  geom_point(size = 3, color= 'grey40', shape = 21) + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax, y = NULL), data = paper_sigs_mean, width = 0.1, color = 'grey70')  + 
  geom_point(size = 6, aes(x = group, y = u_t), data = paper_sigs_mean, shape = '_', color = 'grey40')  + 
  theme_sox9() +
  scale_fill_manual(values = sox9_cols)  + 
  facet_wrap(~signature, nrow = 1) 

# Stats
paper_sigs_test <- paper_sigs %>%
  group_by(signature) %>%
  do(tidy(pairwise.t.test(x = .$TotalScore, g = .$group, p.adjust.method = 'BH')))
write_csv(paper_sigs_test, path = 'results/signature_pvals.csv' )

paper_sigs_plot
ggsave('results/plots/paper_signatures_point.pdf', height = 6, width = 12)
