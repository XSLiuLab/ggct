library(UCSCXenaTools)
library(dplyr)

mut_query = XenaData %>% 
  filter(XenaDatasets=='TCGA.LUAD.sampleMap/mutation_broad_gene') %>% 
  XenaGenerate() %>% XenaQuery() %>% XenaDownload(destdir = "data")

luad_mut = XenaPrepare(mut_query)
luad_mut$sample[grepl("KRAS", luad_mut$sample)]        

kras_status = luad_mut %>% 
  filter(sample == "KRAS") %>% 
  tidyr::gather(barcodes, status, -sample) %>% 
  select(-sample)

luad_expr = XenaPrepare('HiSeqV2')
ggct_expr = luad_expr %>% 
  filter(sample == "GGCT") %>% 
  tidyr::gather(barcodes, expr, -sample) %>% 
  select(-sample)

merged_df = kras_status %>% 
  left_join(ggct_expr) %>% 
  mutate(KRAS_Status = ifelse(status == 1, 'Mut', 'Non-Mut')) %>% 
  filter(!is.na(expr))

#library(ggplot2)
library(ggpubr)
library(cowplot)

ggboxplot(merged_df, x = 'KRAS_Status', y = 'expr', add = 'mean_se') +
  stat_compare_means(aes(group = KRAS_Status), label.x = 1.5, method = "t.test") + 
  labs(x = 'KRAS Status', y = "GGCT expression (log2(RSEM_normalized_count+1))")

# ggplot(merged_df, aes(x = KRAS_Status, y = expr)) +
#   geom_bar(stat = "identity") + 
#   ggpubr::stat_compare_means(method = "t.test")

# Try gdcHub
expr_query = XenaData %>% 
  filter(XenaDatasets%in%c("TCGA-LUAD/Xena_Matrices/TCGA-LUAD.htseq_fpkm.tsv",
                           "TCGA-LUAD/Xena_Matrices/TCGA-LUAD.htseq_counts.tsv")) %>% 
  XenaGenerate() %>% XenaQuery() %>% XenaDownload(destdir = "data", download_probeMap = TRUE, trans_slash = TRUE)

df_list = XenaPrepare(expr_query)
df_list$probeMaps__gencode.v22.annotation.gene.probeMap.gz = df_list$probeMaps__gencode.v22.annotation.gene.probeMap.gz %>%
  dplyr::select(id, gene)

ggct_id = df_list$probeMaps__gencode.v22.annotation.gene.probeMap.gz %>%
  dplyr::filter(gene == "GGCT") %>% pull(id)
 

ggct_expr_count = df_list$TCGA.LUAD__Xena_Matrices__TCGA.LUAD.htseq_counts.tsv.gz %>%
  rename(sample = Ensembl_ID) %>%
  filter(sample == ggct_id) %>% 
  tidyr::gather(barcodes, expr, -sample) %>% 
  select(-sample) %>% 
  mutate(barcodes = substr(barcodes, 1, 15))

ggct_expr_fpkm = df_list$TCGA.LUAD__Xena_Matrices__TCGA.LUAD.htseq_fpkm.tsv.gz %>%
  rename(sample = Ensembl_ID) %>%
  filter(sample == ggct_id) %>% 
  tidyr::gather(barcodes, expr, -sample) %>% 
  select(-sample) %>% 
  mutate(barcodes = substr(barcodes, 1, 15))

merged_df2 = kras_status %>% 
  left_join(ggct_expr_count) %>% 
  mutate(KRAS_Status = ifelse(status == 1, 'Mut', 'Non-Mut')) %>% 
  filter(!is.na(expr))

merged_df3 = kras_status %>% 
  left_join(ggct_expr_fpkm) %>% 
  mutate(KRAS_Status = ifelse(status == 1, 'Mut', 'Non-Mut')) %>% 
  filter(!is.na(expr))

# Plots
ggboxplot(merged_df, x = 'KRAS_Status', y = 'expr', add = 'mean_se') +
  stat_compare_means(aes(group = KRAS_Status), label.x = 1.5, method = "t.test") + 
  labs(x = 'KRAS Status', y = "GGCT expression (log2(RSEM_normalized_count+1))")

ggboxplot(merged_df2, x = 'KRAS_Status', y = 'expr', add = 'mean_se') +
  stat_compare_means(aes(group = KRAS_Status), label.x = 1.5, method = "t.test") + 
  labs(x = 'KRAS Status', y = "GGCT expression (log2(count+1))")

ggboxplot(merged_df3, x = 'KRAS_Status', y = 'expr', add = 'mean_se') +
  stat_compare_means(aes(group = KRAS_Status), label.x = 1.5, method = "t.test") + 
  labs(x = 'KRAS Status', y = "GGCT expression (log2(fpkm+1))")

summary(merged_df$expr[merged_df$KRAS_Status=='Mut'])
summary(merged_df$expr[merged_df$KRAS_Status=='Non-Mut'])

summary(merged_df3$expr[merged_df3$KRAS_Status=='Mut'])
summary(merged_df3$expr[merged_df3$KRAS_Status=='Non-Mut'])
