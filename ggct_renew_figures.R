library(TCGAbiolinks)
library(tidyverse)
query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, directory = "data/")
data <- GDCprepare(query, directory = "data/")
data = data %>% 
  select(-GDC_Aliquot) %>% 
  select(Sample, Chromosome, Start, End, Num_Probes, Segment_Mean)
write_tsv(data, path="data/TCGA_LUAD.seg")

data$arm = "q"
data = data %>% 
  select(Sample, Chromosome, arm, Start, End, Num_Probes, Segment_Mean)
colnames(data) = c("sampleID", "chrom", "arm", "start.pos", "end.pos", "n.probes", "mean")

data$chrom = gsub(pattern = "X", replacement = "23", 
                  x = data$chrom, fixed = TRUE)
data$chrom = gsub(pattern = "Y", replacement = "24", 
                  x = data$chrom, fixed = TRUE)

centerloc=sigminer::centromeres.hg38

centerloc$chrom = gsub(pattern = "chr", replacement = "", 
                  x = centerloc$chrom, fixed = TRUE)
centerloc$chrom = gsub(pattern = "X", replacement = "23", 
                      x = centerloc$chrom, fixed = TRUE)
centerloc$chrom = gsub(pattern = "Y", replacement = "24", 
                      x = centerloc$chrom, fixed = TRUE)

data_m = left_join(data, centerloc, by="chrom")


data_res = data_m %>% 
  dplyr::mutate(
    arm = ifelse(end.pos < right.base, "p", "q"),
    chrom = as.integer(chrom),
    start.pos = as.integer(start.pos),
    end.pos = as.integer(end.pos)
  ) %>% 
  dplyr::select(-left.base, -right.base) %>% 
  dplyr::filter(chrom != 24) %>% 
  as.data.frame()

plotFreq(segments=data_res, thres.gain=1)

library(maftools)
all.lesions = "data/TCGA-LUAD-TP-CopyNumber_Gistic2_Level_4-20160128/all_lesions.conf_99.txt"
amp.genes = "data/TCGA-LUAD-TP-CopyNumber_Gistic2_Level_4-20160128/amp_genes.conf_99.txt"
del.genes = "data/TCGA-LUAD-TP-CopyNumber_Gistic2_Level_4-20160128/del_genes.conf_99.txt"
scores.gis = "data/TCGA-LUAD-TP-CopyNumber_Gistic2_Level_4-20160128/scores.gistic"

gistic = readGistic(
  gisticAllLesionsFile = all.lesions,
  gisticAmpGenesFile = amp.genes,
  gisticDelGenesFile = del.genes,
  gisticScoresFile = scores.gis
)

gisticChromPlot(gistic, markBands = NA, fdrCutOff = 1)

plotCBSsegments(cbsFile = "data/TCGA_LUAD.seg")


#------------- CCLE
library(UCSCXenaTools)

host = "https://ucscpublic.xenahubs.net"
dataset_cp = "ccle/CCLE_copynumber_byGene_2013-12-03"
dataset_rpkm = "ccle/CCLE_DepMap_18Q2_RNAseq_RPKM_20180502"
#dataset_rppa = "ccle/CCLE_RPPA_20180123"

ggct_cp = fetch_dense_values(host, dataset_cp, "GGCT")
ggct_rpkm = fetch_dense_values(host, dataset_rpkm, "GGCT", use_probeMap = TRUE)
#ggct_rppa = fetch_dense_values(host, dataset_rppa, "GGCT", use_probeMap = TRUE)

ggct_cp = ggct_cp[, grep("LUNG",colnames(ggct_cp),value = T)]
ggct_rpkm = ggct_rpkm[, grep("LUNG",colnames(ggct_rpkm),value = T)]

ggct.cp = dplyr::tibble(
  name = names(ggct_cp),
  copynumber = as.numeric(ggct_cp),
  Status = ifelse(copynumber>0, "Amplification", "Deletion")
) %>% 
  dplyr::arrange(desc(copynumber))

ggct.rpkm = dplyr::tibble(
  name = names(ggct_rpkm),
  rpkm = as.numeric(ggct_rpkm)
)

library(ggpubr)
ggpubr::ggbarplot(data=ggct.cp, x="name", y="copynumber", 
                  fill="Status",color="Status",
                  xlab="185 Lung cancer cell lines",
                  ylab="GGCT CNV", 
                  font.x=c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"), 
                  font.ytickslab = c(14, "plain", "black")) + 
  rremove("x.text") + rremove("x.ticks")  


ggct = dplyr::inner_join(ggct.cp,
                        ggct.rpkm, by="name")

ggpubr::ggscatter(data = ggct, x="copynumber", y="rpkm",
                  add="reg.line",
                  xlab = "GGCT CNV",
                  ylab = "GGCT mRNA expression (RPKM)",
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.x = -0.5, label.sep = "\n")
                  )
