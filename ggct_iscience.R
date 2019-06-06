# GGCT modification analysis for iScience


# Figure 2: GGCT TCGA Distribution -----------------------------------
library(UCSCXenaTools)
library(dplyr)

tcga = XenaData %>% filter(XenaHostNames == "tcgaHub")

# TCGA clinical info download and clean
xe_cli = tcga %>% 
  XenaGenerate() %>% 
  XenaFilter(filterDatasets = "clinical") %>% 
  XenaQuery() %>% XenaDownload(destdir = "data/TCGA_Clinical", trans_slash = TRUE)
cli_files = dir("data/TCGA_Clinical/")
projects = setdiff(sub("TCGA\\.(.*)\\.sampleMap.*", "\\1", cli_files), c("FPPP", "COADREAD", "GBMLGG", "LUNG"))

# keep valuable columns of clinical datasets
select_cols = c("sampleID", "OS", "OS.time", "OS.unit", "RFS", "RFS.time", "RFS.unit", 
                "age_at_initial_pathologic_diagnosis", "gender", "tobacco_smoking_history",
                "tobacco_smoking_history_indicator", "sample_type", "pathologic_M",
                "pathologic_N", "pathologic_T", "pathologic_stage")

cli_list = XenaPrepare(sapply(projects, function(x) {
  file.path("data/TCGA_Clinical/", grep(paste0("__", x, "_"), cli_files, value = TRUE))
}))

TCGA_Clinical = tibble()
for (i in 1:length(projects)){
  clinical = names(cli_list)[i]
  project = projects[i]
  df = cli_list[[clinical]]
  col_exist = select_cols %in% colnames(df)
  res = tibble()
  if(!all(col_exist)){
    res = df[, select_cols[col_exist]]
    res[, select_cols[!col_exist]] = NA
  }else{
    res = df[, select_cols]
  }
  res$Project = project
  res %>% select(Project, select_cols) -> res
  TCGA_Clinical = bind_rows(TCGA_Clinical, res)
}

length(unique(TCGA_Clinical$sampleID))
table(TCGA_Clinical$sample_type)

# Filter sample type
TCGA_Clinical = filter(TCGA_Clinical, sample_type %in% c("Primary Tumor", "Solid Tissue Normal"))
table(TCGA_Clinical$sample_type)

# Obtain proble values
host = "https://pancanatlas.xenahubs.net"
host2 = names(UCSCXenaTools:::.xena_hosts[2])
dataset.expr = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
dataset.cnv = "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes"
dataset.cnv_gene = "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"

all_samples.expr =  .p_dataset_samples(host, dataset.expr, NULL)
ggct.expr = .p_dataset_probe_values(host, dataset.expr, all_samples.expr, as.list("GGCT"))
ggct.expr = as.numeric(ggct.expr[[2]])

all_samples.cnv = .p_dataset_samples(host2, dataset.cnv, NULL)
ggct.cnv = .p_dataset_fetch(host2, dataset.cnv, all_samples.cnv, as.list("GGCT"))
ggct.cnv = as.numeric(ggct.cnv)

all_samples.cnv_gene = .p_dataset_samples(host2, dataset.cnv_gene, NULL)
ggct.cnv_gene = .p_dataset_fetch(host2, dataset.cnv_gene, all_samples.cnv, as.list("GGCT"))
ggct.cnv_gene = as.numeric(ggct.cnv_gene)

identical(all_samples.cnv, all_samples.cnv_gene)

df.ggct = purrr::reduce(list(
  tibble(
    sampleID=all_samples.expr,
    expr=ggct.expr
  ),
  tibble(
    sampleID=all_samples.cnv,
    copynumber=ggct.cnv,
    copynumber_threshold=ggct.cnv_gene
  ),
  TCGA_Clinical[, c("Project", "sampleID", "sample_type")]
), dplyr::full_join)

library(ggpubr)

ggboxplot(filter(df.ggct, !is.na(Project), !is.na(expr)), x = "Project", y = "expr", 
          color="sample_type", xlab = FALSE, ylab = "GGCT expression (log2 based)") +
  rotate_x_text(angle = 45)

# Filter project which has no enough x observations

valid_projects = 
  df.ggct %>% 
  filter(!is.na(expr)) %>% 
  group_by(Project) %>% 
  summarise(N=length(unique(sample_type)), N_normal=sum(sample_type=="Solid Tissue Normal")) %>% 
  filter(N>1 & N_normal>10) %>% pull(Project)

# Remove KIRC: which GGCT significantly downgraduated
valid_projects = setdiff(valid_projects, "KIRC")

# check
setdiff(unique(df.ggct$Project), valid_projects)
table(df.ggct[!is.na(df.ggct$expr), ]$sample_type, df.ggct[!is.na(df.ggct$expr), ]$Project)


df.ggct2 = df.ggct %>% filter(Project %in% valid_projects)

# for (i in valid_projects) {
#   print(i)
#   t.test(expr ~ sample_type, data = df.ggct2[df.ggct2$Project==i,])
# }

compare_means(expr ~ sample_type, data = df.ggct2, group.by = "Project")

p_expr= ggboxplot(df.ggct2, x = "Project", y = "expr", 
          color="sample_type", xlab = FALSE, ylab = "GGCT expression (log2 based)") +
  rotate_x_text(angle = 45) +
  stat_compare_means(aes(group=sample_type), label="p.signif", method="t.test")

p_cnv = ggdotplot(df.ggct2, x = "Project", y = "copynumber", 
          color="Project", xlab = FALSE,
          ylab = "GGCT copy number ratio (log2 based)", legend="none", binwidth = 0.005) +
  rotate_x_text(angle = 45)  + 
  geom_hline(size=0.2, yintercept = 0, linetype=2)


library(patchwork)
library(cowplot)
p_expr+p_cnv+plot_layout(ncol=1)

## GGCT analysis in LUAD

load("./ggct_survival_data_on_metastasis_and_stages.Rdata")
library(survival)
library(survminer)

ggscatter(ggct_on_stages, x = "cnv", y="expression", add="reg.line", 
          xlab="GGCT copy number ratio (log2 based)", ylab = "GGCT expression (log2 based)") +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.2, size= 6, label.y = 13.5
  ) + scale_y_continuous()


df.early = ggct_on_stages[ggct_on_stages$classes == "early_stage", ]
df.late = ggct_on_stages[ggct_on_stages$classes == "late_stage", ]

fit = survfit(Surv(OS, OS_IND) ~ group_by_exp, data = ggct_on_stages)
ggsurvplot_facet(fit, ggct_on_stages, facet.by = "classes", palette = "jco", pval = TRUE)


ggscatter(ggct_on_stages, x = "cnv", y="expression", add="reg.line", facet.by = "classes",
          xlab="GGCT copy number ratio (log2 based)", ylab = "GGCT expression (log2 based)") +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.2, size= 6, label.y = 13.5
  ) + scale_y_continuous()


ggerrorplot(ggct_on_stages, x="classes", y="expression", 
            desc_stat = "mean_sd", color = "black", xlab=FALSE, ylab="GGCT expression (log2 based)",
            add = "violin", add.params = list(color = "darkgray")) + 
  stat_compare_means(comparisons = list(c("early_stage", "late_stage")), method = "t.test", label = "p.signif")

ggerrorplot(ggct_on_stages, x="classes", y="cnv", 
            desc_stat = "mean_sd", color = "black", xlab=FALSE, ylab="GGCT expression (log2 based)",
            add = "violin", add.params = list(color = "darkgray")) + 
  stat_compare_means(comparisons = list(c("early_stage", "late_stage")), method = "t.test", label = "p.signif")

