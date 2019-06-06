require(tidyverse)
setwd("/Volumes/data/Wangshixiang/GGCT")
chr_cnv = read.table("data/TCGA-LUAD-TP-CopyNumber_Gistic2_Level_4-20160128/broad_values_by_arm.txt", header=T, sep="\t")
head(chr_cnv)
df <- chr_cnv
rownames(df) <- df$Chromosome.Arm
head(df)
df <- df[,-1]
head(df)
df_7p <- df["7p",]
head(df_7p)

summary(as.numeric(df_7p))
no_amp <- as.numeric(df_7p) == 0
id_no_amp <- colnames(df_7p[, no_amp])
length(id_no_amp)

id.no_amp <- gsub(".", "-",substr(id_no_amp, 1, 15), fixed=TRUE)
head(id.no_amp)
write.table(id.no_amp, file="./ids_with_no_7p_amplification.txt", sep="\t", col.names=F, quote=F,row.names=F)

amp <- as.numeric(df_7p) > 0.4
id_amp <- colnames(df_7p[, amp])
length(id_amp)
id.amp <- gsub(".", "-", substr(id_amp, 1, 15), fixed=TRUE)
head(id.amp)

write.table(id.amp, file="./ids_with_7p_amplification_threshold_0.4.txt", sep="\t", col.names=F, quote=F,row.names=F)


cli   <- read_tsv(file="./LUAD_clinicalMatrix")
head(cli)


nchar(cli$sampleID[1])

table(cli$sample_type)
table(cli$new_neoplasm_event_type) # This is the wrong cloumn to identify metastasis status
                                   # We will use pathologic_M instead

table(cli$pathologic_M)

cli <- subset(cli, sample_type=="Primary Tumor")
dim(cli)
table(cli$pathologic_M)

# extract sample IDs from clinical data
cli_non_metastasis_samples <- cli[which(cli$pathologic_M=="M0"),]$sampleID
cli_metastasis_samples     <- cli[which(cli$pathologic_M%in%c("M1", "M1a", "M1b")), ]$sampleID

#check
#head(cli_metastasis_samples)

# patients with non metastasis and 7p amplification
non_metas_7p_amp <- intersect(cli_non_metastasis_samples, id.amp)
length(non_metas_7p_amp)
# patients with non metastasis and 7p no amplification
non_metas_7p_no_amp <- intersect(cli_non_metastasis_samples, id.no_amp)
length(non_metas_7p_no_amp)

# patients with metastasis and 7p amplification
metas_7p_amp <- intersect(cli_metastasis_sample, id.amp)
length(metas_7p_amp)

# patients with metastasis and 7p no amplification
metas_7p_no_amp <- intersect(cli_metastasis_sample, id.no_amp)
length(metas_7p_no_amp)

# generate survival data for non metastasis patients with 7p amplification (value > 0.4) or normal
non_meta.res <- c(non_metas_7p_amp, non_metas_7p_no_amp)
non_meta.res_os <- cli[match(non_meta.res, cli$sampleID),'_OS']
non_meta.res_os_ind <- cli[match(non_meta.res, cli$sampleID),'_OS_IND']
non_meta.df <- data.frame(sampleID=non_meta.res, OS=non_meta.res_os, Event=non_meta.res_os_ind)
non_meta.df$amp_7p <- c(non_meta.df$X_OS_IND[1:length(non_metas_7p_amp)], rep(NA,times=length(non_metas_7p_no_amp)))
non_meta.df$norm_7p <- c(rep(NA,times=length(non_metas_7p_amp)),non_meta.df$X_OS_IND[(length(non_metas_7p_amp)+1):length(non_meta.res)])

dim(non_meta.df)
head(non_meta.df)
res_non_meta <- non_meta.df[!is.na(non_meta.df$X_OS),]
dim(res_non_meta)

write.table(res_non_meta, file="./ids_survival_data_non_metastasis.txt", sep="\t", col.names=T, quote=F,row.names=F)


# generate survival data for distant metastasis patients with 7p amplification or normal
meta.res <- c(metas_7p_amp, metas_7p_no_amp)
meta.res_os <- cli[match(meta.res, cli$sampleID),'_OS']
meta.res_os_ind <- cli[match(meta.res, cli$sampleID),'_OS_IND']
meta.df <- data.frame(sampleID=meta.res, OS=meta.res_os, Event=meta.res_os_ind)
meta.df$amp_7p <- c(meta.df$X_OS_IND[1:length(metas_7p_amp)], rep(NA,times=length(metas_7p_no_amp)))
meta.df$norm_7p <- c(rep(NA,times=length(metas_7p_amp)),meta.df$X_OS_IND[(length(metas_7p_amp)+1):length(meta.res)])

dim(meta.df)
head(meta.df)
res_meta <- meta.df[!is.na(meta.df$X_OS),]
dim(res_meta)

write.table(res_meta, file="./ids_survival_data_metastasis.txt", sep="\t", col.names=T, quote=F,row.names=F)



# load expression and cnv data
exprs <- as.data.frame(read_tsv(file="./HiSeqV2"))
cnv   <- as.data.frame(read_tsv(file="./Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"))

cnv_2   <- as.data.frame(read_tsv(file="./Gistic2_CopyNumber_Gistic2_all_data_by_genes"))

head(exprs)

head(cnv)

head(cnv_2)

# use sample IDs above to filter expression data
rownames(exprs) <- exprs$sample
# check
all(rownames(exprs) == exprs$sample)
# remove sample variable
exprs <- exprs[, -1]
exprs <- exprs[, colnames(exprs)%in%cli$sampleID]
dim(exprs)

rownames(cnv) <- cnv$'Gene Symbol'
head(cnv)
cnv <- cnv[, -1]

rownames(cnv_2) <- cnv_2$'Gene Symbol'
head(cnv_2)
cnv_2 <- cnv_2[, -1]

share_ids <- intersect(colnames(exprs), colnames(cnv))
share_ids_2 <- intersect(colnames(exprs), colnames(cnv_2))

all(share_ids == share_ids_2)
length(share_ids)

ggct_exp <- as.numeric(exprs["GGCT", share_ids])
ggct_cnv <- as.numeric(cnv["GGCT", share_ids])
ggct_cnv_2 <- as.numeric(cnv_2["GGCT", share_ids])
cor(as.numeric(ggct_cnv), as.numeric(ggct_exp))

cor(as.numeric(ggct_cnv_2), as.numeric(ggct_exp))
head(ggct_cnv_2)
head(ggct_exp)
ggct_exp_cnv <- data.frame(ggct_cnv=ggct_cnv_2, ggct_exp=ggct_exp)
head(ggct_exp_cnv)

write.table(ggct_exp_cnv, file="./ggct_expression_and_cnv.txt", sep="\t", col.names=T, quote=F,row.names=F)

length(cli_non_metastasis_samples)
length(cli_metastasis_samples)

all(cli_non_metastasis_samples %in% colnames(cnv_2))

all(cli_metastasis_samples %in% colnames(cnv_2))

#cli_non_metastasis_samples_2 <- cli_non_metastasis_samples[cli_non_metastasis_samples %in% colnames(exprs)]

cli_non_metastasis_samples_2 <- cli_non_metastasis_samples[cli_non_metastasis_samples %in% intersect(colnames(exprs), colnames(cnv_2))]

length(cli_non_metastasis_samples_2)


stasis_samples <- c(cli_non_metastasis_samples_2, cli_metastasis_samples)
ggct_exp2 <- as.numeric(exprs["GGCT", stasis_samples])
ggct_cnv2 <- as.numeric(cnv_2["GGCT", stasis_samples])

ggct_on_metastasis <- data.frame(samples=stasis_samples, expression=ggct_exp2, cnv=ggct_cnv2)

head(ggct_on_metastasis)

ggct_on_metastasis$group_by_exp <- NA
ggct_on_metastasis$group_by_cnv <- NA
ggct_on_metastasis$pathologic_M <- NA


# pathologic_M M0/M1
ggct_on_metastasis[1:length(cli_non_metastasis_samples_2), "pathologic_M"] <- "M0"
ggct_on_metastasis[(length(cli_non_metastasis_samples_2)+1):nrow(ggct_on_metastasis), "pathologic_M"] <- "M1"

tail(ggct_on_metastasis)

# split high/low expression/cnv groups
exp_sm <- summary(ggct_on_metastasis$expression)
cnv_sm <- summary(ggct_on_metastasis$cnv)
ggct_on_metastasis$group_by_exp[ggct_on_metastasis$expression<=exp_sm[2]] <- "Low"
ggct_on_metastasis$group_by_exp[ggct_on_metastasis$expression>=exp_sm[5]] <- "High"

ggct_on_metastasis$group_by_cnv[ggct_on_metastasis$cnv <= cnv_sm[2]] <- "Low"
ggct_on_metastasis$group_by_cnv[ggct_on_metastasis$cnv >= cnv_sm[5]] <- "High"


head(ggct_on_metastasis)



cli <- as.data.frame(cli)
rownames(cli) =  cli$sampleID

ggct_on_metastasis$samples <- as.character(ggct_on_metastasis$samples)
#cli[ggct_on_metastasis$samples, "_OS"]
ggct_on_metastasis$OS <- cli[ggct_on_metastasis$samples, "_OS"]
ggct_on_metastasis$OS_IND <- cli[ggct_on_metastasis$samples, "_OS_IND"]

str(ggct_on_metastasis)

write.table(ggct_on_metastasis, file="./ggct_ge_cnv_values_on_metastasis.txt", row.names=F, col.names=T, quote=F, sep="\t")
#summary(ggct_on_metastasis$expression)


# survial analysis on different TNM stages
# >>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<

dim(cli)

table(cli$pathologic_stage)
# extract sample IDs from clinical data
cli_early_stage_samples <- cli[which(cli$pathologic_stage%in%c("Stage I", "Stage IA", "Stage IB")),]$sampleID
length(cli_early_stage_samples)
cli_late_stage_samples     <- cli[which(cli$pathologic_stage%in%c("Stage II", "Stage IIA", "Stage IIB", "Stage IIIA", "Stage IIIB", "Stage IV")), ]$sampleID

length(cli_late_stage_samples)


# we only use sample which have cnv and gene expression data both.
cli_early_stage_samples_2 <- cli_early_stage_samples[cli_early_stage_samples %in% intersect(colnames(exprs), colnames(cnv_2))]
length(cli_early_stage_samples_2)

cli_late_stage_samples_2 <- cli_late_stage_samples[cli_late_stage_samples %in% intersect(colnames(exprs), colnames(cnv_2))]
length(cli_late_stage_samples_2)


stage_samples <- c(cli_early_stage_samples_2, cli_late_stage_samples_2)
length(stage_samples)

ggct_exp2 <- as.numeric(exprs["GGCT", stage_samples])
ggct_cnv2 <- as.numeric(cnv_2["GGCT", stage_samples])
length(ggct_cnv2)
length(ggct_exp2)


ggct_on_stages <- data.frame(samples=stage_samples, expression=ggct_exp2, cnv=ggct_cnv2)

head(ggct_on_stages)

ggct_on_stages$group_by_exp <- NA
ggct_on_stages$group_by_cnv <- NA
ggct_on_stages$classes <- NA


# pathologic_M M0/M1
ggct_on_stages[1:length(cli_early_stage_samples_2), "classes"] <- "early_stage"
ggct_on_stages[(length(cli_early_stage_samples_2)+1):nrow(ggct_on_stages), "classes"] <- "late_stage"

table(ggct_on_stages$classes)

# split high/low expression/cnv groups
exp_sm2 <- summary(ggct_on_stages$expression)
cnv_sm2 <- summary(ggct_on_stages$cnv)
ggct_on_stages$group_by_exp[ggct_on_stages$expression<=exp_sm2[2]] <- "Low"
ggct_on_stages$group_by_exp[ggct_on_stages$expression>=exp_sm2[5]] <- "High"

ggct_on_stages$group_by_cnv[ggct_on_stages$cnv <= cnv_sm2[2]] <- "Low"
ggct_on_stages$group_by_cnv[ggct_on_stages$cnv >= cnv_sm2[5]] <- "High"

table(ggct_on_stages$group_by_exp)



#ggct_on_metastasis$samples <- as.character(ggct_on_metastasis$samples)
#cli[ggct_on_metastasis$samples, "_OS"]
ggct_on_stages$OS <- cli[ggct_on_stages$samples, "_OS"]
ggct_on_stages$OS_IND <- cli[ggct_on_stages$samples, "_OS_IND"]

str(ggct_on_stages)

write.table(ggct_on_stages, file="./ggct_ge_cnv_values_on_TNM_stages.txt", row.names=F, col.names=T, quote=F, sep="\t")


save(ggct_on_metastasis, ggct_on_stages, exp_sm, exp_sm2, cnv_sm, cnv_sm2, file="./ggct_survival_data_on_metastasis_and_stages.Rdata")
