# Pathway analysis for GGCT high-relation Gene sets (Person coeffiency |r| > 0.4)
require(tidyverse)
setwd("F:/Wangshixiang/GGCT")

exprs <- as.data.frame(read_tsv(file="./HiSeqV2"))
cli   <- read_tsv(file="./LUAD_clinicalMatrix")
cli

nchar(colnames(exprs)[2])
nchar(cli$sampleID[1])

table(cli$sample_type)
table(cli$new_neoplasm_event_type)

# extract sample IDs from clinical data
cli_non_metastasis_samples <- cli[is.na(cli$new_neoplasm_event_type),]$sampleID
cli_metastasis_sample      <- cli[which(cli$new_neoplasm_event_type == "Distant Metastasis"), ]$sampleID

# use sample IDs above to filter expression data
rownames(exprs) <- exprs$sample
# check
all(rownames(exprs) == exprs$sample)
# remove sample variable
exprs <- exprs[, -1]
# only keep non-metastasis and metastasis samples
exprs <- exprs[, colnames(exprs)%in%c(cli_metastasis_sample, cli_non_metastasis_samples)]

# next, we select GGCT high-relation Gene sets
# first, get expression matrix respectively
exprs_non_meta <- exprs[, colnames(exprs)%in%cli_non_metastasis_samples]
exprs_meta     <- exprs[, colnames(exprs)%in%cli_metastasis_sample]

ggct_exp_non_meta <- as.numeric(exprs_non_meta["GGCT", ])
ggct_exp_meta     <- as.numeric(exprs_meta["GGCT", ])

exprs_non_meta <- exprs_non_meta[rownames(exprs_non_meta)!="GGCT", ]
exprs_meta     <- exprs_meta[rownames(exprs_meta)!="GGCT", ]

# calculate correlation coeffiency "r"
calc_r <- function(x, gene_expr){
    cor(x, gene_expr)
}

r_non_meta <- apply(exprs_non_meta, 1, calc_r, gene_expr=ggct_exp_non_meta)
# which(is.na(r_non_meta))
# any(is.na(exprs_non_meta["OR8J1",]))
# cor(as.numeric(exprs_non_meta["OR8J1",]), ggct_exp_non_meta)

r_meta <- apply(exprs_meta, 1, calc_r, gene_expr=ggct_exp_meta)

summary(r_non_meta)
summary(r_meta)

geneset_non_meta <- names(r_non_meta[which(abs(r_non_meta) > 0.4)])
geneset_meta <- names(r_meta[which(abs(r_meta) > 0.4)])


# pathway analysis
require(clusterProfiler)

length(geneset_non_meta)
grep('/', x=geneset_non_meta, value=TRUE)

non_meta.id <- bitr(geneset_non_meta, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

meta.id <- bitr(geneset_meta, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

geneList <- rownames(exprs)
genelist.id <- bitr(geneList, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = "org.Hs.eg.db")


doGO <- function(id, ont=c("BP", "CC", "MF")){
    res <- list()
    res$bp <-   enrichGO(gene          = id,
                         keyType       = "SYMBOL",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
    res$cc <-   enrichGO(gene          = id,
                         keyType       = "SYMBOL",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
    res$mf <-   enrichGO(gene          = id,
                         keyType       = "SYMBOL",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
    return(res)
}

go.non_meta <- doGO(geneset_non_meta)
go.meta     <- doGO(geneset_meta)

k1 <- enrichKEGG(gene  = non_meta.id$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05, pAdjustMethod="BH",
                 qvalueCutoff=0.1)
k2 <- enrichKEGG(gene  = meta.id$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05, pAdjustMethod="BH",
                 qvalueCutoff=0.1)

# k1 <- enrichKEGG(gene  = non_meta.id$ENTREZID,
#                  universe = genelist.id$ENTREZID,
#                  organism     = 'hsa',
#                  pvalueCutoff = 0.05, pAdjustMethod="BH",
#                  qvalueCutoff=0.1)


dotplot(k1)
dotplot(k2)

barplot(go.non_meta$bp, showCategory = 10)
barplot(go.non_meta$cc, showCategory = 10)
barplot(go.non_meta$mf, showCategory = 10)

barplot(go.meta$bp, showCategory = 10)
barplot(go.meta$cc, showCategory = 10)
barplot(go.meta$mf, showCategory = 10)

dotplot(go.non_meta$bp, showCategory = 10)
dotplot(go.non_meta$cc, showCategory = 10)
dotplot(go.non_meta$mf, showCategory = 10)

dotplot(go.meta$bp, showCategory = 10)
dotplot(go.meta$cc, showCategory = 10)
dotplot(go.meta$mf, showCategory = 10)


plotGOgraph(go.non_meta$bp)
plotGOgraph(go.meta$bp)

# GO
# go.non_meta <- enrichGO(gene          = non_meta.id$ENTREZID,
#                         keyType       = "SYMBOL",
#                         OrgDb         = org.Hs.eg.db,
#                         ont           = "BP",
#                         pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.01,
#                         qvalueCutoff  = 0.05)
#
# go.meta     <- enrichGO(gene          = meta.id$ENTREZID,
#                         OrgDb         = org.Hs.eg.db,
#                         ont           = "BP",
#                         pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.01,
#                         qvalueCutoff  = 0.05)

# barplot(go.non_meta,showCategory=10)
# dotplot(go.non_meta, showCategory=10)
#
# barplot(go.meta,showCategory=10)
# dotplot(go.meta, showCategory=10)
#
# plotGOgraph(go.non_meta)
enrichMap(go.non_meta$bp, n=9)
enrichMap(go.meta$bp, n=9)
#
# # cnetplot(go.non_meta, showCategory = 1)
#
# k1 <- enrichKEGG(gene  = geneset_non_meta,
#                         keyType = "SYMBOL",
#                        organism     = 'hsa',
#                        pvalueCutoff = 0.1)
# k2 <- enrichKEGG(gene  = meta.id$ENTREZID,
#                  organism     = 'hsa',
#                  pvalueCutoff = 0.1)
# barplot(k1, showCategory = 10)
# barplot(k2, showCategory = 10)
