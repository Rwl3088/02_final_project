rm(list=ls())

library(liger)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

data <- read.table(
    sep=',',
    header=TRUE,
    check.names=FALSE,
    '/B_ALL/results/v00/B/B.qc.csv',
    row.names=1,
)

meta.data <- read.table(
    sep=',',
    header=TRUE,
    check.names=FALSE,
    '/B_ALL/results/v00/B/B.metadata.csv',
    row.names=1,
)

sdata <- CreateSeuratObject(
    t(data),
    project='FULL',
    assay='RNA',
    meta.data=meta.data,
)

sdata <- NormalizeData(sdata)
sdata <- FindVariableFeatures(sdata)
sdata <- ScaleData(sdata, split.by='batch', do.center=FALSE)
sdata <- RunOptimizeALS(sdata, k=20, lambda=5, split.by='batch')
sdata <- RunQuantileAlignSNF(sdata, split.by='batch')
sdata <- RunUMAP(sdata, dims = 1:ncol(sdata[["iNMF"]]), reduction = "iNMF")

pdf("liger_umap.pdf")
DimPlot(sdata, group.by = c("batch"), ncol = 1)
dev.off()

saveRDS(sdata, file='/B_ALL/results/v00/B/B.liger.rds')
sdata <- readRDS('/B_ALL/results/v00/B/B.liger.rds')

write.table(
    as.data.frame(sdata@reductions$iNMF@cell.embeddings),
    file='/B_ALL/results/v00/B/B.liger.csv',
    sep=',',
    quote=FALSE,
)
