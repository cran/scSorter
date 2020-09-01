## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(scSorter)

## ------------------------------------------------------------------------
load(url('https://github.com/hyguo2/scSorter/blob/master/inst/extdata/TMpancreas.RData?raw=true'))

## ------------------------------------------------------------------------
expr[1:5, 1:5]
dim(expr)

## ------------------------------------------------------------------------
head(anno)

## ------------------------------------------------------------------------
library(Seurat)
expr_obj = CreateSeuratObject(expr)
expr_obj <- NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)

## ------------------------------------------------------------------------
expr_obj <- FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
topgenes <- head(VariableFeatures(expr_obj), 2000)

expr = GetAssayData(expr_obj)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

## ------------------------------------------------------------------------
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

## ------------------------------------------------------------------------
rts <- scSorter(expr, anno)

## ------------------------------------------------------------------------
print(table(rts$Pred_Type))

## ------------------------------------------------------------------------
mis_rate = 1 - mean(rts$Pred_Type == true_type)
round(mis_rate, 4)

## ------------------------------------------------------------------------
table(true_type, rts$Pred_Type)

