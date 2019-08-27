library(dplyr)
library(Seurat)
library(purrr)
library(cowplot)
library(parallel)
library(roxygen2)
library(reshape2)
library(tibble)



Output.Meta = function(seurat.ob,filt,annotations,multi = NA) {
  if (!is.na(multi)) {
    Idents(seurat.ob) = seurat.ob$stim
    seurat.obs = lapply(1:length(multi), function (x) SubsetData(seurat.ob,ident.use = multi[x]))
    ret = lapply(1:length(multi), function (x) Output.Meta(seurat.obs[[x]],filt,annotations))
    names(ret) = multi
    return(ret)
  }
  
  cells = names(seurat.ob$seurat_clusters)
  cells.filtered = names(seurat.ob$seurat_clusters[as.numeric(seurat.ob$seurat_clusters) %in% filt])
  cluster.filtered = as.numeric(seurat.ob$seurat_clusters[as.numeric(seurat.ob$seurat_clusters) %in% filt])
  
  cluster.filtered.ind = cluster.filtered + 1
  cluster.anno = annotations[cluster.filtered.ind]
  
  names(cluster.anno) = cells.filtered
  anno.out = rownames_to_column(as.data.frame(cluster.anno,stringsAsFactors = F)) 
  colnames(anno.out) = c('Cell',"cell_type")
  return(anno.out)
}



Output.Count = function(seurat.ob,filt,multi = NA) {
  if (!is.na(multi)) {
    Idents(seurat.ob) = seurat.ob$stim
    seurat.obs = lapply(1:length(multi), function (x) SubsetData(seurat.ob,ident.use = multi[x]))
    ret = lapply(1:length(multi), function (x) Output.Count(seurat.obs[[x]],filt))
    names(ret) = multi
    return(ret)
  }
  
  cells = names(seurat.ob$seurat_clusters)
  cells.filtered = names(seurat.ob$seurat_clusters[as.numeric(seurat.ob$seurat_clusters) %in% filt])
  
  normalized.count = as.matrix(seurat.ob@assays$RNA@data)
  normalized.count.filtered = normalized.count[,cells.filtered]
  
  normalized.count.filtered = as.data.frame(normalized.count.filtered)
  normalized.count.filtered = rownames_to_column(normalized.count.filtered,var = "Gene")
  normalized.count.filtered$Gene = toupper(normalized.count.filtered$Gene)
  return(normalized.count.filtered)
}

Filter.NA.ROW = function(mat,sel = 2:5) {
  row.na = apply(mat[,sel],MARGIN = 1, function(x) all(is.na(x)))
  mat.filtered = mat[row.na == F,]
  return(mat.filtered)
}



