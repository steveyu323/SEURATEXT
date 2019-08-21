library(dplyr)
library(Seurat)
library(purrr)
library(cowplot)
library(parallel)
library(roxygen2)
library(reshape2)
library(tibble)

#'  Build.Seurat.Cluster
#' @description A wrpper function for Seurat cluster building process
#' @param seurat.ob The seurat object that does not have UMAP run and does not
#' hold identified clusters
#' @return A seurat object with UMAP reduction and clusters annotated
Build.Seurat.Cluster = function(seurat.ob) {
  seurat.ob <- ScaleData(seurat.ob, verbose = T)
  seurat.ob <- RunPCA(seurat.ob, npcs = 30, verbose = T)
  # t-SNE and Clustering
  seurat.ob <- RunUMAP(seurat.ob, reduction = "pca", dims = 1:30)
  seurat.ob <-
    FindNeighbors(seurat.ob, reduction = "pca", dims = 1:30)
  seurat.ob <- FindClusters(seurat.ob, resolution = 0.6)
}



#'  Build.ConserveMarkers.All
#' @description A wrapper function to simultaneously find conserve markers for
#' all the clusters
#' @param seurat.ob The seurat object with clustering information
#' @return A list object with each cluster's conserved element as a marker
Build.ConserveMarkers.All = function(seurat.ob) {
  Idents(seurat.ob) = seurat.ob$seurat_clusters
  max.clust.num = length(levels(seurat.ob$seurat_clusters)) - 1
  markers = lapply(0:max.clust.num, function(x)
    FindConservedMarkers(
      seurat.ob,
      ident.1 = x,
      group = "stim",
      verbose = T
    ))
  return(markers)
}

# FindConservedMarkers.flat = function(seurat.ob, ident.1, group , verbose ) {
#   markers = FindConservedMarkers(seurat.ob, ident.1, grouping.var = group, verbose )
#   markers$cluster = ident.1
#   return(markers)
# }

#'  FeaturePlot.All
#' @description A wrapper function to plot the top-three condident markers
#' from each cluster
#' @param seurat.ob The seurat object with clustering information
#' @param markers The marker list with each element as a dataframe from
#' FindMarkers output
#' @return feature plots
FeaturePlot.All = function(seurat.ob, markers) {
  max.clust.num = length(levels(seurat.ob$seurat_clusters)) - 1
  p = lapply(0:max.clust.num, function(x)
    FeaturePlot(
      seurat.ob,
      features = row.names(markers[[x + 1]])[1:3],
      split.by = "stim",
      max.cutoff = 3,
      cols = c("grey", "red")
    ))
  return(p)
}



#'  Find.Markers.Each
#' @description A wrapper function to output: For each condition in the "stim"
#' field, each cluster's marker that are significantly overexpressed in a
#' given cluster compared to the rest of the clusters
#' @param seurat.ob The seurat object with clustering information
#' @param min.pct as in Seurat::FindAllMarkers function
#' @param logfc.threshold the minimum log fold change threshold to
#' detect differentially expressed change
#' @param multi if "stim" field is not "CTRL" VS "STIM", then the multi field
#' would take in a vector with the name of the stim conditions
#' @return a object for each stim condition, stored in a ret object list
Find.Markers.Each = function(seurat.ob,
                             min.pct = 0.1,
                             logfc.threshold = 0.3,
                             multi = NA) {
  if (!is.na(multi)) {
    Idents(seurat.ob) = seurat.ob$stim
    seurat.obs = lapply(1:length(multi), function (x)
      SubsetData(seurat.ob, ident.use = multi[x]))
    for (x in 1:length(multi)) {
      Idents(seurat.obs[[x]]) = seurat.obs[[x]]$seurat_clusters
    }
    ret = lapply(1:length(multi), function (x)
      FindAllMarkers(
        seurat.obs[[x]],
        only.pos = TRUE,
        min.pct = min.pct ,
        logfc.threshold = logfc.threshold
      ))
    names(ret) = multi
  } else {
    seurat.ob.wt = subset(seurat.ob, subset = stim == 'CTRL')
    seurat.ob.ko = subset(seurat.ob, subset = stim == 'STIM')
    ctrl.markers <-
      FindAllMarkers(
        seurat.ob.wt,
        only.pos = TRUE,
        min.pct = min.pct ,
        logfc.threshold = logfc.threshold
      )
    stim.markers <-
      FindAllMarkers(
        seurat.ob.ko,
        only.pos = TRUE,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
      )
    ret = c()
    ret$ctrl.markers = ctrl.markers
    ret$stim.markers = stim.markers
  }

  return(ret)
}



#'  DE.Each.Cluster
#' @description A wrapper function to output: For each anchored cluster
#' output the differentially expressed gene between a pair
#' @param seurat.ob The seurat object with clustering information
#' @param pair a pair of conditions to be differentially analyzed in the
#' corresponding "stim" field
#' @return flattened list of differentially expressed genes table
DE.Each.Cluster = function(seurat.ob, pair = NA) {
  if (!is.na(pair)) {
    Idents(seurat.ob) = seurat.ob$stim
    seurat.ob = SubsetData(seurat.ob, ident.use = pair)
  }

  Idents(seurat.ob) = seurat.ob$seurat_clusters
  seurat.ob$celltype.stim <-
    paste(Idents(seurat.ob), seurat.ob$stim, sep = "_")
  seurat.ob$celltype <- Idents(seurat.ob)
  Idents(seurat.ob) <- "celltype.stim"

  max.clust.num = length(levels(seurat.ob$seurat_clusters)) - 1
  diff.out = lapply(0:max.clust.num, function(x)
    DiffOutput(seurat.ob, x, pair))
  return(rbind_list(diff.out))
}


#'  DE.Each.Cluster
#' @description Helper function for DE.Each.Cluster
DiffOutput = function(metadata, clust.name, pair) {
  if (is.na(pair)) {
    pair = c("CTRL", "STIM")
  }
  id.1 = paste0(clust.name, '_', pair[1])
  id.2 = paste0(clust.name, '_', pair[2])
  proceed.1 = tryCatch(
    WhichCells(metadata, idents = id.1),
    error = function(e)
    {
      NA
    }
  )
  proceed.2 = tryCatch(
    WhichCells(metadata, idents = id.2),
    error = function(e)
    {
      NA
    }
  )
  if ((!is.na(proceed.1)) & (!is.na(proceed.2))) {
    clust.diff <-
      FindMarkers(
        metadata,
        ident.1 = id.1,
        ident.2 = id.2,
        verbose = T
      )
    clust.diff = rownames_to_column(clust.diff, "gene")
    clust.diff$cluster = clust.name
    clust.diff = clust.diff %>% filter(pct.1 > 0.2  & pct.2 > 0.2)
    return(clust.diff)
  }
}

#'  DE.Each.Cluster
#' @description A wrapper function to output: For each anchored cluster
#' output the differentially expressed gene between a pair
#' @param seurat.ob The seurat object with clustering information
#' @param pair a pair of conditions to be differentially analyzed in the
#' corresponding "stim" field
#' @param multiple a logical boolean on with there are more than 2 conditions,
#' default to FALSE. If FALSE, the condition must be indicated as "CTRL" and
#' "STIM"
#' @return A seurat out object that can be used for input into app.R for
#' ShinyApp visualization
Shine.Out = function(ob, diff, markers.each, markers.conserved,multiple = F) {
  if (multiple) {
    out = c()
    out$markers.each = markers.each
    out$conserved.markers = markers.conserved
    out$diff = diff
    out$ob = ob
    return(out)
  } else {
    out = c()
    out$stim.markers = markers.each$stim.markers
    out$ctrl.markers = markers.each$ctrl.markers
    out$conserved.markers = markers.conserved
    out$diff = diff
    out$ob = ob
    return(out)
  }
}


