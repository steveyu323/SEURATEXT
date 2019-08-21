# Library dependency
library(dplyr)
library(Seurat)
library(purrr)
library(cowplot)
library(parallel)
library(roxygen2)
library(reshape2)
library(tibble)


#'  Build.Ref.Markers
#' @description takes in a .csv table and build up the ref.marker object for assignment of cell types. The list should be in the following format gene | cell.type1,cell.type2,....
#' @example atlas_markers = Build.Ref.Markers(path =  "../../data/jay_cell_markers.txt",sep = '\t',splitter = c("[[:punct:]]","[[:space:]]"))
#' @param path the path to the csv file with reference marker
#' @param sep the delimiter of the scv file
#' @param splitter the seperator between markers
#' @return a ref_markers object
Build.Ref.Markers = function(path,sep = ',',splitter = ',' ) {
  ref_markers = read.csv(file =  path,sep = sep,stringsAsFactors = F)

  # re-structure the file
  ref_markers.reframe = data.frame("hi","bye")
  names(ref_markers.reframe)<-c("gene","celltype")

  for (i in 1:nrow(ref_markers)) {
    curr.genes = unlist(strsplit(ref_markers[i,2],split = splitter))
    for (j in 1:length(curr.genes)) {
      temp = data.frame(curr.genes[j],ref_markers[i,1] )
      names(temp)<-c("gene","celltype")
      ref_markers.reframe = rbind(ref_markers.reframe,temp)
    }
  }
  ref_markers.reframe$celltype = as.character(ref_markers.reframe$celltype)
  ref_markers.reframe = ref_markers.reframe[2:nrow(ref_markers.reframe),]
  ref_markers = ref_markers.reframe

  # change the label all to lower case
  ref_markers$gene = tolower(ref_markers$gene)
  # sanity filter to get rid of NA row due to format issue
  ref_markers = ref_markers %>% filter(gene != 'NA' & celltype != 'NA') %>% distinct()
  return(ref_markers)
}


#'  Build.Exp.Markers
#' @description From the Seurat FindAllMarkers Output to a Input for Assign.Cell.Type
#' @param multiple if true, will take the list object from the output of lapply and flatten the returned markers with responsive cluster as the new column
#' @param markers the marker table output from FindAllMarkers function
#' @example
#' markers <- FindAllMarkers(object = data.0531.copy, only.pos = T, min.pct = 0.25, logfc.threshold = 1)
#' markers = Build.Exp.Markers(markers)
#' @return a exp_marker object
Build.Exp.Markers = function(markers,multiple = F) {
  if (multiple) {
    markers = lapply(1:length(markers), function(x) cbind(markers[[x]], sample = x))
    markers.reframe = rbind_list(markers)
    markers.reframe$cluster = as.integer(markers.reframe$cluster)
    markers = markers.reframe
  } else {
    markers = cbind(markers, sample = 1)
    markers.reframe = markers
    markers.reframe$cluster = as.integer(as.character(markers.reframe$cluster))
    markers = markers.reframe
  }
  markers$gene = tolower(markers$gene)
  return(markers)
}



#'  Assign.Cell.Type
#' @param ref.markers Markers | Cell Type dataframe from outside reference
#' @param exp.markers each sample's cluster has a set of overlapped markers with ref.markers
#' @description The function takes top 10 (if more than) of the exp.markers to decide the top three cell type a cluster
#' could belong to. use the logFC fold as the weight for the matrix, and output the cell type and corresponding marker
#' @return A data frame with sample|cluster|type1|marker1|type2|marker2|type3|marker3
Assign.Cell.Type = function(ref.markers,exp.markers) {
  # for debug input
  # exp.markers = markers
  # ref.markers = atlas_markers
  # data.0531  = readRDS(file = '../temp/0531_seurat.rds')
  # levels(Idents(data.0531))
  # data.0531.copy = data.0531
  # new.cluster.ids <- as.character( seq(1,12,1))
  # names(new.cluster.ids) <- levels(data.0531.copy)
  # data.0531.copy <- RenameIdents(data.0531.copy, new.cluster.ids)


  # filter to keep only the genes overlapping with the atlas_markers set
  sample.clust = as.matrix(exp.markers %>% group_by(sample,cluster) %>% summarise())
  exp.markers = exp.markers %>% filter(gene %in% ref.markers$gene)
  cell.types = unique(ref.markers$celltype)

  result <- data.frame(sample=integer(),
                       cluster=integer(),
                       type1.name=character(),
                       type1.gene=character(),
                       type1.fc=double(),
                       type2.name=character(),
                       type2.gene=character(),
                       type2.fc=double(),
                       type3.name=character(),
                       type3.gene=character(),
                       type3.fc=double(),
                       stringsAsFactors=FALSE)


  for (i in 1:nrow(sample.clust)) {
    # create a list for all different cell types
    score = rep(0,length(cell.types))
    gene.cache = rep('',length(cell.types))
    names(score) = cell.types
    names(gene.cache) = cell.types

    # get the current cluster and sample id
    curr.sample = sample.clust[i,1]
    curr.clust = sample.clust[i,2]
    # get the current marker genes and the corresponding logFC
    curr.exp.markers = as.data.frame(exp.markers %>% filter(sample == curr.sample & cluster == curr.clust) %>% select(c("gene","avg_logFC")))
    for (j in 1:nrow(curr.exp.markers)) {
      print(curr.exp.markers$gene[j])
      temp.cell.type = ref.markers[ref.markers$gene == curr.exp.markers$gene[j],2]
      for (type in temp.cell.type) {
        score[type] = score[type] + curr.exp.markers$avg_logFC[j]
        gene.cache[type] = paste(gene.cache[type],paste(curr.exp.markers$gene[j],round(curr.exp.markers$avg_logFC[j],digits = 2),sep = '-'),sep = ' ')
      }
    }

    # output a report dataframe
    fc =  score[order(-score)][1:3]
    gene = gene.cache[order(-score)][1:3]
    name = names(fc)
    name[fc == 0] = 'unknown'
    # bind to the final result
    temp.result = data.frame(sample=curr.sample,
                             cluster=curr.clust,
                             type1.name=name[1],
                             type1.gene=gene[1],
                             type1.fc=fc[1],
                             type2.name= name[2],
                             type2.gene= gene[2],
                             type2.fc=fc[2],
                             type3.name=name[3],
                             type3.gene=gene[3],
                             type3.fc=fc[3],
                             stringsAsFactors=FALSE)
    result = rbind(result,temp.result)

  }

  return(result)
}


#'  Filter.Gene.Mat
#' @param count.matrix a matrix type that have a count for value and proper column name as cells and rowname as the gene
#' @param gene the gene gonna impose the threshold filter on
#' @param threshold the threshold to filter from count.matrix
#' @description keep only the cells that have expression of the gene higher than the threshold, and output the new matrix object
#' @return the filtered matrix
Filter.Gene.Mat = function(count.matrix,gene,threshold) {
  # filter by row 'Cd8a'
  pos.names = colnames(count.matrix)[count.matrix[gene,] >= threshold]
  pos.matrix = count.matrix[,pos.names]
  return(pos.matrix)
}




# CreateSeurat = function(count.matrix,project.name = 'project'){
#   ## Build a new Seurat object
#   seurat.ob = CreateSeuratObject(counts = count.matrix, project = project.name, min.cells = 3, min.features = 200)
#
#   ### Normalization
#   seurat.ob[["percent.mt"]] <- PercentageFeatureSet(seurat.ob, pattern = "^MT-|^mt")
#   head(seurat.ob@meta.data, 5)
#   # Visualize QC metrics as a violin plot
#   VlnPlot(seurat.ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#   seurat.ob <- subset(seurat.ob, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 5)
#   seurat.ob <- NormalizeData(seurat.ob, normalization.method = "LogNormalize", scale.factor = 10000)
#
#   ### Find Variable Features
#   seurat.ob <- FindVariableFeatures(seurat.ob, selection.method = "vst", nfeatures = 2000)
#   # Identify the 10 most highly variable genes
#   top10 <- head(VariableFeatures(seurat.ob), 10)
#
#   # plot variable features with and without labels
#   plot1 <- VariableFeaturePlot(seurat.ob)
#   plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#   CombinePlots(plots = list(plot1, plot2))
#
#   ### ScaleData
#   all.genes <- rownames(seurat.ob)
#   seurat.ob <- ScaleData(seurat.ob, features = all.genes)
#
#   ### Run PCA
#   seurat.ob <- RunPCA(seurat.ob, features = VariableFeatures(object = seurat.ob))
#   VizDimLoadings(seurat.ob, dims = 1:2, reduction = "pca")
#   DimPlot(seurat.ob, reduction = "pca")
#   DimHeatmap(seurat.ob, dims = 1, cells = 500, balanced = TRUE)
#   DimHeatmap(seurat.ob, dims = 1:15, cells = 500, balanced = TRUE)
#
#   seurat.ob <- JackStraw(seurat.ob, num.replicate = 100)
#   seurat.ob <- ScoreJackStraw(seurat.ob, dims = 1:20)
#   JackStrawPlot(seurat.ob, dims = 1:15)
#   ElbowPlot(seurat.ob)
#
#
#   ### Linear Clustering
#   seurat.ob <- FindNeighbors(seurat.ob, dims = 1:20)
#   seurat.ob <- FindClusters(seurat.ob, resolution = 0.6)
#   head(Idents(seurat.ob), 5)
#
#   ### Non-linear UMAP
#   seurat.ob <- RunUMAP(seurat.ob, dims = 1:20)
#   DimPlot(seurat.ob, reduction = "umap")
#   # saveRDS(seurat.ob, file = "../temp/seurat.ob.rds")
#   return(seurat.ob)
# }
#
# Annotate.DimPlot = function(seurat.ob,assignment) {
#   new.cluster.ids <- assignment$type1.name
#   names(new.cluster.ids) <- levels(seurat.ob)
#   seurat.ob <- RenameIdents(seurat.ob, new.cluster.ids)
#   return(seurat.ob)
# }
#
# CreateSeurat.Compare = function(ctrl.matrix,stim.matrix) {
#   ctrl = CreateSeuratObject(counts = ctrl.matrix, project = 'CTRL', min.cells = 3, min.features = 200)
#   ctrl$stim = "CTRL"
#   ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-|^mt")
#   ctrl <- subset(ctrl, subset = nFeature_RNA > 500  & percent.mt < 5)
#   print(ctrl)
#   ctrl <- NormalizeData(ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
#   ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
#
#   stim = CreateSeuratObject(counts = stim.matrix, project = 'STIM', min.cells = 3, min.features = 200)
#   stim$stim = "STIM"
#   stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-|^mt")
#   stim <- subset(stim, subset = nFeature_RNA > 500  & percent.mt < 5)
#   print(stim)
#   stim <- NormalizeData(stim, normalization.method = "LogNormalize", scale.factor = 10000)
#   stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)
#
#   immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim),k.anchor = 3)
#   immune.combined <- IntegrateData(anchorset = immune.anchors,k.weight = 30)
#   #
#   # immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:25)
#   # immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:25)
#   #
#   DefaultAssay(immune.combined) <- "integrated"
#
#   # Run the standard workflow for visualization and clustering
#   immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#   immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
#   # t-SNE and Clustering
#   immune.combined <- RunUMAP(immune.combined, reduction = "pca",dims = 1:18)
#   immune.combined <- FindNeighbors(immune.combined, reduction = "pca",dims = 1:18)
#   immune.combined <- FindClusters(immune.combined, resolution = 0.6)
#   # saveRDS(immune.combined,file = "../temp/injured-compare.rds")
#   # Visualization
#   # immune.combined = readRDS("../temp/injured-compare.rds")
#   p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
#   p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
#   plot_grid(p1, p2)
#   DimPlot(immune.combined, reduction = "umap", split.by = "stim",label = T)
#
#   return(immune.combined)
#   # saveRDS(object = immune.markers,file = "../temp/immune-markers.rds")
# }

