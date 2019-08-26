#' Filter.Gene
#' @description keep only the cells barcode that have expression of gene higher
#' than the threshold for now use seurat.ob->assays$RNA@data field,
#' but may change the assay type
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param gene the gene gonna impose the threshold filter on
#' @param threshold the count number (exclusively greater) to define the
#' presence of a certain gene in the cell. Default to 0
#' @param reverse whether we would select the cells absent with the gene
#' count smaller or equal to the threshold
#' @return  return the remained cell's barcode
Filter.Gene = function(seurat.ob,gene,threshold = 0,reverse = F) {
  # filter by row 'Cd8a'
  barcodes = colnames(seurat.ob@assays$integrated@data)
  if(reverse){
    barcodes.sel = seurat.ob@assays$RNA@data[gene,] <= threshold
  } else {
    barcodes.sel = seurat.ob@assays$RNA@data[gene,] > threshold
  }
  barcodes.filtered = barcodes[barcodes.sel]
  return(barcodes.filtered)
}


#' BarcodeToObject
#' @description Given the selected barcodes, the function will create a new a
#' new seurat object that contains only the cell of interest. Other fields will
#' remain unchanged
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param barcodes.filtered The barcode of the cells of interests to be
#' extracted out from the object
#' @return  The filtered seurat object
BarcodeToObject = function(seurat.ob,barcodes.filtered) {
  filt = names(seurat.ob$orig.ident) %in% barcodes.filtered
  seurat.ob$gene.filt = F

  filt = names(seurat.ob$orig.ident) %in% barcodes.filtered
  seurat.ob$gene.filt[filt] = T

  filt.cells = subset(seurat.ob, subset = gene.filt == T)
  filt.cells$gene.filt = NULL
  return(filt.cells)
}



#' Plot.FeatureScatter
#' @description Wrapper function for the Seurat.FeatureScatter for 2 different
#' conditions
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param x The gene name of the first cell marker
#' @param y The gene name of the second cell marker
#' @return  a grid-arranged FeatureScatter Seurat plot
Plot.FeatureScatter = function(seurat.ob,x,y){
  # seurat.ob = macrophage.combined
  # x = 'Irf5'
  # y = 'Irf4'
  DefaultAssay(seurat.ob) <- "RNA"
  ctrl = subset(seurat.ob,subset = stim =='CTRL')
  stim = subset(seurat.ob,subset = stim =='STIM')
  p1 = FeatureScatter(ctrl, feature1 = x, feature2 = y,slot = 'data')
  p2 = FeatureScatter(stim, feature1 = x, feature2 = y,slot = 'data')
  p = plot_grid(p1,p2,labels = c('CTRL','STIM'))
  return(p)
}


#' Build.Cyto
#' @description The four quadrant: x+/y+ = 1 x+/y- = 2  x-/y- = 3   x-/y+ = 4
#' The function will add a $cyto field into the current seurat object to
#' specify the quadant a cell belongs to
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param x The gene name of the first cell marker
#' @param y The gene name of the second cell marker
#' @param x.thresh threshold to deem positive for the first cell marker
#' @param y.thresh threshold to deem positive for the second cell marker
#' @return  filtered seurat object with $cyto field added used for further
Build.Cyto = function(seurat.ob,x,y,x.thresh = 0,y.thresh = 0) {
  x.pos = Filter.Gene(seurat.ob  ,x,x.thresh)
  x.neg = Filter.Gene(seurat.ob  ,x,x.thresh,reverse = T)
  y.pos = Filter.Gene(seurat.ob ,y,y.thresh)
  y.neg = Filter.Gene(seurat.ob ,y,y.thresh,reverse = T)
  seurat.ob$cyto = 0
  names(seurat.ob$cyto) = names(seurat.ob$orig.ident)
  seurat.ob$cyto[intersect(x.pos,y.pos)] = 1
  seurat.ob$cyto[intersect(x.pos,y.neg)] = 2
  seurat.ob$cyto[intersect(x.neg,y.neg)] = 3
  seurat.ob$cyto[intersect(x.neg,y.pos)] = 4
  return(seurat.ob)
}

#' Plot.Cyto.Count
#' @description plot the quadrant grapsh with
#' x+/y+ = 1 x+/y- = 2  x-/y- = 3   x-/y+ = 4
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param x The gene name of the first cell marker
#' @param y The gene name of the second cell marker
#' @param split whether there are two conditions as CTRL and STIM
#' @param return.stats whether to return the stats default to FALSE
#' @return Plot that Count the number of cells in each quadrant
#' of the seurat.cyto object
Plot.Cyto.Count = function(seurat.ob,x,y,split = F,return.stats = F) {
  if (split) {
    ctrl = subset(seurat.ob,subset = stim =='CTRL')
    stim = subset(seurat.ob,subset = stim =='STIM')
    p.ctrl = Plot.Cyto.Count(ctrl,x,y,split = F)
    p.stim = Plot.Cyto.Count(stim,x,y,split = F)
    p = plot_grid(p.ctrl,p.stim,labels = c('CTRL','STIM'))
  } else {
    num.1 =  sum(as.numeric(seurat.ob$cyto == 1))
    num.2 =  sum(as.numeric(seurat.ob$cyto == 2))
    num.3 =  sum(as.numeric(seurat.ob$cyto == 3))
    num.4 =  sum(as.numeric(seurat.ob$cyto == 4))
    total = num.1 + num.2 + num.3 + num.4
    x.label = c('pos','pos','neg','neg')
    y.label = c('pos','neg','neg','pos')
    count = c(num.1,num.2,num.3,num.4)
    dat = data.frame(x.label,y.label,count)
    dat.melted = melt(data = dat)
    dat.melted$percent = dat.melted$value/total *100
    p = ggplot(data = dat.melted,aes(x = x.label,y = y.label,fill = -value)) + geom_tile() +
      labs(title = "", x = x, y = y)  +
      geom_text(aes(label = value),colour = 'white') +
      geom_text(aes(label = paste(round(percent,digits = 3),'%',sep = '')),colour = 'white',vjust = 2)
  }
  return (p)
}

#' Plot.Cyto.Cluster
#' @description Map back the cytometry assignment to the DimPlot of the
#' Seurat Objects clutsers
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param x The gene name of the first cell marker
#' @param y The gene name of the second cell marker
#' @param split whether there are two conditions as CTRL and STIM
#' @return Plot that map red dots on the original DimPlot to indicate the
#' presence of the two marker genes with four quadrants
Plot.Cyto.Cluster = function(seurat.ob,x,y,split = F) {
  if(split){
    ctrl = subset(seurat.ob,subset = stim =='CTRL')
    stim = subset(seurat.ob,subset = stim =='STIM')
    p.ctrl = Plot.Cyto.Cluster(ctrl,x,y,split = F)
    p.stim = Plot.Cyto.Cluster(stim,x,y,split = F)
    p = c()
    p$ctrl = p.ctrl
    p$stim = p.stim
  } else {
    # x= "Cd3e"
    # y = "Adgre1"
    # seurat.ob = immune.combined.cyto
    xlabs = paste(x,c('-','+','-','+'))
    ylabs = paste(y,c('+','+','-','-'))
    labs = paste(xlabs,ylabs,sep = '/')
    seurat.ob[[labs[1]]] = seurat.ob$cyto == 4
    seurat.ob[[labs[2]]] = seurat.ob$cyto == 1
    seurat.ob[[labs[3]]] = seurat.ob$cyto == 3
    seurat.ob[[labs[4]]] = seurat.ob$cyto == 2
    p4 = FeaturePlot(seurat.ob, features = labs[1], cols = c("grey", "red")) + NoLegend() + NoAxes()
    p1 = FeaturePlot(seurat.ob, features = labs[2], cols = c("grey", "red")) + NoLegend() + NoAxes()
    p3 = FeaturePlot(seurat.ob, features = labs[3], cols = c("grey", "red")) + NoLegend() + NoAxes()
    p2 = FeaturePlot(seurat.ob, features = labs[4], cols = c("grey", "red")) + NoLegend() + NoAxes()
    p = plot_grid(p4,p1,p3,p2,ncol = 2)
  }
  return(p)
}



#' Recluster.Quadrant
#' @description with n.quadrant label:
#' x+/y+ = 1 x+/y- = 2  x-/y- = 3   x-/y+ = 4
#' return the subset of seurat ob cells that lies in a given quadrant
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param n.quadrant the numbering of the selected quadrant that would like to
#' be reclustered
#' @return a new seurat object only within a specified quadrant
Recluster.Quadrant = function(seurat.ob,n.quadrant = 1) {
  # seurat.ob = immune.combined.cyto.Adgre1
  # n.quadrant = 2
  quad.barcode = names(seurat.ob$cyto)[(seurat.ob$cyto == n.quadrant)]
  seurat.ob.filtered = BarcodeToObject(seurat.ob,quad.barcode)
  return(seurat.ob.filtered)
}


#' Plot.Cluster.Percentage
#' @param seurat.ob a seurat object type that have a count a count table inside
#' with proper column name as cells and rowname as the gene
#' @param grouped a boolean to specify if the seurat object is grouped as
#' 'CTRL' and 'STIM', default to false
#' @param multi a boolean to specify if the seurat object is grouped as
#' a specified vector of group names, default to NA
#' @return return barplot with the number of cells of each cluster
#' within the sample/population
Plot.Cluster.Percentage = function(seurat.ob,grouped = F,multi = NA) {
  if (grouped) {
    seurat.ob.wt = subset(seurat.ob,subset = stim =='CTRL')
    seurat.ob.ko = subset(seurat.ob,subset = stim =='STIM')
    p1 = Plot.Cluster.Percentage(seurat.ob.wt)
    p2 = Plot.Cluster.Percentage(seurat.ob.ko)
    p = plot_grid(p1,p2,labels = c('CTRL','STIM'))
  } else if (!is.na(multi)){
    Idents(seurat.ob) = seurat.ob$stim
    seurat.obs = lapply(1:length(multi), function (x) SubsetData(seurat.ob,ident.use = multi[x]))
    p = lapply(1:length(multi), function (x) Plot.Cluster.Percentage(seurat.obs[[x]]))
    p = plot_grid(plotlist = p,labels = multi)
  } else {
    dat.sum = rownames_to_column(as.data.frame(summary(seurat.ob$seurat_clusters)))
    colnames(dat.sum) = c('name','value')
    dat.sum$value = dat.sum$value/sum(dat.sum$value)
    p = ggplot(dat = dat.sum,aes(x = name,y = value)) + geom_bar(stat ='identity') + ylim(0, 1)
  }
  return(p)
}
