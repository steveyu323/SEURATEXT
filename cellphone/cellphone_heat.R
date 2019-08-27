library(pheatmap)

heatmaps_plot = function(meta_file, pvalues_file, count_filename, log_filename, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  
  
  #######   Network
  print(meta_file)
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  
  dim(all_intr)
  
  # filter out the collagen related interactions
  all_intr = all_intr %>% filter(!str_detect(all_intr$interacting_pair,"^COL"))
  
  dim(all_intr)
  
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  

  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = sort(unique(meta[,2]))
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  # here we disable this file generation
  # all_count = all_count[-1,]
  # all_count = as.data.frame(all_count)
  # all_count[,"count"] = as.numeric(all_count[,"count"])
  # all_count = all_count %>% arrange(-count)
  
  # Changhua Yu: for selecting the top interacting groups
  # col.sel = paste0(all_count[1:50,1], '|',all_count[1:50,2])
  
  # write.table(all_count, 'count_network.txt', sep='\t', quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
  }
  
  
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= sort(unique(meta[,2]))
    colnames(count_matrix)= sort(unique(meta[,2]))
    
    #here we disable this file generation
    # all_sum = rowSums(count_matrix)
    # all_sum = cbind(names(all_sum), all_sum)
    # write.table(all_sum, file='interactions_sum.txt', quote=F, sep='\t', row.names=F)
    
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    
    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
    
    pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}






Output.Cytoscape = function(meta_file, pvalues_file, meta_sep='\t', pvalues_sep='\t', pvalue=0.05,cluster.sel = NA,network_filename,intr_filename){
  
  #######   Network
  
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  
  dim(all_intr)
  
  all_intr = all_intr %>% filter(!str_detect(all_intr$interacting_pair,"^COL"))
  
  dim(all_intr)
  
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  ## Changhua Yu: A small module for cluster selection
  if (!is.na(cluster.sel)) { 
    pairs1_all = cluster.sel
  }
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  pairs2 = c()
  for (i in 1:length(pairs1_all))
    for (j in i:length(pairs1_all))
      pairs2 = c(pairs2,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs2))
  {
    p1 = strsplit(pairs2[i], split_sep)[[1]][1]
    p2 = strsplit(pairs2[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs2[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  # here we disable this file generation
  all_count = all_count[-1,]
  all_count = as.data.frame(all_count)
  all_count[,"count"] = as.numeric(all_count[,"count"])
  all_count = all_count %>% arrange(-count)
  all_count = all_count[seq(2,nrow(all_count),2),]
  # Changhua Yu: for selecting the top interacting groups
  # col.sel = paste0(all_count[1:50,1], '|',all_count[1:50,2])
  
  write.table(all_count, network_filename, sep='\t', quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(pairs1_all), ncol=length(pairs1_all))
    rownames(count_matrix)= pairs1_all
    colnames(count_matrix)= pairs1_all
    
    #here we disable this file generation
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    write.table(all_sum, file=intr_filename, quote=F, sep='\t', row.names=F)
    
  }
}

