## Functions for Fig6 and Supp
library(fgsea)
library(RColorBrewer)
library(reactome.db)
library(STRINGdb)
library(topGO)
library(org.Hs.eg.db)
library(data.table)
library(stringr)
library(EnhancedVolcano)
library(data.table)
library(dplyr)
library(reshape)

set.seed(4)

sampleType_infectionGroup_cols2 = c('SAT:Mock' = "#339CFF", 'SAT:SARS-CoV-2' = "#201A80", 'VAT:Mock' = "#FFA833", 'VAT:SARS-CoV-2' = "#801E1A")

### functions below
fgseaOuts = function(path, prefix) {
  allfiles <- list.files(path, pattern = ".csv")
  df_files <- list()
  for(x in allfiles){
    onefile <- fread(paste0(path, x))
    onefile$clust <- sapply(strsplit(x, "[_]"), "[", 1)
    df_files[[x]] <- onefile
  } 
  all_cluster <- rbindlist(df_files)
  all_cluster <- all_cluster[which(all_cluster$p_val_adj < 0.1),]
  all_cluster_comb <- split(as.data.table(all_cluster), by = "clust")
  pathways = lapply(all_cluster_comb, runReactomePathway)
  
  pathways_collapsed = lapply(pathways, function(obj) {
    return(obj$collapsed_Path)
  })
  pathways_all = lapply(pathways, function(obj) {
    return(obj$fgseaMultilevelRes)
  })
  
  pathways_collapsed = add_clust(pathways_collapsed, prefix)
  pathways_all = add_clust(pathways_all, prefix)
  
  outputList = list()
  outputList[["pathways_collapsed"]] = pathways_collapsed
  outputList[["pathways_all"]] = pathways_all
  return(outputList)
}

getClusterPaths = function(inputDat, prefix) {
  pathways_collapsed = inputDat$pathways_collapsed
  pathways_all = inputDat$pathways_all
  
  pathways_collapsed_sig = subset_imp(pathways_collapsed)
  pathways_all_sig = subset_imp(pathways_all)
  
  ## collapsed
  collapsed_out = rbindlist(pathways_collapsed_sig)
  collapsed_out <- cast(collapsed_out, pathway~clust)
  collapsed_out[is.na(collapsed_out)] = 0
  rownames(collapsed_out) = collapsed_out$pathway
  collapsed_out$pathway = NULL
  
  donornames = lapply(ordered_adiposeclust, function(vec) {return(paste0(prefix,vec))}) %>% unlist()
  donornames = donornames[donornames %in% colnames(collapsed_out)]
  collapsed_out = collapsed_out[donornames]
  ## allPaths
  all_sig_out = rbindlist(pathways_all_sig)
  all_sig_out <- cast(all_sig_out, pathway~clust)
  all_sig_out[is.na(all_sig_out)] = 0
  rownames(all_sig_out) = all_sig_out$pathway
  all_sig_out$pathway = NULL
  
  donornames = lapply(ordered_adiposeclust, function(vec) {return(paste0(prefix,vec))}) %>% unlist()
  donornames = donornames[donornames %in% colnames(all_sig_out)]
  all_sig_out = all_sig_out[donornames]
  
  outputList = list()
  outputList[["collapsed_paths"]] = collapsed_out
  outputList[["all_paths"]] = all_sig_out
  return(outputList)
}


## add list name to each table (add cluster info)
add_clust = function(obj, prefix) {
  for (elem in names(obj)) {
    obj[[elem]]$clust = paste0(prefix,elem)
  }
  return(obj)
}


## get only significant pathways and maintain 3 columns
subset_imp = function(obj) {
  obj = lapply(obj, function(obj_one) {
    obj_one <- obj_one[obj_one$padj < 0.05,]
    obj_one = obj_one[,c("clust","pathway","NES")]
    return(obj_one)
  })
  return(obj)
}

### patwhay function that runs per dataset
runReactomePathway = function(DEG_table) {
  output = list()
  DEG_table$gene = DEG_table$V1
  Symbol2entrezIDs <- AnnotationDbi::select(x = org.Hs.eg.db, keys = DEG_table$gene, keytype = "SYMBOL", columns = "ENTREZID") %>% as.data.table
  Symbol2entrezIDs_unique = Symbol2entrezIDs[!duplicated(Symbol2entrezIDs$SYMBOL), ]
  
  ## have both dataframes in the same order by symbol
  DEG_table = DEG_table[order(DEG_table$gene), ]
  Symbol2entrezIDs = Symbol2entrezIDs_unique[order(Symbol2entrezIDs_unique$SYMBOL),]
  
  # check that gene symbols are in the same order
  if (sum(Symbol2entrezIDs$SYMBOL != DEG_table$gene) > 0) {
    print("ERROR!")
    invokeRestart("abort")
  }
  
  # grab ranks (logfc) & name it by geneid
  ranks <- DEG_table$avg_log2FC
  names(ranks) = Symbol2entrezIDs$ENTREZID
  ranks = sort(ranks, decreasing = TRUE)
  ranks = ranks[!is.na(names(ranks))]
  
  ## if not enough DEGS
  if (length(ranks)<5) {
    emptydf = setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("pathway", "pval", "padj","log2err","ES","NES","size","leadingEdge"))
    output[["fgseaMultilevelRes"]] = emptydf %>% as.data.table()
    output[["collapsed_Path"]] = emptydf %>% as.data.table()
    return(output)
  }
  
  # get reactome pathways
  pathways <- reactomePathways(names(ranks))
  
  # get enriched pathways
  if (length(ranks[ranks<0]) > 0) {
    fgseaMultilevelRes <- fgseaMultilevel(pathways = pathways, 
                                          stats = ranks,
                                          minSize=15,
                                          maxSize=500)
  }
  if (length(ranks[ranks<0]) == 0) {
    fgseaMultilevelRes <- fgseaMultilevel(pathways = pathways, 
                                          stats = ranks,
                                          minSize=15,
                                          maxSize=500, scoreType = "pos")
  }
  
  
  # grab significant only pathways
  fgseaMultilevelRes_sig = fgseaMultilevelRes[fgseaMultilevelRes$padj<0.001,]
  
  # collapse pathways to independent, unrelated pathways (~optional~ might dump this)
  collapsedPathways <- collapsePathways(fgseaMultilevelRes_sig[order(pval)][padj < 0.01], 
                                        pathways, ranks)
  
  # final path list
  finalPaths <- fgseaMultilevelRes_sig[pathway %in% collapsedPathways$mainPathways][order(-NES),]
  output[["fgseaMultilevelRes"]] = fgseaMultilevelRes
  output[["collapsed_Path"]] = finalPaths
  return(output)
}

## combine outputs
combineThree_collapsed = function(obj1, obj2, obj3){
  allObj= merge(obj1$collapsed_paths, obj2$collapsed_paths, by = 0, all = TRUE)
  rownames(allObj) = allObj$Row.names
  allObj$Row.names = NULL
  allObj = merge(allObj, obj3$collapsed_paths, by  = 0, all = TRUE)
  allObj[is.na(allObj)] = 0
  rownames(allObj) = allObj$Row.names
  allObj$Row.names = NULL
  return(allObj)
}

combineThree_all = function(obj1, obj2, obj3){
  allObj= merge(obj1$all_paths, obj2$all_paths, by = 0, all = TRUE)
  rownames(allObj) = allObj$Row.names
  allObj$Row.names = NULL
  allObj = merge(allObj, obj3$all_paths, by  = 0, all = TRUE)
  allObj[is.na(allObj)] = 0
  rownames(allObj) = allObj$Row.names
  allObj$Row.names = NULL
  return(allObj)
}