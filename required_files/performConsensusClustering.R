# Consensus clustering alleviates common issues that arise in most clustering methods, 
# such as random initialization, choosing K, intuitive visualization and assessing stability of clusters. 
# This method gathers a consensus of cluster assignments based on sub-sampling the dataset and 
# running a chosen clustering algorithm multiple times. 

performConsensusClustering <- function(exprFile, # expression data file with sample rows and feature columns
                                       exprRowsSamples = TRUE, # if expression data rows are features, change to FALSE
                                       outputFolder, # provide project name
                                       traitFile = NULL, # optional traits/meta data, for heatmap top annotation, rows as samples
                                       discreteTraits = NULL, # c("gender", "disease") the discrete traits from the meta data to include 
                                       continuousTraits = NULL, # c("age", "BMI") the continuous traits from the meta data to include
                                       na.list = NULL,
                                       scaleExprData = TRUE, # column-wise scale expression data 
                                       nVarexpr = NULL, # optional number of highly variable expresion data features to include ex. 2000
                                       annotColor = NULL, # colors for annotation of discrete traits
                                       # clustering options for rows/samples
                                       maxK.sample = 8, # number of sample clusters to be evaluated, from 2 to maxK.sample
                                       maxK.expr = 8, # number of expr clusters to be evaluated, from 2 to maxK.expr
                                       nReps = 500, # number of subsamples
                                       pItem = 0.9, # proportion of items to subsample
                                       pFeature = 1, # proportion of features to subsample
                                       clusterAlg = "hc", # hc = hierarchical, "pam" paritioning around medoids, "km" k-means, or custom function
                                       sample.dist.method = "euclidean", # "pearson" = (1-Pearson cor), "spearman" (1-Spearman cor), "euclidean", "binary", "maximum", "canberra", "minkowski" or custom distance function
                                       expr.dist.method = "spearman",
                                       hclust.inner.linkage = "ward.D", # hierarchical linkage method for subsampling; "single", "complete", "average", "centroid"
                                       hclust.final.linkage = "ward.D", # hierarchical linkage method for consensus matrix
                                       corUse = "pairwise.complete.obs", # how to handle missing data in correlation distances
                                       seed = 7654321, # seed for generating resampling indices
                                       ...
                                       
){ # function code begins below
  
  # install the R packages not already installed and load all the required R packages
  list.of.packages <- c("ConsensusClusterPlus", "ggplot2", "clustree", "circlize",
                        "dendextend", "ggstatsplot", "gginnards")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  library(ComplexHeatmap)
  # library(tidyverse) check if required
  # Create output directory and start output log
  if (!dir.exists(outputFolder)) dir.create(path = outputFolder) # output folder
  #sink(file = paste(outputFolder, "consoleLog.Rout", sep = "/"), type = c("output", "message"))
  
  # Add default na.list to the user-specified NA strings
  na.list <- union(na.list,c("NA", "", "Don't know", "don't know", "unknown", "Unknown"))
  
  # Read data
  if (exprRowsSamples) data <- data.frame(read.csv(exprFile, na.strings = na.list, row.names = 1), check.names = F, check.rows = F) else
    data <- data.frame(t(read.csv(exprFile, na.strings = na.list, row.names = 1)), check.names = F, check.rows = F)
  
  cat(getTime(), " Number of samples (raw): ", nrow(data),"\n", sep = "")
  cat(getTime(), " Number of features (raw): ", ncol(data),"\n", sep = "")
  
  if(!is.null(traitFile)) {
    if(is.null(discreteTraits) & is.null(continuousTraits)){
      stop("For trait/meta data, both discrete and continuous traits are NULL. \n
         Please check the arguments discreteTraits and continuousTraits. \n
         At least one of these two must be non-null. \n")
    } else {
      traitData <- read.csv(traitFile, na.strings = na.list, row.names = 1)
      traitData <- subset(traitData, select = c(discreteTraits, continuousTraits))
      if (!is.null(discreteTraits)){ 
        for (disc in discreteTraits) traitData[, disc] <- factor(traitData[, disc])
      }
      cat(getTime(), " Number of samples (raw) in trait/meta data: ", nrow(traitData),"\n", sep = "")
      cat(getTime(), " Number of annotation columns in trait/meta data: ", ncol(traitData),"\n", sep = "")
    } # end of else: is.null(discreteTraits) & is.null(continuousTraits)
  } # end of if: !is.null(traitFile)
  
  
  # preprocessing
  
  # Remove features unexpressed across all samples
  keep <- colSums(data, na.rm = TRUE) != 0
  data <- data[,keep]
  # Remove samples with no expression across all features
  keep <- rowSums(data, na.rm = TRUE) != 0
  data <- data[keep,]
  
  # Remove features with zero variance
  keep <- apply(data, 2, var, na.rm = TRUE) != 0
  data <- data[,keep]
  
  # If meta data is available, restrict expression and meta data to common samples
  commonSamples <- intersect(rownames(data), rownames(traitData))
  data <- data[commonSamples, ]
  traitData <- traitData[commonSamples, ,drop=F]
  
  
  # check expression and clinical data samples are in the same order
  #identical(rownames(traitData), rownames(data))
  cat(getTime(), " Number of samples (after preprocessing): ",nrow(data), "\n", sep = "")
  cat(getTime(), " Number of features (after preprocessing): ",ncol(data), "\n", sep = "")
  
  # include the nVarexpr highly variable expression features
  if (is.null(nVarexpr)) {nVarexpr <- ncol(data)} else {
    nVarexpr <- min(nVarexpr, ncol(data)) # all features will be included if the nVarexpr is greater than the number of features in the data
  }
  
  sdExpr <- apply(data, 2, sd)
  sdExpr <- sort(sdExpr, decreasing = TRUE)
  hvExpr <- names(sdExpr)[1:nVarexpr]
  data <- data[, hvExpr] # Get the expression data for the genes of interest (top 2000, 3000, etc.)
  cat(getTime(), ncol(data), " Highly variable features (based on sd) included in analysis.", "\n", sep = "")
  
  # change output folder to the highly variable features
  # outputFolder <- paste(outputFolder, paste0("nFeatures_", nVarexpr, "_nReps_", nReps), sep = "/")
  outputFolder <- paste(outputFolder, paste0("nFeatures_", nVarexpr, "_nReps_", nReps, "_maxK", maxK.sample, "_", maxK.expr), sep = "/")
  
  if (!dir.exists(outputFolder)) dir.create(outputFolder) # this will be the default folder from now on
  
  # save data files
  if (!dir.exists(paste(outputFolder, "data", sep = "/"))) dir.create(paste(outputFolder, "data", sep = "/")) # to save data file
  write.csv(data, file = paste(outputFolder, "data/expr_data_preprocessed.csv", sep = "/"), row.names = TRUE)
  write.csv(traitData, file = paste(outputFolder, "data/trait_data_preprocessed.csv", sep = "/"), row.names = TRUE)
  file.copy(from = exprFile, to = paste(outputFolder, "data/expr_data_raw.csv", sep = "/"))
  file.copy(from = traitFile, to = paste(outputFolder, "data/trait_data_raw.csv", sep = "/"))
  # scale data if scaleExprData = TRUE
  if (scaleExprData) {
    data <- data.frame(scale(data), check.names = F)
    write.csv(data, paste(outputFolder, "data/expr_data_scaled.csv",sep = "/"), row.names = TRUE)
  } 
  
  # consensus clustering - samples
  cat(getTime(), " Sample consensus clustering in progress..", "\n", sep = "")
  sampleConsensus <- getConsensusClusters(d = t(data), dist.method = sample.dist.method, clusterAlg = clusterAlg, 
                                          innerLinkage = hclust.inner.linkage, finalLinkage = hclust.final.linkage,
                                          maxK = maxK.sample, nReps = nReps, pItem = pItem, pFeature = pFeature, 
                                          corUse = corUse, seed = seed,
                                          folder = paste(outputFolder, "sampleConsensus", sep = "/"))
  cat(getTime(), " Sample consensus clustering DONE", "\n", sep = "")
  
  # consensus clustering - features
  cat(getTime(), " Expression consensus clustering in progress..", "\n", sep = "")
  exprConsensus <- getConsensusClusters(d = data, dist.method = expr.dist.method, clusterAlg = clusterAlg, 
                                        innerLinkage = hclust.inner.linkage, finalLinkage = hclust.final.linkage,
                                        maxK = maxK.expr, nReps = nReps, pItem = pItem, pFeature = pFeature,
                                        corUse = corUse, seed = seed,
                                        folder = paste(outputFolder, "exprConsensus", sep = "/"))
  cat(getTime(), " Sample consensus clustering DONE", "\n", sep = "")
  
  # heatmaps
  cat(getTime(), " Generating heatmaps..", "\n", sep = "")
  dir.create(paste(outputFolder, "heatmaps", sep = "/"))
  heat <- list()
  
  # depending on the scaling, take care of the color scheme
  #cellColor <- colorRamp2(c(min(data), median(data), max(data)), c("blue", "white", "red")) # all heatmaps will have the same color limits
  
  cellColor <- colorRamp2(c(min(data), 0, max(data)), c("blue", "white", "red")) # all heatmaps will have the same color limits

    
  for (kRow in 2:maxK.expr){
    for (kCol in 2:maxK.sample){ 
      heatmapName <- paste("sampleK_", kCol, "_exprK_", kRow, sep = "")
      heatfolder <- paste(outputFolder, "heatmaps",  heatmapName, sep = "/")
      dir.create(heatfolder)
      
      hcRow <- exprConsensus[[kRow]]$consensusTree
      hcCol <- sampleConsensus[[kCol]]$consensusTree
      
      
      # discrete and continuous trait differences across clusters
      traitData$cluster <- sampleConsensus[[kCol]]$consensusClass
      traitData$cluster <- factor(traitData$cluster)
      pdf(paste(heatfolder, paste(heatmapName, "_traitData_Violin_and_Bar_plots.pdf", sep =""), sep = "/"))
      # violin plots for continuous traits
      if(!is.null(continuousTraits)){
        lapply(continuousTraits, FUN = function(cont){
          temp.aov <- summary(aov(traitData[, cont] ~ traitData[, "cluster"] , data = traitData))
          print(
          ggplot(traitData, aes(x = traitData[, "cluster"], y = traitData[, cont])) +
            geom_violin(trim = FALSE, aes(fill = cluster)) +
            xlab(paste("ANOVA Pr(>F):", format(temp.aov[[1]]$`Pr(>F)`[1], scientific = TRUE, digits = 3))) + ylab(cont) +
            #          labs(title = paste("ANOVA Pr(>F):", format(temp.aov[[1]]$`Pr(>F)`[1], scientific = TRUE, digits = 3))) +
            scale_fill_manual(values = sampleConsensus[[kCol]]$clrs[[3]]) +
            theme_classic(20))
        }) # end of lapply
      } # end of continuousTraits
      
      # bars for discrete traits
      suppressMessages( # to suppress the scale fill message
        if (!is.null(discreteTraits)){
          lapply(discreteTraits, FUN = function(disc){
            temp.f.test <- fisher.test(table(traitData[, disc], traitData[, "cluster"]), simulate.p.value = TRUE)
            print(
            ggbarstats(data = traitData, x = !!disc, y = cluster,
                       results.subtitle = FALSE,
                       subtitle = paste0("Fisher's exact test p-value: ", format(temp.f.test$p.value, scientific = TRUE, digits = 3))) +
              theme(axis.text.x = element_text(angle = 90)) +
              scale_fill_manual(values = annotColor[[disc]], #  sampleConsensus[[kCol]]$clrs[[3]]
                                guide = guide_legend(reverse = TRUE)) +
              theme_classic(20))
          }) # end of lapply
        } # end of discreteTraits
      )
      
      dev.off()
      
      clust.col <- sampleConsensus[[kCol]]$clrs[[3]]
      names(clust.col) <- 1:length(clust.col)
      if (is.null(get0("annotColor"))) annotColor <- list("cluster" = clust.col) else annotColor$cluster <- clust.col
      
      annotBar <- HeatmapAnnotation(df = traitData, which = "column", col = annotColor,
                                    simple_anno_size = unit(0.3, "cm")) #, col = annotColor
      
      if (scaleExprData) dataRangeName <- "scaled (mean=0, var=1)" else dataRangeName <- "data range"
      pdf(paste(heatfolder, paste0(heatmapName, "_heatmap.pdf"), sep = "/"), width = 12)
      heat[[heatmapName]] <- draw(Heatmap(as.matrix(t(data)), 
                                          cluster_rows = hcRow, cluster_columns = hcCol,
                                          show_row_dend = TRUE, show_column_dend = TRUE,
                                          show_row_names = TRUE, show_column_names = TRUE,
                                          row_names_gp = gpar(fontsize = 1), column_names_gp = gpar(fontsize = 5),
                                          column_title = paste("Heatmap: consensus clustering (samples = ", nrow(data),", expr = ", ncol(data), ")", sep = ""),
                                          column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                                          col = cellColor, #colGradient(c("blue","white","red"),length=15),
                                          top_annotation = annotBar,
                                          show_heatmap_legend = TRUE,
                                          name = dataRangeName,
                                          use_raster = FALSE,
                                          row_split = kRow, column_split = kCol, 
      ), merge_legend = TRUE)
      dev.off()
      
      # save results 
      sampleClusters <- data.frame("consensusCluster" = sampleConsensus[[kCol]]$consensusClass)
      exprClusters <- data.frame("consensusCluster" = exprConsensus[[kRow]]$consensusClass)
      write.csv(sampleClusters, file = paste(heatfolder, "sampleConsensusClusters.csv", sep = "/"), row.names = TRUE)
      write.csv(exprClusters, file = paste(heatfolder, "exprConsensusClusters.csv", sep = "/"), row.names = TRUE)
      write.csv(rownames(data)[sampleConsensus[[kCol]]$consensusTree$order], file = paste(heatfolder, "heatColOrder.csv", sep = "/"), row.names = TRUE)
      write.csv(colnames(data)[exprConsensus[[kRow]]$consensusTree$order], file = paste(heatfolder, "heatRowOrder.csv", sep = "/"), row.names = TRUE)
      
      write.csv(data.frame("clusterColor" = sampleConsensus[[kCol]]$clrs[[3]]), file = paste(heatfolder, "sampleClusterColors.csv", sep = "/"), row.names = TRUE)
      write.csv(data.frame("clusterColor" = exprConsensus[[kRow]]$clrs[[3]]), file = paste(heatfolder, "exprClusterColors.csv", sep = "/"), row.names = TRUE)
    } # end of for in kCol
  } # end of for in kRow
  cat(getTime(), " Heatmaps DONE", "\n", sep = "")
  
  #sink()
} # end of function performConsensusClustering

# function to get consensus clusters of the columns of a dataframe
getConsensusClusters <- function(d, folder, dist.method = "spearman", clusterAlg = "hc", 
                                 innerLinkage = "ward.D", finalLinkage = "ward.D",
                                 maxK = 8, nReps = 500, pItem = 0.9, pFeature = 1, corUse = "pairwise.complete.obs", 
                                 seed = 7654321){
  library(ConsensusClusterPlus)
  library(clustree)
  library(ggplot2)
  
  if (!dir.exists(folder)){ dir.create(folder, recursive = TRUE) }
  
  if (dist.method == "spearman" | dist.method == "pearson")
    distMat <- as.dist(1-cor(d, use = corUse, method = dist.method)) else # cor() returns cor between cols of dataframe
      distMat <- as.dist(dist(t(d), method = dist.method)) # dist() returns dist between rows of the dataframe so t(d)
    
    write.csv(as.matrix(distMat), file = paste(folder, paste0("distMat_", dist.method, ".csv"), sep = "/"), row.names = TRUE)
    
    consensusRes <- ConsensusClusterPlus(distMat, maxK = maxK, reps = nReps, pItem = pItem, pFeature = pFeature,
                                         innerLinkage = innerLinkage, finalLinkage = finalLinkage, clusterAlg = clusterAlg, 
                                         title = folder, seed = seed, plot = "png")
    # cluster tree for expression/features
    df_clusters <- data.frame(row.names = colnames(d))
    df_clusters[, 1] <- rep(1, ncol(d))
    for (k in 2:maxK){ df_clusters[, k] <- consensusRes[[k]][["consensusClass"]] }
    colnames(df_clusters) <- paste("K", 1:maxK, sep = "")
    write.csv(df_clusters, paste(folder, paste0("clustree_", dist.method, ".csv"), sep = "/"), row.names = TRUE)
    
    #df_clusters <- df_clusters[, c(1:12, 18, 24)]
    clusterTree <- clustree(df_clusters, prefix = "K", node_label = "size", prop_filter = 0, node_text_size = 5, node_label_size = 5,
                            node_label_nudge = -0.4) + 
      guides(fill="none", x = "none", y = "none") + labs(x = NULL, y = NULL) +
      theme_classic(15) # use_core_edges = FALSE, edge_arrow = TRUE, edge_arrow_ends = "both"
    
    pdf(paste(folder, paste0("clustree_", dist.method, ".pdf"), sep = "/"), height = 10, width = 10)
    print(clusterTree) 
    dev.off()
    
    return(consensusRes)
    
} # end of function getConsensusClusters

# get time
getTime <- function(){paste0("[", format(Sys.time(), usetz = T), "]")}
