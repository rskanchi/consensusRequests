# the Rscript to initiate consensus clustering
source("performConsensusClustering.R")
# provide expr file name (.csv) using the argument exprFile
performConsensusClustering(
  exprFile, # expression data file with sample rows and feature columns
  exprRowsSamples = TRUE, # if expression data rows are features, change to FALSE
  outputFolder = "output", # output folder name
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
  sample.dist.method = "spearman", # "pearson" = (1-Pearson cor), "spearman" (1-Spearman cor), "euclidean", "binary", "maximum", "canberra", "minkowski" or custom distance function
  expr.dist.method = "spearman",
  hclust.inner.linkage = "ward.D", # hierarchical linkage method for subsampling; "single", "complete", "average", "centroid"
  hclust.final.linkage = "ward.D", # hierarchical linkage method for consensus matrix
  corUse = "pairwise.complete.obs", # how to handle missing data in correlation distances
  seed = 7654321, # seed for generating resampling indices
)