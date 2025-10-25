consensusClustR

[Consensus clustering](https://link.springer.com/article/10.1023/A:1023949509487) is a resampling-based method for discovering robust sample or feature clusters like stable patient subtypes and their (molecular) signatures. It addresses challenges in traditional clustering such as determining the correct number of clusters and assessing their stability.  

**Consensus Matrix, CDF, and Delta Area Plots**: To identify the optimal number of clusters **K**, consensus clustering evaluates the stability of sample groupings across repeated subsampling. The consensus matrix shows how frequently each pair of samples/features is clustered together across iterations, with values close to 1 indicating highly stable associations. From this matrix, a Cumulative Distribution Function (CDF) plot is generated, summarizing the distribution of consensus values for each **K**. The delta area plot calculates the proportional increase in the area under the CDF curve as **K** increases. Together, these plots help find the “elbow point”, the value of **K** where adding more clusters yields diminishing returns, indicating that further clusters mostly capture noise rather than meaningful structure.

---  

- This repo is an R-based workflow for performing consensus clustering using the package `ConsensusClusterPlus`:
  - locally, for smaller expression matrices, or
  - on an HPC cluster for large (genomic) data.  

**Use the `consensusWorkflow.rmd` file to for instructions on running consensus clustering analyses.**

