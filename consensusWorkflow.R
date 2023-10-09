# Processing consensus clustering requests
# You are in the folder for consensusRequests with a folder named required files

##### Step 1: provide project name 
project_name <- "Lastname_ProjectId" # the folder name for the project

##### RUN the lines below until Step 2
dir.create(project_name)
file.copy("required_files/template.slurm.sh", to = project_name)
file.copy("required_files/template.R", to = project_name)
file.copy("required_files/performConsensusClustering.R", to = project_name) 
# the functions required are getTime(), getConsensusClusters(), performConsensusClustering()

setwd(project_name)
shFile <- readLines("template.slurm.sh", warn = FALSE)
shFile2 <- gsub(pattern = "project_name", replace = project_name, x = shFile)
writeLines(shFile2, con = paste0(project_name, ".slurm.sh"))
file.rename(from = "template.R", to = paste0(project_name, ".R"))

file.remove("template.slurm.sh")
rm(shFile, shFile2)

##### Step 2: Open project_name.R file from the project folder and provide values for arguments to perform consensus clustering
##### Step 3: login to the cluster
  # Copy the project folder to consensus requests folder, say /mount/a_cluster_folder/yourfolder/consensusRequests
  # cd to the project folder
  # sh the slurm file



