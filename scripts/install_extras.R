pkgs_cran <- c("Seurat", "BiocManager")
to_install <- pkgs_cran[!pkgs_cran %in% rownames(installed.packages())]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")

# CellChat from GitHub
if(!"CellChat" %in% rownames(installed.packages())){
  remotes::install_github("sqjin/CellChat")
}

# if(!"SCENIC" %in% rownames(installed.packages())){
#   BiocManager::install(c("RcisTarget","GENIE3","AUCell"))
#   remotes::install_github("aertslab/SCENIC")
# }
# if(!"SCENT" %in% rownames(installed.packages())){
#   remotes::install_github("aet21/SCENT")
# }
