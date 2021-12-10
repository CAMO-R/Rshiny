rm(list=ls())
#download Rshiny and save it as a directory to the following path
setwd("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/CAMO (WEZ97@pitt.edu)/Github_CAMO-R/")

installed <- installed.packages()[,"Package"]
# for (package in c("utils", "devtools", "shinyBS", "cluster", "parallel", "igraph",
#                   "Rcpp", "RcppArmadillo", "RcppGSL", "snowfall", "ggplot2",
#                   "gplots", "shinyjs","dplyr", "plotly","gridExtra", "reticulate")) {
#   if (!(package %in% installed)) {
#     install.packages(package, repos='http://cran.us.r-project.org')
#   }
# }
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# if (!("cvTools" %in% installed)) {
#   BiocManager::install("cvTools")
# }
# if (!("samr" %in% installed)) {
#   BiocManager::install("samr")
# }
# if (!("ggrepel" %in% installed)) {
#   BiocManager::install("ggrepel")
# }
# if (!("biomaRt" %in% installed)) {
#   BiocManager::install("biomaRt")
# }
# if (!("RcppGSL" %in% installed)) {
#   BiocManager::install("RcppGSL")
# }
# if (!("AnnotationDbi" %in% installed)) {
#   BiocManager::install("AnnotationDbi")
# }
# if (!("limma" %in% installed)) {
#   BiocManager::install("limma")
# }
# if (!("DESeq2" %in% installed)) {
#   BiocManager::install("DESeq2")
# }
# if (!("KEGGREST" %in% installed)) {
#   BiocManager::install("KEGGREST")
# }
# if (!("org.Hs.eg.db" %in% installed)) {
#   BiocManager::install("org.Hs.eg.db")
# }
# if (!("pathview" %in% installed)) {
#   BiocManager::install("pathview")
# }
# if (!("ConsensusClusterPlus" %in% installed)) {
#   BiocManager::install("ConsensusClusterPlus")
# }
# if (!("KEGGgraph" %in% installed)) {
#   BiocManager::install("KEGGgraph")
# }
# if (!("reactome.db" %in% installed)) {
#   BiocManager::install("reactome.db")
# }
# if (!("DT" %in% installed)) {
#   BiocManager::install("DT")
# }

# # install source packages
# if(!("DMwR" %in% installed)){# (deprecated and found in https://cran.r-project.org/src/contrib/Archive/DMwR/)
#   install.packages(c("zoo","xts","quantmod", "abind", "ROCR"))
#   install.packages("DMwR_0.4.1.tar.gz",repos=NULL, type="source")
# }
# 
# if (!("preproc" %in% installed)) {
#   install_github("metaOmic/preproc")
# }
# 
# if (!("KEGG.db" %in% installed)) {
#   install.packages("KEGG.db_3.2.3.tar.gz",repos=NULL, type="source")
# }

if(!("CAMO" %in% installed)){
  install.packages("CAMO_1.0.tar.gz",repos=NULL, type="source")
}
library(CAMO)
# library(DMwR)
# library(devtools)
# library(dplyr)
# library(shiny)
# library(shinyjs)
# library(shinyBS)
# library(cluster)
# library("org.Hs.eg.db")
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(CAMO)
# library(KEGGgraph)
# library(KEGGREST)
# library(igraph)
# library(parallel)
# library(reticulate)

# Settings maybe required for mac system:
#Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#Sys.setenv("RETICULATE_PYTHON" = "/home/rstudio/.local/share/r-miniconda/envs/r-reticulate/bin/python")
#use_python(RETICULATE_PYTHON)
#path_to_python <- "/home/rstudio/.local/share/r-miniconda/envs/r-reticulate/bin/python"
#use_python(path_to_python)

### create a virtual env for python
#conda_create("r-reticulate")
#use_condaenv("r-reticulate")
if(!py_module_available("pillow")) {
  py_install("pillow", pip = TRUE)
}
fh <- import("PIL")

if(!py_module_available("pandas")) {
  py_install("pandas", pip = TRUE)
}
fh <- import("pandas")

shiny::runApp('Rshiny-main', port=9987, launch.browser=T)


