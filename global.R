# This file will be executed prior to app startup to setup the necessary environment
installed <- installed.packages()[,"Package"]

#CRAN packages:
# for (package in c("utils", "devtools", "shinyBS", "cluster", "parallel", "igraph",
#                   "Rcpp", "RcppArmadillo", "RcppGSL", "ggplot2",
#                   "gplots", "shinyjs","dplyr","gridExtra", "reticulate")) {
#   if (!(package %in% installed)) {
#     install.packages(package, repos='http://cran.us.r-project.org')
#   }
# }
for (package in c("utils", "devtools", "DT", "shinyBS","ggrepel")) {
  if (!(package %in% installed)) {
    install.packages(package, repos='http://cran.us.r-project.org')
  }
}

#Bioconductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# for (package in c("utils", "devtools", "shinyBS", "cluster", "parallel", "igraph",
#                   "Rcpp", "RcppArmadillo", "RcppGSL", "ggplot2",
#                   "gplots", "shinyjs","dplyr","gridExtra", "reticulate","ggrepel")) {
#   if (!(package %in% installed)) {
#     BiocManager::install(package)
#   }
# }
# if(!("DMwR" %in% installed)){# (deprecated and found in https://cran.r-project.org/src/contrib/Archive/DMwR/)
#   remotes::install_url("http://cran.nexr.com/src/contrib/DMwR_0.4.1.tar.gz")
# }
# if (!("preproc" %in% installed)) {
#   devtools::install_github("metaOmic/preproc")
# }
if(!("CAMO" %in% installed)){
  stop("Please install CAMO R package first")
}

# library(preproc)
library(shiny)
library(shinyBS)
library(shinyjs)
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
# library(snowfall)
# library(cvTools)
# library(samr)
library(ggrepel)
# library(limma) 
# library(MASS)
# library(ggrepel)
# library(ggplot2)
# library(gplots)
# library(WGCNA)
# library(tightClust)
# library(ConsensusClusterPlus)
# library(biomaRt)
# library(KEGG.db)
library(KEGGREST)
library(KEGGgraph)
# library(org.Hs.eg.db)
# library(pathview)
# library(cluster)
# library(gridExtra)
# library(grid)
# library(pathview)
# library(plotly)
library(igraph)
# library(parallel)

#library(DMwR)
library(devtools)
library(dplyr)
library(cluster)
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(devtools)
# library(KEGGgraph)
# library(KEGGREST)
# library(igraph)
# library(parallel)
library(reticulate) 

if(!py_module_available("pillow")) {
  py_install("pillow", pip = TRUE)
}
fh <- import("PIL")

if(!py_module_available("pandas")) {
  py_install("pandas", pip = TRUE)
}
fh <- import("pandas")

###############
library(CAMO)  # main pkg

#Include all global functions
dir <- "global"
for (f in list.files(path=dir)) {
  source(paste(dir, f, sep="/"))
}
# for (f in list.files(path=dir)) {
#   sourceCpp(paste(dir, f, sep="/"))
# }

# Create the directory for database prior to application startup
db <- new("Database", name="studies")

# Include all server modules
dir <- "server"
for (f in list.files(path=dir)) {
  source(paste(dir, f, sep="/"))
}

# Include all UI modules
dir <- "ui"
for (f in list.files(path=dir)) {
  source(paste(dir, f, sep="/"))
}

# Setting default working directory
tryCatch({
  DB.load.working.dir(db)
}, error=function(error){
  DB.set.working.dir(db, paste(getwd(), "data", sep="/"))
})

