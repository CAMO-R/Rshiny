FROM rocker/shiny:latest

MAINTAINER Wei Zong "wez97@pitt.edu"

# Install dependencies and Download and install shiny server
RUN R -e "install.packages(c('shinyBS'))" && \ 
    R -e "install.packages(c('BiocManager','devtools','utils','shinyjs','RcppGSL','RcppArmadillo', 'Rcpp', 'MASS', 'parallel', 'methods', 'igraph', 'gridExtra', 'grid', 'ggplot2', 'gplots', 'reticulate', 'DT'), repos='https://cran.rstudio.com/')" && \ 
    R -e "BiocManager::install(c('DESeq2','limma','ConsensusClusterPlus','pathview','KEGGgraph','KEGGREST','org.Hs.eg.db','org.Mm.eg.db','org.Rn.eg.db','org.Dm.eg.db','reactome.db'), dependencies = TRUE)" && \
    R -e "devtools::install_github('https://github.com/CAMO-R/Rpackage')"

# Copy app and change directory permissions
COPY . /srv/shiny-server/CAMO
RUN chown -R shiny:shiny /srv/shiny-server/CAMO && \
    chown -R shiny:shiny /usr/local/lib/R/site-library && \
    mkdir /srv/shiny-server/CAMO/.Database && \
    chown -R shiny:shiny /srv/shiny-server/CAMO/.Database