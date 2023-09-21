FROM rocker/shiny-verse:latest

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev\
    libbz2-dev \
    libhdf5-dev \
    build-essential \
    libgd-dev \
    libudunits2-dev \
    libproj-dev \
    libgdal-dev \
    libglpk-dev \
    libnetcdf-dev


#Needed for clusterProfiler
#RUN apt-get install libglpk-dev
#Needed for DEP
#RUN apt-get install libnetcdf-dev

RUN Rscript -e 'install.packages(c("devtools", "tidyverse", "ggrepel", "httr", "rjson", "mvtnorm", "tmvtnorm", \
"DT", "BiocManager", "shiny", "shinyalert", "shinydashboard", \
"shinyjs", "svglite", "shinyBS", "shinyWidgets", "shinycssloaders", \
"shiny.info", "circlize", "data.table", "scales", "gridtext", "waiter", "shinythemes", "RColorBrewer", "readxl"), dependencies=TRUE)'

#FROM bioconductor/bioconductor_docker:RELEASE_3_15
RUN Rscript -e 'BiocManager::install(pkgs=c("DEP", "ComplexHeatmap","EnhancedVolcano","org.Hs.eg.db","HDO.db"), ask=F, dependencies=TRUE)'
RUN R -e "BiocManager::install('clusterProfiler')"
RUN sudo rm /srv/shiny-server/index.html

# Make the ShinyApp available at port 80
EXPOSE 80
COPY ./ /srv/shiny-server/shiny-dep
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
#RUN rm -f /srv/shiny-server/lfq-analyst/.Rprofile
RUN chmod -R +r /srv/shiny-server/shiny-dep
#To fix "cannot open file 'rplots.pdf'" error
RUN chown -R shiny:shiny /srv/shiny-server