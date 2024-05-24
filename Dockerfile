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
"shiny.info", "circlize", "data.table", "scales", "gridtext", "waiter", "shinythemes", "RColorBrewer", "readxl", \
"rmarkdown", "tinytex"), dependencies=TRUE)'

#FROM bioconductor/bioconductor_docker
RUN Rscript -e 'BiocManager::install(pkgs=c("DEP", "ComplexHeatmap","EnhancedVolcano","org.Hs.eg.db","HDO.db"), ask=F, dependencies=TRUE)'
RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('org.Mm.eg.db')"
RUN R -e "BiocManager::install('VennDetail')"
RUN apt-get install perl wget -y  
#fetch TinyTeX installer and run it
RUN wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh 
#add TinyTeX to PATH variable
ENV PATH="$PATH:/root/.TinyTeX/bin/x86_64-linux"
#Install extra TinyTEX packages
RUN tlmgr install sourcesanspro
RUN tlmgr install ly1
RUN tlmgr install bookmark
#Copy TinyTEx directory from root to /home/shiny
RUN cp -r /root/.TinyTeX /home/shiny
RUN sudo rm /srv/shiny-server/index.html

# Make the ShinyApp available at port 80
EXPOSE 80
COPY ./ /srv/shiny-server/shiny-saint-dep
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

RUN chmod -R +r /srv/shiny-server/shiny-saint-dep

RUN chown -R shiny:shiny /srv/shiny-server