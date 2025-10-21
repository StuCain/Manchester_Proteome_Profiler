
# Base image with R, Shiny Server, and tidyverse preinstalled
FROM rocker/shiny-verse:4.5.0

# ---- Metadata Labels ----
LABEL org.opencontainers.image.title="Manchester Proteome Profiler"
LABEL org.opencontainers.image.description="R Shiny app for downstream analysis of quantitative proteomics data. Supports MaxQuant, FragPipe, Proteome Discoverer, and Spectronaut. Includes normalization, imputation, clustering, GO enrichment, and protein interaction analysis."
LABEL org.opencontainers.image.version="1.0"
LABEL org.opencontainers.image.url="https://github.com/StuCain/Manchester_Proteome_Profiler"
LABEL org.opencontainers.image.source="https://github.com/StuCain/Manchester_Proteome_Profiler"
LABEL org.opencontainers.image.documentation="https://github.com/StuCain/Manchester_Proteome_Profiler"
LABEL org.opencontainers.image.authors="Stuart Cain"

# ---- System dependencies ----
RUN apt-get update && apt-get install -y --no-install-recommends     sudo     pandoc     libcurl4-gnutls-dev     libcairo2-dev     libxt-dev     libssl-dev     libssh2-1-dev     libbz2-dev     libhdf5-dev     build-essential     libgd-dev     libudunits2-dev     libproj-dev     libgdal-dev     libglpk-dev     libnetcdf-dev

# ---- CRAN package installation ----
RUN Rscript -e 'install.packages(c("devtools", "tidyverse", "ggrepel", "httr", "rjson", "mvtnorm", "tmvtnorm", "DT", "BiocManager", "shiny", "shinyalert", "shinydashboard", "shinyjs", "svglite", "shinyBS", "shinyWidgets", "shinycssloaders", "shiny.info", "circlize", "data.table", "scales", "gridtext", "waiter", "shinythemes", "RColorBrewer", "readxl", "rmarkdown", "tinytex", "assertthat", "fdrtool", "promises", "future", "ggtangle", "tidygraph", "ggraph"), dependencies=TRUE)'

# ---- Bioconductor package installation ----
RUN Rscript -e "BiocManager::install(c('DEP', 'ComplexHeatmap', 'EnhancedVolcano', 'org.Hs.eg.db', 'org.Mm.eg.db', 'HDO.db', 'VennDetail', 'clusterProfiler', 'QFeatures', 'org.At.tair.db', 'org.Xl.eg.db', 'STRINGdb'), ask=FALSE, dependencies=TRUE)"
RUN Rscript -e "BiocManager::install('preprocessCore', configure.args='--disable-threading', force=TRUE, update=TRUE, type='source')"

# ---- GitHub installs ----
RUN Rscript -e 'devtools::install_github("npjc/biogridr")'
RUN Rscript -e 'devtools::install_github("guokai8/VennDetail")'

# ---- TinyTeX installation ----
RUN Rscript -e "install.packages('tinytex', repos='https://cloud.r-project.org'); tinytex::install_tinytex(force = TRUE)"
ENV PATH="/root/.TinyTeX/bin/x86_64-linux:$PATH"
RUN tlmgr update --self
RUN Rscript -e "tinytex::is_tinytex()"
RUN tlmgr install sourcesanspro ly1 bookmark
RUN cp -r /root/.TinyTeX /home/shiny

# ---- Clean up default index and copy app ----
RUN sudo rm /srv/shiny-server/index.html

# Expose Shiny app on port 3838
EXPOSE 3838

# Copy app code and custom Shiny server configuration
COPY ./ /srv/shiny-server/shiny-mpp
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Set permissions for the Shiny app
RUN chmod -R +r /srv/shiny-server/shiny-mpp
RUN chown -R shiny:shiny /srv/shiny-server
