# Manchester Proteome Profiler
Manchester Proteome Profiler (MPP) is an open-source R Shiny application that streamlines downstream analysis of quantitative proteomics data. Compatible with grouped protein intensities tables from MaxQuant, FragPipe, Proteome Discoverer, or other custom layouts, MPP provides an integrated platform for filtering, normalization, imputation, differential expression analysis, and cluster analysis across user chosen experimental conditions. It supports both single- and dual-dataset comparisons, incorporates SAINTexpress for affinity purification and BioID experiments, and downstream analysis of the significant protein list clusters to functional enrichment and interaction networks via Gene Ontology, BioGRID, and STRING. 
Manchester Proteome Profiler is a web application which is easy to use and allows users to upload grouped protein intensity data from MaxQuant, FragPipe or Proteome Discoverer. Following filtering, normalisation and imputation of missing values, Manchester Proteome Profiler uses a number of R packages to perform cluster analysis to allow data visualization and discovery of significantly enriched proteins. The suite also performs correlation and Principal Component Analysis (PCA) on the sample runs to evaluate replicate quality. Complex Heatmaps can be interactively generated from chosen condition comparisons, fold change and p.adj values. K means clustering allows proteins of similar fold change values to be grouped and subjected to GO ontology enrichment analysis. Enhanced Volcano plots and box plots are also generated to validate possible proteins of interest.  The suite also has built in SAINTexpress analysis, which allows identification proteins that have a higher abundance compared to the control samples.

<img
  src="app/www/Introduction.png"
  alt="Alt text"
  title="Optional title"
  style="display: inline-block; margin: 0 auto; max-width: 100px">
  

## Table of Contents
1. [Quick Instructions](#quick-instructions)
2. [Installation](#installation)


# Quick Instructions
#### 1. Load data (Pick Data File type first)
#### 2. Save Exp. Template
#### 3. Modify Exp. Template (Use a spreadsheet editor)
#### 4. Load modified Exp. Template
#### 5. Pick Comparisons or Generate comparison file from Comp. template and upload
#### 6. Press "Analyse"

There is currently a public server hosted by the University of Manchester available at https://mpp.sbs.manchester.ac.uk/ where more detailed instructions can also be found.
------------

## Data Pipeline
<img
  src="app/www/data_flow.jpg"
  alt="Alt text"
  title="Optional title"
  style="display: inline-block; margin: 0 auto; max-width: 100px">



## Installation
Once all the prerequisites are installed, follow steps below to build and run the server locally.

You can build it natively using the following code or use the renv_MPP script:

``` sh
# Clone the repository
git clone https://github.com/StuCain/Manchester_Proteome_Profiler

# Move to the folder
cd Manchester_Proteome_Profiler

# Install renv and BiocManager from CRAN
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Initialize renv project (creates renv.lock and renv library folder)
renv::init(bare = TRUE)

# Install CRAN packages using renv
cran_pkgs <- c(
  # Core development and utilities
  "devtools", "assertthat", "fdrtool",

  # Shiny framework & UI
  "shiny", "shinydashboard", "shinyWidgets", "shinythemes",
  "shinyalert", "shinyjs", "shinyBS", "shinycssloaders",
  "waiter", "DT",

  # Data handling and plotting
  "tidyverse", "data.table", "scales", "RColorBrewer",
  "readxl", "rmarkdown", "svglite", "gridtext", "ggrepel",
  "circlize", "httr", "rjson", "mvtnorm", "tmvtnorm", "ggtangle",
  "tidygraph", "ggraph",

  # Asynchronous and performance tools
  "promises", "future",

  # Reporting and LaTeX
  "tinytex"
)

renv::install(cran_pkgs)

# Install Bioconductor packages using renv
bioc_pkgs <- c(
  "DEP", "ComplexHeatmap", "clusterProfiler", "QFeatures",
  "EnhancedVolcano", "HDO.db", "org.Hs.eg.db", "org.Mm.eg.db",
  "VennDetail", "STRINGdb",
  "org.At.tair.db", "org.Xl.eg.db"
)
for (pkg in bioc_pkgs) renv::install(paste0("bioc::", pkg))

# Install GitHub Packages packages using renv
if (!requireNamespace("devtools", quietly = TRUE)) renv::install("devtools")
renv::install("npjc/biogridr")
renv::install("guokai8/VennDetail")


# ---- Snapshot the environment ----
# This records exact package versions into renv.lock
renv::snapshot()

# ---- Optional: TinyTeX setup ----
# Ensures LaTeX-based PDF reporting works as in Docker
if (!tinytex::is_tinytex()) {
  tinytex::install_tinytex(force = TRUE)
  system("tlmgr install sourcesanspro ly1 bookmark")
}
# Execute
library("shiny")
shiny::runApp("app/app.R")
```

Or build a Docker image

``` sh
# Clone the repository
git clone https://github.com/StuCain/Manchester_Proteome_Profiler
# Move to the folder
cd Manchester_Proteome_Profiler

# Build Manchester_Proteome_Profiler (Any name after -t)
docker buildx build -f Dockerfile -t manchester_proteome_profiler  --output=type=docker --platform=linux/amd64 .

# Run Manchester Proteome Profiler
docker run -it --platform=linux/amd64 -d -p 3838:3838 manchester_proteome_profiler

# Open local interface
http://localhost:3838/
```

The Manchester Proteome Profiler (MPP) Docker image is available on Docker Hub.
To download (pull) the latest version:

``` sh
# Pull the repository
docker pull stucain/mpp

# Run Manchester Proteome Profiler
docker run -d -p 3838:3838 stucain/mpp

# Open local interface
http://localhost:3838/
```