# ============================================================
# Manchester Proteome Profiler (MPP)
# ============================================================

# ---- Dependency managers ----
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

renv::init(bare = TRUE)

# ---- CRAN Packages ----
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

# ---- Bioconductor Packages ----
bioc_pkgs <- c(
  "DEP", "ComplexHeatmap", "clusterProfiler", "QFeatures",
  "EnhancedVolcano", "HDO.db", "org.Hs.eg.db", "org.Mm.eg.db",
  "VennDetail", "STRINGdb",
  "org.At.tair.db", "org.Xl.eg.db"
)
for (pkg in bioc_pkgs) renv::install(paste0("bioc::", pkg))

# ---- GitHub Packages ----
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

shiny::runApp("app/app.R")