#Shiny_DEP_Heatmaps  by Stuart Cain
#The University of Manchester
#This script imports Maxquant output file proteinGroups.txt and combines it with a experimental
#design file. The DEP package is used to calculate foldchange values, which can be filtered using Cutoff values.
#The complexHeatmap package is used to render heatmaps, and the clusters are then analysed by clusterprofiler

#Use to install packages not listed on CRAN
#install.packages("BiocManager")
#BiocManager::install(c("EnhancedVolcano", "enrichplot",  "clusterProfiler",
#"ComplexHeatmap", "SummarizedExperiment", "DEP", "org.Hs.eg.db", "HDO.db", "VennDetail", "QFeatures"))
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.At.tair.db")
#BiocManager::install("org.Xl.eg.db")
#install.packages("devtools")  # if not already installed
#library(devtools)
#devtools::install_github("npjc/biogridr")


# Set a custom temp directory for Shiny Server
if (.Platform$OS.type == "unix") {
Sys.setenv(TMPDIR = "/srv/shiny-server/tmp")
}
library(DEP)
library(tidyverse) #contains: ggplot2, dplyr, tidyr, 
# readr, purr, tibble, stringr, forcats
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(scales)
library(gridtext)
library(EnhancedVolcano)
library(waiter)
library(shinythemes)
library("RColorBrewer")
library(readxl)
library(rmarkdown)
library(SummarizedExperiment)
library(QFeatures)
library(promises)
library(future)
library(ggtangle)
library(VennDetail)
library(ragg)
library(biogridr)
library(tidygraph)
library(ggraph)

#library(limma)
#library(MSnbase)
#library(assertthat)
source("R/functions.R")
#source("R/functions_DEP.R")
Version = "0.16.2-beta"

set.seed(123) #random seed to give consistent kmeans clustering and imputation

#Check TinyTex
options(shiny.sanitize.errors = FALSE)
tt_check<-tinytex:::is_tinytex()
if(tt_check != "TRUE"){tinytex::install_tinytex()}
data_other <- reactiveVal(NULL)

# Use multicore processing for Linux
if (.Platform$OS.type == "unix") {
  plan(multicore)
} else if (.Platform$OS.type == "windows") {
  plan(multisession)
}


####UI####
ui <- fluidPage(
  tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=Montserrat:wght@300;600&family=Open+Sans:wght@300;400&display=swap'); /* Example Google Font */
    body {
      font-family: 'Open Sans', sans-serif; /* Apply the font to the entire body */
      /* Other body styles if needed */
    }
    h1, h2, h3 { /* Apply to headings */
      font-family: 'Montserrat', sans-serif;
      font-weight: 600;
    }
    .custom-class { /* Apply to specific elements */
      font-family: 'Open Sans', sans-serif;
    }
  ")), # Close tags$style
  tags$style(HTML("
    #conversion_status {
      background-color: white;
      color: grey;
      font-size: 14px;
      font-family: monospace;
      #white-space: pre-wrap;
      #padding: 10px;
      #border: 1px solid #ccc;
      #border-radius: 5px;
    }
  ")),
  tags$head(includeHTML(("google-analytics.html"))),
  theme = shinytheme("spacelab"),
  tags$head(tags$style(HTML('.progress-bar {background-color: #2297E6;}'))),
  #tags$head(tags$style(HTML('.navbar-default {background-color: #2297E6; border-color: #dddddd}'))),
  tags$head(tags$link(rel = "icon", type = "image/gif/png", href = "app_icon.png")), 
  tags$head(tags$title("Manchester Proteome Profiler"),
            tags$meta(name = "description", content = "Manchester Proteome Profiler is an esy to use tool to allow analysis of protein intensity data from MaxQuant, FragPipe,  Proteome Discoverer or any grouped protein intensity data. MPP can then normalise, impute missing values and carry out pairwise comparisons. You can then check the data quality with PCA and correlation plots, generate heatmaps, volcano plots and box plots. The resulting significant proteins can be then clustered and subjected to Gene Ontology analysis.")
            ),
  tags$head(tags$style(".main-header .logo {height: 70px}, {background-color: #FFFFFF;}")),
  autoWaiter(html = spin_dots(), color = "white"),
  dashboardPage(
    dashboardHeader( # Set height of dashboardHeader
      title = tags$img(src='app_logo v7.jpg', width = "100%", height = "100%"),
      tags$li(class="dropdown", tags$img(src='CELL MATRIX RGB2.jpg', width = "40%", height = "40%", style="float:right")),
      titleWidth = 300
    ),
    dashboardSidebar(width=300,#changed with of sidebar
                     sidebarMenu(
                       tags$hr(), #adds a line
                       helpText("Two Dataset comparison mode"),
                       switchInput(
                         inputId = "data_mode",
                         size = "mini",
                         onStatus = "success", 
                         offStatus = "danger"
                       ),
                       
                       radioButtons("data_type",
                                    "Data File type",
                                    choices = c("MaxQuant", "Proteome Discoverer", "FragPipe", "Other (Experimental)"),
                                    #selected = "MaxQuant"
                       ),
                       conditionalPanel(
                         condition = "input.data_type == 'MaxQuant'",
                         fileInput('file1',
                                   'ProteinGroups.txt',
                                   accept=c('text/csv',
                                            'text/comma-separated-values,text/plain',
                                            '.csv'))
                       ),
                       conditionalPanel(
                         condition = "input.data_type == 'Proteome Discoverer'",
                         fileInput('filePD',
                                   'PD Output .xlsx',
                                   accept=c(".xlsx"))
                       ), 
                       conditionalPanel(
                         condition = "input.data_type == 'FragPipe'",
                         selectInput("exp", "Data type:", choices = c("LFQ", "DIA"), selected = "LFQ"),
                         conditionalPanel(
                           condition = "input.exp == 'LFQ'",
                           fileInput('lfq_expr', 'Upload FragPipe combined_protein.tsv',
                                     accept=c('text/tsv',
                                              'text/tab-separated-values,text/plain',
                                              '.tsv')),
                           radioButtons("lfq_type",
                                        "Intensity Type",
                                        choices = c("Intensity"="Intensity",
                                                    "MaxLFQ Intensity"="MaxLFQ"),
                                        #                                                            "Spectral Count"="Spectral Count"),
                                        selected = "MaxLFQ"),
                         ),
                         conditionalPanel(
                           condition = "input.exp == 'DIA'",
                           fileInput('dia_expr', 'Upload protein group (PG) matrix *.tsv',
                                     accept=c('text/tsv',
                                              'text/tab-separated-values,text/plain',
                                              '.tsv')),
                           radioButtons("dia_species",
                                        "Choose Species",
                                        choices = c("Human", "Mouse"),
                                        selected = "Human"),
                         ),
                       ),
                       conditionalPanel(
                         condition = "input.data_type == 'Other (Experimental)'",
                         fileInput('file_other',
                                   'Tab seperated text file',
                                   accept=c('text/tsv',
                                            'text/tab-separated-values,text/plain',
                                            '.tsv')),
                         menuItemOutput("columns")
                       ),
                       
                       
                       #                                menuItemOutput("file_input"),
                       #                                infoBoxOutput("ibox"), #Info box for testing
                       fileInput('file2',
                                 'ExperimentalTemplate.txt',
                                 accept=c('text/csv',
                                          'text/comma-separated-values,text/plain',
                                          '.csv')),
                       helpText("Generate Experimental template"),
                       #(=h3("Only one set of comparisons selected", style = 'font-size:16px;color:blue;')
                       downloadButton('downloadexp_template', 'Save Exp. template', class = "btn-primary"),
                       
                       #menuItemOutput("columns"),
                       conditionalPanel(
                         condition = "input.data_mode == false",
                         radioButtons("comp_type",
                                      "Comparisons",
                                      choices = c("1 Round of comparisons", "2 Rounds of comparisons", "Use SAINTexpress"),
                                      #selected = "1 Round of comparisons"
                         ),
                       ),
                       conditionalPanel(
                         condition = "input.data_mode == true",
                         radioButtons("data_type2",
                                      "Data File type 2",
                                      choices = c("MaxQuant", "Proteome Discoverer", "FragPipe", "Other (Experimental)"),
                         ),
                       ),
                       conditionalPanel(
                         condition = "input.data_mode == true & input.data_type2 == 'MaxQuant'",
                         fileInput('fileMax2',
                                   'ProteinGroups.txt',
                                   accept=c('text/csv',
                                            'text/comma-separated-values,text/plain',
                                            '.csv'))
                       ),
                       conditionalPanel(
                         condition = "input.data_mode == true &input.data_type2 == 'Proteome Discoverer'",
                         fileInput('filePD2',
                                   'PD Output .xlsx',
                                   accept=c(".xlsx"))
                       ), 
                       conditionalPanel(
                         condition = "input.data_mode == true & input.data_type2 == 'FragPipe'",
                         selectInput("exp2", "Data type:", choices = c("LFQ", "DIA"), selected = "LFQ"),
                         conditionalPanel(
                           condition = "input.exp2 == 'LFQ'",
                           fileInput('lfq_expr2', 'Upload FragPipe combined_protein.tsv',
                                     accept=c('text/tsv',
                                              'text/tab-separated-values,text/plain',
                                              '.tsv')),
                           radioButtons("lfq_type2",
                                        "Intensity Type",
                                        choices = c("Intensity"="Intensity",
                                                    "MaxLFQ Intensity"="MaxLFQ"),
                                        #                                                            "Spectral Count"="Spectral Count"),
                                        selected = "MaxLFQ"),
                         ),
                         conditionalPanel(
                           condition = "input.exp2 == 'DIA'",
                           fileInput('dia_expr2', 'Upload protein group (PG) matrix *.tsv',
                                     accept=c('text/tsv',
                                              'text/tab-separated-values,text/plain',
                                              '.tsv')),
                           radioButtons("dia_species2",
                                        "Choose Species",
                                        choices = c("Human", "Mouse"),
                                        selected = "Human"),
                         ),
                       ),
                       conditionalPanel(
                         condition = "input.data_type2 == 'Other (Experimental)'",
                         fileInput('file_other2',
                                   'Tab seperated text file',
                                   accept=c('text/tsv',
                                            'text/tab-separated-values,text/plain',
                                            '.tsv')),
                         menuItemOutput("columns2")
                       ),
                       conditionalPanel(
                         condition = "input.data_mode == true",
                         fileInput('file22',
                                   'ExperimentalTemplate.txt',
                                   accept=c('text/csv',
                                            'text/comma-separated-values,text/plain',
                                            '.csv')),
                         helpText("Generate Experimental template"),
                         #(=h3("Only one set of comparisons selected", style = 'font-size:16px;color:blue;')
                         downloadButton('downloadexp_template2', 'Save Exp. template', class = "btn-primary")
                         
                       ),
                       
                       menuItemOutput("comparisons"),
                       actionButton("analyze", "Analyze", icon = icon("magnifying-glass-chart"), class="btn btn-primary btn-lg"),
                       menuItemOutput("imputation_options"),
                       menuItemOutput("normalisations"),
                       menuItemOutput("padj_options"),
                       tags$hr(),
                       uiOutput("downloadTable"),
                       uiOutput("downloadButton"),
                       tags$hr(), #adds a line
                       helpText("Supported by"),
                       imageOutput("logo_img", height = "100px"),
                       helpText("Hosted by"),
                       imageOutput("logo_wccr",  height = "100px")
                       
                     )
    ),
    
    #Body
    dashboardBody(
      #      tags$head(tags$style(HTML('.box {margin: 5px;}' ))),
      #      tags$head(tags$style(HTML('.box {margin-left: 0px;}' ))),
      #      tags$head(tags$style(HTML('.box {margin-right: 0px;}' ))),
      helpText(paste("By Stuart Cain, Version:", Version)),
      #set infobox height
      tags$head(tags$style(HTML('.info-box {min-height: 120px;} .info-box-icon {height: 120px; line-height: 120px;} .info-box-content {padding-top: 20px; padding-bottom: 20px;}'))),
      #Google Analytics
      tags$head(includeHTML(("google-analytics.html"))),
      #Top row of boxes
      fluidRow(
        box(numericInput("thr",
                         "Missing value threshold",
                         min = 0, max = 3, value = 1),
            width = 1, height = 140 ),
        box(numericInput("p",
                         "adj. P value",
                         min = 0.0001, max = 0.1, value = 0.05),
            checkboxInput("pval_use", "Use p.Val", FALSE),
            width = 1, height = 140),
        box(numericInput("lfc",
                         "Fold change",
                         min = 0, max = 32, value = 2),
            width = 1, height = 140),
        infoBoxOutput("significantBox"), #shows number of Sig. Proteins
        conditionalPanel(
          condition = "input.comp_type == 'Use SAINTexpress'",
          box(numericInput("bfdr",
                           "BFDR Cut-Off",
                           min = 0, max = 1, value = 0.05, step = 0.05),
              width = 1, height = 140),
          box(numericInput("bfdr_freq",
                           "BFDR Freq. Min",
                           min = 0, max = 10, value = 1),
              width = 1, height = 140),
        ),
        conditionalPanel(
          condition = "input.data_mode == true || input.comp_type == '2 Rounds of comparisons'",
          box(numericInput("p2",
                           "adj. P value (Comp2)",
                           min = 0.0001, max = 0.1, value = 0.05),
              width = 1, height = 140),
          box(numericInput("lfc2",
                           "Fold change (Comp2)",
                           min = 0, max = 32, value = 2),
              width = 1, height = 140),
        ),
        conditionalPanel(
          condition = "input.analyze >= 1",
          box(
            radioButtons("report_type",
                         "Report type",
                         choices = c("Short", "with Box Plots", "Long"),
                         selected = "Short"),
            width = 2, height = 140),
        ),
        conditionalPanel(
          condition = "input.data_mode == true || input.comp_type == '2 Rounds of comparisons'",
          infoBoxOutput("significantBox2"), #shows number of Sig. Proteins
          infoBoxOutput("significantBox3"), #shows number of Sig. Proteins
        ),
        conditionalPanel(
          condition = "input.comp_type == 'Use SAINTexpress'",
          infoBoxOutput("significantBox_bfdr"), #shows number of Sig. Proteins (BFDR)
          infoBoxOutput("significantBox_bfdr_filtered") #shows number of Sig. Proteins (BFDR)
        ),
        conditionalPanel(
          condition = "input.analyze >= 1",
          box(
            uiOutput("downloadreport"),
            helpText("Long report takes several minutes to generate, check GO plots have completed before downloading report"),
            #helpText("Check Heatmap has no single protein clusters"),
            width = 2, height = 140)
        ),
        box(
          radioButtons("species",
                       "GO species",
                       choices = c("Human", "Mouse", "Arabidopsis", "Xenopus"),
                       selected = "Human"),
          width = 1, height = 140),
      ),
      uiOutput("conversion_status"),
      
      #Main body of
      tabsetPanel(id = "main_body_tab",
                  tabPanel(title = "Introduction",
                           imageOutput("intro_img", fill = FALSE)
                  ),
                  tabPanel("DEP", 
                           fluidRow(
                             column(width = 6,
                                    box(title = "Chosen Comparisons",
                                        DT::dataTableOutput("table_comp"), width = 12),
                                    box(title = "Top Table",
                                        checkboxInput("top_table_opt", "Show all proteins", FALSE),
                                        DT::dataTableOutput("table"), width = 12),
                                    box(title = "Exp Table",
                                        DT::dataTableOutput("exp_des"), width = 12), #exp_design
                             ),
                             column(width = 6,
                                    tabBox(id = "result_plots",
                                           title = "Result Plots", width = 12,
                                           tabPanel(title = "Selected Protein",
                                                    fluidRow(
                                                      box(radioButtons("protein_plot",
                                                                       "Choose plot type",
                                                                       c("contrast", "centered", "normalised", "imputed"),
                                                                       selected = "contrast", inline=T,
                                                      ),
                                                      width = 12),
                                                    ),
                                                    fluidRow(
                                                      plotOutput("selected_plot"),
                                                    ),
                                                    downloadButton('downloadPlot', 'Save plot', class = "btn-primary")),
                                           tabPanel(title = "Selected Protein (Dataset 2)",
                                                    fluidRow(
                                                      box(radioButtons("protein_plot2",
                                                                       "Choose plot type",
                                                                       c("contrast", "centered", "normalised", "imputed"),
                                                                       selected = "contrast", inline=T,
                                                      ),
                                                      width = 12),
                                                    ),
                                                    fluidRow(
                                                      plotOutput("selected_plot2"),
                                                    ),
                                                    downloadButton('downloadPlot2', 'Save plot', class = "btn-primary")),
                                           tabPanel(title = "Abundance Heatmap",
                                                    fluidRow(
                                                      box(numericInput("size",
                                                                       "Number of Most Abundant Proteins",
                                                                       min = 0, max = 1000, value = 100, step = 50),
                                                          width = 6),
                                                      box(radioButtons("abund_choice",
                                                                       "Data Input",
                                                                       c("normalised", "imputed"),
                                                                       selected = "imputed"),
                                                          width = 3),
                                                      box(checkboxInput("signif_choice",
                                                                        "Significant Proteins Only",
                                                                        value = FALSE),
                                                          width = 3)
                                                    ),
                                                    fluidRow(
                                                      uiOutput("plot"),
                                                      downloadButton('downloadHeatmap', 'Save heatmap', class = "btn-primary")),
                                                    fluidRow(
                                                      DT::dataTableOutput("table_abund"), #add table for abundance
                                                      downloadButton('downloadabundancetable', 'Save Abundance Table', class = "btn-primary"),
                                                    ),
                                           ),
                                           tabPanel(title = "Abundance Heatmap (Dataset 2)",
                                                    fluidRow(
                                                      box(numericInput("size2",
                                                                       "Number of Most Abundant Proteins",
                                                                       min = 0, max = 1000, value = 100, step = 50),
                                                          width = 6),
                                                      box(radioButtons("abund_choice2",
                                                                       "Data Input",
                                                                       c("normalised", "imputed"),
                                                                       selected = "imputed"),
                                                          width = 3),
                                                      box(checkboxInput("signif_choice2",
                                                                        "Significant Proteins Only",
                                                                        value = FALSE),
                                                          width = 3)
                                                    ),
                                                    fluidRow(
                                                      uiOutput("plot2"),
                                                      downloadButton('downloadHeatmap2', 'Save heatmap', class = "btn-primary")),
                                                    fluidRow(
                                                      DT::dataTableOutput("table_abund2"), #add table for abundance
                                                      downloadButton('downloadabundancetable2', 'Save Abundance Table', class = "btn-primary"),
                                                    ),
                                           ),
                                           tabPanel(title = "Volcano plot",
                                                    fluidRow(
                                                      box(uiOutput("volcano_cntrst"), width = 6),
                                                      box(numericInput("fontsize",
                                                                       "Font size",
                                                                       min = 0, max = 8, value = 3),
                                                          checkboxInput("check_names",
                                                                        "Selected names only",
                                                                        value = FALSE),
                                                          width = 3),
                                                      box(radioButtons("p_adj",
                                                                       "Label type",
                                                                       choices = c("p-adj sig.","p-Val sig.", "p-adj values")
                                                      ),
                                                          width = 3)
                                                    ),
                                                    fluidRow(
                                                      plotOutput("volcano", height = 700),
                                                      downloadButton('downloadVolcano', 'Save volcano', class = "btn-primary")
                                                    )
                                           ),
                                           tabPanel(title = "Volcano plot (Comp2)",
                                                    fluidRow(
                                                      box(uiOutput("volcano_cntrst_c2"), width = 6),
                                                      box(numericInput("fontsize_c2",
                                                                       "Font size",
                                                                       min = 0, max = 8, value = 3),
                                                          checkboxInput("check_names_c2",
                                                                        "Selected names only",
                                                                        value = FALSE),
                                                          width = 3),
                                                      box(radioButtons("p_adj_c2",
                                                                       "Label type",
                                                                       choices = c("p-adj sig.","p-Val sig.", "p-adj values")
                                                      ),
                                                      width = 3)
                                                    ),
                                                    fluidRow(
                                                      plotOutput("volcano_c2", height = 700),
                                                      downloadButton('downloadVolcano_c2', 'Save volcano', class = "btn-primary")
                                                    )
                                           ),
                                           tabPanel(title = "Significance Counts",
                                                    fluidRow(
                                                      box(numericInput("max_cutoff",
                                                                       "Maximum Cut-off",
                                                                       min = 0, max = 40, value = 20, step = 5),
                                                          width = 6),
                                                    ),
                                                    fluidRow(
                                                      plotOutput("signif_plot", height = 600),
                                                      downloadButton('downloadsignif_plot', 'Save', class = "btn-primary")
                                                    )
                                           ),
                                           tabPanel(title = "Significance Counts (Comp2)",
                                                    fluidRow(
                                                      box(numericInput("max_cutoff2",
                                                                       "Maximum Cut-off",
                                                                       min = 0, max = 40, value = 20, step = 5),
                                                          width = 6),
                                                    ),
                                                    fluidRow(
                                                      plotOutput("signif_plot2", height = 600),
                                                      downloadButton('downloadsignif_plot2', 'Save', class = "btn-primary")
                                                    )
                                           )
                                           
                                    ),
                                    fluidRow(
                                      box(numericInput("pdfwidth",
                                                       "PDF Width",
                                                       min = 3, max = 30, value = 7),
                                          width = 6),
                                      box(numericInput("pdfheight",
                                                       "PDF Height",
                                                       min = 3, max = 30, value = 7),
                                          width = 6)
                                      
                                    ),
                                    tabBox(id = "correlation_plots",
                                           title = "Correlation Plots", width = 12,
                                           tabPanel(title = "PCA",
                                                    fluidRow(
                                                      box(numericInput("pca_num",
                                                                       "Protein number",
                                                                       min = 0, max = 5000, value = 500, step = 100),
                                                          width = 3
                                                      ),
                                                      box(numericInput("pcx_num",
                                                                       "x = PC",
                                                                       min = 1, max = 10, value = 1),
                                                          width = 3
                                                      ),
                                                      box(numericInput("pcy_num",
                                                                       "y = PC",
                                                                       min = 1, max = 10, value = 2),
                                                          width = 3
                                                      ),
                                                      box(checkboxInput("pca_split",
                                                                        "Split by condition_type",
                                                                        value = FALSE),
                                                          width = 3)
                                                    ),
                                                    plotOutput("pca", height = 600),
                                                    downloadButton('downloadpca', 'Save Plot', class = "btn-primary"),
                                                    downloadButton('downloadpca_txt', 'Save Values', class = "btn-primary")
                                           ),
                                           tabPanel(title = "PCA Scree Plot",
                                                    plotOutput("pca_scree", height = 600),
                                                    downloadButton('downloadpca_scree', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Correlation",
                                                    plotOutput("cor", height = 600),
                                                    downloadButton('downloadcor', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Freq. Plot",
                                                    plotOutput("cond", height = 600),
                                                    downloadButton('downloadcond', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "PCA (Comp2)",
                                                    fluidRow(
                                                      box(numericInput("pca_num2",
                                                                       "Protein number",
                                                                       min = 0, max = 5000, value = 500, step = 100),
                                                          width = 3
                                                      ),
                                                      box(numericInput("pcx_num2",
                                                                       "x = PC",
                                                                       min = 1, max = 10, value = 1),
                                                          width = 3
                                                      ),
                                                      box(numericInput("pcy_num2",
                                                                       "y = PC",
                                                                       min = 1, max = 10, value = 2),
                                                          width = 3
                                                      ),
                                                      
                                                      box(checkboxInput("pca_split2",
                                                                        "Split by condition_type",
                                                                        value = FALSE),
                                                          width = 3)
                                                    ),
                                                    plotOutput("pca_c2", height = 600),
                                                    downloadButton('downloadpca_c2', 'Save Plot', class = "btn-primary"),
                                                    downloadButton('downloadpca_c2_txt', 'Save Values', class = "btn-primary")
                                           ),
                                           tabPanel(title = "PCA Scree Plot (Comp2)",
                                                    plotOutput("pca_scree2", height = 600),
                                                    downloadButton('downloadpca_scree2', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Correlation (Comp2)",
                                                    plotOutput("cor_c2", height = 600),
                                                    downloadButton('downloadcor_c2', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Freq. Plot (Comp2)",
                                                    plotOutput("cond_c2", height = 600),
                                                    downloadButton('downloadcond_c2', 'Save', class = "btn-primary")
                                           )
                                           
                                    ),
                                    tabBox(id = "QC_plots",
                                           title = "QC Plots", width = 12,
                                           tabPanel(title = "Protein Numbers",
                                                    plotOutput("numbers", height = 600),
                                                    downloadButton('downloadNumbers', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Protein Numbers (Dataset 2)",
                                                    plotOutput("numbers2", height = 600),
                                                    downloadButton('downloadNumbers2', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Sample coverage",
                                                    plotOutput("coverage", height = 600),
                                                    downloadButton('downloadCoverage', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Normalization",
                                                    fluidRow(
                                                      box(checkboxInput("norm_orient",
                                                                        "Horizontal plots",
                                                                        value = FALSE),
                                                          width = 12)
                                                    ),
                                                    plotOutput("norm", height = 600),
                                                    downloadButton('downloadNorm', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Normalization (Dataset 2)",
                                                    fluidRow(
                                                      box(checkboxInput("norm_orient2",
                                                                        "Horizontal plots",
                                                                        value = FALSE),
                                                          width = 12)
                                                    ),
                                                    plotOutput("norm2", height = 600),
                                                    downloadButton('downloadNorm2', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Missing values - Quant",
                                                    plotOutput("detect", height = 600),
                                                    downloadButton('downloadDetect', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Missing values - Heatmap",
                                                    plotOutput("missval", height = 600),
                                                    downloadButton('downloadMissval', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Imputation",
                                                    plotOutput("imputation", height = 600),
                                                    downloadButton('downloadImputation', 'Save', class = "btn-primary")
                                           ),
                                           tabPanel(title = "Imputation (Dataset 2)",
                                                    plotOutput("imputation2", height = 600),
                                                    downloadButton('downloadImputation2', 'Save', class = "btn-primary")
                                           ),
                                    )
                             )
                           )
                  ),
                  ####UI for Complex Heatmap tab####
                  tabPanel("Cluster Analysis",
                           fluidRow(
                             column(width = 4,
                                    #colour functions
                                    fluidRow(
                                      #Kmeans Cluster box
                                      box(numericInput("k1",
                                                       "Kmeans clusters",
                                                       min = 1, max = 10, value = 4),
                                          checkboxInput("cluster_cols",
                                                        "Cluster columns",
                                                        value = TRUE),
                                          width = 6, height = 120),
                                      box(
                                        selectizeInput("high_col",
                                                       "High Colour",
                                                       choices=c("red","orangered","yellow","green", "darkgreen", "cornflowerblue","blue", "brown", "black"),
                                                       selected = "red"),
                                        width = 3, height = 120),
                                      box(
                                        selectizeInput("low_col",
                                                       "Low Colour",
                                                       choices=c("red","orangered","yellow","green","darkgreen","cornflowerblue","blue", "brown", "black"),
                                                       selected = "blue"),
                                        width = 3, height = 120),
                                    ),
                                    fluidRow(
                                      box(
                                        selectizeInput("palette_choice",
                                                       "Choose palette",
                                                       choices=c("default", "nejm", "Set1","Set2","Set3","Pastel1","Pastel2", "Paired", "Dark2", "Accent", "Reds", "Blues", "Greens", "Greys"),
                                                       selected = "default"),
                                        width = 3, height = 120),
                                      box(numericInput("heatmap_pad",
                                                       "Increase height",
                                                       min = -200, max = 400, value = 0, step = 10),
                                          width = 3, height = 120),
                                      box(numericInput("alpha_input",
                                                       "Alpha",
                                                       min = 0, max = 10, value = 10),
                                          width = 3, height = 120),
                                      box(numericInput("limit_comp",
                                                       "Colour limit (log2)",
                                                       min = 0, max = 16, value = 6),
                                          width = 3, height = 120)
                                    ),
                                    
                                    #Heatmap in variable box
                                    #fluidRow(
                                    #infoBoxOutput("ibox"), #Info box for testing
                                    #),
                                    fluidRow(
                                      uiOutput("plot_comp"),
                                      box(numericInput("pdfwidth_ht",
                                                       "PDF Width",
                                                       min = 3, max = 30, value = 5),
                                          width = 6),
                                      box(numericInput("pdfheight_ht",
                                                       "PDF Height",
                                                       min = 3, max = 30, value = 14),
                                          width = 6),
                                      downloadButton('downloadcomplexheatmap', 'Save ComplexHeatmap', class = "btn-primary")
                                    ),
                                    fluidRow(
                                      plotOutput("swatch", height = 600),
                                      downloadButton('downloadswatch', 'Save Swatch', class = "btn-primary")
                                    )
                             ),
                             
                             #Second Column
                             
                             column(width = 8,
                                    tabBox(id = "comp_heat_tab",
                                           title = "Cluster Analysis Plots", width = 12,
                                           #                                 tabPanel(title = "SAINT Input",
                                           #                                          #DT::dataTableOutput("tablesaint"), width = 12
                                           #                                 ),
                                           tabPanel(title = "SAINT Output",
                                                    DT::dataTableOutput("tablesaint_output1"), width = 12,
                                                    downloadButton('download_saint1', 'Save SAINTexpress Table 1', class = "btn-primary"),
                                                    DT::dataTableOutput("tablesaint_output2"), width = 12,
                                                    downloadButton('download_saint2', 'Save SAINTexpress Table 2', class="btn btn-success"),
                                                    DT::dataTableOutput("tablesaint_output3"), width = 12,
                                                    downloadButton('download_saint3', 'Save SAINTexpress Table 3', class="btn btn-info"),
                                                    DT::dataTableOutput("tablesaint_output4"), width = 12,
                                                    downloadButton('download_saint4', 'Save SAINTexpress Table 4', class ="btn btn-danger")
                                           ),
                                           tabPanel(title = "BFDR Filtering",
                                                    DT::dataTableOutput("tablesaint"), width = 12,
                                                    downloadButton('download_tablesaint', 'Save BFDR Filtered table', class = "btn-primary")
                                                    
                                           ),
                                           tabPanel(title = "Results",
                                                    DT::dataTableOutput("table2"), width = 12
                                           ),
                                           tabPanel(title = "Comparison Heatmaps",
                                                    fluidRow(
                                                      box(title = "Overlap Details",
                                                          column(width=6,
                                                                 fluidRow(
                                                                   plotOutput("venn", height = 400),
                                                                   downloadButton('download_venn', 'Save Venn Diagram', class = "btn-primary"),
                                                                 ),
                                                          ),
                                                          column(width=6,
                                                                 fluidRow(
                                                                   DT::dataTableOutput("venn_table"), width = 12,
                                                                   downloadButton('download_venntable', 'Save Venn table', class = "btn-primary")
                                                                 ),
                                                          ),
                                                          width = 12),
                                                    ),
                                                    fluidRow(
                                                      column(width = 6,
                                                             fluidRow(
                                                               #Kmeans Cluster box
                                                               box(numericInput("k1_c1",
                                                                                "Kmeans clusters",
                                                                                min = 1, max = 10, value = 4),
                                                                   width = 3),
                                                             ),
                                                             #Heatmap in variable box
                                                             fluidRow(
                                                               uiOutput("plot_comp_c1"),
                                                               box(numericInput("pdfwidth_ht_c1",
                                                                                "PDF Width",
                                                                                min = 3, max = 30, value = 5),
                                                                   width = 6),
                                                               box(numericInput("pdfheight_ht_c1",
                                                                                "PDF Height",
                                                                                min = 3, max = 30, value = 14),
                                                                   width = 6),
                                                               downloadButton('downloadcomplexheatmap_c1', 'Save ComplexHeatmap', class = "btn-primary")
                                                             ),
                                                      ),
                                                      #Second Column
                                                      column(width = 6,
                                                             fluidRow(
                                                               #Kmeans Cluster box
                                                               box(numericInput("k1_c2",
                                                                                "Kmeans clusters",
                                                                                min = 1, max = 10, value = 4),
                                                                   width = 3),
                                                             ),
                                                             #Heatmap in variable box
                                                             fluidRow(
                                                               uiOutput("plot_comp_c2"),
                                                               box(numericInput("pdfwidth_ht_c2",
                                                                                "PDF Width",
                                                                                min = 3, max = 30, value = 5),
                                                                   width = 6),
                                                               box(numericInput("pdfheight_ht_c2",
                                                                                "PDF Height",
                                                                                min = 3, max = 30, value = 14),
                                                                   width = 6),
                                                               downloadButton('downloadcomplexheatmap_c2', 'Save ComplexHeatmap', class = "btn-primary")
                                                             ),
                                                      ),
                                                    ),
                                                    
                                           ),
                                           
                                           tabPanel(title = "Cluster",
                                                    fluidRow(
                                                      column(width = 8,
                                                             fluidRow(
                                                               box(numericInput("cluster_no",
                                                                                "Cluster No.",
                                                                                min = 1, max = 10, value = 1),
                                                                   width = 4),
                                                               box(numericInput("cluster_pad",
                                                                                "Increase height",
                                                                                min = -200, max = 400, value = 0, step = 10),
                                                                   width = 4),
                                                             ),
                                                             fluidRow(
                                                               uiOutput("plot_clust"),
                                                               box(numericInput("pdfwidth_htc",
                                                                                "PDF Width",
                                                                                min = 3, max = 30, value = 5),
                                                                   width = 6),
                                                               box(numericInput("pdfheight_htc",
                                                                                "PDF Height",
                                                                                min = 3, max = 30, value = 6),
                                                                   width = 6),
                                                               downloadButton('downloadcomplexheatmapcluster', 'Save ComplexHeatmap', class = "btn-primary")
                                                             )
                                                      ),
                                                      column(width = 4,
                                                             title = "Results",
                                                             DT::dataTableOutput("table3"), #add table for cluster
                                                             downloadButton('downloadclustertable', 'Save Cluster Table', class = "btn-primary"),
                                                             downloadButton('downloadclusterlist', 'Save All Clusters', class = "btn-primary")
                                                      )
                                                    ),
                                           ),
                                           tabPanel(title = "Enhanced Volcano Plots",
                                                    fluidRow(
                                                      column(width = 6,
                                                             fluidRow(
                                                               box(uiOutput("volcano_cntrst_enh1"), width = 8),
                                                               box( radioButtons("p_adj_enh1",
                                                                                 "Label type",
                                                                                 choices = c("p-adj sig.","p-Val sig.", "p-adj values")
                                                               ),
                                                               width = 4)
                                                             ),
                                                             fluidRow(
                                                               plotOutput("enh_volcano_1", height = 600),
                                                               box(checkboxInput("enh_vol_data",
                                                                                 "Significant proteins only",
                                                                                 value = TRUE),
                                                                   width = 3, height = 90),
                                                               box(checkboxInput("enh_vol_names",
                                                                                 "Selected names only",
                                                                                 value = FALSE),
                                                                   width = 3, height = 90),
                                                               box(checkboxInput("enh_vol_connect",
                                                                                 "Draw Arrows",
                                                                                 value = TRUE),
                                                                   width = 3, height = 90),
                                                               box(numericInput("pdfwidth_enh",
                                                                                "PDF Width",
                                                                                min = 3, max = 30, value = 6),
                                                                   width = 3, height = 90)
                                                             ),
                                                             fluidRow(
                                                               downloadButton('downloadVolcano_enh1', 'Save volcano', class = "btn-primary")
                                                             )
                                                             
                                                      ),
                                                      
                                                      #Second Column
                                                      
                                                      column(width = 6,
                                                             fluidRow(
                                                               box(uiOutput("volcano_cntrst_enh2"), width = 8),
                                                               box( radioButtons("p_adj_enh2",
                                                                                 "Label type",
                                                                                 choices = c("p-adj sig.","p-Val sig.", "p-adj values")
                                                               ),
                                                               width = 4)
                                                             ),
                                                             fluidRow(
                                                               plotOutput("enh_volcano_2", height = 600),
                                                               box(numericInput("pdfheight_enh",
                                                                                "PDF Height",
                                                                                min = 3, max = 30, value = 8),
                                                                   width = 3, height = 90),
                                                               box(checkboxInput("enh_vol_data2",
                                                                                 "Significant proteins only",
                                                                                 value = TRUE),
                                                                   width = 3, height = 90),
                                                               box(checkboxInput("enh_vol_names2",
                                                                                 "Selected names only",
                                                                                 value = FALSE),
                                                                   width = 3, height = 90),
                                                               box(checkboxInput("enh_vol_connect2",
                                                                                 "Draw Arrows",
                                                                                 value = TRUE),
                                                                   width = 3, height = 90),
                                                             ),
                                                             fluidRow(
                                                               downloadButton('downloadVolcano_enh2', 'Save volcano', class = "btn-primary")
                                                             )
                                                      ),
                                                    ),
                                                    fluidRow(
                                                      helpText("Select Proteins from the list below. You can search for name or cluster number"),
                                                      DT::dataTableOutput("tablecludf"), width = 12,
                                                    ),
                                           ),
                                           tabPanel(title = "Box Plots",
                                                    fluidRow(
                                                      box(numericInput("cluster_no_box",
                                                                       "Cluster No.",
                                                                       min = 1, max = 10, value = 1),
                                                          width = 3),
                                                      box(radioButtons("pres_box",
                                                                       "Data Type",
                                                                       c("normalised", "imputed"),
                                                                       selected = "imputed"),
                                                          width = 3),
                                                      box(numericInput("pdfwidth_box",
                                                                       "PDF Width",
                                                                       min = 3, max = 30, value = 12),
                                                          width = 3),
                                                      box(numericInput("pdfheight_box",
                                                                       "PDF Height",
                                                                       min = 3, max = 30, value = 12),
                                                          width = 3)
                                                    ),
                                                    fluidRow(
                                                      plotOutput("box_plot", height = 800),
                                                      downloadButton('download_box_plot', 'Save Plot', class = "btn-primary")
                                                    ),
                                           ),
                                           tabPanel(title = "Molecular Function",
                                                    fluidRow(
                                                      box(
                                                      actionButton("analyze_MF", "Analyze Molecular Function", icon = icon("magnifying-glass-chart"), class = "btn-success"),
                                                      helpText("Warning: This can take a couple of minutes"),
                                                      width = 4, height = 90),
                                                      box(numericInput("gomf_level",
                                                                       "GO Level",
                                                                       min = 1, max = 10, value = 5),
                                                          width = 2, height = 90),
                                                      box(numericInput("gomf_node",
                                                                       "Node Size",
                                                                       min = 1, max = 5, value = 2, step = 0.5),
                                                          width = 2, height = 90),
                                                      box(numericInput("pdfwidth_MF",
                                                                       "PDF Width",
                                                                       min = 3, max = 30, value = 8),
                                                          width = 2, height = 90),
                                                      box(numericInput("pdfheight_MF",
                                                                       "PDF Height",
                                                                       min = 3, max = 30, value = 8),
                                                          width = 2, height = 90),
                                                    ),
                                                    fluidRow(
                                                      column(width = 5,
                                                             fluidRow(
                                                               plotOutput("GOMF_plot", height = 800),
                                                               downloadButton('download_MF_plot', 'Save Plot', class = "btn-success")
                                                             )
                                                      ),
                                                      
                                                      #Second Column
                                                      
                                                      column(width = 7,
                                                             fluidRow(
                                                               plotOutput("GOMF_cnet", height = 800),
                                                               downloadButton('download_MF_cnet', 'Save Plot', class = "btn-success")
                                                             )
                                                      ),
                                                    ),
                                                    fluidRow(
                                                      downloadButton('download_MF_all', 'Download MF Data', class = "btn-success")
                                                    ),
                                           ),
                                           tabPanel(title = "Cellular Components",
                                                    fluidRow(
                                                      box(
                                                        actionButton("analyze_CC", "Analyze Cellular Components", icon = icon("magnifying-glass-chart"), class = "btn-info"),
                                                        helpText("Warning: This can take a couple of minutes"),
                                                        width = 4, height = 90),
                                                      box(numericInput("gocc_level",
                                                                       "GO Level",
                                                                       min = 1, max = 10, value = 5),
                                                          width = 2, height = 90),
                                                      box(numericInput("gocc_node",
                                                                       "Node Size",
                                                                       min = 1, max = 5, value = 2, step = 0.5),
                                                          width = 2, height = 90),
                                                      box(numericInput("pdfwidth_CC",
                                                                       "PDF Width",
                                                                       min = 3, max = 30, value = 8),
                                                          width = 2, height = 90),
                                                      box(numericInput("pdfheight_CC",
                                                                       "PDF Height",
                                                                       min = 3, max = 30, value = 8),
                                                          width = 2, height = 90),
                                                    ),
                                                    fluidRow(
                                                      column(width = 5,
                                                             fluidRow(
                                                               plotOutput("GOCC_plot", height = 800),
                                                               downloadButton('download_CC_plot', 'Save Plot', class = "btn-info")
                                                             )
                                                      ),
                                                      
                                                      #Second Column
                                                      
                                                      column(width = 7,
                                                             fluidRow(
                                                               plotOutput("GOCC_cnet", height = 800),
                                                               downloadButton('download_CC_cnet', 'Save Plot', class = "btn-info")
                                                             )
                                                      ),
                                                    ),
                                                    fluidRow(
                                                      downloadButton('download_CC_all', 'Download CC Data', class = "btn-info")
                                                    ),
                                                    
                                           ),
                                           tabPanel(title = "Biological Process",
                                                    fluidRow(
                                                      box(
                                                        actionButton("analyze_BP", "Analyze Biological Process", icon = icon("magnifying-glass-chart"), class = "btn-warning"),
                                                        helpText("Warning: This can take a couple of minutes"),
                                                        width = 4, height = 90),
                                                      box(numericInput("gobp_level",
                                                                       "GO Level",
                                                                       min = 1, max = 10, value = 5),
                                                          width = 2, height = 90),
                                                      box(numericInput("gobp_node",
                                                                       "Node Size",
                                                                       min = 1, max = 5, value = 2, step = 0.5),
                                                          width = 2, height = 90),
                                                      box(numericInput("pdfwidth_BP",
                                                                       "PDF Width",
                                                                       min = 3, max = 30, value = 8),
                                                          width = 2, height = 90),
                                                      box(numericInput("pdfheight_BP",
                                                                       "PDF Height",
                                                                       min = 3, max = 30, value = 8),
                                                          width = 2, height = 90),
                                                      
                                                    ),
                                                    fluidRow(
                                                      column(width = 5,
                                                             fluidRow(
                                                               plotOutput("GOBP_plot", height = 800),
                                                               downloadButton('download_BP_plot', 'Save Plot', class = "btn-warning")
                                                             )
                                                      ),
                                                      
                                                      #Second Column
                                                      
                                                      column(width = 7,
                                                             fluidRow(
                                                               plotOutput("GOBP_cnet", height = 800),
                                                               downloadButton('download_BP_cnet', 'Save Plot', class = "btn-warning")
                                                             )
                                                      ),
                                                    ),
                                                    fluidRow(
                                                      downloadButton('download_BP_all', 'Download BP Data', class = "btn-warning")
                                                    ),
                                           ),
                                           tabPanel(title = "BioGrid Interactions",
                                                    fluidRow(
                                                      box(
                                                        actionButton("analyze_BG", "Analyze BioGrid Interactions", icon = icon("magnifying-glass-chart"), class = "btn-primary"),
                                                        helpText("Warning: This can take a couple of minutes"),
                                                        helpText("Max limit of 160 significant proteins"),
                                                        width = 4, height = 110),
                                                      box(numericInput("BG_level",
                                                                       "Interaction Count <=",
                                                                       min = 1, max = 10, value = 2),
                                                          width = 2, height = 110),
                                                      box(numericInput("BG_node",
                                                                       "Node Size",
                                                                       min = 1, max = 5, value = 1, step = 0.5),
                                                          width = 2, height = 110),
                                                      box(checkboxInput("use_human",
                                                                        "Use Human Interaction Network Data ",
                                                                        value = FALSE),
                                                          helpText("This will convert Gene names to uppercase"),
                                                          width = 4, height = 110),
                                                      
                                                    ),
                                                    fluidRow(
                                                      plotOutput("BG_cnet", height = 800),
                                                    ),
                                                    fluidRow(
                                                      box(downloadButton('download_BG_cnet', 'Save Plot', class = "btn-primary"),
                                                          width = 2, height = 80),
                                                      box(downloadButton('download_BG_data', 'Download BioGrid Search', class = "btn-primary"),
                                                          width = 3, height = 80),
                                                      box(numericInput("pdfwidth_BG",
                                                                       "PDF Width",
                                                                       min = 3, max = 30, value = 10),
                                                          width = 3, height = 80),
                                                      box(numericInput("pdfheight_BG",
                                                                       "PDF Height",
                                                                       min = 3, max = 30, value = 10),
                                                          width = 3, height = 80),
                                                    ),
                                                    
                                           )
                                    )
                             )
                           )
                  ),
                  #Third Tab (Instructions)
                  tabPanel("Instructions",
                           mainPanel(
                             style = "background-color: #ffffff;",
                             width = 12,
                             htmlOutput("help")
                           )
                  )
      )
    )
  )
)
#End of UI

####Server####

server <- shinyServer(function(input, output, session) {
  #options(shiny.usecairo = TRUE)  # Force Cairo rendering in Shiny
  options(shiny.useragg = TRUE)
  options(shiny.graphics.device = ragg::agg_png)
  conversion_messages <- reactiveVal(NULL)
  
  
  
  # Allow reconnect if the session is interrupted
  session$allowReconnect("force")
  
  options(shiny.maxRequestSize=200*1024^2)
  
  # Reactive trigger for programmatic start
  trigger_GO_MF <- reactiveVal(FALSE)
  trigger_GO_CC <- reactiveVal(FALSE)
  trigger_GO_BP <- reactiveVal(FALSE)
  trigger_BG <- reactiveVal(FALSE)
  
  
  
  # Initialise the function variable as NULL
  GO_MF_report_input <- reactiveVal(NULL)  # Store async result
  GO_CC_report_input <- reactiveVal(NULL)  # Store async result
  GO_BP_report_input <- reactiveVal(NULL)  # Store async result
  BG_report_input <- reactiveVal(NULL)  # Store async result
  
  
  
  #### UI functions --------------------------------------------------------
  #### Side bar menu items ####
  #CRUK Logo
  output$logo_img <- renderImage({
    
    list(src = "www/logo.jpg",
         width = "100%",
         height = 100)
    
  }, deleteFile = F)
  
  output$logo_wccr <- renderImage({
    
    list(src = "www/CELL MATRIX RGB2.jpg",
         width = "100%",
         height = 100)
    
  }, deleteFile = F)
  
  output$intro_img <- renderImage({
    width_img  <- session$output_plot2_width
    height_img <- session$output_plot2_height
    list(src = "www/Introduction.svg",
         width = width_img,
         height = height_img)
    
  }, deleteFile = F)
  
  output$conversion_status <- renderUI({
    msgs <- conversion_messages()
    
    # Filter specific error-related messages
    #error_msgs <- msgs[grepl("None of the keys entered are valid|No enrichment found in any of gene cluster", msgs)]
    error_msgs <- msgs[grepl("None of the keys entered are valid", msgs)]
    
    if (length(error_msgs) > 0) {
      print("Error in asigning GO terms, please check the correct G0 species is selected")
    }
  })
  
  
  
  # getPage<-function() {
  #    return(includeHTML("www/index.html"))
  #  }
  #  output$inc<-renderUI({getPage()})
  
  addResourcePath("src", "./www")
  output$help <- renderUI({
    tags$iframe(
      seamless="seamless",
      src="index.html",
      width="100%",
      height="800px",
      style="border: none")
  })
  
  #Side menu for picking other data catergories
  output$columns <- renderMenu({
    if (!is.null(data_other())) {
      col_names<-c("",colnames(data_other()))
      menuItem("Columns",
               startExpanded = TRUE,
               selectizeInput("name_other",
                              "Name column",
                              choices=col_names
               ),
               selectizeInput("id_other",
                              "ID column",
                              choices=col_names
               ),
               selectizeInput("description_other",
                              "Description",
                              choices=col_names
               ),
               selectizeInput("aa_other",
                              "aa length",
                              choices=col_names
               ),
               selectizeInput("int_col_start",
                              "1st Intensity Column",
                              choices=col_names
               ),
               selectizeInput("int_col_last",
                              "Last Intensity Column",
                              choices=col_names
               ),
               radioButtons("other_species",
                            "Choose Species",
                            choices = c("Human", "Mouse"),
                            selected = "Human"),
               
               helpText("Select matching column headings from the"),
               helpText("data table, none selected columns will be"),
               helpText("filled from Uniprot data of the chosen species")
               
      )}
  })
  
  output$columns2 <- renderMenu({
    if (!is.null(data_other2())) {
      col_names<-c("",colnames(data_other2()))
      menuItem("Columns",
               startExpanded = TRUE,
               selectizeInput("name_other2",
                              "Name column",
                              choices=col_names
               ),
               selectizeInput("id_other2",
                              "ID column",
                              choices=col_names
               ),
               selectizeInput("description_other2",
                              "Description",
                              choices=col_names
               ),
               selectizeInput("aa_other2",
                              "aa length",
                              choices=col_names
               ),
               selectizeInput("int_col_start2",
                              "1st Intensity Column",
                              choices=col_names
               ),
               selectizeInput("int_col_last2",
                              "Last Intensity Column",
                              choices=col_names
               ),
               radioButtons("other_species2",
                            "Choose Species",
                            choices = c("Human", "Mouse"),
                            selected = "Human"),
               
               helpText("Select matching column headings from the"),
               helpText("data table, none selected columns will be"),
               helpText("filled from Uniprot data of the chosen species")
               
      )}
  })
  
  
  
  output$comparisons <- renderMenu({
    menuItem("Comparisons",
             if (!is.null(expdesign())) {
               combinations<-combn(unique(expdesign()$condition),2)
               comp1_f<-paste0(combinations[2,],"_vs_",combinations[1,])
               comp1_r<-paste0(combinations[1,],"_vs_",combinations[2,])
               comp1_all<-c(comp1_f,comp1_r)
               selectizeInput("comp1f",
                              "1: Select 1st round comparisons",
                              choices=make.names(comp1_all),
                              multiple = TRUE)
             },
             
             if (!is.null(expdesign()) && input$comp_type == "2 Rounds of comparisons") {
               combinations<-combn(unique(expdesign()$condition),2)
               comp2_f<-paste0(combinations[2,],"_vs_",combinations[1,])
               comp2_r<-paste0(combinations[1,],"_vs_",combinations[2,])
               comp2_all<-c(comp2_f,comp2_r)
               selectizeInput("comp2f",
                              "2: Select 2nd round comparisons",
                              choices=make.names(comp2_all),
                              multiple = TRUE)
             },
             
             if (!is.null(expdesign2()) && input$data_mode == TRUE) {
               combinations<-combn(unique(expdesign2()$condition),2)
               comp2_f<-paste0(combinations[2,],"_vs_",combinations[1,])
               comp2_r<-paste0(combinations[1,],"_vs_",combinations[2,])
               comp2_all<-c(comp2_f,comp2_r)
               selectizeInput("comp2f",
                              "2: Select 2nd Dataset comparisons",
                              choices=make.names(comp2_all),
                              multiple = TRUE)
             },
             
             
             helpText("Load Comparisons file"),
             fileInput('comp_file',
                       'comparisons.txt',
                       accept=c('text/csv',
                                'text/comma-separated-values,text/plain',
                                '.csv')),
             downloadButton('downloadcomp_template', 'Save Comp. template', class="btn btn-primary")
    )
  })
  
  output$imputation_options <- renderMenu({
    menuItem("Imputation options",
             radioButtons("imputation",
                          "Imputation type (Dataset 1)",
                          choices = c("man", "NIPALS", "knn", "QRILC", "MinDet", "MinProb",  "min", "zero", "bpca"),
                          #choices = c("man", MsCoreUtils::imputeMethods())[1:12],
                          selected = "man"),
             helpText("Manual (man) Imputation settings"),
             numericInput("shift",
                          "Left-shifted Gaussian distribution",
                          min = 1.1, max = 2.9, value = 1.8, step = 0.1),
             numericInput("width",
                          "Width",
                          min = 0, max = 1, value = 0.3, step = 0.1),
             p(a("Detailed information link",
                 href = "https://rformassspectrometry.github.io/MsCoreUtils/reference/imputation.html",
                 target="https://rformassspectrometry.github.io/MsCoreUtils/reference/imputation.html")),
             checkboxInput("impute_ave", "Average Imputation of 3 different random seeds (1st Dataset)", FALSE),
             if (input$data_mode == TRUE) {
               radioButtons("imputation2",
                            "Imputation type (Dataset 2)",
                            choices = c("man", "NIPALS", "knn", "QRILC", "MinDet", "MinProb",  "min", "zero", "bpca"),
                            #choices = c("man", MsCoreUtils::imputeMethods())[1:12],
                            selected = "man")
             },
             if (input$data_mode == TRUE) {
               checkboxInput("impute_ave2", "Average Imputation of 3 different random seeds (2nd Dataset)", FALSE)
             }
    )
  })
  
  output$normalisations <- renderMenu({
    menuItem("Normalisation options",
             radioButtons("normalisation_choice",
                          "Normalisation Method",
                          choices = c("vsn", "median center", "quantile", "none"),
                          selected = "vsn"),
             if (input$data_mode == TRUE) {
               radioButtons("normalisation_choice2",
                            "Normalisation Method Dataset 2)",
                            choices = c("vsn", "median center", "quantile", "none"),
                            selected = "vsn")
             }
             
    )
  })
  
  output$padj_options <- renderMenu({
    menuItem("p.adj method options",
             radioButtons("padjust_method",
                          "p.adj (FDR) Methods",#should be one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
                          choices = c("default (DEP)" = "default",
                                      "Benjamini-Hochberg" = "BH",
                                      "Benjamini-Yekutieli" = "BY",
                                      "Bonferroni" = "bonferroni",
                                      "Hochberg" = "hochberg",
                                      #"Holm" = "holm", # Hochberg and Holm give same result
                                      "Hommel" = "hommel"),
                          selected = "default"),
             if (input$data_mode == TRUE) {
               radioButtons("padjust_method2",
                            "p.adj (FDR) Methods (Dataset 2",#should be one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
                            choices = c("default (DEP)" = "default",
                                        "Benjamini-Hochberg" = "BH",
                                        "Benjamini-Yekutieli" = "BY",
                                        "Bonferroni" = "bonferroni",
                                        "Hochberg" = "hochberg",
                                        #"Holm" = "holm",
                                        "Hommel" = "hommel"),
                            selected = "default")
             }
             
    )
  })
  
  
  #Hide/show Introduction Tab
  observeEvent(input$analyze, {
    if (input$analyze > 0) {
      hideTab(inputId =  "main_body_tab", target = "Introduction")
    }else{showTab(inputId =  "main_body_tab", target = "Introduction")}
  })
  
  #Hide/show comparison heatmaps
  observeEvent(ignoreInit = FALSE, list(input$comp_type, input$data_mode), {
    if (input$comp_type != "1 Round of comparisons" | input$data_mode == TRUE) {
      showTab(inputId =  "comp_heat_tab", target = "Comparison Heatmaps")
    }else{hideTab(inputId =  "comp_heat_tab", target = "Comparison Heatmaps")}
    if (input$comp_type == "1 Round of comparisons"|input$comp_type == "2 Rounds of comparisons") {
      hideTab(inputId =  "comp_heat_tab", target = "BFDR Filtering")
      #hideTab(inputId =  "comp_heat_tab", target = "SAINT Input")
      hideTab(inputId =  "comp_heat_tab", target = "SAINT Output")
    }else{showTab(inputId =  "comp_heat_tab", target = "BFDR Filtering")
      #showTab(inputId =  "comp_heat_tab", target = "SAINT Input")
      showTab(inputId =  "comp_heat_tab", target = "SAINT Output")
    }
  })
  
  #Hide/show Results tab Plots
  observeEvent(ignoreInit = FALSE, list(input$comp_type, input$data_mode), {
    if (input$comp_type == "2 Rounds of comparisons" | input$data_mode == TRUE) {
      showTab(inputId =  "result_plots", target = "Volcano plot (Comp2)")
    }else{hideTab(inputId =  "result_plots", target = "Volcano plot (Comp2)")}
    if (input$data_mode == TRUE) {
      showTab(inputId =  "result_plots", target = "Abundance Heatmap (Dataset 2)")
    }else{hideTab(inputId =  "result_plots", target = "Abundance Heatmap (Dataset 2)")}
    if (input$comp_type == "2 Rounds of comparisons" | input$data_mode == TRUE) {
      showTab(inputId =  "result_plots", target = "Significance Counts (Comp2)")
    }else{hideTab(inputId =  "result_plots", target = "Significance Counts (Comp2)")}
    if (input$data_mode == TRUE) {
      showTab(inputId =  "result_plots", target = "Selected Protein (Dataset 2)")
    }else{hideTab(inputId =  "result_plots", target = "Selected Protein (Dataset 2)")}
  })
  
  #Hide/show Correlation plots/PCA Plots/Frequency Plots
  observeEvent(ignoreInit = FALSE, list(input$comp_type, input$data_mode), {
    if (input$data_mode == TRUE) {
      showTab(inputId =  "correlation_plots", target = "PCA (Comp2)")
    }else{hideTab(inputId =  "correlation_plots", target = "PCA (Comp2)")}
    if (input$data_mode == TRUE) {
      showTab(inputId =  "correlation_plots", target = "PCA Scree Plot (Comp2)")
    }else{hideTab(inputId =  "correlation_plots", target = "PCA Scree Plot (Comp2)")}
    if (input$comp_type == "2 Rounds of comparisons" | input$data_mode == TRUE) {
      showTab(inputId =  "correlation_plots", target = "Correlation (Comp2)")
    }else{hideTab(inputId =  "correlation_plots", target = "Correlation (Comp2)")}
    if (input$comp_type == "2 Rounds of comparisons" | input$data_mode == TRUE) {
      showTab(inputId =  "correlation_plots", target = "Freq. Plot (Comp2)")
    }else{hideTab(inputId =  "correlation_plots", target = "Freq. Plot (Comp2)")}
  })
  
  #Hide/show Numbers/Normalisation/Imputation Plots
  observeEvent(ignoreInit = FALSE, list(input$comp_type, input$data_mode), {
    if (input$data_mode == TRUE) {
      showTab(inputId =  "QC_plots", target = "Protein Numbers (Dataset 2)")
    }else{hideTab(inputId =  "QC_plots", target = "Protein Numbers (Dataset 2)")}
    if (input$data_mode == TRUE) {
      showTab(inputId =  "QC_plots", target = "Normalization (Dataset 2)")
    }else{hideTab(inputId =  "QC_plots", target = "Normalization (Dataset 2)")}
    if (input$data_mode == TRUE) {
      showTab(inputId =  "QC_plots", target = "Imputation (Dataset 2)")
    }else{hideTab(inputId =  "QC_plots", target = "Imputation (Dataset 2)")}
  })
  
  
  #### Reactive functions  --------------------------------------------------
  set.seed(123) #random seed to give consistent kmeans clustering
  
  #Output information for testing
  #output$ibox <- renderInfoBox({
  #infoBox(
  #"Title",
  #Check for single protein clusters
  #validate(need(try(cludf_check<-cluster_heatmap_input()[[19]] != "try-error"), "One of your Clusters contains a single protein, 
  #  please alter the Kmeans number"),
  #         print("All checks out")
  #)
  #icon = icon("credit-card")
  #)
  #})
  
  data_other <- reactive({
    inFile <- input$file_other
    if (is.null(inFile))
      return(NULL)
    line<-paste(Sys.time(),"Other")
    write(line,file="count.txt",append=TRUE)
    data<-read.csv(inFile$datapath, header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE) %>%
      mutate(id = row_number())
    return(data)
  })
  
  data_other2 <- reactive({
    inFile <- input$file_other2
    if (is.null(inFile))
      return(NULL)
    line<-paste(Sys.time(),"Other (Dataset 2")
    write(line,file="count.txt",append=TRUE)
    data<-read.csv(inFile$datapath, header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE) %>%
      mutate(id = row_number())
    return(data)
  })
  
  
  data_other_choice <- reactive({
    if (is.null(data_other()))
      return(NULL)
    tempdata1<-data_other()
    tempdata2<-tempdata1 %>% dplyr::rename(Gene.names=input$name_other)
    tempdata2<- tempdata2[tempdata2$Gene.names != "", ]
    tempdata2$Gene.names  <- sub('\\;.*', '', tempdata2$Gene.names)
    
    if(input$other_species == "Human"){uniprot<-read.csv("./uniprot/human_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
    if(input$other_species == "Mouse"){uniprot<-read.csv("./uniprot/mouse_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
    tempdata3<-merge(tempdata2, uniprot, by.x = "Gene.names", by.y = "uniprot_gene", all.x = TRUE, all.y = FALSE)
    
    if (input$id_other!= ""){tempdata4<-tempdata3 %>% dplyr::rename(Protein.IDs=input$id_other)}else{tempdata4<-tempdata3 %>% dplyr::rename(Protein.IDs=uniprot_entry)}
    if (input$description_other!= ""){tempdata5<-tempdata4 %>% dplyr::rename(Protein.names=input$description_other)}else{tempdata5<-tempdata4 %>% dplyr::rename(Protein.names=uniprot_description)}
    if (input$aa_other!= ""){tempdata6<-tempdata5 %>% dplyr::rename(Sequence.length=input$aa_other)}else{tempdata6<-tempdata5 %>% dplyr::rename(Sequence.length=uniprot_length)}
    
    col_num_start <- which(names(tempdata6) == input$int_col_start)
    col_num_end <- which(names(tempdata6) == input$int_col_last)
    lfq_columns <- col_num_start:col_num_end
    #if (input$id_other)
    data1<-tempdata6[, c(grep("Gene.names", colnames(tempdata6)), grep("Protein.IDs", colnames(tempdata6)), grep("Protein.names", colnames(tempdata6)), grep("\\Sequence.length\\b", colnames(tempdata6)),lfq_columns)]
    
    return(data1)
  })
  
  data_other_choice2 <- reactive({
    if (is.null(data_other2()))
      return(NULL)
    tempdata1<-data_other2()
    tempdata2<-tempdata1 %>% dplyr::rename(Gene.names=input$name_other2)
    tempdata2<- tempdata2[tempdata2$Gene.names != "", ]
    tempdata2$Gene.names  <- sub('\\;.*', '', tempdata2$Gene.names)
    
    if(input$other_species2 == "Human"){uniprot<-read.csv("./uniprot/human_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
    if(input$other_species2 == "Mouse"){uniprot<-read.csv("./uniprot/mouse_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
    tempdata3<-merge(tempdata2, uniprot, by.x = "Gene.names", by.y = "uniprot_gene", all.x = TRUE, all.y = FALSE)
    
    if (input$id_other2!= ""){tempdata4<-tempdata3 %>% dplyr::rename(Protein.IDs=input$id_other2)}else{tempdata4<-tempdata3 %>% dplyr::rename(Protein.IDs=uniprot_entry)}
    if (input$description_other2!= ""){tempdata5<-tempdata4 %>% dplyr::rename(Protein.names=input$description_other2)}else{tempdata5<-tempdata4 %>% dplyr::rename(Protein.names=uniprot_description)}
    if (input$aa_other2!= ""){tempdata6<-tempdata5 %>% dplyr::rename(Sequence.length=input$aa_other2)}else{tempdata6<-tempdata5 %>% dplyr::rename(Sequence.length=uniprot_length)}
    
    col_num_start <- which(names(tempdata6) == input$int_col_start2)
    col_num_end <- which(names(tempdata6) == input$int_col_last2)
    lfq_columns <- col_num_start:col_num_end
    #if (input$id_other)
    data1<-tempdata6[, c(grep("Gene.names", colnames(tempdata6)), grep("Protein.IDs", colnames(tempdata6)), grep("Protein.names", colnames(tempdata6)), grep("\\Sequence.length\\b", colnames(tempdata6)),lfq_columns)]
    
    return(data1)
  })
  
  
  data <- reactive({
    if (input$data_type == "MaxQuant"){
      inFile <- input$file1
      if (is.null(inFile))
        return(NULL)
      line<-paste(Sys.time(),"MaxQuant")
      write(line,file="count.txt",append=TRUE)
      maxdata<-read.csv(inFile$datapath, header = TRUE,
                        sep = "\t", stringsAsFactors = FALSE) %>%
        mutate(id = row_number())
      data1<-maxdata[, c(grep("Protein.IDs", colnames(maxdata)),grep("Protein.names", colnames(maxdata)),grep("Gene.names", colnames(maxdata)),  grep("LFQ.",colnames(maxdata)), grep ("\\Sequence.length\\b", colnames(maxdata)), grep("Reverse",colnames(maxdata)), grep("Potential.contaminant",colnames(maxdata)))]
    }
    
    else if (input$data_type == "Proteome Discoverer"){
      inFile <- input$filePD
      if (is.null(inFile))
        return(NULL)
      line<-paste(Sys.time(),"Proteome_Discover")
      write(line,file="count.txt",append=TRUE)
      pd<-read_excel(path = inFile$datapath, sheet = 1, .name_repair = "universal") %>%
        mutate(id = row_number())
      data<-pd %>% 
        dplyr::rename(
          Gene.names=Gene.Symbol,
          Protein.IDs=Accession,
          Protein.names=Description,
          Sequence.length=..AAs
        ) 
      #      data1<-data[, c(1,2,3,  grep("Abundance..*Sample..", colnames(data)))]
      data1<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("Protein.names", colnames(data)),  grep("Abundance..*Sample..", colnames(data)), grep("Abundance..*Control..", colnames(data)), grep ("Sequence.length", colnames(data)))]
    }
    
    else if (input$data_type == "FragPipe"){
      if (input$exp == "LFQ") {
        inFile <- input$lfq_expr
      } else if (input$exp == "TMT") {
        inFile <- input$tmt_expr
      } else if (input$exp == "DIA") {
        inFile <- input$dia_expr
      }
      if(is.null(inFile))
        return(NULL)
      if (input$exp == "LFQ") {
        line<-paste(Sys.time(),"FragPipe_DDA")
        write(line,file="count.txt",append=TRUE)
        temp_data <- read.csv(inFile$datapath,
                              header = TRUE,
                              fill= TRUE, # to fill any missing data
                              sep = "\t",
                              quote = "",
                              comment.char = "",
                              blank.lines.skip = F,
                              check.names = T)
        colnames(temp_data) <- make.unique(colnames(temp_data), "_")
        if ("Protein" %in% colnames(temp_data)) {
          temp_data <- temp_data[!grepl("contam", temp_data$Protein),]
        }
        data<-temp_data %>% 
          dplyr::rename(
            Gene.names=Gene,
            Protein.IDs=Protein.ID,
            Protein.names=Description,
            Sequence.length=Protein.Length
          )
        
        if (input$lfq_type == "Intensity") {
          lfq_columns <- setdiff(grep("Intensity", colnames(data)),
                                 grep("MaxLFQ", colnames(data)))
          lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data)))
          lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data)))
          data1<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("Protein.names", colnames(data)), lfq_columns,  grep ("Sequence.length", colnames(data)))]
          #        data1<-data1 %>% 
          #          rename_with(~ gsub(".Intensity", "", .x, fixed = TRUE))
          
        } else if (input$lfq_type == "MaxLFQ") {
          lfq_columns<-grep("MaxLFQ", colnames(data))
          if (length(lfq_columns) == 0) {
            stop(safeError("No MaxLFQ column available. Please make sure your files have MaxLFQ intensity columns."))
          }
          data1<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("Protein.names", colnames(data)), lfq_columns,  grep ("Sequence.length", colnames(data)))]
        }
      }else if(input$exp == "DIA"){
        line<-paste(Sys.time(),"FragPipe_DIA")
        write(line,file="count.txt",append=TRUE)
        temp_data <- read.table(inFile$datapath,
                                header = T,
                                sep="\t",
                                check.names = T,
                                stringsAsFactors = FALSE,
                                fill= TRUE, # to fill any missing data
                                quote = "",
                                comment.char = "",
                                blank.lines.skip = F,
        )
        temp_data<-temp_data %>% dplyr::filter(!grepl('contam_sp', Protein.Group))
        temp_data<-temp_data %>% dplyr::filter(!grepl('rev_sp', Protein.Group)) 
        
        data<-temp_data %>% 
          dplyr::rename(
            Gene.names=Genes,
            Protein.IDs=Protein.Group,
            Protein.names=First.Protein.Description,
            #  Sequence.length=Protein.Length
          )
        data$Protein.IDs  <- sub('\\;.*', '', data$Protein.IDs)
        lfq_columns<- 6:ncol(data)
        if(input$dia_species == "Human"){uniprot<-read.csv("./uniprot/human_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
        if(input$dia_species == "Mouse"){uniprot<-read.csv("./uniprot/mouse_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
        data_m<-merge(data, uniprot, by.x = "Protein.IDs", by.y = "uniprot_entry", all.x = TRUE, all.y = FALSE)
        data<-data_m %>% 
          dplyr::rename(
            Sequence.length=uniprot_length
          )
        data_2<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("uniprot_description", colnames(data)), lfq_columns,  grep ("Sequence.length", colnames(data)))]
        data1<-data_2 %>% 
          dplyr::rename(
            Protein.names=uniprot_description)
        data$Protein.names  <- gsub("\\(.*?\\)",".", data$Protein.names)
        
        
        
      } # end of DIA
    }
    
    else if (input$data_type == "Other (Experimental)"){
      if (is.null(data_other()))
        return(NULL)
      data1<-data_other_choice()
    }
    
    return(data1)
  })
  
  data2 <- reactive({
    if (input$data_type2 == "MaxQuant"){
      inFile <- input$fileMax2
      if (is.null(inFile))
        return(NULL)
      line<-paste(Sys.time(),"MaxQuant_data2")
      write(line,file="count.txt",append=TRUE)
      maxdata<-read.csv(inFile$datapath, header = TRUE,
                        sep = "\t", stringsAsFactors = FALSE) %>%
        mutate(id = row_number())
      data1<-maxdata[, c(grep("Protein.IDs", colnames(maxdata)),grep("Protein.names", colnames(maxdata)),grep("Gene.names", colnames(maxdata)),  grep("LFQ.",colnames(maxdata)), grep ("\\Sequence.length\\b", colnames(maxdata)), grep("Reverse",colnames(maxdata)), grep("Potential.contaminant",colnames(maxdata)))]
    }
    
    else if (input$data_type2 == "Proteome Discoverer"){
      inFile <- input$filePD2
      if (is.null(inFile))
        return(NULL)
      line<-paste(Sys.time(),"Proteome_Discover_data2")
      write(line,file="count.txt",append=TRUE)
      pd<-read_excel(path = inFile$datapath, sheet = 1, .name_repair = "universal") %>%
        mutate(id = row_number())
      data<-pd %>% 
        dplyr::rename(
          Gene.names=Gene.Symbol,
          Protein.IDs=Accession,
          Protein.names=Description,
          Sequence.length=..AAs
        ) 
      #      data1<-data[, c(1,2,3,  grep("Abundance..*Sample..", colnames(data)))]
      data1<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("Protein.names", colnames(data)),  grep("Abundance..*Sample..", colnames(data)), grep("Abundance..*Control..", colnames(data)), grep ("Sequence.length", colnames(data)))]
    }
    
    else if (input$data_type2 == "FragPipe"){
      if (input$exp2 == "LFQ") {
        inFile <- input$lfq_expr2
      } else if (input$exp2 == "TMT") {
        inFile <- input$tmt_expr2
      } else if (input$exp2 == "DIA") {
        inFile <- input$dia_expr2
      }
      if(is.null(inFile))
        return(NULL)
      if (input$exp2 == "LFQ") {
        line<-paste(Sys.time(),"FragPipe_DDA_data2")
        write(line,file="count.txt",append=TRUE)
        temp_data <- read.csv(inFile$datapath,
                              header = TRUE,
                              fill= TRUE, # to fill any missing data
                              sep = "\t",
                              quote = "",
                              comment.char = "",
                              blank.lines.skip = F,
                              check.names = T)
        colnames(temp_data) <- make.unique(colnames(temp_data), "_")
        if ("Protein" %in% colnames(temp_data)) {
          temp_data <- temp_data[!grepl("contam", temp_data$Protein),]
        }
        data<-temp_data %>% 
          dplyr::rename(
            Gene.names=Gene,
            Protein.IDs=Protein.ID,
            Protein.names=Description,
            Sequence.length=Protein.Length
          )
        
        if (input$lfq_type2 == "Intensity") {
          lfq_columns <- setdiff(grep("Intensity", colnames(data)),
                                 grep("MaxLFQ", colnames(data)))
          lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data)))
          lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data)))
          data1<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("Protein.names", colnames(data)), lfq_columns,  grep ("Sequence.length", colnames(data)))]
          #        data1<-data1 %>% 
          #          rename_with(~ gsub(".Intensity", "", .x, fixed = TRUE))
          
        } else if (input$lfq_type2 == "MaxLFQ") {
          lfq_columns<-grep("MaxLFQ", colnames(data))
          if (length(lfq_columns) == 0) {
            stop(safeError("No MaxLFQ column available. Please make sure your files have MaxLFQ intensity columns."))
          }
          data1<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("Protein.names", colnames(data)), lfq_columns,  grep ("Sequence.length", colnames(data)))]
        }
      }else if(input$exp2 == "DIA"){
        line<-paste(Sys.time(),"FragPipe_DIA_data2")
        write(line,file="count.txt",append=TRUE)
        temp_data <- read.table(inFile$datapath,
                                header = T,
                                sep="\t",
                                check.names = T,
                                stringsAsFactors = FALSE,
                                fill= TRUE, # to fill any missing data
                                quote = "",
                                comment.char = "",
                                blank.lines.skip = F,
        )
        temp_data<-temp_data %>% dplyr::filter(!grepl('contam_sp', Protein.Group))
        temp_data<-temp_data %>% dplyr::filter(!grepl('rev_sp', Protein.Group)) 
        
        data<-temp_data %>% 
          dplyr::rename(
            Gene.names=Genes,
            Protein.IDs=Protein.Group,
            Protein.names=First.Protein.Description,
            #  Sequence.length=Protein.Length
          )
        data$Protein.IDs  <- sub('\\;.*', '', data$Protein.IDs)
        lfq_columns<- 6:ncol(data)
        if(input$dia_species2 == "Human"){uniprot<-read.csv("./uniprot/human_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
        if(input$dia_species2 == "Mouse"){uniprot<-read.csv("./uniprot/mouse_uniprot.tsv", head=T, sep="\t", stringsAsFactors = FALSE)}
        data_m<-merge(data, uniprot, by.x = "Protein.IDs", by.y = "uniprot_entry", all.x = TRUE, all.y = FALSE)
        data<-data_m %>% 
          dplyr::rename(
            Sequence.length=uniprot_length
          )
        data_2<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("uniprot_description", colnames(data)), lfq_columns,  grep ("Sequence.length", colnames(data)))]
        data1<-data_2 %>% 
          dplyr::rename(
            Protein.names=uniprot_description)
        data$Protein.names  <- gsub("\\(.*?\\)",".", data$Protein.names)
        
      } # end of DIA
    }
    
    else if (input$data_type2 == "Other (Experimental)"){
      if (is.null(data_other2()))
        return(NULL)
      data1<-data_other_choice2()
    }
    
    return(data1)
  })
  
  
  
  #make experimental template
  exp_template<- reactive({
    if (input$data_type == "MaxQuant"){
      LFQ_columns<- grep("LFQ.",colnames(data())) 
      label<-colnames(data()[LFQ_columns])
      label <- gsub('LFQ.intensity.', '', label)
      condition<-""
      substrRight <- function(x, n){
        substr(x, nchar(x)-n+1, nchar(x))
      }
      replicate<-substrRight(label, 1)
      condition_type<-"1"
      bait_id<-"T"
      exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
    }
    
    else if (input$data_type == "Proteome Discoverer"){
      Sample_columns<- grep("Abundance..*Sample..", colnames(data()))
      Control_columns<- grep("Abundance..*Control..", colnames(data()))
      LFQ_columns<-c(Sample_columns, Control_columns)
      label<-colnames(data()[LFQ_columns])
      condition<-""
      replicate<-""
      condition_type<-"1"
      bait_id<-"T"
      exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
    }
    
    else if (input$data_type == "FragPipe"){
      if(input$exp == "LFQ")
      {
        LFQ_columns<- grep(".Intensity", colnames(data()))
        label<-colnames(data()[LFQ_columns])
        label <- gsub('.MaxLFQ.Intensity', '', label)
        label <- gsub('.Intensity', '', label)
        condition<-""
        substrRight <- function(x, n){
          substr(x, nchar(x)-n+1, nchar(x))
        }
        replicate<-substrRight(label, 1)
        condition_type<-"1"
        bait_id<-"T"
        exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
      }else if (input$exp == "DIA")
      {
        end<-ncol(data())
        end_1<-end-1
        LFQ_columns<- 4:end_1
        label<-colnames(data()[LFQ_columns])
        #label <- gsub('.MaxLFQ.Intensity', '', label)
        #label <- gsub('.Intensity', '', label)
        condition<-""
        #substrRight <- function(x, n){
        #  substr(x, nchar(x)-n+1, nchar(x))
        #}
        replicate<-""
        condition_type<-"1"
        bait_id<-"T"
        exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
        
      }
    }else if (input$data_type == "Other (Experimental)")
    {
      end<-ncol(data())
      LFQ_columns<- 5:end
      label<-colnames(data()[LFQ_columns])
      #label <- gsub('.MaxLFQ.Intensity', '', label)
      #label <- gsub('.Intensity', '', label)
      condition<-""
      #substrRight <- function(x, n){
      #  substr(x, nchar(x)-n+1, nchar(x))
      #}
      replicate<-""
      condition_type<-"1"
      bait_id<-"T"
      exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
    }
    
    
    return(exp_template_new)
  })
  
  exp_template2<- reactive({
    if (input$data_type2 == "MaxQuant"){
      LFQ_columns<- grep("LFQ.",colnames(data2())) 
      label<-colnames(data2()[LFQ_columns])
      label <- gsub('LFQ.intensity.', '', label)
      condition<-""
      substrRight <- function(x, n){
        substr(x, nchar(x)-n+1, nchar(x))
      }
      replicate<-substrRight(label, 1)
      condition_type<-"1"
      bait_id<-"T"
      exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
    }
    
    else if (input$data_type2 == "Proteome Discoverer"){
      Sample_columns<- grep("Abundance..*Sample..", colnames(data2()))
      Control_columns<- grep("Abundance..*Control..", colnames(data2()))
      LFQ_columns<-c(Sample_columns, Control_columns)
      label<-colnames(data2()[LFQ_columns])
      condition<-""
      replicate<-""
      condition_type<-"1"
      bait_id<-"T"
      exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
    }
    
    else if (input$data_type2 == "FragPipe"){
      if(input$exp2 == "LFQ")
      {
        LFQ_columns<- grep(".Intensity", colnames(data2()))
        label<-colnames(data2()[LFQ_columns])
        label <- gsub('.MaxLFQ.Intensity', '', label)
        label <- gsub('.Intensity', '', label)
        condition<-""
        substrRight <- function(x, n){
          substr(x, nchar(x)-n+1, nchar(x))
        }
        replicate<-substrRight(label, 1)
        condition_type<-"1"
        bait_id<-"T"
        exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
      }else if (input$exp == "DIA")
      {
        end<-ncol(data2())
        end_1<-end-1
        LFQ_columns<- 4:end_1
        label<-colnames(data2()[LFQ_columns])
        #label <- gsub('.MaxLFQ.Intensity', '', label)
        #label <- gsub('.Intensity', '', label)
        condition<-""
        #substrRight <- function(x, n){
        #  substr(x, nchar(x)-n+1, nchar(x))
        #}
        replicate<-""
        condition_type<-"1"
        bait_id<-"T"
        exp_template_new<- data.frame(label, condition, replicate, condition_type, bait_id)
        
      }
    }
    
    return(exp_template_new)
  })
  
  
  #Experimental table download
  output$downloadexp_template <- downloadHandler(
    filename = "ExperimentalTemplate.txt",
    content = function(file) {
      write.table(exp_template(),
                  file,
                  col.names =TRUE,
                  row.names=FALSE,
                  quote=FALSE,
                  sep ="\t") }
  )
  
  output$downloadexp_template2 <- downloadHandler(
    filename = "ExperimentalTemplate2.txt",
    content = function(file) {
      write.table(exp_template2(),
                  file,
                  col.names =TRUE,
                  row.names=FALSE,
                  quote=FALSE,
                  sep ="\t") }
  )
  
  
  #Comparisons Template
  comp_template<- reactive({
    condition<-""
    condition2<-""
    comp1_A<-""
    comp1_B<-""
    comp2_A<-""
    comp2_B<-""
    inFile <- input$file2
    inFile2 <- input$file22
    expfile<-expdesign()
    expfile2<-expdesign2()
    if (!is.null(inFile)){condition<-unique(expfile$condition)}
    
    if (!is.null(inFile2)){condition2<-unique(expfile2$condition)
    if (length(condition)>length(condition2)){length(condition2)<-length(condition)}
    if (length(condition2)>length(condition)){length(condition)<-length(condition2)}
    }
    comp_template_new<- data.frame(condition, comp1_A, comp1_B, comp2_A, comp2_B, condition2)
    
    return(comp_template_new)
  })
  
  #Comparison table download
  output$downloadcomp_template <- downloadHandler(
    filename = "ComparisonsTemplate.txt",
    content = function(file) {
      write.table(comp_template(),
                  file,
                  col.names =TRUE,
                  row.names=FALSE,
                  quote=FALSE,
                  sep ="\t") }
  )
  
  #data annotations for clusters
  data_defn <- reactive({
    data_unique <- make_unique(data(), "Gene.names", "Protein.IDs", delim = ";")
    data_unique_defn <-data_unique[, c("name", "Protein.names", "ID")]
    data_unique_defn$Protein.names  <- sub('\\;.*', '', data_unique_defn$Protein.names)
    data_unique_defn<-as.data.frame(data_unique_defn)
  })
  
  data_defn2 <- reactive({
    data_unique <- make_unique(data2(), "Gene.names", "Protein.IDs", delim = ";")
    data_unique_defn <-data_unique[, c("name", "Protein.names", "ID")]
    data_unique_defn$Protein.names  <- sub('\\;.*', '', data_unique_defn$Protein.names)
    data_unique_defn<-as.data.frame(data_unique_defn)
  })
  
  #Experimental design files input
  expdesign <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE,
             sep = "\t", stringsAsFactors = FALSE, colClasses = c("character")) %>%
      mutate(id = row_number())
  })
  
  expdesign2 <- reactive({
    inFile <- input$file22
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE,
             sep = "\t", stringsAsFactors = FALSE, colClasses = c("character")) %>%
      mutate(id = row_number())
  })
  
  
  #comparisons file input
  compfile <- reactive({
    inFile <- input$comp_file
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE,
             sep = "\t", stringsAsFactors = FALSE) %>%
      mutate(id = row_number())
  })
  
  #Filter the Data using missing value threshold
  filt <- reactive({
    data <- data()
    
    if (input$data_type == "MaxQuant"){
      cols <- grep("^LFQ", colnames(data))
      if('Reverse' %in% names(data)==TRUE && 'Potential.contaminant' %in% names(data)==TRUE){
        filtered <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
      }else{filtered<-data}
      unique_names <- make_unique(filtered, "Gene.names", "Protein.IDs", delim = ";")
    }
    
    else if (input$data_type == "Proteome Discoverer"){
      Sample_columns<- grep("Abundance..*Sample..", colnames(data()))
      Control_columns<- grep("Abundance..*Control..", colnames(data()))
      cols<-c(Sample_columns, Control_columns)
      unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
    }
    
    else if (input$data_type == "FragPipe"){
      if(input$exp == "LFQ")
      {
        cols<- grep(".Intensity", colnames(data))
        data<-data %>% 
          rename_with(~ gsub(".MaxLFQ.Intensity", "", .x, fixed = TRUE))
        data<-data %>% 
          rename_with(~ gsub(".Intensity", "", .x, fixed = TRUE))
        unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
      }else if(input$exp == "DIA")
      {
        end<-ncol(data())
        end_1<-end-1
        cols<- 4:end_1
        unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
      }
      
    }
    
    else if (input$data_type == "Other (Experimental)"){
      end<-ncol(data())
      cols<- 5:end
      unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
    }
    
    #Check for Experimental Design File
    validate(need(!is.null(expdesign()), "You need to upload an Experimental Design File.
    Use the Experimental template to assist you. ")
    )
    
    se <- make_se(unique_names, cols, expdesign())
    filter_missval(se, input$thr) #filtering set less stringent
  })
  
  filt2 <- reactive({
    data <- data2()
    
    if (input$data_type2 == "MaxQuant"){
      cols <- grep("^LFQ", colnames(data))
      #filtered <- DEP:::filter_MaxQuant(data, input$filt)
      filtered <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
      #unique_names <- make_unique(filtered, input$name, input$id)
      unique_names <- make_unique(filtered, "Gene.names", "Protein.IDs", delim = ";")
    }
    
    else if (input$data_type2 == "Proteome Discoverer"){
      Sample_columns<- grep("Abundance..*Sample..", colnames(data2()))
      Control_columns<- grep("Abundance..*Control..", colnames(data2()))
      cols<-c(Sample_columns, Control_columns)
      unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
    }
    
    else if (input$data_type2 == "FragPipe"){
      if(input$exp2 == "LFQ")
      {
        cols<- grep(".Intensity", colnames(data))
        data<-data %>% 
          rename_with(~ gsub(".MaxLFQ.Intensity", "", .x, fixed = TRUE))
        data<-data %>% 
          rename_with(~ gsub(".Intensity", "", .x, fixed = TRUE))
        unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
      }else if(input$exp2 == "DIA")
      {
        end<-ncol(data2())
        end_1<-end-1
        cols<- 4:end_1
        unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
      }
      
    }
    
    else if (input$data_type2 == "Other (Experimental)"){
      end<-ncol(data2())
      cols<- 5:end
      unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
    }
    
    #Check for Experimental Design File
    validate(need(!is.null(expdesign2()), "You need to upload an Experimental Design File for 2nd dataset.
    Use the Experimental template to assist you. ")
    )
    
    se <- make_se(unique_names, cols, expdesign2())
    filter_missval(se, input$thr) #filtering set less stringent
  })
  
  norm <- reactive({
    normalize_custom(filt(), method = input$normalisation_choice)
#    if (input$normalisation_choice=="vsn")
#    {normalize_vsn(filt())}
#    else if(input$normalisation_choice=="quantile"){
#      QFeatures::normalize(filt(), method = "quantiles")}
#    else if (input$normalisation_choice=="none")
#    {norm<-filt()}
  })
  
  norm2 <- reactive({
    normalize_custom(filt2(), method = input$normalisation_choice2)
    
#    if (input$normalisation_choice2=="vsn")
#    {normalize_vsn(filt2())}
#    else if(input$normalisation_choice2=="quantile"){
#      QFeatures::normalize(filt2(), method = "quantiles")}
#    else if (input$normalisation_choice2=="none")
#    {norm<-filt2()}
  })
  
  #Get Normalized data
  norm_data <- reactive({
    data_norm<-SummarizedExperiment::assay(norm())
    data_norm<-as.data.frame(data_norm)
    return(data_norm)
  })
  
  norm_data_table <- reactive({
    rownames_to_column(norm_data(), var = "name")
  })
  
  
  norm_data2 <- reactive({
    data_norm<-SummarizedExperiment::assay(norm2())
    data_norm<-as.data.frame(data_norm)
    return(data_norm)
  })
  
  
  imp <- reactive({
    if(input$impute_ave == TRUE){
      set.seed(123)
      if (input$imputation == "man"){
        imp1<-impute_custom(norm(), fun = "man", shift = input$shift, scale = input$width)}else{
          imp1<-impute_custom(norm(), input$imputation)}
      imp1_data<-SummarizedExperiment::assay(imp1)
      
      set.seed(321)
      if (input$imputation == "man"){
        imp2<-impute_custom(norm(), fun = "man", shift = input$shift, scale = input$width)}else{
          imp2<-impute_custom(norm(), input$imputation)}
      imp2_data<-SummarizedExperiment::assay(imp2)
      
      set.seed(456)
      if (input$imputation == "man"){
        imp3<-impute_custom(norm(), fun = "man", shift = input$shift, scale = input$width)}else{
          imp3<-impute_custom(norm(), input$imputation)}
      imp3_data<-SummarizedExperiment::assay(imp3)
      
      imp_ave<-(imp1_data + imp2_data + imp3_data)/3
      SummarizedExperiment::assay(imp1)<-imp_ave
      #set.seed(123)
      return(imp1)
    }else{set.seed(123)
      if (input$imputation == "man"){
        impute_custom(norm(), fun = "man", shift = input$shift, scale = input$width)}else{
          impute_custom(norm(), input$imputation)}
    }
  })
  
  imp2 <- reactive({
    if(input$impute_ave2 == TRUE){
      set.seed(123)
      if (input$imputation2 == "man"){
        imp1<-impute_custom(norm2(), fun = "man", shift = input$shift, scale = input$width)}else{
          imp1<-impute_custom(norm2(), input$imputation)}
      imp1_data<-SummarizedExperiment::assay(imp1)
      
      set.seed(321)
      if (input$imputation2 == "man"){
        imp2<-impute_custom(norm2(), fun = "man", shift = input$shift, scale = input$width)}else{
          imp2<-impute_custom(norm2(), input$imputation)}
      imp2_data<-SummarizedExperiment::assay(imp2)
      
      set.seed(456)
      if (input$imputation2 == "man"){
        imp3<-impute_custom(norm2(), fun = "man", shift = input$shift, scale = input$width)}else{
          imp3<-impute_custom(norm2(), input$imputation)}
      imp3_data<-SummarizedExperiment::assay(imp3)
      
      imp_ave<-(imp1_data + imp2_data + imp3_data)/3
      SummarizedExperiment::assay(imp1)<-imp_ave
      #set.seed(123)
      return(imp1)
    }else{set.seed(123)
      if (input$imputation2 == "man"){
        impute_custom(norm2(), fun = "man", shift = input$shift, scale = input$width)}else{
          impute_custom(norm2(), input$imputation)}
    }
  })
  
  
  #Get imputed data table
  imp_data <- reactive({
    data_imp<-SummarizedExperiment::assay(imp())
    return(as.data.frame(data_imp))
  })
  
  imp_data_table <- reactive({
    rownames_to_column(imp_data(), var = "name")
  })
  
  
  #Get imputed data table
  imp_data2 <- reactive({
    data_imp2<-SummarizedExperiment::assay(imp2())
    return(as.data.frame(data_imp2))
  })
  
  
  comp <- reactive({
    inFile <- input$comp_file
    if (is.null(inFile)){
      #      comp<-c(input$comp1f,input$comp1r)
      comp<-c(input$comp1f)
    }else{
      comp_input<-data.frame(compfile())
      comp<-paste0(comp_input$comp1_A, "_vs_", comp_input$comp1_B)
      comp<-comp[comp != "_vs_"]
      comp<-comp[comp != ""]
    }
  })
  
  comp2 <- reactive({
    if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons")
      return(NULL)
    inFile <- input$comp_file
    if (is.null(inFile)){
      #      comp2<-c(input$comp2f,input$comp2r)
      comp2<-c(input$comp2f)
      
    }else{
      comp_input<-data.frame(compfile())
      comp2<-paste0(comp_input$comp2_A, "_vs_", comp_input$comp2_B)
      comp2<-comp2[comp2 != "_vs_"]
      comp2<-comp2[comp2 != ""]
    }
    
  })
  
  #Comparisons Table
  comp_table<- reactive({
    R1 <- data.frame(comp())
    R2 <- data.frame(comp2())
    comp_table_new<- R1 %>%
      rownames_to_column() %>% 
      left_join(R2 %>% rownames_to_column()) %>%
      select(-rowname)
    
    #     return(comp_template_new)
  })
  
  
  df <- reactive({
    if (input$padjust_method == "default"){
      DEP::test_diff(imp(), type = "manual", test = comp())}
    else{test_diff_custom(imp(), type = "manual", test = comp(), method = input$padjust_method)}
    
    
  })
  
  df2 <- reactive({
    if (input$comp_type == "2 Rounds of comparisons"){
      if (input$padjust_method == "default"){
        DEP::test_diff(imp(), type = "manual", test = comp2())}
      else{test_diff_custom(imp(), type = "manual", test = comp2(), method = input$padjust_method)}
    }
    else{
      if (input$padjust_method2 == "default"){
        DEP::test_diff(imp2(), type = "manual", test = comp2())}
      else{test_diff_custom(imp2(), type = "manual", test = comp2(), method = input$padjust_method2)}
    }
  })
  
  dep <- reactive({
    #Check for Experimental Design File
    validate(need(!is.null(expdesign()), "You need to upload an Experimental Design File.
    Use the Experimental template to assist you. ")
    )
    
    #Check for Comparisons
    validate(need(comp() != "", "You need to select some comparisons, or import a comparisons file.
          Use the Comp. template to assist you. ")
    )
    #validate(need(try(df() != "try-error"), base::stop(shiny::safeError(warning(test_diff(imp(), type = "manual", test = comp()))))))
    if(input$pval_use == FALSE){
      DEP::add_rejections(df(), input$p, log2(input$lfc))
    }else{add_rejections_pval(df(), input$p, log2(input$lfc))}
  }) 
  
  dep2 <- reactive({
    if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons")
      return(NULL)
    #Check for Comparisons
    validate(need(comp2() != "", "You need to select some 2nd round comparisons, or import a comparisons file. Use the Comp. template to assist you. ")
    )
    if(input$pval_use == FALSE){
      DEP::add_rejections(df2(), input$p2, log2(input$lfc2))
    }else{add_rejections_pval(df2(), input$p2, log2(input$lfc2))}
    
  }) 
  
  
  dep_num2 <- reactive({
    if (input$comp_type == "1 Round of comparisons")
      return(NULL)
    num_signif <- dep2() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  })
  
  df_long_input <- reactive({
    norm_long <-DEP::get_df_long(norm())
  })
  
  df_long_input2 <- reactive({
    norm_long2 <-DEP::get_df_long(norm2())
  })
  
  
  
  #SAINTexpress Processing
  saint_output <- reactive({
    
    experimental_design<-data.frame(expdesign())
    #Check for C's
    validate(need('C' %in% experimental_design$bait_id, "For SAINTexpress processing, you must indicate the control samples with a 'C' in the 'bait_id' column of the Experimental Design file. Use the ExperimentalTemplate to assist you. ")
    )
    SysName <- Sys.info()['sysname']
    norm_long <-norm_output()
    saint_list<-list()
    
    for(i in sort(unique(experimental_design$condition_type))){
      #norm_df<-pg_filt[, c(1, grep("Norm",colnames(pg_filt)))]
      
      #bait
      bait<-experimental_design[experimental_design$condition_type==i, c(1,2,grep("bait_id",colnames(experimental_design)))]
      baitfile<-paste("data/bait_", i, ".txt",sep="")
      write.table(bait, baitfile, sep="\t", quote = F, row.names = F, col.names = F)
      
      #prey
      prey<-data.frame(PROT_ID=norm_long$ID,
                       LENGTH=norm_long$Sequence.length ,
                       GENE_ID=norm_long$name)
      prey$PROT_ID  <- sub('\\;.*', '', prey$PROT_ID)
      prey$GENE_ID  <- sub('\\;.*', '', prey$GENE_ID)
      prey<-prey %>% 
        distinct(across(-PROT_ID), # reduces data frame to only unique rows (keeps first one of any duplicates)
                 .keep_all = TRUE) 
      
      preyfile<-paste("data/prey_", i, ".txt",sep="")
      write.table(prey, preyfile, sep="\t", quote = F, row.names = F, col.names = F)
      
      #inter
      inter<-norm_long[norm_long$condition_type==i, c(grep("label", colnames(norm_long)),grep("\\condition\\b", colnames(norm_long)),grep("\\ID\\b", colnames(norm_long)), grep ("intensity", colnames(norm_long)) )]
      inter$ID  <- sub('\\;.*', '', inter$ID)
      
      #remove NA's and reverse log2 transformation
      inter<-inter[!is.na(inter$intensity),]
      inter$intensity<-2^(inter$intensity)
      
      ##reorder columns for ease
      #inter<-inter[, c("EXP", "Bait", "PROT_ID", "INT")]
      
      interfile<-paste("data/inter_", i, ".txt",sep="")
      write.table(inter, interfile, sep="\t", quote = F, row.names = F, col.names = F)
      
      #run SAINT (-int for intensities and -spc for spectral counts)
      saintfile<-paste("SAINT_", i, ".txt",sep="")
      if(SysName=="Linux"){saint_location<-"data/SAINTexpress-int"}
      #if(SysName=="Windows"){saint_location<-"data/SAINTexpress-spc.exe"}
      if(SysName=="Windows"){saint_location<-"data/SAINTexpress-int.exe"}
      
      system(paste(saint_location,interfile, preyfile, baitfile))
      saint_df<-read.csv("list.txt", head=T, sep="\t", stringsAsFactors = FALSE)
      if(SysName=="Linux"){system("rm -r data/*.txt")}
      #cmd<-"del data\ *.txt"
      #if(SysName=="Windows"){system(paste0(cmd))}
      
      saint_list[[i]]<-saint_df
      if(SysName=="Linux"){system("rm -r list.txt")}
    }
    return(saint_list)
    
  })
  
  norm_output <- reactive({
    norm_output <-get_df_long(norm())
  })
  
  saintlist <- reactive({
    experimental_design<-data.frame(expdesign())
    cond_num<-length(unique(experimental_design$condition_type))
    
    if(cond_num==1){unfilt_saint<-saint_output()[[1]]}
    if(cond_num==2){unfilt_saint<-rbind(saint_output()[[1]],saint_output()[[2]])}
    if(cond_num==3){unfilt_saint<-rbind(saint_output()[[1]],saint_output()[[2]],saint_output()[[3]])}
    if(cond_num==4){unfilt_saint<-rbind(saint_output()[[1]],saint_output()[[2]],saint_output()[[3]],saint_output()[[4]])}
    
    
    #unfilt_saint<- c(saintfile1()) #Make single list
    #if (!is.null(saintfile2())){unfilt_saint<- rbind(saintfile1(),saintfile2())}
    #if (!is.null(saintfile3())){unfilt_saint<- rbind(saintfile1(),saintfile2(),saintfile3())}
    unfilt_saint<-data.frame(unfilt_saint)
    unfilt_saint<-unfilt_saint %>% dplyr::filter(BFDR <= input$bfdr)
    #    return(unfilt_saint)
    
    saint_count<-as.data.frame(table(unfilt_saint$PreyGene)) #Make into dataframe
    saint_count<-saint_count %>%dplyr::rename(PreyGene=Var1)%>% dplyr::rename(BFDR_Freq=Freq) %>% filter(PreyGene !=(""))
    saint_count <- saint_count[saint_count$BFDR_Freq>=input$bfdr_freq, ]#filter by BFDR_Freq
    saint_count_no <-nrow(saint_count)
    saintlist_out<-list(unfilt_saint,saint_count,saint_count_no)
    return(saintlist_out)
    
  })
  
  
  ### All object and functions upon 'Analyze' input ----------------------
  
  observeEvent(input$analyze, {
    
    #### Interactive UI functions  ------------------------------------------
    output$downloadTable <- renderUI({
      selectizeInput("dataset",
                     "Choose a dataset to save" ,
                     c("results","significant_proteins",
                       "displayed_subset","full_dataset","Cluster_list", "normalised_values", "imputed_values", "results_long", "results_long_signif"))
    })
    
    output$downloadButton <- renderUI({
      downloadButton('downloadData', 'Save', class="btn btn-primary")
    })
    
    output$significantBox <- renderInfoBox({
      num_total <- dep() %>%
        nrow()
      num_signif <- dep() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()
      frac <- num_signif / num_total
      
      if(frac > 0.2) {
        info_box <- infoBox("Significant proteins (Comp1)",
                            paste0(num_signif,
                                   " out of ",
                                   num_total),
                            paste0("Too large fraction (",
                                   signif(frac * 100, digits = 3),
                                   "%) of proteins differentially expressed"),
                            icon = icon("minus", lib = "glyphicon"),
                            color = "orange",
                            width = 2)
      }
      if(frac == 0) {
        info_box <- infoBox("Significant proteins (Comp1)",
                            paste0(num_signif,
                                   " out of ",
                                   num_total),
                            "No proteins differentially expressed",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 2)
      }
      if(frac > 0 & frac <= 0.2) {
        info_box <- 		infoBox("Significant proteins (Comp1)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0(signif(frac * 100, digits = 3),
                                     "% of proteins differentially expressed"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 2)
      }
      info_box
    })
    
    output$significantBox2 <- renderInfoBox({
      if (input$comp_type == "2 Rounds of comparisons" | input$data_mode == TRUE)
      {
        num_total <- dep2() %>%
          nrow()
        num_signif <- dep2() %>%
          .[SummarizedExperiment::rowData(.)$significant, ] %>%
          nrow()
        frac <- num_signif / num_total
        
        if(frac > 0.2) {
          info_box <- infoBox("Significant proteins (Comp2)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0("Too large fraction (",
                                     signif(frac * 100, digits = 3),
                                     "%) of proteins differentially expressed"),
                              icon = icon("minus", lib = "glyphicon"),
                              color = "orange",
                              width = 2)
        }
        if(frac == 0) {
          info_box <- infoBox("Significant proteins (Comp2)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              "No proteins differentially expressed",
                              icon = icon("thumbs-down", lib = "glyphicon"),
                              color = "red",
                              width = 2)
        }
        if(frac > 0 & frac <= 0.2) {
          info_box <- 		infoBox("Significant proteins (Comp2)",
                                paste0(num_signif,
                                       " out of ",
                                       num_total),
                                paste0(signif(frac * 100, digits = 3),
                                       "% of proteins differentially expressed"),
                                icon = icon("thumbs-up", lib = "glyphicon"),
                                color = "green",
                                width = 2)
        }
      }
      info_box
    })
    
    output$significantBox3 <- renderInfoBox({
      if (input$comp_type == "2 Rounds of comparisons" | input$data_mode == TRUE)
      {
        num_total <- dep() %>%
          nrow()
        num_signif <- data_results3() %>%
          nrow()
        frac <- num_signif / num_total
        
        if(frac > 0.2) {
          info_box <- infoBox("Significant proteins (Overlap)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0("Too large fraction (",
                                     signif(frac * 100, digits = 3),
                                     "%) of proteins differentially expressed"),
                              icon = icon("minus", lib = "glyphicon"),
                              color = "orange",
                              width = 2)
        }
        if(frac == 0) {
          info_box <- infoBox("Significant proteins (Overlap)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              "No proteins differentially expressed",
                              icon = icon("thumbs-down", lib = "glyphicon"),
                              color = "red",
                              width = 2)
        }
        if(frac > 0 & frac <= 0.2) {
          info_box <- 		infoBox("Significant proteins (Overlap)",
                                paste0(num_signif,
                                       " out of ",
                                       num_total),
                                paste0(signif(frac * 100, digits = 3),
                                       "% of proteins differentially expressed"),
                                icon = icon("thumbs-up", lib = "glyphicon"),
                                color = "green",
                                width = 2)
        }
      }
      info_box
    })
    
    output$significantBox_bfdr <- renderInfoBox({
      if (input$comp_type == "Use SAINTexpress")
      {
        num_total <- dep() %>%
          nrow()
        num_signif <- saintlist()[[2]] %>%
          nrow()
        frac <- num_signif / num_total
        
        if(frac > 0.2) {
          info_box <- infoBox("Significant proteins (BFDR)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0("Too large fraction (",
                                     signif(frac * 100, digits = 3),
                                     "%) of proteins differentially expressed"),
                              icon = icon("minus", lib = "glyphicon"),
                              color = "orange",
                              width = 2)
        }
        if(frac == 0) {
          info_box <- infoBox("Significant proteins (BFDR)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              "No proteins differentially expressed",
                              icon = icon("thumbs-down", lib = "glyphicon"),
                              color = "red",
                              width = 2)
        }
        if(frac > 0 & frac <= 0.2) {
          info_box <- 		infoBox("Significant proteins (BFDR)",
                                paste0(num_signif,
                                       " out of ",
                                       num_total),
                                paste0(signif(frac * 100, digits = 3),
                                       "% of proteins differentially expressed"),
                                icon = icon("thumbs-up", lib = "glyphicon"),
                                color = "green",
                                width = 2)
        }
      }
      info_box
    })
    
    output$significantBox_bfdr_filtered <- renderInfoBox({
      if (input$comp_type == "Use SAINTexpress")
      {
        num_total <- dep() %>%
          nrow()
        num_signif <- data_results5() %>%
          nrow()
        frac <- num_signif / num_total
        
        if(frac > 0.2) {
          info_box <- infoBox("Significant proteins (BFDR Overlap)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0("Too large fraction (",
                                     signif(frac * 100, digits = 3),
                                     "%) of proteins differentially expressed"),
                              icon = icon("minus", lib = "glyphicon"),
                              color = "orange",
                              width = 2)
        }
        if(frac == 0) {
          info_box <- infoBox("Significant proteins (BFDR Overlap)",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              "No proteins differentially expressed",
                              icon = icon("thumbs-down", lib = "glyphicon"),
                              color = "red",
                              width = 2)
        }
        if(frac > 0 & frac <= 0.2) {
          info_box <- 		infoBox("Significant proteins (BFDR Overlap)",
                                paste0(num_signif,
                                       " out of ",
                                       num_total),
                                paste0(signif(frac * 100, digits = 3),
                                       "% of proteins differentially expressed"),
                                icon = icon("thumbs-up", lib = "glyphicon"),
                                color = "green",
                                width = 2)
        }
      }
      info_box
    })
    
    
    #  output$ibox <- renderInfoBox({
    #    infoBox(
    #      "Title",
    #      input$file_type,
    #      #icon = icon("credit-card")
    #    )
    #  })
    
    
    output$volcano_cntrst <- renderUI({
      if (!is.null(dep())) {
        df <- SummarizedExperiment::rowData(dep())
        cols <- grep("_significant$",colnames(df))
        selectizeInput("volcano_cntrst",
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]))
      }
    })
    
    output$volcano_cntrst_c2 <- renderUI({
      if (!is.null(dep())) {
        selectizeInput("volcano_cntrst_c2",
                       "Contrast",
                       choices = comp2())
      }
    })
    
    
    #Selection for Enhanced volcano1
    output$volcano_cntrst_enh1 <- renderUI({
      df_choice <- comp()
      if (input$comp_type == "1 Round of comparisons")
      {df_choice<-comp()}
      if (input$comp_type == "2 Rounds of comparisons")
      {df_choice<-comp2()}
      
      
      if (!is.null(dep())) {
        selectizeInput("volcano_cntrst_enh1",
                       "Contrast",
                       choices = df_choice)
      }
    })
    
    #Selection for Enhanced volcano2
    output$volcano_cntrst_enh2 <- renderUI({
      df_choice <- comp()
      if (input$comp_type == "1 Round of comparisons")
      {df_choice<-comp()}
      if (input$comp_type == "2 Rounds of comparisons")
      {df_choice<-comp2()}
      
      if (!is.null(dep())) {
        selectizeInput("volcano_cntrst_enh2",
                       "Contrast",
                       choices = df_choice)
      }
    })
    
    
    
    #### Reactive functions  (Generate plots and tables)------------------------
    
    
    #Make results table 
    data_dep <- reactive({
      get_results(dep()) %>% filter(significant) %>% select(-significant) 
    })
    
    data_dep2 <- reactive({
      if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons")
        return(NULL)
      if (input$data_mode == FALSE && input$comp_type == "Use SAINTexpress")
        return(NULL)
      get_results(dep2()) %>% filter(significant) %>% select(-significant) 
    })
    
    data_results1<- reactive({
      data.frame(data_dep())
    })
    
    data_results2<- reactive({
      if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons")
        return(NULL)
      if (input$data_mode == FALSE && input$comp_type == "Use SAINTexpress")
        return(NULL)
      data.frame(data_dep2())
    })
    
    data_results3<- reactive({
      #      if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons")
      #        return(NULL)
      #      if (input$data_mode == FALSE && input$comp_type == "Use SAINTexpress")
      #        return(NULL)
      if (input$comp_type == "2 Rounds of comparisons"){
        data_results2() [data_results2()$name %in% data_results1()$name, ]
      }
      else{
        data_results1() [data_results1()$name %in% data_results2()$name, ]
      }
    })
    
    #Unfiltered BFDR results
    data_results4<- reactive({
      if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons" | input$comp_type == "2 Rounds of comparisons")
        return(NULL)
      data<-get_results(dep())
      saint_data<-data.frame(saintlist()[[2]])
      data<-data[data$name %in% saint_data$PreyGene, ]
      data1 <- data %>% arrange(name)
      return(data1)
    })
    
    data_results5<- reactive({
      if (input$comp_type == "1 Round of comparisons" | input$comp_type == "2 Rounds of comparisons")
        return(NULL)
      data_results4() [data_results4()$name %in% data_results1()$name, ]
    })
    
    
    data_results<- reactive({
      if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons")
      {data<-data_results1()}
      if (input$data_mode == FALSE && input$comp_type == "2 Rounds of comparisons"){
        data<-data_results3()}
      if (input$data_mode == FALSE && input$comp_type == "Use SAINTexpress"){
        data<-data_results5()}
      if (input$data_mode == TRUE){
        data<-data_results3()}
      return(data)
    })
    
    dep_num <- reactive({
      num_signif <- data_results() %>% nrow()
    })
    
    dep_num_c1 <- reactive({
      num_signif <- data_results1() %>% nrow()
    })
    
    dep_num_c2 <- reactive({
      num_signif <- data_results2() %>% nrow()
    })
    
    dep_num_bfdr <- reactive({
      num_signif <- data_results4() %>% nrow()
    })
    
    #Selected Protein Plots
    selected_plot_input <- reactive ({
      if(!is.null(input$table_rows_selected)) {
        if (input$protein_plot == "contrast" | input$protein_plot == "centered"){
          if (input$protein_plot == "contrast"){validate(need(length(comp()) != 1, "   You need to select a second comparison for this plot to generate"))}
          if (input$comp_type == "2 Rounds of comparisons"){
            if (input$top_table_opt == TRUE){selected_id <- get_results(dep2())[input$table_rows_selected,1]}
            else{selected_id <- data_results()[input$table_rows_selected,1]}
            plot_single(dep2(), selected_id, input$protein_plot)}
          else if (input$top_table_opt == TRUE){selected_id <- get_results(dep())[input$table_rows_selected,1]
          plot_single(dep(), selected_id, input$protein_plot)}
          else{selected_id <- data_results()[input$table_rows_selected,1]
          plot_single(dep(), selected_id, input$protein_plot)}
        }else{
          if (input$comp_type == "2 Rounds of comparisons"){
            if (input$top_table_opt == TRUE){GeneID <- get_results(dep2())[input$table_rows_selected,1]}
            else{GeneID <- data_results()[input$table_rows_selected,1]}
          }
          else if (input$top_table_opt == TRUE){GeneID <- get_results(dep())[input$table_rows_selected,1]
          }
          else{GeneID <- data_results()[input$table_rows_selected,1]
          }
          Cluster <- "cluster0"
          clu_df_temp<-data.frame(GeneID, Cluster)
          plot_box_plot (norm_data(),
                         imp_data(),
                         clu_df=clu_df_temp,
                         cluster_no=0,
                         input$palette_choice
          )
        }
      }
    })
    
    selected_plot_input2 <- reactive ({
      if(!is.null(input$table_rows_selected)) {
        if (input$protein_plot2 == "contrast" | input$protein_plot2 == "centered"){
          if (input$protein_plot2 == "contrast"){validate(need(length(comp2()) != 1, "   You need to select a second comparison for this plot to generate"))}
          if (input$top_table_opt == TRUE){selected_id <- get_results(dep())[input$table_rows_selected,1]}
          else{selected_id <- data_results()[input$table_rows_selected,1]}
          plot_single(dep2(), selected_id, input$protein_plot2)
        }else{
          if (input$top_table_opt == TRUE){GeneID <- get_results(dep())[input$table_rows_selected,1]}
          else{GeneID <- data_results()[input$table_rows_selected,1]}
          Cluster <- "cluster0"
          clu_df_temp<-data.frame(GeneID, Cluster)
          plot_box_plot (norm_data2(),
                         imp_data2(),
                         clu_df=clu_df_temp,
                         cluster_no=0,
                         input$palette_choice
          )
        }
      }
    })
    
    
    #Abundance Heatmap
    heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if(input$abund_choice == "imputed"){
          abund_matrix<-imp_data()
          if(input$signif_choice==TRUE){abund_matrix<- abund_matrix[rownames(abund_matrix) %in% data_results()$name, ]}
          protein_no<-abund_matrix %>% nrow()
          if (protein_no>input$size){protein_no<-input$size}
          abund_matrix<- abund_matrix[ order(rowMeans(abund_matrix), decreasing = T), ]
          abund_matrix<-abund_matrix[1:protein_no,]
          abund_matrix<-as.matrix(abund_matrix)
          switch_no<-1
        }
        
        if(input$abund_choice == "normalised"){
          abund_matrix<-norm_data()
          if(input$signif_choice==TRUE){abund_matrix<- abund_matrix[rownames(abund_matrix) %in% data_results()$name, ]}
          protein_no<-abund_matrix %>% nrow()
          if (protein_no>input$size){protein_no<-input$size}
          #abund_matrix[is.na(abund_matrix)] = 0
          abund_matrix<- abund_matrix[ order(rowMeans(abund_matrix,na.rm=TRUE), decreasing = T), ]
          abund_matrix[abund_matrix==0]<-NA
          abund_matrix<-abund_matrix[1:protein_no,]
          abund_matrix<-as.matrix(abund_matrix)
          switch_no<-2
        }
        
        if(input$signif_choice==TRUE){sign_txt<-"significant proteins only"}else{sign_txt<-"all proteins"}
        ht_title<-paste(input$data_type, "<br>", protein_no, " most abundant proteins", "<br>", input$abund_choice, " intensity values", "<br>", sign_txt, sep = "")
        
        ht<-plot_complex_heatmap(abund_matrix,
                                 #                     input$pres,
                                 #                     kmeans = TRUE,
                                 kmeans=1,
                                 ht_title,
                                 colour_functions_input()[[1]],
                                 colour_functions_input()[[2]],
                                 switch = switch_no,
                                 font_size=8,
                                 cluster_cols=FALSE,
                                 comp_list=comp()
        )
      })
    })
    
    heatmap_input2 <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if(input$abund_choice2 == "imputed"){
          abund_matrix<-imp_data2()
          if(input$signif_choice2==TRUE){abund_matrix<- abund_matrix[rownames(abund_matrix) %in% data_results()$name, ]}
          protein_no<-abund_matrix %>% nrow()
          if (protein_no>input$size2){protein_no<-input$size2}
          abund_matrix<- abund_matrix[ order(rowMeans(abund_matrix), decreasing = T), ]
          abund_matrix<-abund_matrix[1:protein_no,]
          abund_matrix<-as.matrix(abund_matrix)
          switch_no<-1
        }
        
        if(input$abund_choice2 == "normalised"){
          abund_matrix<-norm_data2()
          if(input$signif_choice2==TRUE){abund_matrix<- abund_matrix[rownames(abund_matrix) %in% data_results()$name, ]}
          protein_no<-abund_matrix %>% nrow()
          if (protein_no>input$size2){protein_no<-input$size2}
          #abund_matrix[is.na(abund_matrix)] = 0
          abund_matrix<- abund_matrix[ order(rowMeans(abund_matrix,na.rm=TRUE), decreasing = T), ]
          abund_matrix[abund_matrix==0]<-NA
          abund_matrix<-abund_matrix[1:protein_no,]
          abund_matrix<-as.matrix(abund_matrix)
          switch_no<-2
        }
        
        if(input$signif_choice2==TRUE){sign_txt<-"significant proteins only"}else{sign_txt<-"all proteins"}
        ht_title<-paste(input$data_type2, "<br>", protein_no, " most abundant proteins", "<br>", input$abund_choice, " intensity values", "<br>", sign_txt, sep = "")
        
        ht<-plot_complex_heatmap(abund_matrix,
                                 #                     input$pres,
                                 #                     kmeans = TRUE,
                                 kmeans=1,
                                 ht_title,
                                 colour_functions_input()[[1]],
                                 colour_functions_input()[[2]],
                                 switch = switch_no,
                                 font_size=8,
                                 cluster_cols=FALSE,
                                 comp_list=comp()
        )
      })
    })
    
    
    #data annotations for abundance table
    data_defn_abundance <- reactive({
      abund_table_matrix<-heatmap_input()[[2]]
      abund_table <- data.frame(name = row.names(abund_table_matrix), abund_table_matrix)
      abund_table_merge <- merge(data_defn(), abund_table, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)
      abund_table_merge[match(abund_table$name, abund_table_merge$name),]
    })
    
    data_defn_abundance2 <- reactive({
      abund_table_matrix<-heatmap_input2()[[2]]
      abund_table <- data.frame(name = row.names(abund_table_matrix), abund_table_matrix)
      abund_table_merge <- merge(data_defn2(), abund_table, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)
      abund_table_merge[match(abund_table$name, abund_table_merge$name),]
    })
    
    volcano_input <- reactive({
      if(!is.null(input$volcano_cntrst)) {
        if (input$comp_type == "2 Rounds of comparisons"){
          if (input$top_table_opt == TRUE){GeneID <- get_results(dep2())[input$table_rows_selected,1]}
          else{GeneID <- data_results()[input$table_rows_selected,1]}
        }
        else if (input$top_table_opt == TRUE){GeneID <- get_results(dep())[input$table_rows_selected,1]
        }
        else{GeneID <- data_results()[input$table_rows_selected,1]
        }
        
        plot_enhanced_volcano_dep (data_results = get_results(dep()),
                                   #get_results(dep()),
                                   #cludf_check_input(),
                                   comp_pick = input$volcano_cntrst,
                                   FC=input$lfc,
                                   pVal=input$p,
                                   padj_check=input$p_adj,
                                   mypal = colour_functions_input()[[2]],
                                   colour_functions_input()[[3]],
                                   #connect=input$enh_vol_connect,
                                   selected=input$check_names, 
                                   selected_names = GeneID,
                                   label_size = input$fontsize
        )
      }
    })
    
    volcano_c2_input <- reactive({
      if(!is.null(input$volcano_cntrst_c2)) {
        if (input$comp_type == "2 Rounds of comparisons"){
          if (input$top_table_opt == TRUE){GeneID <- get_results(dep2())[input$table_rows_selected,1]}
          else{GeneID <- data_results()[input$table_rows_selected,1]}
        }
        else if (input$top_table_opt == TRUE){GeneID <- get_results(dep())[input$table_rows_selected,1]
        }
        else{GeneID <- data_results()[input$table_rows_selected,1]
        }
        
        plot_enhanced_volcano_dep (data_results = get_results(dep2()),
                                   #get_results(dep()),
                                   #cludf_check_input(),
                                   comp_pick = input$volcano_cntrst_c2,
                                   FC=input$lfc2,
                                   pVal=input$p2,
                                   padj_check=input$p_adj_c2,
                                   mypal = colour_functions_input()[[2]],
                                   colour_functions_input()[[3]],
                                   #connect=input$enh_vol_connect,
                                   selected=input$check_names_c2, 
                                   selected_names = GeneID,
                                   label_size = input$fontsize_c2
        )
      }
    })
    
    #Significance Plot
    signif_plot_input <- reactive({
      if(input$pval_use == FALSE){
        signif_count(df(), max_cutoff = input$max_cutoff)
      }else{signif_count_pval(df(), max_cutoff = input$max_cutoff)}
    })
    
    signif_plot_input2 <- reactive({
      if(input$pval_use == FALSE){
        signif_count(df2(), input$max_cutoff2)
      }else{signif_count_pval(df2(), input$max_cutoff2)}
    })
    
    
    norm_input <- reactive({
      if(input$norm_orient == TRUE){plot_normalization_hor(filt(), norm())}else{
        plot_normalization_vert(filt(), norm())}
    })
    
    norm_input2 <- reactive({
      if(input$norm_orient2 == TRUE){plot_normalization_hor(filt2(), norm2())}else{
        plot_normalization_vert(filt2(), norm2())}
    })
    
    
    missval_input <- reactive({
      plot_missval(norm())
    })
    
    detect_input <- reactive({
      plot_detect(norm())
    })
    
    imputation_input <- reactive({
      plot_imputation(norm(),
                      df())
    })
    
    imputation_input2 <- reactive({
      plot_imputation(norm2(),
                      df2())
    })
    
    numbers_input <- reactive({
      plot_numbers(norm())
    })
    
    numbers_input2 <- reactive({
      plot_numbers(norm2())
    })
    
    coverage_input <- reactive({
      plot_coverage(norm())
    })
    
    #PCA Input
    pca_input <- reactive({
      num_total <- dep() %>% nrow()
      if(num_total >= input$pca_num){
        if(input$pca_split==FALSE){
          plot_pca_txt(dep(), x = input$pcx_num, y = input$pcy_num, n=input$pca_num, point_size = 4)
        }else{plot_pca_txt(dep(), x = input$pcx_num, y = input$pcy_num, indicate = c("condition", "replicate", "condition_type"), n=input$pca_num, point_size = 4)}
      }else{
        plot_pca_txt(dep(), x = input$pcx_num, y = input$pcy_num, n=num_total, point_size = 4)
      }
    })
    
    #PCA Input (Comp2)
    pca_c2_input <- reactive({
      num_total <- dep2() %>% nrow()
      if(num_total >= input$pca_num2){
        if(input$pca_split2==FALSE){
          plot_pca_txt(dep2(), x = input$pcx_num2, y = input$pcy_num2, n=input$pca_num2, point_size = 4)
        }else{plot_pca_txt(dep2(), x = input$pcx_num2, y = input$pcy_num2, indicate = c("condition", "replicate", "condition_type"), n=input$pca_num2, point_size = 4)}
      }else{
        plot_pca_txt(dep2(), x = input$pcx_num2, y = input$pcy_num2, n=num_total, point_size = 4)
      }
    })
    
    #Correlation Input
    cor_input <- reactive({
      expfile<-expdesign()
      condition_no<-length(unique(expfile$condition))
      if(condition_no<=12){  
        plot_cor(dep(), significant = TRUE, lower = 0, upper = 1, pal = "Reds", indicate = c("condition", "replicate"), font_size = 10)
      }else{plot_cor(dep(), significant = TRUE, lower = 0, upper = 1, pal = "Reds", font_size = 10)}
    })
    
    #Correlation Input (Comp2)
    cor_c2_input <- reactive({
      expfile<-expdesign()
      condition_no<-length(unique(expfile$condition))
      if(condition_no<=12){  
        plot_cor(dep2(), significant = TRUE, lower = 0, upper = 1, pal = "Reds", indicate = c("condition", "replicate"), font_size = 10)
      }else{plot_cor(dep2(), significant = TRUE, lower = 0, upper = 1, pal = "Reds", font_size = 10)}
    })
    
    #Frequency plot Input
    cond_input <- reactive({
      plot_cond(dep())
    })
    
    #Frequency plot Input
    cond_c2_input <- reactive({
      plot_cond(dep2())
    })
    
    #Colour Functions
    colour_functions_input <- reactive({
      colour_functions(input$limit_comp, input$high_col, input$low_col, input$palette_choice, input$alpha_input)
    }) %>% debounce(3000)
    
    #complex heatmap Input
    comp_heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
        if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons"){
          ht_title<-paste(input$data_type, "<br>", "FC=", input$lfc, p_txt, input$p, "<br>", "(",dep_num_c1(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<-plot_complex_heatmap(data_results(),
                                   #                     input$pres,
                                   #                     kmeans = TRUE,
                                   input$k1,
                                   ht_title,
                                   colour_functions_input()[[1]],
                                   colour_functions_input()[[2]],
                                   font_size=8,
                                   cluster_cols=input$cluster_cols,
                                   comp_list=comp()
          )
        }
        
        if (input$data_mode == FALSE && input$comp_type == "Use SAINTexpress"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          ht_title<-paste(input$data_type, "<br>", "SAINTexpress Filtered", "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "BFDR_Cutoff=", input$bfdr,  "<br>", "BFDR_Freq=", input$bfdr_freq,  "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_BFDR(data_results(),
                                        #                     input$pres,
                                        #                     kmeans = TRUE,
                                        input$k1,
                                        ht_title,
                                        colour_functions_input()[[1]],
                                        colour_functions_input()[[2]],
                                        font_size=8,
                                        cluster_cols=input$cluster_cols,
                                        comp_list=comp(),
                                        data_results1=data_results1(),
                                        saintlist=saintlist()[[2]])
        }
        
        if (input$data_mode == FALSE && input$comp_type == "2 Rounds of comparisons"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          ht_title<-paste(input$data_type, "<br>", "Comparison 1+2", "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "Comp 2 FC=", input$lfc2, p_txt, input$p2, "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<- plot_complex_heatmap_2xcomp(data_results(),
                                           #                     input$pres,
                                           #                     kmeans = TRUE,
                                           input$k1,
                                           ht_title,
                                           colour_functions_input()[[1]],
                                           colour_functions_input()[[2]],
                                           font_size=8,
                                           cluster_cols=input$cluster_cols,
                                           comp_list=comp2(),
                                           data_results1=data_results1(),
                                           data_results2=data_results2())
        }
        
        if (input$data_mode == TRUE){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          if(input$data_type2 == "FragPipe"){data_type2 <- paste(input$data_type2, " ", input$exp2 , sep="")}else{data_type2<-input$data_type2}
          ht_title<-paste(data_type, " and ",  data_type2, "<br>", "Dataset 1+2", "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "Comp 2 FC=", input$lfc2, p_txt, input$p2, "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<- plot_complex_heatmap_2xcomp(data_results(),
                                           #                     input$pres,
                                           #                     kmeans = TRUE,
                                           input$k1,
                                           ht_title,
                                           colour_functions_input()[[1]],
                                           colour_functions_input()[[2]],
                                           font_size=8,
                                           cluster_cols=input$cluster_cols,
                                           comp_list=comp(),
                                           data_results1=data_results1(),
                                           data_results2=data_results2())
        }
        
        return(ht)
      })
    })
    
    comp_heatmap_input_pdf <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          ht_title<-paste(input$data_type, "<br>", "FC=", input$lfc, p_txt, input$p, "<br>", "(",dep_num_c1(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<- plot_complex_heatmap(data_results(),
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1,
                                    ht_title,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=6,
                                    cluster_cols=input$cluster_cols,
                                    comp_list=comp()
          )
        }
        
        if (input$comp_type == "Use SAINTexpress"){
          if(input$data_mode == FALSE && input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          ht_title<-paste(input$data_type, "<br>", "SAINTexpress Filtered", "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "BFDR_Cutoff=", input$bfdr,  "<br>", "BFDR_Freq=", input$bfdr_freq, "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<- plot_complex_heatmap_BFDR(data_results(),
                                         #                     input$pres,
                                         #                     kmeans = TRUE,
                                         input$k1,
                                         ht_title,
                                         colour_functions_input()[[1]],
                                         colour_functions_input()[[2]],
                                         font_size=6,
                                         cluster_cols=input$cluster_cols,
                                         comp_list=comp(),
                                         data_results1=data_results1(),
                                         saintlist=saintlist()[[2]])
        }
        
        if (input$data_mode == TRUE | input$comp_type == "2 Rounds of comparisons"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          ht_title<-paste(input$data_type, "<br>", "Comparison 1+2", "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "Comp 2 FC=", input$lfc2, p_txt, input$p2, "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<- plot_complex_heatmap_2xcomp(data_results(),
                                           #                     input$pres,
                                           #                     kmeans = TRUE,
                                           input$k1,
                                           ht_title,
                                           colour_functions_input()[[1]],
                                           colour_functions_input()[[2]],
                                           font_size=6,
                                           cluster_cols=input$cluster_cols,
                                           comp_list=comp2(),
                                           data_results1=data_results1(),
                                           data_results2=data_results2())
        }
        
        if (input$data_mode == TRUE){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          if(input$data_type2 == "FragPipe"){data_type2 <- paste(input$data_type2, " ", input$exp2 , sep="")}else{data_type2<-input$data_type2}
          ht_title<-paste(data_type, " and ",  data_type2, "<br>", "Dataset 1+2", "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "Comp 2 FC=", input$lfc2, p_txt, input$p2, "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ", sep="")
          
          ht<- plot_complex_heatmap_2xcomp(data_results(),
                                           #                     input$pres,
                                           #                     kmeans = TRUE,
                                           input$k1,
                                           ht_title,
                                           colour_functions_input()[[1]],
                                           colour_functions_input()[[2]],
                                           font_size=6,
                                           cluster_cols=input$cluster_cols,
                                           comp_list=comp(),
                                           data_results1=data_results1(),
                                           data_results2=data_results2())
        }
        
        return(ht)
      })
    })
    
    comp_heatmap_input_c1 <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$data_mode == TRUE | input$comp_type == "2 Rounds of comparisons"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          ht_title<-paste(data_type, "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>",  "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c1," Clusters ", sep="")
          ht<-plot_complex_heatmap_2xcomp(data_results1(),
                                          #                     input$pres,
                                          #                     kmeans = TRUE,
                                          input$k1_c1,
                                          ht_title,
                                          colour_functions_input()[[1]],
                                          colour_functions_input()[[2]],
                                          font_size=8,
                                          cluster_cols=input$cluster_cols,
                                          comp_list=comp(),
                                          data_results1=data_results1(),
                                          data_results2=data_results2())
        }
        
        if (input$comp_type == "Use SAINTexpress"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          ht_title<-paste(data_type, "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c1," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_BFDR(data_results1(),
                                        #                     input$pres,
                                        #                     kmeans = TRUE,
                                        input$k1_c1,
                                        ht_title,
                                        colour_functions_input()[[1]],
                                        colour_functions_input()[[2]],
                                        font_size=8,
                                        cluster_cols=input$cluster_cols,
                                        comp_list=comp(),
                                        data_results1=data_results1(),
                                        saintlist=saintlist()[[2]])
        }
        return(ht)
      })
    })
    
    comp_heatmap_input_c1_pdf <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$data_mode == TRUE | input$comp_type == "2 Rounds of comparisons"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          ht_title<-paste(data_type, "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>",  "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c1," Clusters ", sep="")
          ht<-plot_complex_heatmap_2xcomp(data_results1(),
                                          #                     input$pres,
                                          #                     kmeans = TRUE,
                                          input$k1_c1,
                                          ht_title,
                                          colour_functions_input()[[1]],
                                          colour_functions_input()[[2]],
                                          font_size=6,
                                          cluster_cols=input$cluster_cols,
                                          comp_list=comp(),
                                          data_results1=data_results1(),
                                          data_results2=data_results2())
        }
        
        if (input$comp_type == "Use SAINTexpress"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          ht_title<-paste(data_type, "<br>", "Comp 1 FC=", input$lfc, p_txt, input$p, "<br>", "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c1," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_BFDR(data_results1(),
                                        #                     input$pres,
                                        #                     kmeans = TRUE,
                                        input$k1_c1,
                                        ht_title,
                                        colour_functions_input()[[1]],
                                        colour_functions_input()[[2]],
                                        font_size=6,
                                        cluster_cols=input$cluster_cols,
                                        comp_list=comp(),
                                        data_results1=data_results1(),
                                        saintlist=saintlist()[[2]])
        }
      })
    })
    
    comp_heatmap_input_c2 <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$comp_type == "Use SAINTexpress"){
          if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          ht_title<-paste(data_type, "<br>", "SAINTexpress List", "<br>", "BFDR_Cutoff=", input$bfdr,  "<br>", "BFDR_Freq=", input$bfdr_freq,  "<br>", "(", dep_num_bfdr(), " Proteins)", "<br>",input$k1_c2," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_BFDR(data_results4(),
                                        #                     input$pres,
                                        #                     kmeans = TRUE,
                                        input$k1_c2,
                                        ht_title,
                                        colour_functions_input()[[1]],
                                        colour_functions_input()[[2]],
                                        font_size=8,
                                        cluster_cols=input$cluster_cols,
                                        comp_list=comp(),
                                        data_results1=data_results1(),
                                        saintlist=saintlist()[[2]])
        }
        else{
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          if(input$data_mode == FALSE){
            if(input$data_type == "FragPipe"){data_type <- paste(input$data_type, " ", input$exp , sep="")}else{data_type<-input$data_type}
          }else{if(input$data_type2 == "FragPipe"){data_type <- paste(input$data_type2, " ", input$exp2 , sep="")}else{data_type<-input$data_type2}}
          ht_title<-paste(data_type, "<br>", "Comp 2 FC=", input$lfc2, p_txt, input$p2, "<br>",  "(", dep_num_c2(), " Proteins)", "<br>",input$k1_c2," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_2xcomp(data_results2(),
                                          #                     input$pres,
                                          #                     kmeans = TRUE,
                                          input$k1_c2,
                                          ht_title,
                                          colour_functions_input()[[1]],
                                          colour_functions_input()[[2]],
                                          font_size=8,
                                          cluster_cols=input$cluster_cols,
                                          comp_list=comp2(),
                                          data_results1=data_results1(),
                                          data_results2=data_results2())
        }
        return(ht)
      })
    })
    
    comp_heatmap_input_c2_pdf <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        
        if (input$data_mode == FALSE && input$comp_type == "2 Rounds of comparisons"){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          ht_title<-paste(input$data_type, "<br>", "Comp 2 FC=", input$lfc2, p_txt, input$p2, "<br>",  "(", dep_num_c2(), " Proteins)", "<br>",input$k1_c2," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_2xcomp(data_results2(),
                                          #                     input$pres,
                                          #                     kmeans = TRUE,
                                          input$k1_c2,
                                          ht_title,
                                          colour_functions_input()[[1]],
                                          colour_functions_input()[[2]],
                                          font_size=6,
                                          cluster_cols=input$cluster_cols,
                                          comp_list=comp2(),
                                          data_results1=data_results1(),
                                          data_results2=data_results2())
        }
        
        if (input$data_mode == FALSE && input$comp_type == "Use SAINTexpress"){
          ht_title<-paste(input$data_type, "<br>", "SAINTexpress List", "<br>", "BFDR_Cutoff=", input$bfdr,  "<br>", "BFDR_Freq=", input$bfdr_freq,  "<br>", "(", dep_num_bfdr(), " Proteins)", "<br>",input$k1_c2," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_BFDR(data_results4(),
                                        #                     input$pres,
                                        #                     kmeans = TRUE,
                                        input$k1_c2,
                                        ht_title,
                                        colour_functions_input()[[1]],
                                        colour_functions_input()[[2]],
                                        font_size=6,
                                        cluster_cols=input$cluster_cols,
                                        comp_list=comp(),
                                        data_results1=data_results1(),
                                        saintlist=saintlist()[[2]])
        }
        
        if (input$data_mode == TRUE){
          if(input$pval_use == FALSE){p_txt<-" p.adj="}else{p_txt<-" p.Val="}
          ht_title<-paste(input$data_type2, "<br>", "Comp 2 FC=", input$lfc2, p_txt, input$p2, "<br>",  "(", dep_num_c2(), " Proteins)", "<br>",input$k1_c2," Clusters ", sep="")
          
          ht<-plot_complex_heatmap_2xcomp(data_results2(),
                                          #                     input$pres,
                                          #                     kmeans = TRUE,
                                          input$k1_c2,
                                          ht_title,
                                          colour_functions_input()[[1]],
                                          colour_functions_input()[[2]],
                                          font_size=6,
                                          cluster_cols=input$cluster_cols,
                                          comp_list=comp2(),
                                          data_results1=data_results1(),
                                          data_results2=data_results2())
        }
        
        return(ht)
      })
    })
    
    #Venn Diagram
    venn_input<- reactive({
      if (input$comp_type == "2 Rounds of comparisons"){
        comparison_1 <- data_results1()$name
        comparison_2 <- data_results2()$name
        ven <- venndetail(list(Round1 = comparison_1, Round2 = comparison_2))
        plot_venn<-plot(ven)
        venn_list<-result(ven)
        ven_all<-list(plot_venn, venn_list, ven)
        
      }
      if (input$comp_type == "Use SAINTexpress"){
        comparison_1 <- data_results1()$name
        comparison_2 <- data_results4()$name
        ven <- venndetail(list(Round1 = comparison_1, SAINTexpress = comparison_2))
        plot_venn<-plot(ven)
        venn_list<-result(ven)
        ven_all<-list(plot_venn, venn_list, ven)
        
      }
      if (input$data_mode == TRUE){
        
        comparison_1 <- data_results1()$name
        comparison_2 <- data_results2()$name
        ven <- venndetail(list(Dataset1 = comparison_1, Dataset2 = comparison_2))
        plot_venn<-plot(ven)
        venn_list<-result(ven)
        ven_all<-list(plot_venn, venn_list, ven)
        
      }
      
      return(ven_all)  
    })
    
    
    #Return my_mat from heatmap (object [[2]])
    data_results_mat<- reactive({
      data.frame(comp_heatmap_input()[[2]])
    })
    
    #    ht<-comp_heatmap_input()[[1]]
    #    my_mat<-comp_heatmap_input()[[2]]
    #    dfdata_sig2<-comp_heatmap_input()[[3]]
    
    cludf_check_input<- reactive({
      cludf<- cluster_check(
        comp_heatmap_input()[[1]],
        comp_heatmap_input()[[2]],
        input$k1
        #        colour_functions_input()[[1]]
      )
      data.frame(cludf)
    })
    
    #Heatmap Clusters Input
    cluster_heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        #Check for single protein clusters
        validate(need(try(cludf_check_input() != "try-error"), "One of your Clusters contains a single protein, 
          please alter the Kmeans number")
        )
        
        if (input$data_mode == TRUE | input$comp_type == "2 Rounds of comparisons"){
          ht<-plot_heatmap_clusters_2xcomp(comp_heatmap_input()[[1]],
                                           comp_heatmap_input()[[2]],
                                           comp_heatmap_input()[[3]],
                                           #                     input$pres,
                                           #                     kmeans = TRUE,
                                           input$k1,
                                           input$lfc,
                                           input$p,
                                           colour_functions_input()[[1]],
                                           colour_functions_input()[[2]],
                                           font_size=8,
                                           data_unique_defn = data_defn()
          )
        }
        
        if (input$data_mode == FALSE && input$comp_type == "1 Round of comparisons"){
          ht<-plot_heatmap_clusters(comp_heatmap_input()[[1]],
                                    comp_heatmap_input()[[2]],
                                    comp_heatmap_input()[[3]],
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1,
                                    input$lfc,
                                    input$p,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=8,
                                    data_unique_defn = data_defn()
          )
        }
        
        if (input$data_mode == FALSE && input$comp_type == "Use SAINTexpress"){
          ht<-plot_heatmap_clusters_2xBFDR(comp_heatmap_input()[[1]],
                                           comp_heatmap_input()[[2]],
                                           comp_heatmap_input()[[3]],
                                           #                     input$pres,
                                           #                     kmeans = TRUE,
                                           input$k1,
                                           input$lfc,
                                           input$p,
                                           colour_functions_input()[[1]],
                                           colour_functions_input()[[2]],
                                           font_size=8,
                                           data_unique_defn = data_defn()
          )
        }
        return(ht) 
      })
    }) 
    
    
    
    #Enhanced Volcano Input 1
    enhanced_volcano_1_input <- reactive({
      if (input$comp_type == "2 Rounds of comparisons"){
        lfc<-input$lfc2
        p<-input$p2
        if(input$enh_vol_data==FALSE){vol_data<-get_results(dep2())}
        else{vol_data<-data_results()}
      }else{lfc<-input$lfc
      p<-input$p
      if(input$enh_vol_data==FALSE){vol_data<-get_results(dep())}
      else{vol_data<-data_results()}
      }
      #Check for single protein clusters
      validate(need(try(cludf_check_input() != "try-error"), "One of your Clusters contains a single protein, 
          please alter the Kmeans number")
      )
      input_volcano_cntrst_enh1 <- req(input$volcano_cntrst_enh1)
      plot_enhanced_volcano (vol_data,
                             #get_results(dep()),
                             cludf_check_input(),
                             comp_pick = input_volcano_cntrst_enh1,
                             FC=lfc,
                             pVal=p,
                             padj_check=input$p_adj_enh1,
                             colour_functions_input()[[2]],
                             colour_functions_input()[[3]],
                             connect=input$enh_vol_connect,
                             selected=input$enh_vol_names, 
                             selected_names = cludf_check_input()[input$tablecludf_rows_selected,1]
      )
    })
    
    #Enhanced Volcano Input 2
    enhanced_volcano_2_input <- reactive({
      if (input$comp_type == "2 Rounds of comparisons"){
        lfc<-input$lfc2
        p<-input$p2
        if(input$enh_vol_data2==FALSE){vol_data<-get_results(dep2())}
        else{vol_data<-data_results()}
      }else{lfc<-input$lfc
      p<-input$p
      if(input$enh_vol_data2==FALSE){vol_data<-get_results(dep())}
      else{vol_data<-data_results()}
      }
      #Check for single protein clusters
      validate(need(try(cludf_check_input() != "try-error"), "One of your Clusters contains a single protein, 
          please alter the Kmeans number")
      )
      input_volcano_cntrst_enh2 <- req(input$volcano_cntrst_enh2)
      plot_enhanced_volcano (vol_data,
                             cludf_check_input(),
                             comp_pick = input_volcano_cntrst_enh2,
                             FC=lfc,
                             pVal=p,
                             padj_check=input$p_adj_enh2,
                             colour_functions_input()[[2]],
                             colour_functions_input()[[3]],
                             connect=input$enh_vol_connect2,
                             selected=input$enh_vol_names2, 
                             selected_names = cludf_check_input()[input$tablecludf_rows_selected,1]
                             
                             
      )
    })
    
    #Plot_Box_Plots
    box_plot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        #Check for single protein clusters
        validate(need(try(cludf_check_input() != "try-error"), "One of your Clusters contains a single protein, 
          please alter the Kmeans number")
        )
        
        plot_box_plot (norm_data(),
                       imp_data(),
                       clu_df=cludf_check_input(),
                       cluster_no=input$cluster_no_box,
                       input$palette_choice
        )
      })
    })
    
    #Plot_Box_Plots_reports
    box_plot_report_input <- reactive({
      withProgress(message = 'Generating Box Plots', value = 0.66, {
        #Check for single protein clusters
        validate(need(try(cludf_check_input() != "try-error"), "One of your Clusters contains a single protein, 
          please alter the Kmeans number")
        )
        
        box_plot_list<-list()
        for(i in 1:input$k1){
          box_plot<- plot_box_plot (
            norm_data(),
            imp_data(),
            clu_df=cludf_check_input(),
            cluster_no=i,
            input$palette_choice
          )
          box_plot_list[[i]]<-box_plot  
        }
        box_plot_list
      })
    })
    
    #Make entrezlist
    entrezlist_input <- reactive({
      validate(
        need(try(cludf_check_input() != "try-error"), "One of your Clusters contains a single protein, please alter the Kmeans number")
      )
      
      messages <- character()
      
      result <- withCallingHandlers(
        get_entrezlist(cludf_check_input(), input$species),
        message = function(m) {
          messages <<- c(messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      
      conversion_messages(messages)  # Store messages for later display
      
      msgs <- conversion_messages()
      error_msgs <- msgs[grepl("None of the keys entered are valid", msgs)]
      
      if (length(error_msgs) > 0) {
        print("Error in asigning GO terms, please check the correct G0 species is selected")
        showNotification("Error in asigning GO terms, please check the correct G0 species is selected", type = "error", duration = 20)
      }
      
      
      return(result)
    })
    
    GOMF_cnet_input <- reactive({
      if (is.null(GO_MF_report_input())) {
        # Show a text message while waiting for the plot data
        plot.new()  # Create an empty plot
        text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
      } else {
        GOtitle <- paste(" Molecular Function up to level ", input$gomf_level, sep="")

        cnet_MF<-cnetplot_custom(cnet_data = GO_MF_report_input(), 
                        mypal = colour_functions_input()[[2]], 
                        GOtitle = GOtitle, 
                        showCategory = input$gomf_level,
                        node_scale = input$gomf_node)
        
        print(cnet_MF)
      }
    })
    
    GOCC_cnet_input <- reactive({
      if (is.null(GO_CC_report_input())) {
        # Show a text message while waiting for the plot data
        plot.new()  # Create an empty plot
        text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
      } else {
        GOtitle <- paste(" Cellular Component up to level ", input$gocc_level, sep="")
        
        cnet_CC<-cnetplot_custom(cnet_data = GO_CC_report_input(), 
                                 mypal = colour_functions_input()[[2]], 
                                 GOtitle = GOtitle, 
                                 showCategory = input$gocc_level,
                                 node_scale = input$gocc_node)
        
        print(cnet_CC)
      }
    })
    
    GOBP_cnet_input <- reactive({
      if (is.null(GO_BP_report_input())) {
        # Show a text message while waiting for the plot data
        plot.new()  # Create an empty plot
        text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
      } else {
        GOtitle <- paste(" Biological Process up to level ", input$gobp_level, sep="")
        
        cnet_BP<-cnetplot_custom(cnet_data = GO_BP_report_input(), 
                                 mypal = colour_functions_input()[[2]], 
                                 GOtitle = GOtitle, 
                                 showCategory = input$gobp_level,
                                 node_scale = input$gobp_node)
        
        print(cnet_BP)
      }
    })
    
    #BioGrid Plot
    BG_cnet_input <- reactive({
      if (is.null(BG_report_input())) {
        # Show a text message while waiting for the plot data
        plot.new()  # Create an empty plot
        if(nrow(cludf_check_input()) > 160){text(0.5, 0.5, "Lower the number of signifcant protein to below 160", cex = 1.5)}else
          text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
      } else {
        clu_df<-cludf_check_input()
        if(input$use_human == TRUE){
          clu_df$GeneID <- toupper(clu_df$GeneID)}# Convert GeneID column to uppercase
        biogrid_network_input <- biogrid_network(biogrid_list = BG_report_input(),
                                                 cludf =clu_df,
                                                 mypal = colour_functions_input()[[2]],
                                                 inter_count = input$BG_level,
                                                 node_scale = input$BG_node)
        
        print (biogrid_network_input)
      }
    })
    
    

    #### Output functions  --------------------------------------------------
    output$table <- DT::renderDataTable({
      if (input$top_table_opt == TRUE){
        if (input$comp_type == "2 Rounds of comparisons"){get_results(dep2())}else{get_results(dep())}
      }else{
        data_results() #should be default
        #pca_input()[[2]]
        #heatmap_input()[[2]]
        #norm_data()
        #data_other_choice() #to test data input
        #data2()
      }
    }, options = list(pageLength = 25, scrollX = T),
    selection = list(selected = c(1)))
    
    #Output of table on 2nd tab
    output$table2 <- DT::renderDataTable({
      data_results_mat()
    }, options = list(pageLength = 25, scrollX = T),
    )
    
    #Output of table on Cluster tab
    output$table3 <- DT::renderDataTable({
      validate(need(!(input$cluster_no > input$k1), "The chosen Cluster No. is higher than the number of clusters"))
      cluster_heatmap_input()[[3]][[input$cluster_no]]
    }, options = list(pageLength = 25, scrollX = T),
    )
    
    #settings table
    output$table4 <- DT::renderDataTable({
      settings_input()
    })
    
    #Exp Design Table 
    output$exp_des <- DT::renderDataTable({
      expdesign()
    })
    
    #Comparisons Table
    output$table_comp <- DT::renderDataTable({
      comp_table()
    }) 
    
    #Test_table
    output$test <- DT::renderDataTable({
      df_long_input()
      #norm_data()
    })
    
    #cludf table
    output$tablecludf <- DT::renderDataTable({
      cludf_check_input()
    }) 
    
    
    #Abundance Table
    output$table_abund <- DT::renderDataTable({
      data_defn_abundance()
    })
    
    output$table_abund2 <- DT::renderDataTable({
      data_defn_abundance2()
    })
    
    
    #Saint Outputs
    output$tablesaint_output1 <- DT::renderDataTable({
      saint_output()[[1]]
    })
    output$tablesaint_output2 <- DT::renderDataTable({
      # validate(need(saint_output()[[2]] != "", " "))
      
      saint_output()[[2]]
    }) 
    output$tablesaint_output3 <- DT::renderDataTable({
      #validate(need(saint_output()[[3]] != "", " "))
      
      saint_output()[[3]]
    }) 
    output$tablesaint_output4 <- DT::renderDataTable({
      #validate(need(saint_output()[[4]] != "", " "))
      
      saint_output()[[4]]
    }) 
    
    
    #Saint table
    output$tablesaint <- DT::renderDataTable({
      saintlist()[[2]]
    }) 
    
    #Selcted Protein plots
    output$selected_plot <- renderPlot({
      if (input$protein_plot == "contrast" | input$protein_plot == "centered"){
        selected_plot_input()
      }else{
        if (input$protein_plot=="normalised"){choice=1}else{choice=2}
        selected_plot_input()[[choice]]
      }
    })
    
    output$selected_plot2 <- renderPlot({
      if (input$protein_plot2 == "contrast" | input$protein_plot2 == "centered"){
        selected_plot_input2()
      }else{
        if (input$protein_plot2=="normalised"){choice=1}else{choice=2}
        selected_plot_input2()[[choice]]
      }
    })
    
    
    output$heatmap <- renderPlot({
      heatmap_input()[[1]]
    })
    
    output$heatmap2 <- renderPlot({
      heatmap_input2()[[1]]
    })
    
    
    output$volcano <- renderPlot({
      volcano_input()
    })
    
    output$volcano_c2 <- renderPlot({
      volcano_c2_input()
    })
    
    #signif plot output
    output$signif_plot <- renderPlot({
      signif_plot_input()
    })
    
    output$signif_plot2 <- renderPlot({
      signif_plot_input2()
    })
    
    output$norm <- renderPlot({
      norm_input()
    })
    
    output$norm2 <- renderPlot({
      norm_input2()
    })
    
    
    output$missval <- renderPlot({
      missval_input()
    })
    
    output$detect <- renderPlot({
      detect_input()
    })
    
    output$imputation <- renderPlot({
      imputation_input()
    })
    
    output$imputation2 <- renderPlot({
      imputation_input2()
    })
    
    output$numbers <- renderPlot({
      numbers_input()
    })
    
    output$numbers2 <- renderPlot({
      numbers_input2()
    })
    
    output$coverage <- renderPlot({
      coverage_input()
    })
    
    #PCA output
    output$pca <- renderPlot({
      pca_input()[[1]]
    })
    
    #PCA output scree
    output$pca_scree <- renderPlot({
      pca_input()[[3]]
    })
    
    
    #PCA output (Comp2)
    output$pca_c2 <- renderPlot({
      pca_c2_input()[[1]]
    })
    
    #PCA output scree
    output$pca_scree2 <- renderPlot({
      pca_c2_input()[[3]]
    })
    
    
    #Correlation output
    output$cor <- renderPlot({
      cor_input()
    })
    
    #Correlation output (Comp2)
    output$cor_c2 <- renderPlot({
      cor_c2_input()
    })
    
    #Frequency plot output
    output$cond <- renderPlot({
      cond_input()
    })
    
    #Frequency plot output (Comp2)
    output$cond_c2 <- renderPlot({
      cond_c2_input()
    })
    
    #swatch
    output$swatch <- renderPlot({
      scales::show_col(colour_functions_input()[[2]])
    })
    
    #complex heatmap output
    output$comp_heatmap <- renderPlot({
      comp_heatmap_input()[[1]]
    })
    
    output$comp_heatmap_c1 <- renderPlot({
      comp_heatmap_input_c1()[[1]]
    })
    
    output$comp_heatmap_c2 <- renderPlot({
      comp_heatmap_input_c2()[[1]]
    })
    
    #Venn Diagram
    output$venn <- renderPlot({
      plot(venn_input()[[3]])
      #venn_input()[[1]]
    })
    
    #Venn table
    output$venn_table <- DT::renderDataTable({
      DT::datatable(venn_input()[[2]], options = list(lengthMenu = c(5, 30, 50), pageLength = 7))
    })
    
    #cluster heatmap output
    output$cluster_heatmap <- renderPlot({
      cluster_heatmap_input()[[1]][[input$cluster_no]]
    })
    
    #Box Plots
    output$box_plot <- renderPlot({
      validate(need(!(input$cluster_no_box > input$k1), "The chosen Cluster No. is higher than the number of clusters"))
      if (input$pres_box=="normalised"){choice=1}else{choice=2}
      box_plot_input()[[choice]]
    })
    
    #Enhanced Volcano plot 1
    output$enh_volcano_1 <- renderPlot({
      enhanced_volcano_1_input()
    })
    
    #Enhanced Volcano plot 2
    output$enh_volcano_2 <- renderPlot({
      enhanced_volcano_2_input()
    })
    
    #Report
    output$downloadreport <- renderUI({
      if(!is.null(dep())){
        #Check for single protein clusters
        #        validate(need(try(cludf_check_input() != "try-error"), "One of your Clusters contains a single protein, please alter the Kmeans number")
        #       )
        if(input$report_type=="Short"){
        }
        
        if(input$report_type=="with Box Plots"){
          box<-box_plot_report_input()
          
        }
        
        if(input$report_type=="Long"){
          box<-box_plot_report_input()
        }
        downloadButton('downloadReport', 'Download Report', )
        
      }
    })
    
    
    
    ####Observe####
    # Auto-trigger only when 'long' is selected
    observeEvent(input$report_type, {
      if (input$report_type == "Long") {
        #trigger_GO_MF <- reactiveVal(TRUE)
        trigger_GO_MF(TRUE)
        trigger_GO_CC(TRUE)
        trigger_GO_BP(TRUE)
        trigger_BG(TRUE)

      }
    })
    
    observeEvent(input$analyze_MF, {
      trigger_GO_MF(TRUE)
    })
    
    observeEvent(input$analyze_CC, {
      trigger_GO_CC(TRUE)
    })
    
    observeEvent(input$analyze_BP, {
      trigger_GO_BP(TRUE)
    })
    
    observeEvent(input$analyze_BG, {
      trigger_BG(TRUE)
    })
    
    
    observe({
      output$plot <- renderUI({
        ht_height<-((input$size*12)+270)
        if (ht_height < 420){ht_height<-ht_height+80}
        plotOutput("heatmap", height = ht_height)
        #plotOutput("heatmap", height = (100 * as.numeric(input$size)))
      })
    })
    
    observe({
      output$plot2 <- renderUI({
        ht_height<-((input$size2*12)+270)
        if (ht_height < 420){ht_height<-ht_height+80}
        plotOutput("heatmap2", height = ht_height)
        #plotOutput("heatmap", height = (100 * as.numeric(input$size)))
      })
    })
    
    
    observe({
      output$plot_comp <- renderUI({
        ht_height<-((dep_num()*12)+270+input$heatmap_pad)
        if (ht_height < 420){ht_height<-ht_height+80}
        plotOutput("comp_heatmap", height = ht_height)
      })
    })
    
    observe({
      output$plot_comp_c1 <- renderUI({
        plotOutput("comp_heatmap_c1", height = ((dep_num_c1()*12)+270))
      })
    })
    
    observe({
      output$plot_comp_c2 <- renderUI({
        if (input$data_mode == TRUE | input$comp_type == "2 Rounds of comparisons"){
          plotOutput("comp_heatmap_c2", height = ((dep_num_c2()*12)+270))
        }
        else if (input$comp_type == "Use SAINTexpress"){
          plotOutput("comp_heatmap_c2", height = ((dep_num_bfdr()*12)+270))
        }
      })
    })
    
    observe({
      output$plot_clust <- renderUI({
        validate(need(!(input$cluster_no > input$k1), "The chosen Cluster No. is higher than the number of clusters"))
        plotOutput("cluster_heatmap", height = (cluster_heatmap_input()[[2]][[input$cluster_no]]*12)+270+input$cluster_pad)
      })
    })
    
    
    #Observe GO Plots
    #Observe MF
    observeEvent(trigger_GO_MF(),{
      if (trigger_GO_MF() == TRUE) {
        sendSweetAlert(session, title = "Plotting Molecular Function", html = TRUE, type = NULL, text = HTML('<img src="app_icon_anim.gif" width="100"><br><br><strong>This may take a while...</strong>'), timer = 10000)
        
        # ✅ Correct way to get reactive values
        entrez_list <- isolate(entrezlist_input())  # Extract reactive value
        species <- isolate(input$species)  # Capture input value
        
        future_promise({
          # Validate input inside async process
          validate(need(try(entrez_list != "try-error"), 
                        "None of the Gene_IDs matched to the species database, make sure you have selected the correct GO species"))
          
          plot_GO_MF(entrez_list, species)  # Call function with extracted values
        }, seed = 123) %...>% GO_MF_report_input %...!% (function(e) {
          GO_MF_report_input(NULL)
          showNotification(paste("Error:", e$message), type = "error")
        })
        
        GOMF_plot_input <- reactive({
          if (is.null(GO_MF_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
                 GOMF_dotplot<- dotplot(GO_MF_report_input(), showCategory = input$gomf_level)
          }
        })


        #Plot output
        output$GOMF_plot <- renderPlot({
          if (is.null(GO_MF_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
            GOMFtitle<-paste("Molecular Function up to level ", input$gomf_level, sep="")
            print(GOMF_plot_input() + ggtitle(GOMFtitle))
          }
        })
        
        output$GOMF_cnet <- renderPlot({
          if (is.null(GO_MF_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
            print(GOMF_cnet_input())
          }
        })
        
        #PDF Output
        output$download_MF_plot <- downloadHandler(
          filename = function() { paste("Molecular_Function_Plot_level_",input$gomf_level, ".pdf", sep = "") },
          content = function(file) {
            cairo_pdf(file, height = input$pdfheight_MF, width = input$pdfwidth_MF)
            GOMFtitle<-paste("Molecular Function up to level ", input$gomf_level, sep="")
            print(GOMF_dotplot<- dotplot(GO_MF_report_input(), showCategory = input$gomf_level) + ggtitle(GOMFtitle))
            dev.off()
          }
        )
        output$download_MF_cnet <- downloadHandler(
          filename = function() { paste("Molecular_Function_Map_level_",input$gomf_level, ".pdf", sep = "") },
          content = function(file) {
            cairo_pdf(file, height = input$pdfheight_MF, width = input$pdfwidth_MF)
            print(GOMF_cnet_input())
            dev.off()
          }
        )
        output$download_MF_all <- downloadHandler(
          filename = function() { paste("Molecular_Function_Data.txt", sep = "") },
          content = function(file) {
            write.table(GO_MF_report_input(),
                        file,
                        col.names = TRUE,
                        row.names = FALSE,
                        sep ="\t") 
          }
        )
        
        trigger_GO_MF(FALSE)
      }
    })
    
    #Observe CC
    observeEvent(trigger_GO_CC(),{
      if (trigger_GO_CC() == TRUE) {
        sendSweetAlert(session, title = "Plotting Cellular Component", html = TRUE, type = NULL, text = HTML('<img src="app_icon_anim.gif" width="100"><br><br><strong>This may take a while...</strong>'), timer = 10000)
        
        # ✅ Correct way to get reactive values
        entrez_list <- isolate(entrezlist_input())  # Extract reactive value
        species <- isolate(input$species)  # Capture input value
        
        future_promise({
          # Validate input inside async process
          validate(need(try(entrez_list != "try-error"), 
                        "None of the Gene_IDs matched to the species database, make sure you have selected the correct GO species"))
          
          plot_GO_CC(entrez_list, species)  # Call function with extracted values
        }, seed = 123) %...>% GO_CC_report_input %...!% (function(e) {
          GO_CC_report_input(NULL)
          showNotification(paste("Error:", e$message), type = "error")
        })
        
        #Plot output
        output$GOCC_plot <- renderPlot({
          if (is.null(GO_CC_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
            GOCCtitle<-paste("Cellular Component up to level ", input$gocc_level, sep="")
            print(GOCC_dotplot<- dotplot(GO_CC_report_input(), showCategory = input$gocc_level) + ggtitle(GOCCtitle))
          }
        })
        output$GOCC_cnet <- renderPlot({
          if (is.null(GO_CC_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
            print(GOCC_cnet_input())
          }
        })
        #PDF Output
        output$download_CC_plot <- downloadHandler(
          filename = function() { paste("Cellular_Component_Plot_level_",input$gocc_level, ".pdf", sep = "") },
          content = function(file) {
            cairo_pdf(file, height = input$pdfheight_CC, width = input$pdfwidth_CC)
            GOCCtitle<-paste("Cellular Component up to level ", input$gocc_level, sep="")
            print(GOCC_dotplot<- dotplot(GO_CC_report_input(), showCategory = input$gocc_level) + ggtitle(GOCCtitle))
            dev.off()
          }
        )
        output$download_CC_cnet <- downloadHandler(
          filename = function() { paste("Cellular_Component_Map_level_",input$gocc_level, ".pdf", sep = "") },
          content = function(file) {
            cairo_pdf(file, height = input$pdfheight_CC, width = input$pdfwidth_CC)
            print(GOCC_cnet_input())
            dev.off()
          }
        )
        output$download_CC_all <- downloadHandler(
          filename = function() { paste("Cellular_Component_Data.txt", sep = "") },
          content = function(file) {
            write.table(GO_CC_report_input(),
                        file,
                        col.names = TRUE,
                        row.names = FALSE,
                        sep ="\t") 
          }
        )
        trigger_GO_CC(FALSE)
      }
    })
    
    #Observe BP
    observeEvent(trigger_GO_BP(),{
      if (trigger_GO_BP() == TRUE) {
        sendSweetAlert(session, title = "Plotting Biological Function", html = TRUE, type = NULL, text = HTML('<img src="app_icon_anim.gif" width="100"><br><br><strong>This may take a while...</strong>'), timer = 10000)
        
        # ✅ Correct way to get reactive values
        entrez_list <- isolate(entrezlist_input())  # Extract reactive value
        species <- isolate(input$species)  # Capture input value
        
        future_promise({
          # Validate input inside async process
          validate(need(try(entrez_list != "try-error"), 
                        "None of the Gene_IDs matched to the species database, make sure you have selected the correct GO species"))
          
          plot_GO_BP(entrez_list, species)  # Call function with extracted values
        }, seed = 123) %...>% GO_BP_report_input %...!% (function(e) {
          GO_BP_report_input(NULL)
          showNotification(paste("Error:", e$message), type = "error")
        })
        
        
        #Plot output
        output$GOBP_plot <- renderPlot({
          if (is.null(GO_BP_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
            GOBPtitle<-paste("Biological Process up to level ", input$gobp_level, sep="")
            print(GOBP_dotplot<- dotplot(GO_BP_report_input(), showCategory = input$gobp_level) + ggtitle(GOBPtitle))
          }
        })
        
        output$GOBP_cnet <- renderPlot({
          if (is.null(GO_BP_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
            print(GOBP_cnet_input())
          }
        })
        
        #PDF Output
        output$download_BP_plot <- downloadHandler(
          filename = function() { paste("Biological_Process_Plot_level_",input$gobp_level, ".pdf", sep = "") },
          content = function(file) {
            cairo_pdf(file, height = input$pdfheight_BP, width = input$pdfwidth_BP)
            GOBPtitle<-paste("Biological Process up to level ", input$gobp_level, sep="")
            print(GOBP_dotplot<- dotplot(GO_BP_report_input(), showCategory = input$gobp_level) + ggtitle(GOBPtitle))
            dev.off()
          }
        )
        
        output$download_BP_cnet <- downloadHandler(
          filename = function() { paste("Biological_Process_Map_level_",input$gobp_level, ".pdf", sep = "") },
          content = function(file) {
            cairo_pdf(file, height = input$pdfheight_BP, width = input$pdfwidth_BP)
            print(GOBP_cnet_input())
            dev.off()
          }
        )
        
        output$download_BP_all <- downloadHandler(
          filename = function() { paste("Biological_Process_Data.txt", sep = "") },
          content = function(file) {
            write.table(GO_BP_report_input(),
                        file,
                        col.names = TRUE,
                        row.names = FALSE,
                        sep ="\t") 
          }
        )
        trigger_GO_BP(FALSE)
      }
    })
    
    #Observe BioGrid
    observeEvent(trigger_BG(),{
      if (trigger_BG() == TRUE) {
        sendSweetAlert(session, title = "Plotting BioGrid Interactions", html = TRUE, type = NULL, text = HTML('<img src="app_icon_anim.gif" width="100"><br><br><strong>This may take a while...</strong>'), timer = 10000)
        
        # ✅ Correct way to get reactive values
        clu_df <- isolate(cludf_check_input())  # Extract reactive value
        #species <- 9606
        species_choice <- isolate(input$species)  # Capture input value
        if(species_choice =="Human"){species <- 9606}
        if(species_choice =="Mouse"){species <- 10090}
        if(species_choice =="Arabidopsis"){species <- 3702}
        if(species_choice =="Xenopus"){species <- 8355}
        
        if(input$use_human == TRUE){
          clu_df$GeneID <- toupper(clu_df$GeneID)# Convert GeneID column to uppercase
          species <- 9606 } 
        

        future_promise({
          validate(
            need(nrow(clu_df) < 160, 
                 "BioGRID searches are limited to a significant protein count below 160")
          )
          biogrid_query(cludf = clu_df, speciesID = species)
          #plot_GO_BP(entrez_list, species)  # Call function with extracted values
        }) %...>% BG_report_input %...!% (function(e) {
          BG_report_input(NULL)
          showNotification(paste("Error:", e$message), type = "error")
        })
        
        

        output$BG_cnet <- renderPlot({
          if (is.null(BG_report_input())) {
            # Show a text message while waiting for the plot data
            plot.new()  # Create an empty plot
            if(nrow(clu_df) > 160){text(0.5, 0.5, "Lower the number of signifcant protein to below 160", cex = 1.5)}else
            text(0.5, 0.5, "Generating plot. Please wait...", cex = 1.5)
          } else {
            print(BG_cnet_input())
          }
        })
        

        output$download_BG_cnet <- downloadHandler(
          filename = function() { paste("BioGrid_Map_level_",input$BG_level, ".pdf", sep = "") },
          content = function(file) {
            cairo_pdf(file, height = input$pdfheight_BG, width = input$pdfwidth_BG)
            print(BG_cnet_input())
            dev.off()
          }
        )
        
        output$download_BG_data <- downloadHandler(
          filename = function() { paste("BioGRid_Search_Data.txt", sep = "") },
          content = function(file) {
            write.table(BG_report_input(),
                        file,
                        col.names = TRUE,
                        row.names = FALSE,
                        sep ="\t") 
          }
        )
        trigger_BG(FALSE)
      }
    })
    
        
    #### Download objects and functions  ------------------------------------
    datasetInput <- reactive({
      switch(input$dataset,
             "results" = get_results(dep()),
             "significant_proteins" = get_results(dep()) %>%
               filter(significant) %>%
               select(-significant),
             "displayed_subset" = data_results(), #%>%
             #filter(significant) %>%
             #select(-significant),
             "full_dataset" = get_df_wide(dep()),
             "Cluster_list" = cludf_check_input(),
             "normalised_values"= norm_data_table(),
             "imputed_values"= imp_data_table(),
             "results_long"=get_df_long(dep()),
             "results_long_signif" = get_df_long(dep())%>% filter(significant) %>% select(-significant)
             
      )
    }
    )
    
    output$downloadData <- downloadHandler(
      filename = function() { paste(input$dataset, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
      
    )
    
    #Abundance table
    output$downloadabundancetable <- downloadHandler(
      filename = function() { paste("Top_", input$size, "_Abundant_Proteins_", input$abund_choice, ".txt", sep = "") },
      content = function(file) {
        write.table(data_defn_abundance(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    #Abundance table2
    output$downloadabundancetable2 <- downloadHandler(
      filename = function() { paste("Top_", input$size, "_Abundant_Proteins_dataset2", input$abund_choice, ".txt", sep = "") },
      content = function(file) {
        write.table(data_defn_abundance2(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    
    #Cluster table
    output$downloadclustertable <- downloadHandler(
      filename = function() { paste("Cluster ", input$cluster_no, ".txt", sep = "") },
      content = function(file) {
        write.table(cluster_heatmap_input()[[3]][[input$cluster_no]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    #cludf table
    output$downloadclusterlist <- downloadHandler(
      filename = function() { paste("Cluster_list.txt", sep = "") },
      content = function(file) {
        write.table(cludf_check_input(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    
    #SAINTtables table
    output$download_saint1 <- downloadHandler(
      filename = "SAINTexpress_1.txt",
      content = function(file) {
        write.table(saint_output()[[1]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    output$download_saint2 <- downloadHandler(
      filename = "SAINTexpress_2.txt",
      content = function(file) {
        write.table(saint_output()[[2]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    output$download_saint3 <- downloadHandler(
      filename = "SAINTexpress_3.txt",
      content = function(file) {
        write.table(saint_output()[[3]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    output$download_saint4 <- downloadHandler(
      filename = "SAINTexpress_4.txt",
      content = function(file) {
        write.table(saint_output()[[4]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    #BFDR Table
    output$download_tablesaint <- downloadHandler(
      filename = "BFDR_filtered_list.txt",
      content = function(file) {
        write.table(saintlist()[[2]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    #Variable table
    output$downloadsettings <- downloadHandler(
      filename = "Settings.txt",
      content = function(file) {
        write.table(settings_input(),
                    file,
                    col.names = TRUE,
                    row.names = TRUE,
                    sep ="\t") 
      }
    )
    
    #Venn Table
    output$download_venntable <- downloadHandler(
      filename = "Venn_list.txt",
      content = function(file) {
        write.table(venn_input()[[2]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    #Venn Diagram
    output$download_venn <- downloadHandler(
      filename = 'Venn.pdf',
      content = function(file) {
        cairo_pdf(file)
        plot(venn_input()[[3]])
        dev.off()
      }
    )
    
    
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0("Barplot_", data_results()[input$table_rows_selected,1], ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        if (input$protein_plot == "contrast" | input$protein_plot == "centered"){
          print(selected_plot_input())
        }else{
          if (input$protein_plot=="normalised"){choice=1}else{choice=2}
          print(selected_plot_input()[[choice]])
        }
        dev.off()
      }
    )
    
    output$downloadPlot2 <- downloadHandler(
      filename = function() {
        paste0("Barplot_data2_", data_results()[input$table_rows_selected,1], ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        if (input$protein_plot2 == "contrast" | input$protein_plot2 == "centered"){
          print(selected_plot_input2())
        }else{
          if (input$protein_plot2=="normalised"){choice=1}else{choice=2}
          print(selected_plot_input2()[[choice]])
        }
        dev.off()
      }
    )
    
    output$downloadHeatmap <- downloadHandler(
      filename = 'Heatmap.pdf',
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(heatmap_input())
        dev.off()
      }
    )
    
    output$downloadHeatmap2 <- downloadHandler(
      filename = 'Heatmap_dataset2.pdf',
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(heatmap_input2())
        dev.off()
      }
    )
    
    
    output$downloadVolcano <- downloadHandler(
      filename = function() {
        paste0("Volcano_", input$volcano_cntrst, ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(volcano_input())
        dev.off()
      }
    )
    
    output$downloadVolcano_c2 <- downloadHandler(
      filename = function() {
        paste0("Volcano_", input$volcano_cntrst_c2, ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(volcano_c2_input())
        dev.off()
      }
    )
    
    
    output$downloadsignif_plot <- downloadHandler(
      if (input$pval_use == FALSE){pval_txt<-"padj"}else{pval_txt<-"pVal"},
      filename = function() {paste0("Signif_plot_", pval_txt, ".pdf")},
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(signif_plot_input())
        dev.off()
      }
    )
    
    output$downloadsignif_plot2 <- downloadHandler(
      if (input$pval_use == FALSE){pval_txt<-"padj"}else{pval_txt<-"pVal"},
      filename = function() {paste0("Signif_plot2_", pval_txt, ".pdf")},
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(signif_plot_input2())
        dev.off()
      }
    )
    
    
    output$downloadNorm <- downloadHandler(
      filename = "normalization.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(norm_input())
        dev.off()
      }
    )
    
    output$downloadNorm2 <- downloadHandler(
      filename = "normalization (dataset 2).pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(norm_input2())
        dev.off()
      }
    )
    
    
    output$downloadMissval <- downloadHandler(
      filename = "missing_values_heatmap.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(missval_input())
        dev.off()
      }
    )
    
    output$downloadDetect <- downloadHandler(
      filename = "missing_values_quant.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        gridExtra::grid.arrange(detect_input())
        dev.off()
      }
    )
    
    output$downloadImputation <- downloadHandler(
      filename = "imputation.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(imputation_input())
        dev.off()
      }
    )
    
    output$downloadImputation2 <- downloadHandler(
      filename = "imputation (dataset 2).pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(imputation_input2())
        dev.off()
      }
    )
    
    
    #PCA Plot Download
    output$downloadpca <- downloadHandler(
      filename = "PCA.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        plot(pca_input()[[1]])
        dev.off()
      }
    )
    
    output$downloadpca_c2 <- downloadHandler(
      filename = "PCA2.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        plot(pca_c2_input()[[1]])
        dev.off()
      }
    )
    
    #PCA Scree Download
    output$downloadpca_scree <- downloadHandler(
      filename = "PCA_scree.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        plot(pca_input()[[3]])
        dev.off()
      }
    )
    
    output$downloadpca_scree2 <- downloadHandler(
      filename = "PCA_scree2.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        plot(pca_c2_input()[[3]])
        dev.off()
      }
    )
    
    #PCA table 1
    output$downloadpca_txt <- downloadHandler(
      filename = function() { paste("PCA1_values.txt", sep = "") },
      content = function(file) {
        write.table(pca_input()[[2]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    output$downloadpca_c2_txt <- downloadHandler(
      filename = function() { paste("PCA2_values.txt", sep = "") },
      content = function(file) {
        write.table(pca_c2_input()[[2]],
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") 
      }
    )
    
    #Correlation Plot Download
    output$downloadcor <- downloadHandler(
      filename = "Correlation.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        plot(cor_input())
        dev.off()
      }
    )
    
    #Correlation Plot Download (Comp2)
    output$downloadcor_c2 <- downloadHandler(
      filename = "Correlation (Comp2).pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        plot(cor_c2_input())
        dev.off()
      }
    )
    
    #Frequency Plot Download
    output$downloadcond <- downloadHandler(
      filename = "Overlap.pdf",
      content = function(file) {
        cairo_pdf(file, onefile = TRUE, bg = "white", height = input$pdfheight, width = input$pdfwidth)
        plot(cond_input())
        dev.off()
      }
    )
    
    #Frequency Plot Download (Comp2)
    output$downloadcond_c2 <- downloadHandler(
      filename = "Overlap (Comp2).pdf",
      content = function(file) {
        cairo_pdf(file, onefile = TRUE, bg = "white", height = input$pdfheight, width = input$pdfwidth)
        plot(cond_c2_input())
        dev.off()
      }
    )
    
    
    output$downloadNumbers <- downloadHandler(
      filename = "numbers.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(numbers_input())
        dev.off()
      }
    )
    
    output$downloadNumbers2 <- downloadHandler(
      filename = "numbers (dataset 2).pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(numbers_input2())
        dev.off()
      }
    )
    
    output$downloadCoverage <- downloadHandler(
      filename = "coverage.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(coverage_input())
        dev.off()
      }
    )
    
    output$downloadcomplexheatmap <- downloadHandler(
      filename = "comp_heatmap.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight_ht, width = input$pdfwidth_ht)
        print(comp_heatmap_input_pdf()[[1]])
        dev.off()
      }
    )
    
    output$downloadcomplexheatmap_c1 <- downloadHandler(
      filename = "comp_heatmap_c1.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight_ht_c1, width = input$pdfwidth_ht_c1)
        print(comp_heatmap_input_c1_pdf()[[1]])
        dev.off()
      }
    )
    
    output$downloadcomplexheatmap_c2 <- downloadHandler(
      filename = "comp_heatmap_c2.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight_ht_c2, width = input$pdfwidth_ht_c2)
        print(comp_heatmap_input_c2_pdf()[[1]])
        dev.off()
      }
    )
    
    #swatch
    output$downloadswatch <- downloadHandler(
      filename = "Swatch.pdf",
      content = function(file) {
        cairo_pdf(file, height = 6, width = 6)
        print(scales::show_col(colour_functions_input()[[2]]))
        dev.off()
      }
    )
    
    output$downloadcomplexheatmapcluster <- downloadHandler(
      filename = function() {
        paste0("Heatmap_Cluster_", input$cluster_no, ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight_htc, width = input$pdfwidth_htc)
        print(cluster_heatmap_input()[[1]][[input$cluster_no]])
        dev.off()
      }
    )
    
    output$downloadVolcano_enh1 <- downloadHandler(
      filename = function() {
        paste0("Enh_Volcano_", input$volcano_cntrst_enh1, "_", input$p_adj_enh1, ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight_enh, width = input$pdfwidth_enh)
        print(enhanced_volcano_1_input())
        dev.off()
      }
    )
    
    output$downloadVolcano_enh2 <- downloadHandler(
      filename = function() {
        paste0("Enh_Volcano_", input$volcano_cntrst_enh2, "_", input$p_adj_enh2,".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight_enh, width = input$pdfwidth_enh)
        print(enhanced_volcano_2_input())
        dev.off()
      }
    )
    
    output$download_box_plot <- downloadHandler(
      filename = function() {
        paste0("Box_Plot_Cluster_", input$cluster_no_box, "_", input$pres_box, ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight_box, width = input$pdfwidth_box)
        if (input$pres_box=="normalised"){choice=1}else{choice=2}
        print(box_plot_input()[[choice]])
        dev.off()
      }
    )
    
    #####===== Download Report =====#####
    output$downloadReport <-downloadHandler(
      # Create a Progress object
      # For PDF output, change this to "report.pdf"
      filename = "MPP_report.pdf",
      if(input$report_type=="Short"){message_text<-"Generating Report, please wait!"}else{message_text<-"Generating Report, please wait! (If fails retry generation)"},
      content = function(file) {
        withProgress(message = message_text, {
          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          if (input$pval_use == TRUE){padj_method <- "pval"}else{padj_method <- input$padjust_method}
          if(input$report_type=="Short"){
            tempReport <- file.path(tempdir(), "SAINT_DEP_short_report.Rmd")
            file.copy(paste0("SAINT_DEP_short_report.Rmd"), tempReport, overwrite = TRUE)
            line<-paste(Sys.time(), input$comp_type, input$data_mode, padj_method, input$report_type,"",input$imputation, sep = ";")
            write(line,file="report.txt",append=TRUE)
          }
          if(input$report_type=="with Box Plots"){
            tempReport <- file.path(tempdir(), "SAINT_DEP_med_report.Rmd")
            file.copy(paste0("SAINT_DEP_med_report.Rmd"), tempReport, overwrite = TRUE)
            line<-paste(Sys.time(), input$comp_type, input$data_mode, padj_method, input$report_type,"",input$imputation, sep = ";")
            write(line,file="report.txt",append=TRUE)
          }
          if(input$report_type=="Long"){
            if (is.null(GO_BP_report_input())) {
              showNotification("GO plots not ready, please check before dowloading", type = "message", duration = 20)
            } else {
              tempReport <- file.path(tempdir(), "SAINT_DEP_long_report.Rmd")
              file.copy(paste0("SAINT_DEP_long_report.Rmd"), tempReport, overwrite = TRUE)
              line<-paste(Sys.time(), input$comp_type, input$data_mode, padj_method, input$report_type, input$species,input$imputation, sep = ";")
              write(line,file="report.txt",append=TRUE)
            }
          }
          templogo <- file.path(tempdir(), "logo.png")
          file.copy(paste0("logo.png"), templogo, overwrite = TRUE)
          
          templogo_app <- file.path(tempdir(), "app_logo.png")
          file.copy(paste0("app_logo.png"), templogo_app, overwrite = TRUE)
          
          updateRadioButtons(session, "report_type", choices = c("Short", "with Box Plots", "Long"), selected = "Short")
          
          sig_proteins<-dep() %>%
            .[replace_na(SummarizedExperiment::rowData(.)$significant, FALSE), ] %>%
            nrow()
          
          tested_contrasts<- gsub("_p.adj", "", 
                                  colnames(SummarizedExperiment::rowData(dep()))[grep("p.adj", 
                                                                                      colnames(SummarizedExperiment::rowData(dep())))])
          selected_names <- cludf_check_input()[input$tablecludf_rows_selected,1]
          
          #pg_width<- ncol(imp_data()) / 2.5
          clust_height<-(dep_num()/8)+4
          # Set up parameters to pass to Rmd document
          params <- list(data = data,
                         data2 = data2,
                         expdesign = expdesign,
                         expdesign2 = expdesign2,
                         data_mode = input$data_mode,
                         data_type = input$data_type,
                         exp = input$exp,
                         lfq_type = input$lfq_type,
                         data_type2 = input$data_type2,
                         exp2 = input$exp2,
                         lfq_type2 = input$lfq_type2,
                         thr = input$thr,
                         pval_use = input$pval_use,
                         p = input$p,
                         lfc = input$lfc,
                         p2 = input$p2,
                         padjust_method = input$padjust_method,
                         padjust_method2 = input$padjust_method2,
                         lfc2 = input$lfc2,
                         bfdr = input$bfdr,
                         bfdr_freq = input$bfdr_freq,
                         k1 = input$k1,
                         comp_type = input$comp_type,
                         comp= comp,
                         comp2= comp2,
                         imputation = input$imputation,
                         imputation2 = input$imputation2,
                         shift = input$shift,
                         width = input$width,
                         #missval_input = missval_input,
                         # p_hist_input = p_hist_input,
                         pca_input = pca_input,
                         pca_c2_input = pca_c2_input,
                         norm_input = norm_input,
                         norm_input2 = norm_input2,
                         normalisation_choice = input$normalisation_choice,
                         normalisation_choice2 = input$normalisation_choice2,
                         numbers_input = numbers_input,
                         numbers_input2 = numbers_input2, 
                         #coverage_input = coverage_input,
                         cor_input = cor_input,
                         cor_c2_input = cor_c2_input,
                         cond_input = cond_input,
                         cond_c2_input = cond_c2_input,
                         imputation_input = imputation_input,
                         imputation_input2 = imputation_input2,
                         signif_plot_input = signif_plot_input,
                         signif_plot_input2 = signif_plot_input2,
                         heatmap_input = heatmap_input,
                         heatmap_input2 = heatmap_input2,
                         comp_heatmap_input = comp_heatmap_input,
                         comp_heatmap_input_c1 = comp_heatmap_input_c1,
                         comp_heatmap_input_c2 = comp_heatmap_input_c2,
                         venn_input = venn_input,
                         clust_height = clust_height,
                         cluster_heatmap_input = cluster_heatmap_input,
                         colour_functions_input = colour_functions_input,
                         data_results = data_results,
                         p_adj_enh1 = input$p_adj_enh1,
                         dep = dep,
                         dep2 = dep2,
                         enh_vol_data = input$enh_vol_data,
                         enh_vol_connect = input$enh_vol_connect,
                         enh_vol_names = input$enh_vol_names,
                         selected_names = selected_names,
                         dep_num = dep_num,
                         dep_num_c1 = dep_num_c1,
                         dep_num_c2 = dep_num_c2,
                         dep_num_bfdr = dep_num_bfdr,
                         norm_data = norm_data,
                         imp_data = imp_data,
                         palette_choice = input$palette_choice,
                         GO_MF_report_input = GO_MF_report_input,
                         GOMF_cnet_input = GOMF_cnet_input,
                         GO_CC_report_input = GO_CC_report_input,
                         GOCC_cnet_input = GOCC_cnet_input,
                         GO_BP_report_input = GO_BP_report_input,
                         GOBP_cnet_input = GOBP_cnet_input,
                         gomf_level = input$gomf_level,
                         gocc_level = input$gocc_level,
                         gobp_level = input$gobp_level,
                         box_plot_report_input = box_plot_report_input,
                         rendered_by_shiny = TRUE
          )
          
          # Knit the document, passing in the `params` list
          rmarkdown::render(tempReport, output_file = file,
                            params = params,
                            envir = new.env(parent = globalenv())
          )
        })
      }
    )
    
    
    
    
  }) #End of Analyse button event
  
  
})
#End of Server

shinyApp(ui = ui, server = server)


