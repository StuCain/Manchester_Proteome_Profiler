#Shiny_DEP_Heatmaps  by Stuart Cain
#The University of Manchester
#This script imports Maxquant output file proteinGroups.txt and combines it with a experimental
#design file. The DEP package is used to calculate foldchange values, which can be filtered using Cutoff values.
#The complexHeatmap package is used to render heatmaps, and the clusters are then analysed by clusterprofiler

#Use to install packages not listed on CRAN
#install.packages("BiocManager")
#BiocManager::install(c("EnhancedVolcano", "enrichplot",  "clusterProfiler",
#"ComplexHeatmap", "SummarizedExperiment", "DEP", "org.Hs.eg.db", "HDO.db"))

library(DEP)
library(tidyverse) #contains: ggplot2, dplyr, tidyr, 
# readr, purr, tibble, stringr, forcats
library(shiny)
library(shinydashboard)
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

source("R/functions.R")

Version = "0.1.1-beta"


####UI####
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  #autoWaiter(html = spin_loaders(id = 12, color = "#2fa3e6", style = NULL), color = "white"),
  autoWaiter(html = spin_dots(), color = "white"),
  dashboardPage(
    dashboardHeader(title = "DEP + Cluster"),
    dashboardSidebar(width=250,#changed with of sidebar
                     sidebarMenu(
                       menuItem(
                         radioButtons("data_type",
                                      "Data File type",
                                      choices = c("MaxQuant", "Proteome Discoverer"),
                                      #selected = "MaxQuant"
                         )
                       ),
                       menuItem("Files", selected = TRUE,
                                menuItemOutput("file_input"),
                                #                                infoBoxOutput("ibox"), #Info box for testing
                                fileInput('file2',
                                          'ExperimentalDesign.txt',
                                          accept=c('text/csv',
                                                   'text/comma-separated-values,text/plain',
                                                   '.csv')),
                                helpText("Generate Experimental template"),
                                #(=h3("Only one set of comparisons selected", style = 'font-size:16px;color:blue;')
                                downloadButton('downloadexp_template', 'Save Exp. template')
                       ),
                       menuItem("SAINT Files", selected = FALSE,
                                fileInput('saintfile1',
                                          'saint_file1.txt',
                                          accept=c('text/csv',
                                                   'text/comma-separated-values,text/plain',
                                                   '.csv')),
                                fileInput('saintfile2',
                                          'saint_file2.txt',
                                          accept=c('text/csv',
                                                   'text/comma-separated-values,text/plain',
                                                   '.csv')),
                                fileInput('saintfile3',
                                          'saint_file3.txt',
                                          accept=c('text/csv',
                                                   'text/comma-separated-values,text/plain',
                                                   '.csv'))
                       ),
                     
                     menuItemOutput("columns"),
                     menuItem(
                       radioButtons("comp_type",
                                    "Comparisons",
                                    choices = c("1 Round of comparisons", "2 Rounds of comparisons", "Use SaintExpress"),
                                    #selected = "1 Round of comparisons"
                       )
                     ),
                     menuItemOutput("comparisons"),
                     menuItem("Imputation options",
                              radioButtons("imputation",
                                           "Imputation type",
                                           choices = c("man", MsCoreUtils::imputeMethods())[1:9],
                                           selected = "man"),
                              p(a("Detailed information link",
                                  href = "https://rformassspectrometry.github.io/MsCoreUtils/reference/imputation.html",
                                  target="https://rformassspectrometry.github.io/MsCoreUtils/reference/imputation.html"))
                     ),
                     actionButton("analyze", "Analyze"),
                     tags$hr(),
                     uiOutput("downloadTable"),
                     uiOutput("downloadButton"),
                     tags$hr(), #adds a line
                     helpText("Supported by"),
                     imageOutput("logo_img", fill = TRUE)
    )
  ),
  
  #Body
  dashboardBody(
    helpText(paste("By Stuart Cain, Version:", Version)),
    #Top row of boxes
    fluidRow(
      box(numericInput("thr",
                       "Missing value threshold",
                       min = 0, max = 3, value = 1),
          width = 1),
      box(numericInput("bfdr",
                       "BFDR Cut-Off",
                       min = 0, max = 1, value = 0.1),
          width = 1),
      box(numericInput("p",
                       "adj. P value",
                       min = 0.0001, max = 0.1, value = 0.05),
          width = 1),
      box(numericInput("lfc",
                       "Fold change",
                       min = 0, max = 32, value = 2),
          width = 1),
      infoBoxOutput("significantBox"), #shows number of Sig. Proteins
      box(numericInput("p2",
                       "adj. P value (Comp2)",
                       min = 0.0001, max = 0.1, value = 0.05),
          width = 1),
      box(numericInput("lfc2",
                       "Fold change (Comp2)",
                       min = 0, max = 32, value = 2),
          width = 1),
      box(radioButtons("pres",
                       "Data presentation",
                       c("contrast", "centered"),
                       selected = "contrast"),
          width = 1),
      infoBoxOutput("significantBox2"), #shows number of Sig. Proteins
      infoBoxOutput("significantBox3"), #shows number of Sig. Proteins
      infoBoxOutput("significantBox_bfdr"), #shows number of Sig. Proteins (BFDR)
      infoBoxOutput("significantBox_bfdr_filtered") #shows number of Sig. Proteins (BFDR)
    ),
    
    #Main body of
    tabsetPanel(
      tabPanel("DEP", 
               fluidRow(
                 column(width = 6,
                        box(title = "Top Table",
                            DT::dataTableOutput("table"), width = 12)
                 ),
                 column(width = 6,
                        tabBox(id = "result_plots",
                               title = "Result Plots", width = 12,
                               tabPanel(title = "Selected Protein",
                                        plotOutput("selected_plot"),
                                        downloadButton('downloadPlot', 'Save plot')),
                               tabPanel(title = "Heatmap",
                                        fluidRow(
                                          box(numericInput("limit",
                                                           "Color limit (log2)",
                                                           min = 0, max = 16, value = 6),
                                              width = 3),
                                          box(numericInput("size",
                                                           "Heatmap size (4-30)",
                                                           min = 4, max = 30, value = 10),
                                              width = 3),
                                          #Kmeans Cluster box
                                          box(numericInput("k",
                                                           "Kmeans clusters",
                                                           min = 1, max = 6, value = 4),
                                              width = 3)
                                          
                                        ),
                                        fluidRow(
                                          uiOutput("plot"),
                                          downloadButton('downloadHeatmap', 'Save heatmap'))
                               ),
                               tabPanel(title = "Volcano plot",
                                        fluidRow(
                                          box(uiOutput("volcano_cntrst"), width = 6),
                                          box(numericInput("fontsize",
                                                           "Font size",
                                                           min = 0, max = 8, value = 4),
                                              width = 3),
                                          box(checkboxInput("check_names",
                                                            "Display names",
                                                            value = TRUE),
                                              checkboxInput("p_adj",
                                                            "Adjusted p values",
                                                            value = FALSE),
                                              width = 3)
                                        ),
                                        fluidRow(
                                          plotOutput("volcano", height = 600),
                                          downloadButton('downloadVolcano', 'Save volcano')
                                        )
                               ),
                               tabPanel(title = "Volcano plot (Comp2)",
                                        fluidRow(
                                          box(uiOutput("volcano_cntrst_c2"), width = 6),
                                          box(numericInput("fontsize_c2",
                                                           "Font size",
                                                           min = 0, max = 8, value = 4),
                                              width = 3),
                                          box(checkboxInput("check_names_c2",
                                                            "Display names",
                                                            value = TRUE),
                                              checkboxInput("p_adj_c2",
                                                            "Adjusted p values",
                                                            value = FALSE),
                                              width = 3)
                                        ),
                                        fluidRow(
                                          plotOutput("volcano_c2", height = 600),
                                          downloadButton('downloadVolcano_c2', 'Save volcano')
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
                                        plotOutput("pca", height = 600),
                                        downloadButton('downloadpca', 'Save')
                               ),
                               tabPanel(title = "Correlation",
                                        plotOutput("cor", height = 600),
                                        downloadButton('downloadcor', 'Save')
                               ),
                               tabPanel(title = "Freq. Plot",
                                        plotOutput("cond", height = 600),
                                        downloadButton('downloadcond', 'Save')
                               ),
                               #                                 tabPanel(title = "PCA (Comp2)",
                               #                                          plotOutput("pca_c2", height = 600),
                               #                                          downloadButton('downloadpca_c2', 'Save')
                               #                                 ),
                               tabPanel(title = "Correlation (Comp2)",
                                        plotOutput("cor_c2", height = 600),
                                        downloadButton('downloadcor_c2', 'Save')
                               ),
                               tabPanel(title = "Freq. Plot (Comp2)",
                                        plotOutput("cond_c2", height = 600),
                                        downloadButton('downloadcond_c2', 'Save')
                               )
                               
                        ),
                        tabBox(title = "QC Plots", width = 12,
                               tabPanel(title = "Protein Numbers",
                                        plotOutput("numbers", height = 600),
                                        downloadButton('downloadNumbers', 'Save')
                               ),
                               tabPanel(title = "Sample coverage",
                                        plotOutput("coverage", height = 600),
                                        downloadButton('downloadCoverage', 'Save')
                               ),
                               tabPanel(title = "Normalization",
                                        plotOutput("norm", height = 600),
                                        downloadButton('downloadNorm', 'Save')
                               ),
                               tabPanel(title = "Missing values - Quant",
                                        plotOutput("detect", height = 600),
                                        downloadButton('downloadDetect', 'Save')
                               ),
                               tabPanel(title = "Missing values - Heatmap",
                                        plotOutput("missval", height = 600),
                                        downloadButton('downloadMissval', 'Save')
                               ),
                               tabPanel(title = "Imputation",
                                        plotOutput("imputation", height = 600),
                                        downloadButton('downloadImputation', 'Save')
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
                                           min = 1, max = 6, value = 4),
                              checkboxInput("cluster_cols",
                                            "Cluster columns",
                                            value = TRUE),
                              width = 6),
                          box(
                            selectizeInput("high_col",
                                           "High Colour",
                                           choices=c("red","orangered","yellow","green", "darkgreen", "cornflowerblue","blue", "brown", "black"),
                                           selected = "red"),
                            width = 3),
                          box(
                            selectizeInput("low_col",
                                           "Low Colour",
                                           choices=c("red","orangered","yellow","green","darkgreen","cornflowerblue","blue", "brown", "black"),
                                           selected = "blue"),
                            width = 3),
                        ),
                        fluidRow(
                          box(
                            selectizeInput("palette_choice",
                                           "Choose palette",
                                           choices=c("default","Set1","Set2","Set3","Pastel1","Pastel2", "Paired", "Dark2", "Accent", "Reds", "Blues", "Greens", "Greys"),
                                           selected = "default"),
                            width = 6),
                          box(numericInput("alpha_input",
                                           "Alpha",
                                           min = 0, max = 10, value = 10),
                              width = 3),
                          box(numericInput("limit_comp",
                                           "Color limit (log2)",
                                           min = 0, max = 16, value = 6),
                              width = 3)
                        ),
                        
                        #Heatmap in variable box
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
                          downloadButton('downloadcomplexheatmap', 'Save ComplexHeatmap')
                        ),
                        fluidRow(
                          plotOutput("swatch", height = 600),
                          downloadButton('downloadswatch', 'Save Swatch')
                        ),
                        
                 ),
                 
                 #Second Column
                 
                 column(width = 8,
                        tabBox(id = "comp_heat_tab",
                               title = "Cluster Analysis Plots", width = 12,
                               tabPanel(title = "SAINT Analysis",
                                        DT::dataTableOutput("tablesaint"), width = 12
                               ),
                               tabPanel(title = "Results",
                                        DT::dataTableOutput("table2"), width = 12
                               ),
                               tabPanel(title = "Comparison Heatmaps",
                                        fluidRow(
                                          column(width = 6,
                                                 fluidRow(
                                                   #Kmeans Cluster box
                                                   box(numericInput("k1_c1",
                                                                    "Kmeans clusters",
                                                                    min = 1, max = 6, value = 4),
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
                                                   downloadButton('downloadcomplexheatmap_c1', 'Save ComplexHeatmap')
                                                 ),
                                          ),
                                          #Second Column
                                          column(width = 6,
                                                 fluidRow(
                                                   #Kmeans Cluster box
                                                   box(numericInput("k1_c2",
                                                                    "Kmeans clusters",
                                                                    min = 1, max = 6, value = 4),
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
                                                   downloadButton('downloadcomplexheatmap_c2', 'Save ComplexHeatmap')
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
                                                                    min = 1, max = 6, value = 1),
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
                                                   downloadButton('downloadcomplexheatmapcluster', 'Save ComplexHeatmap')
                                                 )
                                          ),
                                          column(width = 4,
                                                 title = "Results",
                                                 DT::dataTableOutput("table3"), #add table for cluster
                                                 downloadButton('downloadclustertable', 'Save Cluster Table')
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
                                                   downloadButton('downloadVolcano_enh1', 'Save volcano')
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
                                                   downloadButton('downloadVolcano_enh2', 'Save volcano')
                                                 )
                                          ),
                                        )
                                        
                               ),
                               tabPanel(title = "Box Plots",
                                        fluidRow(
                                          box(numericInput("cluster_no_box",
                                                           "Cluster No.",
                                                           min = 1, max = 6, value = 1),
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
                                          downloadButton('download_box_plot', 'Save Plot')
                                        ),
                               ),
                               tabPanel(title = "Molecular Function",
                                        fluidRow(
                                          actionButton("analyze_MF", "Analyze Molecular Function"),
                                          helpText("Warning: This can take a couple of minutes"),
                                        ),
                                        fluidRow(
                                          column(width = 5,
                                                 fluidRow(
                                                   plotOutput("GOMF_plot", height = 800),
                                                   downloadButton('download_MF_plot', 'Save Plot')
                                                 )
                                          ),
                                          
                                          #Second Column
                                          
                                          column(width = 7,
                                                 fluidRow(
                                                   plotOutput("GOMF_cnet", height = 800),
                                                   downloadButton('download_MF_cnet', 'Save Plot')
                                                 )
                                          ),
                                        ),
                                        fluidRow(
                                          box(numericInput("pdfwidth_MF",
                                                           "PDF Width",
                                                           min = 3, max = 30, value = 6),
                                              width = 3),
                                          box(numericInput("pdfheight_MF",
                                                           "PDF Height",
                                                           min = 3, max = 30, value = 8),
                                              width = 3)
                                        ),
                               ),
                               tabPanel(title = "Cellular Components",
                                        fluidRow(
                                          actionButton("analyze_CC", "Analyze Cellular Components"),
                                          helpText("Warning: This can take a couple of minutes"),
                                        ),
                                        fluidRow(
                                          column(width = 5,
                                                 fluidRow(
                                                   plotOutput("GOCC_plot", height = 800),
                                                   downloadButton('download_CC_plot', 'Save Plot')
                                                 )
                                          ),
                                          
                                          #Second Column
                                          
                                          column(width = 7,
                                                 fluidRow(
                                                   plotOutput("GOCC_cnet", height = 800),
                                                   downloadButton('download_CC_cnet', 'Save Plot')
                                                 )
                                          ),
                                        ),
                                        fluidRow(
                                          box(numericInput("pdfwidth_CC",
                                                           "PDF Width",
                                                           min = 3, max = 30, value = 6),
                                              width = 3),
                                          box(numericInput("pdfheight_CC",
                                                           "PDF Height",
                                                           min = 3, max = 30, value = 8),
                                              width = 3)
                                        ),
                                        
                               ),
                               tabPanel(title = "Biological Process",
                                        fluidRow(
                                          actionButton("analyze_BP", "Analyze Biological Process"),
                                          helpText("Warning: This can take a couple of minutes"),
                                        ),
                                        fluidRow(
                                          column(width = 5,
                                                 fluidRow(
                                                   plotOutput("GOBP_plot", height = 800),
                                                   downloadButton('download_BP_plot', 'Save Plot')
                                                 )
                                          ),
                                          
                                          #Second Column
                                          
                                          column(width = 7,
                                                 fluidRow(
                                                   plotOutput("GOBP_cnet", height = 800),
                                                   downloadButton('download_BP_cnet', 'Save Plot')
                                                 )
                                          ),
                                        ),
                                        fluidRow(
                                          box(numericInput("pdfwidth_BP",
                                                           "PDF Width",
                                                           min = 3, max = 30, value = 6),
                                              width = 3),
                                          box(numericInput("pdfheight_BP",
                                                           "PDF Height",
                                                           min = 3, max = 30, value = 8),
                                              width = 3)
                                        ),
                               ),
                               tabPanel(title = "App Variables",
                                        fluidRow(
                                          helpText("Save settings to reproduce results in future sessions"),
                                        ),
                                        fluidRow(
                                          DT::dataTableOutput("table4"), #Setting Table
                                          downloadButton('downloadsettings', 'Save Settings')
                                        )
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
                 htmlOutput("inc")
               )
      )
    )
  )
)
)
#End of UI

####Server####

server <- shinyServer(function(input, output) {
  options(shiny.maxRequestSize=200*1024^2)
  
  #### UI functions --------------------------------------------------------
  #### Side bar menu items ####
  #CRUK Logo
  output$logo_img <- renderImage({
    
    list(src = "www/logo.jpg",
         width = "100%",
         height = 100)
    
  }, deleteFile = F)
  
  getPage<-function() {
    return(includeHTML("www/index.html"))
  }
  output$inc<-renderUI({getPage()})
  
  #Side bar menu items
  output$file_input <- renderMenu({
    if (input$data_type == "Proteome Discoverer"){
      fileInput('filePD',
                'PD Output .xlsx',
                accept=c(".xlsx"))}else{
                  fileInput('file1',
                            'ProteinGroups.txt',
                            accept=c('text/csv',
                                     'text/comma-separated-values,text/plain',
                                     '.csv'))
                }
  })
  #Side menu for picking comparisons
  output$columns <- renderMenu({
    menuItem("Columns",
             selectizeInput("name",
                            "Name column",
                            choices=colnames(data()),
                            selected = "Gene.names"),
             selectizeInput("id",
                            "ID column",
                            choices=colnames(data()),
                            selected = "Protein.IDs"),
             if (input$data_type == "MaxQuant"){
               selectizeInput("filt",
                              "Filter on columns" ,
                              colnames(data()),
                              multiple = TRUE,
                              selected = c("Reverse","Potential.contaminant"))
             }
    )
  })
  
  output$comparisons <- renderMenu({
    menuItem("Comparisons",
             if (!is.null(expdesign())) {
               combinations<-combn(unique(expdesign()$condition),2)
               comp1_f<-paste0(combinations[2,],"_vs_",combinations[1,])
               selectizeInput("comp1f",
                              "1: Select forward comparisons",
                              choices=make.names(comp1_f),
                              multiple = TRUE)
             },
             
             if (!is.null(expdesign())) {
               combinations<-combn(unique(expdesign()$condition),2)
               comp1_r<-paste0(combinations[1,],"_vs_",combinations[2,])
               selectizeInput("comp1r",
                              "1: Select reverse comparisons",
                              choices=make.names(comp1_r),
                              multiple = TRUE)
             },
             
             if (!is.null(expdesign()) && input$comp_type == "2 Rounds of comparisons") {
               combinations<-combn(unique(expdesign()$condition),2)
               comp2_f<-paste0(combinations[2,],"_vs_",combinations[1,])
               selectizeInput("comp2f",
                              "2: Select forward comparisons",
                              choices=make.names(comp2_f),
                              multiple = TRUE)
             },
             
             if (!is.null(expdesign()) && input$comp_type == "2 Rounds of comparisons") {
               combinations<-combn(unique(expdesign()$condition),2)
               comp2_r<-paste0(combinations[1,],"_vs_",combinations[2,])
               selectizeInput("comp2r",
                              "2: Select forward comparisons",
                              choices=make.names(comp2_r),
                              multiple = TRUE)
             }
    )
  })
  #Hide/show comparison heatmaps
  observeEvent(input$comp_type, {
    if (input$comp_type == "1 Round of comparisons") {
      hideTab(inputId =  "comp_heat_tab", target = "Comparison Heatmaps")
    }else{showTab(inputId =  "comp_heat_tab", target = "Comparison Heatmaps")}
    if (input$comp_type == "1 Round of comparisons"|input$comp_type == "2 Rounds of comparisons") {
      hideTab(inputId =  "comp_heat_tab", target = "SAINT Analysis")
    }else{showTab(inputId =  "comp_heat_tab", target = "SAINT Analysis")}
    
  })
  
  #Hide/show Volcano Plts
  observeEvent(input$comp_type, {
    if (input$comp_type == "1 Round of comparisons"|input$comp_type == "Use SaintExpress") {
      hideTab(inputId =  "result_plots", target = "Volcano plot (Comp2)")
    }else{showTab(inputId =  "result_plots", target = "Volcano plot (Comp2)")}
  })
  
  #Hide/show Correlation plots
  #Hide/show Volcano Plts
  observeEvent(input$comp_type, {
    #    if (input$comp_type == "1 Round of comparisons") {
    #      hideTab(inputId =  "correlation_plots", target = "PCA (Comp2)")
    #    }else{showTab(inputId =  "correlation_plots", target = "PCA (Comp2)")}
    if (input$comp_type == "1 Round of comparisons") {
      hideTab(inputId =  "correlation_plots", target = "Correlation (Comp2)")
    }else{showTab(inputId =  "correlation_plots", target = "Correlation (Comp2)")}
    if (input$comp_type == "1 Round of comparisons") {
      hideTab(inputId =  "correlation_plots", target = "Freq. Plot (Comp2)")
    }else{showTab(inputId =  "correlation_plots", target = "Freq. Plot (Comp2)")}
  })
  
  #### Reactive functions  --------------------------------------------------
  set.seed(123) #random seed to give consistent kmeans clustering
  
  #Output information for testing
  #  output$ibox <- renderInfoBox({
  #    infoBox(
  #      "Title",
  #      input$data_type,
  #      #icon = icon("credit-card")
  #    )
  #  })
  
  
  data <- reactive({
    if (input$data_type == "MaxQuant"){
      inFile <- input$file1
      if (is.null(inFile))
        return(NULL)
      maxdata<-read.csv(inFile$datapath, header = TRUE,
                        sep = "\t", stringsAsFactors = FALSE) %>%
        mutate(id = row_number())
      data1<-maxdata[, c(grep("Protein.IDs", colnames(maxdata)),grep("Protein.names", colnames(maxdata)),grep("Gene.names", colnames(maxdata)),  grep("LFQ.",colnames(maxdata)), grep("Reverse",colnames(maxdata)), grep("Potential.contaminant",colnames(maxdata)))]
    }
    
    if (input$data_type == "Proteome Discoverer"){
      inFile <- input$filePD
      if (is.null(inFile))
        return(NULL)
      pd<-read_excel(path = inFile$datapath, sheet = 1, .name_repair = "universal") %>%
        mutate(id = row_number())
      data<-pd %>% 
        dplyr::rename(
          Gene.names=Gene.Symbol,
          Protein.IDs=Accession,
          Protein.names=Description  
        ) 
      #      data1<-data[, c(1,2,3,  grep("Abundance..*Sample..", colnames(data)))]
      data1<-data[, c(grep("Gene.names", colnames(data)), grep("Protein.IDs", colnames(data)), grep("Protein.names", colnames(data)),  grep("Abundance..*Sample..", colnames(data)), grep("Abundance..*Control..", colnames(data)))]
      
    }
    return(data1)
  })
  
  #SAINT File imports
  saintfile1 <- reactive({
    inFile <- input$saintfile1
    if (is.null(inFile))
      return(NULL)
    data<-read.csv(inFile$datapath, header = TRUE,
             sep = "\t", stringsAsFactors = FALSE)
    data<-data.frame(data)
    data<-data %>% dplyr::filter(BFDR <= input$bfdr)
    saint_BFDR <- data[, c(grep("Bait", colnames(data)),grep("PreyGene", colnames(data)),grep("BFDR", colnames(data)))]
    return(saint_BFDR)
  })
  
  saintfile2 <- reactive({
    inFile <- input$saintfile2
    if (is.null(inFile))
      return(NULL)
    data<-read.csv(inFile$datapath, header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
    data<-data.frame(data)
    data<-data %>% dplyr::filter(BFDR <= input$bfdr)
    saint_BFDR <- data[, c(grep("Bait", colnames(data)),grep("PreyGene", colnames(data)),grep("BFDR", colnames(data)))]
    return(saint_BFDR)
  })
  
  saintfile3 <- reactive({
    inFile <- input$saintfile3
    if (is.null(inFile))
      return(NULL)
    data<-read.csv(inFile$datapath, header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
    data<-data.frame(data)
    data<-data %>% dplyr::filter(BFDR <= input$bfdr)
    saint_BFDR <- data[, c(grep("Bait", colnames(data)),grep("PreyGene", colnames(data)),grep("BFDR", colnames(data)))]
    return(saint_BFDR)
  })
  
  saintlist <- reactive({
    unfilt_saint<- c(saintfile1()) #Make single list
    if (!is.null(saintfile2())){unfilt_saint<- rbind(saintfile1(),saintfile2())}
    if (!is.null(saintfile3())){unfilt_saint<- rbind(saintfile1(),saintfile2(),saintfile3())}
    unfilt_saint<-data.frame(unfilt_saint)
#    return(unfilt_saint)
    
    saint_count<-as.data.frame(table(unfilt_saint$PreyGene)) #Make into dataframe
    saint_count<-saint_count %>%dplyr::rename(PreyGene=Var1)%>% dplyr::rename(BFDR_Freq=Freq) %>% filter(PreyGene !=(""))
    saint_count_no <-nrow(saint_count)
    saintlist_out<-list(unfilt_saint,saint_count,saint_count_no)
    return(saintlist_out)
    
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
      exp_template_new<- data.frame(label, condition, replicate)
    }
    
    if (input$data_type == "Proteome Discoverer"){
      Sample_columns<- grep("Abundance..*Sample..", colnames(data()))
      Control_columns<- grep("Abundance..*Control..", colnames(data()))
      LFQ_columns<-c(Sample_columns, Control_columns)
      label<-colnames(data()[LFQ_columns])
      condition<-""
      replicate<-""
      exp_template_new<- data.frame(label, condition, replicate)
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
  
  
  #data annotations for clusters
  data_defn <- reactive({
    data_unique <- make_unique(data(), "Gene.names", "Protein.IDs", delim = ";")
    data_unique_defn <-data_unique[, c("name", "Protein.names", "ID")]
    data_unique_defn$Protein.names  <- sub('\\;.*', '', data_unique_defn$Protein.names)
    data_unique_defn<-as.data.frame(data_unique_defn)
  })
  
  expdesign <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE,
             sep = "\t", stringsAsFactors = FALSE) %>%
      mutate(id = row_number())
  })
  
  filt <- reactive({
    data <- data()
    
    if (input$data_type == "MaxQuant"){
      cols <- grep("^LFQ", colnames(data))
      filtered <- DEP:::filter_MaxQuant(data, input$filt)
      unique_names <- make_unique(filtered, input$name, input$id)
    }
    
    if (input$data_type == "Proteome Discoverer"){
      Sample_columns<- grep("Abundance..*Sample..", colnames(data()))
      Control_columns<- grep("Abundance..*Control..", colnames(data()))
      cols<-c(Sample_columns, Control_columns)
      unique_names <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
    }
    
    se <- make_se(unique_names, cols, expdesign())
    filter_missval(se, input$thr) #filtering set less stringent
  })
  
  norm <- reactive({
    normalize_vsn(filt())
  })
  
  #Get Normalized data
  norm_data <- reactive({
    data_norm<-SummarizedExperiment::assay(norm())
    data_norm<-as.data.frame(data_norm)
    return(data_norm)
  })
  
  imp <- reactive({
    DEP::impute(norm(), input$imputation)
  })
  
  #Get imputed data table
  imp_data <- reactive({
    data_imp<-SummarizedExperiment::assay(imp())
    return(as.data.frame(data_imp))
  })
  
  comp <- reactive({
    comp<-c(input$comp1f,input$comp1r)
  })
  
  comp2 <- reactive({
    if (input$comp_type == "1 Round of comparisons")
      return(NULL)
    comp2<-c(input$comp2f,input$comp2r)
  })
  
  
  df <- reactive({
    test_diff(imp(), type = "manual", test = comp())
  })
  
  df2 <- reactive({
    if (input$comp_type == "1 Round of comparisons")
      return(NULL)
    if (input$comp_type == "2 Rounds of comparisons"){
      test_diff(imp(), type = "manual", test = comp2())
    }
  })
  
  dep <- reactive({
    add_rejections(df(), input$p, log2(input$lfc))
  }) 
  
  dep2 <- reactive({
    if (input$comp_type == "1 Round of comparisons")
      return(NULL)
    add_rejections(df2(), input$p2, log2(input$lfc2))
  }) 
  
  
  dep_num2 <- reactive({
    if (input$comp_type == "1 Round of comparisons")
      return(NULL)
    num_signif <- dep2() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  })
  
  ### All object and functions upon 'Analyze' input ----------------------
  
  observeEvent(input$analyze, {
    
    #### Interactive UI functions  ------------------------------------------
    output$downloadTable <- renderUI({
      selectizeInput("dataset",
                     "Choose a dataset to save" ,
                     c("results","significant_proteins",
                       "displayed_subset","full_dataset"))
    })
    
    output$downloadButton <- renderUI({
      downloadButton('downloadData', 'Save')
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
      info_box<-box(width = 1, height = 0)
      
      if (input$comp_type == "2 Rounds of comparisons")
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
      info_box<-box(width = 1, height = 0)
      
      if (input$comp_type == "2 Rounds of comparisons")
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
      info_box<-box(width = 1, height = 0)
      
      if (input$comp_type == "Use SaintExpress")
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
      info_box<-box(width = 1, height = 0)
      
      if (input$comp_type == "Use SaintExpress")
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
    
    
    output$ibox <- renderInfoBox({
      infoBox(
        "Title",
        input$file_type,
        #icon = icon("credit-card")
      )
    })
    
    
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
    
    
    #res <- reactive({
    #  get_results(dep())
    #})
    
    #    table <- reactive({
    #      DEP:::get_table(res(), input$pres)
    #    })
    
    #Make results table 
    data_dep <- reactive({
      get_results(dep()) %>% filter(significant) %>% select(-significant) 
    })
    
    data_dep2 <- reactive({
      if (input$comp_type == "1 Round of comparisons")
        return(NULL)
      get_results(dep2()) %>% filter(significant) %>% select(-significant) 
    })
    
    data_results1<- reactive({
      data.frame(data_dep())
    })
    
    data_results2<- reactive({
      if (input$comp_type == "1 Round of comparisons")
        return(NULL)
      data.frame(data_dep2())
    })
    
    data_results3<- reactive({
      if (input$comp_type == "1 Round of comparisons")
        return(NULL)
      data_results2() [data_results2()$name %in% data_results1()$name, ]
    })
    
    #Unfiltered BFDR results
    data_results4<- reactive({
      if (input$comp_type == "1 Round of comparisons" | input$comp_type == "2 Rounds of comparisons")
        return(NULL)
      data<-get_results(dep())
      saint_data<-data.frame(saintlist()[[2]])
      data<-data[data$name %in% saint_data$PreyGene, ]
    })
    
    data_results5<- reactive({
      if (input$comp_type == "1 Round of comparisons" | input$comp_type == "2 Rounds of comparisons")
        return(NULL)
      data_results4() [data_results4()$name %in% data_results1()$name, ]
    })
    
    
    
    data_results<- reactive({
      if (input$comp_type == "2 Rounds of comparisons"){
        data_results3()}
      if (input$comp_type == "Use SaintExpress"){
        data_results5()
      }
      else{data_results1()}
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
    
    selected_plot_input <- reactive ({
      if(!is.null(input$table_rows_selected)) {
        selected_id <- data_results()[input$table_rows_selected,1]
        plot_single(dep(), selected_id, input$pres)
      }
    })
    
    heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        plot_heatmap(dep(),
                     input$pres,
                     kmeans = TRUE,
                     input$k,
                     input$limit)
      })
    })
    
    volcano_input <- reactive({
      if(!is.null(input$volcano_cntrst)) {
        plot_volcano(dep(),
                     input$volcano_cntrst,
                     input$fontsize,
                     input$check_names,
                     input$p_adj)
      }
    })
    
    volcano_c2_input <- reactive({
      if(!is.null(input$volcano_cntrst_c2)) {
        plot_volcano(dep2(),
                     contrast = input$volcano_cntrst_c2,
                     label_size = input$fontsize_c2,
                     add_names = input$check_names_c2,
                     input$p_adj_c2)
      }
    })
    
    
    norm_input <- reactive({
      plot_normalization(filt(),
                         norm())
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
    
    numbers_input <- reactive({
      plot_numbers(norm())
    })
    
    coverage_input <- reactive({
      plot_coverage(norm())
    })
    
    #PCA Input
    pca_input <- reactive({
      plot_pca(dep(), x = 1, y = 2, point_size = 4)
    })
    
    #    #PCA Input (Comp2)
    #    pca_c2_input <- reactive({
    #      plot_pca(dep2(), x = 1, y = 2, point_size = 4)
    #    })
    
    #Correlation Input
    cor_input <- reactive({
      plot_cor(dep(), significant = TRUE, lower = 0, upper = 1, pal = "Reds", indicate = c("condition", "replicate"))
    })
    
    #Correlation Input (Comp2)
    cor_c2_input <- reactive({
      plot_cor(dep2(), significant = TRUE, lower = 0, upper = 1, pal = "Reds", indicate = c("condition", "replicate"))
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
        if (input$comp_type == "1 Round of comparisons"){
          ht_title<-paste("FC=", input$lfc, " pVal=", input$p, "<br>", "(",dep_num_c1(), " Proteins)", "<br>",input$k1," Clusters ",input$pres , sep="")
          
          plot_complex_heatmap(data_results(),
                               #                     input$pres,
                               #                     kmeans = TRUE,
                               input$k1,
                               ht_title,
                               colour_functions_input()[[1]],
                               colour_functions_input()[[2]],
                               font_size=8,
                               cluster_cols=input$cluster_cols
          )
        }

        if (input$comp_type == "Use SaintExpress"){
          ht_title<-paste("SAINTexpress Filtered", "<br>", "Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>", "BDFR_Cutoff=", input$bdfr,  "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ",input$pres, sep="")
          
          plot_complex_heatmap_BDFR(data_results(),
                                      #                     input$pres,
                                      #                     kmeans = TRUE,
                                      input$k1,
                                      ht_title,
                                      colour_functions_input()[[1]],
                                      colour_functions_input()[[2]],
                                      font_size=8,
                                      cluster_cols=input$cluster_cols,
                                      data_results1=data_results1(),
                                      saintlist=saintlist()[[2]])
        }
        
        if (input$comp_type == "2 Rounds of comparisons"){
          ht_title<-paste("Comparison 1+2", "<br>", "Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>", "Comp 2 FC=", input$lfc2, " pVal=", input$p2, "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ",input$pres, sep="")
          
          plot_complex_heatmap_2xcomp(data_results(),
                                      #                     input$pres,
                                      #                     kmeans = TRUE,
                                      input$k1,
                                      ht_title,
                                      colour_functions_input()[[1]],
                                      colour_functions_input()[[2]],
                                      font_size=8,
                                      cluster_cols=input$cluster_cols,
                                      data_results1=data_results1(),
                                      data_results2=data_results2())
        }
        
      })
    })
    
    comp_heatmap_input_pdf <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$comp_type == "1 Round of comparisons"){
          ht_title<-paste("FC=", input$lfc, " pVal=", input$p, "<br>", "(",dep_num_c1(), " Proteins)", "<br>",input$k1," Clusters ",input$pres , sep="")
          
          plot_complex_heatmap(data_results(),
                               #                     input$pres,
                               #                     kmeans = TRUE,
                               input$k1,
                               ht_title,
                               colour_functions_input()[[1]],
                               colour_functions_input()[[2]],
                               font_size=6,
                               cluster_cols=input$cluster_cols
          )
        }
        
        if (input$comp_type == "Use SaintExpress"){
          ht_title<-paste("SAINTexpress Filtered", "<br>", "Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>", "BDFR_Cutoff=", input$bdfr,  "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ",input$pres, sep="")
          
          plot_complex_heatmap_BDFR(data_results(),
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1,
                                    ht_title,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=6,
                                    cluster_cols=input$cluster_cols,
                                    data_results1=data_results1(),
                                    saintlist=saintlist()[[2]])
        }
        
        if (input$comp_type == "2 Rounds of comparisons"){
          ht_title<-paste("Comparison 1+2", "<br>", "Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>", "Comp 2 FC=", input$lfc2, " pVal=", input$p2, "<br>", "(", dep_num(), " Proteins)", "<br>",input$k1," Clusters ",input$pres, sep="")
          
          plot_complex_heatmap_2xcomp(data_results(),
                                      #                     input$pres,
                                      #                     kmeans = TRUE,
                                      input$k1,
                                      ht_title,
                                      colour_functions_input()[[1]],
                                      colour_functions_input()[[2]],
                                      font_size=6,
                                      cluster_cols=input$cluster_cols,
                                      data_results1=data_results1(),
                                      data_results2=data_results2())
        }
        
      })
    })
    
    comp_heatmap_input_c1 <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$comp_type == "2 Rounds of comparisons"){
        ht_title<-paste("Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>",  "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c1," Clusters ",input$pres, sep="")
        plot_complex_heatmap_2xcomp(data_results1(),
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1_c1,
                                    ht_title,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=8,
                                    cluster_cols=input$cluster_cols,
                                    data_results1=data_results1(),
                                    data_results2=data_results2())
        }
        
        if (input$comp_type == "Use SaintExpress"){
          ht_title<-paste("Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>", "BFDR_Cutoff=", input$bfdr,  "<br>", "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c1," Clusters ",input$pres, sep="")
          
          plot_complex_heatmap_BDFR(data_results1(),
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1_c1,
                                    ht_title,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=8,
                                    cluster_cols=input$cluster_cols,
                                    data_results1=data_results1(),
                                    saintlist=saintlist()[[2]])
        }
        
      })
    })
    
    comp_heatmap_input_c1_pdf <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        ht_title<-paste("Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>",  "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c1," Clusters ",input$pres, sep="")
        
        plot_complex_heatmap_2xcomp(data_results1(),
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1_c1,
                                    ht_title,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=6,
                                    cluster_cols=input$cluster_cols,
                                    data_results1=data_results1(),
                                    data_results2=data_results2())
      })
    })
    
    comp_heatmap_input_c2 <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$comp_type == "2 Rounds of comparisons"){
          ht_title<-paste("Comp 2 FC=", input$lfc2, " pVal=", input$p2, "<br>",  "(", dep_num_c2(), " Proteins)", "<br>",input$k1_c2," Clusters ",input$pres, sep="")
          
          plot_complex_heatmap_2xcomp(data_results2(),
                                      #                     input$pres,
                                      #                     kmeans = TRUE,
                                      input$k1_c2,
                                      ht_title,
                                      colour_functions_input()[[1]],
                                      colour_functions_input()[[2]],
                                      font_size=8,
                                      cluster_cols=input$cluster_cols,
                                      data_results1=data_results1(),
                                      data_results2=data_results2())
        }
        
        if (input$comp_type == "Use SaintExpress"){
          ht_title<-paste("SAINTexpress List", "<br>", "Comp 1 FC=", input$lfc, " pVal=", input$p, "<br>", "BFDR_Cutoff=", input$bfdr,  "<br>", "(", dep_num_c1(), " Proteins)", "<br>",input$k1_c2," Clusters ",input$pres, sep="")
          data<-data.frame(data_results4())
          plot_complex_heatmap_BDFR(data,
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1_c2,
                                    ht_title,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=8,
                                    cluster_cols=input$cluster_cols,
                                    data_results1=data_results1(),
                                    saintlist=saintlist()[[2]])
        }
      })
    })
    
    comp_heatmap_input_c2_pdf <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        ht_title<-paste("Comp 2 FC=", input$lfc2, " pVal=", input$p2, "<br>",  "(", dep_num_c2(), " Proteins)", "<br>",input$k1_c2," Clusters ",input$pres, sep="")
        
        plot_complex_heatmap_2xcomp(data_results2(),
                                    #                     input$pres,
                                    #                     kmeans = TRUE,
                                    input$k1_c2,
                                    ht_title,
                                    colour_functions_input()[[1]],
                                    colour_functions_input()[[2]],
                                    font_size=6,
                                    cluster_cols=input$cluster_cols,
                                    data_results1=data_results1(),
                                    data_results2=data_results2())
      })
    })
    
    #Return my_mat from heatmap (object [[2]])
    data_results_mat<- reactive({
      data.frame(comp_heatmap_input()[[2]])
    })
    
    #    ht<-comp_heatmap_input()[[1]]
    #    my_mat<-comp_heatmap_input()[[2]]
    #    dfdata_sig2<-comp_heatmap_input()[[3]]
    
    #Heatmap Clusters Input
    cluster_heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        if (input$comp_type == "2 Rounds of comparisons"){
          plot_heatmap_clusters_2xcomp(comp_heatmap_input()[[1]],
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
        
        if (input$comp_type == "Use SaintExpress"){
          plot_heatmap_clusters_2xcomp(comp_heatmap_input()[[1]],
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
        
        if (input$comp_type == "1 Round of comparisons"){
          plot_heatmap_clusters(comp_heatmap_input()[[1]],
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
      })
    }) 
    
    #Enhanced Volcano Input 1
    enhanced_volcano_1_input <- reactive({
      if (input$comp_type == "2 Rounds of comparisons"){
        lfc<-input$lfc2
        p<-input$p2
      }else{lfc<-input$lfc
      p<-input$p
      }
      input_volcano_cntrst_enh1 <- req(input$volcano_cntrst_enh1)
      plot_enhanced_volcano (data_results(),
                             cluster_heatmap_input()[[19]],
                             comp_pick = input_volcano_cntrst_enh1,
                             FC=lfc,
                             pVal=p,
                             padj_check=input$p_adj_enh1,
                             colour_functions_input()[[2]],
                             colour_functions_input()[[3]]
      )
    })
    
    #Enhanced Volcano Input 2
    enhanced_volcano_2_input <- reactive({
      if (input$comp_type == "2 Rounds of comparisons"){
        lfc<-input$lfc2
        p<-input$p2
      }else{lfc<-input$lfc
      p<-input$p
      }
      input_volcano_cntrst_enh2 <- req(input$volcano_cntrst_enh2)
      plot_enhanced_volcano (data_results(),
                             cluster_heatmap_input()[[19]],
                             comp_pick = input_volcano_cntrst_enh2,
                             FC=lfc,
                             pVal=p,
                             padj_check=input$p_adj_enh2,
                             colour_functions_input()[[2]],
                             colour_functions_input()[[3]]
                             
      )
    })
    
    #Plot_Box_Plots
    box_plot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        plot_box_plot (norm_data(),
                       imp_data(),
                       clu_df=cluster_heatmap_input()[[19]],
                       cluster_no=input$cluster_no_box,
                       input$palette_choice
        )
      })
    })
    
    #Make entrezlist
    entrezlist_input <- reactive({
      get_entrezlist(cluster_heatmap_input()[[19]])
    })
    
    #Make Variables Table
    settings_input <- reactive({
      file1<-input$file1[1]
      file2<-input$file2[1]
      filePD<-input$filePD[1]
      comp_1<-paste(comp(),collapse=',')
      comp_2<-paste(comp2(),collapse=',')
      if(is.null(file1) == TRUE){file1<-"NA"}
      if(is.null(file2) == TRUE){file2<-"NA"}
      if(is.null(filePD) == TRUE){filePD<-"NA"}
      values_all<-rbind (input$data_type, file1, filePD, file2,  input$comp_type, comp_1, comp_2)
      variables_all <- c("Data_type", "MaxQuant_File", "PD_File", "Experimental_Settings_File", "Comparison_No.", "Comp1", "Comp2")
      all_settings<-data.frame(
        Variables=variables_all,
        Set_Values=values_all)
      #all_settings
      
      myvalues <- NULL
      newvalues <- NULL
      for(i in c(1:29,31:56)){
        newvalues <- paste(names(input)[i], input[[names(input)[i]]], sep=":")
        myvalues <- append(myvalues, newvalues)
      }
      DF <- data.frame(do.call(rbind, strsplit(myvalues, ":", fixed=TRUE)))
      DF <- tibble::rownames_to_column(DF, "variable")
      DF<-data.frame(DF[-1])
      colnames(DF) <- colnames(all_settings)
      total_settings<-rbind(all_settings,DF)
      total_settings
      
    })
    
    
    #### Output functions  --------------------------------------------------
    output$table <- DT::renderDataTable({
      data_results4()
      #data_results()
      #     data() to test data input
    }, options = list(pageLength = 25, scrollX = T),
    selection = list(selected = c(1)))
    
    #Output of table on 2nd tab
    output$table2 <- DT::renderDataTable({
      #data_results_mat()
      comp_heatmap_input()[[3]]
    }, options = list(pageLength = 25, scrollX = T),
    )
    
    #Output of table on Cluster tab
    output$table3 <- DT::renderDataTable({
      cluster_heatmap_input()[[input$cluster_no+12]]
    }, options = list(pageLength = 25, scrollX = T),
    )
    
    #settings table
    output$table4 <- DT::renderDataTable({
      settings_input()
    }) 
    #Saint table
    output$tablesaint <- DT::renderDataTable({
      saintlist()[[2]]
    }) 

    output$selected_plot <- renderPlot({
      selected_plot_input()
    })
    
    output$heatmap <- renderPlot({
      heatmap_input()
    })
    
    output$volcano <- renderPlot({
      volcano_input()
    })
    
    output$volcano_c2 <- renderPlot({
      volcano_c2_input()
    })
    
    
    output$norm <- renderPlot({
      norm_input()
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
    
    output$numbers <- renderPlot({
      numbers_input()
    })
    
    output$coverage <- renderPlot({
      coverage_input()
    })
    
    #PCA output
    output$pca <- renderPlot({
      pca_input()
    })
    
    #    #PCA output (Comp2)
    #    output$pca_c2 <- renderPlot({
    #      pca_c2_input()
    #    })
    
    
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
    
    #cluster heatmap output
    output$cluster_heatmap <- renderPlot({
      cluster_heatmap_input()[[input$cluster_no]]
    })
    
    #Box Plots
    output$box_plot <- renderPlot({
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
    
    
    ####Observe####
    observe({
      output$plot <- renderUI({
        plotOutput("heatmap", height = (100 * as.numeric(input$size)))
      })
    })
    
    observe({
      output$plot_comp <- renderUI({
        plotOutput("comp_heatmap", height = ((dep_num()*12)+270))
      })
    })
    
    observe({
      output$plot_comp_c1 <- renderUI({
        plotOutput("comp_heatmap_c1", height = ((dep_num_c1()*12)+270))
      })
    })
    
    observe({
      output$plot_comp_c2 <- renderUI({
        plotOutput("comp_heatmap_c2", height = ((dep_num_c2()*12)+270))
      })
    })
    
    observe({
      output$plot_clust <- renderUI({
        plotOutput("cluster_heatmap", height = (cluster_heatmap_input()[[input$cluster_no+6]]*12)+270)
      })
    })
    
    
    #Observe GO Plots
    #Observe MF
    observeEvent(input$analyze_MF,{
      GO_MF_input <- reactive({
        withProgress(message = 'Plotting', detail = 'This may take a while...', max=10, value = 6.6,{
          for (i in 1:100) {
            incProgress(1/100)
            Sys.sleep(0.5)
          }
          plot_GO_MF(entrezlist_input())
        })
      })
      #Plot output
      output$GOMF_plot <- renderPlot({
        GO_MF_input()[[1]]
      })
      output$GOMF_cnet <- renderPlot({
        GO_MF_input()[[2]]
      })
      #PDF Output
      output$download_MF_plot <- downloadHandler(
        filename = "Molecular_Function_Plot.pdf",
        content = function(file) {
          cairo_pdf(file, height = input$pdfheight_MF, width = input$pdfwidth_MF)
          print(GO_MF_input()[[1]])
          dev.off()
        }
      )
      output$download_MF_cnet <- downloadHandler(
        filename = "Molecular_Function_Map.pdf",
        content = function(file) {
          cairo_pdf(file, height = input$pdfheight_MF, width = input$pdfwidth_MF)
          print(GO_MF_input()[[2]])
          dev.off()
        }
      )
    })
    
    #Observe CC
    observeEvent(input$analyze_CC,{
      GO_CC_input <- reactive({
        withProgress(message = 'Plotting', detail = 'This may take a while...', max=10, value = 6.6,{
          for (i in 1:100) {
            incProgress(1/100)
            Sys.sleep(0.5)
          }
          plot_GO_CC(entrezlist_input())
        })
      })
      #Plot output
      output$GOCC_plot <- renderPlot({
        GO_CC_input()[[1]]
      })
      output$GOCC_cnet <- renderPlot({
        GO_CC_input()[[2]]
      })
      #PDF Output
      output$download_CC_plot <- downloadHandler(
        filename = "Cellular_Components_Plot.pdf",
        content = function(file) {
          cairo_pdf(file, height = input$pdfheight_CC, width = input$pdfwidth_CC)
          print(GO_CC_input()[[1]])
          dev.off()
        }
      )
      output$download_CC_cnet <- downloadHandler(
        filename = "Cellular_Components_Map.pdf",
        content = function(file) {
          cairo_pdf(file, height = input$pdfheight_CC, width = input$pdfwidth_CC)
          print(GO_CC_input()[[2]])
          dev.off()
        }
      )
    })
    
    #Observe BP
    observeEvent(input$analyze_BP,{
      GO_BP_input <- reactive({
        withProgress(message = 'Plotting', detail = 'This may take a while...', max=10, value = 6.6,{
          for (i in 1:100) {
            incProgress(1/100)
            Sys.sleep(0.5)
          }
          plot_GO_BP(entrezlist_input())
        })
      })
      #Plot output
      output$GOBP_plot <- renderPlot({
        GO_BP_input()[[1]]
      })
      output$GOBP_cnet <- renderPlot({
        GO_BP_input()[[2]]
      })
      #PDF Output
      output$download_BP_plot <- downloadHandler(
        filename = "Biological_Process_Plot.pdf",
        content = function(file) {
          cairo_pdf(file, height = input$pdfheight_BP, width = input$pdfwidth_BP)
          print(GO_BP_input()[[1]])
          dev.off()
        }
      )
      output$download_BP_cnet <- downloadHandler(
        filename = "Biological_Process_Map.pdf",
        content = function(file) {
          cairo_pdf(file, height = input$pdfheight_BP, width = input$pdfwidth_BP)
          print(GO_BP_input()[[2]])
          dev.off()
        }
      )
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
             "full_dataset" = get_df_wide(dep()))
    })
    
    output$downloadData <- downloadHandler(
      filename = function() { paste(input$dataset, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )
    
    #Cluster table
    output$downloadclustertable <- downloadHandler(
      filename = function() { paste("Cluster ", input$cluster_no, ".txt", sep = "") },
      content = function(file) {
        write.table(cluster_heatmap_input()[[input$cluster_no+12]],
                    file,
                    col.names = TRUE,
                    row.names = TRUE,
                    sep ="\t") }
    )
    
    #Variable table
    output$downloadsettings <- downloadHandler(
      filename = "Settings.txt",
      content = function(file) {
        write.table(settings_input(),
                    file,
                    col.names = TRUE,
                    row.names = TRUE,
                    sep ="\t") }
    )
    
    #   #Variable table V2
    #    output$downloadsettingsV2 <- downloadHandler(
    #      filename = "Settings.txt",
    #     content = function(file) {
    ##       write.table(AllInputs(),
    #                  file,
    ##                  col.names = TRUE,
    #                 row.names = TRUE,
    #                  sep ="\t") }
    #  )
    
    
    
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0("Barplot_", data_results()[input$table_rows_selected,1], ".pdf")
      },
      content = function(file) {
        cairo_pdf(file)
        print(selected_plot_input())
        dev.off()
      }
    )
    
    output$downloadHeatmap <- downloadHandler(
      filename = 'Heatmap.pdf',
      content = function(file) {
        cairo_pdf(file, height = 10, paper = "a4")
        print(heatmap_input())
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
    
    output$downloadNorm <- downloadHandler(
      filename = "normalization.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        print(norm_input())
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
    
    #PCA Plot Download
    output$downloadpca <- downloadHandler(
      filename = "PCA.pdf",
      content = function(file) {
        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
        plot(pca_input())
        dev.off()
      }
    )
    
    #    #PCA Plot Download (comp2)
    #    output$downloadpca <- downloadHandler(
    #      filename = "PCA (Comp2).pdf",
    #      content = function(file) {
    #        cairo_pdf(file, height = input$pdfheight, width = input$pdfwidth)
    #        print(pca_c2_input())
    #        dev.off()
    #      }
    #    )
    
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
        print(cluster_heatmap_input()[[input$cluster_no]])
        dev.off()
      }
    )
    
    output$downloadVolcano_enh1 <- downloadHandler(
      filename = function() {
        paste0("Enh_Volcano_", input$volcano_cntrst_enh1, "_", input$p_adj_enh1, ".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = 8, width = 6)
        print(enhanced_volcano_1_input())
        dev.off()
      }
    )
    
    output$downloadVolcano_enh2 <- downloadHandler(
      filename = function() {
        paste0("Enh_Volcano_", input$volcano_cntrst_enh2, "_", input$p_adj_enh2,".pdf")
      },
      content = function(file) {
        cairo_pdf(file, height = 8, width = 6)
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
    
    
    
  }) #End of Analyse button event
  
  
})
#End of Server

shinyApp(ui = ui, server = server)
