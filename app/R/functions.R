#DEP+Heatmaps_all.R by Stuart Cain
#The University of Manchester
#This script imports Maxquant output file proteinGroups.txt and combines it with a experimental
#design file. The DEP package is used to calculate foldchange values, which can be filtered using Cutoff values.
#The complexHeatmap package is used to render heatmaps, and the clusters are then analysed by clusterprofiler

Version = 2.4

#library(DEP)
#library(SummarizedExperiment)
#library(tidyverse) #contains: ggplot2, dplyr, tidyr, 
## readr, purr, tibble, stringr, forcats
#library(ComplexHeatmap)
#library(circlize)
#library(data.table)
#library(clusterProfiler)
#library(enrichplot)
#library(scales)
#library(gridtext)
#library(EnhancedVolcano)
#library("RColorBrewer")


set.seed(123) #random seed to give consistent kmeans clustering

####Colour functions####
colour_functions <- function( limit_comp=6, high_col="red", low_col="blue", palette_choice="default", alpha_input=10 ) {
  col_fun = colorRamp2(c(-limit_comp, 0, limit_comp), c(low_col, "white", high_col))
  col_fun(seq(-limit_comp, limit_comp))
  
  #Default Colours
  c6h <- c(
    "1"="#DF536B", # Maroon
    "2"="#61D04F", # Green
    "3"="#2297E6", # Dark Blue
    "4"="#28E2E5", # Light Blue
    "5"="#CD0BBC", # Purple
    "6"="#F5C710", # Yellow
    "7"="#9E9E9E", # grey
    "8"="#539002"  # dark green 
  )
  
  if (palette_choice=="default"){
    mypal<-c6h
  }else{
    mypal = brewer.pal(n = 8, name = palette_choice)
  }
  
  if (palette_choice=="Blues"){
    mypal = brewer.pal(n = 12, name = palette_choice)
    mypal <- mypal[3:10]
  }
  if (palette_choice=="Reds"){
    mypal = brewer.pal(n = 12, name = palette_choice)
    mypal <- mypal[3:10]
  }
  if (palette_choice=="Greens"){
    mypal = brewer.pal(n = 12, name = palette_choice)
    mypal <- mypal[3:10]
  }
  if (palette_choice=="Greys"){
    mypal = brewer.pal(n = 12, name = palette_choice)
    mypal <- mypal[3:10]
  }
  if(alpha_input<10){
    mypal<-adjustcolor(mypal, alpha.f = alpha_input/10)
  }
  
  col_fun_out<-list(col_fun, mypal, alpha_input)
  
  return(col_fun_out)
} 


####Complex Heatmap function####
plot_complex_heatmap <- function(data_results, kmeans, ht_title, col_fun, mypal, switch=3, font_size=6, cluster_cols=TRUE, comp_list) {
  set.seed(123)
  FC_data <- data_results[, c(1,grep("_ratio",colnames(data_results)))]
  #Change column names
  for ( col in 1:ncol(FC_data)){
    colnames(FC_data)[col] <-  sub("_ratio", "", colnames(FC_data)[col])
  }
  FC_matrix <- FC_data %>% remove_rownames %>% column_to_rownames(var="name")
  
 
  
  #Set up label data
  if (switch == 1){
    label = "Imputed_Int_Sig"
  }
  if (switch == 2){
    label = "Orig_Int_Sig"
  }
  if (switch == 3) {
    label = "Contrast"
  } 
  

  ####Set up Heatmaps####
  #Choose which set ofdata to use (Intensities or Log2FC)
  if (switch == 3) {
    #  my_mat <-
    #    as.matrix(dfdata2) #converts table into matrix called my_mat
    my_mat<-FC_matrix
    my_mat <- as.matrix(my_mat)
    if (ncol(my_mat)>1){my_mat<-my_mat[,comp_list]}
    CS = colMeans(FC_matrix) #creates a mean of each column
    col_fun = col_fun
    kmeans<-kmeans
  } 
  if (switch == 1) {
    my_mat <- as.matrix(imp_matrix_filt)
    #my_mat <- log2(my_mat)
    CSr = rowMeans(my_mat)
    my_mat <- my_mat - CSr
    CSi = colMeans(my_mat)
    col_fun = col_int
    kmeans<-kmeans
  }
  
  if (switch == 2) {
    my_mat <- as.matrix(norm_matrix_filt)
    #my_mat <- log2(my_mat)
    CSr = rowMeans(my_mat,na.rm=TRUE)
    my_mat <- my_mat - CSr
    my_mat[is.na(my_mat)] = 0
    CSi = colMeans(my_mat)
    col_fun = col_int
    kmeans<-kmeans
  }
  
  #Set Proteins count
  LFQ_GL_Count<-nrow(my_mat)
  
  
  ####Generate Annotated Heatmap####

  CS= colMeans(my_mat) #creates a mean of each column
  
  #Set up Heatmap annotations
  if (switch == 3) {
    ha = HeatmapAnnotation(
      log2FC = anno_density(
        my_mat,
        type = "violin",
        height = unit(2, "cm"),
        gp = gpar(fill = mypal[c(1:8)])
      ),
      FC_mean = CS,
      col = list(FC_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend=FALSE
    )
    htname = "log2FC" 
  }else{
    ha = HeatmapAnnotation(
      log2Int = anno_density(
        my_mat,
        type = "violin",
        height = unit(2, "cm"),
        gp = gpar(fill = mypal[c(1:8)])
      ),
      Int_mean = CS,
      col = list(Int_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend=FALSE
    )
    htname = "Log2 Int."
  }
  
  #Right Annotation
  #For RA
  dfdata_sig<-data_results[, c(1,grep("_significant",colnames(data_results)))]
  dfdata_sig2 <- dfdata_sig%>% mutate(sigcount = rowSums(dfdata_sig[-1], na.rm = TRUE))
  #saint_count3 = setDT(dfdata_sig2)[name %chin% LFQ_GL]
  sc3max<- max(dfdata_sig2$sigcount)
  
  ra = rowAnnotation("Sig.
Freq" = anno_barplot(
  dfdata_sig2$sigcount,
  ylim = c(0 , sc3max),
  gp = gpar(fill = mypal[1], fontsize = 8)))
  
  
  
  #Set up left annotation to match kmeans number
  if (kmeans == 1) {
    la = rowAnnotation(foo = anno_block(
      gp = gpar(fill = mypal[c(1:8)]),
      labels = c("Cluster 1"),
      labels_gp = gpar(col = "white", fontsize = font_size+4)
    ))
  } else{
    if (kmeans == 2) {
      la = rowAnnotation(foo = anno_block(
        gp = gpar(fill = mypal[c(1:8)]),
        labels = c("Cluster 1", "Cluster 2"),
        labels_gp = gpar(col = "white", fontsize = font_size+4)
      ))
    } else{
      if (kmeans == 3) {
        la = rowAnnotation(foo = anno_block(
          gp = gpar(fill = mypal[c(1:8)]),
          labels = c("Cluster 1", "Cluster 2" , "Cluster 3"),
          labels_gp = gpar(col = "white", fontsize = font_size+4)
        ))
      } else{
        if (kmeans == 4) {
          la = rowAnnotation(foo = anno_block(
            gp = gpar(fill = mypal[c(1:8)]),
            labels = c("Cluster 1", "Cluster 2" , "Cluster 3", "Cluster 4"),
            labels_gp = gpar(col = "white", fontsize = font_size+4)
          ))
        } else{
          if (kmeans == 5) {
            la = rowAnnotation(foo = anno_block(
              gp = gpar(fill = mypal[c(1:8)]),
              labels = c(
                "Cluster 1",
                "Cluster 2",
                "Cluster 3",
                "Cluster 4",
                "Cluster 5"
              ),
              labels_gp = gpar(col = "white", fontsize = font_size+4)
            ))
          } else{
            if (kmeans == 6) {
              la = rowAnnotation(foo = anno_block(
                gp = gpar(fill = mypal[c(1:8)]),
                labels = c(
                  "Cluster 1",
                  "Cluster 2",
                  "Cluster 3",
                  "Cluster 4",
                  "Cluster 5",
                  "Cluster 6"
                ),
                labels_gp = gpar(col = "white", fontsize = font_size+4)
              ))
            } else{
              if (kmeans > 6) {
                la = rowAnnotation(" " = anno_empty(
                  border = FALSE,
                  zoom = FALSE,
                  width = unit(0.01, "cm")
                ))
              }
            }
          }
        }
      }
    }
  }
  #ht_title<-paste("FC=", FC, " pVal=", pVal, "<br>", "(",LFQ_GL_Count, " Proteins)", "<br>",kmeans," Clusters ",label , sep="")
  
  ht = Heatmap(my_mat,
               name = htname,
               cluster_rows = TRUE, 
               cluster_columns = cluster_cols,
               column_dend_side = "bottom",
               top_annotation = ha,
               left_annotation = la,
               right_annotation = ra,
               #bottom_annotation = ba,
               row_km = kmeans, row_km_repeats = 2000, # kmeans clustering
               cluster_row_slices = FALSE,
               #column_km = 5, column_km_repeats = 2000, # kmeans clustering
               col = col_fun,
               #column_split = column_split,
               border = TRUE,
               row_names_gp = gpar(fontsize = font_size),
               column_names_gp = gpar(fontsize = font_size+2),
               column_title = gt_render (ht_title))
  #ht = draw(ht)
  
  compht<-list(ht, my_mat, dfdata_sig2)
  
  return(compht)
  
}


plot_complex_heatmap_2xcomp <- function(data_results, kmeans, ht_title, col_fun, mypal, switch=3, font_size=6, cluster_cols=TRUE, comp_list, data_results1, data_results2) {
  set.seed(123)
  FC_data <- data_results[, c(1,grep("_ratio",colnames(data_results)))]
  #Change column names
  for ( col in 1:ncol(FC_data)){
    colnames(FC_data)[col] <-  sub("_ratio", "", colnames(FC_data)[col])
  }
  FC_matrix <- FC_data %>% remove_rownames %>% column_to_rownames(var="name")
  
  
  
  #Set up label data
  if (switch == 1){
    label = "Imputed_Int_Sig"
  }
  if (switch == 2){
    label = "Orig_Int_Sig"
  }
  if (switch == 3) {
    label = "Contrast"
  } 
  
  
  ####Set up Heatmaps####
  #Choose which set ofdata to use (Intensities or Log2FC)
  if (switch == 3) {
    #  my_mat <-
    #    as.matrix(dfdata2) #converts table into matrix called my_mat
    my_mat<-FC_matrix
    my_mat <- as.matrix(my_mat)
    if (ncol(my_mat)>1){my_mat<-my_mat[,comp_list]}
    CS = colMeans(FC_matrix) #creates a mean of each column
    col_fun = col_fun
    kmeans<-kmeans
  } 
  if (switch == 1) {
    my_mat <- as.matrix(imp_matrix_filt)
    #my_mat <- log2(my_mat)
    CSr = rowMeans(my_mat)
    my_mat <- my_mat - CSr
    CSi = colMeans(my_mat)
    col_fun = col_int
    kmeans<-kmeans
  }
  
  if (switch == 2) {
    my_mat <- as.matrix(norm_matrix_filt)
    #my_mat <- log2(my_mat)
    CSr = rowMeans(my_mat,na.rm=TRUE)
    my_mat <- my_mat - CSr
    my_mat[is.na(my_mat)] = 0
    CSi = colMeans(my_mat)
    col_fun = col_int
    kmeans<-kmeans
  }
  
  #Set Proteins count
  LFQ_GL_Count<-nrow(my_mat)
  
  
  ####Generate Annotated Heatmap####
  
  CS= colMeans(my_mat) #creates a mean of each column
  
  #Set up Heatmap annotations
  if (switch == 3) {
    ha = HeatmapAnnotation(
      log2FC = anno_density(
        my_mat,
        type = "violin",
        height = unit(2, "cm"),
        gp = gpar(fill = mypal[c(1:8)])
      ),
      FC_mean = CS,
      col = list(FC_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend=FALSE
    )
    htname = "log2FC" 
  }else{
    ha = HeatmapAnnotation(
      log2Int = anno_density(
        my_mat,
        type = "violin",
        height = unit(2, "cm"),
        gp = gpar(fill = mypal[c(1:8)])
      ),
      Int_mean = CS,
      col = list(Int_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend=FALSE
    )
    htname = "Log2 Int."
  }
  
  #Right Annotation
  #For RA
  dfdata_sig<-data_results[, c(1,grep("_significant",colnames(data_results)))]
  dfdata_sig_c1<-data_results1[, c(1,grep("_significant",colnames(data_results1)))]
  dfdata_sig_c2<-data_results2[, c(1,grep("_significant",colnames(data_results2)))]
  
  #dfdata_sig2 <- dfdata_sig%>% mutate(sigcount = rowSums(dfdata_sig[-1], na.rm = TRUE))
  dfdata_sig2_c1 <- dfdata_sig_c1%>% mutate(sigcount_c1 = rowSums(dfdata_sig_c1[-1], na.rm = TRUE))
  dfdata_sig2_c2<- dfdata_sig_c2%>% mutate(sigcount_c2 = rowSums(dfdata_sig_c2[-1], na.rm = TRUE))
  
  dfdata_sig2=merge(dfdata_sig, dfdata_sig2_c1, by.x = "name", by.y = "name", all.x = TRUE, all.y = FALSE)
  dfdata_sig_c1m<-dfdata_sig2
  dfdata_sig_c1m=merge(dfdata_sig_c1m, dfdata_sig2_c2, by.x = "name", by.y = "name", all.x = TRUE, all.y = FALSE)
  dfdata_sig2<-dfdata_sig_c1m[, c("name", "sigcount_c1", "sigcount_c2")]
  dfdata_sig2[is.na(dfdata_sig2)] <- 0
  
  #saint_count3 = setDT(dfdata_sig2)[name %chin% LFQ_GL]
  sc3max_c1<- max(dfdata_sig2$sigcount_c1)
  sc3max_c2<- max(dfdata_sig2$sigcount_c2)
  
  ra_2 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  dfdata_sig2$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 8), 
)) +
    rowAnnotation("Sig.
Freq2" = anno_barplot(
  dfdata_sig2$sigcount_c2,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  

  
  
  #Set up left annotation to match kmeans number
  if (kmeans == 1) {
    la = rowAnnotation(foo = anno_block(
      gp = gpar(fill = mypal[c(1:8)]),
      labels = c("Cluster 1"),
      labels_gp = gpar(col = "white", fontsize = font_size+4)
    ))
  } else{
    if (kmeans == 2) {
      la = rowAnnotation(foo = anno_block(
        gp = gpar(fill = mypal[c(1:8)]),
        labels = c("Cluster 1", "Cluster 2"),
        labels_gp = gpar(col = "white", fontsize = font_size+4)
      ))
    } else{
      if (kmeans == 3) {
        la = rowAnnotation(foo = anno_block(
          gp = gpar(fill = mypal[c(1:8)]),
          labels = c("Cluster 1", "Cluster 2" , "Cluster 3"),
          labels_gp = gpar(col = "white", fontsize = font_size+4)
        ))
      } else{
        if (kmeans == 4) {
          la = rowAnnotation(foo = anno_block(
            gp = gpar(fill = mypal[c(1:8)]),
            labels = c("Cluster 1", "Cluster 2" , "Cluster 3", "Cluster 4"),
            labels_gp = gpar(col = "white", fontsize = font_size+4)
          ))
        } else{
          if (kmeans == 5) {
            la = rowAnnotation(foo = anno_block(
              gp = gpar(fill = mypal[c(1:8)]),
              labels = c(
                "Cluster 1",
                "Cluster 2",
                "Cluster 3",
                "Cluster 4",
                "Cluster 5"
              ),
              labels_gp = gpar(col = "white", fontsize = font_size+4)
            ))
          } else{
            if (kmeans == 6) {
              la = rowAnnotation(foo = anno_block(
                gp = gpar(fill = mypal[c(1:8)]),
                labels = c(
                  "Cluster 1",
                  "Cluster 2",
                  "Cluster 3",
                  "Cluster 4",
                  "Cluster 5",
                  "Cluster 6"
                ),
                labels_gp = gpar(col = "white", fontsize = font_size+4)
              ))
            } else{
              if (kmeans > 6) {
                la = rowAnnotation(" " = anno_empty(
                  border = FALSE,
                  zoom = FALSE,
                  width = unit(0.01, "cm")
                ))
              }
            }
          }
        }
      }
    }
  }
  #ht_title<-paste("FC=", FC, " pVal=", pVal, "<br>", "(",LFQ_GL_Count, " Proteins)", "<br>",kmeans," Clusters ",label , sep="")
  
  ht = Heatmap(my_mat,
               name = htname,
               cluster_rows = TRUE, 
               cluster_columns = cluster_cols,
               column_dend_side = "bottom",
               top_annotation = ha,
               left_annotation = la,
               #right_annotation = ra,
               #bottom_annotation = ba,
               row_km = kmeans, row_km_repeats = 2000, # kmeans clustering
               cluster_row_slices = FALSE,
               #column_km = 5, column_km_repeats = 2000, # kmeans clustering
               col = col_fun,
               #column_split = column_split,
               border = TRUE,
               row_names_gp = gpar(fontsize = font_size),
               column_names_gp = gpar(fontsize = font_size+2),
               column_title = gt_render (ht_title))
  
  ht = draw(ht+ra_2, auto_adjust = FALSE)
  
  
  compht<-list(ht, my_mat, dfdata_sig2)
  
  return(compht)
  
}

#BFDR Heatmap
plot_complex_heatmap_BFDR <- function(data_results, kmeans, ht_title, col_fun, mypal, switch=3, font_size=6, cluster_cols=TRUE, comp_list, data_results1, saintlist) {
  set.seed(123)
  FC_data <- data_results[, c(1,grep("_ratio",colnames(data_results)))]
  #Change column names
  for ( col in 1:ncol(FC_data)){
    colnames(FC_data)[col] <-  sub("_ratio", "", colnames(FC_data)[col])
  }
  FC_matrix <- FC_data %>% remove_rownames %>% column_to_rownames(var="name")
  
  
  
  #Set up label data
  if (switch == 1){
    label = "Imputed_Int_Sig"
  }
  if (switch == 2){
    label = "Orig_Int_Sig"
  }
  if (switch == 3) {
    label = "Contrast"
  } 
  
  
  ####Set up Heatmaps####
  #Choose which set ofdata to use (Intensities or Log2FC)
  if (switch == 3) {
    #  my_mat <-
    #    as.matrix(dfdata2) #converts table into matrix called my_mat
    my_mat<-FC_matrix
    my_mat <- as.matrix(my_mat)
    if (ncol(my_mat)>1){my_mat<-my_mat[,comp_list]}
    CS = colMeans(FC_matrix) #creates a mean of each column
    col_fun = col_fun
    kmeans<-kmeans
  } 
  if (switch == 1) {
    my_mat <- as.matrix(imp_matrix_filt)
    #my_mat <- log2(my_mat)
    CSr = rowMeans(my_mat)
    my_mat <- my_mat - CSr
    CSi = colMeans(my_mat)
    col_fun = col_int
    kmeans<-kmeans
  }
  
  if (switch == 2) {
    my_mat <- as.matrix(norm_matrix_filt)
    #my_mat <- log2(my_mat)
    CSr = rowMeans(my_mat,na.rm=TRUE)
    my_mat <- my_mat - CSr
    my_mat[is.na(my_mat)] = 0
    CSi = colMeans(my_mat)
    col_fun = col_int
    kmeans<-kmeans
  }
  
  #Set Proteins count
  LFQ_GL_Count<-nrow(my_mat)
  
  
  ####Generate Annotated Heatmap####
  
  CS= colMeans(my_mat) #creates a mean of each column
  
  #Set up Heatmap annotations
  if (switch == 3) {
    ha = HeatmapAnnotation(
      log2FC = anno_density(
        my_mat,
        type = "violin",
        height = unit(2, "cm"),
        gp = gpar(fill = mypal[c(1:8)])
      ),
      FC_mean = CS,
      col = list(FC_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend=FALSE
    )
    htname = "log2FC" 
  }else{
    ha = HeatmapAnnotation(
      log2Int = anno_density(
        my_mat,
        type = "violin",
        height = unit(2, "cm"),
        gp = gpar(fill = mypal[c(1:8)])
      ),
      Int_mean = CS,
      col = list(Int_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend=FALSE
    )
    htname = "Log2 Int."
  }
  
  #Right Annotation
  #For RA
  dfdata_sig<-data_results[, c(1,grep("_significant",colnames(data_results)))]
  dfdata_sig_c1<-data_results1[, c(1,grep("_significant",colnames(data_results1)))]
  #dfdata_sig_c2<-saintlist
  
  #dfdata_sig2 <- dfdata_sig%>% mutate(sigcount = rowSums(dfdata_sig[-1], na.rm = TRUE))
  dfdata_sig2_c1 <- dfdata_sig_c1%>% mutate(sigcount_c1 = rowSums(dfdata_sig_c1[-1], na.rm = TRUE))
  dfdata_sig2_c2<-saintlist 
  
  dfdata_sig2=merge(dfdata_sig, dfdata_sig2_c1, by.x = "name", by.y = "name", all.x = TRUE, all.y = FALSE)
  dfdata_sig_c1m<-dfdata_sig2
  dfdata_sig_c1m=merge(dfdata_sig_c1m, dfdata_sig2_c2, by.x = "name", by.y = "PreyGene", all.x = TRUE, all.y = FALSE)
  dfdata_sig2<-dfdata_sig_c1m[, c("name", "sigcount_c1", "BFDR_Freq")]
  dfdata_sig2[is.na(dfdata_sig2)] <- 0
  
  #saint_count3 = setDT(dfdata_sig2)[name %chin% LFQ_GL]
  sc3max_c1<- max(dfdata_sig2$sigcount_c1)
  sc3max_c2<- max(dfdata_sig2$BFDR_Freq)
  
  ra_2 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  dfdata_sig2$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 8), 
)) +
    rowAnnotation("BFDR
Freq" = anno_barplot(
  dfdata_sig2$BFDR_Freq,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  
  
  
  
  #Set up left annotation to match kmeans number
  if (kmeans == 1) {
    la = rowAnnotation(foo = anno_block(
      gp = gpar(fill = mypal[c(1:8)]),
      labels = c("Cluster 1"),
      labels_gp = gpar(col = "white", fontsize = font_size+4)
    ))
  } else{
    if (kmeans == 2) {
      la = rowAnnotation(foo = anno_block(
        gp = gpar(fill = mypal[c(1:8)]),
        labels = c("Cluster 1", "Cluster 2"),
        labels_gp = gpar(col = "white", fontsize = font_size+4)
      ))
    } else{
      if (kmeans == 3) {
        la = rowAnnotation(foo = anno_block(
          gp = gpar(fill = mypal[c(1:8)]),
          labels = c("Cluster 1", "Cluster 2" , "Cluster 3"),
          labels_gp = gpar(col = "white", fontsize = font_size+4)
        ))
      } else{
        if (kmeans == 4) {
          la = rowAnnotation(foo = anno_block(
            gp = gpar(fill = mypal[c(1:8)]),
            labels = c("Cluster 1", "Cluster 2" , "Cluster 3", "Cluster 4"),
            labels_gp = gpar(col = "white", fontsize = font_size+4)
          ))
        } else{
          if (kmeans == 5) {
            la = rowAnnotation(foo = anno_block(
              gp = gpar(fill = mypal[c(1:8)]),
              labels = c(
                "Cluster 1",
                "Cluster 2",
                "Cluster 3",
                "Cluster 4",
                "Cluster 5"
              ),
              labels_gp = gpar(col = "white", fontsize = font_size+4)
            ))
          } else{
            if (kmeans == 6) {
              la = rowAnnotation(foo = anno_block(
                gp = gpar(fill = mypal[c(1:8)]),
                labels = c(
                  "Cluster 1",
                  "Cluster 2",
                  "Cluster 3",
                  "Cluster 4",
                  "Cluster 5",
                  "Cluster 6"
                ),
                labels_gp = gpar(col = "white", fontsize = font_size+4)
              ))
            } else{
              if (kmeans > 6) {
                la = rowAnnotation(" " = anno_empty(
                  border = FALSE,
                  zoom = FALSE,
                  width = unit(0.01, "cm")
                ))
              }
            }
          }
        }
      }
    }
  }
  #ht_title<-paste("FC=", FC, " pVal=", pVal, "<br>", "(",LFQ_GL_Count, " Proteins)", "<br>",kmeans," Clusters ",label , sep="")
  
  ht = Heatmap(my_mat,
               name = htname,
               cluster_rows = TRUE, 
               cluster_columns = cluster_cols,
               column_dend_side = "bottom",
               top_annotation = ha,
               left_annotation = la,
               #right_annotation = ra,
               #bottom_annotation = ba,
               row_km = kmeans, row_km_repeats = 2000, # kmeans clustering
               cluster_row_slices = FALSE,
               #column_km = 5, column_km_repeats = 2000, # kmeans clustering
               col = col_fun,
               #column_split = column_split,
               border = TRUE,
               row_names_gp = gpar(fontsize = font_size),
               column_names_gp = gpar(fontsize = font_size+2),
               column_title = gt_render (ht_title))
  
  ht = draw(ht+ra_2, auto_adjust = FALSE)
  
  
  compht<-list(ht, my_mat, dfdata_sig2)
  
  return(compht)
  
}


####Heatmap Clusters Function####
plot_heatmap_clusters <- function(ht, my_mat, dfdata_sig2, kmeans, FC=2, pVal=0.05, col_fun, mypal, switch=3, font_size=6, data_unique_defn){
  set.seed(123)
  #Set up label data
  if (switch == 1){
    label = "Imputed_Int_Sig"
  }
  if (switch == 2){
    label = "Orig_Int_Sig"
  }
  if (switch == 3) {
    label = "Contrast"
  } 
  
  ht_C1<-"No Heatmap"
  ht_C2<-"No Heatmap"
  ht_C3<-"No Heatmap"
  ht_C4<-"No Heatmap"
  ht_C5<-"No Heatmap"
  ht_C6<-"No Heatmap"
  
  comp_length<-ncol(my_mat)
  
  #Reorder columns to match clustering
  if(comp_length>1){
    co = column_order(ht)
    my_mat<-my_mat[,co]
    }
  
  #Split into the Clusters
  
  #Obtain Gene names from each kmeans cluster
  if (kmeans>1){
    ro = row_order(ht)
    
    if (comp_length>1)  {clu_df <- lapply(names(ro), function(i) {
      out <- data.frame(
        GeneID = rownames(my_mat[ro[[i]], ]),
        Cluster = paste0("cluster", i),
        stringsAsFactors = FALSE
      )
      return(out)
    }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
      do.call(rbind, .)
    }
    
    if (comp_length== 1){
      my_mat_r <-cbind(my_mat, rep(my_mat[,1],1))
      
      ht_r = Heatmap(my_mat_r,
                     #name = htname,
                     cluster_rows = TRUE, 
                     cluster_columns = TRUE,
                     column_dend_side = "bottom",
                     #top_annotation = ha,
                     #left_annotation = la,
                     #right_annotation = ra,
                     #bottom_annotation = ba,
                     row_km = kmeans, row_km_repeats = 2000, # kmeans clustering
                     cluster_row_slices = FALSE,
                     #column_km = 5, column_km_repeats = 2000, # kmeans clustering
                     col = col_fun,
                     #column_split = column_split,
                     border = TRUE,
                     row_names_gp = gpar(fontsize = 5),
                     column_names_gp = gpar(fontsize = 8)
                     #column_title = gt_render (hti_title)
      )
      
      ht_r = draw(ht_r)  
      
      ro = row_order(ht_r)    
      clu_df <- lapply(names(ro), function(i) {
        out <- data.frame(
          GeneID = rownames(my_mat_r[ro[[i]], ]),
          Cluster = paste0("cluster", i),
          stringsAsFactors = FALSE
        )
        return(out)
      }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
        do.call(rbind, .)
    }
  }
  
  if(kmeans==1){
    clu_df_new <- data.frame(GeneID = row.names(my_mat), my_mat)
    clu_df_new$Cluster <- "cluster1"
    clu_df <- clu_df_new[, c("GeneID", "Cluster")]
  }
  
  #check
  clu_df
  

  #Create filtered lists
  Cluster1 =clu_df[clu_df$Cluster == "cluster1", ]
  Cluster2 =clu_df[clu_df$Cluster == "cluster2", ]
  Cluster3 =clu_df[clu_df$Cluster == "cluster3", ]
  Cluster4 =clu_df[clu_df$Cluster == "cluster4", ]
  Cluster5 =clu_df[clu_df$Cluster == "cluster5", ]
  Cluster6 =clu_df[clu_df$Cluster == "cluster6", ]
  
  #Count each cluster
  C1_Count<-nrow(Cluster1)
  C2_Count<-nrow(Cluster2)
  C3_Count<-nrow(Cluster3)
  C4_Count<-nrow(Cluster4)
  C5_Count<-nrow(Cluster5)
  C6_Count<-nrow(Cluster6)
  
  #Generate separate cluster matrices
  if(comp_length==1){
    if (switch ==3){
      my_mat_r <-cbind(my_mat, rep(my_mat[,1],1))
      
      my_mat_c1 <- my_mat_r [rownames(my_mat_r) %in% Cluster1$GeneID, ]
      my_mat_c2 <- my_mat_r [rownames(my_mat_r) %in% Cluster2$GeneID, ]
      my_mat_c3 <- my_mat_r [rownames(my_mat_r) %in% Cluster3$GeneID, ]
      my_mat_c4 <- my_mat_r [rownames(my_mat_r) %in% Cluster4$GeneID, ]
      my_mat_c5 <- my_mat_r [rownames(my_mat_r) %in% Cluster5$GeneID, ]
      my_mat_c6 <- my_mat_r [rownames(my_mat_r) %in% Cluster6$GeneID, ]
      my_mat_c1 <- data.frame(name = row.names(my_mat_c1), my_mat_c1)
      my_mat_c2 <- data.frame(name = row.names(my_mat_c2), my_mat_c2)
      my_mat_c3 <- data.frame(name = row.names(my_mat_c3), my_mat_c3)
      my_mat_c4 <- data.frame(name = row.names(my_mat_c4), my_mat_c4)
      my_mat_c5 <- data.frame(name = row.names(my_mat_c5), my_mat_c5)
      my_mat_c6 <- data.frame(name = row.names(my_mat_c6), my_mat_c6)
      my_mat_c1 <- my_mat_c1[,-3]
      my_mat_c2 <- my_mat_c2[,-3]
      my_mat_c3 <- my_mat_c3[,-3]
      my_mat_c4 <- my_mat_c4[,-3]
      my_mat_c5 <- my_mat_c5[,-3]
      my_mat_c6 <- my_mat_c6[,-3]
      
      if (kmeans>0){my_mat_c1 <- my_mat_c1 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>1){my_mat_c2 <- my_mat_c2 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>2){my_mat_c3 <- my_mat_c3 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>3){my_mat_c4 <- my_mat_c4 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>4){my_mat_c5 <- my_mat_c5 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>5){my_mat_c6 <- my_mat_c6 %>% remove_rownames %>% column_to_rownames(var="name")}
      
    }else{
      my_mat_c1 <- my_mat [rownames(my_mat) %in% Cluster1$GeneID, ]
      my_mat_c2 <- my_mat [rownames(my_mat) %in% Cluster2$GeneID, ]
      my_mat_c3 <- my_mat [rownames(my_mat) %in% Cluster3$GeneID, ]
      my_mat_c4 <- my_mat [rownames(my_mat) %in% Cluster4$GeneID, ]
      my_mat_c5 <- my_mat [rownames(my_mat) %in% Cluster5$GeneID, ]
      my_mat_c6 <- my_mat [rownames(my_mat) %in% Cluster6$GeneID, ]
    }
    
  }
  
  
  if(comp_length>1){
    my_mat_c1 <- my_mat [rownames(my_mat) %in% Cluster1$GeneID, ]
    my_mat_c2 <- my_mat [rownames(my_mat) %in% Cluster2$GeneID, ]
    my_mat_c3 <- my_mat [rownames(my_mat) %in% Cluster3$GeneID, ]
    my_mat_c4 <- my_mat [rownames(my_mat) %in% Cluster4$GeneID, ]
    my_mat_c5 <- my_mat [rownames(my_mat) %in% Cluster5$GeneID, ]
    my_mat_c6 <- my_mat [rownames(my_mat) %in% Cluster6$GeneID, ]
  }
  
  CS1= colMeans(my_mat_c1)
  CS2= colMeans(my_mat_c2)
  CS3= colMeans(my_mat_c3)
  CS4= colMeans(my_mat_c4)
  CS5= colMeans(my_mat_c5)
  CS6= colMeans(my_mat_c6)
  
  #Generate Cluster txt files
  Cluster1FC <- data.frame(name = row.names(my_mat_c1), my_mat_c1)
  Cluster2FC <- data.frame(name = row.names(my_mat_c2), my_mat_c2)
  Cluster3FC <- data.frame(name = row.names(my_mat_c3), my_mat_c3)
  Cluster4FC <- data.frame(name = row.names(my_mat_c4), my_mat_c4)
  Cluster5FC <- data.frame(name = row.names(my_mat_c5), my_mat_c5)
  Cluster6FC <- data.frame(name = row.names(my_mat_c6), my_mat_c6)
  
  if (kmeans>0){Cluster1FC <- merge(data_unique_defn, Cluster1FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>1){Cluster2FC <- merge(data_unique_defn, Cluster2FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>2){Cluster3FC <- merge(data_unique_defn, Cluster3FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>3){Cluster4FC <- merge(data_unique_defn, Cluster4FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>4){Cluster5FC <- merge(data_unique_defn, Cluster5FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>5){Cluster6FC <- merge(data_unique_defn, Cluster6FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  
  
  # Converts each cluster into a list of values
  dfc1 <- Cluster1$GeneID 
  dfc2 <- Cluster2$GeneID
  dfc3 <- Cluster3$GeneID
  dfc4 <- Cluster4$GeneID
  dfc5 <- Cluster5$GeneID
  dfc6 <- Cluster6$GeneID
  
  
  ####Generate Separate Cluster Heatmaps####
  #Generate Annotations
  #Top heatmap annotation
  if (switch == 3) {
    if (kmeans > 0) {
      ha1 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c1,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS1,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 1) {
      ha2 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c2,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS2,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 2) {
      ha3 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c3,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS3,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 3) {
      ha4 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c4,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS4,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 4) {
      ha5 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c5,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS5,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 5) {
      ha6 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c6,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS6,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    htname = "log2FC"
  } else{
    if (kmeans > 0) {
      ha1 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c1,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS1,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 1) {
      ha2 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c2,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS2,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 2) {
      ha3 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c3,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS3,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 3) {
      ha4 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c4,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS4,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 4) {
      ha5 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c5,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS5,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 5) {
      ha6 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c6,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS6,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    htname = "Log2 Int."
  }
  
  
  la1 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[1]),
                                       labels = c("Cluster 1"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la2 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[2]),
                                       labels = c("Cluster 2"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la3 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[3]),
                                       labels = c("Cluster 3"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la4 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[4]),
                                       labels = c("Cluster 4"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la5 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[5]),
                                       labels = c("Cluster 5"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la6 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[6]),
                                       labels = c("Cluster 6"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  
  #ba = HeatmapAnnotation(foo2 = anno_block(gp = gpar(fill = 5:7),
  #                                         labels = c("2kPa vs 25kPa", "2kPa vs 8kPa", "8kPa vs 25kPa"),
  #                                         labels_gp = gpar(col = "white", fontsize = 10)))
  #
  
  #Calculate Right annotations
  #if (switch==FALSE) {
  saint_count3 <- dfdata_sig2
  sc3max<- max(dfdata_sig2$sigcount)
  saint_count_cl1 <- saint_count3 %>% filter(saint_count3$name %in% Cluster1$GeneID )
  saint_count_cl2 <- saint_count3 %>% filter(saint_count3$name %in% Cluster2$GeneID )
  saint_count_cl3 <- saint_count3 %>% filter(saint_count3$name %in% Cluster3$GeneID )
  saint_count_cl4 <- saint_count3 %>% filter(saint_count3$name %in% Cluster4$GeneID )
  saint_count_cl5 <- saint_count3 %>% filter(saint_count3$name %in% Cluster5$GeneID )
  saint_count_cl6 <- saint_count3 %>% filter(saint_count3$name %in% Cluster6$GeneID )
  #
  #}
  
  #if (switch==FALSE) {
  ra1=rowAnnotation("Sig.
Freq" = anno_barplot(saint_count_cl1$sigcount, ylim = c( 0 , sc3max ), gp = gpar(fill = mypal[1], fontsize=8)))
  ra2=rowAnnotation("Sig.
Freq" = anno_barplot(saint_count_cl2$sigcount, ylim = c( 0 , sc3max ), gp = gpar(fill = mypal[1], fontsize=8)))
  ra3=rowAnnotation("Sig.
Freq" = anno_barplot(saint_count_cl3$sigcount, ylim = c( 0 , sc3max ), gp = gpar(fill = mypal[1], fontsize=8)))
  ra4=rowAnnotation("Sig.
Freq" = anno_barplot(saint_count_cl4$sigcount, ylim = c( 0 , sc3max ), gp = gpar(fill = mypal[1], fontsize=8)))
  ra5=rowAnnotation("Sig.
Freq" = anno_barplot(saint_count_cl5$sigcount, ylim = c( 0 , sc3max ), gp = gpar(fill = mypal[1], fontsize=8)))
  ra6=rowAnnotation("Sig.
Freq" = anno_barplot(saint_count_cl6$sigcount, ylim = c( 0 , sc3max ), gp = gpar(fill = mypal[1], fontsize=8)))
  
  
  if (kmeans > 0) {
    ht_C1 = Heatmap(
      my_mat_c1,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha1,
      left_annotation = la1,
      #bottom_annotation = ba,
      right_annotation = ra1,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (
        paste(
          label,
          " Cluster 1",
          " (",
          C1_Count,
          ")",
          sep = ""
        )
      )
    )
  }
  if (kmeans > 1) {
    ht_C2 = Heatmap(
      my_mat_c2,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha2,
      left_annotation = la2,
      #bottom_annotation = ba,
      right_annotation = ra2,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster 2", " (", C2_Count, ")", sep =
                              ""))
    )
  }
  if (kmeans > 2) {
    ht_C3 = Heatmap(
      my_mat_c3,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha3,
      left_annotation = la3,
      #bottom_annotation = ba,
      right_annotation = ra3,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster 3", " (", C3_Count, ")", sep =
                              ""))
    )
  }
  if (kmeans > 3) {
    ht_C4 = Heatmap(
      my_mat_c4,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha4,
      left_annotation = la4,
      #bottom_annotation = ba,
      right_annotation = ra4,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster 4", " (", C4_Count, ")", sep = ""))
    )
  }
  if (kmeans > 4) {
    ht_C5 = Heatmap(
      my_mat_c5,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha5,
      left_annotation = la5,
      #bottom_annotation = ba,
      right_annotation = ra5,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster 5", " (", C5_Count, ")", sep = ""))
    )
  } 
  if (kmeans > 5) {
    ht_C6 = Heatmap(
      my_mat_c6,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha6,
      left_annotation = la6,
      #bottom_annotation = ba,
      right_annotation = ra6,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster 6", " (", C6_Count, ")", sep = ""))
    )
  } 
  
  
  comphtclust<-list(ht_C1, ht_C2, ht_C3, ht_C4, ht_C5, ht_C6, C1_Count, C2_Count, C3_Count, C4_Count, C5_Count, C6_Count, 
                    Cluster1FC, Cluster2FC, Cluster3FC, Cluster4FC, Cluster5FC, Cluster6FC, clu_df)
  return(comphtclust)
}


plot_heatmap_clusters_2xcomp <- function(ht, my_mat, dfdata_sig2, kmeans, FC=2, pVal=0.05, col_fun, mypal, switch=3, font_size=6, data_unique_defn){
  set.seed(123)
  #Set up label data
  if (switch == 1){
    label = "Imputed_Int_Sig"
  }
  if (switch == 2){
    label = "Orig_Int_Sig"
  }
  if (switch == 3) {
    label = "Contrast"
  } 
  
  ht_C1<-"No Heatmap"
  ht_C2<-"No Heatmap"
  ht_C3<-"No Heatmap"
  ht_C4<-"No Heatmap"
  ht_C5<-"No Heatmap"
  ht_C6<-"No Heatmap"
  
  comp_length<-ncol(my_mat)
  
  #Reorder columns to match clustering
  if(comp_length>1){
    co = column_order(ht)
    my_mat<-my_mat[,co]
  }
  
  #Split into the Clusters
  
  #Obtain Gene names from each kmeans cluster
  if (kmeans>1){
    ro = row_order(ht)
    
    if (comp_length>1)  {clu_df <- lapply(names(ro), function(i) {
      out <- data.frame(
        GeneID = rownames(my_mat[ro[[i]], ]),
        Cluster = paste0("cluster", i),
        stringsAsFactors = FALSE
      )
      return(out)
    }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
      do.call(rbind, .)
    }
    
    if (comp_length== 1){
      my_mat_r <-cbind(my_mat, rep(my_mat[,1],1))
      
      ht_r = Heatmap(my_mat_r,
                     #name = htname,
                     cluster_rows = TRUE, 
                     cluster_columns = TRUE,
                     column_dend_side = "bottom",
                     #top_annotation = ha,
                     #left_annotation = la,
                     #right_annotation = ra,
                     #bottom_annotation = ba,
                     row_km = kmeans, row_km_repeats = 2000, # kmeans clustering
                     cluster_row_slices = FALSE,
                     #column_km = 5, column_km_repeats = 2000, # kmeans clustering
                     col = col_fun,
                     #column_split = column_split,
                     border = TRUE,
                     row_names_gp = gpar(fontsize = 5),
                     column_names_gp = gpar(fontsize = 8)
                     #column_title = gt_render (hti_title)
      )
      
      ht_r = draw(ht_r)  
      
      ro = row_order(ht_r)    
      clu_df <- lapply(names(ro), function(i) {
        out <- data.frame(
          GeneID = rownames(my_mat_r[ro[[i]], ]),
          Cluster = paste0("cluster", i),
          stringsAsFactors = FALSE
        )
        return(out)
      }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
        do.call(rbind, .)
    }
  }
  
  if(kmeans==1){
    clu_df_new <- data.frame(GeneID = row.names(my_mat), my_mat)
    clu_df_new$Cluster <- "cluster1"
    clu_df <- clu_df_new[, c("GeneID", "Cluster")]
  }
  
  #check
  clu_df
  
  
  #Create filtered lists
  Cluster1 =clu_df[clu_df$Cluster == "cluster1", ]
  Cluster2 =clu_df[clu_df$Cluster == "cluster2", ]
  Cluster3 =clu_df[clu_df$Cluster == "cluster3", ]
  Cluster4 =clu_df[clu_df$Cluster == "cluster4", ]
  Cluster5 =clu_df[clu_df$Cluster == "cluster5", ]
  Cluster6 =clu_df[clu_df$Cluster == "cluster6", ]
  
  #Count each cluster
  C1_Count<-nrow(Cluster1)
  C2_Count<-nrow(Cluster2)
  C3_Count<-nrow(Cluster3)
  C4_Count<-nrow(Cluster4)
  C5_Count<-nrow(Cluster5)
  C6_Count<-nrow(Cluster6)
  
  #Generate separate cluster matrices
  if(comp_length==1){
    if (switch ==3){
      my_mat_r <-cbind(my_mat, rep(my_mat[,1],1))
      
      my_mat_c1 <- my_mat_r [rownames(my_mat_r) %in% Cluster1$GeneID, ]
      my_mat_c2 <- my_mat_r [rownames(my_mat_r) %in% Cluster2$GeneID, ]
      my_mat_c3 <- my_mat_r [rownames(my_mat_r) %in% Cluster3$GeneID, ]
      my_mat_c4 <- my_mat_r [rownames(my_mat_r) %in% Cluster4$GeneID, ]
      my_mat_c5 <- my_mat_r [rownames(my_mat_r) %in% Cluster5$GeneID, ]
      my_mat_c6 <- my_mat_r [rownames(my_mat_r) %in% Cluster6$GeneID, ]
      my_mat_c1 <- data.frame(name = row.names(my_mat_c1), my_mat_c1)
      my_mat_c2 <- data.frame(name = row.names(my_mat_c2), my_mat_c2)
      my_mat_c3 <- data.frame(name = row.names(my_mat_c3), my_mat_c3)
      my_mat_c4 <- data.frame(name = row.names(my_mat_c4), my_mat_c4)
      my_mat_c5 <- data.frame(name = row.names(my_mat_c5), my_mat_c5)
      my_mat_c6 <- data.frame(name = row.names(my_mat_c6), my_mat_c6)
      my_mat_c1 <- my_mat_c1[,-3]
      my_mat_c2 <- my_mat_c2[,-3]
      my_mat_c3 <- my_mat_c3[,-3]
      my_mat_c4 <- my_mat_c4[,-3]
      my_mat_c5 <- my_mat_c5[,-3]
      my_mat_c6 <- my_mat_c6[,-3]
      
      if (kmeans>0){my_mat_c1 <- my_mat_c1 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>1){my_mat_c2 <- my_mat_c2 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>2){my_mat_c3 <- my_mat_c3 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>3){my_mat_c4 <- my_mat_c4 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>4){my_mat_c5 <- my_mat_c5 %>% remove_rownames %>% column_to_rownames(var="name")}
      if (kmeans>5){my_mat_c6 <- my_mat_c6 %>% remove_rownames %>% column_to_rownames(var="name")}
      
    }else{
      my_mat_c1 <- my_mat [rownames(my_mat) %in% Cluster1$GeneID, ]
      my_mat_c2 <- my_mat [rownames(my_mat) %in% Cluster2$GeneID, ]
      my_mat_c3 <- my_mat [rownames(my_mat) %in% Cluster3$GeneID, ]
      my_mat_c4 <- my_mat [rownames(my_mat) %in% Cluster4$GeneID, ]
      my_mat_c5 <- my_mat [rownames(my_mat) %in% Cluster5$GeneID, ]
      my_mat_c6 <- my_mat [rownames(my_mat) %in% Cluster6$GeneID, ]
    }
    
  }
  
  
  if(comp_length>1){
    my_mat_c1 <- my_mat [rownames(my_mat) %in% Cluster1$GeneID, ]
    my_mat_c2 <- my_mat [rownames(my_mat) %in% Cluster2$GeneID, ]
    my_mat_c3 <- my_mat [rownames(my_mat) %in% Cluster3$GeneID, ]
    my_mat_c4 <- my_mat [rownames(my_mat) %in% Cluster4$GeneID, ]
    my_mat_c5 <- my_mat [rownames(my_mat) %in% Cluster5$GeneID, ]
    my_mat_c6 <- my_mat [rownames(my_mat) %in% Cluster6$GeneID, ]
  }
  
  CS1= colMeans(my_mat_c1)
  CS2= colMeans(my_mat_c2)
  CS3= colMeans(my_mat_c3)
  CS4= colMeans(my_mat_c4)
  CS5= colMeans(my_mat_c5)
  CS6= colMeans(my_mat_c6)
  
  #Generate Cluster txt files
  Cluster1FC <- data.frame(name = row.names(my_mat_c1), my_mat_c1)
  Cluster2FC <- data.frame(name = row.names(my_mat_c2), my_mat_c2)
  Cluster3FC <- data.frame(name = row.names(my_mat_c3), my_mat_c3)
  Cluster4FC <- data.frame(name = row.names(my_mat_c4), my_mat_c4)
  Cluster5FC <- data.frame(name = row.names(my_mat_c5), my_mat_c5)
  Cluster6FC <- data.frame(name = row.names(my_mat_c6), my_mat_c6)
  
  if (kmeans>0){Cluster1FC <- merge(data_unique_defn, Cluster1FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>1){Cluster2FC <- merge(data_unique_defn, Cluster2FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>2){Cluster3FC <- merge(data_unique_defn, Cluster3FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>3){Cluster4FC <- merge(data_unique_defn, Cluster4FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>4){Cluster5FC <- merge(data_unique_defn, Cluster5FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  if (kmeans>5){Cluster6FC <- merge(data_unique_defn, Cluster6FC, by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)}
  
  
  # Converts each cluster into a list of values
  dfc1 <- Cluster1$GeneID 
  dfc2 <- Cluster2$GeneID
  dfc3 <- Cluster3$GeneID
  dfc4 <- Cluster4$GeneID
  dfc5 <- Cluster5$GeneID
  dfc6 <- Cluster6$GeneID
  
  
  ####Generate Separate Cluster Heatmaps####
  #Generate Annotations
  #Top heatmap annotation
  if (switch == 3) {
    if (kmeans > 0) {
      ha1 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c1,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS1,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 1) {
      ha2 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c2,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS2,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 2) {
      ha3 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c3,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS3,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 3) {
      ha4 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c4,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS4,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 4) {
      ha5 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c5,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS5,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 5) {
      ha6 = HeatmapAnnotation(
        log2FC = anno_density(
          my_mat_c6,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        FC_mean = CS6,
        col = list(FC_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    htname = "log2FC"
  } else{
    if (kmeans > 0) {
      ha1 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c1,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS1,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 1) {
      ha2 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c2,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS2,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 2) {
      ha3 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c3,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS3,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 3) {
      ha4 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c4,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS4,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 4) {
      ha5 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c5,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS5,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    if (kmeans > 5) {
      ha6 = HeatmapAnnotation(
        log2Int = anno_density(
          my_mat_c6,
          type = "violin",
          height = unit(3, "cm"),
          gp = gpar(fill = mypal[c(1:8)])
        ),
        Int_mean = CS6,
        col = list(Int_mean = col_fun),
        gp = gpar(col = "black"),
        show_legend = FALSE
      )
    }
    htname = "Log2 Int."
  }
  
  
  la1 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[1]),
                                       labels = c("Cluster 1"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la2 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[2]),
                                       labels = c("Cluster 2"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la3 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[3]),
                                       labels = c("Cluster 3"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la4 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[4]),
                                       labels = c("Cluster 4"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la5 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[5]),
                                       labels = c("Cluster 5"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  la6 = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[6]),
                                       labels = c("Cluster 6"), 
                                       labels_gp = gpar(col = "white", fontsize = font_size+4)))
  
  #ba = HeatmapAnnotation(foo2 = anno_block(gp = gpar(fill = 5:7),
  #                                         labels = c("2kPa vs 25kPa", "2kPa vs 8kPa", "8kPa vs 25kPa"),
  #                                         labels_gp = gpar(col = "white", fontsize = 10)))
  #
  
  #Calculate Right annotations
  #if (switch==FALSE) {
  saint_count3 <- dfdata_sig2
  sc3max_c1<- max(dfdata_sig2$sigcount_c1)
  sc3max_c2<- max(dfdata_sig2$sigcount_c2)
  saint_count_cl1 <- saint_count3 %>% filter(saint_count3$name %in% Cluster1$GeneID )
  saint_count_cl2 <- saint_count3 %>% filter(saint_count3$name %in% Cluster2$GeneID )
  saint_count_cl3 <- saint_count3 %>% filter(saint_count3$name %in% Cluster3$GeneID )
  saint_count_cl4 <- saint_count3 %>% filter(saint_count3$name %in% Cluster4$GeneID )
  saint_count_cl5 <- saint_count3 %>% filter(saint_count3$name %in% Cluster5$GeneID )
  saint_count_cl6 <- saint_count3 %>% filter(saint_count3$name %in% Cluster6$GeneID )
  #
  #}
  
  ra1 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  saint_count_cl1$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 6), 
)) +
    rowAnnotation("Sig.
Freq2" = anno_barplot(
  saint_count_cl1$sigcount_c2,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  ra2 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  saint_count_cl2$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 6), 
)) +
    rowAnnotation("Sig.
Freq2" = anno_barplot(
  saint_count_cl2$sigcount_c2,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  ra3 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  saint_count_cl3$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 6), 
)) +
    rowAnnotation("Sig.
Freq2" = anno_barplot(
  saint_count_cl3$sigcount_c2,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  ra4 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  saint_count_cl4$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 6), 
)) +
    rowAnnotation("Sig.
Freq2" = anno_barplot(
  saint_count_cl4$sigcount_c2,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  ra5 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  saint_count_cl5$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 6), 
)) +
    rowAnnotation("Sig.
Freq2" = anno_barplot(
  saint_count_cl5$sigcount_c2,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  
  ra6 =  rowAnnotation("Sig.
Freq1" = anno_barplot(
  saint_count_cl6$sigcount_c1,
  ylim = c(0 , sc3max_c1),
  gp = gpar(fill = mypal[1], fontsize = 6), 
)) +
    rowAnnotation("Sig.
Freq2" = anno_barplot(
  saint_count_cl6$sigcount_c2,
  ylim = c(0 , sc3max_c2),
  gp = gpar(fill = mypal[3], fontsize = 6), 
))
  
  
  if (kmeans > 0) {
    ht_C1 = Heatmap(
      my_mat_c1,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha1,
      left_annotation = la1,
      #bottom_annotation = ba,
      #right_annotation = ra1,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = (
        paste(
          label,
          " Cluster 1",
          " (",
          C1_Count,
          ")",
          sep = ""
        )
      )
    )
    ht_C1 = draw(ht_C1+ra1, auto_adjust = FALSE)
  }
  if (kmeans > 1) {
    ht_C2 = Heatmap(
      my_mat_c2,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha2,
      left_annotation = la2,
      #bottom_annotation = ba,
      #right_annotation = ra2,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = (paste("Cluster 2", " (", C2_Count, ")", sep =
                              ""))
    )
    ht_C2 = draw(ht_C2+ra2, auto_adjust = FALSE)
  }
  if (kmeans > 2) {
    ht_C3 = Heatmap(
      my_mat_c3,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha3,
      left_annotation = la3,
      #bottom_annotation = ba,
      #right_annotation = ra3,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = (paste("Cluster 3", " (", C3_Count, ")", sep =
                              ""))
    )
    ht_C3 = draw(ht_C3+ra3, auto_adjust = FALSE)
  }
  if (kmeans > 3) {
    ht_C4 = Heatmap(
      my_mat_c4,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha4,
      left_annotation = la4,
      #bottom_annotation = ba,
      #right_annotation = ra4,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = (paste("Cluster 4", " (", C4_Count, ")", sep = ""))
    )
    ht_C4 = draw(ht_C4+ra4, auto_adjust = FALSE)
  }
  if (kmeans > 4) {
    ht_C5 = Heatmap(
      my_mat_c5,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha5,
      left_annotation = la5,
      #bottom_annotation = ba,
      #right_annotation = ra5,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = (paste("Cluster 5", " (", C5_Count, ")", sep = ""))
    )
    ht_C5 = draw(ht_C5+ra5, auto_adjust = FALSE)
  } 
  if (kmeans > 5) {
    ht_C6 = Heatmap(
      my_mat_c6,
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha6,
      left_annotation = la6,
      #bottom_annotation = ba,
      #right_annotation = ra6,
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = (paste("Cluster 6", " (", C6_Count, ")", sep = ""))
    )
    ht_C6 = draw(ht_C6+ra6, auto_adjust = FALSE)
  } 
  
  
  comphtclust<-list(ht_C1, ht_C2, ht_C3, ht_C4, ht_C5, ht_C6, C1_Count, C2_Count, C3_Count, C4_Count, C5_Count, C6_Count, 
                    Cluster1FC, Cluster2FC, Cluster3FC, Cluster4FC, Cluster5FC, Cluster6FC, clu_df)
  return(comphtclust)
}



####Enhanced Volcano Function####
#choices ("p-adj sig.","p-Val sig.", "p-adj values")
plot_enhanced_volcano <- function(data_results, clu_df, comp_pick,  FC=2, pVal=0.05, padj_check="p-adj sig.", mypal, alpha_input=10){
  
  
  Cluster1 =clu_df[clu_df$Cluster == "cluster1", ]
  Cluster2 =clu_df[clu_df$Cluster == "cluster2", ]
  Cluster3 =clu_df[clu_df$Cluster == "cluster3", ]
  Cluster4 =clu_df[clu_df$Cluster == "cluster4", ]
  Cluster5 =clu_df[clu_df$Cluster == "cluster5", ]
  Cluster6 =clu_df[clu_df$Cluster == "cluster6", ]
  
  # Converts each cluster into a list of values
  dfc1 <- Cluster1$GeneID 
  dfc2 <- Cluster2$GeneID
  dfc3 <- Cluster3$GeneID
  dfc4 <- Cluster4$GeneID
  dfc5 <- Cluster5$GeneID
  dfc6 <- Cluster6$GeneID
  
  FC_data <- data_results[, c(1,grep(paste(comp_pick, sep = ""),colnames(data_results)))]
  colnames(FC_data) = gsub(paste(comp_pick,"_", sep = ""),"", colnames(FC_data))
  FC_data<-data.frame(FC_data)
  FC_data_sig<-FC_data %>% filter(significant)
  
  #Set up colour values for each cluster
  keyvals.colour <- ifelse(data_results$name %in% dfc1, mypal[1],
                           ifelse(data_results$name %in% dfc2, mypal[2],
                                  ifelse(data_results$name %in% dfc3, mypal[3],
                                         ifelse(data_results$name %in% dfc4, mypal[4],
                                                ifelse(data_results$name %in% dfc5, mypal[5],
                                                       ifelse(data_results$name %in% dfc6, mypal[6],
                                                              "#8F8F8F"))))))
  
  keyvals.colour[is.na(keyvals.colour)] <- "#8F8F8F"
  names(keyvals.colour)[keyvals.colour == mypal[1]] <- 'Cluster1'
  names(keyvals.colour)[keyvals.colour == mypal[2]] <- 'Cluster2'
  names(keyvals.colour)[keyvals.colour == mypal[3]] <- 'Cluster3'
  names(keyvals.colour)[keyvals.colour == mypal[4]] <- 'Cluster4'
  names(keyvals.colour)[keyvals.colour == mypal[5]] <- 'Cluster5'
  names(keyvals.colour)[keyvals.colour == mypal[6]] <- 'Cluster6'
  names(keyvals.colour)[keyvals.colour == '#8F8F8F'] <- 'NA'
  
  if (padj_check =="p-adj sig."){
    enh_vol<-plot(EnhancedVolcano(
      FC_data,
      lab = FC_data$name,
      x = "ratio",
      y = "p.val",
      title = comp_pick,
      subtitle = paste("Fold change cutoff = ", FC ,", p-adj cutoff =  ", pVal, sep = ""),
      pCutoff = pVal,
      FCcutoff = log2(FC),
      pointSize = 3.0,
      labSize = 3.0,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      #maxoverlapsConnectors = Inf,
      colCustom = keyvals.colour,
      colAlpha = (alpha_input/10)*0.5,
      selectLab = FC_data_sig$name #show labels for significant padj proteins
    ))
  }  
  #Volcano plots with pVal significant labels
  if (padj_check =="p-Val sig."){
    enh_vol<-plot(EnhancedVolcano(
      FC_data,
      lab = FC_data$name,
      x = "ratio",
      y = "p.val",
      title = comp_pick,
      subtitle = paste("Fold change cutoff = ", FC ,", p-Val cutoff =  ", pVal, sep = ""),
      pCutoff = pVal,
      FCcutoff = log2(FC),
      pointSize = 3.0,
      labSize = 3.0,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      #maxoverlapsConnectors = Inf,
      colCustom = keyvals.colour,
      colAlpha = (alpha_input/10)*0.5
      #selectLab = FC_data_sig$name #show labels for significant padj proteins
    ))
  } 
  
  #Volcano plot with padj values and padj labels
  
  if (padj_check =="p-adj values"){
    enh_vol<-plot(EnhancedVolcano(
      FC_data,
      lab = FC_data$name,
      x = "ratio",
      y = "p.adj",
      title = comp_pick,
      subtitle = paste("Fold change cutoff = ", FC ,", p-adj cutoff =  ", pVal, sep = ""),
      ylab = bquote(~-Log[10] ~ italic(P) ~ adj), #shows padj
      pCutoff = pVal,
      FCcutoff = log2(FC),
      pointSize = 3.0,
      labSize = 3.0,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      #maxoverlapsConnectors = Inf,
      colCustom = keyvals.colour,
      colAlpha = (alpha_input/10)*0.5
      #selectLab = FC_data_sig$name #show labels for significant padj proteins
    ))
  }  
  return(enh_vol)
}

####Gene Ontology Plots####

####Cellualr Compartment####
get_entrezlist <- function(clu_df){
  Cluster1 =clu_df[clu_df$Cluster == "cluster1", ]
  Cluster2 =clu_df[clu_df$Cluster == "cluster2", ]
  Cluster3 =clu_df[clu_df$Cluster == "cluster3", ]
  Cluster4 =clu_df[clu_df$Cluster == "cluster4", ]
  Cluster5 =clu_df[clu_df$Cluster == "cluster5", ]
  Cluster6 =clu_df[clu_df$Cluster == "cluster6", ]
  
  # Converts each cluster into a list of values
  dfc1 <- Cluster1$GeneID 
  dfc2 <- Cluster2$GeneID
  dfc3 <- Cluster3$GeneID
  dfc4 <- Cluster4$GeneID
  dfc5 <- Cluster5$GeneID
  dfc6 <- Cluster6$GeneID
  
  
  ####To do gene Ontology analysis
  #Convert gene lists (from "SYMBOL", for example) into EntrezIDs.
  # Generate lists of EntrezIDs
  idCluster1 <- c(bitr(
    dfc1,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ))
  idCluster2 <- c(bitr(
    dfc2,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ))
  idCluster3 <- c(bitr(
    dfc3,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ))
  idCluster4 <- c(bitr(
    dfc4,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ))
  idCluster5 <- c(bitr(
    dfc5,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ))
  idCluster6 <- c(bitr(
    dfc6,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ))
  
  # Assemble converted BaitID vectors into a list object
  entrezlist <- list(
    idCluster1$ENTREZID,
    idCluster2$ENTREZID,
    idCluster3$ENTREZID,
    idCluster4$ENTREZID,
    idCluster5$ENTREZID,
    idCluster6$ENTREZID
  )
  names(entrezlist) <-
    c(
      "Cluster 1",
      "Cluster 2",
      "Cluster 3",
      "Cluster 4",
      "Cluster 5",
      "Cluster 6"
    )
return(entrezlist)
}

plot_GO_CC <- function(entrezlist){
  
  #Perform functional enrichment analysis using compareCluster
  #ont = "CC" (cellular compartment),"BP" (Biological process) or "MF" Molecular Function
  #Cellular Compartment
  GOCC <- compareCluster(
    geneClusters = entrezlist,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    ont = "CC",
    readable = T
  )
  #GOCCfile<-paste("Clusters_CC_", label, ".txt",sep="")
  GOCCtitle<-paste("Cellular Component ",sep="")
  #write.table(GOCC, file = GOCCfile, sep = "\t",
  #            row.names = TRUE, col.names = NA)
  
  GOCC_dotplot<- dotplot(GOCC, showCategory = 8) + ggtitle(GOCCtitle)
  
  GOCC_cnetplot<- cnetplot(GOCC, showCategory = 5) + ggtitle(GOCCtitle) 
  
  GO_CC<-list(GOCC_dotplot, GOCC_cnetplot)
  
  return(GO_CC)
  
}

#####Molecular Function####

plot_GO_MF <- function(entrezlist){

  GOMF <- compareCluster(
    geneClusters = entrezlist,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    ont = "MF",
    readable = T
  )
  #GOMFfile<-paste("Clusters_MF_", label, ".txt",sep="")
  GOMFtitle<-paste("Molecular Function ", sep="")
  #write.table(GOMF, file = GOMFfile, sep = "\t",
  #            row.names = TRUE, col.names = NA)
  
  GOMF_dotplot<- dotplot(GOMF, showCategory = 8) + ggtitle(GOMFtitle)
  
  GOMF_cnetplot<- cnetplot(GOMF, showCategory = 5) + ggtitle(GOMFtitle) 
  
  GO_MF<-list(GOMF_dotplot, GOMF_cnetplot)
  
  return(GO_MF)
  
}

####Biological Process####

plot_GO_BP <- function(entrezlist){
  

  GOBP <- compareCluster(
    geneClusters = entrezlist,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    ont = "BP",
    readable = T
  )
  GOBPtitle<-paste("Biological Process ",sep="")
  
  GOBP_dotplot<- dotplot(GOBP, showCategory = 8) + ggtitle(GOBPtitle)
  
  GOBP_cnetplot<- cnetplot(GOBP, showCategory = 5) + ggtitle(GOBPtitle) 
  
  GO_BP<-list(GOBP_dotplot, GOBP_cnetplot)
  
  return(GO_BP)
  
}

####BOX Plot Generator####

#Box Plot Function
plot_box_plot <- function(normdata, intdata, clu_df, cluster_no, palette_choice="default") {
  
  if(palette_choice=="default"){palette_choice<-"Dark2"}
  
  Cluster1 =clu_df[clu_df$Cluster == "cluster1", ]
  Cluster2 =clu_df[clu_df$Cluster == "cluster2", ]
  Cluster3 =clu_df[clu_df$Cluster == "cluster3", ]
  Cluster4 =clu_df[clu_df$Cluster == "cluster4", ]
  Cluster5 =clu_df[clu_df$Cluster == "cluster5", ]
  Cluster6 =clu_df[clu_df$Cluster == "cluster6", ]
  
  # Converts each cluster into a list of values
  dfc1 <- Cluster1$GeneID 
  dfc2 <- Cluster2$GeneID
  dfc3 <- Cluster3$GeneID
  dfc4 <- Cluster4$GeneID
  dfc5 <- Cluster5$GeneID
  dfc6 <- Cluster6$GeneID

  clusterID <- paste("dfc", cluster_no,  sep = "")
  clusterID <- eval(parse(text = clusterID))#converts text to variable
  protein_list <- clusterID
  
  #origdata2 <- normdata %>% remove_rownames %>% column_to_rownames(var="name")
  
  #orig_log2<-log2(origdata2)
  orig_melt<-normdata%>%rownames_to_column()%>%gather(ID, val, -rowname)%>%mutate(Condition=ID)%>%mutate(Rep=str_sub(ID,-1))%>%arrange(rowname)
  orig_melt$Condition = substr(orig_melt$Condition,1,nchar(orig_melt$Condition)-2)
  orig_melt_filt = setDT(orig_melt)[rowname %chin% protein_list]
  
  origtitle<-paste("Original values ", " Cluster ", cluster_no,sep="")
  
  orig<-ggplot(orig_melt_filt, aes(Condition, val))+
    geom_boxplot()+
    geom_jitter(aes(color = factor(Rep)),
                size = 3, position = position_dodge(width=0.3)) +
    labs(
      y = expression(log[2]~"Intensity"),
      col = "Replicates") +
    facet_wrap(~rowname) +
    scale_color_brewer(palette = palette_choice)+
    #theme_linedraw() +
    theme (axis.text.x = element_text(angle=90, hjust = 1, vjust = 1))+
    ggtitle(origtitle)
  #theme(axis.title.x = element_blank())
  
  
  
  ##Imputed box plot
  #intdata2 <- intdata %>% remove_rownames %>% column_to_rownames(var="name")
  #int_log2<-log2(intdata2)
  int_melt<-intdata%>%rownames_to_column()%>%gather(ID, val, -rowname)%>%mutate(Condition=ID)%>%mutate(Rep=str_sub(ID,-1))%>%arrange(rowname)
  int_melt$Condition = substr(int_melt$Condition,1,nchar(int_melt$Condition)-2)
  int_melt_filt = setDT(int_melt)[rowname %chin% protein_list]
  
  Imptitle<-paste("Imputed values ",  " Cluster ", cluster_no,sep="")
  
  Imp<-ggplot(int_melt_filt, aes(Condition, val))+
    geom_boxplot()+
    geom_jitter(aes(color = factor(Rep)),
                size = 3, position = position_dodge(width=0.3)) +
    labs(
      y = expression(log[2]~"Intensity"),
      col = "Replicates") +
    facet_wrap(~rowname) +
    scale_color_brewer(palette = palette_choice)+
    #theme_linedraw() +
    theme (axis.text.x = element_text(angle=90, hjust = 1, vjust = 1))+
    ggtitle(Imptitle)
  #theme(axis.title.x = element_blank())
  
  
  box_plot<-list(orig, Imp)
  
  return(box_plot)
  
}
