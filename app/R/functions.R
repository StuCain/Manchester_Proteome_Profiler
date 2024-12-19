#DEP+Heatmaps_all.R by Stuart Cain
#The University of Manchester
#This script imports Maxquant output file proteinGroups.txt and combines it with a experimental
#design file. The DEP package is used to calculate foldchange values, which can be filtered using Cutoff values.
#The complexHeatmap package is used to render heatmaps, and the clusters are then analysed by clusterprofiler

Version = 4.7


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
  
  nejm <- c(
    "1"="#BC3C29", 
    "2"="#0072B5", 
    "5"="#E18727", 
    "4"="#20854E", 
    "6"="#7876B1", 
    "3"="#6F99AD", 
    "7"="#FFDC91", 
    "8"="#EE4C97"
  )
  
  
  if (palette_choice=="default"){
    mypal<-c6h
  }
  
  if (palette_choice=="nejm"){
    mypal<-nejm
  } 
  
  if (palette_choice!="default"&&palette_choice!="nejm"){
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

####p.Val function####
#add rejections function using pVal instead of padj
add_rejections_pval<-function (diff, alpha = 0.05, lfc = 1) 
{
  if (is.integer(alpha)) 
    alpha <- as.numeric(alpha)
  if (is.integer(lfc)) 
    lfc <- as.numeric(lfc)
  
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"), 
                          is.numeric(alpha), length(alpha) == 1, is.numeric(lfc), 
                          length(lfc) == 1)
  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(diff)), "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  if (length(grep("_p.val|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.val' columns are not present in '", 
         deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns", 
         call. = FALSE)
  }
  cols_p <- grep("_p.val", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  if (length(cols_p) == 1) {
    rowData(diff)$significant <- row_data[, cols_p] <= alpha & 
      abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub("p.val", 
                                                                            "significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df, significant = apply(sign_df, 
                                                  1, function(x) any(x)))
    colnames(sign_df) <- gsub("_p.val", "_significant", colnames(sign_df))
    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    rowData(diff) <- merge(rowData(diff, use.names = FALSE), 
                           sign_df, by = "name")
  }
  return(diff)
}

####Significance Count Functions####
#Using p.adj
signif_count <- function( df, max_cutoff=20 ) {
  sig_count_0.05<-vector(mode="numeric")
  sig_count_0.1<-vector(mode="numeric")
  sig_count_0.5<-vector(mode="numeric")
  
  for (i in 1:20) {
    dep_i<-add_rejections(df, 0.05, log2(i))
    sig_count_0.05[i]<-dep_i %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  }
  for (i in 1:20) {
    dep_i<-add_rejections(df, 0.1, log2(i))
    sig_count_0.1[i]<-dep_i %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  }
  for (i in 1:20) {
    dep_i<-add_rejections(df, 0.5, log2(i))
    sig_count_0.5[i]<-dep_i %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  }
  
  fc_val<-1:20
  sig_plot_data <-data.frame(fc_val,sig_count_0.05, sig_count_0.1, sig_count_0.5)
  sig_plot_graph<- sig_plot_data %>% gather(key = p.adj, value = value, sig_count_0.05:sig_count_0.5)
  sig_plot<-ggplot(sig_plot_graph, aes(x = fc_val, y = value, group = p.adj, colour = p.adj)) +
    labs(title = "Significance Count Plot", 
         x = "Fold Change Cut-off", 
         y = "Number of Significant Proteins" 
    ) + 
    theme_DEP1()+
    geom_line() +
    geom_point()+
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_discrete (name = "p.adj cut-off", labels = c("0.05", "0.1", "0.5"))
  
  return(sig_plot)
}

#Using p.adj
signif_count_pval <- function( df, max_cutoff=20 ) {
  sig_count_0.05<-vector(mode="numeric")
  sig_count_0.1<-vector(mode="numeric")
  sig_count_0.5<-vector(mode="numeric")
  
  for (i in 1:20) {
    dep_i<-add_rejections_pval(df, 0.05, log2(i))
    sig_count_0.05[i]<-dep_i %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  }
  for (i in 1:20) {
    dep_i<-add_rejections_pval(df, 0.1, log2(i))
    sig_count_0.1[i]<-dep_i %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  }
  for (i in 1:20) {
    dep_i<-add_rejections_pval(df, 0.5, log2(i))
    sig_count_0.5[i]<-dep_i %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
  }
  
  fc_val<-1:20
  sig_plot_data <-data.frame(fc_val,sig_count_0.05, sig_count_0.1, sig_count_0.5)
  sig_plot_graph<- sig_plot_data %>% gather(key = p.adj, value = value, sig_count_0.05:sig_count_0.5)
  sig_plot<-ggplot(sig_plot_graph, aes(x = fc_val, y = value, group = p.adj, colour = p.adj)) +
    labs(title = "Significance Count Plot", 
         x = "Fold Change Cut-off", 
         y = "Number of Significant Proteins" 
    ) + 
    theme_DEP1()+
    geom_line() +
    geom_point()+
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_discrete (name = "p.Val cut-off", labels = c("0.05", "0.1", "0.5"))
  
  return(sig_plot)
}


####PCA with txt output####
plot_pca_txt<-function (dep, x = 1, y = 2, indicate = c("condition", "replicate"), 
                        label = FALSE, n = 500, point_size = 4, label_size = 3, plot = TRUE) 
{
  if (is.integer(x)) 
    x <- as.numeric(x)
  if (is.integer(y)) 
    y <- as.numeric(y)
  if (is.integer(n)) 
    n <- as.numeric(n)
  if (is.integer(point_size)) 
    point_size <- as.numeric(point_size)
  if (is.integer(label_size)) 
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.numeric(x), length(x) == 1, is.numeric(y), length(y) == 
                            1, is.numeric(n), length(n) == 1, is.character(indicate), 
                          is.logical(label), is.numeric(point_size), length(point_size) == 
                            1, is.numeric(label_size), length(label_size) == 
                            1, is.logical(plot), length(plot) == 1)
  if (x > ncol(dep) | y > ncol(dep)) {
    stop(paste0("'x' and/or 'y' arguments are not valid\n", 
                "Run plot_pca() with 'x' and 'y' <= ", ncol(dep), 
                "."), call. = FALSE)
  }
  if (n > nrow(dep)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_pca() with 'n' <= ", 
                nrow(dep), "."), call. = FALSE)
  }
  columns <- colnames(colData(dep))
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_pca() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"), 
                  "' column(s) is/are not present in ", deparse(substitute(dep)), 
                  ".\nValid columns are: '", paste(columns, collapse = "', '"), 
                  "'."), call. = FALSE)
    }
  }
  var <- apply(assay(dep), 1, sd)
  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)], 
  ]
  pca <- prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame() %>% rownames_to_column() %>% 
    left_join(., data.frame(colData(dep)), by = c(rowname = "ID"))
  percent <- round(100 * pca$sdev^2/sum(pca$sdev^2), 1)
  for (feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }
  scree_data <- data.frame(
    Component = 1:length(pca$sdev),
    Variance = round(100 * pca$sdev^2/sum(pca$sdev^2), 1)
  )
  
  scree_plot<-ggplot(scree_data, aes(x = Component, y = Variance)) +
    geom_bar(stat = "identity", fill = "#F8766D") +
    geom_line() +
    geom_text(data=scree_data, aes(label=paste0(Variance,"%"),vjust=-1), size =5) +
    geom_point() +
    theme_DEP1() +
    xlab("Principal Components") +
    ylab("Percentage of Variance Explained") +
    ggtitle("Scree Plot")
  
  p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC", y)))) + 
    labs(title = paste0("PCA plot - top ", n, " variable proteins"), 
         x = paste0("PC", x, ": ", percent[x], "%"), 
         y = paste0("PC", y, ": ", percent[y], "%")) + 
    coord_fixed() + theme_DEP1()
  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]]), 
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], 
                            shape = pca_df[[indicate[2]]]), size = point_size) + 
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], 
                            shape = pca_df[[indicate[2]]]), size = point_size) + 
      facet_wrap(~pca_df[[indicate[3]]]) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    pca_out<-list(p, pca_df, scree_plot)
    return(pca_out)
  }
  else {
    df <- pca_df %>% select(rowname, paste0("PC", c(x, y)), 
                            match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    pca_out<-list(df, pca_df, scree_data)
    return(pca_out)
  }
}

####Custom Normalisation Plot####
plot_normalization_vert <-function (se,  ...) 
{
  call <- match.call()
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x, "SummarizedExperiment"), 
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "ID", "condition", "replicate") %in% 
            colnames(colData(x)))) {
      stop("'label', 'ID', 'condition' and/or 'replicate' ", 
           "columns are not present in (one of) the input object(s)", 
           "\nRun make_se() or make_se_parse() to obtain the required columns", 
           call. = FALSE)
    }
  })
  gather_join <- function(se) {
    assay(se) %>% data.frame() %>% gather(ID, val) %>% left_join(., 
                                                                 data.frame(colData(se)), by = "ID")
  }
  df <- map_df(arglist, gather_join, .id = "var") %>% mutate(var = factor(var, 
                                                                          levels = names(arglist)))
  # New facet label names
  df$var <- factor(df$var, labels = c("Filtered Data", "Normalized Data"))
  ggplot(df, aes(x = ID, y = val, fill = condition)) + geom_boxplot(notch = TRUE, 
                                                                    na.rm = TRUE) + coord_flip() + facet_wrap(~var, ncol = 1) + 
    labs(x = "", y = expression(log[2] ~ "Intensity")) + 
    theme_DEP1()
}

plot_normalization_hor <-function (se,  ...) 
{
  call <- match.call()
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x, "SummarizedExperiment"), 
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "ID", "condition", "replicate") %in% 
            colnames(colData(x)))) {
      stop("'label', 'ID', 'condition' and/or 'replicate' ", 
           "columns are not present in (one of) the input object(s)", 
           "\nRun make_se() or make_se_parse() to obtain the required columns", 
           call. = FALSE)
    }
  })
  gather_join <- function(se) {
    assay(se) %>% data.frame() %>% gather(ID, val) %>% left_join(., 
                                                                 data.frame(colData(se)), by = "ID")
  }
  df <- map_df(arglist, gather_join, .id = "var") %>% mutate(var = factor(var, 
                                                                          levels = names(arglist)))
  # New facet label names
  df$var <- factor(df$var, labels = c("Filtered Data", "Normalized Data"))
  ggplot(df, aes(x = ID, y = val, fill = condition)) + geom_boxplot(notch = TRUE, 
                                                                    na.rm = TRUE) + coord_flip() + facet_wrap(~var, ncol = 2) + 
    labs(x = "", y = expression(log[2] ~ "Intensity")) + 
    theme_DEP1()
}


####Complex Heatmap function####
plot_complex_heatmap <- function(data_results, kmeans, ht_title, col_fun, mypal, switch=3, font_size=6, cluster_cols=TRUE, comp_list) {
  mypal2<-adjustcolor(mypal, alpha.f = 0.5)
  mypal<-c(mypal,mypal2)
  
  set.seed(123)
  if (switch==3){
    FC_data <- data_results[, c(1,grep("_ratio",colnames(data_results)))]
    #Change column names
    for ( col in 1:ncol(FC_data)){
      colnames(FC_data)[col] <-  sub("_ratio", "", colnames(FC_data)[col])
    }
    FC_matrix <- FC_data %>% remove_rownames %>% column_to_rownames(var="name")
  }
  
  if (switch==1){imp_matrix<-data_results}
  if (switch==2){norm_matrix<-data_results}
  
  
  ####Set up Heatmaps####
  #Choose which set ofdata to use (Intensities or Log2FC)
  if (switch == 3) {
    #  my_mat <-
    #    as.matrix(dfdata2) #converts table into matrix called my_mat
    my_mat<-FC_matrix
    my_mat <- as.matrix(my_mat)
    if (ncol(my_mat)>1){my_mat<-my_mat[,comp_list]}
    #CS = colMeans(FC_matrix) #creates a mean of each column
    col_fun = col_fun
    kmeans<-kmeans
  } 
  if (switch == 1) {
    my_mat <- as.matrix(imp_matrix)
    #my_mat<- my_mat[ order(rowMeans(my_mat), decreasing = T), ]
    #my_mat<-my_mat[1:100,]
    #my_mat_a <- log2(my_mat)
    CSr = rowMeans(my_mat)
    CSrmax<-max(CSr)
    CSrmin<-min(CSr)
    col_fun = colorRamp2(c(CSrmax+1, CSrmin-1), c("blue", "yellow"))
    col_fun(seq(CSrmax+1, CSrmin-1))
    
    #my_mat <- my_mat - CSr
    CSi = colMeans(my_mat)
    col_fun = col_fun
    kmeans<-kmeans
  }
  
  if (switch == 2) {
    my_mat <- as.matrix(norm_matrix)
    #my_mat <- log2(my_mat)
    CSr = rowMeans(my_mat,na.rm=TRUE)
    CSrmax<-max(CSr)
    CSrmin<-min(CSr)
    col_fun = colorRamp2(c(CSrmax+1, CSrmin-1), c("blue", "yellow"))
    col_fun(seq(CSrmax+1, CSrmin-1))
    #my_mat <- my_mat - CSr
    #my_mat[is.na(my_mat)] = 0
    CSi = colMeans(my_mat,na.rm=TRUE)
    col_fun = col_fun
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
        gp = gpar(fill = mypal[c(1:16)])
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
        gp = gpar(fill = mypal[c(1:16)])
      ),
      Int_mean = CSi,
      col = list(Int_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend=FALSE
    )
    htname = "Log2 Int."
  }
  
  #Right Annotation
  #For RA
  if(switch==3){
    dfdata_sig<-data_results[, c(1,grep("_significant",colnames(data_results)))]
    dfdata_sig2 <- dfdata_sig%>% mutate(sigcount = rowSums(dfdata_sig[-1], na.rm = TRUE))
    #saint_count3 = setDT(dfdata_sig2)[name %chin% LFQ_GL]
    sc3max<- max(dfdata_sig2$sigcount)
    
    ra = rowAnnotation("Sig.
Freq" = anno_barplot(
  dfdata_sig2$sigcount,
  ylim = c(0 , sc3max),
  gp = gpar(fill = mypal[1], fontsize = 8)))
  }else{
    ra = rowAnnotation("Log2Int" = anno_points(
      CSr,
      ylim = c(CSrmin,CSrmax),
      gp = gpar(col = mypal[1], fontsize = 8),
      width = unit(2, "cm")
    ))
  }
  
  
  #Set up left annotation to match kmeans number
  if(switch==3){
    label_txt<-vector()
    for (i in 1:kmeans) {
      label_txt[i]<-paste0("Cluster ", i)
    }
    if (kmeans <= 16) {
      la = rowAnnotation(foo = anno_block(
        gp = gpar(fill = mypal[c(1:16)]),
        labels = label_txt,
        labels_gp = gpar(col = "white", fontsize = font_size+4)
      ))
    } else{
      if (kmeans > 16) {
        la = rowAnnotation(" " = anno_empty(
          border = FALSE,
          zoom = FALSE,
          width = unit(0.01, "cm")
        ))
      }
    }
  }else{la = rowAnnotation(" " = anno_empty(
    border = FALSE,
    zoom = FALSE,
    width = unit(0.01, "cm")
  ))
  }
  
  #ht_title<-paste("FC=", FC, " pVal=", pVal, "<br>", "(",LFQ_GL_Count, " Proteins)", "<br>",kmeans," Clusters ",label , sep="")
  if(switch==3){
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
                 col = col_fun,
                 #column_split = column_split,
                 border = TRUE,
                 row_names_gp = gpar(fontsize = font_size),
                 column_names_gp = gpar(fontsize = font_size+2),
                 column_title = gt_render (ht_title))
  }else{
    
    
    ht = Heatmap(my_mat,
                 name = htname,
                 cluster_rows = FALSE, 
                 cluster_columns = cluster_cols,
                 column_dend_side = "bottom",
                 top_annotation = ha,
                 left_annotation = la,
                 right_annotation = ra,
                 ##bottom_annotation = ba,
                 #row_km = kmeans, row_km_repeats = 2000, # kmeans clustering
                 #cluster_row_slices = FALSE,
                 col = col_fun,
                 na_col = "grey",
                 #column_split = column_split,
                 border = TRUE,
                 row_names_gp = gpar(fontsize = font_size),
                 column_names_gp = gpar(fontsize = font_size+2),
                 column_title = gt_render (ht_title))
  }
  ht = draw(ht) #Important to draw ht, to fix variations in heatmap
  if(switch==3){
    compht<-list(ht, my_mat, dfdata_sig2)}else{compht<-list(ht, my_mat)}
  
  return(compht)
  
}

####2 Round ComplexHeatmap####
plot_complex_heatmap_2xcomp <- function(data_results, kmeans, ht_title, col_fun, mypal, switch=3, font_size=6, cluster_cols=TRUE, comp_list, data_results1, data_results2) {
  mypal2<-adjustcolor(mypal, alpha.f = 0.5)
  mypal<-c(mypal,mypal2)
  
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
        gp = gpar(fill = mypal[c(1:16)])
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
        gp = gpar(fill = mypal[c(1:16)])
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
  
  
  
  
  if(switch==3){
    label_txt<-vector()
    for (i in 1:kmeans) {
      label_txt[i]<-paste0("Cluster ", i)
    }
    if (kmeans <= 16) {
      la = rowAnnotation(foo = anno_block(
        gp = gpar(fill = mypal[c(1:16)]),
        labels = label_txt,
        labels_gp = gpar(col = "white", fontsize = font_size+4)
      ))
    } else{
      if (kmeans > 16) {
        la = rowAnnotation(" " = anno_empty(
          border = FALSE,
          zoom = FALSE,
          width = unit(0.01, "cm")
        ))
      }
    }
  }else{la = rowAnnotation(" " = anno_empty(
    border = FALSE,
    zoom = FALSE,
    width = unit(0.01, "cm")
  ))
  }
  
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

####BFDR Heatmap####
plot_complex_heatmap_BFDR <- function(data_results, kmeans, ht_title, col_fun, mypal, switch=3, font_size=6, cluster_cols=TRUE, comp_list, data_results1, saintlist) {
  mypal2<-adjustcolor(mypal, alpha.f = 0.5)
  mypal<-c(mypal,mypal2)
  
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
        gp = gpar(fill = mypal[c(1:16)])
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
        gp = gpar(fill = mypal[c(1:16)])
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
  
  
  
  
  if(switch==3){
    label_txt<-vector()
    for (i in 1:kmeans) {
      label_txt[i]<-paste0("Cluster ", i)
    }
    if (kmeans <= 16) {
      la = rowAnnotation(foo = anno_block(
        gp = gpar(fill = mypal[c(1:16)]),
        labels = label_txt,
        labels_gp = gpar(col = "white", fontsize = font_size+4)
      ))
    } else{
      if (kmeans > 16) {
        la = rowAnnotation(" " = anno_empty(
          border = FALSE,
          zoom = FALSE,
          width = unit(0.01, "cm")
        ))
      }
    }
  }else{la = rowAnnotation(" " = anno_empty(
    border = FALSE,
    zoom = FALSE,
    width = unit(0.01, "cm")
  ))
  }
  
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

####Cluster Check####
cluster_check <- function(ht, my_mat, kmeans){
  set.seed(123)
  col_fun = colorRamp2(c(-5, 0, 5), c("red", "white", "blue"))
  col_fun(seq(-5, 5))
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
  
  #checkclust<-list(ht_C1, ht_C2, ht_C3, ht_C4, ht_C5, ht_C6, C1_Count, C2_Count, C3_Count, C4_Count, C5_Count, C6_Count, 
  #                  Cluster1FC, Cluster2FC, Cluster3FC, Cluster4FC, Cluster5FC, Cluster6FC, clu_df)
  return(clu_df)
}

####Heatmap Clusters Function####
plot_heatmap_clusters <- function(ht, my_mat, dfdata_sig2, kmeans, FC=2, pVal=0.05, col_fun, mypal, switch=3, font_size=6, data_unique_defn){
  mypal2<-adjustcolor(mypal, alpha.f = 0.5)
  mypal<-c(mypal,mypal2)
  
  set.seed(123)
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
  
  #Generate Lists
  Cluster_list<- vector(mode='list', length=kmeans)
  C_Count_list<- vector(mode='list', length=kmeans)
  my_mat_list <- vector(mode='list', length=kmeans)
  ha_list<- vector(mode='list', length=kmeans)
  la_list<- vector(mode='list', length=kmeans)
  ra_list<- vector(mode='list', length=kmeans)
  ht_list<- vector(mode='list', length=kmeans)
  CS_list<- vector(mode='list', length=kmeans)
  ClusterFC_list<- vector(mode='list', length=kmeans)
  saint_count_list<- vector(mode='list', length=kmeans)
  
  #Create filtered lists
  for (i in 1:kmeans){
    cluster_no<-paste0("cluster", i)
    Cluster_list[[i]] =clu_df[clu_df$Cluster == cluster_no, ]
    C_Count_list[[i]]<-nrow( Cluster_list[[i]])
  }
  
  #Generate separate cluster matrices
  if(comp_length==1){
    my_mat_r <-cbind(my_mat, rep(my_mat[,1],1))
    
    for (i in 1:kmeans){
      my_mat_list[[i]] <- my_mat_r [rownames(my_mat_r) %in% Cluster_list[[i]]$GeneID, ]
      my_mat_list[[i]] <- data.frame(name = row.names(my_mat_list[[i]]), my_mat_list[[i]])
      my_mat_list[[i]] <- my_mat_list[[i]][,-3]
      my_mat_list[[i]] <- my_mat_list[[i]] %>% remove_rownames %>% column_to_rownames(var="name")
    }
  }
  
  
  if(comp_length>1){
    for (i in 1:kmeans){
      my_mat_list[[i]] <- my_mat [rownames(my_mat) %in% Cluster_list[[i]]$GeneID, ]
    }
  }
  
  for (i in 1:kmeans){
    CS_list[[i]]= colMeans(my_mat_list[[i]])
  }
  
  #Generate Cluster txt files
  for (i in 1:kmeans){
    ClusterFC_list[[i]] <- data.frame(name = row.names(my_mat_list[[i]]), my_mat_list[[i]])
    ClusterFC_list[[i]] <- merge(data_unique_defn, ClusterFC_list[[i]], by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)
  }
  
  
  #Calculate Right annotations
  #saint_count3 <- dfdata_sig2
  sc3max<- max(dfdata_sig2$sigcount)
  
  for (i in 1:kmeans){
    saint_count_list[[i]] <-  dfdata_sig2 %>% filter( dfdata_sig2$name %in% Cluster_list[[i]]$GeneID ) 
  }
  
  
  ####Generate Separate Cluster Heatmaps####
  #Generate Annotations
  #Top heatmap annotation
  for (i in 1:kmeans){
    ha_list[[i]] = HeatmapAnnotation(
      log2FC = anno_density(
        my_mat_list[[i]],
        type = "violin",
        height = unit(3, "cm"),
        gp = gpar(fill = mypal[c(1:16)])
      ),
      FC_mean = CS_list[[i]],
      col = list(FC_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend = FALSE
    )
    htname = "log2FC"
  }
  
  #Left heatmap annotation
  for (i in 1:kmeans){
    la_list[[i]] = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[i]),
                                                  labels = paste0("Cluster ", i), 
                                                  labels_gp = gpar(col = "white", fontsize = font_size+4)))
  }
  
  #Right heatmap annotation
  for (i in 1:kmeans){
    ra_list[[i]]=rowAnnotation("Sig.
Freq" = anno_barplot(saint_count_list[[i]]$sigcount, ylim = c( 0 , sc3max ), gp = gpar(fill = mypal[1], fontsize=8)))
  }
  
  for (i in 1:kmeans){
    ht_c = Heatmap(
      my_mat_list[[i]],
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha_list[[i]],
      left_annotation = la_list[[i]],
      #bottom_annotation = ba,
      right_annotation = ra_list[[i]],
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster ",i, " (", C_Count_list[[i]], ")", sep =""))
    )
    ht_list[[i]] = draw(ht_c)
  }
  
  comphtclust<-list(ht_list, C_Count_list, ClusterFC_list, clu_df)
  
  return(comphtclust)
}


plot_heatmap_clusters_2xcomp <- function(ht, my_mat, dfdata_sig2, kmeans, FC=2, pVal=0.05, col_fun, mypal, switch=3, font_size=6, data_unique_defn){
  mypal2<-adjustcolor(mypal, alpha.f = 0.5)
  mypal<-c(mypal,mypal2)
  
  set.seed(123)
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
  
  #Generate Lists
  Cluster_list<- vector(mode='list', length=kmeans)
  C_Count_list<- vector(mode='list', length=kmeans)
  my_mat_list <- vector(mode='list', length=kmeans)
  ha_list<- vector(mode='list', length=kmeans)
  la_list<- vector(mode='list', length=kmeans)
  ra_list<- vector(mode='list', length=kmeans)
  ht_list<- vector(mode='list', length=kmeans)
  CS_list<- vector(mode='list', length=kmeans)
  ClusterFC_list<- vector(mode='list', length=kmeans)
  saint_count_list<- vector(mode='list', length=kmeans)
  
  #Create filtered lists
  for (i in 1:kmeans){
    cluster_no<-paste0("cluster", i)
    Cluster_list[[i]] =clu_df[clu_df$Cluster == cluster_no, ]
    C_Count_list[[i]]<-nrow( Cluster_list[[i]])
  }
  
  #Generate separate cluster matrices
  if(comp_length==1){
    my_mat_r <-cbind(my_mat, rep(my_mat[,1],1))
    
    for (i in 1:kmeans){
      my_mat_list[[i]] <- my_mat_r [rownames(my_mat_r) %in% Cluster_list[[i]]$GeneID, ]
      my_mat_list[[i]] <- data.frame(name = row.names(my_mat_list[[i]]), my_mat_list[[i]])
      my_mat_list[[i]] <- my_mat_list[[i]][,-3]
      my_mat_list[[i]] <- my_mat_list[[i]] %>% remove_rownames %>% column_to_rownames(var="name")
    }
  }
  
  
  if(comp_length>1){
    for (i in 1:kmeans){
      my_mat_list[[i]] <- my_mat [rownames(my_mat) %in% Cluster_list[[i]]$GeneID, ]
    }
  }
  
  for (i in 1:kmeans){
    CS_list[[i]]= colMeans(my_mat_list[[i]])
  }
  
  #Generate Cluster txt files
  for (i in 1:kmeans){
    ClusterFC_list[[i]] <- data.frame(name = row.names(my_mat_list[[i]]), my_mat_list[[i]])
    ClusterFC_list[[i]] <- merge(data_unique_defn, ClusterFC_list[[i]], by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)
  }
  
  
  #Calculate Right annotations
  #saint_count3 <- dfdata_sig2
  sc3max_c1<- max(dfdata_sig2$sigcount_c1)
  sc3max_c2<- max(dfdata_sig2$sigcount_c2)
  
  for (i in 1:kmeans){
    saint_count_list[[i]] <-  dfdata_sig2 %>% filter( dfdata_sig2$name %in% Cluster_list[[i]]$GeneID ) 
  }
  
  #Right heatmap annotation
  for (i in 1:kmeans){
    ra_list[[i]] =  rowAnnotation("Sig.
Freq1" = anno_barplot(saint_count_list[[i]]$sigcount_c1, ylim = c(0 , sc3max_c1), gp = gpar(fill = mypal[1], fontsize = 6))) +
      rowAnnotation("Sig.
Freq2" = anno_barplot(saint_count_list[[i]]$sigcount_c2,ylim = c(0 , sc3max_c2), gp = gpar(fill = mypal[3], fontsize = 6)))
  }
  
  ####Generate Separate Cluster Heatmaps####
  #Generate Annotations
  #Top heatmap annotation
  for (i in 1:kmeans){
    ha_list[[i]] = HeatmapAnnotation(
      log2FC = anno_density(
        my_mat_list[[i]],
        type = "violin",
        height = unit(3, "cm"),
        gp = gpar(fill = mypal[c(1:16)])
      ),
      FC_mean = CS_list[[i]],
      col = list(FC_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend = FALSE
    )
    htname = "log2FC"
  }
  
  #Left heatmap annotation
  for (i in 1:kmeans){
    la_list[[i]] = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[i]),
                                                  labels = paste0("Cluster ", i), 
                                                  labels_gp = gpar(col = "white", fontsize = font_size+4)))
  }
  
  #Generate heatmaps
  for (i in 1:kmeans){
    ht_list[[i]] = Heatmap(
      my_mat_list[[i]],
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha_list[[i]],
      left_annotation = la_list[[i]],
      #bottom_annotation = ba,
      #right_annotation = ra_list[[i]],
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster ",i, " (", C_Count_list[[i]], ")", sep ="")))
    
    ht_list[[i]] = draw(ht_list[[i]]+ra_list[[i]], auto_adjust = FALSE)
  }
  
  comphtclust<-list(ht_list, C_Count_list, ClusterFC_list, clu_df)
  
  return(comphtclust)
}


plot_heatmap_clusters_2xBFDR <- function(ht, my_mat, dfdata_sig2, kmeans, FC=2, pVal=0.05, col_fun, mypal, switch=3, font_size=6, data_unique_defn){
  mypal2<-adjustcolor(mypal, alpha.f = 0.5)
  mypal<-c(mypal,mypal2)
  
  set.seed(123)
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
  
  #Generate Lists
  Cluster_list<- vector(mode='list', length=kmeans)
  C_Count_list<- vector(mode='list', length=kmeans)
  my_mat_list <- vector(mode='list', length=kmeans)
  ha_list<- vector(mode='list', length=kmeans)
  la_list<- vector(mode='list', length=kmeans)
  ra_list<- vector(mode='list', length=kmeans)
  ht_list<- vector(mode='list', length=kmeans)
  CS_list<- vector(mode='list', length=kmeans)
  ClusterFC_list<- vector(mode='list', length=kmeans)
  saint_count_list<- vector(mode='list', length=kmeans)
  
  #Create filtered lists
  for (i in 1:kmeans){
    cluster_no<-paste0("cluster", i)
    Cluster_list[[i]] =clu_df[clu_df$Cluster == cluster_no, ]
    C_Count_list[[i]]<-nrow( Cluster_list[[i]])
  }
  
  #Generate separate cluster matrices
  if(comp_length==1){
    my_mat_r <-cbind(my_mat, rep(my_mat[,1],1))
    
    for (i in 1:kmeans){
      my_mat_list[[i]] <- my_mat_r [rownames(my_mat_r) %in% Cluster_list[[i]]$GeneID, ]
      my_mat_list[[i]] <- data.frame(name = row.names(my_mat_list[[i]]), my_mat_list[[i]])
      my_mat_list[[i]] <- my_mat_list[[i]][,-3]
      my_mat_list[[i]] <- my_mat_list[[i]] %>% remove_rownames %>% column_to_rownames(var="name")
    }
  }
  
  
  if(comp_length>1){
    for (i in 1:kmeans){
      my_mat_list[[i]] <- my_mat [rownames(my_mat) %in% Cluster_list[[i]]$GeneID, ]
    }
  }
  
  for (i in 1:kmeans){
    CS_list[[i]]= colMeans(my_mat_list[[i]])
  }
  
  #Generate Cluster txt files
  for (i in 1:kmeans){
    ClusterFC_list[[i]] <- data.frame(name = row.names(my_mat_list[[i]]), my_mat_list[[i]])
    ClusterFC_list[[i]] <- merge(data_unique_defn, ClusterFC_list[[i]], by.x = "name", by.y = "name", all.x = FALSE, all.y = TRUE)
  }
  
  ####Generate Separate Cluster Heatmaps####
  #Generate Annotations
  
  #Calculate Right annotations
  #saint_count3 <- dfdata_sig2
  sc3max_c1<- max(dfdata_sig2$sigcount_c1)
  sc3max_c2<- max(dfdata_sig2$BFDR_Freq)
  
  for (i in 1:kmeans){
    saint_count_list[[i]] <-  dfdata_sig2 %>% filter( dfdata_sig2$name %in% Cluster_list[[i]]$GeneID ) 
  }
  
  #Right heatmap annotation
  for (i in 1:kmeans){
    ra_list[[i]] =  rowAnnotation("Sig.
Freq1" = anno_barplot(saint_count_list[[i]]$sigcount_c1, ylim = c(0 , sc3max_c1), gp = gpar(fill = mypal[1], fontsize = 6))) +
      rowAnnotation("BFDR
Freq" = anno_barplot(saint_count_list[[i]]$BFDR_Freq,ylim = c(0 , sc3max_c2), gp = gpar(fill = mypal[3], fontsize = 6)))
  }
  
  #Top heatmap annotation
  for (i in 1:kmeans){
    ha_list[[i]] = HeatmapAnnotation(
      log2FC = anno_density(
        my_mat_list[[i]],
        type = "violin",
        height = unit(3, "cm"),
        gp = gpar(fill = mypal[c(1:16)])
      ),
      FC_mean = CS_list[[i]],
      col = list(FC_mean = col_fun),
      gp = gpar(col = "black"),
      show_legend = FALSE
    )
    htname = "log2FC"
  }
  
  #Left heatmap annotation
  for (i in 1:kmeans){
    la_list[[i]] = rowAnnotation(foo = anno_block(gp = gpar(fill = mypal[i]),
                                                  labels = paste0("Cluster ", i), 
                                                  labels_gp = gpar(col = "white", fontsize = font_size+4)))
  }
  
  #Generate heatmaps
  for (i in 1:kmeans){
    ht_list[[i]] = Heatmap(
      my_mat_list[[i]],
      name = htname,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      #column_split = column_split,
      top_annotation = ha_list[[i]],
      left_annotation = la_list[[i]],
      #bottom_annotation = ba,
      #right_annotation = ra_list[[i]],
      #row_km = 3, row_km_repeats = 2000, # kmeans clustering
      col = col_fun,
      #border_gp = gpar(col = "black", lty = 2),
      border = TRUE,
      #rect_gp = gpar(col = "white", lwd = 1),
      row_names_gp = gpar(fontsize = font_size),
      column_names_gp = gpar(fontsize = font_size+2),
      column_title = (paste("Cluster ",i, " (", C_Count_list[[i]], ")", sep ="")))
    
    ht_list[[i]] = draw(ht_list[[i]]+ra_list[[i]], auto_adjust = FALSE)
  }
  
  comphtclust<-list(ht_list, C_Count_list, ClusterFC_list, clu_df)
  
  return(comphtclust)
  
}

#####Enhanced Volcano Function####
#choices ("p-adj sig.","p-Val sig.", "p-adj values")
plot_enhanced_volcano <- function(data_results, clu_df, comp_pick,  FC=2, pVal=0.05, padj_check="p-adj sig.", mypal, alpha_input=10, connect=TRUE, selected=FALSE, selected_names="")
{
  mypal2<-adjustcolor(mypal, alpha.f = 0.5)
  mypal<-c(mypal,mypal2)
  
  Cluster1 =clu_df[clu_df$Cluster == "cluster1", ]
  Cluster2 =clu_df[clu_df$Cluster == "cluster2", ]
  Cluster3 =clu_df[clu_df$Cluster == "cluster3", ]
  Cluster4 =clu_df[clu_df$Cluster == "cluster4", ]
  Cluster5 =clu_df[clu_df$Cluster == "cluster5", ]
  Cluster6 =clu_df[clu_df$Cluster == "cluster6", ]
  Cluster7 =clu_df[clu_df$Cluster == "cluster7", ]
  Cluster8 =clu_df[clu_df$Cluster == "cluster8", ]
  Cluster9 =clu_df[clu_df$Cluster == "cluster9", ]
  Cluster10 =clu_df[clu_df$Cluster == "cluster10", ]
  
  
  
  # Converts each cluster into a list of values
  dfc1 <- Cluster1$GeneID 
  dfc2 <- Cluster2$GeneID
  dfc3 <- Cluster3$GeneID
  dfc4 <- Cluster4$GeneID
  dfc5 <- Cluster5$GeneID
  dfc6 <- Cluster6$GeneID
  dfc7 <- Cluster7$GeneID
  dfc8 <- Cluster8$GeneID
  dfc9 <- Cluster9$GeneID
  dfc10 <- Cluster10$GeneID
  
  
  FC_data <- data_results[, c(1,grep(paste(comp_pick, sep = ""),colnames(data_results)))]
  colnames(FC_data) = gsub(paste(comp_pick,"_", sep = ""),"", colnames(FC_data))
  FC_data<-data.frame(FC_data)
  FC_data_sig<-FC_data %>% filter(significant)
  
  if (selected == TRUE){vol_names <- selected_names}else{vol_names<-FC_data_sig$name}
  
  #obtain minimum p_val number
  p_val_min<-min(FC_data$p.val)
  p_Val_min<- -log10(p_val_min)
  p_adj_min<-min(FC_data$p.adj)
  p_adj_min<- -log10(p_adj_min)
  
  
  #Set up colour values for each cluster
  keyvals.colour <- ifelse(data_results$name %in% dfc1, mypal[1],
                           ifelse(data_results$name %in% dfc2, mypal[2],
                                  ifelse(data_results$name %in% dfc3, mypal[3],
                                         ifelse(data_results$name %in% dfc4, mypal[4],
                                                ifelse(data_results$name %in% dfc5, mypal[5],
                                                       ifelse(data_results$name %in% dfc6, mypal[6],
                                                              ifelse(data_results$name %in% dfc7, mypal[7],
                                                                     ifelse(data_results$name %in% dfc8, mypal[8],
                                                                            ifelse(data_results$name %in% dfc9, mypal[9],
                                                                                   ifelse(data_results$name %in% dfc10, mypal[10],
                                                                                          "#8F8F8F"))))))))))
  
  keyvals.colour[is.na(keyvals.colour)] <- "#8F8F8F"
  names(keyvals.colour)[keyvals.colour == mypal[1]] <- 'Cluster1'
  names(keyvals.colour)[keyvals.colour == mypal[2]] <- 'Cluster2'
  names(keyvals.colour)[keyvals.colour == mypal[3]] <- 'Cluster3'
  names(keyvals.colour)[keyvals.colour == mypal[4]] <- 'Cluster4'
  names(keyvals.colour)[keyvals.colour == mypal[5]] <- 'Cluster5'
  names(keyvals.colour)[keyvals.colour == mypal[6]] <- 'Cluster6'
  names(keyvals.colour)[keyvals.colour == mypal[7]] <- 'Cluster7'
  names(keyvals.colour)[keyvals.colour == mypal[8]] <- 'Cluster8'
  names(keyvals.colour)[keyvals.colour == mypal[9]] <- 'Cluster9'
  names(keyvals.colour)[keyvals.colour == mypal[10]] <- 'Cluster10'
  names(keyvals.colour)[keyvals.colour == '#8F8F8F'] <- 'ns'
  
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
      drawConnectors = connect,
      widthConnectors = 0.5,
      ylim = c(0, p_Val_min),
      #maxoverlapsConnectors = Inf,
      colCustom = keyvals.colour,
      colAlpha = (alpha_input/10)*0.5,
      selectLab = vol_names #show labels for significant padj proteins
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
      drawConnectors = connect,
      widthConnectors = 0.5,
      ylim = c(0, p_Val_min),
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
      drawConnectors = connect,
      widthConnectors = 0.5,
      ylim = c(0, p_adj_min),
      #maxoverlapsConnectors = Inf,
      colCustom = keyvals.colour,
      colAlpha = (alpha_input/10)*0.5,
      selectLab = vol_names #show labels for significant padj proteins
    ))
  }  
  return(enh_vol)
}

####Gene Ontology Plots####

####Entrelist####
get_entrezlist <- function(clu_df, species="Human"){
  if (species=="Mouse"){species_db<-"org.Mm.eg.db"}
  if (species=="Human"){species_db<-"org.Hs.eg.db"}
  if (species=="Arabidopsis"){species_db<-"org.At.tair.db"}
  
  clu_df_no<-length(unique(clu_df$Cluster))
  clu_df_val<-unique(clu_df$Cluster)
  Cluster_list<- vector(mode='list', length=clu_df_no)
  entrezlist <- vector(mode='list', length=clu_df_no)
  entrezlist_names <- vector()
  
  for (i in 1:clu_df_no){
    Cluster_list[[i]] =clu_df[clu_df$Cluster == clu_df_val[i], ]
  } 
  
  
  ####To do gene Ontology analysis
  #Convert gene lists (from "SYMBOL", for example) into EntrezIDs.
  # Generate lists of EntrezIDs
  
  for (i in 1:clu_df_no){
    entrezlist[[i]] <- c(bitr(
      Cluster_list[[i]]$GeneID,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = species_db
    ))
  }
  
  for (i in 1:clu_df_no){
    entrezlist[[i]] <-entrezlist[[i]]$ENTREZID
  } 
  
  for (i in 1:clu_df_no){
    entrezlist_names[i]<-paste0("Cluster ", i)
  } 
  
  names(entrezlist)<-c(entrezlist_names)
  return(entrezlist)
}

####Cellular Compartment####

plot_GO_CC <- function(entrezlist, species="Human"){
  if (species=="Mouse"){species_db<-"org.Mm.eg.db"}
  if (species=="Human"){species_db<-"org.Hs.eg.db"}
  if (species=="Arabidopsis"){species_db<-"org.At.tair.db"}
  
  GOCC <- compareCluster(
    geneClusters = entrezlist,
    fun = "enrichGO",
    OrgDb = species_db,
    pvalueCutoff = 0.05,
    ont = "CC",
    readable = T
  )

  return(GOCC)
  
}

#####Molecular Function####

plot_GO_MF <- function(entrezlist, species="Human"){
  if (species=="Mouse"){species_db<-"org.Mm.eg.db"}
  if (species=="Human"){species_db<-"org.Hs.eg.db"}
  if (species=="Arabidopsis"){species_db<-"org.At.tair.db"}
  
  GOMF <- compareCluster(
    geneClusters = entrezlist,
    fun = "enrichGO",
    OrgDb = species_db,
    pvalueCutoff = 0.05,
    ont = "MF",
    readable = T
  )
  return(GOMF)
  
}

####Biological Process####

plot_GO_BP <- function(entrezlist, species="Human"){
  if (species=="Mouse"){species_db<-"org.Mm.eg.db"}
  if (species=="Human"){species_db<-"org.Hs.eg.db"}
  if (species=="Arabidopsis"){species_db<-"org.At.tair.db"}
  
  GOBP <- compareCluster(
    geneClusters = entrezlist,
    fun = "enrichGO",
    OrgDb = species_db,
    pvalueCutoff = 0.05,
    ont = "BP",
    readable = T
  )

  return(GOBP)
  
}

####BOX Plot Generator####

#Box Plot Function
plot_box_plot <- function(normdata, intdata, clu_df, cluster_no, palette_choice="default") {
  clu_df_no<-length(unique(clu_df$Cluster))
  clu_df_val<-unique(clu_df$Cluster)
  Cluster_list<- vector(mode='list', length=clu_df_no)
  
  if(palette_choice=="default"){palette_choice<-"Dark2"}
  
  if(cluster_no>0){
  for (i in 1:clu_df_no){
    Cluster_list[[i]] =clu_df[clu_df$Cluster == clu_df_val[i], ]
  } }
  
  if(cluster_no>0){
  protein_list <- Cluster_list[[cluster_no]]$GeneID
  }else{protein_list <- clu_df$GeneID}
  
  #origdata2 <- normdata %>% remove_rownames %>% column_to_rownames(var="name")
  
  #orig_log2<-log2(origdata2)
  orig_melt<-normdata%>%rownames_to_column()%>%gather(ID, val, -rowname)%>%mutate(Condition=ID)%>%mutate(Rep=str_sub(ID,-1))%>%arrange(rowname)
  orig_melt$Condition = substr(orig_melt$Condition,1,nchar(orig_melt$Condition)-2)
  orig_melt_filt = setDT(orig_melt)[rowname %chin% protein_list]
  
  if(cluster_no > 0){
    origtitle<-paste("Original values ", " Cluster ", cluster_no,sep="")
  }else{origtitle<-paste("Original values", sep="")}
  
  orig<-ggplot(orig_melt_filt, aes(Condition, val))+
    geom_boxplot()+
    theme_DEP1()+
    geom_jitter(aes(color = factor(Rep)),
                size = 3, position = position_dodge(width=0.3)) +
    labs(
      y = expression(log[2]~"Intensity"),
      col = "Replicates") +
    facet_wrap(~rowname) +
    #scale_color_brewer(palette = palette_choice)+
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
  
  if(cluster_no > 0){
    Imptitle<-paste("Imputed values ",  " Cluster ", cluster_no,sep="")
  }else{Imptitle<-paste("Imputed values ", sep="")}
  
  
  Imp<-ggplot(int_melt_filt, aes(Condition, val))+
    geom_boxplot()+
    theme_DEP1()+
    geom_jitter(aes(color = factor(Rep)),
                size = 3, position = position_dodge(width=0.3)) +
    labs(
      y = expression(log[2]~"Intensity"),
      col = "Replicates") +
    facet_wrap(~rowname) +
    #scale_color_brewer(palette = palette_choice)+
    #theme_linedraw() +
    theme (axis.text.x = element_text(angle=90, hjust = 1, vjust = 1))+
    ggtitle(Imptitle)
  #theme(axis.title.x = element_blank())
  
  
  box_plot<-list(orig, Imp)
  
  return(box_plot)
  
}
