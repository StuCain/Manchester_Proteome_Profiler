Version = 1.0

make_unique <- function (proteins, names, ids, delim = ";") 
{
  assertthat::assert_that(is.data.frame(proteins), is.character(names), 
                          length(names) == 1, is.character(ids), length(ids) == 
                            1, is.character(delim), length(delim) == 1)
  col_names <- colnames(proteins)
  if (!names %in% col_names) {
    stop("'", names, "' is not a column in '", deparse(substitute(proteins)), 
         "'", call. = FALSE)
  }
  if (!ids %in% col_names) {
    stop("'", ids, "' is not a column in '", deparse(substitute(proteins)), 
         "'", call. = FALSE)
  }
  if (tibble::is_tibble(proteins)) 
    proteins <- as.data.frame(proteins)
  double_NAs <- apply(proteins[, c(names, ids)], 1, function(x) all(is.na(x)))
  if (any(double_NAs)) {
    stop("NAs in both the 'names' and 'ids' columns")
  }
  proteins_unique <- proteins %>% mutate(name = gsub(paste0(delim, 
                                                            ".*"), "", get(names)), ID = gsub(paste0(delim, ".*"), 
                                                                                              "", get(ids)), name = make.unique(ifelse(name == "" | 
                                                                                                                                         is.na(name), ID, name)))
  return(proteins_unique)
}

delete_prefix <- function (words)
{
  prefix <- get_prefix(words)
  gsub(paste0("^", prefix), "", words)
}

get_prefix <- function (words)
{
  assertthat::assert_that(is.character(words))
  if (length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  if (any(is.na(words))) {
    stop("'words' contains NAs")
  }
  minlen <- min(nchar(words))
  truncated <- substr(words, 1, minlen)
  if (minlen < 1) {
    stop("At least one of the elements is too short")
  }
  mat <- data.frame(strsplit(truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) ==
                       1)
  prefix <- as.logical(cumprod(identical))
  paste(mat[prefix, 1], collapse = "")
}

get_suffix<- function (words) 
{
  assertthat::assert_that(is.character(words))
  if (length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  if (any(is.na(words))) {
    stop("'words' contains NAs")
  }
  minlen <- min(nchar(words))
  truncated <- substr(words, nchar(words) - minlen + 1, nchar(words))
  if (minlen < 1) {
    stop("At least one of the elements is too short")
  }
  rev_string <- function(str) {
    paste(rev(strsplit(str, "")[[1]]), collapse = "")
  }
  rev_truncated <- vapply(truncated, rev_string, character(1))
  mat <- data.frame(strsplit(rev_truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 
                       1)
  prefix <- as.logical(cumprod(identical))
  rev_string(paste(mat[prefix, 1], collapse = ""))
}


delete_suffix <- function (words)
{
  suffix <- get_suffix(words)
  gsub(paste0("^", suffix), "", words)
}


make_se <- function (proteins_unique, columns, expdesign) 
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.integer(columns), 
                          is.data.frame(expdesign))
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns", 
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns", 
         "are not present in the experimental design", call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument", 
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique)) 
    proteins_unique <- as.data.frame(proteins_unique)
  if (tibble::is_tibble(expdesign)) 
    expdesign <- as.data.frame(expdesign)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>% 
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(delete_prefix(expdesign$label)), 
                   make.names(delete_prefix(colnames(raw))))
  if (any(is.na(matched))) {
    stop("None of the labels in the experimental design match ", 
         "with column names in 'proteins_unique'", "\nRun make_se() with the correct labels in the experimental design", 
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  se <- SummarizedExperiment(assays = as.matrix(raw), colData = expdesign, 
                             rowData = row_data)
  return(se)
}

filter_missval <-function (se, thr = 0) 
{
  if (is.integer(thr)) 
    thr <- as.numeric(thr)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), 
                          is.numeric(thr), length(thr) == 1)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '", 
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", 
         call. = FALSE)
  }
  max_repl <- max(colData(se)$replicate)
  if (thr < 0 | thr > max_repl) {
    stop("invalid filter threshold applied", "\nRun filter_missval() with a threshold ranging from 0 to ", 
         max_repl)
  }
  bin_data <- assay(se)
  idx <- is.na(assay(se))
  bin_data[!idx] <- 1
  bin_data[idx] <- 0
  keep <- bin_data %>% data.frame() %>% rownames_to_column() %>% 
    gather(ID, value, -rowname) %>% left_join(., data.frame(colData(se)), 
                                              by = "ID") %>% group_by(rowname, condition) %>% summarize(miss_val = n() - 
                                                                                                          sum(value)) %>% filter(miss_val <= thr) %>% spread(condition, 
                                                                                                                                                             miss_val)
  se_fltrd <- se[keep$rowname, ]
  return(se_fltrd)
}

normalize_vsn <- function (se) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_vsn <- se
  vsn.fit <- vsn::vsnMatrix(2^assay(se_vsn))
  assay(se_vsn) <- vsn::predict(vsn.fit, 2^assay(se_vsn))
  return(se_vsn)
}

impute <-function (se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet", 
                      "MinProb", "man", "min", "zero", "mixed", "nbavg"), ...) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), 
                          is.character(fun))
  fun <- match.arg(fun)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  if (!any(is.na(assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)), 
            "'. ", "Returning the unchanged object.", call. = FALSE)
    return(se)
  }
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
  if (fun == "man") {
    se <- manual_impute(se, ...)
  }
  else {
    MSnSet_data <- as(se, "MSnSet")
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, 
                                      ...)
    assay(se) <- MSnbase::exprs(MSnSet_imputed)
  }
  return(se)
}

manual_impute <- function(se, scale = 0.3, shift = 1.8, seed=NULL) {
  if (is.null(seed)) seed <- 123
  if (is.integer(scale)) scale <- is.numeric(scale)
  if (is.integer(shift)) shift <- is.numeric(shift)
  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(se, "SummarizedExperiment"),
    is.numeric(scale),
    length(scale) == 1,
    is.numeric(shift),
    length(shift) == 1
  )
  
  se_assay <- assay(se)
  
  # Show error if there are no missing values
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), "'",
         call. = FALSE
    )
  }
  
  # Get descriptive parameters of the current sample distributions
  stat <- se_assay %>%
    data.frame(check.names = F) %>%
    rownames_to_column() %>%
    gather(samples, value, -rowname) %>%
    filter(!is.na(value)) %>%
    group_by(samples) %>%
    summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value),
      n = n(),
      infin = nrow(se_assay) - n
    )
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
  set.seed(seed)
  for (a in seq_len(nrow(stat))) {
    assay(se)[is.na(assay(se)[, stat$samples[a]]), stat$samples[a]] <-
      rnorm(stat$infin[a],
            mean = stat$median[a] - shift * stat$sd[a],
            sd = stat$sd[a] * scale
      )
  }
  return(se)
}


test_diff <-function (se, type = c("control", "all", "manual"), control = NULL, 
          test = NULL, design_formula = formula(~0 + condition)) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), 
                          is.character(type), class(design_formula) == "formula")
  type <- match.arg(type)
  col_data <- colData(se)
  raw <- assay(se)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '", 
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", 
         call. = FALSE)
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), 
            "'")
  }
  if (!is.null(control)) {
    assertthat::assert_that(is.character(control), length(control) == 
                              1)
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '", 
           paste0(unique(col_data$condition), collapse = "', '"), 
           "'", call. = FALSE)
    }
  }
  variables <- terms.formula(design_formula) %>% attr(., "variables") %>% 
    as.character() %>% .[-1]
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, 
                    collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control, 
                                                    "- ", sep = " "), "", .) %>% paste(" - ", control, 
                                                                                       sep = "")
      }
    }
  }
  if (type == "control") {
    if (is.null(control)) 
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control], 
                    control, sep = " - ")
  }
  if (type == "manual") {
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))
    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'", 
           ".\nValid contrasts should contain combinations of: '", 
           paste0(conditions, collapse = "', '"), "', for example '", 
           paste0(conditions[1], "_vs_", conditions[2]), 
           "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", test)
  }
  message("Tested contrasts: ", paste(gsub(" - ", "_vs_", cntrst), 
                                      collapse = ", "))
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  if (any(is.na(raw))) {
    for (i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- makeContrasts(contrasts = i, levels = design[, 
                                                                      covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates], 
                                           single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 
                                                                         1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 
                                                                             1]
    }
  }
  eB_fit <- eBayes(contrast_fit)
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf, 
                    confint = TRUE)
    res <- res[!is.na(res$t), ]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }
  limma_res <- map_df(cntrst, retrieve_fun)
  table <- limma_res %>% select(rowname, logFC, CI.L, CI.R, 
                                P.Value, qval, comparison) %>% mutate(comparison = gsub(" - ", 
                                                                                        "_vs_", comparison)) %>% gather(variable, value, -c(rowname, 
                                                                                                                                            comparison)) %>% mutate(variable = recode(variable, logFC = "diff", 
                                                                                                                                                                                      P.Value = "p.val", qval = "p.adj")) %>% unite(temp, comparison, 
                                                                                                                                                                                                                                    variable) %>% spread(temp, value)
  rowData(se) <- merge(rowData(se, use.names = FALSE), table, 
                       by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)
  return(se)
}

add_rejections<-function (diff, alpha = 0.05, lfc = 1) 
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
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '", 
         deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns", 
         call. = FALSE)
  }
  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  if (length(cols_p) == 1) {
    rowData(diff)$significant <- row_data[, cols_p] <= alpha & 
      abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub("p.adj", 
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
    colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    rowData(diff) <- merge(rowData(diff, use.names = FALSE), 
                           sign_df, by = "name")
  }
  return(diff)
}

get_results<-function (dep) 
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(dep)), "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '", 
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns", 
         call. = FALSE)
  }
  row_data$mean <- rowMeans(assay(dep), na.rm = TRUE)
  centered <- assay(dep) - row_data$mean
  centered <- data.frame(centered) %>% rownames_to_column() %>% 
    gather(ID, val, -rowname) %>% left_join(., data.frame(colData(dep)), 
                                            by = "ID")
  centered <- group_by(centered, rowname, condition) %>% summarize(val = mean(val, 
                                                                              na.rm = TRUE)) %>% mutate(val = signif(val, digits = 3)) %>% 
    spread(condition, val)
  colnames(centered)[2:ncol(centered)] <- paste(colnames(centered)[2:ncol(centered)], 
                                                "_centered", sep = "")
  ratio <- as.data.frame(row_data) %>% column_to_rownames("name") %>% 
    select(ends_with("diff")) %>% signif(., digits = 3) %>% 
    rownames_to_column()
  colnames(ratio)[2:ncol(ratio)] <- gsub("_diff", "_ratio", 
                                         colnames(ratio)[2:ncol(ratio)])
  df <- left_join(ratio, centered, by = "rowname")
  pval <- as.data.frame(row_data) %>% column_to_rownames("name") %>% 
    select(ends_with("p.val"), ends_with("p.adj"), ends_with("significant")) %>% 
    rownames_to_column()
  pval[, grep("p.adj", colnames(pval))] <- pval[, grep("p.adj", 
                                                       colnames(pval))] %>% signif(digits = 3)
  ids <- as.data.frame(row_data) %>% select(name, ID)
  table <- left_join(ids, pval, by = c(name = "rowname"))
  table <- left_join(table, df, by = c(name = "rowname")) %>% 
    arrange(desc(significant))
  return(table)
}

get_df_long<-function (se) 
{
  assert_that(inherits(se, "SummarizedExperiment"))
  if (!"name" %in% colnames(rowData(se, use.names = FALSE))) {
    stop("'name' column is not present in '", deparse(substitute(se)), 
         "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  col_data <- colData(se) %>% data.frame() %>% rownames_to_column() %>% 
    select(-ID)
  row_data <- rowData(se, use.names = FALSE) %>% data.frame()
  assay_data <- assay(se) %>% data.frame() %>% rownames_to_column()
  colnames(assay_data)[1] <- "name"
  long_assay <- assay_data %>% gather("rowname", "intensity", 
                                      -name)
  long <- long_assay %>% full_join(col_data, ., by = "rowname") %>% 
    select(-rowname) %>% full_join(., row_data, by = "name")
  return(long)
}

get_suffix<- function (words) 
{
  assertthat::assert_that(is.character(words))
  if (length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  if (any(is.na(words))) {
    stop("'words' contains NAs")
  }
  minlen <- min(nchar(words))
  truncated <- substr(words, nchar(words) - minlen + 1, nchar(words))
  if (minlen < 1) {
    stop("At least one of the elements is too short")
  }
  rev_string <- function(str) {
    paste(rev(strsplit(str, "")[[1]]), collapse = "")
  }
  rev_truncated <- vapply(truncated, rev_string, character(1))
  mat <- data.frame(strsplit(rev_truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 
                       1)
  prefix <- as.logical(cumprod(identical))
  rev_string(paste(mat[prefix, 1], collapse = ""))
}


theme_DEP1<-function () 
{
  basesize <- 12
  theme <- ggplot2::theme_bw(base_size = basesize)
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 2
  theme$plot.title$hjust <- 0.5
  theme$axis.title.x$size <- basesize + 2
  theme$axis.title.y$size <- basesize + 2
  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"
  theme$legend.title$size <- basesize + 2
  theme$legend.text$size <- basesize
  theme$strip.text$face <- "bold"
  theme$strip.text$size <- basesize + 2
  theme$strip.text$colour <- "black"
  return(theme)
}

theme_DEP2<-function () 
{
  theme <- theme_DEP1()
  theme$axis.text.x$angle <- 90
  theme$axis.text.x$hjust <- 1
  theme$axis.text.x$vjust <- 0.5
  return(theme)
}

plot_normalization <-function (se, ...) 
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
  ggplot(df, aes(x = ID, y = val, fill = condition)) + geom_boxplot(notch = TRUE, 
                                                                    na.rm = TRUE) + coord_flip() + facet_wrap(~var, ncol = 1) + 
    labs(x = "", y = expression(log[2] ~ "Intensity")) + 
    theme_DEP1()
}

plot_detect<-function (se) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), 
         "'", call. = FALSE)
  }
  df <- se_assay %>% data.frame() %>% rownames_to_column() %>% 
    gather(ID, val, -rowname)
  stat <- df %>% group_by(rowname) %>% summarize(mean = mean(val, 
                                                             na.rm = TRUE), missval = any(is.na(val)))
  cumsum <- stat %>% group_by(missval) %>% arrange(mean) %>% 
    mutate(num = 1, cs = cumsum(num), cs_frac = cs/n())
  p1 <- ggplot(stat, aes(mean, col = missval)) + geom_density(na.rm = TRUE) + 
    labs(x = expression(log[2] ~ "Intensity"), y = "Density") + 
    guides(col = guide_legend(title = "Missing values")) + 
    theme_DEP1()
  p2 <- ggplot(cumsum, aes(mean, cs_frac, col = missval)) + 
    geom_line() + labs(x = expression(log[2] ~ "Intensity"), 
                       y = "Cumulative fraction") + guides(col = guide_legend(title = "Missing values")) + 
    theme_DEP1()
  gridExtra::grid.arrange(p1, p2, ncol = 1)
}

plot_imputation <-function (se, ...) 
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
  ggplot(df, aes(val, col = condition)) + geom_density(na.rm = TRUE) + 
    facet_wrap(~var, ncol = 1) + labs(x = expression(log[2] ~ 
                                                       "Intensity"), y = "Density") + theme_DEP1()
}

plot_numbers <- function (se, plot = TRUE) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), 
                          is.logical(plot), length(plot) == 1)
  df <- assay(se) %>% data.frame() %>% rownames_to_column() %>% 
    gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin), 
                                                      0, 1))
  stat <- df %>% group_by(ID) %>% summarize(n = n(), sum = sum(bin)) %>% 
    left_join(., data.frame(colData(se)), by = "ID")
  p <- ggplot(stat, aes(x = ID, y = sum, fill = condition)) + 
    geom_col() + geom_hline(yintercept = unique(stat$n)) + 
    labs(title = "Proteins per sample", x = "", y = "Number of proteins") + 
    theme_DEP2()
  if (plot) {
    return(p)
  }
  else {
    df <- as.data.frame(stat)
    colnames(df)[seq_len(3)] <- c("sample", "total_proteins", 
                                  "proteins_in_sample")
    return(df)
  }
}

plot_coverage <- function (se, plot = TRUE) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), 
                          is.logical(plot), length(plot) == 1)
  df <- assay(se) %>% data.frame() %>% rownames_to_column() %>% 
    gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin), 
                                                      0, 1))
  stat <- df %>% group_by(rowname) %>% summarize(sum = sum(bin))
  table <- table(stat$sum) %>% data.frame()
  p <- ggplot(table, aes(x = "all", y = Freq, fill = Var1)) + 
    geom_col(col = "white") + scale_fill_grey(start = 0.8, 
                                              end = 0.2) + labs(title = "Protein coverage", x = "", 
                                                                y = "Number of proteins", fill = "Samples") + theme_DEP1()
  if (plot) {
    return(p)
  }
  else {
    df <- as.data.frame(table)
    colnames(df) <- c("samples", "proteins")
    return(df)
  }
}

plot_cor <- function (dep, significant = TRUE, lower = -1, upper = 1, pal = "PRGn", 
          pal_rev = FALSE, indicate = NULL, font_size = 12, plot = TRUE, 
          ...) 
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.logical(significant), length(significant) == 1, is.numeric(lower), 
                          length(lower) == 1, is.numeric(upper), length(upper) == 
                            1, is.character(pal), length(pal) == 1, is.logical(pal_rev), 
                          length(pal_rev) == 1, is.numeric(font_size), length(font_size) == 
                            1, is.logical(plot), length(plot) == 1)
  if (!(lower >= -1 & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid\n         Run plot_pca() with 'lower' and 'upper' between -1 and 1", 
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>% 
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel", " (qualitative panels also not allowed)\n", 
         "Run plot_pca() with one of the following 'pal' options: ", 
         paste(pals$rowname, collapse = "', '"), "'", call. = FALSE)
  }
  if (any(is.na(assay(dep)))) {
    stop("Missing values in '", deparse(substitute(dep)), 
         "'. Use plot_dist() instead")
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- colData(dep) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"), 
           "' column(s) is/are not present in ", deparse(substitute(dep)), 
           ".\nValid columns are: '", paste(columns, collapse = "', '"), 
           "'.", call. = FALSE)
    }
    anno <- colData(dep) %>% data.frame() %>% select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1) 
        cols <- c("black")
      if (length(var) == 2) 
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2) 
        cols <- RColorBrewer::brewer.pal(length(var), 
                                         "Pastel1")
      if (length(var) >= 7) 
        cols <- RColorBrewer::brewer.pal(length(var), 
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }
  else {
    ha1 <- NULL
  }
  if (significant) {
    if (!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
      stop("'significant' column is not present in '", 
           deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required column", 
           call. = FALSE)
    }
    dep <- dep[rowData(dep, use.names = FALSE)$significant, 
    ]
  }
  cor_mat <- cor(assay(dep))
  ht1 = Heatmap(cor_mat, col = circlize::colorRamp2(seq(lower, 
                                                        upper, ((upper - lower)/7)), if (pal_rev) {
                                                          rev(RColorBrewer::brewer.pal(8, pal))
                                                        }
                                                    else {
                                                      RColorBrewer::brewer.pal(8, pal)
                                                    }), heatmap_legend_param = list(color_bar = "continuous", 
                                                                                    legend_direction = "horizontal", legend_width = unit(5, 
                                                                                                                                         "cm"), title_position = "topcenter"), name = "Pearson correlation", 
                column_names_gp = gpar(fontsize = font_size), row_names_gp = gpar(fontsize = font_size), 
                top_annotation = ha1, ...)
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(cor_mat)
    return(df)
  }
}


plot_cond <-function (dep, plot = TRUE) 
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.logical(plot), length(plot) == 1)
  row_data <- rowData(dep, use.names = FALSE)
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '", 
         deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required columns", 
         call. = FALSE)
  }
  if (!"significant" %in% colnames(row_data)) {
    stop("'significant' column is not present in '", deparse(substitute(dep)), 
         "'\nRun add_rejections() to obtain the required column", 
         call. = FALSE)
  }
  significant <- dep[row_data$significant, ]
  row_data <- rowData(significant, use.names = FALSE) %>% data.frame()
  cols <- grep("_significant", colnames(row_data))
  colnames(row_data)[cols] <- gsub("_significant", "", colnames(row_data)[cols])
  row_data_renamed <- row_data
  colnames(row_data_renamed)[cols] <- LETTERS[seq(to = length(cols))]
  legend <- data.frame(symbol = colnames(row_data_renamed)[cols], 
                       names = colnames(row_data)[cols])
  df <- select(row_data_renamed, name, ID, cols)
  counts <- table(df[, 3:(length(cols) + 2)]) %>% as.data.frame() %>% 
    filter(Freq > 0)
  mat <- counts[, seq_len(ncol(counts) - 1)] == TRUE
  counts$conditions <- apply(mat, 1, function(x) {
    paste0(colnames(counts[, seq_len(ncol(counts) - 1)])[x], 
           collapse = " ")
  })
  counts$n_con <- apply(mat, 1, function(x) length(which(x)))
  counts <- counts %>% arrange(n_con, conditions)
  counts$conditions <- parse_factor(counts$conditions, levels = counts$conditions)
  counts <- mutate(counts, ID = paste0(conditions, ": ", Freq))
  if (nrow(counts) <= 40) {
    labelsize = 12
  }
  else {
    labelsize = 12/(nrow(counts)/40)
  }
  p1 <- ggplot(counts, aes(x = n_con, y = Freq, fill = conditions)) + 
    geom_col() + geom_text(aes(label = ID), size = 3, hjust = 0.5, 
                           position = position_stack(vjust = 0.5)) + labs(x = "Number of significant conditions", 
                                                                          y = "Number of proteins", title = "Overlap between conditions") + 
    theme_DEP1()
  ttheme <- gridExtra::ttheme_minimal(core = list(fg_params = list(hjust = 0, 
                                                                   x = 0.1)))
  p2 <- gridExtra::tableGrob(legend, theme = ttheme, rows = NULL, 
                             cols = NULL)
  if (!plot) {
    df <- counts %>% select(conditions, Freq, n_con)
    colnames(df) <- c("conditions", "proteins", "number_of_conditions")
    list <- list(counts = df, legend = legend)
    return(list)
  }
  else {
    gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(0.8, 
                                                         0.2))
  }
}

plot_volcano<-function (dep, contrast, label_size = 3, add_names = TRUE, adjusted = FALSE, 
          plot = TRUE) 
{
  if (is.integer(label_size)) 
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.character(contrast), length(contrast) == 1, is.numeric(label_size), 
                          length(label_size) == 1, is.logical(add_names), length(add_names) == 
                            1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot), 
                          length(plot) == 1)
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."), 
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), 
         call. = FALSE)
  }
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."), 
         call. = FALSE)
  }
  if (length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) == 
      0) {
    valid_cntrsts <- row_data %>% data.frame() %>% select(ends_with("_diff")) %>% 
      colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '", 
                                paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n", 
         valid_cntrsts_msg, call. = FALSE)
  }
  diff <- grep(paste(contrast, "_diff", sep = ""), colnames(row_data))
  if (adjusted) {
    p_values <- grep(paste(contrast, "_p.adj", sep = ""), 
                     colnames(row_data))
  }
  else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""), 
                     colnames(row_data))
  }
  signif <- grep(paste(contrast, "_significant", sep = ""), 
                 colnames(row_data))
  df <- data.frame(x = row_data[, diff], y = -log10(row_data[, 
                                                             p_values]), significant = row_data[, signif], name = row_data$name) %>% 
    filter(!is.na(significant)) %>% arrange(significant)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  p <- ggplot(df, aes(x, y)) + geom_vline(xintercept = 0) + 
    geom_point(aes(col = significant)) + geom_text(data = data.frame(), 
                                                   aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1, 
                                                                                                      0), vjust = c(-1, -1), label = c(name1, name2), size = 5, 
                                                       fontface = "bold")) + labs(title = contrast, x = expression(log[2] ~ 
                                                                                                                     "Fold change")) + theme_DEP1() + theme(legend.position = "none") + 
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey"))
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = filter(df, significant), 
                                      aes(label = name), size = label_size, box.padding = unit(0.1, 
                                                                                               "lines"), point.padding = unit(0.1, "lines"), 
                                      segment.size = 0.5)
  }
  if (adjusted) {
    p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  }
  else {
    p <- p + labs(y = expression(-log[10] ~ "P-value"))
  }
  if (plot) {
    return(p)
  }
  else {
    df <- df %>% select(name, x, y, significant) %>% arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", 
                                  "p_value_-log10")
    if (adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}



plot_single<-function (dep, proteins, type = c("contrast", "centered"), plot = TRUE) 
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.character(proteins), is.character(type), is.logical(plot), 
                          length(plot) == 1)
  type <- match.arg(type)
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("label", "condition", "replicate") %in% colnames(colData(dep)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '", 
         deparse(substitute(dep)), "'\nRun make_se() or make_se_parse() to obtain the required columns", 
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", 
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns", 
         call. = FALSE)
  }
  if (!"name" %in% colnames(row_data)) {
    stop("'name' column not present in '", deparse(substitute(dep)), 
         "'\nRun make_se() or make_se_parse() to obtain the required columns", 
         call. = FALSE)
  }
  if (all(!proteins %in% row_data$name)) {
    if (length(proteins) == 1) {
      rows <- grep(substr(proteins, 1, nchar(proteins) - 
                            1), row_data$name)
      possibilities <- row_data$name[rows]
    }
    else {
      rows <- lapply(proteins, function(x) grep(substr(x, 
                                                       1, nchar(x) - 1), row_data$name))
      possibilities <- row_data$name[unlist(rows)]
    }
    if (length(possibilities) > 0) {
      possibilities_msg <- paste0("Do you mean: '", paste0(possibilities, 
                                                           collapse = "', '"), "'")
    }
    else {
      possibilities_msg <- NULL
    }
    stop("please run `plot_single()` with a valid protein names in the 'proteins' argument\n", 
         possibilities_msg, call. = FALSE)
  }
  if (any(!proteins %in% row_data$name)) {
    proteins <- proteins[proteins %in% row_data$name]
    warning("Only used the following protein(s): '", paste0(proteins, 
                                                            collapse = "', '"), "'")
  }
  subset <- dep[proteins]
  if (type == "centered") {
    means <- rowMeans(assay(subset), na.rm = TRUE)
    df_reps <- data.frame(assay(subset) - means) %>% rownames_to_column() %>% 
      gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)), 
                                              by = "ID")
    df_reps$replicate <- as.factor(df_reps$replicate)
    df <- df_reps %>% group_by(condition, rowname) %>% summarize(mean = mean(val, 
                                                                             na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
      mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean - 
               error, CI.R = mean + error) %>% as.data.frame()
    df$rowname <- parse_factor(df$rowname, levels = proteins)
    p <- ggplot(df, aes(condition, mean)) + geom_hline(yintercept = 0) + 
      geom_col(colour = "black", fill = "grey") + geom_point(data = df_reps, 
                                                             aes(condition, val, col = replicate), shape = 18, 
                                                             size = 5, position = position_dodge(width = 0.3)) + 
      geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) + 
      labs(x = "Baits", y = expression(log[2] ~ "Centered intensity" ~ 
                                         "(±95% CI)"), col = "Rep") + facet_wrap(~rowname) + 
      theme_DEP2()
  }
  if (type == "contrast") {
    df <- rowData(subset, use.names = FALSE) %>% data.frame() %>% 
      select(name, ends_with("_diff"), ends_with("_CI.L"), 
             ends_with("_CI.R")) %>% gather(var, val, -name) %>% 
      mutate(contrast = gsub("_diff|_CI.L|_CI.R", "", var), 
             var = gsub(".*_", "", var)) %>% spread(var, val)
    df$name <- parse_factor(df$name, levels = proteins)
    suffix <- get_suffix(df$contrast)
    if (length(suffix)) {
      df$contrast <- delete_suffix(df$contrast)
    }
    p <- ggplot(df, aes(contrast, diff)) + geom_hline(yintercept = 0) + 
      geom_col(colour = "black", fill = "grey") + geom_errorbar(aes(ymin = CI.L, 
                                                                    ymax = CI.R), width = 0.3) + labs(x = suffix, y = expression(log[2] ~ 
                                                                                                                                   "Fold change" ~ "(±95% CI)")) + facet_wrap(~name) + 
      theme_DEP2()
  }
  if (plot) {
    return(p)
  }
  else {
    if (type == "centered") {
      df <- df %>% select(rowname, condition, mean, CI.L, 
                          CI.R)
      colnames(df) <- c("protein", "condition", "log2_intensity", 
                        "CI.L", "CI.R")
    }
    if (type == "contrast") {
      df <- df %>% select(name, contrast, diff, CI.L, CI.R) %>% 
        mutate(contrast = paste0(contrast, suffix))
      colnames(df) <- c("protein", "contrast", "log2_fold_change", 
                        "CI.L", "CI.R")
    }
    return(df)
  }
}


plot_missval <-function (se) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), 
         "'", call. = FALSE)
  }
  df <- se_assay %>% data.frame(.)
  missval <- df[apply(df, 1, function(x) any(is.na(x))), ]
  missval <- ifelse(is.na(missval), 0, 1)
  ht2 = Heatmap(missval, col = c("white", "black"), column_names_side = "top", 
                show_row_names = FALSE, show_column_names = TRUE, name = "Missing values pattern", 
                column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(at = c(0, 
                                                                                          1), labels = c("Missing value", "Valid value")))
  draw(ht2, heatmap_legend_side = "top")
}