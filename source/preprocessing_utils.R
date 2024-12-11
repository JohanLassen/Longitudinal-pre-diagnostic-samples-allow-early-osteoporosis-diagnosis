



#' Transform
#'
#' Wrapper for all the transformation methods. Only used in the app.
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param method Which transformation method to be used (includes fourth root and log10)
#'
#' @return ms object
transform <- function(ms, method, print_methods = FALSE){
  
  methods <- list("log10"=transform_log, "fourth root" = transform_fourth_root)
  
  if (print_methods){
    print(names(methods))
    return()
  }
  
  msX <- ms
  msX$values[is.na(msX$values)] <- 0
  
  msX <- methods[[method]](msX)
  return(msX)
}


#' Impute
#'
#' Imputes missing values. Only used for the app
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param method Which imputation method to be used. Includes knn or 0-imputation
#'
#' @return ms object
impute <- function(ms, method, print_methods=FALSE){
  if (method == "default") return(ms)
  methods <- list("knn" = impute_knn)
  if (print_methods){
    print(c(names(methods), "default"))
  }
  ms$values[ms$values == 0] <- NA # we imputed by 0 and now we reverse to make a more advanced imputation
  ms <- methods[[method]](ms)
  return(ms)
}


#' Normalize
#'
#' Normalizes the data. Only used in the app.
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param method Which normalization method to be used.
#'
#' @return ms
normalize <- function(ms, method){
  
  methods <- list("PQN" = normalize_pqn,
                  "Rownorm"=normalize_robust_row,
                  "Quantile normalization"=normalize_limma_quantile,
                  "WaveICA2"=normalize_waveICA2,
                  "Combat" = normalize_combat)
  
  if (print_methods){
    print(names(methods))
    return()
  }
  
  ms <- methods[[method]](ms)
  return(ms)
}


#' Fourth root transform data
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- transform_fourth_root(ms)
#' head(ms$values)
#' head(ms$rowinfo)

transform_fourth_root <- function(ms){
  ms$values     <- ms$values^0.25 %>% tibble::as_tibble()
  return(ms)
}


#' Log transform data
#' log(x+1)
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- transform_log(ms)
#' head(ms$values)
#' head(ms$rowinfo)
transform_log <- function(ms){
  ms$values     <- log(ms$values+1) %>% tibble::as_tibble()
  return(ms)
}

#' Missing value imputation: Zero replacement
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' head(ms$values)
#' head(ms$rowinfo)
impute_zero <- function(ms){
  ms$values[is.na(ms$values)] <- 0
  return(ms)
}


#' Missing value imputation: knn replacement
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_knn(ms)
#' head(ms$values)
#' head(ms$rowinfo)
impute_knn <- function(ms){
  tmp1          <- ms$values %>% as.matrix() %>% t()
  tmp1[tmp1==0] <- NA
  tmp1_imputed  <- impute::impute.knn(tmp1 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
  tmp1          <- tmp1_imputed$data %>% t() %>% tibble::as_tibble()
  ms$values     <- tmp1
  return(ms)
}


# impute_svd <- function(ms, k = 20){
#   feature_names <- colnames(ms$values)
#   tmp1          <- ms$values %>% as.matrix()
#   tmp1[tmp1==0] <- NA
#   tmp2          <- bcv::impute.svd(tmp1, k=k)
#   ms$values     <- tmp2$x %>% tibble::as_tibble()
#   colnames(ms$values) <- feature_names
#   return(ms)
# }

#
# #' Missing value imputation: Random Forest replacement
# #'
# #' Beware! This algorithm might run longer than your lifespan.
# #'
# #' @param ms List object containing value data and rowinfo data (meta data)
# #'
# #' @return ms
# #' @export
# #'
# #' @examples
# #' # First convert the pneumonia object to a tibble.
# #' pneumonia<- tibble::tibble(pneumonia)
# #'
# #' # Generate list object
# #' ms <- list()
# #'
# #' # Assign feature values to ms$values
# #' start_column <- 8 # The first column with feature values
# #' end_column <- ncol(pneumonia) # The last column of the dataset
# #' ms$values <- pneumonia[1:10, start_column:end_column]
# #'
# #' # Assign metadata to ms$rowinfo
# #' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
# #' ms <- impute_zero(ms)
# #' ms <- impute_RF(ms)
# #' head(ms$values)
# #' head(ms$rowinfo)
# impute_RF <- function(ms){
#   tmp1          <- ms$values %>% as.matrix()
#   tmp1[tmp1==0] <- NA
#   doParallel::registerDoParallel(cores=8)
#   tmp1_imputed  <-
#     missForest::missForest(
#       tmp1, maxiter = 10, ntree = 100, variablewise = FALSE, verbose = TRUE, parallelize = "variables"
#     )
#
#   tmp1          <- tmp1_imputed$data %>% t() %>% tibble::as_tibble()
#   ms$values     <- tmp1
#   return(ms)
# }



#' Remove features inflated with zeros (batch-wise)
#'
#' This might help avoiding batch effect and identification of batches through signals with binary (on/off) behaviour
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param batch column name of batch. If none provided the dataset is considered as one large batch.
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#' # Generate list object
#' ms <- list()
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- rm_feature_batch_inflated_zeros(ms)
#' head(ms$values)
#' head(ms$rowinfo)

rm_feature_batch_inflated_zeros <- function(ms, batch = NULL){
  
  raw     <- ms$values
  rowinfo <- ms$rowinfo
  
  if (is.null(batch)){
    BATCH = 1
  } else{
    BATCH = rowinfo[[batch]]
  }
  
  tmp1 <- rowinfo %>%
    dplyr::mutate(BATCH = BATCH) %>%
    dplyr::select(BATCH) %>%
    dplyr::bind_cols(raw) %>%
    tidyr::pivot_longer(cols = starts_with("M")) %>%
    dplyr::group_by(BATCH, name) %>%
    dplyr::summarise(zeros = mean(value==0))
  
  bad_features <- tmp1[tmp1$zeros>0.2,] %>% pull(name) %>% table %>%  names() %>%  unique()
  
  raw <- raw %>% dplyr::select(-dplyr::any_of(bad_features))
  
  ms$values <- raw
  ms$rowinfo <- ms$rowinfo
  
  return(ms)
}


#' Removes observations if their batches are below a given threshold
#'
#' Is usable if preceeded by batch normalization
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param min_occurence the minimum number of samples per batch (to avoid batch effect being predictive for the outcome)
#' @param batch Column name of the batch info column
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- rm_sample_min_batch_occurence(ms)
#' head(ms$values)
#' head(ms$rowinfo)


rm_sample_min_batch_occurence <- function(ms, min_occurence = 5, batch = NULL){
  
  rowinfo <- ms$rowinfo
  values <- ms$values
  
  
  if (is.null(batch)){
    BATCH = 1
  } else{
    BATCH = rowinfo[[batch]]
  }
  
  
  small_batches <- table(rowinfo$BATCH)[table(rowinfo$BATCH) < min_occurence] %>% names
  rowinfo <- rowinfo %>% dplyr::filter(!BATCH %in% small_batches)
  values <- values[rowinfo$rowid,]
  rowinfo <- rowinfo %>% dplyr::mutate(rowid = dplyr::row_number())
  
  ms$rowinfo <- rowinfo
  ms$values <- values
  
  return(ms)
}


#' Removes outliers by PCA
#'
#' Samples are removed if: abs(x - median(x)) > (1.5 * quantile(x,0.95)-quantile(x,0.05)), where x is the PCA scores.
#' Scores from the first 12 components are used.
#'
#' May be applied before row normalization but after removal of extreme features.
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param plot set to TRUE to visualize the removed samples
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- rm_sample_pca_outliers(ms)
rm_sample_pca_outliers <- function(ms, plot = FALSE) {
  print(ms)
  
  
  raw     <- ms$values
  rowinfo <- ms$rowinfo
  
  
  tmp1 <- rowinfo %>% dplyr::mutate(rowid2 = dplyr::row_number()) #
  tmp2 <- raw[tmp1$rowid,] %>% as.matrix()
  
  r  <- prcomp(x = tmp2, retx = T, center=T, scale. = T, rank. = 12)
  
  bad_rows<- tibble::tibble(rowid2=apply(r$x, 2, function(x) {
    which(abs(x - median(x)) > (1.5 * quantile(x,0.95)-quantile(x,0.05)))
  }) %>%
    unlist() %>%
    as.vector()) %>%
    dplyr::count(rowid2)
  
  tmp1 <- tmp1 %>%
    dplyr::left_join(bad_rows) %>%
    dplyr::mutate(n=ifelse(is.na(n), 0,n))
  
  pd <- r$x %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(tmp1) %>%
    {.}
  
  pd <- pd %>%
    dplyr::mutate(response = ifelse(n>0,"Outlier", "Not outlier")) %>%
    dplyr::mutate(response = factor(response))
  
  plotlist <- list()
  
  for(i in 1:(ncol(r$x)/2)) {
    xvar <- names(pd)[2*i-1]
    yvar <- names(pd)[2*i]
    p1 <- ggplot2::ggplot(pd,ggplot2::aes_string(x=xvar, y=yvar,
                                                 fill="response"))+
      ggplot2::geom_point(shape=21, color="#FFFFFFFF", size=3) +
      ggplot2::scale_fill_manual(values = c("#D0D0D0", "#D04040")) +
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.position="none") +
      NULL
    
    plotlist[[length(plotlist)+1]] <- p1
    rm(p1)
  }
  
  outlier_plot <- cowplot::plot_grid(plotlist = plotlist,nrow=1)
  
  bad_rows <- tmp1 %>% dplyr::filter(n>0)
  if (nrow(bad_rows) > 0) {
    ms$values  <- raw[-bad_rows$rowid2,] %>% tibble::as_tibble()
    ms$rowinfo <- rowinfo[-bad_rows$rowid2,]
  }
  
  
  ms$rowinfo <- ms$rowinfo %>%
    dplyr::mutate(rowid = dplyr::row_number())
  
  if (plot) return(list("ms"=ms, "outlier_plot"=outlier_plot))
  return(ms)
}


#' Remove samples if internal standards are outliers
#'
#'
#' @param quantiles {How extreme should an IS value be before being assigned as an outlier.}
#' @param ms {List object containing value data and rowinfo data (meta data)}
#' @param standards {Vector containing column names of the internal standards e.g., M215T123}
#' @param tolerance {Number of how many internal standards that must lie outside 95 percent of the distribution before a sample is excluded.
#' Max number of standards
#' Min 1
#' The smaller the more restrictive. The tolerance reflects the number of outlier ISs allowed per sample.
#' Outlier is given by lying in the 5 percent most extreme values. This means that if a threshold of 1 is used 5 percent of samples are removed.
#' This should be used carefully.}
#' @param batch column name of batch info
#'
#' @return ms
#' @export
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia <- tibble::tibble(pneumonia)
#' # Generate list object
#' ms <- list()
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia %>%
#'     dplyr::select(id, group, age, gender) %>%
#'     dplyr::slice(1:10) %>%
#'     dplyr::mutate(rowid = dplyr::row_number())
#' ms <- impute_zero(ms)
#' ms <- rm_IS_outliers(ms, standards=c("M363T419","M512T603"), tolerance=4, quantiles=c(0.01,0.99))
#'
rm_IS_outliers <- function(ms, standards, tolerance = 0, quantiles = c(0.025, 0.975), batch = NULL){
  tmp1 <- ms$rowinfo
  tmp2 <- ms$values
  
  if (is.null(batch)){
    batch = rep(1, length(tmp1$rowid))
  } else{
    batch = tmp1[[batch]]
  }
  tmp1 <- tmp1 %>% dplyr::mutate(BATCH = batch)
  good_data <-
    tmp2 %>%
    dplyr::select(all_of(standards)) %>%
    dplyr::bind_cols(tmp1) %>%
    tidyr::pivot_longer(cols = all_of(standards)) %>%
    dplyr::mutate(batch = as.factor(BATCH))
  
  good_data <-
    good_data %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(outlier = value >= quantile(value, quantiles[2]) | value <= quantile(value, quantiles[1])) %>%
    dplyr::ungroup()
  
  a <-
    ggplot2::ggplot(good_data) +
    ggplot2::geom_point(ggplot2::aes(y=value, x=rowid, color = batch), show.legend = F)+
    ggplot2::geom_point(data = good_data %>% dplyr::filter(outlier), ggplot2::aes(y=value, x=rowid), color = "gray40")+
    ggplot2::facet_wrap(~name, ncol = 1) +
    NULL
  
  bad_samples <-
    good_data %>%
    dplyr::filter(outlier == TRUE) %>%
    dplyr::count(rowid) %>%
    dplyr::filter(n >= tolerance) %>%
    dplyr::pull(rowid)
  
  ms$rowinfo <- ms$rowinfo %>% dplyr::filter(!(rowid %in% bad_samples))
  ms$values <- ms$values[ms$rowinfo$rowid,]
  ms$rowinfo <- ms$rowinfo %>% dplyr::mutate(rowid = dplyr::row_number())
  
  return(ms)
}


#' Removes outliers by PCA
#'
#' Samples are removed if: abs(x - median(x)) > (1.5 * quantile(x,0.95)-quantile(x,0.05)), where x is the PCA scores.
#' Scores from the first 12 components are used.
#'
#' May be applied before row normalization but after removal of extreme features.
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @examples
#' # First convert the pneumonia object to a tibble. 
#' pneumonia<- tibble::tibble(pneumonia)
#' # Generate list object
#' ms <- list()
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- rm_feature_extreme_values(ms)
rm_feature_extreme_values <- function(ms) {
  
  raw     <- ms$values
  rowinfo <- ms$rowinfo
  
  tmp1 <-
    tibble::tibble(rowid = rowinfo$rowid) %>%
    dplyr::bind_cols(tibble::as_tibble(raw))
  
  tmp1 <-
    tmp1 %>%
    tidyr::pivot_longer(names_to = "compound", values_to = "value",  cols= c(-rowid))
  
  tmp2 <- tmp1 %>%
    dplyr::group_by(compound) %>%
    dplyr::summarise(n_bad = sum(value > median(value)+2*quantile(value,0.95))) %>%
    {.}
  
  bad_features <- tmp2 %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_bad > 0) %>%
    dplyr::select(compound) %>%
    dplyr::distinct()
  
  ms$values <- raw %>% dplyr::select(-dplyr::any_of(bad_features$compound))
  ms$rowinfo <- rowinfo
  
  return(ms)
}



#' Robust row normalization
#'
#' Only stable features used to normalize the rows (each sample)
#' Corrects matrix effects
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' library(dplyr)
#' pneumonia<- tibble::tibble(pneumonia)
#' # Generate list object
#' ms <- list()
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:100, start_column:end_column]
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:100,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- normalize_robust_row(ms)
normalize_robust_row <- function(ms) {
  
  target_info   <- ms$rowinfo
  target_values <- ms$values %>% tibble::as_tibble()
  
  tmp1 <- target_info
  tmp2 <- target_values[tmp1$rowid,]
  
  stable_features <-
    tmp1 %>%
    dplyr::bind_cols(tmp2) %>%
    tidyr::pivot_longer(dplyr::starts_with("M")) %>%
    dplyr::group_by(rowid) %>%
    dplyr::mutate(rank=rank(value)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(median = median(rank),
                     range = max(rank)-min(rank)) %>%
    dplyr::ungroup() %>%
    dplyr::slice_min(order_by = median, prop = 0.8) %>%
    dplyr::slice_max(order_by = median, prop = 0.8) %>%
    dplyr::slice_min(order_by = range, prop = 0.8)
  
  raw    <- target_values
  data.x <- raw
  tmp    <- rowSums(target_values %>% dplyr::select(dplyr::any_of(stable_features$name)))
  raw    <- max(raw)*raw / tmp
  
  ms$values  <- as_tibble(raw)
  ms$rowinfo <- target_info
  
  return(ms)
}



#' Probabilistic quotient normalization
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param n Summary function
#' @param QC Available QC column?
#'
#' @return ms
#' @export
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- normalize_pqn(ms)
#' head(ms$values)
#' head(ms$rowinfo)
normalize_pqn <- function(ms, n = "median", QC = NULL) {
  X <- ms$values %>% as.matrix()
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  
  if (!is.null(QC)) {
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }
  
  # do the actual normalisation
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ] / median(as.numeric(X[a, ] / mX)))
  }
  
  ms$values <- X.norm %>% tibble::as_tibble()
  return(ms)
}


#' Normalize using waveICA2
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- normalize_waveICA2(ms)
#' head(ms$values)
#' head(ms$rowinfo)
#'
normalize_waveICA2 <- function(ms){
  tmp1 <- ms$rowinfo
  tmp2 <- ms$values %>% as.matrix()
  tmp3 <- WaveICA_2.0(tmp2, Injection_Order = 1:nrow(tmp2), Cutoff = 0.1, wf = "haar", K = 20, alpha = 0)
  tmp3 <- tmp3$data_wave
  ms$values <- tmp3 %>% tibble::as_tibble()
  return(ms)
}

#' Normalize using combat
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#'

normalize_combat <- function(ms){
  ms$values <-
    ComBat(t(ms$values), batch = as.numeric(as.factor(ms$rowinfo$BATCH))) %>%
    t() %>%
    tibble::as_tibble()
  return(ms)
}

#' Quantile normalize using limma (quantile normalization)
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- normalize_limma_quantile(ms)
#' head(ms$values)
#' head(ms$rowinfo)
#'
normalize_limma_quantile <- function(ms){
  ms$values <- limma::normalizeBetweenArrays(as.matrix(ms$values), method = "quantile") %>% tibble::as_tibble()
  return(ms)
}


#' Cyclic loss normalize using limma
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#'
#' @return ms
#' @export
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- normalize_limma_cyclicloess(ms)
#' head(ms$values)
#' head(ms$rowinfo)
#'
normalize_limma_cyclicloess <- function(ms){
  ms$values <- limma::normalizeBetweenArrays(ms$values, method = "cyclicloess")
  return(ms)
}


#' Select most variable features (beta version)
#' Only works with features starting with "M" as in M314T182 (untargeted data)
#' Effective for model screening to avoid heavy work loads. Often top 500 contains most of the signal.
#'
#' @param ms List object containing value data and rowinfo data (meta data)
#' @param n the number of most variable features selected
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <-select_most_variable_features(ms, n=10)
#' head(ms$values)
#' head(ms$rowinfo)
#'
select_most_variable_features <- function(ms, n=500) {
  raw <- ms$values
  rowinfo <- ms$rowinfo
  
  tmp1 <- rowinfo
  tmp2 <- raw[tmp1$rowid,]
  
  good_features <- tmp2 %>%
    dplyr::bind_cols(tmp1) %>%
    tidyr::pivot_longer(dplyr::starts_with("M")) %>% # Convert to tidy format
    dplyr::group_by(name) %>%
    dplyr::summarise(mean = mean(value),
                     variation = var(value)/mean(value)) %>%
    dplyr::arrange(-variation) %>%
    dplyr::slice(1:n) %>%
    pull(name)
  
  ms$rowinfo <- rowinfo
  ms$values  <- raw %>% dplyr::select(dplyr::any_of(good_features))
  
  return(ms)
}



#' Z-standardization of batches
#'
#'
#' @param ms List object containing value and meta data
#' @param batch_column Column name of batch info
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' ms <- standardize_z_batch(ms)
#' head(ms$values)
#' head(ms$rowinfo)
#'
standardize_z_batch <- function(ms, batch_column = NULL) {
  raw <- ms$values
  target_info <- ms$rowinfo
  
  if (is.null(batch_column)){
    batch = rep(1, length(target_info$rowid))
  } else{
    batch = target_info[[batch_column]]
  }
  
  tmp1<- tibble::tibble(batch = batch, rowid = target_info$rowid) %>%
    dplyr::bind_cols(tibble::as_tibble(raw)) %>%
    tidyr::pivot_longer(names_to = "compound", values_to = "value",  cols= c(-batch, -rowid))
  
  tmp2 <- tmp1 %>%
    dplyr::group_by(batch, compound) %>%
    dplyr::summarise(value.batch_mean = mean(value),
                     value.batch_sd   = sd(value))
  
  
  tmp1 <- dplyr::full_join(x = tmp1, y = tmp2, all.x = TRUE) %>%
    dplyr::mutate(value2 = (value- value.batch_mean)/(value.batch_sd)) %>% #
    dplyr::mutate(value2 = ifelse((
      value.batch_mean==0 &
        value.batch_sd==0) |
        value == 0, 0, value2
    )) %>%
    dplyr::select(batch, rowid, compound,value2) %>%
    tidyr::pivot_wider(names_from = compound, values_from = value2)
  
  raw <- tmp1 %>%
    dplyr::arrange(rowid) %>%
    dplyr::select(-batch, -rowid)
  
  ms$rowinfo <- target_info
  ms$values  <- raw
  return(ms)
}


#' Plot intensity vs injection order
#'
#' @param ms List object of value and meta data
#' @param standards vector containing column names of internal standards
#' Feature names of standards
#' @param batch column name of batch info
#'
#' @return ggplot2 plot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#'
#'
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:10, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:10,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' plot_standards(ms, standards = c("M363T419","M512T603","M364T419","M365T392", "M143T177"))
plot_standards <- function(ms, standards, batch = NULL) {
  values <- ms$values
  info   <- ms$rowinfo
  
  if (is.null(batch)){
    BATCH = 1
  } else{
    BATCH = info[[batch]]
  }
  
  plot_data <-
    values %>%
    dplyr::select(dplyr::all_of(standards)) %>%
    dplyr::bind_cols(info) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(standards)) %>%
    dplyr::mutate(batch = as.factor(as.numeric(as.factor(BATCH)) %% 2))
  
  p <-
    ggplot2::ggplot(plot_data) +
    ggplot2::geom_point(ggplot2::aes(y=value, x=rowid, fill = batch, color = batch), shape = 21, show.legend = F)+
    ggplot2::facet_wrap(~name, ncol = 1) +
    ggplot2::scale_fill_manual(values = c("#ffffcc", "#41b6c4"))+
    ggplot2::scale_color_manual(values = c("#7fcdbb","#225ea8")) +
    ggplot2::theme_classic()
  NULL
  
  return(p)
}


#' UMAP of data (beta)
#'
#' @param tmp1 metadata
#' @param tmp2 valuedata
#' @param color_label coloring variable
#'
#' @return ggplot2 plot
#' @export
#' @import ggplot2
#' @import dplyr
#' @examples
#' # First convert the pneumonia object to a tibble.
#' pneumonia<- tibble::tibble(pneumonia)
#'
#' # Generate list object
#' ms <- list()
#'
#' # Assign feature values to ms$values
#' start_column <- 8 # The first column with feature values
#' end_column <- ncol(pneumonia) # The last column of the dataset
#' ms$values <- pneumonia[1:100, start_column:end_column]
#'
#' # Assign metadata to ms$rowinfo
#' ms$rowinfo <- pneumonia[1:100,] %>% dplyr::select(id, group, age, gender)
#' ms <- impute_zero(ms)
#' tmp1 <- ms$rowinfo
#' tmp2 <- ms$values
#' make_umap(tmp1, tmp2, color_label = "group")
make_umap <- function(tmp1, tmp2, color_label) {
  r <- umap::umap(tmp2)
  o<- tibble::tibble("umap1"=r$layout[,1],
                     "umap2"=r$layout[,2])
  pd <- o %>%
    dplyr::bind_cols(tmp1) %>%
    {.}
  
  p <- ggplot2::ggplot(pd,ggplot2::aes(x=umap1, y=umap2,
                                       fill=.data[[color_label]]))+
    ggplot2::geom_point(shape=21, color="#FFFFFFFF", size=3, show.legend = F) +
    ggplot2::labs(fill = color_label)+
    ggplot2::theme() +
    NULL
  return(p)
}


#' Importance from caret model
#'
#' @param model the caret model trained on the ms list
#'
#' @return ggplot2 plot
#' @export
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import ggplot2
#' @import caret
plot_importance <- function(model){
  
  # FEATURE IMPORTANCE - BEST MODEL
  importance <- caret::varImp(model) %>%
    {
      tibble::tibble(Feature = rownames(.$importance),
                     Importance = .$importance %>% unlist %>% round(2)) %>%
        dplyr::arrange(-Importance)
    }
  
  
  importance[1:10,] %>%
    ggplot2::ggplot(ggplot2::aes(y=factor(Feature, levels = rev(importance$Feature)), x=Importance)) +
    ggplot2::geom_col()+
    ggplot2::labs(y = "Features")+
    ggplot2::theme_minimal()
  
  important <- importance %>% dplyr::filter(Importance>5) %>% dplyr::arrange(-Importance) %>% dplyr::slice(1:25)
  
  
  
  outcome <- model %>% dplyr::select(obs, rowIndex) %>% dplyr::distinct() %>% dplyr::arrange(rowIndex)
  
  outcome %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(model$trainingData) %>%
    dplyr::select(obs, all_of(important$Feature)) %>%
    dplyr::filter(!is.na(obs)) %>%
    tidyr::pivot_longer(names_to = "Feature", cols=-obs) %>%
    ggplot2::ggplot() +
    ggplot2::geom_boxplot(ggplot2::aes(x=obs, y=value), fill = "#6DA398", show.legend = F)+
    ggplot2::facet_wrap(~Feature, scales = "free_y")+
    ggplot2::theme_bw()
}


#' Scatter plot of probabilities
#'
#' Used when classifying over a continuous variable
#'
#' @param model the caret model trained on the ms list
#' @param target_info outcome variable
#' @param var_name outcome variable
#' @import dplyr
#' @import ggplot2
#' @return NULL
#' @export
prob_vs_continuous <- function(model, target_info, var_name){
  pdata <-
    model$pred %>%
    dplyr::arrange(.data$rowIndex) %>%
    dplyr::mutate(rowid = .data$rowIndex) %>%
    dplyr::inner_join(target_info %>% dplyr::drop_na() %>% dplyr::mutate(rowid = dplyr::row_number()), by = "rowid")
  
  classes <- unique(model_cv_r$obs)
  
  ggplot2::ggplot(pdata, ggplot2::aes_string(y=classes[1], x = var_name, fill = "obs")) +
    ggplot2::geom_point(shape = 21, color = "white") +
    ggplot2::scale_fill_manual(values = c("#d53e4f", "#3288bd"))+
    ggplot2::geom_segment(aes(x = min(age), xend = max(age), y=0.5, yend = 0.5), color = "black", linetype = "dashed")
}

