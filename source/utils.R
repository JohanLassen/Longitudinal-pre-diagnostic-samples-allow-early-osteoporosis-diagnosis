

plot_pca <- function(feature_matrix, metadata, color_labels="rowid", palette = "Spectral", tech_rep = "rowid", panels = 3) {
  
  # Feature matrix --> features on columns and samples on rows. Only features, no metadata
  # metadata --> only metadata. used for coloring
  # palette --> palette used for coloring
  # tech_rep --> if the data contain techical replicates. If the replicates are duplicates: c(1, 1, 2, 2, etc,) - i.e., one value per group
  
  cowplot_labels <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)")
  
  tmp1 <- metadata
  tmp2 <- feature_matrix
  if (!"rowid"%in%colnames(tmp1)){tmp1 <- tmp1 %>% tibble::rowid_to_column()}
  if (length(palette)==1) {palette <- rep(palette, length(color_labels))}
  
  r    <- prcomp(x = tmp2, retx = TRUE, center = T, scale = T, rank. = panels*2)
  
  variance_explained <- summary(r)$importance[2,1:(panels*2)]
  variance_explained <- round(variance_explained, 3)*100
  
  variable_plotlist <- list()
  
  for (color_idx in seq_along(color_labels)){
    color_label = color_labels[color_idx]
    
    pd <- r$x %>%
      tibble::as_tibble() %>%
      dplyr::bind_cols(tmp1 %>% dplyr::select(dplyr::all_of(c(color_label, tech_rep)))) %>%
      {.}
    
    plotlist <- list()
    for(i in 1:(ncol(r$x)/2)) {
      
      xvar <- names(pd)[2*i-1]
      yvar <- names(pd)[2*i]
      p1 <-
        ggplot2::ggplot(pd, ggplot2::aes(x=.data[[xvar]], y=.data[[yvar]], fill=.data[[color_label]], group=.data[[tech_rep]]))+
        ggplot2::geom_polygon(aes(fill=.data[[color_label]]), show.legend = F, alpha = 0.6)+
        #ggplot2::geom_line(aes(color=.data[[color_label]]), show.legend = F)+
        ggplot2::geom_point(shape=21, color = "white", size=2, show.legend = F) +#, color="black"
        ggplot2::labs(fill = fill, 
                      x = paste0(xvar, " (", variance_explained[xvar], " %)"), 
                      y = paste0(yvar, " (", variance_explained[yvar], " %)"))+
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.title.x = element_text(size = 8),
                       axis.title.y = element_text(size = 8),
                       axis.text.x = element_text(size = 8),
                       axis.text.y = element_text(size = 8))
      
      if (!is.numeric(tmp1[[color_label]])){
        p1 <- p1+scale_fill_brewer(palette = palette[color_idx])+scale_color_brewer(palette = palette[color_idx])
      } else{
        p1 <- p1+scale_fill_distiller(palette = palette[color_idx])+scale_color_distiller(palette = palette[color_idx])
      }
      plotlist[[length(plotlist)+1]] <- p1
      
      if (i == 1){
        p1 <-
          ggplot2::ggplot(pd, ggplot2::aes_string(x=xvar, y=yvar, fill=color_label))+
          ggplot2::geom_point(shape=21, color="#FFFFFFFF", size=3)+
          ggplot2::theme(legend.text = element_text(size=8), legend.title = element_text(size=8), legend.key=element_blank()) +
          #ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 2)))+
          NULL
        
        if (!is.numeric(tmp1[[color_label]])){
          p1 <- p1+scale_fill_brewer(palette = palette[color_idx])
        } else{
          p1 <- p1+scale_fill_distiller(palette = palette[color_idx])
        }
        legend <- 
          cowplot::get_legend(
            p1 + 
              ggplot2::theme(legend.box.margin = ggplot2::margin(0, 0, 0, 0)))
      }
    }
    
    p1 <- cowplot::plot_grid(plotlist = plotlist, nrow=1)
    p1 <- cowplot::plot_grid(p1, legend, rel_widths = c(3, .4))
    if (color_idx==1){
      title <- cowplot::ggdraw() + cowplot::draw_label(paste("PCA based on", ncol(tmp2), "Features"), size = 8)
    }
    
    variable_plotlist[[color_idx]] <- p1
    
  }
  p1    <- 
    cowplot::plot_grid(
      plotlist = variable_plotlist, 
      labels = cowplot_labels[1:length(color_labels)], 
      label_size = 10,
      label_fontface = "plain",
      ncol=1)
  final <- cowplot::plot_grid(title, p1, nrow=2, rel_heights = c(1, 10*length(variable_plotlist)))
  
  return(final)
}


make_slopes <- function(ms) {
  rowinfo <- ms %>% select(!starts_with("M"))
  values <- ms %>% select(starts_with("M"), -MCS)
  
  tmp1 <- 
    rowinfo %>% 
    filter(!is.na(group)) %>% 
    group_by(group) %>% 
    mutate(n=n()) %>% 
    ungroup() %>% 
    filter(n>4) %>% 
    select(-n)
  
  tmp2 <- values[tmp1$rowid,]
  tmp2 <- tmp2[,tmp2 %>% map_lgl(~sum(is.na(.x))==0)]
  
  try2 <- 
    tmp1 %>% 
    select(group, time_point = age_at_donation, rowid, status, batch) %>% 
    bind_cols(tmp2) %>%
    select(-rowid) %>% 
    pivot_longer(cols = -c(group, time_point, status, batch)) %>%
    group_by(group, status, name) %>% 
    mutate(value = scale(value)) %>% 
    summarise(
      batch = batch[1],
      metric = sum((time_point-mean(time_point, na.rm = T))*(value-mean(value, na.rm = T)))/(sum((time_point-mean(time_point, na.rm = T))^2, na.rm = T))
    ) %>% ungroup() %>% 
    pivot_wider(names_from = name, values_from = metric) %>% 
    ungroup() %>% 
    select(-group)
  
  ms <- list()
  x <- try2 %>% select(starts_with("M"))
  x <- x[(x %>% map(is.na) %>% map_dbl(sum))==0]
  ms$values <- x
  ms$rowinfo <- try2 %>% select(-starts_with("M"))
  
  return(ms)
}


generate_xy <- function(objective, ms, exclude_samples="none", exclude_features="none"){
  
  if (objective == "simple"){
    ms$rowid <- 1:length(ms$rowid)
    tmp2     <- 
      ms %>% 
      select(!starts_with("M")) %>% 
      filter(type1=="sample") %>% 
      filter(time_point == 1)
    
    x          <- ms %>% select(starts_with("M"), -MCS) %>% .[tmp2$rowid, ]
    x          <- x[((x %>% map(is.na) %>% map_dbl(sum))==0)]
    y          <- tmp2$status
    return(list(x, y, tmp2))
  }
  
  if (objective == "slope"){
    rowinfo <- ms %>% select(!starts_with("M"))
    values <- ms %>% select(starts_with("M"), -MCS)%>% select(-any_of(exclude_features))
    
    tmp1 <- 
      rowinfo %>% 
      filter(!is.na(group)) %>% 
      group_by(group) %>% 
      mutate(n=n()) %>% 
      ungroup() %>% 
      mutate(identifyer = paste0(group, status)) %>% 
      filter(n>4, (!identifyer %in% c(exclude_samples))) %>% 
      select(-n)
    
    tmp2 <- values[tmp1$rowid,]
    tmp2 <- tmp2[,tmp2 %>% map_lgl(~sum(is.na(.x))==0)] 
    
    try2 <- 
      tmp1 %>% 
      select(group, time_point = age_at_donation, rowid, status, batch) %>% 
      bind_cols(tmp2) %>%
      select(-rowid) %>% 
      pivot_longer(cols = -c(group, time_point, status, batch)) %>%
      group_by(group, status, name) %>% 
      mutate(value = scale(value)) %>% 
      summarise(
        batch = batch[1],
        metric = 
          sum((time_point-mean(time_point, na.rm = T))*
                (value-mean(value, na.rm = T)))/
          (sum((time_point-mean(time_point, na.rm = T))^2, na.rm = T))
      ) %>% ungroup()%>% 
        pivot_wider(names_from = name, values_from = metric) %>%
        ungroup() %>%
        select(-group)
    
    excl_samples <- paste0(tmp1$group, tmp1$status)

    x <- try2 %>% select(starts_with("M"))
    x <- x[(x %>% map(is.na) %>% map_dbl(sum))==0]
    
    y <- try2 %>% pull(status)
    xy <- list(x, y, try2)
    return(xy)
  }
  return("Objective not found")
}

generate_folds <- function(info, variable){
  folds <- as.list(unique(as.character(info[[variable]]))) %>% map(~{which(.x != as.character(info[[variable]]))})
  names(folds) <- paste0("Fold", 1:length(folds))
  return(folds)
}

fit <- function(ms, 
                permute = F, 
                objective = "simple", 
                cv_method = "loocv", 
                method = "glmnet", 
                tuneLength = 5, 
                permutations = 1000,
                exclude_samples = "none",
                exclude_features = "none"){
  # cv_method: loocv or leave-one-group-out where cv_method must be a group variable in ms
  
  set.seed(1)
  
  xy <- generate_xy(objective = objective, ms, exclude_samples = exclude_samples, exclude_features = exclude_features)  
  
  if (cv_method != "loocv") {
    folds  <- generate_folds(info = xy[[3]], variable = cv_method)
    trctrl <- 
      caret::trainControl(
        #summaryFunction = twoClassSummary,
        method      = "cv",
        index = folds,
        verboseIter = !permute | exclude_samples != "none", 
        savePredictions = T, 
        classProbs = T
      )
  } else {
    trctrl <- 
      caret::trainControl(
        #summaryFunction = twoClassSummary,
        method      = "loocv",
        verboseIter = !permute | exclude_samples != "none", 
        savePredictions = T, 
        classProbs = T
      )
  }
  
  if(permute) {
    
    fit <- 
      as.list(1:permutations) %>% 
      map(
        ~{
          xy[[2]] <- sample(xy[[2]]) 
          
          caret::train(
              x = xy[[1]], 
              metric = "ROC",
              y = xy[[2]],
              tuneLength = tuneLength,
              trControl = trctrl,
              method = method
            )
        },
        .progress = TRUE
      )

    
    
  } else {
    fit <- 
        caret::train(
          x = xy[[1]], 
          metric = "ROC",
          y = xy[[2]],
          tuneLength = tuneLength,
          trControl = trctrl,
          method = method
        )
    }
  return(fit)
}


# Plot model result including permutation

model_result <- function(prediction_both_models, permutation){
  
  
  HA_ROC <- 
    prediction_both_models %>% 
    map(~pROC::roc(.x$obs, .x$case, )) %>% 
    map(~{pROC:::get.coords.for.ggplot(.x) %>% as_tibble()}) %>% 
    bind_rows(.id = "dataset") %>% 
    mutate(Ionization = strsplit(dataset, "_") %>% map_chr(1),
           method = gsub("[.]", " ", strsplit(dataset, "_") %>% map_chr(2)))
  
  
  H0_ROC <- 
    permutation[[1]] %>% 
    get_best_glmnet_preds %>% 
    map(~arrange(.x, rowIndex)) %>% 
    map(~pROC::roc(.x$obs, .x$case, )) %>% 
    map_dfr(pROC:::get.coords.for.ggplot) %>% 
    as_tibble() %>% 
    group_by(`1-specificity`) %>% 
    summarise(
      sens_mean = mean(sensitivity),
      upper = sens_mean+1.96*sqrt((sens_mean*(1-sens_mean))/length(permutation[[1]])),
      lower = sens_mean-1.96*sqrt((sens_mean*(1-sens_mean))/length(permutation[[1]])),
      sens_lower = quantile(sensitivity, 0.05),
      sensitivity = quantile(sensitivity, 0.95))
  
  ggplot()+
    geom_line(data = HA_ROC, aes(x = `1-specificity`, y=sensitivity, color = Ionization))+
    geom_ribbon(data = H0_ROC, aes(x = `1-specificity`, ymin = lower, ymax=upper), alpha = 0.1)+
    geom_segment(data = H0_ROC, aes(x = 0, y = 0, xend = 1, yend = 1), color = "gray60", linetype = 2)+
    scale_color_brewer(palette = "Set1", direction = -1)+
    theme_minimal()+
    facet_wrap(~method)+
    theme(panel.grid = element_blank(),
          strip.text = element_text(size = 7),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
  
}



# Get the predictions of the hyperparameter set that yields the best AUC ROC. 
# Takes a list of caret models as input
# Returns a list of prediction data frames
get_best_glmnet_preds <- function(glmnet_caret_outputs) {
  glmnet_caret_outputs %>% 
    map(~{
      best_params <- 
        .x$pred %>% 
        group_split(alpha, lambda) %>%
        map_dfr(~tibble(alpha = .x$alpha[1], 
                        lambda = .x$lambda[1], 
                        ROC = twoClassSummary(as.data.frame(.x), lev = levels(.x$obs))["ROC"])) %>% 
        arrange(-ROC)
      .x$pred %>% 
        filter(alpha==best_params$alpha[1], lambda==best_params$lambda[1])
    }) %>% 
    map(~arrange(.x, rowIndex)) 
}


get_best_ranger_preds <- function(ranger_caret_outputs) {
  ranger_caret_outputs %>% 
    map(~{
      best_params <- 
        .x$pred %>% 
        group_split(mtry, splitrule) %>%
        map_dfr(~tibble(mtry = .x$mtry[1], 
                        splitrule = .x$splitrule[1], 
                        ROC = twoClassSummary(as.data.frame(.x), lev = levels(.x$obs))["ROC"])) %>% 
        arrange(-ROC)
      .x$pred %>% 
        filter(mtry==best_params$mtry[1], splitrule==best_params$splitrule[1])
    }) %>% 
    map(~arrange(.x, rowIndex)) 
}


