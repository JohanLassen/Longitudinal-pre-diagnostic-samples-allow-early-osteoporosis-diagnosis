
library(tidyverse)
library(xcms)
source("/data/projects/osteoporosis_metabolomics/johan/scripts/markdown/utils.R")
source("/data/projects/osteoporosis_metabolomics/johan/scripts/markdown/preprocessing_utils.R")
library(lubridate)

all_text_same_size <- function(s = 10){
  theme(plot.title = element_text(face = "plain", size = s), axis.text = element_text(face = "plain", size = s), axis.title = element_text(face = "plain", size = s), legend.title = element_text(face = "plain", size = s), legend.text = element_text(face = "plain", size = s), plot.subtitle = element_text(face = "plain", size = s))
}

filename = "msobject_neg.rds" 
# filename = "msobject.rds"

ms <- readr::read_rds(paste0("/data/projects/osteoporosis_metabolomics/johan/data/", filename))
ms$rowinfo$storage_time <- ms$rowinfo$age - ms$rowinfo$age_at_donation
ms$rowinfo$smoker <- ms$rowinfo$smoker>0
ms$rowinfo$batch <- as.factor(ms$rowinfo$batch) 

  
## Visualize transformation effect
(feature_distribution_log2 <-
    ms$values %>%
    pivot_longer(cols = everything()) %>%
    ggplot(aes(x=log2(value+1)))+
    geom_histogram(color = "white")+
    theme_classic()+
    labs(x="log2(intensity+1)"))

(feature_distribution <-
    ms$values %>%
    pivot_longer(cols = everything()) %>%
    ggplot(aes(x=value))+
    geom_histogram(color = "white", breaks = seq(0, 20000, length.out = 30))+
    coord_cartesian(xlim = c(0,20000))+
    theme_classic()+
    labs(x="intensity")
)

library(patchwork)

transformation_plot <- feature_distribution+feature_distribution_log2

ggsave(transformation_plot, 
       filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS1_transformed_features_", filename, ".pdf"), 
       width = 183, height = 100, unit = "mm")


## Impute NAs with 0 and log transform
ms$values0  <- ms$values
a <- ms$values %>% map(is.na) %>% map_dbl(mean)
ms$values0  <- ms$values0[,a<0.1]
ms$values0[is.na(ms$values0)] <- 0
ms$values0  <- log2(ms$values0+1) %>% as_tibble()
ms$rowinfo0 <- ms$rowinfo

tmp1 <- ms$rowinfo0
tmp1 <- tmp1 %>% filter(grepl("replicate", type2)) %>% mutate(rep = gsub(".*_", "", sample_id)) %>% select(rowid, sample_origin, rep)
tmp2 <- ms$values0[tmp1$rowid,]
df <- 
  tmp1 %>% 
  bind_cols(tmp2) %>% 
  pivot_longer(cols = -c(sample_origin, rowid, rep)) %>% 
  select(-rowid) %>% 
  pivot_wider(names_from = "rep", values_from = "value") %>% 
  group_by(sample_origin) %>% 
  mutate(peak_number = row_number())

# Plot PCA

rm_blinds <- ms$rowinfo0$type1 != "blind"

(ms0_pca <- 
    plot_pca(
      feature_matrix = ms$values0, 
      metadata = ms$rowinfo, 
      color_labels = c("batch", "storage_time", "region"), 
      palette = c("Set1", "YlOrRd", "Set1"))
)

ggsave(plot = ms0_pca, 
       filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS2_transformed_pca_", filename, ".pdf"),
       bg = "white",width = 183, height = 150, units = "mm")


tmp2 <- ms$rowinfo0 %>% filter(grepl("sample_rep", type2))
tmp1 <- ms$values0[tmp2$rowid,]
correlations <- 
  tmp2 %>% 
  select(sample_origin, sample_id) %>% 
  bind_cols(tmp1) %>% 
  mutate(sample_id = gsub(".*_", "", sample_id)) %>% 
  pivot_longer(cols = starts_with("M")) %>% 
  pivot_wider(names_from = sample_id, values_from = value) %>% 
  group_by(name) %>% 
  summarise(correlation = cor(rep1, rep2))

all_replicate_correlations <-
  ggplot(correlations, aes(x=correlation))+
  geom_histogram(color = "white", breaks = seq(-0.5, 1, length.out = 100))+
  theme_bw()+
  coord_cartesian(xlim=c(0.5, 1))

ggsave(
  filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS3_replicate_cor_before_norm_", filename, ".pdf"), 
  plot = all_replicate_correlations, 
  bg = "white", 
  width = 89, height = 40, units = "mm")

stable_features <- 
  correlations %>% 
  filter(correlation > 0.95) %>% 
  pull(name)


(cor_plots <- 
    tmp2 %>% 
    select(sample_origin, sample_id, batch) %>% 
    bind_cols(tmp1 %>% select(all_of(stable_features))) %>% 
    mutate(sample_id = gsub(".*_", "", sample_id)) %>% 
    pivot_longer(cols = starts_with("M")) %>% 
    pivot_wider(names_from = sample_id, values_from = value) %>% 
    filter(name %in% sample(stable_features, size=12)) %>% 
    ggplot(aes(x=rep1, y=rep2))+
    geom_point(aes(fill = batch), shape=21, color="white")+
    facet_wrap(~name, ncol=4, scales="free")+
    theme_minimal()+
    scale_fill_brewer(palette = "Set1")+
    labs(x="Replicate 1", y="Replicate 2"))

ggsave(
  filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS4_sampled_correlations_", filename, ".pdf"), 
  plot = cor_plots, 
  bg = "white",width = 183, height = 150, units = "mm")


## Row normalization 
target_info   <- ms$rowinfo0
target_values <- ms$values0 %>% tibble::as_tibble()

tmp1 <- target_info %>% filter(type1 == "sample") %>% select(rowid)
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
  dplyr::slice_min(order_by = range, prop = 0.8) %>% 
  pull(name)

tmp    <- rowSums(
  ms$values0[stable_features])
raw    <- max(tmp)* (ms$values0 / tmp)

ms$values1  <- as_tibble(raw)
ms$rowinfo1 <- target_info

# Plot rownorm result
rm_blinds_and_batch = !(ms$rowinfo1$type1 %in% c("QC", "blind"))
(pca_rownorm_noqc <- 
    plot_pca(
      feature_matrix = ms$values1[rm_blinds_and_batch,], 
      metadata = ms$rowinfo1[rm_blinds_and_batch,], 
      color_label = c("batch", "storage_time", "region"),
      palette = c("Set1", "YlOrRd", "Set1"),
      panels = 3
    )
)

ggsave(plot = pca_rownorm_noqc, 
       filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS5_rownorm_pca_", filename, ".pdf"), 
       bg = "white",width = 183, height = 150, units = "mm")


# Defining data and confounding variables
ms$i  <-ms$rowinfo1 %>% filter(type1 == "sample")
ms$v <- ms$values1[ms$i$rowid,]
ms$i$rowid <- 1:nrow(ms$i)

injection_order <- ms$i$injection_order 
batch           <- factor(ms$i$batch)
storage_time    <- ms$i$storage_time
region          <- factor(ms$i$region)


# Normalization - simultaneous removal of confounder effects
ms$v1 <- 
  ms$v %>% 
  map_dfc(~(.x-predict(lm(.x~region+batch*injection_order+batch*storage_time),
                       tibble(batch, injection_order, storage_time, region))))

(pca_storage_scaled <- 
    plot_pca(
      feature_matrix = ms$v1, 
      metadata = ms$i, 
      color_labels = c("batch", "storage_time", "region"),
      palette = c("Set1", "YlOrRd", "Set1"),
      panels = 3
    )
)


(pca_status <- 
    plot_pca(
      feature_matrix = ms$v1, 
      metadata = ms$i, 
      color_labels = c("status"),
      palette = c("Set1"),
      panels = 4
    )
)


ggsave(plot = pca_storage_scaled, 
       filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS6_pca_corrected_", filename, ".pdf"), 
       bg = "white",width = 183, height = 150, units = "mm")


## Remove features if QC relative standard deviation is more than 30 % in the non-batch normalized data
tmp1 <- ms$values1
tmp2 <- ms$rowinfo1

# Features that only exist in biological samples
QC_CV <- 
  tmp1 %>% 
  mutate(QC = ifelse(grepl("QC", tmp2$type1), "QC", "sample")) %>% 
  pivot_longer(cols = -QC) %>% 
  group_by(QC, name) %>% 
  summarise(
    RSD = sd(value)/mean(value)
  ) %>% 
  pivot_wider(names_from = QC, values_from = RSD)


# Correlations on normalized data

tmp2 <- ms$i %>% mutate(rowid=row_number()) %>% filter(grepl("sample_rep", type2))
tmp1 <- ms$v1[tmp2$rowid,]

correlations <- 
  tmp2 %>% 
  select(sample_origin, sample_id) %>% 
  bind_cols(tmp1) %>% 
  mutate(sample_id = gsub(".*_", "", sample_id)) %>% 
  pivot_longer(cols = starts_with("M")) %>% 
  pivot_wider(names_from = sample_id, values_from = value) %>% 
  group_by(name) %>% 
  summarise(correlation = cor(rep1, rep2))

all_replicate_correlations_after_norm <- 
  ggplot(correlations, aes(x = correlation))+
  geom_histogram(color ="white", show.legend = FALSE)+
  theme_bw()#+scale_fill_brewer(palette = "Set1")

ggsave(
  filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS3b_replicate_cor_after_norm_", filename, ".pdf"), 
  plot = all_replicate_correlations, 
  bg = "white", 
  width = 89, height = 40, units = "mm")


(reproducibility_plot <- 
    QC_CV %>% 
    ungroup() %>% 
    left_join(correlations %>% ungroup()) %>% 
    ggplot(aes(x=correlation, y=QC))+
    geom_hex(bins = 70)+
    scale_fill_viridis_c(trans = "log10")+
    theme_minimal()+
    labs(x="Replicate correlation (Pearson)", y="QC (coefficient of variation)", title = "Reproducibility")+
    geom_vline(xintercept = 0.8, color="red")+
    geom_hline(yintercept = 0.1, color = "red")+
    all_text_same_size()
)

ggsave(plot=reproducibility_plot, 
       filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS7_reproducibility_plot_", filename, ".pdf"),
       bg = "white",width = 183, height = 100, units = "mm")


# Remove bad features

### QC and replicate filtering
stable_features <-
  QC_CV %>% 
  ungroup() %>% 
  left_join(correlations %>% ungroup()) %>% 
  filter(correlation > 0.80 & QC < 0.1) %>% 
  pull(name)

ms$v2 <- 
  ms$v1 %>% 
  select(all_of(stable_features))


# Plot final dataset af feature filtering
(
  pca_preprocessed_storage <- 
    plot_pca(
      feature_matrix = ms$v2, 
      metadata = ms$i, 
      color_labels = c("batch", "storage_time", "region"),
      palette = c("Set1", "YlOrRd", "Set1"),
      panels = 3
    )
)


ggsave(plot = pca_preprocessed_storage, 
       filename = paste0("/data/projects/osteoporosis_metabolomics/johan/figures/figureS8_pca_filtered_final_", filename, ".pdf"), 
       bg = "white",width = 183, height = 150, units = "mm")


if (filename != "msobject_neg.rds"){
# ANnotation dataset
annotations <- read_csv("/data/projects/osteoporosis_metabolomics/johan/annotations.csv")
ms$annotated <- ms$v2 %>% select(any_of(annotations$name))
ms$values %>% select(any_of(annotations$name))
colnames(ms$annotated) <- annotations %>% filter(name %in% colnames(ms$annotated)) %>% pull(annotation)

not_annotated <- ms$v2 %>% select(-any_of(annotations$name))
annotated <- ms$v2 %>% select(any_of(annotations$name))
ms$v2 <- bind_cols(annotated, not_annotated)

# Plot final annotation
(
  pca_preprocessed_batch_annotated <- 
    plot_pca(
      feature_matrix = ms$annotated, 
      metadata = ms$i, 
      color_labels = c("batch", "storage_time", "status"),
      palette = c("Set1", "Spectral", "Set1"),
      panels = 3
    )
)
}

write_rds(ms, file = paste0("ml_preprocessed_", filename))

