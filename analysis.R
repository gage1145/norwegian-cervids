library(tidyverse)
library(pROC)
library(lme4)
library(emmeans)
library(broom)
library(scales)

main_theme <- theme(
	plot.title = element_text(size=30),
	axis.title = element_text(size=24),
	axis.text = element_text(size=12),
	legend.title = element_text(size=24),
	legend.text = element_text(size=20),
	strip.text = element_text(size=24)
)

df_clean <- read.csv("data/clean_data.csv", check.names=FALSE, row.names = 1) %>%
  filter(Dilutions > -6) %>%
  filter(!(Tissue == "LN" & Dilutions > -2))


# Step 1: Animal-level aggregation ---------------------------------------


df_animal <- df_clean %>%
  mutate(
		Tissue = ifelse(str_detect(Tissue, "SK") | Tissue == "Ear", "skin", Tissue),
		cutoff = as.factor(cutoff),
    .by = c(species, Animal, Assay, Tissue, cutoff, Dilutions)
	)


# PCA --------------------------------------------------------------------


pca <- df_animal %>%
  select(MPR, MS, AUC) %>%
  prcomp(scale. = TRUE)
summary(pca)

df_long <- df_animal %>%
  mutate(
    pc1 = pca$x[,1], 1, 0, 
    pc2 = pca$x[,2], 1, 0, 
    pc3 = pca$x[,3], 1, 0
  )

# PCA Visualization
df_long %>%
  mutate(cutoff = paste(cutoff, "hr")) %>%
  ggplot(aes(pc1, pc2, color = ELISA)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.2) +
  scale_color_manual(values = c("darkblue", "maroon")) +
  facet_grid(cols = vars(cutoff), rows = vars(Tissue)) +
  # coord_fixed() +
  main_theme +
  theme(
    legend.position = "bottom"
  )
ggsave("pca.png", path = "figures", height=6, width = 12)

distinct_combos <- df_long %>%
  summarize(
    # value = mean(value),
    .by = c(species, Tissue, Assay, Dilutions, cutoff)
  )

get_roc <- function(species, Tissue, Assay, Dilutions, cutoff) {
  temp_df <- df_long %>%
    filter(
      species   == .env$species,
      Tissue    == .env$Tissue,
      Assay     == .env$Assay,
      Dilutions == .env$Dilutions,
      cutoff    == .env$cutoff
    )
  
  if (n_distinct(temp_df$ELISA) < 2) {
    return(NULL)
  }

  temp_roc <- roc(temp_df, ELISA, pc1, direction = ">")
  temp_roc$species   <- species
  temp_roc$Tissue    <- Tissue
  temp_roc$Assay     <- Assay
  temp_roc$Dilutions <- Dilutions
  temp_roc$cutoff    <- cutoff

  temp_roc
}

roc_list <- distinct_combos %>%
  pmap(get_roc) %>%
  compact()

get_roc_auc <- function(x) {
  data.frame(
    species   = x$species,
    Tissue    = x$Tissue,
    Assay     = x$Assay,
    Dilutions = x$Dilutions,
    cutoff    = x$cutoff,
    auc       = as.numeric(auc(x))
  )
}

get_roc_info <- function(x) {
  temp_coords <- coords(x) 
    # mutate(youden = sensitivity + specificity - 1)
  # optimal_threshold <- temp_coords$threshold[which(temp_coords$youden == max(temp_coords$youden))]
  tibble(
    species     = x$species,
    Tissue      = x$Tissue,
    Assay       = x$Assay,
    Dilutions   = x$Dilutions,
    cutoff      = x$cutoff,
    metric      = x$metric,
    threshold   = list(temp_coords$threshold),
    sensitivity = list(temp_coords$sensitivity),
    specificity = list(temp_coords$specificity),
    # youden      = list(temp_coords$youden),
    auc = as.numeric(auc(x)),
    # optimal_threshold = optimal_threshold
  )
}

roc_aucs <- map_dfr(roc_list, get_roc_auc) %>%
  mutate(
    auc = round(auc, 3),
    y = ifelse(Assay == "NQ", 0.25, 0.15)
  )

roc_coords <- map_dfr(roc_list, get_roc_info) %>%
  unnest(c(sensitivity, specificity, threshold)) %>%
  group_by(species, Tissue, Assay, Dilutions, cutoff) %>%
  arrange(
    sensitivity,
    .by_group = TRUE
  ) %>%
  ungroup()

# Find top combinations of variables.
top_combos <- roc_aucs %>% 
  group_by(species, Tissue, Assay, Dilutions) %>%
  slice_max(auc, with_ties = TRUE) %>%
  slice_min(as.numeric(cutoff))

make_roc_plot <- function(x) {
  x %>%
		ggplot(aes(specificity, sensitivity, color = Assay, group = Assay)) +
		geom_step(linewidth=1.2) +
		geom_abline(
      slope = 1, intercept = 1, linetype = "dashed", 
      color = "gray50", inherit.aes = FALSE
    ) +
		facet_grid(cols = vars(round(Dilutions, 3)), rows = vars(Tissue)) +
    scale_color_manual(values=c("#ff9c59", "purple")) +
    scale_x_reverse(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    coord_fixed() +
		labs(
			x = "Specificity", y = "Sensitivity"
		) +
		main_theme
}

# Moose
top_moose <- top_combos %>%
  filter(species == "M")

roc_coords %>%
  filter(species == "M", 
    Dilutions %in% top_moose$Dilutions, cutoff %in% top_moose$cutoff
  ) %>%
  mutate(cutoff = as.numeric(cutoff)) %>%
  filter(cutoff == min(cutoff)) %>%
	make_roc_plot() +
  geom_text(
    aes(x = 0.35, y = y, color = Assay, label = sprintf("%s, %shr AUC = %s", Assay, cutoff, format(round(auc, 3), nsmall=3)), group = Assay),
    data = top_moose, inherit.aes = FALSE, size = 5
  ) +
  ggtitle("Moose ROC Curves")

ggsave("figures/moose_roc_curves.png", width = 16, height = 10)


# Reindeer
top_reindeer <- top_combos %>%
  filter(species == "R", Tissue != "skin")

roc_coords %>%
  filter(species == "R", Dilutions %in% top_reindeer$Dilutions, cutoff %in% top_reindeer$cutoff) %>%
  mutate(cutoff = as.numeric(cutoff)) %>%
  filter(cutoff == min(cutoff)) %>%
	make_roc_plot() +
  geom_text(
    aes(x = 0.35, y = y, color = Assay, label = sprintf("%s, %shr AUC = %s", Assay, cutoff, format(round(auc, 3), nsmall=3)), group = Assay),
    data = top_reindeer, inherit.aes = FALSE
  ) +
  ggtitle("Reindeer ROC Curves")


ggsave("figures/reindeer_roc_curves.png", width = 12, height = 8)

# Red deer
top_reddeer <- top_combos %>%
  filter(species == "RD")

roc_coords %>%
  filter(species == "RD", Dilutions %in% top_reddeer$Dilutions, cutoff %in% top_reddeer$cutoff) %>%
  mutate(cutoff = as.numeric(cutoff)) %>%
  filter(cutoff == min(cutoff)) %>%
	make_roc_plot() +
  ggtitle("Red Deer ROC Curves")


ggsave("figures/reddeer_roc_curves.png", width = 12, height = 8)

best <- bind_rows(top_moose, top_reindeer, top_reddeer)
write.csv(best, "data/top_conditions.csv", row.names=FALSE)



# K-means ----------------------------------------------------------------


library(plotly)

df_km <- df_long %>%
  mutate(
    pc2 = pca$x[,2], 
    pc3 = pca$x[,3],
  ) %>%
  filter(if_all(c(MPR, MS, AUC), ~ . > 0)) %>%
  mutate(
    across(c(MPR, MS, AUC), log)
  )

set.seed(102470238)

km.out <- df_km %>%
  select(MPR, MS, AUC) %>%
  kmeans(centers=3, nstart = 20)
km.out

df_km$cluster_id <- factor(km.out$cluster)

df_km %>%
  plot_ly(x = ~MPR, y = ~MS, frame = ~mpi, color = ~ELISA, colors = c('#BF382A', '#0C4B8E'), marker = list(size = 3)) %>%
  add_markers(opacity = 0.1) %>%
  animation_slider()
