library(tidyverse)
library(pROC)
library(ggpubr)
library(ggrepel)
library(janitor)


# Setup ------------------------------------------------------------------

for_manuscript <- TRUE

main_theme <- theme(
  plot.title = element_text(size = 30, hjust=0.5),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 24),
  legend.text = element_text(size = 20),
  strip.text = element_text(size = 24)
)

quic_colors  <- c("Nano-QuIC" = "darkorange", "RT-QuIC" = "purple")
elisa_colors <- c("Negative" = "navy", "Positive" = "red")

# Expand the short codes used in raw filenames/sample IDs into the human-readable
# species/tissue/assay labels used throughout the analysis and figures.
format_labels <- function(df) {
  df <- clean_names(df, replace = c("TtT" = "ttt"))
  species_key <- c(M = "Moose", R = "Reindeer", RD = "Red Deer")
  tissue_key <- c(BR = "Brainstem", LN = "Lymph Node", SK = "Skin", Ear = "Skin")
  assay_key <- c(RT = "RT-QuIC", NQ = "Nano-QuIC")
  if ("species" %in% colnames(df)) df <- mutate(df, species = species_key[str_extract(species, "[A-Z]{1,2}")])
  if ("tissue" %in% colnames(df)) df <- mutate(df, tissue = tissue_key[str_extract(tissue, "[A-Z]{2}|Ear")])
  if ("assay" %in% colnames(df)) df <- mutate(df, assay = assay_key[str_extract(assay, "[A-Z]{2}")])
  mutate(df, across(where(is.character), trimws))
}


# Clean incoming data ----------------------------------------------------


df_clean <- read.csv("data/clean_data.csv", check.names = FALSE, row.names = 1) %>%
  format_labels() %>%
  # Drop very-dilute samples (<= 1e-6) and the higher-concentration lymph node
  # series (> 1e-2), which were only run for a subset of animals.
  filter(dilutions > -6 & !(tissue == "Lymph Node" & dilutions > -2)) %>%
  mutate(
    across(c(assay, animal, animal_id, species, tissue, rxn, id, genotype), as.factor),
    dilutions = round(dilutions, 2),
    # tissue = ifelse(tissue == "Brain", "Brainstem", tissue)
  )


# Metric Correlation -----------------------------------------------------


metric_combos <- as.data.frame(t(combn(c("mpr", "ms", "auc"), 2)))
colnames(metric_combos) <- c("x", "y")

make_cor_plot <- function(x, y) {
  df_clean %>%
    filter(cutoff == 48) %>%
    ggplot(aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(color = "purple", alpha = 0.2) +
    stat_smooth(method = "lm", color = "darkorange") +
    stat_cor(
      aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~','~")),
      label.y.npc = 1, size = 6
    ) +
    labs(
      x = toupper(x),
      y = toupper(y)
    ) +
    main_theme
}

cor_plots <- pmap(metric_combos, make_cor_plot)
ggarrange(plotlist = cor_plots, ncol = 3, align = "hv")
ggsave("corplot_separated.png", path = "figures", width = 16, height = 12)


# Principal Component Analysis -------------------------------------------


# Collapse the three correlated QuIC metrics (MPR, MS, AUC) into orthogonal
# components so downstream ROCs can be driven by a single positivity score.
pca <- df_clean %>%
  select(mpr, ms, auc) %>%
  prcomp(scale. = TRUE)
summary(pca)

df_clean <- mutate(
  df_clean,
  pc1 = pca$x[, 1],
  pc2 = pca$x[, 2],
  pc3 = pca$x[, 3]
)

write.csv(df_clean, "data/pca.csv")

# PCA Visualization
pca_plot <- df_clean %>%
  mutate(
    # Flip pc1 so positives load high, then shift to be strictly positive
    # before log-transform. The shift keeps the smallest value at log(1) = 0.
    pc1 = log(-pc1 + abs(min(-pc1)) + 1),
    elisa = ifelse(elisa, "Positive", "Negative"),
    tissue = as.character(tissue),
    cutoff = paste(cutoff, "hr")
  ) %>%
  ggplot(aes(dilutions, pc1, color = elisa)) +
  geom_point(alpha = 0.1, position = position_jitter(0.2)) +
  scale_color_manual(values = elisa_colors) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  facet_grid(cols = vars(cutoff), rows = vars(tissue)) +
  labs(
    x = "Log Dilution Factors",
    y = sprintf("Log(pc1 + %s)", round(abs(min(-df_clean$pc1)) + 1, 2)),
    color = "ELISA"
  ) +
  main_theme +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
pca_plot
ggsave("pca.png", path = "figures", height = 10, width = 12)


# ROC --------------------------------------------------------------------


## Compute ROC objects ---------------------------------------------------


distinct_combos <- distinct(df_clean, species, tissue, assay, dilutions, cutoff)

# Fit one ROC curve per (species, tissue, assay, dilution, cutoff) cell.
# direction = ">" tells pROC that lower pc1 = ELISA positive (since pc1 was
# built from raw QuIC metrics where positives sit at the negative end of pc1).
# Cells without both ELISA classes are skipped; meta-fields are stashed on
# the returned object so downstream code can recover the cell identity.
get_roc <- function(species, tissue, assay, dilutions, cutoff) {
  meta <- as.list(environment())

  temp_df <- df_clean %>%
    filter(
      species == .env$species,
      tissue == .env$tissue,
      assay == .env$assay,
      dilutions == .env$dilutions,
      cutoff == .env$cutoff
    )

  if (n_distinct(temp_df$elisa) < 2) {
    return(NULL)
  }

  temp_roc <- roc(temp_df, elisa, pc1, direction = ">")
  temp_roc[names(meta)] <- meta
  temp_roc
}

roc_list <- distinct_combos %>%
  pmap(get_roc) %>%
  compact()


## Build coords + AUC tables ---------------------------------------------


# coords() returns variable-length vectors per ROC, so they are stored as
# list-columns and unnested below into a long table of curve points.
get_roc_info <- function(x) {
  temp_coords <- coords(x)
  tibble(
    species     = x$species,
    tissue      = x$tissue,
    assay       = x$assay,
    dilutions   = x$dilutions,
    cutoff      = x$cutoff,
    threshold   = list(temp_coords$threshold),
    sensitivity = list(temp_coords$sensitivity),
    specificity = list(temp_coords$specificity),
    roc_auc     = as.numeric(auc(x))
  )
}

roc_coords <- map_dfr(roc_list, get_roc_info) %>%
  unnest(c(sensitivity, specificity, threshold)) %>%
  group_by(species, tissue, assay, dilutions, cutoff) %>%
  arrange(sensitivity, .by_group = TRUE) %>%
  ungroup()

roc_aucs <- roc_coords %>%
  distinct(species, tissue, assay, dilutions, cutoff, roc_auc) %>%
  mutate(roc_auc = round(roc_auc, 3))

write.csv(roc_aucs, "data/aucs.csv")

# For each (species, tissue, assay, dilution), pick the cutoff with the highest
# AUC; on ties, prefer the earliest cutoff (faster assay readout).
top_combos <- roc_aucs %>%
  group_by(species, tissue, assay, dilutions) %>%
  slice_max(roc_auc, with_ties = TRUE) %>%
  slice_min(cutoff, with_ties = FALSE) %>%
  ungroup()

# Summary table for Zoe
df_sum <- df_clean %>%
  # filter(assay == "Nano-QuIC" & dilutions == -2 | assay == "RT-QuIC" & dilutions == -3) %>%
  group_by(animal, species, tissue, assay, dilutions, cutoff, elisa) %>%
  right_join(top_combos) %>%
  filter(tissue != "Skin") %>%
  group_by(animal, species, tissue, assay, elisa) %>%
  summarize(
    perc_crossed = mean(crossed),
  ) %>%
  pivot_wider(names_from = assay, values_from = perc_crossed) %>%
  pivot_wider(names_from = tissue, values_from = c(`Nano-QuIC`, `RT-QuIC`))
write.csv(df_sum, "sum.csv")

## Plot ROC curves -------------------------------------------------------


make_roc_plot <- function(x, spec = NULL, tissues = NULL, fct_col = "dilutions", fct_row = "tissue") {
  filtered <- x
  top_filtered <- top_combos

  if (!is.null(spec)) {
    filtered <- filter(filtered, species %in% spec)
    top_filtered <- filter(top_filtered, species %in% spec)
  }
  if (!is.null(tissues)) {
    filtered <- filter(filtered, tissue %in% tissues)
    top_filtered <- filter(top_filtered, tissue %in% tissues)
  }

  # Stagger the per-assay AUC labels vertically so they don't overlap.
  top_filtered <- mutate(top_filtered, y = ifelse(assay == "Nano-QuIC", 0.25, 0.1))

  filtered %>%
    inner_join(top_filtered %>% select(-roc_auc, -y)) %>%
    ggplot(aes(specificity, sensitivity, color = assay, group = assay)) +
    geom_step(linewidth = 1.2) +
    geom_abline(
      slope = 1, intercept = 1, linetype = "dashed",
      color = "gray50", inherit.aes = FALSE
    ) +
    geom_label(
      aes(label = paste0(assay, " = ", round(roc_auc, 2)), y = y, color = assay),
      x = -0.5, data = top_filtered, inherit.aes = FALSE, hjust = 0.5, size = 6
    ) +
    facet_grid(cols = vars(.data[[fct_col]]), rows = vars(.data[[fct_row]])) +
    scale_color_manual(values = quic_colors) +
    scale_x_reverse(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~., name = "Log Dilution Factor", breaks = NULL, labels = NULL)
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = "Specificity", y = "Sensitivity"
    ) +
    main_theme +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
    )
}

# roc_plot_specs <- tibble(
#   spec    = list("Moose", "Reindeer", "Red Deer", NULL),
#   tissues = list(c("Brainstem", "Lymph Node"), c("Brainstem", "Lymph Node"), c("Brainstem", "Lymph Node"), c("Skin")),
#   fct_row = c("tissue", "tissue", "tissue", "species"),
#   file    = c("moose_roc.png", "reindeer_roc.png", "reddeer_roc.png", "skin_roc.png"),
#   height  = c(10, 10, 10, 10),
#   width   = c(16, 16, 16, 20)
# )

# pwalk(roc_plot_specs, function(spec, tissues, fct_row, file, height, width) {
#   p <- make_roc_plot(roc_coords, spec = spec, tissues = tissues, fct_row = fct_row)
#   ggsave(file, plot = p, path = "figures", height = height, width = width)
# })


# Cutoff Analysis --------------------------------------------------------


optimal_tissue_params <- top_combos %>%
  mutate(cutoff = as.numeric(as.character(cutoff))) %>%
  summarize(
    across(c(roc_auc, cutoff), mean),
    .by = c(tissue, assay, dilutions)
  ) %>%
  group_by(tissue, assay) %>%
  filter(roc_auc == max(roc_auc))
write.csv(optimal_tissue_params, "data/top_tissue_conditions.csv", row.names = FALSE)


df_opt_dil <- df_clean %>%
  filter(
    !(dilutions %in% c(-5, -4)),
    tissue != "Skin"
  )

cutoff_fig <- function(t) {
  df_opt_dil %>%
    filter(cutoff == 48, tissue == t) %>%
    mutate(elisa = ifelse(elisa, "Positive", "Negative")) %>%
    ggplot(aes(species, ttt, color = elisa)) +
    geom_point(position = position_jitter(0.2), size = 2, alpha = 0.7) +
    geom_hline(yintercept = 32, linetype = "dashed") +
    facet_grid(cols = vars(assay), rows = vars(dilutions)) +
    scale_y_continuous(breaks = seq(0, 48, 6)) +
    scale_color_manual(values = elisa_colors) +
    ggtitle(t) +
    coord_flip() +
    labs(
      title = t,
      y = "Time to Threshold (h)"
    ) +
    main_theme +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank()
    )
}

cutoff_figs <- lapply(unique(df_opt_dil$tissue), cutoff_fig)

ggarrange(
  plotlist = cutoff_figs, font.label = list(size = 30),
  nrow = 1, align = "v", legend = "none"
)
ggsave("cutoffs.png", path = "figures", width = 16, height = 10)


# AUC Line Graph -----------------------------------------------------------


mean_roc_aucs <- roc_aucs %>%
  summarize(
    roc_auc = mean(roc_auc),
    .by = c(tissue, assay, dilutions, cutoff)
  )

best_cutoffs <- mean_roc_aucs %>%
  filter(
    roc_auc == max(roc_auc),
    .by = c(tissue, assay, dilutions)
  ) %>%
  slice_min(cutoff, with_ties = FALSE, by = c(tissue, assay, dilutions))

mean_roc_aucs %>%
  ggplot(aes(cutoff, roc_auc, color = assay, fill = assay, label = round(roc_auc, 2))) +
  geom_area(position = "jitter", alpha = 0.25) +
  geom_line(linewidth = 1.2, alpha = 0.25) +
  geom_hline(
    aes(yintercept = roc_auc, color = assay),
    data = best_cutoffs, alpha = 0.5, linetype = "dashed", linewidth = 1.2, show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = cutoff, color = assay),
    data = best_cutoffs, alpha = 0.5, linetype = "dashed", linewidth = 1.2, show.legend = FALSE
  ) +
  geom_text_repel(
    aes(
      label = paste0(round(roc_auc, 2), ", ", cutoff, "hr"),
      x = 25, y = roc_auc + 0.03, hjust = 0, color = assay
    ),
    data = best_cutoffs, size = 6, fontface = 2, force = 10, direction = "y",
    max.iter = 10000, show.legend = FALSE, inherit.aes = FALSE, point.size = NA,
    min.segment.length = 0, bg.color = "#333", bg.r = 0.03
  ) +
  facet_grid(cols = vars(round(dilutions, 2)), rows = vars(tissue)) +
  scale_y_continuous(breaks = seq(0, 1.1, 0.1), expand = expansion(c(-0.35, 0.05))) +
  scale_x_continuous(breaks = seq(28, 44, 4), limits = c(24, 48), expand = expansion(0.02)) +
  scale_color_manual(values = quic_colors) +
  scale_fill_manual(values = quic_colors) +
  labs(
    x = "Cutoff (hr)",
    y = "Area Under ROC Curve",
    title = "ROC Results by tissue",
    alpha = NULL
  ) +
  main_theme +
  theme(
    legend.position = c(0.8, 0.8),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    legend.background = element_blank(),
    legend.key.spacing.y = unit(0.01, "npc"),
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(fill = guide_legend(byrow = TRUE))

ggsave("roc_results_by_tissue.png", path = "figures", width = 16, height = 12)


# Cutoff vs. Dilution ----------------------------------------------------


# Compute per-tile x-edges so geom_rect can draw heatmap cells of variable
# width along the unevenly spaced dilution axis. Each tile spans halfway to
# its neighbours; endpoints fall back to mirroring the one available side.
df_rect <- expand.grid(
  tissue = sort(unique(roc_aucs$tissue)),
  dilutions = sort(unique(roc_aucs$dilutions))
) %>%
  mutate(
    w_next = abs((dilutions - lag(dilutions)) / 2),
    w_prev = abs((dilutions - lead(dilutions)) / 2),
    dmin = dilutions - w_next,
    dmax = dilutions + w_prev,
    dmin = ifelse(is.na(dmin), dilutions - w_prev, dmin),
    dmax = ifelse(is.na(dmax), dilutions + w_next, dmax),
    .by = tissue
  )

dilution_factors <- round(sort(unique(roc_aucs$dilutions)), 2)

roc_aucs %>%
  summarize(
    roc_auc = mean(roc_auc),
    .by = c(dilutions, assay, tissue, cutoff)
  ) %>%
  left_join(df_rect) %>%
  mutate(
    auc_lab = ifelse(roc_auc >= 0.75, str_remove(as.character(round(roc_auc, 2)), "0"), NA),
    auc_lab_position = (dmin + dmax) / 2
  ) %>%
  ggplot(aes(dilutions, cutoff, z = roc_auc, fill = roc_auc, group = roc_auc)) +
  geom_rect(aes(xmin = dmin, xmax = dmax, ymin = cutoff - 2, ymax = cutoff + 2)) +
  geom_text(aes(label = auc_lab, x = auc_lab_position), color = "white", size = 6) +
  labs(
    x = "Log Dilution Factor",
    y = "Cutoff Time (hr)",
    fill = "AUC of ROC"
  ) +
  facet_grid(tissue ~ assay) +
  scale_fill_gradient2(
    low = "red", mid = "lightgreen", high = "navy", midpoint = mean(roc_aucs$roc_auc),
    breaks = seq(0.4, 1, 0.1), limits = c(0.38, 1.01)
  ) +
  scale_x_continuous(breaks = sort(unique(dilution_factors))) +
  scale_y_continuous(breaks = seq(24, 48, 4)) +
  {if (for_manuscript) ggtitle("Area Under ROC Curve Heatmap") else .} +
  labs(
    fill = ifelse(for_manuscript, "AUC", "AUC of ROC"),
    x = "Log10 Dilution Factor"
  ) +
  main_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = ifelse(for_manuscript, "right", "top"),
    legend.title = element_text(margin = margin(r = 20)),
    legend.key.width =  unit(ifelse(for_manuscript, 0.5, 1), "in"))
  

ggsave(
  sprintf("cutoff_vs_dilution%s.png", ifelse(for_manuscript, "_manuscript", "")), 
  path = "figures", width = 16, height = 16
)



# RAF Graphs -------------------------------------------------------------


make_boxplot <- function(data, y, t) {
  

  cutoffs <- roc_aucs %>%
  filter(tissue == t) %>%
  summarize(
    cutoff = cutoff[which.max(roc_auc)][1],
    .by=c(species, assay, dilutions)
  ) 

  pvals <- data %>%
  filter(tissue == t, dilutions > -5) %>% 
  inner_join(cutoffs) %>%
  summarize(
    p.value = wilcox.test(!!sym(y) ~ elisa)$p.value,
    .by = c(species, assay, dilutions)
  ) %>%
  mutate(
    label = paste("p =", signif(p.value, 2)),
    xmin  = ifelse(assay == "Nano-QuIC", 0.875, 1.125),
    xmax  = ifelse(assay == "Nano-QuIC", 1.875, 2.125),
    y     = ifelse(
      assay == "Nano-QuIC", 
      ifelse(y == "raf", 0.5, 5.5), 
      ifelse(y == "raf", 0.75, 8.25)
    ),
    # color = ifelse(p.value < 0.05, "darkred", "black"),
    species = str_remove_all(species, " "),
    # elisa = ifelse(elisa, "Positive", "Negative"),
    dilutions = paste0("10^", dilutions)
  ) %>%
    filter(p.value < 0.05)

  pvals2 <- data %>%
  filter(tissue == t, dilutions > -5) %>% 
  inner_join(cutoffs) %>%
  summarize(
    p.value = wilcox.test(!!sym(y) ~ assay)$p.value,
    # elisa = sum(elisa) != 0,
    .by = c(species, elisa, dilutions)
  ) %>%
  mutate(
    label = paste("p =", signif(p.value, 2)),
    y     = ifelse(y == "raf", 1, 10),
    # color = ifelse(p.value < 0.05, "darkred", "black"),
    species = str_remove_all(species, " "),
    dilutions = paste0("10^", dilutions),
    elisa = ifelse(elisa, "Positive", "Negative"),
  ) %>%
    filter(p.value < 0.05)
  pvals2

  data %>%
    mutate(raf = ifelse(crossed, raf, 0)) %>%
    right_join(cutoffs) %>%
    filter(tissue == t, dilutions > -5) %>% 
    mutate(
      pc1 = -pc1,
      species = str_remove_all(species, " "),
      elisa = ifelse(elisa, "Positive", "Negative"),
      dilutions = paste0("10^", dilutions)
    ) %>%
    ggplot(aes(elisa, .data[[y]], fill = assay)) +
    geom_boxplot() +
    new_scale_color() +
    geom_text(
      aes(x = elisa, y = y, label = label),
      data = pvals2, inherit.aes = FALSE, size = 5, fontface = "bold"
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = y * 1.1, label = label),
      data = pvals, inherit.aes = FALSE, size = 5, fontface = "bold"
    ) +
    geom_segment(
      aes(x = xmin, xend = xmax, y = y, yend = y),
      data = pvals, inherit.aes = FALSE, linewidth = 0.5
    ) +
    geom_segment(
      aes(x = xmin, xend = xmin, y = y * 0.9, yend = y),
      data = pvals, inherit.aes = FALSE, linewidth = 0.5
    ) +
    geom_segment(
      aes(x = xmax, xend = xmax, y = y * 0.9, yend = y),
      data = pvals, inherit.aes = FALSE, linewidth = 0.5
    ) +
    scale_color_identity() +
    facet_grid(rows = vars(dilutions), cols=vars(species), scale = "free_x", labeller=label_parsed) +
    scale_fill_manual(values = quic_colors) +
    ggtitle(str_wrap(sprintf("RT-QuIC vs Nano-QuIC %s %s", t, toupper(y)), width = 40)) +
    labs(
      y = ifelse(y == "raf", "Rate of amyloid formation (1/h)", "PC1")
    ) +
    main_theme +
    theme(
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}

## Brain RAF Boxplot ---------------------------------------------------------
braf <- df_clean %>%
make_boxplot("raf", "Brainstem")
# ggsave("brain_bp.v3.png", path="figures/raf", width = 12, height = 10)


## Lymph RAF Node Boxplot ----------------------------------------------------
lraf <- make_boxplot(df_clean, "raf", "Lymph Node")
# ggsave("lymph_bp.v3.png", path="figures/raf", width = 12, height = 10)

## Skin RAF Boxplot ----------------------------------------------------------
dils <- sort(unique(df_clean$dilutions))
skraf <- df_clean %>%
  filter(dilutions > -4) %>%
  make_boxplot("raf", "Skin")
# ggsave("skin_bp.v3.png", path="figures/raf", width = 12, height = 10)

## Brain PCA Boxplot ---------------------------------------------------------
bpca <- make_boxplot(df_clean, "pc1", "Brainstem")
# ggsave("brain_bp_pca.v3.png", path="figures/pca", width = 12, height = 10)


## Lymph PCA Node Boxplot ----------------------------------------------------
lpca <- make_boxplot(df_clean, "pc1", "Lymph Node")
# ggsave("lymph_bp_pca.v3.png", path="figures/pca", width = 12, height = 10)

## Skin PCA Boxplot ----------------------------------------------------------
dils <- sort(unique(df_clean$dilutions))
skpca <- df_clean %>%
  filter(dilutions > -4) %>%
  make_boxplot("pc1", "Skin")
# ggsave("skin_bp_pca.v3.png", path="figures/pca", width = 12, height = 10)



# Combined boxplots ------------------------------------------------------


ggarrange(braf, bpca, ncol=2, common.legend=TRUE, legend="bottom", labels=c("A", "B"), font.label=list(size=30))
ggsave("combined_brain_box.png", path="figures", width=20, height=12)

ggarrange(lraf, lpca, ncol=2, common.legend=TRUE, legend="bottom", labels=c("A", "B"), font.label=list(size=30))
ggsave("combined_lymph_box.png", path="figures", width=20, height=12)

ggarrange(skraf, skpca, ncol=2, common.legend=TRUE, legend="bottom", labels=c("A", "B"), font.label=list(size=30))
ggsave("combined_skin_box.png", path="figures", width=20, height=12)


# Tables --------------------------------------------------------------------


top_all <- top_combos %>%
  filter(roc_auc > 0.7) %>%
  arrange(species, tissue, assay, desc(roc_auc), cutoff)
write.csv(top_all, "data/top_all.csv", row.names = FALSE)
