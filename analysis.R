library(tidyverse)
library(pROC)
library(lme4)
library(emmeans)
library(broom)
library(scales)
library(ggpubr)
library(patchwork)

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
pca_plot <- df_long %>%
  # filter(Tissue != "skin") %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Tissue = ifelse(Tissue == "skin", "SK", Tissue),
    cutoff = paste(cutoff, "hr")
  ) %>%
  ggplot(aes(pc1, pc2, color = ELISA)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("darkblue", "maroon")) +
  facet_grid(cols = vars(cutoff), rows = vars(Tissue)) +
  # coord_fixed() +
  ggtitle("Principal Component Analysis") +
  main_theme +
  theme(
    legend.position = "bottom"
  )
pca_plot
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
    y = ifelse(Assay == "NQ", 0.25, 0.1)
  )

roc_coords <- map_dfr(roc_list, get_roc_info) %>%
  unnest(c(sensitivity, specificity, threshold)) %>%
  group_by(species, Tissue, Assay, Dilutions, cutoff) %>%
  arrange(
    sensitivity,
    .by_group = TRUE
  ) %>%
  ungroup() %>%
  mutate(Dilutions = paste0("10^", round(Dilutions, 3)))

# Find top combinations of variables.
top_combos <- roc_aucs %>%
  group_by(species, Tissue, Assay, Dilutions) %>%
  # mutate(cutoff = as.numeric(cutoff)) %>%
  slice_max(auc, with_ties = TRUE) %>%
  slice_min(cutoff, with_ties = FALSE) %>%
  ungroup()

make_roc_plot <- function(x, fct_col = "Dilutions", fct_row = "Tissue") {
  x %>%
    mutate(
      cutoff = as.numeric(cutoff),
      Tissue = ifelse(Tissue == "skin", "SK", Tissue),
      Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC")
    ) %>%
    # filter(cutoff == min(cutoff)) %>%
    # mutate(Dilutions = paste0("10^", round(Dilutions, 3))) %>%
		ggplot(aes(specificity, sensitivity, color = Assay, group = Assay)) +
		geom_step(linewidth=1.2) +
		geom_abline(
      slope = 1, intercept = 1, linetype = "dashed", 
      color = "gray50", inherit.aes = FALSE
    ) +
		facet_grid(cols = vars(.data[[fct_col]]), rows = vars(.data[[fct_row]]), labeller = label_parsed) +
    scale_color_manual(values=c("darkorange", "purple")) +
    scale_x_reverse(limits=c(0,1), breaks = seq(0, 1, 0.2), sec.axis = sec_axis(~ ., name="Dilution Factor", breaks = NULL, labels = NULL)) +
    scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.2)) +
    # coord_fixed() +
		labs(
			x = "Specificity", y = "Sensitivity"
		) +
		main_theme +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
    )
}


# Brain and LN plots -----------------------------------------------------


{
# Moose
top_moose_base <- top_combos %>%
  filter(species == "M", Tissue != "skin") %>%
  mutate(Dilutions = paste0("10^", round(Dilutions, 3)))

top_moose_join <- top_moose_base %>%
  select(species, Tissue, Assay, Dilutions, cutoff)

top_moose <- top_moose_base %>%
  mutate(Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC"))

m_roc_plot <- roc_coords %>%
  inner_join(top_moose_join) %>%
	make_roc_plot() +
  geom_label(
    aes(x = 0.5, y = y, color = Assay, label = sprintf("%shr AUC = %s", cutoff, format(round(auc, 2), nsmall=2)), group = Assay),
    data = top_moose, inherit.aes = FALSE, size = 7, show.legend = FALSE
  ) +
  ggtitle("Moose ROC Curves") +
  theme(
    # axis.title.x = element_blank()
  )
m_roc_plot

# ggsave("figures/moose_roc_curves.png", width = 16, height = 8)


# Reindeer
top_reindeer_base <- top_combos %>%
  filter(species == "R", Tissue != "skin") %>%
  mutate(Dilutions = paste0("10^", round(Dilutions, 3)))

top_reindeer_join <- top_reindeer_base %>%
  select(species, Tissue, Assay, Dilutions, cutoff)

top_reindeer <- top_reindeer_base %>%
  mutate(Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC"))

r_roc_plot <- roc_coords %>%
  inner_join(top_reindeer_join) %>%
	make_roc_plot() +
  geom_label(
    aes(x = 0.5, y = y, color = Assay, label = sprintf("%shr AUC = %s", cutoff, format(round(auc, 2), nsmall=2)), group = Assay),
    data = top_reindeer, inherit.aes = FALSE, size = 7, show.legend = FALSE
  ) +
  ggtitle("Reindeer ROC Curves") +
  theme(
    # axis.title = element_blank()
  )
# r_roc_plot

# ggsave("figures/reindeer_roc_curves.png", width = 10, height = 6)

# Red deer
top_reddeer_base <- top_combos %>%
  filter(species == "RD", Tissue != "skin") %>%
  mutate(Dilutions = paste0("10^", round(Dilutions, 3)))

top_reddeer_join <- top_reddeer_base %>%
  select(species, Tissue, Assay, Dilutions, cutoff)

top_reddeer <- top_reddeer_base %>%
  mutate(Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC"))

rd_roc_plot <- roc_coords %>%
  inner_join(top_reddeer_join) %>%
	make_roc_plot() +
  geom_label(
    aes(x = 0.5, y = y, color = Assay, label = sprintf("%shr AUC = %s", cutoff, format(round(auc, 2), nsmall=2)), group = Assay),
    data = top_reddeer, inherit.aes = FALSE, size = 7, show.legend = FALSE
  ) +
  ggtitle("Red Deer ROC Curves") +
  theme(
    # axis.title = element_blank()
  )
# rd_roc_plot

# ggsave("figures/reddeer_roc_curves.png", width = 12, height = 8)


# Skin Plots -------------------------------------------------------------


top_skin_base <- top_combos %>%
  filter(Tissue == "skin") %>%
  mutate(Dilutions = paste0("10^", round(Dilutions, 3)))

top_skin_join <- top_skin_base %>%
  select(species, Tissue, Assay, Dilutions, cutoff)

top_skin <- top_skin_base %>%
  mutate(
    species = ifelse(species == "M", "Moose", "Red_Deer"),
    Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC")
  )

roc_plot_skin <- roc_coords %>%
  inner_join(top_skin_join) %>%
  mutate(species = ifelse(species == "M", "Moose", "Red_Deer")) %>%
	make_roc_plot(fct_row="species") +
  geom_label(
    aes(x = 0.5, y = y, color = Assay, label = sprintf("%shr AUC = %s", cutoff, format(round(auc, 2), nsmall=2)), group = Assay),
    data = top_skin, inherit.aes = FALSE, size = 7, show.legend = FALSE
  ) +
  ggtitle("Skin ROC Curves")
roc_plot_skin

# ggsave("figures/reindeer_roc_curves.png", width = 10, height = 6)






roc_joined <- ggarrange(m_roc_plot, r_roc_plot, rd_roc_plot, font.label = list(size=30), labels=c("A", "B", "C"), nrow=1, align = "v", legend = "none")
# roc_joined
ggarrange(roc_joined, roc_plot_skin, ncol=1, labels = c("", "D"), font.label = list(size=30), common.legend = TRUE, legend="bottom") +
  theme(
    plot.background = element_rect(fill="white")
  )
}

ggsave("combined.png", path="figures", height=20, width=35)


# Cutoff Analysis --------------------------------------------------------


optimal_tissue_params <- top_combos %>%
  mutate(cutoff = as.numeric(as.character(cutoff))) %>%
  summarize(
    across(c(auc, cutoff), mean),
    .by = c(Tissue, Assay, Dilutions)
  ) %>%
  group_by(Tissue, Assay) %>%
  filter(auc == max(auc))
write.csv(optimal_tissue_params, "data/top_tissue_conditions.csv", row.names = FALSE)


{
  cutoff_fig <- function(t) {
    df_long %>%
      filter(cutoff == 48, Tissue == t) %>%
      ggplot(aes(species, TtT, color=ELISA)) +
      geom_point(position=position_jitter(0.2), size=2, alpha=0.7) + 
      geom_hline(yintercept = 24, linetype = "dashed") +
      # stat_compare_means(
      #   method = "anova",
      #   label = "p.signif",
      #   label.y = 1,
      #   tip.length = 0.01,
      #   size=8,
      #   fontface = "bold",
      #   show.legend = FALSE,
      #   hide.ns = TRUE,
      # ) +
      facet_grid(cols=vars(Assay), rows=vars(Dilutions)) +
      scale_y_continuous(breaks=seq(0, 48, 12)) +
      scale_color_manual(values = c("darkblue", "maroon")) +
      ggtitle(t) +
      coord_flip() +
      # scale_y_log10()+
      main_theme +
      theme(
        axis.title.y = element_blank()
      )
  }
    
  cutoff_figs <- lapply(unique(df_long$Tissue), cutoff_fig)

  cutoffs_arranged <- ggarrange(plotlist = cutoff_figs, font.label = list(size=30), labels=c("A", "B", "C"), nrow=1, align = "v", legend = "none")
  ggarrange(cutoffs_arranged, pca_plot, ncol=1, labels = c("", "D"), font.label = list(size=30), common.legend = TRUE, legend="bottom")
}
