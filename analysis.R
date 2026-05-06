library(tidyverse)
library(pROC)
library(lme4)
library(emmeans)
library(broom)
library(scales)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(skimr)

main_theme <- theme(
	plot.title = element_text(size=30),
	axis.title = element_text(size=24),
	axis.text = element_text(size=20),
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


# mpms <- cor(df_animal$mpr, df_animal$ms)

metric_combos <- combn(c("MPR", "MS", "AUC"), 2) |> t() |> as.data.frame()
colnames(metric_combos) <- c("x", "y")

make_cor_plot <- function(x, y) {
  df_cor <- df_animal %>%
    filter(cutoff == 48) %>%
    mutate(
    Tissue = ifelse(Tissue == "skin", "SK", Tissue),
    cutoff = paste(cutoff, "hr")
    )
  df_cor %>%
    ggplot(aes(x=.data[[x]], y=.data[[y]])) +
    geom_point(color="purple", alpha=0.2) +
    stat_smooth(method = "lm", color="darkorange") +
    # stat_regline_equation(aes(label = ..adj.rr.label..)) +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~','~")), label.y.npc = 1, size=6) +
    # facet_grid(rows=vars(Tissue), cols=vars(cutoff), scales = "fixed") +
    scale_color_gradient(low="navy", high="darkorange") +
    main_theme +
    theme(
      legend.position = "bottom",
      legend.text=element_text(size=12)
    )
}

cor_plots <- pmap(metric_combos, make_cor_plot)
ggarrange(plotlist = cor_plots, ncol=3, align="hv")
ggsave("corplot_separated.png", path="figures", width=16, height=12)
  

cor_plot <- df_animal %>%
  filter(cutoff == 48) %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Tissue = ifelse(Tissue == "skin", "SK", Tissue),
    cutoff = paste(cutoff, "hr")
  ) %>%
  ggplot(aes(x = MPR, y = MS, color = AUC)) +
  geom_point() +
  facet_grid(rows=vars(Tissue), cols=vars(cutoff), scales = "fixed") +
  scale_color_gradient(low="navy", high="darkorange") +
  main_theme +
  labs(
    title = "Metric Correlation",
    x = "Maxpoint Ratio",
    y = "Max Slope"
  ) +
  theme(
    legend.position = "bottom",
    legend.text=element_text(size=12)
  )
cor_plot
ggsave("corplot.png", path="figures", width = 12, height=12)

pca <- df_animal %>%
  select(MPR, MS, AUC) %>%
  prcomp(scale. = TRUE)
summary(pca)

df_long <- df_animal %>%
  mutate(
    pc1 = log(rescale(pca$x[,1], c(0, 1))), 
    pc2 = log(rescale(pca$x[,2], c(0, 1))), 
    pc3 = log(rescale(pca$x[,3], c(0, 1)))
  )

df_long <- df_animal %>%
  mutate(
    pc1 = pca$x[,1], 
    pc2 = pca$x[,2], 
    pc3 = pca$x[,3]
  )

write.csv(df_long, "data/pca.csv")

# PCA Visualization
pca_plot <- df_long %>%
  # filter(Tissue != "skin") %>%
  mutate(
    # pc1 = log(rescale(-pc1, c(0.01, 1))),
    pc1 = log(-pc1 + abs(min(-pc1)) + 1),
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Tissue = ifelse(Tissue == "skin", "SK", Tissue),
    cutoff = paste(cutoff, "hr")
  ) %>%
  ggplot(aes(Dilutions, pc1, color = ELISA)) +
  # geom_hline(yintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.1, position = position_jitter(0.2)) +
  scale_color_manual(values = c("navy", "red")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  # scale_y_log10() +
  facet_grid(cols = vars(cutoff), rows = vars(Tissue)) +
  # coord_fixed() +
  labs(
    x = "Log Dilution Factors",
    y = sprintf("Log(pc1 + %s)", round(abs(min(-df_long$pc1)) + 1, 2))
  ) +
  # ggtitle("Principal Component Analysis") +
  main_theme +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "top",
    legend.text=element_text(size=20),
    axis.text = element_text(size = 16)
  )
pca_plot
ggsave("pca.png", path = "figures", height=10, width = 12)

ggarrange(cor_plot, pca_plot, ncol=2, align="h", widths=c(2, 5))
ggsave("cor_pca_combo.png", path="figures", width=16, height=12)


# ROC --------------------------------------------------------------------


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

write.csv(roc_aucs, "data/aucs.csv")

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
  mutate(Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC")) %>%
  arrange(Tissue, desc(auc), cutoff)

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
  mutate(Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC")) %>%
  arrange(Tissue, desc(auc), cutoff)

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
  mutate(Assay = ifelse(Assay == "RT", "RT-QuIC", "Nano-QuIC")) %>%
  arrange(Tissue, desc(auc), cutoff)

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


df_opt_dil <- df_long %>%
  filter(
      !(Dilutions %in% c(-5, -4)),
      Tissue != "skin"
  ) %>%
  mutate(
      Tissue = factor(Tissue, levels = c("BR", "LN"), labels = c("Brainstem", "Lymph Node")),
      Assay = factor(Assay, levels = c("NQ", "RT"), labels=c("Nano-QuIC", "RT-QuIC")),
      # Dilutions = factor(Dilutions, levels = c(-3, -2), labels = c(expression("10^{-3}"), expression("10^{-2}")))
  )

cutoff_fig <- function(t) {
  df_opt_dil %>%
    filter(cutoff == 48, Tissue == t) %>%
    ggplot(aes(species, TtT, color=ELISA)) +
    geom_point(position=position_jitter(0.2), size=2, alpha=0.7) + 
    geom_hline(yintercept = 32, linetype = "dashed") +
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
    scale_y_continuous(breaks=seq(0, 48, 6)) +
    scale_color_manual(values = c("darkblue", "maroon")) +
    ggtitle(t) +
    coord_flip() +
    labs(
      title = t,
      y = "Time to Threshold (h)"
    ) +
    main_theme +
    theme(
      plot.title = element_text(hjust=0.5),
      axis.title.y = element_blank()
    )
} 
  cutoff_figs <- lapply(unique(df_opt_dil$Tissue), cutoff_fig)

  cutoffs_arranged <- ggarrange(
    plotlist = cutoff_figs, font.label = list(size=30), 
    nrow=1, align = "v", legend = "none"
  )
  cutoffs_arranged
  ggsave("cutoffs.png", path="figures", width=16, height=10)


# AUC Line Graph -----------------------------------------------------------



mean_roc_aucs <- roc_aucs %>%
  mutate(
    cutoff = as.numeric(as.character(cutoff)),
    Assay = factor(Assay, levels = c("NQ", "RT"), labels = c("Nano-QuIC", "RT-QuIC"))
  ) %>%
  summarize(
    auc = mean(auc),
    .by = c(Tissue, Assay, Dilutions, cutoff)
  )

best_cutoffs <- mean_roc_aucs %>%
  filter(
    auc == max(auc),
    .by = c(Tissue, Assay, Dilutions)
  ) %>%
  slice_min(cutoff, with_ties = FALSE, by = c(Tissue, Assay, Dilutions)) %>%
  mutate(
    hjust = ifelse(cutoff > 36, 1, 0),
    offset = ifelse(cutoff > 36, cutoff - 1, cutoff + 1)
  )

mean_roc_aucs %>%
  ggplot(aes(cutoff, auc, color=Assay, fill=Assay, label=round(auc, 2))) +
  # geom_point(size=4) +
  geom_area(position="jitter", alpha = 0.25) +
  geom_line(linewidth=1.2, alpha = 0.25) +
  geom_hline(
    aes(yintercept = auc, color=Assay), 
    data=best_cutoffs, alpha = 0.5, linetype="dashed", linewidth=1.2, show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = cutoff, color = Assay), 
    data=best_cutoffs, alpha = 0.5, linetype="dashed", linewidth=1.2, show.legend = FALSE
  ) +
  geom_text_repel(
    aes(
      label=paste0(round(auc, 2), ", ", cutoff, "hr"), 
      x = 25, y = auc + 0.03, hjust = 0, color = Assay
    ), 
    data = best_cutoffs, size = 6, fontface = 2, force = 10, direction = "y",
    max.iter = 10000, show.legend = FALSE, inherit.aes = FALSE, point.size = NA,
    min.segment.length = 0, bg.color = "#333", bg.r = 0.03
  ) +
  facet_grid(cols=vars(round(Dilutions, 2)), rows=vars(Tissue)) +
  scale_y_continuous(breaks=seq(0, 1.1, 0.1), expand = expansion(c(-0.35, 0.05))) +
  scale_x_continuous(breaks=seq(28,44,4), limits=c(24,48), expand = expansion(0.02)) +
  scale_color_manual(values=c("darkorange", "purple")) +
  scale_fill_manual(values=c("darkorange", "purple")) +
  labs(
    x = "Cutoff (hr)",
    y = "Area Under ROC Curve",
    title = "ROC Results by Tissue",
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
    # legend.spacing.y = unit(10, "npc"),
    plot.title = element_text(hjust=0.5)
  ) +
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))

ggsave("roc_results_by_tissue.png", path="figures", width=16, height=12)



# Cutoff vs. Dilution ----------------------------------------------------


top_conditions <- read.csv("data/top_conditions.csv") %>%
  slice_max(Dilutions, by = c(species, Tissue, Assay))

all_aucs <- read.csv("data/aucs.csv") 


df_pca <- read.csv("data/pca.csv") %>%
  mutate(
    Tissue = ifelse(str_detect(Tissue, "SK"), "skin", Tissue),
    RAF = ifelse(crossed, RAF, 0)
  ) 


df_top <- df_pca %>%
  right_join(top_conditions) 

df_all_aucs <- df_pca %>%
  right_join(all_aucs)

# All
df_top %>%
  # filter(Tissue == "BR" & Dilutions > -6) %>%
  ggplot(aes(ELISA, RAF, fill = Assay, color=Assay)) +
  # geom_point(alpha=0.5, size=1.2, position=position_jitter(0.1))  +
  geom_boxplot(aes(outlier.color = ELISA), color="black") +
  # stat_compare_means(label.y.npc = 0.9)+
  # stat_compare_means(aes(label = after_stat(p.signif)),
  #                   method = "wilcox.test", label.y = 5, size = 7) +
  # facet_grid(rows=vars(Tissue)) +

  facet_grid(rows=vars(Tissue), cols=(vars(species))) +
  geom_label(aes(y = 0.3, label=paste0("Cutoff = ", cutoff, "\nDilution = ", round(Dilutions, 2), "\nAUC = ", round(auc, 3))), fill="white", color="black", inherit.aes = TRUE) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(
    title = "Comparisons of RAF at Optimal Conditions"
  )

df_rect <- expand.grid(
  Tissue = sort(unique(all_aucs$Tissue)),
  Dilutions = sort(unique(all_aucs$Dilutions))
) %>%
  mutate(
    w_next = abs((Dilutions - lag(Dilutions)) / 2),
    w_prev = abs((Dilutions - lead(Dilutions)) / 2),
    dmin = Dilutions - w_next,
    dmax = Dilutions + w_prev,
    dmin = ifelse(is.na(dmin), Dilutions - w_prev, dmin),
    dmax = ifelse(is.na(dmax), Dilutions + w_next, dmax),
    .by = Tissue
  ) 

dilution_factors <- round(sort(unique(all_aucs$Dilutions)), 2)

all_aucs %>%
  summarize(
    auc = mean(auc),
    .by = c(Dilutions, Assay, Tissue, cutoff)
  ) %>%
  left_join(df_rect) %>%
  mutate(
    auc_lab = ifelse(auc >= 0.75, str_remove(as.character(round(auc, 2)), "0"), NA),
    auc_lab_position = (dmin + dmax) / 2
  ) %>%
  ggplot(aes(Dilutions, cutoff, z=auc, fill=auc, group=auc)) +
  geom_rect(aes(xmin=dmin, xmax = dmax, ymin = cutoff - 2, ymax = cutoff + 2)) +
  geom_text(aes(label = auc_lab, x = auc_lab_position), color="white", size = 6) +
  labs(
    x = "Log Dilution Factor",
    y = "Cutoff Time (hr)",
    fill = "AUC of ROC"
  ) +
  facet_grid(Tissue ~ Assay) +
  scale_fill_gradient2(
    low="red", mid="lightgreen", high="navy", midpoint=mean(all_aucs$auc), 
    breaks=seq(0.4, 1, 0.1), limits = c(0.38, 1.01)
  ) +
  scale_x_continuous(breaks = sort(unique(dilution_factors))) +
  scale_y_continuous(breaks = seq(24, 48, 4)) +
  main_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top",
    legend.title = element_text(margin = margin(r = 20)),
    legend.key.width = unit(1, "in")
  )

ggsave("cutoff_vs_dilution.png", path="figures", width=16, height=16)


# Tables --------------------------------------------------------------------


top_all <- bind_rows(top_moose, top_reindeer, top_reddeer) %>%
  select(-y) %>%
  filter(auc > 0.7) %>%
  arrange(species, Tissue, Assay, desc(auc), cutoff)
write.csv(top_all, "data/top_all.csv", row.names=FALSE)
