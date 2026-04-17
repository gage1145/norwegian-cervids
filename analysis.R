library(tidyverse)
library(pROC)
library(lme4)
library(emmeans)
library(broom.mixed)
library(plotROC)

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
  group_by(species, Animal, Assay, Tissue, cutoff, Dilutions) %>%
  mutate(
		Tissue = ifelse(str_detect(Tissue, "SK") | Tissue == "Ear", "skin", Tissue),
		cutoff = as.factor(cutoff)
	)
  # summarise(
  #   across(c(MPR, MS, AUC, TtT, RAF, ELISA), \(x) mean(x, na.rm = TRUE)),
  #   prop_crossed = mean(crossed, na.rm = TRUE),
  #   n_wells = n(),
  #   .groups = "drop"
  # ) 
  # Filter for false ROC signals
  # filter(
    # !(prop_crossed == 0 & ELISA), 
      # (((ELISA == 1 & min(MPR)) > (ELISA == 0 & max(MPR))) | ((ELISA == 0 & min(MPR)) > (ELISA == 1 & max(MPR))))),
    # .by = c(species, Assay, Tissue, cutoff, Dilutions)
  # )


# Step 2: Improved ROC Analysis ------------------------------------------

df_roc_long <- df_animal %>%
  # select(-c(TtT, RAF, prop_crossed, n_wells)) %>%
  mutate(ELISA = as.integer(ELISA)) %>%
  pivot_longer(c(MPR, MS, AUC), names_to = "metric") %>%
  ungroup()

distinct_combos <- df_roc_long %>%
  summarize(
    # value = mean(value),
    .by = c(species, Tissue, Assay, Dilutions, cutoff, metric)
  )

get_roc <- function(species, Tissue, Assay, Dilutions, cutoff, metric) {
  temp_df <- df_roc_long %>%
    filter(
      species   == .env$species,
      Tissue    == .env$Tissue,
      Assay     == .env$Assay,
      Dilutions == .env$Dilutions,
      cutoff    == .env$cutoff,
      metric    == .env$metric
    )
  
  if (n_distinct(temp_df$ELISA) < 2) {
    # print(species, Tissue, Assay, Dilutions, cutoff, metric)
    return(NULL)
  }

  temp_roc <- roc(temp_df, ELISA, value, direction = "<")
  temp_roc$species   <- species
  temp_roc$Tissue    <- Tissue
  temp_roc$Assay     <- Assay
  temp_roc$Dilutions <- Dilutions
  temp_roc$cutoff    <- cutoff
  temp_roc$metric    <- metric
  # temp_roc$value     <- value

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
    metric    = x$metric,
    # value     = x$value,
    auc       = x$auc[1]
  )
}

get_roc_info <- function(x) {
  temp_coords <- coords(x) %>%
    mutate(youden = sensitivity + specificity - 1)
  optimal_threshold <- temp_coords$threshold[which(temp_coords$youden == max(temp_coords$youden))]
  tibble(
    species     = x$species,
    Tissue      = x$Tissue,
    Assay       = x$Assay,
    Dilutions   = x$Dilutions,
    cutoff      = x$cutoff,
    metric      = x$metric,
    # value       = x$value,
    threshold   = list(temp_coords$threshold),
    sensitivity = list(temp_coords$sensitivity),
    specificity = list(temp_coords$specificity),
    youden      = list(temp_coords$youden),
    auc = auc(x)[2],
    optimal_threshold = optimal_threshold
  )
}

roc_aucs <- map_dfr(roc_list, get_roc_auc)

roc_coords <- map_dfr(roc_list, get_roc_info) %>%
  unnest(c(sensitivity, specificity, threshold, youden)) %>%
  arrange(sensitivity)

# roc_coords_combined_metrics <- roc_coords %>%
#   select(species, Tissue, Assay, Dilutions, cutoff, sensitivity, specificity, metric) %>%
#   mutate(idx = row_number(), .by = c(metric)) %>%
#   pivot_wider(names_from = "metric", values_from = c("sensitivity", "specificity"),id_cols = c(species, Tissue, Assay, Dilutions, cutoff, idx)) %>%
#   mutate(
#     sensitivity = rowMeans(pick(sensitivity_MPR, sensitivity_MS, sensitivity_AUC), na.rm = TRUE),
#     specificity = rowMeans(pick(specificity_MPR, specificity_MS, specificity_AUC), na.rm = TRUE)
#   ) %>%
#   select(-c(sensitivity_MPR, sensitivity_MS, sensitivity_AUC, specificity_MPR, specificity_MS, specificity_AUC, idx)) %>%
#   group_by(species, Tissue, Assay, Dilutions, cutoff) %>%
#   arrange(sensitivity, desc(specificity))

# ROC curves for top combinations (AUC > 0.75)
top_combos <- roc_aucs %>% 
  # group_by(metric) %>%
  # mutate(idx = row_number(), .groups = "drop") %>%
  # group_by(species, Assay, Tissue, Dilutions, cutoff, idx) %>%
  # summarize(auc = mean(auc), .groups = "drop") %>%
  # filter(auc > 0.75) %>%
  group_by(species, Tissue, Assay, metric) %>%
  slice_max(auc, with_ties = TRUE) %>%
  slice_min(as.numeric(cutoff))

make_roc_plot <- function(x) {
  x %>%
		ggplot(aes(specificity, sensitivity, color = metric,
									linetype = Assay, group = interaction(metric, Assay))) +
    # geom_line() +
		geom_step() +
    # stat_smooth(method="loess", n=10, se=F, fullrange=TRUE) +
    # geom_roc(n.cuts = 0) +
		geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray50", inherit.aes = FALSE) +
		facet_grid(cols = vars(round(Dilutions, 3)), rows = vars(Tissue)) +
    scale_x_reverse(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    coord_fixed() +
		labs(
			x = "Specificity", y = "Sensitivity"
		) +
		main_theme
}

# Moose
# moose_top_combos <- top_combos %>%
#   filter(species == "M" & !(Tissue == "LN" & Dilutions > -1)) %>%
#   summarize(auc = mean(auc), .by=c(species, Tissue, Assay, Dilutions, cutoff)) %>%
#   arrange(desc(auc)) %>%
  # slice_max(auc, by = Tissue, with_ties = TRUE) %>%
  # slice_min(cutoff, by = Tissue, with_ties = FALSE) %>%
  # head(20)

top_moose <- top_combos %>%
  filter(species == "M")

roc_coords %>%
  filter(species == "M", 
    Dilutions %in% top_moose$Dilutions, cutoff %in% top_moose$cutoff
  ) %>%
  mutate(cutoff = as.numeric(cutoff)) %>%
  filter(cutoff == min(cutoff)) %>%
	make_roc_plot() +
  ggtitle("Moose ROC Curves")

ggsave("figures/moose_roc_curves.png", width = 10, height = 8)

# df_roc_long %>%
#   filter(species == "M", cutoff == 36) %>%
#   ggplot(aes(Tissue, value, color = as.logical(ELISA))) +
#   geom_point(position=position_jitterdodge(0.1, dodge.width = 0.5)) +
#   facet_grid(cols = vars(Dilutions), rows = vars(Assay, metric), scales="free") +
#   scale_y_log10()


# Reindeer
top_reindeer <- top_combos %>%
  filter(species == "R", Tissue != "skin")

roc_coords %>%
  filter(species == "R", Dilutions %in% top_reindeer$Dilutions, cutoff %in% top_reindeer$cutoff) %>%
  mutate(cutoff = as.numeric(cutoff)) %>%
  filter(cutoff == min(cutoff)) %>%
	make_roc_plot() +
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


# Step 3: Species-specific Logistic Regression Models --------------------
# Univariate models per metric to avoid multicollinearity (MPR:MS r=0.99, MPR:AUC r=0.89)
# Average across dilutions first: one row per animal x Assay x Tissue,
# so no repeated structure remains and plain glm() suffices

# df_lr <- df_animal %>%
#   group_by(species, animal_id, Assay, Tissue, ELISA) %>%
#   summarise(across(c(MPR, MS, AUC), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
#   mutate(ELISA = as.integer(ELISA)) %>%
#   pivot_longer(c(MPR, MS, AUC), names_to = "metric")

# lr_models <- df_lr %>%
#   group_by(species, metric) %>%
#   mutate(value = scale(value)) %>%
#   group_modify(~{
#     mod <- tryCatch(
#       glm(ELISA ~ value + Assay * Tissue, data = .x, family = binomial),
#       warning = function(w) NULL,
#       error   = function(e) NULL
#     )
#     if (is.null(mod) || !mod$converged) return(data.frame())
#     broom::tidy(mod, exponentiate = TRUE)
#   })

# if (nrow(lr_models) > 0) {
#   lr_models %>%
#     filter(term != "(Intercept)") %>%
#     ggplot(aes(term, estimate)) +
#     geom_point(position = position_dodge(0.5)) +
#     geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
#                   position = position_dodge(0.5), width = 0.2) +
#     geom_hline(yintercept = 1, linetype = "dashed") +
#     facet_grid(rows=vars(metric), scales = "free_y") +
#     coord_flip() +
#     labs(title = "Odds Ratios from Logistic Regression by Species/Tissue",
#          y = "Odds Ratio (95% CI)") +
#     theme_bw()

#   ggsave("figures/logistic_regression_odds_ratios.png", width = 12, height = 8)
# }


# Step 4: Optimal Threshold Analysis -------------------------------------

make_thresh_plot <- function(x) {
    ggplot(x, aes(Assay, threshold, fill = Tissue, label=round(threshold, 2))) +
    geom_col(position = "dodge", width=0.8/5*n_distinct(x$Tissue)) +
    geom_text(position = "dodge") +
    facet_grid(rows=vars(metric), scales = "free") +
    labs(title = "Optimal Decision Thresholds (Youden's J, AUC > 0.75)",
         x = "Metric", y = "Threshold Value") +
    theme_bw()
}

threshold_results <- df_roc_long %>%
  semi_join(roc_results %>% filter(auc > 0.75),
            by = c("species", "Assay", "Tissue", "metric")) %>%
  group_by(species, Assay, Tissue, metric) %>%
  group_modify(~{
    tryCatch({
      r <- roc(.x, ELISA, value, quiet = TRUE)
      coords_df <- coords(r, "best", best.method = "youden",
                          ret = c("threshold", "sensitivity", "specificity"))
      as.data.frame(coords_df)
    }, error = function(e) data.frame(threshold = NA, sensitivity = NA, specificity = NA))
  }) %>%
  filter(!is.na(threshold))

if (nrow(threshold_results) > 0) {
  # Moose
  threshold_results %>%
    filter(species == "M") %>%
    make_thresh_plot()

  ggsave("figures/moose_optimal_thresholds.png", width = 10, height = 7)
  
  # Moose
  threshold_results %>%
    filter(species == "R") %>%
    make_thresh_plot()

  ggsave("figures/reindeer_optimal_thresholds.png", width = 10, height = 7)
  
  # Moose
  threshold_results %>%
    filter(species == "RD") %>%
    make_thresh_plot()

  ggsave("figures/reddeer_optimal_thresholds.png", width = 10, height = 7)

  print(threshold_results)
}


# Step 5: Per-species Mixed Effects Models --------------------------------

# metrics_list <- c("AUC", "MPR", "MS")

# species_lmer <- map(species_list, function(sp) {
#   df_sp <- df_clean %>% filter(species == sp)

#   map(metrics_list, function(met) {
#     formula <- as.formula(paste0(met, " ~ Assay * Tissue + Assay * Dilutions + (1|animal_id) + (1|rxn)"))
#     tryCatch(
#       list(species = sp, metric = met, model = lmer(formula, df_sp)),
#       error = function(e) NULL
#     )
#   }) %>% compact()
# }) %>% flatten()

# species_emm <- map(species_lmer, function(m) {
#   tryCatch({
#     df_sp <- df_clean %>% filter(species == m$species)
#     emm <- emmeans(m$model, ~ Assay | Tissue,
#                    at = list(Dilutions = unique(df_sp$Dilutions)))
#     list(species = m$species, metric = m$metric, emmeans = emm)
#   }, error = function(e) NULL)
# }) %>% compact()

# species_contrasts <- map(species_emm, function(em) {
#   tryCatch({
#     con <- contrast(em$emmeans, interaction=TRUE)
#     list(species = em$species, metric = em$metric, contrast = con)
#   }, error = function(e) NULL)
# }) %>% compact()

