library(quicR)
library(tidyverse)
library(readxl)
library(pROC)
library(janitor)


# Setup ------------------------------------------------------------------


main_theme <- theme(
  plot.title = element_text(size = 30),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 24),
  legend.text = element_text(size = 20),
  strip.text = element_text(size = 24)
)

quic_colors  <- c("Nano-QuIC" = "darkorange", "RT-QuIC" = "purple")
elisa_colors <- c("Negative" = "navy", "Positive" = "red")

# Load in data -----------------------------------------------------------


files <- list.files("mnquic_raw", full.names=TRUE)

get_raw <- function(file) {
  rxn_str <- str_split_i(file, "/", 2) %>%
    str_remove(".xlsx")

  date_str <- str_match(rxn_str, "\\d{1,2}-\\d{1,2}-\\d{2}")[[1]] %>%
    parse_date(format = "%m-%d-%y")

  file %>%
    organize_tables() %>%
    convert_tables() %>%
    mutate(date = date_str, rxn = rxn_str)
}

df_ <- map_dfr(files, get_raw)

df_meta <- read_xlsx("metadata/allraw.xlsx") %>%
  select(Animal, ID, ELISA, Genotype) %>%
  filter(!is.na(Animal)) %>%
  group_by(Animal) %>%
  reframe(across(everything(), unique))

species_key <- c(M = "Moose", R = "Reindeer", RD = "Red Deer")
tissue_key <- c(BR = "Brain", LN = "Lymph Node")

df_clean <- df_ %>%
  separate_wider_delim(
    "Sample IDs", "-", names=c("Assay", "Animal", "Tissue", "Dilutions"), too_few = "align_start", names_repair = "none"
  ) %>%
  select("Animal", "Assay", "Tissue", "Dilutions", "Ratio580", "RFU", "Wells", "date", "rxn") %>%
  filter(
    Assay != "Empty",
    Animal != "2326",
    !(tolower(Animal) %in% c("neg", "pos"))
  ) %>%
  mutate(
    # Multiply by 10 to account for Peter C's dilution nomenclature.
    Dilutions = -log10(as.numeric(Dilutions) * 10),
    animal_id = str_match(Animal, "\\d+$"),
    species = str_match(Animal, "^[[:alpha:]]+"),
    Ratio580 = as.numeric(Ratio580),
    RFU = as.numeric(RFU),
    species = species_key[species],
    Tissue = tissue_key[Tissue]
  ) %>%
  relocate(species, animal_id, .after = Animal) %>%
  left_join(df_meta) %>%
  clean_names()


# ROC --------------------------------------------------------------------


## Compute ROC objects ---------------------------------------------------


distinct_combos <- distinct(df_clean, species, tissue)

# Fit one ROC curve per (species, tissue, assay, dilution, cutoff) cell.
# direction = ">" tells pROC that lower pc1 = ELISA positive (since pc1 was
# built from raw QuIC metrics where positives sit at the negative end of pc1).
# Cells without both ELISA classes are skipped; meta-fields are stashed on
# the returned object so downstream code can recover the cell identity.
get_roc <- function(species, tissue) {
  meta <- as.list(environment())

  temp_df <- df_clean %>%
    filter(
      species == .env$species,
      tissue == .env$tissue
    ) %>%
    mutate(elisa = as.integer(elisa))

  if (n_distinct(temp_df$elisa) < 2) {
    return(NULL)
  }

  temp_roc <- roc(temp_df, elisa, ratio580, direction = "<")
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
    threshold   = list(temp_coords$threshold),
    sensitivity = list(temp_coords$sensitivity),
    specificity = list(temp_coords$specificity),
    roc_auc     = as.numeric(auc(x))
  )
}

roc_coords <- map_dfr(roc_list, get_roc_info) %>%
  unnest(c(sensitivity, specificity, threshold)) %>%
  group_by(species, tissue) %>%
  arrange(sensitivity, .by_group = TRUE) %>%
  ungroup()

roc_aucs <- roc_coords %>%
  distinct(species, tissue, roc_auc) %>%
  mutate(roc_auc = round(roc_auc, 3))

write.csv(roc_aucs, "data/mnquic_aucs.csv")


# Figures ----------------------------------------------------------------


# 580/ThT Relationship
df_clean %>%
  mutate(elisa = ifelse(elisa, "Positive", "Negative")) %>%
  ggplot(aes(y=ratio580, x=log(rfu), color=elisa)) +
  geom_point(size=2, alpha=0.7) +
  facet_grid(cols = vars(species), rows=(vars(tissue))) +
  scale_color_manual(values = c("blue", "red")) +
  # scale_x_log10() +
  labs(
    x = "Log(RFU)",
    y = "580 Ratio"
  ) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

df_clean %>%
  mutate(elisa = as.integer(elisa)) %>%
  ggplot(aes(ratio580, elisa)) +
  geom_point(size=2) +
  geom_smooth(method = "glm", 
    method.args = list(family = "binomial"), 
    se = T) +
  facet_grid(cols = vars(species), rows = vars(tissue)) 


# ROC Figure
roc_coords %>%
ggplot(aes(specificity, sensitivity)) +
  geom_step(linewidth = 1.2) +
  geom_abline(
    slope = 1, intercept = 1, linetype = "dashed",
    color = "black", inherit.aes = FALSE
  ) +
  geom_label(
    aes(label = paste0("AUC = ", round(roc_auc, 2))),
    x = -0.3, y=0.3, data = roc_aucs, inherit.aes = FALSE, hjust = 0.5, size = 6
  ) +
  facet_grid(cols = vars(species), rows = vars(tissue)) +
  # scale_color_manual(values = quic_colors) +
  scale_x_reverse(
    limits = c(0, 1), breaks = seq(0, 1, 0.2),
    # sec.axis = sec_axis(~., name = "Log Dilution Factor", breaks = NULL, labels = NULL)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Specificity", y = "Sensitivity", title = "MN-QuIC ROC Curves"
  ) +
  main_theme +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "none",
    legend.title = element_blank(),
  )
ggsave("mnquic_roc.png", path = "figures", height = 10, width = 16)
