library(tidyverse)
library(quicR)
library(readxl)
library(skimr)


files <- list.files("formatted_raw", ".xlsx", full.names = TRUE)

get_raw <- function(file) {
  rxn_str <- str_split_i(file, "/", 2) %>%
    str_remove(".xlsx")

  date_str <- str_match(rxn_str, "\\d{1,2}-\\d{1,2}-\\d{2}")[[1]] %>%
    parse_date(format = "%m-%d-%y")
  
  data <- get_quic(file, norm_point = 4, window_size = 3) %>%
    suppressMessages() %>%
    mutate(date = date_str, rxn = rxn_str)

  return(data)
}

df_ <- map_dfr(files, get_raw, .progress=TRUE)

df_formatted <- df_ %>%
  separate(`Sample IDs`, c("Assay", "Animal", "Tissue", "Dilutions"), "-", extra = "merge", fill = "right")

skim(df_formatted)

df_metrics <- df_formatted %>%
  calculate_metrics(c("Assay", "Animal", "Tissue", "Dilutions", "Wells", "date", "rxn"), threshold = 4)

# species_key <- list(M = "Moose", R = "Reindeer", RD = "Red Deer")

df_meta <- read_xlsx("metadata/allraw.xlsx") %>%
  select(Animal, ID, ELISA, Genotype) %>%
  filter(!is.na(Animal)) %>%
  group_by(Animal) %>%
  reframe(across(everything(), unique))

df_clean <- df_metrics %>%
  filter(
    Assay != "Empty",
    Animal != "2326",
    !(tolower(Animal) %in% c("neg", "pos"))
  ) %>%
  mutate(
    Dilutions = -log10(as.numeric(Dilutions)),
    animal_id = str_match(Animal, "\\d+$"),
    species = str_match(Animal, "^[[:alpha:]]+")
  ) %>%
  filter(
    Dilutions > -3
  ) %>%
  relocate(species, animal_id, .after = Animal) %>%
  left_join(df_meta)
skim(df_clean)

write.csv(df_clean, "data/clean_data.csv")


# Analysis ---------------------------------------------------------------


# ROC Analysis -----------------------------------------------------------

library(pROC)

df_roc <- df_clean %>%
  pivot_longer(c(MPR, MS, AUC), names_to = "metric") %>%
  mutate(ELISA = as.integer(ELISA))

roc_list <- list()
names_list <- c()

for (assay in unique(df_roc$Assay)) {
  for (tissue in unique(df_roc$Tissue)) {
    for (dilution in unique(df_roc$Dilutions)) {
      for (spec in unique(df_roc$species)) {
        for (met in unique(df_roc$metric)) {
          roc_name <- paste(assay, tissue, dilution, spec, met, sep="_")
          print(roc_name)
          df_sub <- df_roc %>%
            filter(Assay == assay, Tissue == tissue, Dilutions == dilution, species == spec, metric == met)
          print(nrow(df_sub))
          tryCatch(
            {
              sub_roc <- roc(df_sub, ELISA, value)
              roc_list <- append(roc_list, list(sub_roc))
              names_list <- c(names_list, roc_name)
            }, 
            error = function(e) return()
          )
        }
      }
    }
  }
}
names(roc_list) <- names_list

aucs <- lapply(roc_list, auc) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("roc") %>%
  separate(
    roc, c("assay", "tissue", "dilution", "species", "metric"), "_"
  ) %>%
  mutate(dilution = as.factor(round(as.numeric(str_replace(dilution, ".", "-")), 2)))

# Brain
aucs %>%
  filter(tissue == "BR") %>%
  ggplot(aes(dilution, V1, fill=metric)) +
  geom_col(position = "dodge") +
  facet_grid(species ~ assay, space="free", scales="free")

# Lymph node
aucs %>%
  filter(tissue == "LN") %>%
  ggplot(aes(dilution, V1, fill=metric)) +
  geom_col(position = "dodge") +
  facet_grid(species ~ assay, space="free", scales="free")

# Overall
aucs %>%
  # filter(tissue == "LN") %>%
  ggplot(aes(dilution, V1, fill=metric)) +
  geom_col(position = "dodge") +
  facet_grid(species ~ assay, space="free", scales="free")














library(lme4)
library(emmeans)


# Moose analysis ---------------------------------------------------------

df_moose <- df_clean %>%
  filter(species == "M")


auc_formula <- AUC ~ Assay * Tissue + Assay * Dilutions + (1 | animal_id) + (1 | rxn)

auc_mod <- lmer(auc_formula, df_moose)
summary(auc_mod)

auc_emm <- emmeans(auc_mod, ~ Assay | Tissue + Assay | Dilutions, at = list(Dilutions = unique(df_moose$Dilutions)))
summary(auc_emm)
