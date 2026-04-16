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

get_metrics_w_cutoff <- function(c) {
  df_formatted %>%
    filter(Time <= c) %>%
    mutate(cutoff = c) %>%
    calculate_metrics(c("Assay", "Animal", "Tissue", "Dilutions", "Wells", "cutoff", "date", "rxn"), threshold = 4)
}

cutoffs <- seq(24, 48, 4)

df_metrics <- map_dfr(cutoffs, get_metrics_w_cutoff)

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
  relocate(species, animal_id, .after = Animal) %>%
  left_join(df_meta)
skim(df_clean)

write.csv(df_clean, "data/clean_data.csv")
