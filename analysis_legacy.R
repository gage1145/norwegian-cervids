# Initialize --------------------------------------------------------------


reps <- paste0("rep", seq(1:12))

file <- "allraw.xlsx"

rep_cols <- paste0("rep", seq(1:12))

df_raw <- read_xlsx(file) %>%
  mutate(
    `Tenfold Dil` = -log10(Dilution)
  ) %>%
  pivot_longer(rep_cols, names_to = "rep", values_to = "RAF") %>%
  filter(!is.na(RAF))

df_raf <- read_xlsx(file) %>%
  mutate(metric = "RAF") %>%
  pivot_longer(all_of(reps), names_to = "reps")

df_mpr <- read_xlsx(file, sheet = 2) %>%
  mutate(metric = "MPR") %>%
  pivot_longer(all_of(reps), names_to = "reps")

df_ms <- read_xlsx(file, sheet = 3) %>%
  mutate(metric = "MS") %>%
  pivot_longer(all_of(reps), names_to = "reps")

df_ttt <- read_xlsx(file, sheet = 4) %>%
  mutate(metric = "TtT") %>%
  pivot_longer(all_of(reps), names_to = "reps")

df_ <- df_raf %>%
  bind_rows(df_mpr) %>%
  bind_rows(df_ms) %>%
  bind_rows(df_ttt) %>%
  filter(!is.na(value)) %>%
  mutate(Tissue = ifelse(Tissue %in% c("SK1", "SK2", "SK3", "SK4", "Ear"), "SK", Tissue))


df_skinsep_ <- df_raf %>%
  bind_rows(df_mpr) %>%
  bind_rows(df_ms) %>%
  bind_rows(df_ttt) %>%
  filter(!is.na(value))


# RT-QuIC vs Nano-QuIC Comparison Boxplots ----------------------------------------------------------------


# Brain boxplot - GOOD

df_ %>%
  filter(metric == "RAF") %>%
  filter(ELISA) %>%
  mutate(ELISA = ifelse(ELISA, "Positive", "Negative")) %>%
  filter(Tissue == "BR") %>%
  filter(Dilution < 10000) %>%
  mutate_at("Dilution", ~ as.factor(-log10(as.numeric(.)))) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Species, value, fill = Assay)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 7.8e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE
  ) +
  facet_grid(rows = vars(Dilution), scale = "free_x") +
  scale_y_continuous(limits = c(0, 9.5e-5)) +
  scale_fill_manual(values = c("#ff9c59", "purple")) +
  ggtitle(str_wrap("RT-QuIC vs Nano-QuIC Positive Brain RAFs", width = 40)) +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  main_theme +
  theme(axis.title.x = element_blank())


ggsave("brain_bp.v3.png", width = 18, height = 9)
ggsave("brain_bp_half.v3.png", width = 9, height = 9)


# Lymph Nodes -------------------------------------------------------------


# Lymph Node boxplot - GOOD
df_ %>%
  filter(metric == "RAF") %>%
  filter(ELISA) %>%
  mutate(ELISA = ifelse(ELISA, "Positive", "Negative")) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  filter(Tissue == "LN") %>%
  filter(Dilution < 1000) %>%
  mutate_at("Dilution", ~ as.factor(-log10(as.numeric(.)))) %>%
  ggplot(aes(Species, value, fill = Assay)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.5e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE
  ) +
  facet_grid(rows = vars(Dilution), scale = "free_x") +
  scale_y_continuous(limits = c(0, 9.5e-5)) +
  scale_fill_manual(values = c("#ff9c59", "purple")) +
  ggtitle("RT-QuIC vs Nano-QuIC Positive Lymph Node RAFs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("ln_bp.v3.png", width = 18, height = 9)
ggsave("ln_bp_half.v3.png", width = 9, height = 9)


# Skin --------------------------------------------------------------------


# Combined Moose Skin boxplot - GOOD
df_ %>%
  filter(metric == "RAF") %>%
  # filter(ELISA) %>%
  mutate(ELISA = ifelse(ELISA, "Positive", "Negative")) %>%
  filter(Tissue %in% c("SK")) %>%
  filter(Dilution < 1000) %>%
  filter(Species == "Moose") %>%
  mutate_at("Dilution", ~ as.factor(-log10(as.numeric(.)))) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Dilution, value, fill = Assay)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.5e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(ELISA), scale = "free_x") +
  scale_fill_manual(values = c("#ff9c59", "purple")) +
  scale_y_continuous(limits = c(0, 9.5e-5)) +
  ggtitle("RT-QuIC vs Nano-QuIC Positive Moose Skin RAFs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  main_theme
ggsave("moose_skin_v3.png", width = 18, height = 9)
ggsave("moose_skin_half_v3.png", width = 9, height = 9)


# Individual Skins Moose figure


df_skinsep_ %>%
  filter(metric == "RAF") %>%
  filter(ELISA) %>%
  mutate(ELISA = ifelse(ELISA, "Positive", "Negative")) %>%
  filter(Tissue %in% c("SK1", "SK2", "SK3", "SK4", "Ear")) %>%
  filter(Dilution < 1000) %>%
  mutate_at("Dilution", ~ as.factor(-log10(as.numeric(.)))) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Tissue, value, fill = Assay)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.5e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE
  ) +
  facet_grid(rows = vars(Dilution), scale = "free_x") +
  scale_y_continuous(limits = c(0, 9.5e-5)) +
  scale_fill_manual(values = c("#ff9c59", "purple")) +
  ggtitle("RT-QuIC vs Nano-QuIC Positive Moose Skins RAFs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("moose_indvskin_v3.png", width = 18, height = 9)
ggsave("moose_indvskin_half_v3.png", width = 9, height = 9)


# Positive and Negative Comparison Boxplots ---------------------------------------------------


# Brain
df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution)
  ) %>%
  filter(metric == "RAF", Tissue == "BR", Dilution %in% c(-1, -2, -3)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Brain RAFs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("brain_comps_v3.png", width = 18, height = 9)
ggsave("brain_comps_half_v3.png", width = 9, height = 9)


# Lymph Nodes
df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution)
  ) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  filter(metric == "RAF", Tissue == "LN", Dilution %in% c(-1, -2)) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 7.9e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Lymph Node RAFs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("ln_comps_v3.png", width = 18, height = 9)
ggsave("ln_comps_half_v3.png", width = 9, height = 9)


#   All Species Skin -> doesnt show stats because Reindeer doesnt have neg

df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution)
  ) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  filter(metric == "RAF", Tissue == "SK", Dilution %in% c(-1, -2)) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 7.9e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Skin RAFs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("all_sk_comps_v1.png", width = 18, height = 9)
ggsave("all_sk_comps_half_v1.png", width = 9, height = 9)


# Just Moose Skin for Stats to Add to Combined Graph

df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution)
  ) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  filter(metric == "RAF", Species == "Moose", Tissue == "SK", Dilution %in% c(-1, -2)) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 7.9e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Moose Skin RAFs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("moose_sk_comps_v1.png", width = 18, height = 9)
ggsave("moose_sk_comps_half_v1.png", width = 9, height = 9)


# Combined figure ---------------------------------------------------------


# Combined LN BR figure
df_ %>%
  filter(metric == "RAF") %>%
  filter(ELISA) %>%
  mutate(ELISA = ifelse(ELISA, "Positive", "Negative")) %>%
  filter(Tissue %in% c("LN", "BR")) %>%
  filter(Dilution < 10000) %>%
  mutate_at("Dilution", ~ as.factor(-log10(as.numeric(.)))) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Species, value, fill = Assay)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.5e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  # scale_fill_manual(values=c(""))
  facet_grid(rows = vars(Dilution), cols = vars(Tissue), scale = "free_x") +
  scale_y_continuous(limits = c(0, 9.5e-5)) +
  scale_fill_manual(values = c("#ff9c59", "purple")) +
  ggtitle("RT-QuIC vs Nano-QuIC Positive Brain and Lymph Node RAFs") +
  labs(
    # y = "Maxpoint Ratio"
    y = "Rate of amyloid formation (1/s)"
  ) +
  main_theme +
  theme(axis.title.x = element_blank())
ggsave("combined_lnbr_v3.png", width = 18, height = 9)
ggsave("combined_lnbr_half_v3.png", width = 9, height = 9)


# M10 Graph for Supplementary Figures


df_ %>%
  filter(metric == "RAF") %>%
  filter(Animal == "M10") %>%
  filter(Tissue == "BR") %>%
  filter(Dilution < 1000) %>%
  mutate_at("Dilution", ~ as.factor(log10(as.numeric(.)))) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Dilution, value, fill = Assay)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 1.5e-5,
    tip.length = 0.01,
    size = 15,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  scale_fill_manual(values = c("#ff9c59", "purple")) +
  scale_y_continuous(limits = c(0, 1.7e-5)) +
  ggtitle("RAFs of Moose Subset Brain") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  main_theme
ggsave("m10_NQ_example_v1.png", width = 18, height = 9)
ggsave("m10_NQ_example_half_v1.png", width = 9, height = 9)


# Graphs to determine Cutoff point for false seeding


df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution)
  ) %>%
  filter(metric == "TtT", Tissue == "BR", Dilution %in% c(-1, -2, -3)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Species, value, color = ELISA)) +
  geom_point() +
  geom_jitter(width = 0.2, height = 0) +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Brain RAFs") +
  labs(
    y = "Time to Threshold"
  ) +
  scale_fill_manual(values = c("red", "blue")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("brain_comps_timetothresh.png", width = 18, height = 9)
ggsave("brain_comps_half_timetothresh.png", width = 9, height = 9)


df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution)
  ) %>%
  filter(metric == "TtT", Tissue == "LN", Dilution %in% c(-1, -2)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(Species, value, color = ELISA)) +
  geom_point() +
  geom_jitter(width = 0.2, height = 0) +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive LN RAFs") +
  labs(
    y = "Time to Threshold"
  ) +
  scale_fill_manual(values = c("red", "blue")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("ln_comps_half_timetothresh.png", width = 9, height = 9)


# Different Thresholds Example Graphs of Above Pos vs Neg Graphs -----------------------


# Threshold of 30 hours

df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution),
  ) %>%
  filter(metric == "RAF", Tissue == "BR", Dilution %in% c(-1, -2, -3)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  mutate(
    value = if_else(value < 9.259259e-6, 0, value)
  ) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Brain RAFs 30 hrs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("brain_comps_half_30hr_thresh.png", width = 9, height = 9)


df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution),
  ) %>%
  filter(metric == "RAF", Tissue == "LN", Dilution %in% c(-1, -2)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  mutate(
    value = if_else(value < 9.259259e-6, 0, value)
  ) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive LN RAFs 30 hrs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("ln_comps_half_30hr_thresh.png", width = 9, height = 9)


# Threshold of 35 hours


df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution),
  ) %>%
  filter(metric == "RAF", Tissue == "BR", Dilution %in% c(-1, -2, -3)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  mutate(
    value = if_else(value < 7.9365079e-6, 0, value)
  ) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Brain RAFs 35 hrs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("brain_comps_half_35hr_thresh.png", width = 9, height = 9)


df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution),
  ) %>%
  filter(metric == "RAF", Tissue == "LN", Dilution %in% c(-1, -2)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  mutate(
    value = if_else(value < 7.9365079e-6, 0, value)
  ) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive LN RAFs 35 hrs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("ln_comps_half_35hr_thresh.png", width = 9, height = 9)


# Threshold of 40 hours

df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution),
  ) %>%
  filter(metric == "RAF", Tissue == "BR", Dilution %in% c(-1, -2, -3)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  mutate(
    value = if_else(value < 6.944444444444444444e-6, 0, value)
  ) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Brain RAFs 40 hrs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("brain_comps_half_40hr_thresh.png", width = 9, height = 9)


df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution),
  ) %>%
  filter(metric == "RAF", Tissue == "LN", Dilution %in% c(-1, -2)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  mutate(
    value = if_else(value < 6.944444444444444444e-6, 0, value)
  ) %>%
  ggplot(aes(Species, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive LN RAFs 40 hrs") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.title.x = element_blank())

ggsave("ln_comps_half_40hr_thresh.png", width = 9, height = 9)


## Time to Thresh by Animal ID Graphs-----------------------------------------


# Moose

df_ %>%
  mutate(
    ELISA = ifelse(ELISA, "Positive", "Negative"),
    Dilution = -log10(Dilution),
  ) %>%
  filter(metric == "RAF", Tissue == "BR", Species == "Moose", Dilution %in% c(-1, -2, -3)) %>%
  mutate(
    Assay = case_when(
      Assay == "RT" ~ "RT-QuIC",
      Assay == "Nano" ~ "Nano-QuIC",
      TRUE ~ Assay # keep all other values unchanged
    )
  ) %>%
  ggplot(aes(ID, value, fill = ELISA)) +
  geom_boxplot() +
  stat_compare_means(
    method = "anova",
    label = "p.signif",
    label.y = 8.3e-5,
    tip.length = 0.01,
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    hide.ns = TRUE,
  ) +
  facet_grid(cols = vars(Assay), rows = vars(Dilution)) +
  ggtitle("Negative vs Positive Moose Time to Threshold") +
  labs(
    y = "Rate of amyloid formation (1/s)"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  main_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("br_moose_thresh.png", width = 9, height = 9)
