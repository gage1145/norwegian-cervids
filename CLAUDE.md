# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based biomedical data analysis project for comparing RT-QuIC vs Nano-QuIC prion disease diagnostic assays in wildlife species (moose, reindeer, red deer). It processes raw fluorescence data from QuIC experiments and produces statistical analyses and publication-ready figures.

## Running the Code

```r
# Run the main data processing pipeline (produces clean_data.csv + ROC/mixed model analysis)
source("analysis.R")
```

Required R packages: `tidyverse`, `readxl`, `skimr`, `quicR`, `pROC`, `lme4`, `emmeans`, `ggplot2`

## Architecture & Data Flow

**Input**:
- `formatted_raw/*.xlsx` — 72 raw QuIC run files; filenames encode metadata as `[Date]_ZS_[Animal]_[Tissue]_[Condition].xlsx`
- `metadata/allraw.xlsx` — master metadata with animal ID, ELISA status, genotype

**Pipeline** (`analysis.R`):
1. Reads all xlsx files via `quicR::get_quic()` which extracts fluorescence metrics (MPR, MS, AUC)
2. Parses filename components to extract date, assay type, animal, tissue
3. Splits Sample ID column into Assay, Animal, Tissue, Dilutions fields
4. Converts dilutions to -log10 scale; maps species codes (M=Moose, R=Reindeer, RD=Red Deer)
5. Merges with `metadata/allraw.xlsx`
6. Writes tidy output to `data/clean_data.csv`
7. Runs ROC analysis (via `pROC`) across all assay/tissue/dilution/species/metric combinations
8. Fits linear mixed models (`lmer`) with animal_id and rxn as random effects; uses `emmeans` for pairwise comparisons

**Output**: `data/clean_data.csv` and PNG figures

**`analysis_legacy.R`** (906 lines): Extended exploratory script with boxplots faceted by tissue/dilution/species, threshold sensitivity analysis (30hr, 35hr, 40hr cutoffs for false seeding), and statistical annotations. Color scheme: `#ff9c59` (orange) = RT-QuIC, purple = Nano-QuIC.

## Key Domain Concepts

- **QuIC assays**: RT-QuIC (Real-Time) vs Nano-QuIC — two methods for detecting misfolded prion proteins
- **Metrics**: MPR (Maxpoint Ratio), MS (Maxpoint Slope), AUC (Area Under Curve), RAF (Rate of Amyloid Formation), TtT (Time to Threshold)
- **Samples**: Brain, Lymph Node, Skin tissue types; dilution series tested on each
- **Controls**: `negcontrol` samples are filtered out before analysis; `empty` assays are excluded