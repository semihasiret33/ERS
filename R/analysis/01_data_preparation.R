# ============================================================================
# Script: 01_data_preparation.R
# Purpose: Load and prepare PISA 2018/2022 data for ERS analysis
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(haven)  # For reading SPSS files

# 2. Configuration ----

# Specify which PISA cycle to use
PISA_CYCLE <- "2018"  # or "2022"

# Data file path (update this to your actual data location)
DATA_PATH <- "data/raw/PISA2018_student.sav"

# Countries to include in analysis
# Select countries with known ERS differences
# High ERS: Some Asian and Latin American countries
# Low ERS: Some Western European countries

SELECTED_COUNTRIES <- c(
  "CHN",  # China (typically high ERS)
  "JPN",  # Japan (typically high ERS)
  "MEX",  # Mexico (typically high ERS)
  "DEU",  # Germany (typically low ERS)
  "SWE",  # Sweden (typically low ERS)
  "NLD"   # Netherlands (typically low ERS)
)

# Select scale to analyze
# Example: Student Belonging (ST034) or Learning Goals (ST208)
# Update based on actual PISA questionnaire items

# For PISA 2018, we'll use Sense of Belonging items (ST034Q01-ST034Q06)
SCALE_ITEMS <- c(
  "ST034Q01TA",  # I feel like an outsider (or left out of things)
  "ST034Q02TA",  # I make friends easily
  "ST034Q03TA",  # I feel like I belong
  "ST034Q04TA",  # I feel awkward and out of place
  "ST034Q05TA",  # Other students seem to like me
  "ST034Q06TA"   # I feel lonely
)

# Note: Update item names based on actual PISA 2018/2022 codebook


# 3. Load Data ----

cat("Loading PISA data from:", DATA_PATH, "\n")

# Check if file exists
if (!file.exists(DATA_PATH)) {
  cat("\n*** WARNING: Data file not found! ***\n")
  cat("Please download PISA data from: https://www.oecd.org/pisa/data/\n")
  cat("And update DATA_PATH in this script.\n\n")
  cat("For PISA 2018 Student Questionnaire:\n")
  cat("File: CY07_MSU_STU_QQQ.sav or similar\n\n")
  cat("Creating example placeholder data for demonstration...\n\n")

  # Create placeholder data for demonstration
  set.seed(123)
  n_per_country <- 500

  pisa_data <- expand_grid(
    CNT = SELECTED_COUNTRIES,
    student_id = 1:n_per_country
  ) %>%
    mutate(
      # Simulate responses (1-4 scale)
      across(
        all_of(SCALE_ITEMS),
        ~sample(1:4, n(), replace = TRUE)
      ),
      # Add some realistic patterns
      # High ERS countries tend toward extremes
      .by = CNT
    )

  cat("Using placeholder data for demonstration.\n")
  cat("Replace with actual PISA data for real analysis.\n\n")

} else {

  # Load actual PISA data
  pisa_raw <- read_sav(DATA_PATH)

  cat("Data loaded successfully.\n")
  cat("Total observations:", nrow(pisa_raw), "\n")
  cat("Total variables:", ncol(pisa_raw), "\n\n")

  # Extract relevant variables
  pisa_data <- pisa_raw %>%
    select(
      CNT,              # Country code
      CNTSTUID,         # Student ID
      all_of(SCALE_ITEMS)  # Scale items
    )
}


# 4. Data Cleaning ----

cat("Cleaning data...\n")

# Convert to data frame
pisa_data <- as.data.frame(pisa_data)

# Filter for selected countries
pisa_data <- pisa_data %>%
  filter(CNT %in% SELECTED_COUNTRIES)

cat("Countries included:", paste(SELECTED_COUNTRIES, collapse = ", "), "\n")
cat("Observations after country filter:", nrow(pisa_data), "\n\n")

# Handle missing values
# PISA uses specific codes for missing: often 7, 8, 9 or higher values
# Convert to NA

for (item in SCALE_ITEMS) {
  if (item %in% names(pisa_data)) {
    # Assume valid responses are 1-4 (typical Likert scale)
    pisa_data[[item]][pisa_data[[item]] > 4 | pisa_data[[item]] < 1] <- NA
  }
}

# Remove cases with all items missing
pisa_data <- pisa_data %>%
  filter(rowSums(!is.na(select(., all_of(SCALE_ITEMS)))) > 0)

cat("Observations after removing all-missing cases:", nrow(pisa_data), "\n")

# Optionally: Remove cases with too many missing items
# For example, require at least half of items to be non-missing
min_items_required <- ceiling(length(SCALE_ITEMS) / 2)

pisa_data <- pisa_data %>%
  filter(rowSums(!is.na(select(., all_of(SCALE_ITEMS)))) >= min_items_required)

cat("Observations after requiring", min_items_required, "valid items:", nrow(pisa_data), "\n\n")


# 5. Reverse Coding (if needed) ----

# Check if any items need reverse coding
# For Belonging scale, negatively worded items should be reverse coded
# Example: "I feel like an outsider" should be reverse coded

# Items to reverse code (update based on codebook)
REVERSE_ITEMS <- c("ST034Q01TA", "ST034Q04TA", "ST034Q06TA")

cat("Reverse coding items:", paste(REVERSE_ITEMS, collapse = ", "), "\n")

for (item in REVERSE_ITEMS) {
  if (item %in% names(pisa_data)) {
    # Reverse code: 1->4, 2->3, 3->2, 4->1
    pisa_data[[item]] <- 5 - pisa_data[[item]]
  }
}

cat("Reverse coding completed.\n\n")


# 6. Descriptive Statistics ----

cat("=== Descriptive Statistics ===\n\n")

# Sample sizes by country
sample_sizes <- pisa_data %>%
  count(CNT) %>%
  rename(n_students = n)

cat("Sample sizes by country:\n")
print(sample_sizes)
cat("\n")

# Item correlations
cat("Item intercorrelations:\n")
item_cors <- cor(pisa_data[, SCALE_ITEMS], use = "pairwise.complete.obs")
print(round(item_cors, 2))
cat("\n")

# Overall scale score
pisa_data$scale_score_raw <- rowMeans(pisa_data[, SCALE_ITEMS], na.rm = TRUE)

cat("Overall scale score distribution:\n")
print(summary(pisa_data$scale_score_raw))
cat("\n")

# Country means (raw scores)
country_means_raw <- pisa_data %>%
  group_by(CNT) %>%
  summarise(
    mean_raw = mean(scale_score_raw, na.rm = TRUE),
    sd_raw = sd(scale_score_raw, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_raw))

cat("Country means (raw scores):\n")
print(country_means_raw)
cat("\n")


# 7. Save Prepared Data ----

# Create processed data directory if it doesn't exist
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}

# Save prepared data
saveRDS(pisa_data, file = "data/processed/pisa_clean.rds")
write_csv(pisa_data, file = "data/processed/pisa_clean.csv")

cat("Prepared data saved to:\n")
cat("  - data/processed/pisa_clean.rds\n")
cat("  - data/processed/pisa_clean.csv\n\n")

# Save country means
write_csv(country_means_raw, file = "output/tables/pisa_country_means_raw.csv")

cat("Country means saved to output/tables/pisa_country_means_raw.csv\n\n")


# 8. Export Metadata ----

metadata <- list(
  pisa_cycle = PISA_CYCLE,
  scale_name = "Sense of Belonging",
  items = SCALE_ITEMS,
  reverse_coded = REVERSE_ITEMS,
  countries = SELECTED_COUNTRIES,
  n_total = nrow(pisa_data),
  n_items = length(SCALE_ITEMS),
  date_processed = Sys.time()
)

saveRDS(metadata, file = "data/processed/metadata.rds")

cat("Metadata saved to data/processed/metadata.rds\n")


# 9. Session Info ----
cat("\n=== Session Info ===\n")
cat("Script completed:", as.character(Sys.time()), "\n")
