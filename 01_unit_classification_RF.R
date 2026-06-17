# =============================================================================
# 01_unit_classification_RF.R
#
# Random Forest validation of humpback whale song unit classification
#
# Part of the analysis reported in:
#   Magnúsdóttir et al. (in press). Subarctic feeding grounds as cultural
#   exchange hubs in the North Atlantic humpback whale song network.
#   
#
# Description:
#   Validates human-assigned unit classifications using a Random Forest
#   classifier (ranger package). Unit type is used as the response variable
#   and 11 acoustic parameters as predictors, following Garland et al. (2017)
#   and Dunlop et al. (2007). Out-of-bag (OOB) error rate is reported as
#   the primary validation metric.
#
# Input:
#   A tab-delimited file of unit measurements with columns including:
#     - Unit        : human-assigned unit type label
#     - Freq_Min, Freq_Max, Dur, Peak_freq, StartF, EndF,
#       RangeF, TrendF, Inflections, Pulse_rate  : acoustic parameters
#     - LocYY       : location-year identifier (e.g. ICE11, CAR09)
#     - SongUnitID  : unique identifier per unit
#
# Output:
#   - Printed OOB error rate and confusion matrix
#
# Authors: Edda Elísabet Magnúsdóttir et al.
# Contact: edda@hi.is
# =============================================================================


# -----------------------------------------------------------------------------
# SETTINGS — adjust these for your data
# -----------------------------------------------------------------------------

# Path to unit measurement data
input_file <- "data/Song_Unit_Frequency_data_combined_results.txt"

# Unit IDs to exclude (quality control)
bad_ids <- c("Name of unit to delete from dataset")

# Random Forest parameters
n_trees  <- 10000
mtry_val <- 3
seed_val <- 123

# -----------------------------------------------------------------------------
# PACKAGES
# -----------------------------------------------------------------------------

library(dplyr)
library(stringr)
library(ranger)

# -----------------------------------------------------------------------------
# LOAD AND CLEAN DATA
# -----------------------------------------------------------------------------

x <- read.delim(input_file, header = TRUE)

# Standardise location-year identifier
standardize_locyy <- function(x) {
  x     <- as.character(x)
  Loc   <- str_extract(x, "^[A-Za-z]+")
  Year2 <- dplyr::case_when(
    str_detect(x, "\\d{4}") ~ str_sub(str_extract(x, "\\d{4}"), 3, 4),
    str_detect(x, "\\d{2}") ~ str_extract(x, "\\d{2}"),
    TRUE ~ NA_character_
  )
  paste0(toupper(Loc), Year2)
}

x <- x %>%
  mutate(
    LocYY = standardize_locyy(Loc_Year),
    Unit  = Unit %>% str_trim() %>% str_squish() %>% str_to_lower(),
    Unit  = case_when(
      Unit %in% c("mup_short", "mup3") ~ "mup3",
      Unit %in% c("li1_short", "li1")  ~ "li1",
      TRUE ~ Unit
    ),
    Pulse_rate = ifelse(is.na(Pulse_rate), 0, Pulse_rate)
  )

# Impute missing acoustic values with column means
num_cols <- c("Freq_Min", "Freq_Max", "Dur", "Peak_freq", "StartF", "EndF",
              "RangeF", "TrendF", "Inflections")
x <- x %>%
  mutate(across(all_of(num_cols),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  mutate(Unit = droplevels(as.factor(Unit)))

# Remove flagged bad IDs
units <- x %>% filter(!SongUnitID %in% bad_ids)

cat("Units loaded:", nrow(units), "\n")
cat("Unit types:  ", nlevels(units$Unit), "\n")

# -----------------------------------------------------------------------------
# RANDOM FOREST CLASSIFICATION
# -----------------------------------------------------------------------------

vars <- c("Unit", "Freq_Min", "Freq_Max", "Dur", "Peak_freq",
          "StartF", "EndF", "RangeF", "TrendF", "Inflections", "Pulse_rate")

x_rf <- units %>%
  select(all_of(vars)) %>%
  filter(complete.cases(.)) %>%
  mutate(Unit = droplevels(as.factor(Unit))) %>%
  as.data.frame()

set.seed(seed_val)
train_idx  <- sample(1:nrow(x_rf), 0.7 * nrow(x_rf))
train_data <- x_rf[train_idx, ]
test_data  <- x_rf[-train_idx, ]

RFfit <- ranger(
  Unit ~ Freq_Min + Freq_Max + Dur + Peak_freq +
         StartF + EndF + RangeF + TrendF +
         Inflections + Pulse_rate,
  data      = train_data,
  num.trees = n_trees,
  mtry      = mtry_val
)

cat("\n==== Random Forest results ====\n")
print(RFfit)

cat("\n==== Confusion matrix (full dataset) ====\n")
pred     <- predict(RFfit, data = x_rf)$predictions
conf_mat <- table(observed = x_rf$Unit, predicted = pred)
print(conf_mat)

cat("\nOOB error rate:", round(RFfit$prediction.error * 100, 1), "%\n")
