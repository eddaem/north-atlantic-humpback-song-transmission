# =============================================================================
# 05_transition_entropy.R
#
# Shannon transition entropy of humpback whale theme sequences
#
# Part of the analysis reported in:
#   Magnúsdóttir et al. (in press). Subarctic feeding grounds as cultural
#   exchange hubs in the North Atlantic humpback whale song network.
#   Royal Society Open Science.
#
# Description:
#   Calculates Shannon transition entropy from first-order Markov transition
#   probability matrices. Entropy is computed per theme and then averaged
#   across all themes for a given group (e.g. location-year or time period),
#   providing a summary measure of theme sequence predictability.
#   Higher values = more variable theme ordering.
#   Lower values  = more stereotyped, predictable sequences.
#
#   Users can apply this script to any grouping variable in their transition
#   data (e.g. location_year, season, temporal period) by adjusting the
#   group_var parameter in the SETTINGS section below.
#
# Input:
#   A data frame of first-order Markov transition probabilities with columns:
#     - <group_var>  : grouping variable (e.g. "location_year" or "period")
#     - from         : theme transitioning from
#     - to           : theme transitioning to
#     - probability  : transition probability P(to | from)
#
#   This can be produced from raw theme sequence data using the Markov
#   analysis in 02_markov_sequences.R.
#
# Output:
#   - transition_entropy_phrase_level.csv  : per-theme entropy per group
#   - transition_entropy_summary.csv       : mean entropy summary per group
#
# Authors: Edda Elísabet Magnúsdóttir et al.
# Contact: edda@hi.is
# =============================================================================


# -----------------------------------------------------------------------------
# SETTINGS — adjust these for your data
# -----------------------------------------------------------------------------

# Path to your transition probability data (output of 02_markov_sequences.R)
input_file <- "data/transition_probs.csv"

# Grouping variable name — change to match your data
# Examples: "location_year", "period", "season"
group_var <- "location_year"

# Output folder
output_folder <- "output"

# -----------------------------------------------------------------------------
# PACKAGES
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

transition_probs <- read_csv(input_file)

# Verify required columns are present
required_cols <- c(group_var, "from", "to", "probability")
missing_cols <- setdiff(required_cols, colnames(transition_probs))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------------------------------------------------------
# SHANNON ENTROPY FUNCTION
#
# Computes Shannon entropy (in nats, using natural logarithm) for a vector
# of transition probabilities from a single theme.
# Zero probabilities are excluded before calculation.
#
# Formula: H(i) = -sum( P(j|i) * ln( P(j|i) ) )
# where the sum is over all themes j with non-zero transition probability.
# -----------------------------------------------------------------------------

calc_entropy <- function(probs) {
  probs <- probs[probs > 0]
  -sum(probs * log(probs))
}

# -----------------------------------------------------------------------------
# STEP 1: Per-theme entropy per group
# -----------------------------------------------------------------------------

entropy_phrase_level <- transition_probs %>%
  group_by(across(all_of(c(group_var, "from")))) %>%
  summarise(
    entropy          = calc_entropy(probability),
    max_prob         = max(probability),
    dominant_theme   = to[which.max(probability)],
    n_transitions    = n(),
    .groups          = "drop"
  )

cat("\n==== Per-theme entropy (first rows) ====\n")
print(head(entropy_phrase_level, 20))

# -----------------------------------------------------------------------------
# STEP 2: Summary per group
#
# Mean Shannon entropy averaged equally across all themes in each group.
# Also reports SD, median, and mean maximum transition probability as a
# complementary measure of sequence stereotypy.
# -----------------------------------------------------------------------------

entropy_summary <- entropy_phrase_level %>%
  group_by(across(all_of(group_var))) %>%
  summarise(
    n_themes       = n(),
    mean_entropy   = round(mean(entropy), 3),
    sd_entropy     = round(sd(entropy), 3),
    median_entropy = round(median(entropy), 3),
    mean_max_prob  = round(mean(max_prob), 3),
    .groups        = "drop"
  )

cat("\n==== Shannon entropy summary by", group_var, "====\n")
print(entropy_summary)

# -----------------------------------------------------------------------------
# STEP 3: Save outputs
# -----------------------------------------------------------------------------

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

write_csv(
  entropy_phrase_level,
  file.path(output_folder, "transition_entropy_phrase_level.csv")
)

write_csv(
  entropy_summary,
  file.path(output_folder, "transition_entropy_summary.csv")
)

cat("\n==== Outputs saved to", output_folder, "====\n")
