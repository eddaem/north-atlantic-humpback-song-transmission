# =============================================================================
# 02_markov_sequences.R
#
# First-order Markov chain analysis of humpback whale song theme sequences
#
# Part of the analysis reported in:
#   Magnúsdóttir et al. (in press). Subarctic feeding grounds as cultural
#   exchange hubs in the North Atlantic humpback whale song network.
#   
#
# Description:
#   Reconstructs representative song sequences for each location-year group
#   using first-order Markov chain analysis. Theme transitions are extracted
#   from annotated song recordings, transition counts and probabilities are
#   computed, and the most likely song sequences are reconstructed by
#   following the highest-probability transitions from each theme.
#   Also computes set median theme strings per theme per location-year,
#   which serve as input to the LSI analysis (03_LSI_similarity_dendrograms.R).
#
# Input:
#   A CSV file of annotated theme sequences per recording, with columns:
#     - location_year : location-year identifier (e.g. ICE11, CAR09)
#     - Begin File    : recording file name
#     - Singer        : singer identifier within recording
#     - Phrase        : theme label (e.g. 4b, 13c)
#   (Exported from Raven Pro annotation files)
#
# Output:
#   - transition_counts.csv        : pairwise theme transition counts
#   - transition_probabilities.csv : pairwise theme transition probabilities
#   - most_likely_sequences.csv    : reconstructed song sequences per LocYY
#   - set_median_theme_strings.tsv : set median unit strings per theme/LocYY
#                                    (input to 03_LSI_similarity_dendrograms.R)
#
# Authors: Edda Elísabet Magnúsdóttir et al.
# Contact: edda@hi.is
# =============================================================================


# -----------------------------------------------------------------------------
# SETTINGS — adjust these for your data
# -----------------------------------------------------------------------------

# Path to combined annotated theme sequence data
input_file    <- "data/phrases_combined.csv"

# Output folder
output_folder <- "output"

# Maximum steps for sequence reconstruction
max_steps     <- 15

# -----------------------------------------------------------------------------
# PACKAGES
# -----------------------------------------------------------------------------

library(tidyverse)
library(stringdist)
library(readr)

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# -----------------------------------------------------------------------------
# STEP 1: Load data
# -----------------------------------------------------------------------------

all_data <- read_csv(input_file, show_col_types = FALSE)
loc_years <- sort(unique(all_data$location_year))
cat("Location-years found:", paste(loc_years, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# STEP 2: Extract theme transitions
#
# Transitions are extracted per recording file and singer. Consecutive
# repeated themes are collapsed (a theme is only counted once per run).
# Explicit transition notations (e.g. "4b->5c") in the Phrase column
# are expanded into separate themes before extraction.
# -----------------------------------------------------------------------------

collapse_repeats <- function(phrases) {
  phrases[c(TRUE, phrases[-1] != phrases[-length(phrases)])]
}

extract_transitions <- function(df_group) {
  transitions <- list()
  all_phrases <- c()

  for (i in seq_len(nrow(df_group))) {
    phrase <- df_group$Phrase[i]
    if (str_detect(phrase, "->|-->|-")) {
      parts <- str_split(phrase, "->|-->|-")[[1]] %>% str_trim()
      parts <- parts[parts != ""]
      all_phrases <- c(all_phrases, parts)
    } else {
      all_phrases <- c(all_phrases, phrase)
    }
  }

  all_phrases <- collapse_repeats(all_phrases)

  if (length(all_phrases) >= 2) {
    for (i in seq_len(length(all_phrases) - 1)) {
      transitions <- append(transitions, list(c(all_phrases[i], all_phrases[i + 1])))
    }
  }
  transitions
}

all_transitions <- all_data %>%
  group_by(location_year, `Begin File`, Singer) %>%
  group_split() %>%
  map(~ {
    loc_year    <- .x$location_year[1]
    transitions <- extract_transitions(.x)
    map(transitions, ~ c(location_year = loc_year, from = .x[1], to = .x[2]))
  }) %>%
  flatten() %>%
  map_dfr(~ as_tibble(as.list(.x)))

cat("Total transitions extracted:", nrow(all_transitions), "\n")

# -----------------------------------------------------------------------------
# STEP 3: Transition counts and probabilities
# -----------------------------------------------------------------------------

transition_counts <- all_transitions %>%
  group_by(location_year, from, to) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(location_year, desc(count))

transition_probs <- transition_counts %>%
  group_by(location_year, from) %>%
  mutate(probability = count / sum(count)) %>%
  ungroup()

write_csv(transition_counts, file.path(output_folder, "transition_counts.csv"))
write_csv(transition_probs,  file.path(output_folder, "transition_probabilities.csv"))
cat("Transition counts and probabilities saved\n")

# -----------------------------------------------------------------------------
# STEP 4: Reconstruct most likely song sequences
# -----------------------------------------------------------------------------

reconstruct_sequence <- function(probs_loc, rank = 1, max_steps = 15) {
  ranked_transitions <- probs_loc %>%
    group_by(from) %>%
    slice_max(probability, n = rank, with_ties = FALSE) %>%
    slice_tail(n = 1) %>%
    ungroup()

  start_phrase <- probs_loc %>%
    group_by(from) %>%
    summarise(total = sum(count)) %>%
    slice_max(total, n = 1) %>%
    pull(from)

  sequence <- start_phrase
  visited  <- c(start_phrase)

  for (step in 1:max_steps) {
    next_phrase <- ranked_transitions %>%
      filter(from == tail(sequence, 1)) %>%
      pull(to)
    if (length(next_phrase) == 0 || next_phrase %in% visited) break
    sequence <- c(sequence, next_phrase)
    visited  <- c(visited, next_phrase)
  }
  sequence
}

seq_output <- map_dfr(loc_years, ~ {
  probs_loc <- transition_probs %>% filter(location_year == .x)
  seq1 <- reconstruct_sequence(probs_loc, rank = 1, max_steps)
  seq2 <- reconstruct_sequence(probs_loc, rank = 2, max_steps)
  cat("\n", .x, "\n  Most likely:   ", paste(seq1, collapse = " -> "),
      "\n  Second likely: ", paste(seq2, collapse = " -> "), "\n")
  tibble(
    location_year = .x,
    rank          = c(1, 2),
    sequence      = c(paste(seq1, collapse = " -> "),
                      paste(seq2, collapse = " -> "))
  )
})

write_csv(seq_output, file.path(output_folder, "most_likely_sequences.csv"))
cat("\nMost likely sequences saved\n")

# -----------------------------------------------------------------------------
# STEP 5: Set median theme strings per theme per location-year
#
# For each theme type in each location-year, identifies the most
# representative instance (medoid) as the theme instance with the
# highest summed LSI to all other instances of the same theme type.
# Used as input to LSI similarity analysis (03_LSI_similarity_dendrograms.R).
# -----------------------------------------------------------------------------

# Load phrase-level unit sequence data
# (phrases data frame with columns: PhraseID, LocYY, Phrase, UnitSeqString)
phrases <- read_tsv("data/phrases.tsv")

get_medoid_idx <- function(strings) {
  strings <- as.character(strings)
  if (length(strings) == 1) return(1L)
  D     <- stringdistmatrix(strings, strings, method = "lv")
  LSI   <- 1 - (D / max(D))
  diag(LSI) <- 0
  which.max(rowSums(LSI))
}

median_theme_strings <- phrases %>%
  group_by(LocYY, Phrase) %>%
  reframe(
    n_instances     = n(),
    medoid_idx      = get_medoid_idx(UnitSeqString),
    SetMedianString = UnitSeqString[medoid_idx]
  )

write_tsv(median_theme_strings,
          file.path(output_folder, "set_median_theme_strings.tsv"))

cat("Set median theme strings saved —", nrow(median_theme_strings), "rows\n")
cat("(Use as input to 03_LSI_similarity_dendrograms.R)\n")
