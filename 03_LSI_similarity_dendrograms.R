# =============================================================================
# 03_LSI_similarity_dendrograms.R
#
# Levenshtein Similarity Index (LSI) and pvclust dendrogram analysis
#
# Part of the analysis reported in:
#   Magnúsdóttir et al. (in press). Subarctic feeding grounds as cultural
#   exchange hubs in the North Atlantic humpback whale song network.
#   Royal Society Open Science.
#
# Description:
#   Computes pairwise unweighted LSI between set median theme sequences
#   across location-year groups. All edit operations (insertion, deletion,
#   substitution) are assigned equal cost of 1. LSI is computed as:
#
#     LSI(A, B) = 1 - (Levenshtein distance / max sequence length)
#
#   The unweighted LSI was selected as the primary approach on the basis of
#   best dendrogram fit (CCC = 0.887) and clearest within/between theme
#   separation (within = 0.38, between = 0.885).
#
#   Bootstrapped average-linkage hierarchical clustering (pvclust, 1000
#   resamples) is used to assess cluster support. AU (approximately
#   unbiased) p-values >= 0.95 and BP (bootstrap probability) values >= 0.70
#   are indicated on the dendrogram. The Cophenetic Correlation Coefficient
#   (CCC) is reported as a measure of dendrogram fit; values above 0.8
#   indicate a reliable representation of the underlying similarity matrix.
#
# Input:
#   - set_median_theme_strings.tsv: set median unit strings per theme per
#     location-year. Output of 02_markov_sequences.R.
#   - phrases.tsv: theme sequence data for within/between validation.
#
# Output:
#   - LSI_unweighted_matrix.csv        : pairwise LSI matrix
#   - LSI_unweighted_pvclust.tiff      : pvclust dendrogram
#   - within_between_summary.csv       : validation summary
#
# Authors: Edda Elísabet Magnúsdóttir et al.
# Contact: edda@hi.is
# =============================================================================


# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

output_folder <- "output"
n_bootstrap   <- 1000    # pvclust bootstrap resamples

# -----------------------------------------------------------------------------
# PACKAGES
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(pvclust)
library(dendextend)
library(RColorBrewer)

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# -----------------------------------------------------------------------------
# STEP 1: Load set median theme sequences
#
# S1 must have columns: LocYY, Phrase (theme identifier), SetMedianString
# (space-delimited unit string). Output of 02_markov_sequences.R.
# -----------------------------------------------------------------------------

S1 <- read_tsv("data/set_median_theme_strings.tsv")

S1_units <- S1 %>%
  filter(!is.na(SetMedianString), SetMedianString != "") %>%
  mutate(
    Label   = paste(LocYY, Phrase, sep = " | "),
    UnitSeq = strsplit(SetMedianString, "\\s+")
  )

set_units        <- S1_units$UnitSeq
names(set_units) <- S1_units$Label
n_set            <- length(set_units)

cat("Set median sequences loaded:", n_set, "\n")
cat("Labels (first 5):\n")
print(head(names(set_units), 5))

# -----------------------------------------------------------------------------
# STEP 2: Unweighted Levenshtein distance function
#
# All edit operations (insertion, deletion, substitution) cost 1.
# Normalised by the length of the longer sequence to give LSI in [0, 1].
# -----------------------------------------------------------------------------

lv_unweighted <- function(seq1, seq2) {
  n <- length(seq1); m <- length(seq2)
  if (n == 0 && m == 0) return(0)
  if (n == 0) return(m)
  if (m == 0) return(n)
  D <- matrix(0, n + 1, m + 1)
  for (i in 1:n) D[i + 1, 1] <- i
  for (j in 1:m) D[1, j + 1] <- j
  for (i in 1:n) for (j in 1:m) {
    D[i + 1, j + 1] <- min(
      D[i,     j + 1] + 1,
      D[i + 1, j    ] + 1,
      D[i,     j    ] + if (seq1[i] == seq2[j]) 0 else 1
    )
  }
  D[n + 1, m + 1]
}

# -----------------------------------------------------------------------------
# STEP 3: Compute pairwise unweighted LSI matrix
# -----------------------------------------------------------------------------

seq_lengths <- vapply(set_units, length, integer(1))
Lmax        <- outer(seq_lengths, seq_lengths, pmax)

D_unw <- matrix(0, n_set, n_set,
                dimnames = list(names(set_units), names(set_units)))

cat("Computing pairwise LSI matrix — this may take a few minutes...\n")

for (i in 1:n_set) {
  for (j in i:n_set) {
    if (i == j) next
    d <- lv_unweighted(set_units[[i]], set_units[[j]])
    D_unw[i, j] <- D_unw[j, i] <- d
  }
}

LSI_unw <- 1 - (D_unw / Lmax)

# CCC
hc  <- hclust(as.dist(1 - LSI_unw), method = "average")
CCC <- cor(as.dist(1 - LSI_unw), cophenetic(hc))
cat("CCC:", round(CCC, 3), "\n")

# Save matrix
write_csv(as.data.frame(LSI_unw),
          file.path(output_folder, "LSI_unweighted_matrix.csv"))
cat("LSI matrix saved\n")

# -----------------------------------------------------------------------------
# STEP 4: Within/between theme dissimilarity (validation)
#
# Confirms that mean dissimilarity between theme types substantially
# exceeds that within theme types, validating the classification framework.
# -----------------------------------------------------------------------------

phrases <- read_tsv("data/phrases.tsv")

within_between_summary <- function(lsi_mat) {
  dsi_mat <- 1 - lsi_mat
  as.data.frame(as.table(dsi_mat)) %>%
    rename(id1 = Var1, id2 = Var2, d = Freq) %>%
    mutate(id1 = as.character(id1), id2 = as.character(id2)) %>%
    filter(id1 < id2) %>%
    mutate(
      theme1 = sub(".* \\| ", "", id1),
      theme2 = sub(".* \\| ", "", id2),
      type   = ifelse(theme1 == theme2, "within", "between")
    ) %>%
    group_by(type) %>%
    summarise(
      mean_d = round(mean(d), 3),
      sd_d   = round(sd(d),   3),
      n      = n(),
      .groups = "drop"
    )
}

wb <- within_between_summary(LSI_unw)
cat("\n==== Within vs between theme dissimilarity ====\n")
print(wb)
write_csv(wb, file.path(output_folder, "within_between_summary.csv"))

# -----------------------------------------------------------------------------
# STEP 5: pvclust bootstrapped dendrogram
# -----------------------------------------------------------------------------

build_theme_colours <- function(labels) {
  theme_types  <- sub(".* \\| ", "", labels)
  unique_types <- sort(unique(theme_types))
  n_types      <- length(unique_types)
  pal <- if (n_types <= 12) brewer.pal(max(n_types, 3), "Paired") else
    colorRampPalette(brewer.pal(12, "Paired"))(n_types)
  colour_map <- setNames(pal[seq_len(n_types)], unique_types)
  colour_map[theme_types]
}

cat("\nRunning pvclust (", n_bootstrap, "bootstraps)...\n")

pv <- pvclust(as.matrix(LSI_unw),
              method.hclust = "average",
              method.dist   = function(x) as.dist(1 - cor(x)),
              nboot         = n_bootstrap,
              quiet         = TRUE)

row_labels        <- rownames(LSI_unw)
phrase_colour_map <- build_theme_colours(row_labels)

dend         <- as.dendrogram(pv$hclust)
leaf_order   <- order.dendrogram(dend)
labels(dend) <- row_labels[leaf_order]
dend <- dendextend::set(dend, "branches_lwd", 2.0)

dend <- dendrapply(dend, function(node) {
  if (is.leaf(node)) {
    lbl        <- attr(node, "label")
    theme_type <- sub(".* \\| ", "", lbl)
    col_idx    <- which(sub(".* \\| ", "", names(phrase_colour_map)) == theme_type)
    col        <- if (length(col_idx) > 0) phrase_colour_map[col_idx[1]] else "grey30"
    attr(node, "nodePar") <- list(lab.font = 2, lab.cex = 1.4,
                                   lab.col = col, pch = NA)
  }
  node
})

au_vals <- if (!is.null(pv$edges)) pv$edges$au else pv$au
bp_vals <- if (!is.null(pv$edges)) pv$edges$bp else pv$bp

fname <- file.path(output_folder, "LSI_unweighted_pvclust.tiff")

tiff(fname, width = 20, height = 10, units = "in",
     pointsize = 10, bg = "white", res = 300)
par(mar = c(10, 6, 2, 2))
plot(dend, ylab = "Height (dissimilarity)", xlab = "", sub = "", main = "",
     cex.lab = 1.4, cex.axis = 1.4)

xy       <- dendextend::get_nodes_xy(dend)
xy_df    <- as.data.frame(xy); colnames(xy_df) <- c("x", "y")
internal <- xy_df[xy_df$y > 0, ][order(xy_df[xy_df$y > 0, "y"]), ]
m        <- min(length(au_vals), nrow(internal))
internal <- internal[seq_len(m), ]
y_off    <- 0.012 * max(internal$y, na.rm = TRUE)

idx_bp <- which(bp_vals[seq_len(m)] >= 0.70)
if (length(idx_bp) > 0)
  points(internal$x[idx_bp] - 0.3, internal$y[idx_bp] - y_off,
         pch = 16, col = "green3", cex = 1.6)

idx_au <- which(au_vals[seq_len(m)] >= 0.95)
if (length(idx_au) > 0)
  points(internal$x[idx_au] + 0.3, internal$y[idx_au] + y_off,
         pch = 16, col = "red", cex = 1.6)

legend("topright",
       legend = c("BP \u2265 70%", "AU \u2265 95%",
                  paste0("CCC = ", round(CCC, 3))),
       col = c("green3", "red", NA), pch = c(16, 16, NA),
       pt.cex = 1.6, cex = 1.3, bty = "n")
dev.off()

cat("Dendrogram saved:", fname, "\n")
