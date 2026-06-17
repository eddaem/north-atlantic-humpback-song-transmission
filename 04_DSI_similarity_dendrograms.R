# =============================================================================
# 04_DSI_similarity_dendrograms.R
#
# Dice's Similarity Index (DSI) and pvclust dendrogram analysis
#
# Part of the analysis reported in:
#   Magnúsdóttir et al. (in press). Subarctic feeding grounds as cultural
#   exchange hubs in the North Atlantic humpback whale song network.
#   Royal Society Open Science.
#
# Description:
#   Computes pairwise Dice's Similarity Index (DSI) between location-year
#   groups based on theme presence/absence. DSI quantifies shared theme
#   repertoire independently of sequence order:
#     DSI(A, B) = 2 * |A ∩ B| / (|A| + |B|)
#   where A and B are the sets of themes present in each location-year.
#   Bootstrapped average-linkage hierarchical clustering (pvclust, 1000
#   resamples) is used to assess cluster support (AU and BP values).
#
# Input:
#   - phrases: data frame of theme sequences with columns:
#       LocYY  : location-year identifier
#       Phrase : theme identifier
#     (Output of 02_markov_sequences.R or provided as data/phrases.tsv)
#
# Output:
#   - DSI matrix (printed and saved)
#   - pvclust dendrogram saved as TIFF
#
# Authors: Edda Elísabet Magnúsdóttir et al.
# Contact: edda@hi.is
# =============================================================================


# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

output_folder <- "output"
n_bootstrap   <- 1000

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
# STEP 1: Load theme presence data
# -----------------------------------------------------------------------------

phrases <- read_tsv("data/phrases.tsv")

phrase_presence <- phrases %>%
  filter(!is.na(LocYY), !is.na(Phrase), Phrase != "") %>%
  distinct(LocYY, Phrase)

cat("Location-years:", length(unique(phrase_presence$LocYY)), "\n")
cat("Unique themes: ", length(unique(phrase_presence$Phrase)), "\n")

# -----------------------------------------------------------------------------
# STEP 2: Compute pairwise DSI matrix
#
# DSI(A, B) = 2 * |shared themes| / (|themes in A| + |themes in B|)
# Values range from 0 (no shared themes) to 1 (identical repertoires).
# -----------------------------------------------------------------------------

dice_pair <- function(A, B) {
  A <- unique(A); B <- unique(B)
  shared <- length(intersect(A, B))
  if ((length(A) + length(B)) == 0) return(NA_real_)
  2 * shared / (length(A) + length(B))
}

sets     <- sort(unique(phrase_presence$LocYY))
dice_mat <- matrix(NA_real_, length(sets), length(sets),
                   dimnames = list(sets, sets))

for (i in seq_along(sets)) {
  for (j in seq_along(sets)) {
    A <- phrase_presence %>% filter(LocYY == sets[i]) %>% pull(Phrase)
    B <- phrase_presence %>% filter(LocYY == sets[j]) %>% pull(Phrase)
    dice_mat[i, j] <- dice_pair(A, B)
  }
}

cat("\n==== DSI matrix ====\n")
print(round(dice_mat, 3))

write_csv(as.data.frame(dice_mat),
          file.path(output_folder, "DSI_matrix.csv"))

# -----------------------------------------------------------------------------
# STEP 3: pvclust dendrogram on DSI matrix
# -----------------------------------------------------------------------------

dist_mat <- as.dist(1 - dice_mat)
hc       <- hclust(dist_mat, method = "average")
CCC      <- cor(dist_mat, cophenetic(hc))
cat("\nCCC:", round(CCC, 3), "\n")

cat("Running pvclust (", n_bootstrap, "bootstraps)...\n")
pv <- pvclust(dice_mat,
              method.hclust = "average",
              method.dist   = function(x) as.dist(1 - cor(x)),
              nboot         = n_bootstrap,
              quiet         = TRUE)

dend         <- as.dendrogram(pv$hclust)
leaf_order   <- order.dendrogram(dend)
labels(dend) <- sets[leaf_order]
dend <- dendextend::set(dend, "branches_lwd", 2.0)

au_vals <- if (!is.null(pv$edges)) pv$edges$au else pv$au
bp_vals <- if (!is.null(pv$edges)) pv$edges$bp else pv$bp

fname <- file.path(output_folder, "DSI_pvclust_dendrogram.tiff")

tiff(fname, width = 14, height = 10, units = "in",
     pointsize = 10, bg = "white", res = 300)
par(mar = c(8, 6, 2, 2))
plot(dend,
     ylab     = "Height (dissimilarity)",
     xlab     = "", sub = "", main = "",
     cex.lab  = 1.4, cex.axis = 1.4)

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

cat("Saved:", fname, "\n")
