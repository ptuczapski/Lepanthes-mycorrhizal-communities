# ============================================================
# Script 11. Fungal ASV accumulation curves
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Evaluate fungal ASV accumulation across sampled orchid
# individuals for the General Fungi (GF) or Confirmed
# Mycorrhizal Fungi (CMF) dataset.
#
# ANALYSES:
# 1. Calculate an overall ASV accumulation curve.
# 2. Calculate separate accumulation curves for each of
#    the four focal Lepanthes species.
# 3. Estimate mean cumulative ASV richness and empirical
#    95% sample-order envelopes using 1,000 randomizations.
# 4. Generate a PDF figure of the species-specific curves.
#
# INPUTS:
# - ../allorch2otu.txt
#     GF ASV abundance table.
#     ASVs are rows and orchid samples are columns.
#
# - ../qiime2metadata.tsv
#     Sample metadata containing the orchid.sp. column.
#
# FOR CMF ANALYSIS:
# Replace the GF abundance table with the corresponding
# CMF table, ../allorch2otu_myco.txt.
#
# Load libraries after installation.
library(vegan)
library(ggplot2)
library(tidyverse)

# ------------------------------------------------------------
# 1. Import the fungal ASV abundance table
# ------------------------------------------------------------

# Import the GF abundance table.
#
# The original file contains ASVs as rows and individual
# orchid samples as columns.
#
# For CMF analyses, replace the filename with:
# "../allorch2otu_myco.txt"

otu <- read.table(file = "../allorch2otu.txt", sep = "\t", header = T, row.names = 1, 
                  comment.char = "")
#substitute ".." with "--" in barcode names
otu <- otu %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))
# Transpose the abundance matrix.
otu <- t(as.matrix(otu))

# ------------------------------------------------------------
# 2. Convert read counts to presence/absence
# ------------------------------------------------------------
# This prevents differences in sequencing abundance from
# affecting the cumulative number of detected ASVs.
otu_pa <- 1L * (otu > 0)

# Inspect the number of samples and ASVs.
dim(otu_pa)

# ------------------------------------------------------------
# 3. Overall fungal ASV accumulation
# ------------------------------------------------------------
# Calculate the cumulative number of unique ASVs detected
# as individual orchid samples are progressively added.
#
# method = "random" randomizes sample order.
# permutations = 1000 generates 1,000 randomized orders.
#
# The resulting object contains mean cumulative richness
# and its standard deviation at each sampling effort.
#
# This initial analysis pools all available orchid samples.
set.seed(123)

accum <- specaccum(
  otu_pa,
  method = "random", #The function randomizes the order in which samples are added. It calculates the mean cumulative richness and variability at each sampling effort.
  permutations = 1000
)

# Inspect the overall accumulation results.
accum

# ============================================================
# SPECIES-SPECIFIC ASV ACCUMULATION
# ============================================================

# ------------------------------------------------------------
# 4. Import sample metadata
# ------------------------------------------------------------

metadata <- read.table(file = "../qiime2metadata.tsv", sep = "\t", header = T, row.names = 1, quote = "")

# ------------------------------------------------------------
# 5. Match abundance data and metadata
# ------------------------------------------------------------
# Identify samples represented in both datasets.

common_samples <- intersect(
  rownames(otu),
  rownames(metadata)
)

# Retain only shared samples and arrange them in the
# same order in both objects.
otu <- otu[common_samples, , drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Recalculate presence/absence after sample matching.
otu_pa <- 1L * (otu > 0)

# Check sample alignment
stopifnot(
  identical(rownames(otu_pa), rownames(metadata))
)

# ------------------------------------------------------------
# 6. Calculate species-specific accumulation curves
# ------------------------------------------------------------
# Calculate a separate ASV accumulation curve for each
# orchid species.
#
# ASVs absent from all samples of a particular species
# are removed from that species-specific matrix.
#
# Each curve is based on 1,000 randomized sample orders.
set.seed(123)

species_names <- unique(metadata$orchid.sp.)

accum_list <- lapply(species_names, function(sp) {
  
  # Select samples belonging to focal species
  idx <- metadata$orchid.sp. == sp
  
  # Remove ASVs absent from all samples in this group
  group_otu <- otu_pa[idx, , drop = FALSE]
  group_otu <- group_otu[, colSums(group_otu) > 0, drop = FALSE]
  
  # Calculate the randomized accumulation curve.
  acc <- specaccum(
    group_otu,
    method = "random",
    permutations = 1000
  )
  # Store cumulative richness and standard deviation.
  data.frame(
    samples = acc$sites,
    richness = acc$richness,
    sd = acc$sd,
    species = sp
  )
})
# Combine the results for all four orchid species.
accum_species <- do.call(rbind, accum_list)

# ------------------------------------------------------------
# 7. Plot species-specific accumulation curves
# ------------------------------------------------------------
# This exploratory figure displays mean cumulative ASV
# richness without uncertainty envelopes.
ggplot(
  accum_species,
  aes(x = samples, y = richness, color = species)
) +
  geom_line(linewidth = 1) +
  labs(
    x = "Number of samples",
    y = "Cumulative number of ASVs",
    color = "Orchid species",
    title = "Fungal ASV accumulation by orchid species"
  ) +
  theme_classic(base_size = 14)


# ============================================================
# ACCUMULATION CURVES WITH 95% SAMPLE-ORDER ENVELOPES
# ============================================================
#Confirm sample alignment
metadata <- metadata[
  rownames(otu_pa),
  ,
  drop = FALSE
]

# Inspect dimensions.
# Original GF dataset:
# 387 samples x 6,543 ASVs
# Metadata: 387 samples x 13 variables
#
# These dimensions are dataset-specific and may differ
# when the CMF dataset is analyzed.

dim(otu_pa)       # 387 x 6543
dim(metadata)     # 387 x 13
# Confirm exact sample alignment.
stopifnot(
  identical(rownames(otu_pa), rownames(metadata))
)

# ------------------------------------------------------------
# 8. Calculate accumulation curves with uncertainty
# ------------------------------------------------------------
# Randomly permute the order of orchid samples 1,000 times.
#
# For each permutation:
# - Add samples sequentially.
# - Track the ASVs detected up to each sampling effort.
# - Count the cumulative number of unique ASVs.
#
# Summarize cumulative richness across permutations
# using the mean and empirical 2.5th and 97.5th
# percentiles.
#
# These percentile bounds describe variability among
# randomized sample orders. 

set.seed(123)

B <- 1000  # Number of randomized sample orders
species_names <- unique(metadata$orchid.sp.)

accum_list <- lapply(species_names, function(sp) {
  
  # Select samples belonging to this species
  idx <- metadata$orchid.sp. == sp
  
  group_otu <- otu_pa[idx, , drop = FALSE]
  
  # Remove ASVs absent from all samples in this group
  group_otu <- group_otu[
    , colSums(group_otu) > 0,
    drop = FALSE
  ]
  # Number of orchid samples for this species.
  n <- nrow(group_otu)
  
  # Store cumulative richness for each permutation
  richness_perm <- matrix(
    NA_real_,
    nrow = B,
    ncol = n
  )
  
  for (b in seq_len(B)) {
    
    # Randomize sample order
    perm <- sample.int(n)
    
    # Count ASVs detected across all samples added so far
    cumulative_counts <- integer(ncol(group_otu))
  
  # Calculate cumulative richness correctly for each permutation
  
    for (k in seq_len(n)) {
      cumulative_counts <- cumulative_counts +
        group_otu[perm[k], ]
      
      richness_perm[b, k] <- sum(cumulative_counts > 0)
    }
  }
  
  # Calculate the mean cumulative richness and
  # empirical 95% sample-order envelope at each
  # sampling effort.
  data.frame(
    samples = seq_len(n),
    richness = colMeans(richness_perm),
    lower = apply(
      richness_perm, 2, quantile, probs = 0.025
    ),
    upper = apply(
      richness_perm, 2, quantile, probs = 0.975
    ),
    species = sp
  )
})

# ------------------------------------------------------------
# 9. Combine species-specific accumulation results
# ------------------------------------------------------------

# Combine results into one data frame for plotting.
#
# Each row contains:
# - Number of accumulated samples
# - Mean cumulative ASV richness
# - Lower 2.5th percentile
# - Upper 97.5th percentile
# - Orchid species

accum_species <- do.call(rbind, accum_list)

# ------------------------------------------------------------
# 10. Generate and save accumulation figure
# ------------------------------------------------------------

# For the CMF dataset, change the output filename to:
# "ASV_accumulation_CMF.pdf"
pdf(file = tf <- tempfile(fileext = ".pdf"),width=10,height=8)
ggplot(
  accum_species,
  aes(
    x = samples,
    y = richness,
    color = species,
    fill = species
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  labs(
    x = "Number of samples",
    y = "Cumulative number of ASVs",
    color = "Orchid species",
    fill = "Orchid species",
    title = "Fungal ASV accumulation by orchid species",
    subtitle = "Shaded areas represent 95% sample-order envelopes"
  ) +
  theme_classic(base_size = 14)
dev.off()
shell.exec(tf)

# ============================================================
# End of Script 11
# ============================================================
