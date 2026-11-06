#!/bin/bash

#SBATCH --job-name=diversity
#SBATCH --partition=highmem_p
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=1:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ============================================================
# Script 8. Alpha and beta diversity analyses in QIIME 2
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Evaluate fungal community diversity among the four focal
# Lepanthes species using phylogenetic and non-phylogenetic
# diversity metrics.
#
# ANALYSES:
# 1. Generate alpha rarefaction curves to assess sequencing
#    depth and sampling completeness.
# 2. Calculate alpha diversity, including Faith's
#    phylogenetic diversity (Faith's PD).
# 3. Calculate beta diversity using Jaccard dissimilarity
#    and unweighted UniFrac distances.
# 4. Test for differences in fungal community composition
#    among orchid species using PERMANOVA.
# 5. Test for differences in Faith's PD among orchid spp. using Kruskal-Wallis test
#
# DATASETS:
# GF  = General Fungi
# CMF = Confirmed Mycorrhizal Fungi
#
# The commands below apply to the GF dataset.
# The original workflow repeated these analyses for the
# CMF dataset using a sampling depth of 445 reads.
#
# INPUT FILES:
# - outallorch2otu.qza
#     GF ASV abundance table generated following the
#     species-specific singleton filtering in Scripts 4-5.
#
# - tree_orchidsonly.qza
#     Fungal ASV phylogenetic tree generated in Script 7.
#
# - qiime2metadata.tsv
#     Sample metadata containing the "orchid sp." column.
#
# OUTPUTS:
# - alpha-rarefaction4000.qzv
# - core-metrics-results/
#     Alpha and beta diversity artifacts and visualizations.
# - core-metrics-results/jaccard-orchid-sp.-significance.qzv
# - core-metrics-results/
#     unweighted-unifrac-orchid-sp.-significance.qzv
#
# SAMPLING DEPTHS:
# GF:  679 reads
# CMF: 445 reads
#
# IMPORTANT:
# The original script does not include the commands used
# to import the final GF and CMF abundance tables into
# QIIME 2. These are required before executing this script.
# ============================================================
# ------------------------------------------------------------
# 1. Initialize the computing environment
# ------------------------------------------------------------

# Change to the directory from which the SLURM job
# was submitted.

cd $SLURM_SUBMIT_DIR

# Load the QIIME 2 environment used for filtering.
module load QIIME2/2026.04

# Configure Matplotlib for non-interactive execution.
echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc
# Set locale for UTF-8 character encoding.
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

#load software
module load QIIME2/2026.04
echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

# ------------------------------------------------------------
# 2. Generate alpha rarefaction curves
# ------------------------------------------------------------

# Generate alpha rarefaction curves to evaluate how fungal
# ASV richness and phylogenetic diversity change with
# increasing sequencing depth.
#
# Rarefaction is performed across multiple sequencing
# depths up to a maximum of 4,000 reads per sample.
# The original QIIME 2 implementation uses repeated
# subsampling iterations (default:10) to estimate diversity values.
#
# Sample metadata are included to group the rarefaction
# curves by orchid species.
qiime diversity alpha-rarefaction \
  --i-table outallorch2otu.qza \
  --i-phylogeny tree_orchidsonly.qza \
  --p-max-depth 4000 \
  --m-metadata-file qiime2metadata.tsv \
  --o-visualization alpha-rarefaction4000.qzv
# ------------------------------------------------------------
# 3. Calculate alpha and beta diversity metrics
# ------------------------------------------------------------

# Calculate alpha and beta diversity using the QIIME 2
# core-metrics-phylogenetic pipeline. Kruskal-Wallis test is included.
#
# The feature table is rarefied to a sampling depth of
# 679 reads per sample for the GF dataset. 
# For the CMF dataset, the original analysis used
# a sampling depth of 445 reads.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny tree_orchidsonly.qza \
  --i-table outallorch2otu.qza \
  --p-sampling-depth 679 \
  --m-metadata-file qiime2metadata.tsv \
  --output-dir core-metrics-results 
# ------------------------------------------------------------
# 4. Test differences in Jaccard dissimilarity
# ------------------------------------------------------------

# Test whether fungal community composition differs
# among the four focal Lepanthes species using Jaccard
# dissimilarity and PERMANOVA.  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file qiime2metadata.tsv \
  --m-metadata-column 'orchid sp.' \
  --o-visualization core-metrics-results/jaccard-orchid-sp.-significance.qzv \
  --p-pairwise
# ------------------------------------------------------------
# 5. Test differences in unweighted UniFrac distances
# ------------------------------------------------------------

# Test for differences in fungal community composition
# among orchid species using unweighted UniFrac distances
# and PERMANOVA.
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2metadata.tsv \
  --m-metadata-column 'orchid sp.' \
  --o-visualization core-metrics-results/unweighted-unifrac-orchid-sp.-significance.qzv \
  --p-pairwise
  
# ============================================================
# End of Script 8
# ============================================================
