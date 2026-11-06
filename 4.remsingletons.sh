#!/bin/bash

#SBATCH --job-name=filter-by-orchid
#SBATCH --partition=highmem_p
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=1:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out   #gives output file's name of the job
#SBATCH --error=%x_%j.err

# ============================================================
# Script 4. Species-specific ASV filtering and export
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Generate separate fungal ASV abundance tables for each
# of the four focal Lepanthes species and remove ASVs
# occurring in fewer than two samples within each species.
#
# The resulting tables are exported in TSV format for
# subsequent processing and analysis in R.
#
# INPUT FILES:
# - orchidspfilteredtablenochim.qza
#     QIIME 2 feature table containing fungal ASVs
#     associated with the four focal orchid species.
#
# - qiime2metadata.tsv
#     Sample metadata containing the "orchid sp." column.
#
# OUTPUT FILES:
# - LC2-filtered_table_nochim.tsv
# - LFB2-filtered_table_nochim.tsv
# - LMe2-filtered_table_nochim.tsv
# - LMo2-filtered_table_nochim.tsv
#
# FILTERING CRITERION:
# ASVs must occur in at least two samples within a given
# orchid species to be retained in its species-specific table.
#
# This criterion removes sample-occurrence singletons.
#
# DATASETS:
# The commands below process the General Fungi (GF)
# dataset. The same procedure was applied separately
# to the Confirmed Mycorrhizal Fungi (CMF) dataset.
#
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

# ------------------------------------------------------------
# 2. Lepanthes cribbii
# ------------------------------------------------------------

# Retain samples belonging to Lepanthes cribbii.

qiime feature-table filter-samples \
  --i-table orchidspfilteredtablenochim.qza \
  --m-metadata-file qiime2metadata.tsv \
  --p-where "[orchid sp.] IN ('Lepanthes cribbii')" \
  --o-filtered-table LC-filtered_table_nochim.qza

# Remove ASVs occurring in fewer than two samples
# within L. cribbii.

qiime feature-table filter-features \
  --i-table LC-filtered_table_nochim.qza \
  --p-min-samples 2 \
  --o-filtered-table LC2-filtered_table_nochim.qza

# Export the filtered feature table to BIOM format.

qiime tools export \
  --input-path LC2-filtered_table_nochim.qza \
  --output-path .
# Convert the exported BIOM file to TSV format.

biom convert \
  -i feature-table.biom \
  -o LC2-filtered_table_nochim.tsv \
  --to-tsv

# ------------------------------------------------------------
# 3. Lepanthes falx-bellica
# ------------------------------------------------------------

# Retain samples belonging to Lepanthes falx-bellica.

qiime feature-table filter-samples \
  --i-table orchidspfilteredtablenochim.qza \
  --m-metadata-file qiime2metadata.tsv \
  --p-where "[orchid sp.] IN ('Lepanthes falx-bellica')" \
  --o-filtered-table LFB-filtered_table_nochim.qza

# Remove ASVs occurring in fewer than two samples
# within L. falx-bellica.

qiime feature-table filter-features \
  --i-table LFB-filtered_table_nochim.qza \
  --p-min-samples 2 \
  --o-filtered-table LFB2-filtered_table_nochim.qza
  
# Export the filtered feature table to BIOM format.

qiime tools export \
  --input-path LFB2-filtered_table_nochim.qza \
  --output-path .

# Convert the exported BIOM file to TSV format.

biom convert \
  -i feature-table.biom \
  -o LFB2-filtered_table_nochim.tsv \
  --to-tsv

# ------------------------------------------------------------
# 4. Lepanthes mentosa
# ------------------------------------------------------------

# Retain samples belonging to Lepanthes mentosa.

qiime feature-table filter-samples \
  --i-table orchidspfilteredtablenochim.qza \
  --m-metadata-file qiime2metadata.tsv \
  --p-where "[orchid sp.] IN ('Lepanthes mentosa')" \
  --o-filtered-table LMe-filtered_table_nochim.qza
# Remove ASVs occurring in fewer than two samples
# within L. mentosa.

qiime feature-table filter-features \
  --i-table LMe-filtered_table_nochim.qza \
  --p-min-samples 2 \
  --o-filtered-table LMe2-filtered_table_nochim.qza

# Export the filtered feature table to BIOM format.

qiime tools export \
  --input-path LMe2-filtered_table_nochim.qza \
  --output-path .

# Convert the exported BIOM file to TSV format.

biom convert \
  -i feature-table.biom \
  -o LMe2-filtered_table_nochim.tsv \
  --to-tsv

# ------------------------------------------------------------
# 5. Lepanthes monteverdensis
# ------------------------------------------------------------

# Retain samples belonging to Lepanthes monteverdensis.

qiime feature-table filter-samples \
  --i-table orchidspfilteredtablenochim.qza \
  --m-metadata-file qiime2metadata.tsv \
  --p-where "[orchid sp.] IN ('Lepanthes monteverdensis')" \
  --o-filtered-table LMo-filtered_table_nochim.qza
# Remove ASVs occurring in fewer than two samples
# within L. monteverdensis.

qiime feature-table filter-features \
  --i-table LMo-filtered_table_nochim.qza \
  --p-min-samples 2 \
  --o-filtered-table LMo2-filtered_table_nochim.qza

# Export the filtered feature table to BIOM format.

qiime tools export \
  --input-path LMo2-filtered_table_nochim.qza \
  --output-path .

# Convert the exported BIOM file to TSV format.
biom convert \
  -i feature-table.biom \
  -o LMo2-filtered_table_nochim.tsv \
  --to-tsv
  
# ============================================================
# End of Script 4
# ============================================================
