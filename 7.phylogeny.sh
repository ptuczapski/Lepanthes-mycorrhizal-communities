#!/bin/bash

#SBATCH --job-name=filter_tree
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
# Script 7. Phylogenetic reconstruction of fungal ASVs
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Construct a phylogenetic tree of fungal amplicon sequence
# variants (ASVs) using QIIME 2 and retain the ASVs associated
# with the four focal Lepanthes species.
#
# The resulting phylogeny is exported in Newick format for
# downstream analyses of fungal phylogenetic diversity in R.
#
# WORKFLOW:
# 1. Align representative ASV nucleotide sequences using MAFFT.
# 2. Mask alignment positions using QIIME 2.
# 3. Construct an approximately maximum-likelihood tree
#    using FastTree.
# 4. Generate a rooted phylogenetic tree.
# 5. Filter the tree using the orchid-associated ASV table.
# 6. Export the filtered tree in Newick format.
#
# INPUT FILES:
# - rep-seqs_nochim.qza
#     Representative fungal ASV nucleotide sequences
#     imported into QIIME 2 in Script 3.
#
# - orchidsp-filtered_table_nochim.qza
#     Feature table containing fungal ASVs associated
#     with the four focal Lepanthes species.
#
# OUTPUT FILES:
# - aligned-rep-seqs.qza
# - masked-aligned-rep-seqs.qza
# - unrooted-tree.qza
# - rooted-tree.qza
# - tree_orchidsonly.qza
# - exported-tree/tree.nwk

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
# 2. Align fungal ASV sequences and construct a phylogeny
# ------------------------------------------------------------

# Perform multiple sequence alignment of representative
# fungal ASV nucleotide sequences using MAFFT.
#
# Mask highly variable and poorly conserved alignment
# positions using the default QIIME 2 masking procedure.
#
# Construct an approximately maximum-likelihood phylogeny
# from the masked alignment using FastTree.
#
# Generate a rooted tree using midpoint rooting.
#
# The resulting tree contains the representative ASVs
# provided in rep-seqs_nochim.qza.


qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_nochim.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
# ------------------------------------------------------------
# 3. Filter the phylogenetic tree by orchid-associated ASVs
# ------------------------------------------------------------

# Retain only fungal ASVs represented in the feature table
# containing samples from the four focal Lepanthes species.
#
# This step removes tree tips corresponding to ASVs absent
# from the selected orchid-associated fungal dataset.
#
# The resulting phylogeny is used for subsequent analyses
# of fungal phylogenetic diversity.  
qiime phylogeny filter-tree \
  --i-tree rooted-tree.qza \
  --i-table orchidspfilteredtablenochim.qza \
  --o-filtered-tree tree_orchidsonly.qza
  
# ------------------------------------------------------------
# 4. Export the filtered phylogenetic tree
# ------------------------------------------------------------

# Export the filtered phylogeny from QIIME 2.
#
# The resulting Newick file can be imported into R
# using packages such as ape or phyloseq.
#
# The exported file is expected to be named tree.nwk
# within the specified output directory.
qiime tools export \
  --input-path tree_orchidsonly.qza \
  --output-path exported-tree
  
# ============================================================
# End of Script 7
# ============================================================  