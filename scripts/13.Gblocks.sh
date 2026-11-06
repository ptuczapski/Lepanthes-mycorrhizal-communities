#!/bin/bash

#SBATCH --job-name=Gblocks_sensitivity
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=2:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ============================================================
# Script 13. Gblocks alignment-filtering sensitivity analysis
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Assess the effect of removing poorly aligned positions
# and divergent regions from the fungal ITS alignment
# before phylogenetic reconstruction and UniFrac analysis.
#
# WORKFLOW:
# 1. Align fungal ASV sequences using MAFFT in QIIME 2.
# 2. Export the multiple sequence alignment in FASTA format.
# 3. Filter the alignment using Gblocks v0.91b.
# 4. Import the Gblocks-filtered alignment into QIIME 2.
# 5. Apply QIIME 2 alignment masking.
# 6. Construct and midpoint-root a FastTree phylogeny.
# 7. Calculate diversity metrics using the GF abundance table.
# 8. Test orchid-species differences in unweighted UniFrac
#    distances using QIIME 2 beta-group-significance.
#
# INPUTS:
# - rep-seqs_nochim.qza
#     Representative fungal ASV nucleotide sequences.
#
# - allorch2otu.qza
#     General Fungi (GF) feature table.
#
# - qiime2metadata.tsv
#     Sample metadata containing the "orchid sp." column.
#
# SOFTWARE:
# - QIIME 2 (original script: 2026.04)
# - MAFFT (through QIIME 2 or online tool)
# - Gblocks v0.91b
# - FastTree (through QIIME 2)
#
# REFERENCE:
# Castresana J. (2000).
# Selection of conserved blocks from multiple alignments
# for their use in phylogenetic analysis.
# Molecular Biology and Evolution 17:540-552.
#
# IMPORTANT:
# This is a sensitivity analysis, separate from the
# primary phylogenetic reconstruction in Script 7.
# ============================================================

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
# 2. Align fungal ITS ASV sequences using MAFFT
# ------------------------------------------------------------
# Generate a multiple sequence alignment from the
# representative fungal ASV sequences.
qiime alignment mafft \
  --i-sequences rep-seqs_nochim.qza \
  --o-alignment aligned-rep-seqs.qza \
# ------------------------------------------------------------
# 3. Export the MAFFT alignment
# ------------------------------------------------------------

# Export the aligned nucleotide sequences from QIIME 2
# into a standard FASTA file.
#
# QIIME 2 normally exports the alignment as
# aligned-dna-sequences.fasta.
mkdir fastaMSA
qiime tools export \
  --input-path aligned-rep-seqs.qza \
  --output-path fastaMSA/

# Copy the exported alignment to the filename used
# in the original Gblocks command.
cp \
  fastaMSA/aligned-dna-sequences.fasta \
  outaligneddnasequences.fasta

# ------------------------------------------------------------
# 4. Filter the alignment using Gblocks
# ------------------------------------------------------------
# Load Gblocks version 0.91b.

#load software
module load Gblocks/0.91b
# Remove poorly aligned positions and divergent regions. 
# Gblocks normally creates a filtered FASTA alignment
# with the suffix "-gb".
Gblocks ./outaligneddnasequences.fasta -t=d -b1=7670 -b2=7670 -b3=27886 -b4=2 -b5=a -d=y

# ------------------------------------------------------------
# 5. Import the Gblocks-filtered alignment into QIIME 2
# ------------------------------------------------------------

# Gblocks generates a filtered alignment named:
# outaligneddnasequences.fasta-gb
#
# Copy it to a descriptive FASTA filename.

cp \
  outaligneddnasequences.fasta-gb \
  aligned-dna-sequences-gb.fasta
#load software
module load QIIME2/2026.04
# Import the filtered nucleotide alignment as a
# QIIME 2 FeatureData[AlignedSequence] artifact.
qiime tools import \
  --input-path aligned-dna-sequences-gb.fasta \
  --output-path aligned-sequences-gb.qza \
  --type 'FeatureData[AlignedSequence]'  

# ------------------------------------------------------------
# 6. Apply QIIME 2 alignment masking
# ------------------------------------------------------------
# Mask highly variable alignment positions using
# the default QIIME 2 masking procedure.
qiime alignment mask \
  --i-alignment aligned-sequences-gb.qza \
  --o-masked-alignment masked-aligned-sequences-gb.qza

# ------------------------------------------------------------
# 7. Construct a phylogenetic tree using FastTree
# ------------------------------------------------------------
# Infer an approximately maximum-likelihood phylogeny
# from the masked, Gblocks-filtered alignment.
qiime phylogeny fasttree \
  --i-alignment masked-aligned-sequences-gb.qza \
  --o-tree unrooted-fasttree-gb-tree.qza

# ------------------------------------------------------------
# 8. Midpoint-root the phylogeny
# ------------------------------------------------------------

qiime phylogeny midpoint-root \
  --i-tree unrooted-fasttree-gb-tree.qza \
  --o-rooted-tree rooted-fasttree-gb-tree.qza

# ------------------------------------------------------------
# 9. Calculate phylogenetic diversity metrics
# ------------------------------------------------------------
# Run the QIIME 2 core-metrics-phylogenetic pipeline.
#
# Sampling depth: 679 reads per sample.
#
# This matches the original GF sampling depth.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-fasttree-gb-tree.qza \
  --i-table allorch2otu.qza \
  --p-sampling-depth 679 \
  --m-metadata-file qiime2metadata.tsv \
  --output-dir core-metrics-results  
# ------------------------------------------------------------
# 11. Test orchid-species differences using PERMANOVA
# ------------------------------------------------------------
# Test whether fungal community composition differs
# among the four focal Lepanthes species using
# unweighted UniFrac distances.
#
# The distance matrix is generated from the phylogeny
# constructed using the Gblocks-filtered alignment.
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2metadata.tsv \
  --m-metadata-column 'orchid sp.' \
  --o-visualization core-metrics-results/unweighted-unifrac-orchid-sp.-significance.qzv \
  --p-pairwise
  
# ============================================================
# End of Script 13
# ============================================================
