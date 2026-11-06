# ============================================================
# Script 2. Fungal guild assignment and construction of the
#           putative mycorrhizal fungi (CMF) dataset
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Assign ecological guilds to fungal amplicon sequence variants
# (ASVs) using FUNGuildR and construct a subset containing ASVs
# classified as putative mycorrhizal fungi.
#
# INPUT FILES:
# 1. taxa_names_funguild.csv
#    Taxonomic classifications of ASVs, formatted for FUNGuildR.
#
# 2. seq_list.csv
#    Mapping between ASV identifiers (VCF1, VCF2, etc.) and
#    their corresponding nucleotide sequences.
#
# 3. seqtab_nochim.csv
#    ASV abundance table generated following DADA2 processing
#    and chimera removal.
#
# 4. rep-seqs_nochim.fna.csv
#    Representative ASV sequences exported in CSV format.
#
# OUTPUT FILES:
# 1. funguilds_all.csv
#    FUNGuild assignments for all fungal ASVs.
#
# 2. myco_funguilds.csv
#    FUNGuild assignments for ASVs with guild descriptions
#    containing the term "myco".
#
# 3. seq.list_myco.csv
#    ASV identifiers and nucleotide sequences retained in
#    the putative mycorrhizal dataset.
#
# 4. seqtab_myco.csv
#    Subset of the original ASV abundance table.
#
# 5. rep.seqs_myco.csv
#    Representative sequences of retained ASVs.
#
# NOTE:
# The filtering procedure retains ASVs whose FUNGuild guild
# annotation contains "myco", irrespective of confidence
# ranking or additional ecological guild assignments.
#
# ============================================================

# ------------------------------------------------------------
# 1. Set working directory and load required package
# ------------------------------------------------------------
setwd("../myco_dataset2/")

# Install FUNGuildR from GitHub.
devtools::install_github("brendanf/FUNGuildR") 

library(FUNGuildR)

# ------------------------------------------------------------
# 2. Import ASV taxonomic classifications
# ------------------------------------------------------------
# The input taxonomy file was prepared from assigned_taxa.csv,
# generated during DADA2 processing (Script 1).
#
# The original preparation was performed manually in Excel.
# The required format is:
#
# First column "ASV":
#   Unique ASV labels (VCF1, VCF2, ..., VCFn).
#
# Middle columns: one per taxonomic rank. Last column "Taxonomy":
#   Semicolon-separated taxonomic classifications with
#   prefixes (k_, p_, c_, o_, f_, g_, s_) removed. Between genus and species there is a space instead of semicolon.
#
# Example:
#   Fungi;Basidiomycota;Agaricomycetes;Cantharellales;Tulasnellaceae;Tulasnella sp.
#
# These VCF identifiers must correspond to the identifiers
# in the ASV-to-nucleotide sequence mapping file.DADA2 generates assigned_taxa.csv in that format by default.

tax_df <- 
  read.csv(
    "../FUNguild/taxa_names_funguild.csv") 

# ------------------------------------------------------------
# 3. Assign ecological guilds using FUNGuildR
# ------------------------------------------------------------

# Assign ecological guild classifications to all fungal ASVs
# in the General Fungi (GF) dataset.

myco.tax <- funguild_assign(tax_df)

# Retain the complete FUNGuild assignment table.

funguilds_all <- myco.tax #save data frame

# ------------------------------------------------------------
# 4. Identify putative mycorrhizal ASVs
# ------------------------------------------------------------

# Select ASVs whose assigned ecological guild contains
# the substring "myco" (case-insensitive).
#
# This includes guild descriptions containing terms such
# as "Mycorrhizal" and combinations of ecological guilds
# that include a mycorrhizal designation.
#

myco_funguilds <- myco.tax[
  grepl("myco", myco.tax$guild, ignore.case = TRUE),
]

# ------------------------------------------------------------
# 5. Import ASV identifiers and nucleotide sequences
# ------------------------------------------------------------

# Import the mapping between ASV identifiers (VCF1 to VCFn)
# and their unique nucleotide sequences.
#
# Required columns:
#   VCF_name : ASV identifier VCF1 to VCFn
#   seq      : Corresponding nucleotide sequence
#
# The mapping must remain consistent with the DADA2
# sequence table and representative sequence file.

seq.list <- read.csv("../FUNguild/seq_list.csv")

# ------------------------------------------------------------
# 6. Retain ASVs classified as putative mycorrhizal fungi
# ------------------------------------------------------------

# Match ASV identifiers against the FUNGuild classifications.
# Remove ASVs without a guild annotation containing "myco".
#
# The resulting table contains only ASVs retained for the
# putative mycorrhizal dataset.

seq.list_myco <- seq.list[
  seq.list[,1] %in% myco_funguilds$ASV,
]

# ------------------------------------------------------------
# 7. Import ASV abundance table and representative sequences
# ------------------------------------------------------------

# Import the DADA2-derived ASV abundance table and
# representative sequences in CSV format.

seqtab <- read.csv("../FUNguild/seqtab_nochim.csv")
rep.seqs <- read.csv("../FUNguild/rep-seqs_nochim.fna.csv")

# ------------------------------------------------------------
# 8. Construct the putative mycorrhizal fungi dataset
# ------------------------------------------------------------

# Retain ASVs matching the identifiers in seq.list_myco.
#
# Both the ASV abundance table and representative sequence
# table are filtered using the same ASV identifiers.

seqtab_myco <- subset(
  seqtab,
  VCF_name %in% seq.list_myco$VCF_name
)

rep.seqs_myco <- subset(
  rep.seqs,
  VCF_name %in% seq.list_myco$VCF_name
)

# ------------------------------------------------------------
# 9. Export the putative mycorrhizal dataset
# ------------------------------------------------------------

# Export the filtered ASV abundance table.

write.csv(
  seqtab_myco,
  "./seqtab_myco.csv"
)
# Export representative sequences of retained ASVs
write.csv(
  rep.seqs_myco,
  "./rep-seqs_myco.csv"
)
# Export the ASV identifiers and corresponding sequences.
write.csv(
  seq.list_myco,
  "./seq.list_myco.csv"
)
# ------------------------------------------------------------
# 10. Export FUNGuild classification results
# ------------------------------------------------------------

# Export ecological guild assignments for all fungal ASVs.
write.csv(
  funguilds_all,
  "./funguilds_all.csv"
)

# Export ecological guild assignments for the ASVs
# retained in the putative mycorrhizal dataset.

write.csv(
  myco_funguilds,
  "./myco_funguilds.csv"
)

# ============================================================
# End of Script 2
# ============================================================