# ============================================================
# Script 12. Fungal taxon prevalence and ASV richness
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Calculate the prevalence and ASV richness of fungal
# taxonomic groups in the General Fungi (GF) dataset.
#
# ANALYSES:
# 1. Calculate the percentage of retained orchid individuals
#    in which a fungal taxon or group was detected.
# 2. Calculate the number of distinct ASVs assigned to each
#    fungal group.
# 3. Summarize the most frequently detected fungal orders
#    across all samples and separately by orchid species.
# 4. Examine family-level composition within Cantharellales.
#
# INPUT FILES:
# - ../allorch2otu.txt
# - ../taxonomy.tsv
# - ../qiime2metadata.tsv
# - ../outtreeorchidsonlyperl.nwk
# - ../funguilds_all_fungi_ASVnames.tsv
#
# OUTPUT FILES:
# - ../fungal_group_frequencies.csv
# - ../top10ordersperorchid.csv
# - ../top10orders.csv
#
# DEFINITIONS:
# Prevalence (%) = 100 * (number of retained orchid samples
# containing at least one ASV from a fungal group) /
# (number of retained orchid samples examined).
#
# ASV richness = number of distinct ASVs assigned to a
# particular fungal taxonomic group.
#
# NOTE:
# The calculations use presence/absence rather than
# read abundance.

# Load libraries after installation.
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library("devtools")
library("vegan")
library("tidyverse")

# ------------------------------------------------------------
# 1. Import the GF ASV abundance table
# ------------------------------------------------------------

#Read in ASV table
otu <- read.table(file = "../allorch2otu.txt", sep = "\t", header = T, row.names = 1, 
                  comment.char = "")

#substitute ".." with "--" in barcode names
otu <- otu %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))

### Import taxonomy manually curated from DADA2 assignments ("assigned_taxa.csv")
# to match the phyloseq format.
taxonomy <- read.table(file = "../taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# Clean taxonomy prefixes and fill unclassified ranks. Greengenes format
tax <- taxonomy 

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- "" #remove NAs from the table
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

#Import sample metadata.
metadata <- read.table(file = "../qiime2metadata.tsv", sep = "\t", header = T, row.names = 1, quote = "")


# ------------------------------------------------------------
# 2. Construct phyloseq object and rarefy samples
# ------------------------------------------------------------
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
# Import the fungal phylogeny exported from QIIME 2.
# ASV identifiers must match between the tree, abundance
# table, and taxonomy file.
TREE = read_tree("../outtreeorchidsonlyperl.nwk") #Note that ASV identifiers need to be descriptive and unique, e.g. "Fusarium_babinda_VCF146", not nucleotide sequences. ASV labels need to match ASV table and taxonomy file.
# merge the data
ps <- phyloseq(OTU, SAMPLE, TAX, TREE)

# Rarefy to 679 reads per sample without replacement. Number inferred from QIIME2 rarefaction curve based on alpha diversity iterative estimates.
# The random seed is fixed for reproducibility.
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=679, replace=F)
ps.rarefied_meta = data.frame(sample_data(ps.rarefied))

# Extract the rarefied ASV abundance table.
#
# ASVs are expected as rows and samples as columns.
otu_df <- as.data.frame(otu_table(ps.rarefied))

# ------------------------------------------------------------
# 3. Import FUNGuild annotations
# ------------------------------------------------------------
# ASV identifiers must match those in the phyloseq object.
# Taxonomic and guild annotations are associated with ASVs.
guild <- read.table(file = "../funguilds_all_fungi_ASVnames.tsv", sep = "\t", header = T, row.names = 1, quote = "")

# Retain guild annotations corresponding to ASVs
# in the rarefied abundance table.
#
# match() preserves the ordering of the ASV table.
guilds_orchid <- guild[match(rownames(otu_df), rownames(guild)), ]

# ------------------------------------------------------------
# 4. Convert abundance to presence/absence
# ------------------------------------------------------------
otu_pa <- otu_df
otu_pa[otu_pa > 0] <- 1 #make every ASV occurance count as 1 in each orchid instead of ASV read number

# ------------------------------------------------------------
# 5. Combine ASV occurrence data and annotations
# ------------------------------------------------------------
# Add explicit ASV identifier columns to both tables.
guilds_orchid$ASV <- rownames(guilds_orchid)
otu_pa$ASV <- rownames(otu_pa)
# Join the presence/absence matrix with fungal annotations.
#
# Each row corresponds to an ASV.
# Sample columns contain 0/1 occurrence values.
# Additional columns contain taxonomy and guild information.
otu2 <-
  left_join(
    otu_pa,
    guilds_orchid,
    by="ASV"
  )

# ============================================================
# FUNGAL GROUP PREVALENCE
# ============================================================

# ------------------------------------------------------------
# 6. Define a reusable prevalence calculation function
# ------------------------------------------------------------ 


# This function calculates:
#
# ASVs:
#   Total number of distinct ASVs assigned to a fungal group
#   across the complete retained dataset.
#
# Percent:
#   Percentage of orchid individuals in which at least
#   one ASV assigned to the group was detected.
#
# N:
#   Number of orchid individuals containing the group.
#
# Results are reported:
# - Overall, across all retained orchid individuals.
# - Separately for each orchid species.

calc_frequency <- function(otu2,
                           ps.rarefied_meta,
                           filter_expr,
                           group_name = NULL){
  
  filter_expr <- enquo(filter_expr)
  
  # Sample names
  sample_cols <- rownames(ps.rarefied_meta)
  
  # Filter ASVs
  x <- otu2 %>%
    filter(!!filter_expr)
  
  # Stop if nothing matched
  if(nrow(x) == 0){
    stop("No ASVs matched the filtering expression.")
  }
  
  # Number of ASVs
  n_asvs <- nrow(x)
  
  # Collapse ASVs to presence/absence per sample
  sample_pa <- x %>%
    select(all_of(sample_cols)) %>%
    summarise(across(everything(), max)) %>%
    pivot_longer(
      everything(),
      names_to = "Sample",
      values_to = "Presence"
    )
  
  # Add metadata
  meta2 <- ps.rarefied_meta %>%
    rownames_to_column("Sample")
  
  sample_pa <- sample_pa %>%
    left_join(meta2, by = "Sample")
  
  # Assign the fungal group name.
  if(is.null(group_name)){
    group_name <- as_label(filter_expr)
  }
  
  # Calculate prevalence across all orchid individuals.
  overall <- sample_pa %>%
    summarise(
      Percent = mean(Presence) * 100,
      N = sum(Presence)
    ) %>%
    mutate(
      Group = group_name,
      orchid.sp. = "Overall",
      ASVs = n_asvs
    ) %>%
    select(Group, orchid.sp., ASVs, Percent, N)
  
  # Calculate prevalence separately for each orchid species.
  by_species <- sample_pa %>%
    group_by(orchid.sp.) %>%
    summarise(
      Percent = mean(Presence) * 100,
      N = sum(Presence),
      .groups = "drop"
    ) %>%
    mutate(
      Group = group_name,
      ASVs = n_asvs
    ) %>%
    select(Group, orchid.sp., ASVs, Percent, N)
  # Combine overall and species-specific summaries.
  bind_rows(overall, by_species)
}

# ------------------------------------------------------------
# 7. Calculate prevalence of selected fungal groups
# ------------------------------------------------------------

#Traditional orchid mycorrhizal fungal groups
tulas <- calc_frequency(otu2, ps.rarefied_meta, Family == "Tulasnellaceae", "Tulasnellaceae")
ser <- calc_frequency(otu2, ps.rarefied_meta,Family == "Serendipitaceae", "Serendipitaceae")
cera <- calc_frequency(otu2, ps.rarefied_meta, Family == "Ceratobasidiaceae", "Ceratobasidiaceae")

#Other putatively mycorrhizal fungal groups
atra <- calc_frequency(otu2, ps.rarefied_meta,Order == "Atractiellales", "Atractiellales")
myc <- calc_frequency(otu2, ps.rarefied_meta, Family == "Mycenaceae", "Mycenaceae")
armi <- calc_frequency(otu2, ps.rarefied_meta,Genus == "Armillaria", "Armillaria")
fus <- calc_frequency(otu2, ps.rarefied_meta, Genus == "Fusarium", "Fusarium")

#Major non-mycorrhizal fungal groups
asco <- calc_frequency(otu2, ps.rarefied_meta, Phylum == "Ascomycota", "Ascomycota")
xyl <- calc_frequency(otu2, ps.rarefied_meta,Order == "Xylariales", "Xylariales")
chae <- calc_frequency(otu2, ps.rarefied_meta,Order == "Chaetothyriales", "Chaetothyriales")
hypo <- calc_frequency(otu2, ps.rarefied_meta,Order == "Hypocreales", "Hypocreales")
capno <- calc_frequency(otu2, ps.rarefied_meta,Order == "Capnodiales", "Capnodiales")
coll <- calc_frequency(otu2, ps.rarefied_meta,Genus == "Colletotrichum", "Colletotrichum")
trich <- calc_frequency(otu2, ps.rarefied_meta,Genus == "Trichoderma", "Trichoderma")
helo <- calc_frequency(otu2, ps.rarefied_meta,Order == "Helotiales", "Helotiales")
trechi <- calc_frequency(otu2, ps.rarefied_meta,Order == "Trechisporales", "Trechisporales")
# Combine all selected fungal group prevalence results.
results <- bind_rows(tulas, ser, cera, atra, myc, armi, fus, asco, xyl, chae, hypo, capno, coll, trich, helo, trechi)
write.csv(results, file = "../fungal_group_frequencies.csv",row.names = FALSE)


# ============================================================
# MOST FREQUENT FUNGAL ORDERS BY ORCHID SPECIES
# ============================================================

# ------------------------------------------------------------
# 8. Calculate order prevalence within orchid species
# ------------------------------------------------------------

# For each fungal order, calculate the percentage of orchid
# individuals containing at least one ASV from that order.
#
# The results are calculated separately for each of the
# four focal Lepanthes species.
#
# N represents the number of orchid individuals in which
# the order was detected.

order_frequency_species <- otu2 %>%
  filter(!is.na(Order), Order != "") %>%
  group_by(Order) %>%
  summarise(
    across(
      all_of(rownames(ps.rarefied_meta)),
      ~ max(.x)
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Order,
    names_to = "Sample",
    values_to = "Presence"
  ) %>%
  left_join(
    ps.rarefied_meta %>%
      rownames_to_column("Sample"),
    by = "Sample"
  ) %>%
  group_by(orchid.sp., Order) %>%
  summarise(
    Percent = mean(Presence) * 100,
    N = sum(Presence),
    .groups = "drop"
  ) %>%
  arrange(orchid.sp., desc(Percent))

# ------------------------------------------------------------
# 9. Select and export the most frequent orders by orchid species
# ------------------------------------------------------------
#Top 10 most frequent orders across orchid species
top10_orders_species <- order_frequency_species %>%
  group_by(orchid.sp.) %>%
  slice_head(n = 11) #10 identified orders plus incertae_sedis

write.csv(top10_orders_species, file = "../top10ordersperorchid.csv",row.names = FALSE)
#Top 10 most frequent orders across orchid species

# ============================================================
# MOST FREQUENT FUNGAL ORDERS ACROSS ALL ORCHID SAMPLES
# ============================================================

# ------------------------------------------------------------
# 10. Calculate ASV richness for each fungal order
# ------------------------------------------------------------

sample_cols <- rownames(ps.rarefied_meta)

# Count distinct ASVs assigned to each fungal order.
order_richness <- otu2 %>%
  filter(!is.na(Order), Order != "") %>%
  group_by(Order) %>%
  summarise(
    ASV_richness = n_distinct(ASV),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 11. Calculate overall fungal order prevalence
# ------------------------------------------------------------

order_frequency <- otu2 %>%
  filter(!is.na(Order), Order != "") %>%
  group_by(Order) %>%
  summarise(
    across(
      all_of(sample_cols),
      ~ as.integer(any(.x > 0))
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Order,
    names_to = "Sample",
    values_to = "Presence"
  ) %>%
  group_by(Order) %>%
  summarise(
    Percent = mean(Presence) * 100,
    N = sum(Presence),
    .groups = "drop"
  ) %>%
  left_join(order_richness, by = "Order") %>%
  arrange(desc(Percent))

# ------------------------------------------------------------
# 12. Select and export the most frequent orders
# ------------------------------------------------------------
top10_orders <- order_frequency %>%
  slice_head(n = 11) #10 identified orders plus incertae_sedis

write.csv(top10_orders, file = "../top10orders.csv",row.names = FALSE)

# ============================================================
# FAMILY COMPOSITION WITHIN CANTHARELLALES
# ============================================================
# Identify fungal families assigned to Cantharellales.
#
# For each family:
# - Calculate prevalence across retained orchid samples.
# - Count the number of samples containing the family.
# - Calculate the number of distinct ASVs assigned to it.
#
# This analysis helps distinguish the contribution of
# Tulasnellaceae, Ceratobasidiaceae, and other families
# to the overall prevalence of Cantharellales.

sample_cols <- rownames(ps.rarefied_meta)

canth_families <- otu2 %>%
  filter(Order == "Cantharellales",
         !is.na(Family),
         Family != "") %>%
  group_by(Family) %>%
  summarise(
    across(
      all_of(sample_cols),
      ~ as.integer(any(.x > 0))
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Family,
    names_to = "Sample",
    values_to = "Presence"
  ) %>%
  group_by(Family) %>%
  summarise(
    Percent = mean(Presence) * 100,
    N = sum(Presence),
    .groups = "drop"
  ) %>%
  left_join(
    otu2 %>%
      filter(Order == "Cantharellales",
             !is.na(Family),
             Family != "") %>%
      group_by(Family) %>%
      summarise(
        ASV_richness = n_distinct(ASV),
        .groups = "drop"
      ),
    by = "Family"
  ) %>%
  arrange(desc(Percent))
# Display family-level results.
canth_families


# ============================================================
# End of Script 12
# ============================================================

