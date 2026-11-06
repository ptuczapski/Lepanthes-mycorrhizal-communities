# ============================================================
# Script 6. Construction of the Confirmed Mycorrhizal Fungi
#           (CMF) ASV abundance table after singleton removal
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Combine species-specific CMF abundance tables generated
# following QIIME 2 filtering into a single ASV abundance
# matrix for downstream ecological analyses.
#
# Species-specific filtering retained ASVs occurring in
# at least two orchid individuals within each species.
#
# The original workflow includes manual processing of
# intermediate tables in Microsoft Excel.
#
# INPUT FILES:
# - LMe2-filtered_table_myco-nochim.txt
# - LFB2-filtered_table_myco-nochim.txt
# - LMo2-filtered_table_myco-nochim.txt
# - LC2-filtered_table_myco-nochim.txt
# - orchidsp-filtered_table_myco-nochim.txt
#
# OUTPUT FILE:
# - allorch2otu_myco.txt
#
# ABBREVIATIONS:
# CMF = Confirmed Mycorrhizal Fungi
# LMe = Lepanthes mentosa
# LFB = Lepanthes falx-bellica
# LMo = Lepanthes monteverdensis
# LC  = Lepanthes cribbii
#
# NOTE:
# ASV identifiers must be unique and consistent across
# all input tables. Taxonomic assignments should be
# maintained separately from ASV identifiers.
# ============================================================
# ------------------------------------------------------------
# 1. Load required R packages
# ------------------------------------------------------------
library("data.table")
library("tidyverse")
#####################################################################################################################
# ------------------------------------------------------------
# 2. Import species-specific CMF abundance tables
# ------------------------------------------------------------

# Import the QIIME 2-filtered feature tables.
#
# Each table contains ASVs as rows and orchid individuals
# as columns.
#
# Replace ".." with "--" in column names to restore
# sample identifiers modified during import.


# Lepanthes mentosa
LMe2otu_myco<-read.table(file = "../LMe2-filtered_table_myco-nochim.txt", sep = "\t", header = T, row.names = 1, 
                    comment.char = "")
LMe2otu_myco<- LMe2otu_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))
# Lepanthes falx-bellica
LFB2otu_myco<-read.table(file = "../LFB2-filtered_table_myco-nochim.txt", sep = "\t", header = T, row.names = 1, 
                    comment.char = "")
LFB2otu_myco<- LFB2otu_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))
# Lepanthes monteverdensis
LMo2otu_myco<-read.table(file = "../LMo2-filtered_table_myco-nochim.txt", sep = "\t", header = T, row.names = 1, 
                    comment.char = "")
LMo2otu_myco<- LMo2otu_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))
# Lepanthes cribbii
LC2otu_myco<-read.table(file = "../LC2-filtered_table_myco-nochim.txt", sep = "\t", header = T, row.names = 1, 
                   comment.char = "")
LC2otu_myco<- LC2otu_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))

# ------------------------------------------------------------
# 3. Import the complete orchid CMF abundance table
# ------------------------------------------------------------

# Import the CMF abundance table containing all four
# focal orchid species before species-specific singleton
# filtering.
#
# This table provides the reference set of ASVs used
# to align the species-specific tables.

allorchidsOTUdf_myco <- read.table(file = "../orchidsp-filtered_table_myco-nochim.txt", sep = "\t", header = T, row.names = 1, 
                              comment.char = "")
allorchidsOTUdf_myco<- allorchidsOTUdf_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 

# ------------------------------------------------------------
# 4. Convert ASV row names into explicit columns
# ------------------------------------------------------------
# Convert ASV row names into the column "rowname".
# This column serves as the matching key during full joins.
allorchidsOTU_myco_tibble <- rownames_to_column(allorchidsOTUdf_myco, var = "rowname")
allorchidsOTU_myco_tibble <- as.tibble(allorchidsOTU_myco_tibble)
LMe2otu_myco_tibble <- rownames_to_column(LMe2otu_myco, var = "rowname")
LMe2otu_myco_tibble <- as.tibble(LMe2otu_myco_tibble) 
LFB2otu_myco_tibble <- rownames_to_column(LFB2otu_myco, var = "rowname")
LFB2otu_myco_tibble <- as.tibble(LFB2otu_myco_tibble) 
LMo2otu_myco_tibble <- rownames_to_column(LMo2otu_myco, var = "rowname")
LMo2otu_myco_tibble <- as.tibble(LMo2otu_myco_tibble) 
LC2otu_myco_tibble <- rownames_to_column(LC2otu_myco, var = "rowname")
LC2otu_myco_tibble <- as.tibble(LC2otu_myco_tibble) 
# ------------------------------------------------------------
# 5. Align Lepanthes mentosa CMF ASVs
# ------------------------------------------------------------

# Join the L. mentosa CMF table with the complete orchid
# CMF table using ASV identifiers.
#
# The full join retains ASVs present in either table.
#
# MANUAL PROCESSING:
# Following the full_join() operations, each species-specific
# table was manually processed in Microsoft Excel.
#
# The procedure was performed independently for:
# 1. Lepanthes mentosa
# 2. Lepanthes falx-bellica
# 3. Lepanthes monteverdensis
# 4. Lepanthes cribbii
#
# The manual processing involved:
# - Moving ASV identifiers into the first column.
# - Removing redundant columns introduced by the joins.
# - Retaining the species-specific abundance columns.
# - Saving the modified tables as tab-delimited text files.
#
# The corrected tables were subsequently imported into R.
#
# IMPORTANT:
# The following code does not automatically perform
# the manual Excel modifications.

joinLMe_myco<- full_join(LMe2otu_myco_tibble, allorchidsOTU_myco_tibble, keep=TRUE, by = "rowname") 
joinLMedf_myco = as.data.frame(joinLMe_myco, header = F)
write.table(joinLMedf_myco,"../joinLMedf_myco.txt",sep="\t",row.names=F)
# Read the manually corrected intermediate table.
joinLMedf_myco <- read.table(file = "../joinLMedf_myco.txt", sep = "\t", header = T, row.names = 1, 
                        comment.char = "")
joinLMedf_myco<- joinLMedf_myco %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLMedf_myco<- joinLMedf_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
# Replace missing values with zero.
# Missing abundance values represent
# ASVs absent after species-specific filtering.
joinLMedf_myco <- replace(joinLMedf_myco, is.na(joinLMedf_myco), 0)

# ------------------------------------------------------------
# 6. Align Lepanthes falx-bellica ASVs
# ------------------------------------------------------------
joinLFB_myco<- full_join(LFB2otu_myco_tibble, allorchidsOTU_myco_tibble, keep=TRUE, by = "rowname") 
joinLFBdf_myco = as.data.frame(joinLFB_myco, header = F)
write.table(joinLFBdf_myco,"../joinLFBdf_myco.txt",sep="\t",row.names=F)
# REVIEW REQUIRED:
# Apply the same manual correction to this intermediate
# table before importing it.
joinLFBdf_myco <- read.table(file = "../joinLFBdf_myco.txt", sep = "\t", header = T, row.names = 1, 
                        comment.char = "")
joinLFBdf_myco<- joinLFBdf_myco %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLFBdf_myco<- joinLFBdf_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
joinLFBdf_myco <- replace(joinLFBdf_myco, is.na(joinLFBdf_myco), 0)

# ------------------------------------------------------------
# 7. Align Lepanthes monteverdensis ASVs
# ------------------------------------------------------------
joinLMo_myco<- full_join(LMo2otu_myco_tibble, allorchidsOTU_myco_tibble, keep=TRUE, by = "rowname") 
joinLModf_myco = as.data.frame(joinLMo_myco, header = F)
write.table(joinLModf_myco,"../joinLModf_myco.txt",sep="\t",row.names=F)
# REVIEW REQUIRED:
# Apply the same manual correction to this intermediate
# table before importing it.
joinLModf_myco <- read.table(file = "../joinLModf_myco.txt", sep = "\t", header = T, row.names = 1, 
                        comment.char = "")
joinLModf_myco<- joinLModf_myco %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLModf_myco<- joinLModf_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
joinLModf_myco <- replace(joinLModf_myco, is.na(joinLModf_myco), 0)

# ------------------------------------------------------------
# 8. Align Lepanthes cribbii ASVs
# ------------------------------------------------------------
joinLC_myco<- full_join(LC2otu_myco_tibble, allorchidsOTU_myco_tibble, keep=TRUE, by = "rowname") 
joinLCdf_myco = as.data.frame(joinLC_myco, header = F)
write.table(joinLCdf_myco,"../joinLCdf_myco.txt",sep="\t",row.names=F)
# REVIEW REQUIRED:
# Apply the same manual correction to this intermediate
# table before importing it.
joinLCdf_myco <- read.table(file = "../joinLCdf_myco.txt", sep = "\t", header = T, row.names = 1, 
                       comment.char = "")
joinLCdf_myco<- joinLCdf_myco %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLCdf_myco<- joinLCdf_myco %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
joinLCdf_myco <- replace(joinLCdf_myco, is.na(joinLCdf_myco), 0)

# ------------------------------------------------------------
# 9. Match ASV order across species-specific tables
# ------------------------------------------------------------

# Following manual processing, the four species-specific
# CMF tables contain ASVs in different row orders.
#
# Align ASVs by their unique identifiers before combining
# the abundance columns.
#
# The L. falx-bellica table provides the reference order.
# Align Lepanthes monteverdensis.
sortLMo_myco <- joinLModf_myco[match(rownames(joinLFBdf_myco), rownames(joinLModf_myco)), ]
sortLMo_myco<-as.data.frame(sortLMo_myco)
bind_myco <- cbind(joinLFBdf_myco,sortLMo_myco ) 
# Align Lepanthes mentosa.
sortLMe_myco <- joinLMedf_myco[match(rownames(bind_myco), rownames(joinLMedf_myco)), ]
sortLMe_myco<-as.data.frame(sortLMe_myco)
bind_myco <- cbind(bind_myco,sortLMe_myco ) 
# Align Lepanthes cribbii.
sortLC_myco <- joinLCdf_myco[match(rownames(bind_myco), rownames(joinLCdf_myco)), ]
sortLC_myco<-as.data.frame(sortLC_myco)
bind_myco <- cbind(bind_myco,sortLC_myco ) 

# ------------------------------------------------------------
# 10. Construct the combined CMF abundance table
# ------------------------------------------------------------

# Combine the aligned species-specific CMF abundance tables.
#
# Rows represent fungal ASVs.
# Columns represent individual orchid samples.
#
# Values represent ASV read abundances following
# species-specific singleton filtering.

allorch2otu_myco <- bind_myco

# ------------------------------------------------------------
# 11. Export the combined CMF abundance table
# ------------------------------------------------------------
write.table(allorch2otu_myco,"../allorch2otu_myco.txt",sep="\t",row.names=TRUE)

# ============================================================
# End of Script 6
# ============================================================
