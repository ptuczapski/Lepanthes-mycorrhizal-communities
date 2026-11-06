# ============================================================
# Script 5. Construction of the General Fungi (GF) ASV table
#           following species-specific singleton filtering
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Combine the species-specific ASV abundance tables generated
# in Script 4 into a single GF abundance table.
#
# Species-specific filtering was performed in QIIME 2 using
# a minimum occurrence threshold of two orchid individuals
# per species.
#
# The workflow retains a common set of ASV identifiers across
# the four orchid species, allowing their abundance tables
# to be combined for downstream community analyses.
# INPUT FILES:
# - outLMe2filteredtablenochim.txt
# - outLFB2filteredtablenochim.txt
# - outLMo2filteredtablenochim.txt
# - outLC2filteredtablenochim.txt
# - outorchidspfilteredtablenochim.txt
#
# OUTPUT FILE:
# - allorch2otu.txt
#
# ABBREVIATIONS:
# GF  = General Fungi
# LMe = Lepanthes mentosa
# LFB = Lepanthes falx-bellica
# LMo = Lepanthes monteverdensis
# LC  = Lepanthes cribbii
#
# IMPORTANT:
# The original workflow includes manual modifications of
# intermediate tables. 
#
# ============================================================
# ------------------------------------------------------------
# 1. Load required R packages
# ------------------------------------------------------------
library("data.table")
library("tidyverse")
# ------------------------------------------------------------
# 2. Import species-specific ASV abundance tables
# ------------------------------------------------------------

# Import QIIME 2 feature tables following species-specific
# singleton filtering.
#
# Each table contains ASVs as rows and orchid individuals
# as columns.
#
# ASVs occurring in fewer than two individuals within
# the corresponding orchid species were excluded.
#
# Column names are standardized by replacing ".." with "--".
# This restores characters altered during R file import.

# Lepanthes mentosa
LMe2otu<-read.table(file = "../outLMe2filteredtablenochim.txt", sep = "\t", header = T, row.names = 1, 
                    comment.char = "")
LMe2otu<- LMe2otu %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))

# Lepanthes falx-bellica
LFB2otu<-read.table(file = "../outLFB2filteredtablenochim.txt", sep = "\t", header = T, row.names = 1, 
                    comment.char = "")
LFB2otu<- LFB2otu %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))

# Lepanthes monteverdensis
LMo2otu<-read.table(file = "../outLMo2filteredtablenochim.txt", sep = "\t", header = T, row.names = 1, 
                    comment.char = "")
LMo2otu<- LMo2otu %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))

# Lepanthes cribbii
LC2otu<-read.table(file = "../outLC2filteredtablenochim.txt", sep = "\t", header = T, row.names = 1, 
                   comment.char = "")
LC2otu<- LC2otu %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))

# ------------------------------------------------------------
# 3. Import the complete orchid GF abundance table
# ------------------------------------------------------------

# Import the ASV abundance table containing all four focal
# orchid species before species-specific singleton removal.
#
# This table provides the reference set of ASV identifiers
# used to align the species-specific tables.

allorchidsOTUdf <- read.table(file = "../outorchidspfilteredtablenochim.txt", sep = "\t", header = T, row.names = 1, 
                              comment.char = "")
allorchidsOTUdf<- allorchidsOTUdf %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 

# ------------------------------------------------------------
# 4. Convert row names into explicit ASV identifier columns
# ------------------------------------------------------------

# Convert ASV row names to the column "rowname".
#
# This column is subsequently used to match ASVs during
# full joins between the species-specific and complete tables.

allorchidsOTU_tibble <- rownames_to_column(allorchidsOTUdf, var = "rowname")
allorchidsOTU_tibble <- as.tibble(allorchidsOTU_tibble)
LMe2otu_tibble <- rownames_to_column(LMe2otu, var = "rowname")
LMe2otu_tibble <- as.tibble(LMe2otu_tibble) 
LFB2otu_tibble <- rownames_to_column(LFB2otu, var = "rowname")
LFB2otu_tibble <- as.tibble(LFB2otu_tibble) 
LMo2otu_tibble <- rownames_to_column(LMo2otu, var = "rowname")
LMo2otu_tibble <- as.tibble(LMo2otu_tibble) 
LC2otu_tibble <- rownames_to_column(LC2otu, var = "rowname")
LC2otu_tibble <- as.tibble(LC2otu_tibble)

# ------------------------------------------------------------
# 5. Align Lepanthes mentosa ASVs
# ------------------------------------------------------------

# Join the species-specific L. mentosa table with the
# complete GF table using ASV identifiers.
#
# The full join retains ASVs present in either table.
#
# MANUAL PROCESSING:
#
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

joinLMe<- full_join(LMe2otu_tibble, allorchidsOTU_tibble, keep=TRUE, by = "rowname") 
joinLMedf = as.data.frame(joinLMe, header = F)
write.table(joinLMedf,"../joinLMedf.txt",sep="\t",row.names=F)

# Read the manually corrected intermediate table.
joinLMedf <- read.table(file = "../joinLMedf.txt", sep = "\t", header = T, row.names = 1, 
                        comment.char = "")
joinLMedf<- joinLMedf %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLMedf<- joinLMedf %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
# Replace missing values with zero.
# Missing abundance values represent
# ASVs absent after species-specific filtering.
joinLMedf <- replace(joinLMedf, is.na(joinLMedf), 0)

# ------------------------------------------------------------
# 6. Align Lepanthes falx-bellica ASVs
# ------------------------------------------------------------
joinLFB<- full_join(LFB2otu_tibble, allorchidsOTU_tibble, keep=TRUE, by = "rowname") 
joinLFBdf = as.data.frame(joinLFB, header = F)
write.table(joinLFBdf,"../joinLFBdf.txt",sep="\t",row.names=F)
# REVIEW REQUIRED:
# Apply the same manual correction to this intermediate
# table before importing it.
joinLFBdf <- read.table(file = "../joinLFBdf.txt", sep = "\t", header = T, row.names = 1, 
                        comment.char = "")
joinLFBdf<- joinLFBdf %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLFBdf<- joinLFBdf %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
joinLFBdf <- replace(joinLFBdf, is.na(joinLFBdf), 0)

# ------------------------------------------------------------
# 7. Align Lepanthes monteverdensis ASVs
# ------------------------------------------------------------
joinLMo<- full_join(LMo2otu_tibble, allorchidsOTU_tibble, keep=TRUE, by = "rowname") 
joinLModf = as.data.frame(joinLMo, header = F)
write.table(joinLModf,"../joinLModf.txt",sep="\t",row.names=F)
# REVIEW REQUIRED:
# The original output and input paths differ.
# Confirm that both refer to the intended intermediate file.
joinLModf <- read.table(file = "../OTUtablephyloseq/joinLModf.txt", sep = "\t", header = T, row.names = 1, 
                        comment.char = "")
joinLModf<- joinLModf %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLModf<- joinLModf %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
joinLModf <- replace(joinLModf, is.na(joinLModf), 0)

# ------------------------------------------------------------
# 8. Align Lepanthes cribbii ASVs
# ------------------------------------------------------------
joinLC<- full_join(LC2otu_tibble, allorchidsOTU_tibble, keep=TRUE, by = "rowname") 
joinLCdf = as.data.frame(joinLC, header = F)
write.table(joinLCdf,"../joinLCdf.txt",sep="\t",row.names=F)
# REVIEW REQUIRED:
# Apply the same manual correction to this intermediate
# table before importing it.
joinLCdf <- read.table(file = "../joinLCdf.txt", sep = "\t", header = T, row.names = 1, 
                       comment.char = "")
joinLCdf<- joinLCdf %>%
  rename_with(~ gsub(".x", "", .x, fixed=TRUE)) 
joinLCdf<- joinLCdf %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE)) 
joinLCdf <- replace(joinLCdf, is.na(joinLCdf), 0)

# ------------------------------------------------------------
# 9. Match ASV order across species-specific tables
# ------------------------------------------------------------

# The four species-specific tables may contain ASVs in
# different orders.
#
# Match row names to ensure that corresponding ASVs occupy
# the same row before combining abundance information.
#
# The L. falx-bellica table provides the reference order.

sortLMo <- joinLModf[match(rownames(joinLFBdf), rownames(joinLModf)), ]
sortLMo<-as.data.frame(sortLMo)
bind <- cbind(joinLFBdf,sortLMo ) 
# Align Lepanthes mentosa.
sortLMe <- joinLMedf[match(rownames(bind), rownames(joinLMedf)), ]
sortLMe<-as.data.frame(sortLMe)
bind <- cbind(bind,sortLMe ) 
# Align Lepanthes cribbii.
sortLC <- joinLCdf[match(rownames(bind), rownames(joinLCdf)), ]
sortLC<-as.data.frame(sortLC)
bind <- cbind(bind,sortLC ) 

# ------------------------------------------------------------
# 10. Construct the combined CMF abundance table
# ------------------------------------------------------------

# Combine the aligned species-specific GF abundance tables.
#
# Rows represent fungal ASVs.
# Columns represent individual orchid samples.
#
# Values represent ASV read abundances following
# species-specific singleton filtering.

allorch2otu <- bind

#Save singleton-filtered table
write.table(allorch2otu,"../allorch2otu.txt",sep="\t",row.names=TRUE)

# ============================================================
# End of Script 5
# ============================================================
