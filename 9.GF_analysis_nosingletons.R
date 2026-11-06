# ============================================================
# Script 9. General Fungi (GF): community composition analyses
# ============================================================
# Created by Piotr Tymon Tuczapski
#
# PURPOSE:
# Analyze fungal community composition in the GF dataset.
#
# ANALYSES:
# 1. Construct a phyloseq object.
# 2. Rarefy samples to 679 sequencing reads.
# 3. Generate Jaccard PCoA and NMDS ordinations.
# 4. Test Jaccard multivariate dispersion (PERMDISP).
# 5. Generate unweighted UniFrac PCoA and NMDS ordinations.
# 6. Test unweighted UniFrac multivariate dispersion.
# 7. Perform PERMANOVA to evaluate host-tree and orchid
#    species effects on fungal community composition.
#
# INPUTS:
# ../allorch2otu.txt
# ../taxonomy.tsv
# ../qiime2metadata.tsv
# ../outtreeorchidsonlyperl.nwk
# ============================================================
#
# Below is from "Handoff to phyloseq" by
# https://benjjneb.github.io/dada2/tutorial_1_8.html
# Load libraries after installation.

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library("devtools")
library("vegan")
library("tidyverse")

# ------------------------------------------------------------
# 1. Import ASV abundance table, taxonomy, and metadata
# ------------------------------------------------------------

# Import the ASV abundance table.

otu <- read.table(file = "../allorch2otu.txt", sep = "\t", header = T, row.names = 1, 
                  comment.char = "")

#substitute ".." with "--" in barcode names
otu <- otu %>%
  rename_with(~ gsub("..", "--", .x, fixed=TRUE))

## Import taxonomy manually curated from DADA2 assignments ("assigned_taxa.csv")
# to match the phyloseq format
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
# Inspect the number of samples per orchid species.
ps.rarefied_meta %>% count(orchid.sp.) 

# ============================================================
# JACCARD DISSIMILARITY ANALYSES
# ============================================================
# Calculate presence/absence Jaccard dissimilarity.
jacc <- phyloseq::distance(ps.rarefied, "jaccard", binary = T)

# ------------------------------------------------------------
# 3. Jaccard PCoA
# ------------------------------------------------------------

jacc.ordination = ordinate(ps.rarefied, method="PCoA", distance=jacc) 

p_jacc <- plot_ordination(ps.rarefied, jacc.ordination, color = "orchid.sp.", shape = "orchid.sp.")

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
p_jacc + geom_point(size = 4)  +theme_bw() +theme_classic()+ 
  ggtitle("Jaccard Similarity Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust=2.5)) +
  geom_point(size=4) +
  theme(text = element_text(size = 15)) +
  theme(plot.margin=unit(c(1,1,1,1), 'cm')) +
  labs("orchid.sp.") + scale_color_manual(values = c("salmon", "darkgreen", "mediumorchid1", "dodgerblue2"), name = "Orchid species:") +
  scale_shape_manual(values = c(15, 17, 8, 21), name = "Orchid species:")
dev.off()
shell.exec(tf)

# ------------------------------------------------------------
# 4. Jaccard NMDS
# ------------------------------------------------------------
# Evaluate NMDS stress across 1–10 dimensions.

stress_vec <- numeric(10)
for(i in seq(10)){
  stress_vec[i] <- vegan::metaMDS(jacc, weakties=FALSE, k=i)$stress
}

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
plot(seq(10),stress_vec, type = 'o', ylab = "Stress", xlab = "Number of dimensions",
     col="tomato", pch=19)
abline(h=0.2, lty=2) #adds dotted line
abline(h=0.1, lty=3)
dev.off()
shell.exec(tf)

# Fit an eight-dimensional NMDS solution.
# The resulting figure displays only NMDS1 and NMDS2.
set.seed(1)
orchid.sp._NMDS <- vegan::metaMDS(jacc, weakties=FALSE, trymax=200, k =8)  

# Generate Shepard stress plot.

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
vegan::stressplot(orchid.sp._NMDS)
dev.off()
shell.exec(tf)

# Extract sample ordination scores.
orchid.sp._NMDS2 <- orchid.sp._NMDS %>%  
  vegan::scores() %>%                                         
  as_tibble(rownames = "Sample")

# Convert phyloseq object to a data frame.
ps.rarefied_df = psmelt(ps.rarefied) 
ps.rarefied_tibble <- as.tibble(ps.rarefied_df) 

nmds <- inner_join(ps.rarefied_tibble, orchid.sp._NMDS2, by="Sample") #join tibbles by OTU column so that NMDS scores are added

# Calculate group centroids from the joined coordinates.
centroid <- nmds %>%
  group_by(orchid.sp.) %>%
  summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=orchid.sp.)) +
  geom_point(aes(shape = orchid.sp., color = orchid.sp.),size = 4)  +theme_bw() +theme_classic()+ 
  scale_shape_manual(values = c(15, 17, 8, 21), name = "Orchid species:") +
  ggtitle("Jaccard Similarity Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust=2.5)) +
  theme(text = element_text(size = 15)) + 
  theme(plot.margin=unit(c(1,1,1,1), 'cm')) +
  labs("orchid.sp.") + scale_color_manual(values = c("salmon", "darkgreen", "mediumorchid1", "dodgerblue2"), name = "Orchid species:") +
  geom_point(data=centroid, size=7, shape=21, color="black",
             aes(fill=orchid.sp.), show.legend=FALSE) + scale_fill_manual(values = c("salmon", "darkgreen", "mediumorchid1", "dodgerblue2"))


dev.off()
shell.exec(tf)


# ------------------------------------------------------------
# 5. Jaccard PERMDISP by orchid species
# ------------------------------------------------------------

ps.rarefied_metadata <- as(sample_data(ps.rarefied), "data.frame") #export metadata from rarefied phyloseq object 
groups_orchids <- ps.rarefied_metadata[["orchid.sp."]]
jacc_permdisp <- vegan::betadisper(jacc, groups_orchids, type=c("median"))

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
plot(jacc_permdisp, main = "Jaccard dist. to centroid", pch=c(15, 17, 8, 21), cex = 0.7, label.cex=0.8)+ points(jacc_permdisp$centroids, cex=1.2) 
dev.off()
shell.exec(tf)

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
boxplot(jacc_permdisp, main = "", xlab = "") #shows how much variation there is in distances around centroids
dev.off()
shell.exec(tf)

# Test differences in multivariate dispersion.
anova(jacc_permdisp)
Tukey <- TukeyHSD(jacc_permdisp) #which groups have dif. dispersion?
plot(Tukey)
Tukey <- (as.data.frame(Tukey$group))
write.table(Tukey,"../Tukeyjaccnosingletons.txt",sep="\t",row.names=TRUE)

# ============================================================
# UNWEIGHTED UNIFRAC ANALYSES
# ============================================================
#Calculate unweighted UniFrac distances
unifrac = phyloseq::distance(ps.rarefied, method="unifrac", weighted=F)

# ------------------------------------------------------------
# 6. Unweighted UniFrac PCoA
# ------------------------------------------------------------

unif.ordination = ordinate(ps.rarefied, method="PCoA", distance=unifrac) 

p_unif <- plot_ordination(ps.rarefied, unif.ordination, color = "orchid.sp.", shape = "orchid.sp.")

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
p_unif + geom_point(size = 4)  +theme_bw() +theme_classic()+ 
  ggtitle("Unweighted Unifrac") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust=2.5)) +
  geom_point(size=4) +
  theme(text = element_text(size = 15)) +
  theme(plot.margin=unit(c(1,1,1,1), 'cm')) +
  labs("orchid.sp.") + scale_color_manual(values = c("salmon", "darkgreen", "mediumorchid1", "dodgerblue2"), name = "Orchid species:") +
  scale_shape_manual(values = c(15, 17, 8, 21), name = "Orchid species:")
dev.off()
shell.exec(tf)

# ------------------------------------------------------------
# 7. Unweighted UniFrac NMDS
# ------------------------------------------------------------

# Evaluate stress across 1–10 dimensions.

stress_vec <- numeric(10)
for(i in seq(10)){
  stress_vec[i] <- vegan::metaMDS(unifrac, weakties=FALSE, k=i)$stress
}

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
plot(seq(10),stress_vec, type = 'o', ylab = "Stress", xlab = "Number of dimensions",
     col="tomato", pch=19)
abline(h=0.2, lty=2) #adds dotted line
abline(h=0.1, lty=3)
dev.off()
shell.exec(tf)

# Fit a five-dimensional NMDS solution.
# Only the first two axes are displayed.
set.seed(1)
NMDS_uni <- vegan::metaMDS(unifrac, weakties=FALSE, trymax=200, k =5)

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
vegan::stressplot(NMDS_uni)  #non-linear fit is R2=1-Stress2, meaning the lower stress the lower this value is
dev.off()
shell.exec(tf)


NMDS2_uni <- NMDS_uni %>%  
  vegan::scores() %>%                                         
  as_tibble(rownames = "Sample")

ps.rarefied_tibble <- as.tibble(ps.rarefied_df) 

nmds_uni <- inner_join(ps.rarefied_tibble, NMDS2_uni, by="Sample") #join tiblles by OTU column so that NMDS scores are added

centroid_uni <- nmds_uni %>%
  group_by(orchid.sp.) %>%
  summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
nmds_uni %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=orchid.sp.)) +
  geom_point(aes(shape = orchid.sp., color = orchid.sp.),size = 4)  +theme_bw() +theme_classic()+ 
  scale_shape_manual(values = c(15, 17, 8, 21), name = "Orchid species:") +
  ggtitle("Unweighted Unifrac") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust=2.5)) +
  theme(text = element_text(size = 15)) + 
  theme(plot.margin=unit(c(1,1,1,1), 'cm')) +
  labs("orchid.sp.") + scale_color_manual(values = c("salmon", "darkgreen", "mediumorchid1", "dodgerblue2"), name = "Orchid species:") +
  geom_point(data=centroid_uni, size=7, shape=21, color="black",
             aes(fill=orchid.sp.), show.legend=FALSE) + scale_fill_manual(values = c("salmon", "darkgreen", "mediumorchid1", "dodgerblue2"))


dev.off()
shell.exec(tf)

# ------------------------------------------------------------
# 8. Unweighted UniFrac PERMDISP by orchid species
# ------------------------------------------------------------

ps.rarefied_metadata <- as(sample_data(ps.rarefied), "data.frame") #export metadata from rarefied phyloseq object 
groups_orchids <- ps.rarefied_metadata[["orchid.sp."]]
unif_permdisp <- vegan::betadisper(unifrac, groups_orchids, type=c("median"))

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
plot(unif_permdisp, main = "Unifrac dist. to centroid", pch=c(15, 17, 19), cex = 0.7, label.cex=0.8)+ points(unif_permdisp$centroids, cex=1.2) 
dev.off()
shell.exec(tf)

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
boxplot(unif_permdisp, main = "", xlab = "") #shows how much variation there is in distances around centroids
dev.off()
shell.exec(tf)

# Test differences in multivariate dispersion.
anova(unif_permdisp)
Tukey_unif <- TukeyHSD(unif_permdisp) #which groups have dif. dispersion?
plot(Tukey_unif)
Tukey <- (as.data.frame(Tukey_unif$group))
write.table(Tukey,"../Tukeyunif.txt",sep="\t",row.names=TRUE)

# ------------------------------------------------------------
# 9. Unweighted UniFrac PERMDISP by phorophyte-sharing status
# ------------------------------------------------------------

groups_shared <- ps.rarefied_metadata[["share.status"]]
unif_permdisp_share <- vegan::betadisper(unifrac, groups_shared, type=c("median"))

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
plot(unif_permdisp, main = "Unifrac dist. to centroid", pch=c(15, 17, 19), cex = 0.7, label.cex=0.8)+ points(unif_permdisp$centroids, cex=1.2) 
dev.off()
shell.exec(tf)

pdf(file = tf <- tempfile(fileext = ".pdf"),width=12,height=10)
boxplot(unif_permdisp, main = "", xlab = "") #shows how much variation there is in distances around centroids
dev.off()
shell.exec(tf)

#Do ANOVA of distances between groups to see if there is heterogeneity of variances/dispersion (p<=0.05)
anova(unif_permdisp_share)
Tukey_unif_share <- TukeyHSD(unif_permdisp_share) #which groups have dif. dispersion?
plot(Tukey_unif)
Tukey2 <- (as.data.frame(Tukey_unif_share$group))
write.table(Tukey2,"../Tukeyunif)sharedstatus.txt",sep="\t",row.names=TRUE)

# ------------------------------------------------------------
# 10. PERMANOVA
# ------------------------------------------------------------

# All models use 9,999 permutations and marginal tests.
# Results should be interpreted alongside PERMDISP tests.

#Model 1: Test the association between phorophyte sharing status and fungal community composition while accounting for orchid species
adonis2_unif <- vegan::adonis2(unifrac ~ share.status+orchid.sp., parallel= T, permutation= 9999, data = ps.rarefied_meta, by = "margin")
adonis2_unif
write.table(adonis2_unif,"../adonissharedstatus_orchid.txt",sep="\t",row.names=TRUE)

# Model 2: Host-sharing status × orchid-species interaction.
adonis2_unif2 <- vegan::adonis2(unifrac ~ share.status*orchid.sp., parallel= T, permutation= 9999, data = ps.rarefied_meta, by = "margin")
adonis2_unif2
write.table(adonis2_unif2,"../adonissharedstatus_orchidb.txt",sep="\t",row.names=TRUE)

# Model 3: Test the association between individual phorophyte identity and fungal community composition while accounting for orchid species.
adonis2_unif3 <- vegan::adonis2(unifrac ~ factor(tree.no.)+orchid.sp., parallel= T, permutation= 9999, data = ps.rarefied_meta, by = "margin")
adonis2_unif3
write.table(adonis2_unif3,"../adonistreeind_orchid.txt",sep="\t",row.names=TRUE)

# Model 4: Test the association between phorophyte species and fungal community composition while accounting for orchid species.
adonis2_unif4 <- vegan::adonis2(unifrac ~ tree.sp.+orchid.sp., parallel= T, permutation= 9999, data = ps.rarefied_meta, by = "margin")
adonis2_unif4
write.table(adonis2_unif4,"../adonistreesp_orchid.txt",sep="\t",row.names=TRUE)

# ============================================================
# End of Script 9
# ============================================================
