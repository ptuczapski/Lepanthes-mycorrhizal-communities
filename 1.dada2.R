# ============================================================
# Script 1. Sequence processing, fungal taxonomic assignment, and
#           construction of the root-associated
#           general fungi (GF) dataset 
#           from demultiplexed PacBio CCS reads
# ============================================================
# Created by Piotr Tymon Tuczapski
#
library("devtools")
library(dada2)

# Input: demultiplexed PacBio CCS FASTQ files in processed_FASTQfiles.
# Intermediate input: primer-trimmed FASTQ files generated externally with BBDuk.
# Outputs: chimera tracking CSV, taxonomy CSVs, ASV table, and representative FASTA.
# IMPORTANT: Update local absolute paths before running.
setwd("../processed_FASTQfiles")

#PacBio reads come in both orientations. You'll need to orient all reads 5'->3'. First create a list of files sorted by file name and assign it to an object, comes from https://benjjneb.github.io/dada2/tutorial.html
sequence_files_list <- 
  sort(
    list.files(
      "./",
      pattern = "fq"
    )
  )

#Split the raw sequencing file names and grab the 2nd part, i.e. barcode IDs
sample.names <- 
  sapply(
    strsplit(basename(sequence_files_list), "asymm."),
    `[`, 2
  )

#Split again and grab only barcode names to get rid of "fq"
sample.names <-
  sapply(
    strsplit(sample.names, ".fq"),
    `[`, 1
  )


#Orient PacBio CCS reads using the forward primer without trimming primer sequences (intended for PacBio CCS data specifically); this doesn't remove primers! We tried this function but it doesn't work
#set up directory for oriented reads
oriented_list <- 
  file.path(
    "../ccs_oriented/",
    paste0(sample.names, "_orient.fastq.gz")
  )

#Orient reads by White Forward primer
ITS1F <- "TAGAGGAAGTAAAAGTCGTAA"

out.oriented <-
  removePrimers(
    sequence_files_list,
    oriented_list,
    primer.fwd = ITS1F,
    max.mismatch = 5,       #with 5 mismatches it still works for 99% of reads, with 6 it starts leaving behind unoriented reads
    trim.fwd = FALSE,
    trim.rev = FALSE,
    orient = TRUE,
    verbose = TRUE
  )
head(out.oriented)

# EXTERNAL STEP: Transfer oriented FASTQ files to the cluster and trim primers using BBDuk.
# Primer-presence filtering excludes improperly primed reads (don't have one primer in correct orientation)
# as it screens for, passes, and trims only literal primer sequences; it is not equivalent
# to the DADA2 de novo chimera detection performed below. Transfer trimmed files back.

#Inspect read quality profiles
plotQualityProfile(sequence_files_list[1:2])

#Learn the Error Rates using error model err

# Directory containing individually demultiplexed,
# oriented, and primer-trimmed FASTQ files returned from the cluster.
path <- "~/fungal_fastqs" # CHANGE ME to the directory containing the fastq files.
orientedntrim_dir <- "/Users/macstudio/Desktop/fungal_fastqs"
orientedntrim_list <- sort(list.files(
  orientedntrim_dir,
  pattern = "_trimmed.fastq", #sample names have format: SAMPLENAME_trimmed.fastq
  full.names = TRUE,
  ignore.case = TRUE
))

# Verify that these files correspond to the sample IDs derived above.
# File naming after BBDuk may differ from the original FASTQ naming.
print(basename(orientedntrim_list))
# REVIEW REQUIRED: Confirm that sample.names matches the number and order
# of orientedntrim_list before assigning names to dereplicated samples.
stopifnot(length(sample.names) == length(orientedntrim_list))

#"The main purpose of learnErrors is that it speeds things up by learning the error rates over an appropriate subset of large datasets"
err <- learnErrors(orientedntrim_list, multithread = FALSE) 

plotErrors(err, nominalQ=TRUE) #visualise the estimated error rates. Points are the observed error rates for each consensus quality score. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. Estimated error rates (black line) should fit to the observed rates (points), and the error rates are expected to drop with increased quality.

#Read in and dereplicate a fastq file - dereplicating amplicon sequences from fastq files (if there are any replicate CCS sequences it will get rid of them)

derep <- derepFastq(orientedntrim_list, verbose = TRUE, qualityType = "FastqQuality")

names(derep) <- sample.names

#Sample infererence. Apply inference algorithm to the filtered and trimmed data. How many reads in how many unique sequences/sequence variants?
dada_out <- dada(              #The dada function takes as input dereplicated amplicon sequencing reads and returns the inferred composition of the sample (or samples).  Put another way, dada removes all sequencing errors to reveal the members of the sequenced community
  derep,
  err = err,
  multithread = FALSE
)

#Construct sequence table (feature table in QIIME2) - amplicon sequence variant table ASV, higher-resolution version of OTU table
seqtab <- makeSequenceTable(dada_out)
dim(seqtab)                 #how many samples (rows) and ASVs (columns) in the matrix? (dimension)

#Remove chimeras. In this step algorithm screens for seqs that can be reconstructed by combining a left-segment and a right-segment from two more abundant "parent" ASV seqs
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #Fraction of reads retained after chimera removal (not percent chimeric).
table(nchar(getSequences(seqtab.nochim)))   #Inspect distribution of sequence lengths (how many seqs of given length)

#Track reads through the pipeline. How many input, denoised, chimera reads?
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dada_out, getN),
               rowSums(seqtab.nochim))

colnames(track) <- c("denoised", "nonchim")
rownames(track) <- sample.names
head(track,20)
write.csv(track, "removedchimeras_data.csv")

#Assign taxonomy. assignTaxonomy function takes as input a set of sequences to be classified and a training set of reference sequences with known taxonomy
unite_database <- "../UNITE/sh_general_release_dynamic_s_04.02.2020_dev.fasta"

taxa <- 
  assignTaxonomy(
    seqtab.nochim,
    unite_database,
    tryRC = TRUE,     #tell dada2 to look also at the reverse-complement orientation
    verbose = TRUE
  )
# Removing sequence rownames for display only
taxa.print <- taxa

rownames(taxa.print) <- NULL
head(taxa.print)

# Cleaning data up and getting taxa counts

write.csv(
  taxa.print,
    "./assigned_taxa.csv" )

class_names <-
  sapply(
    strsplit(as.character(taxa.print[,1]), "_"),
    `[`, 3
  )

order_names <-
  sapply(
    strsplit(as.character(taxa.print[,2]), "_"),
    `[`, 3
  )

family_names <-
  sapply(
    strsplit(as.character(taxa.print[,3]), "_"),
    `[`, 3
  )

genus_names <-
  sapply(
    strsplit(as.character(taxa.print[,4]), "_"),
    `[`, 3
  )
species_names <-
  sapply(
    strsplit(as.character(taxa.print[,5]), "_"),
    `[`, 3
  )

# Generate taxonomic labels for descriptive purposes only.
taxa_names <- paste0(class_names, "_", order_names, "_", family_names, "_", genus_names, "_", species_names)

#optional step
#colnames(seqtab.nochim) <- taxa_names  #substitute colnames in seqtab: was nt sequence, now 'class_order_family_genus_species'

# make a loop: for each sample create a list of fungal taxa when in a row for given sample value is not 0
for(sample in 1:length(seqtab.nochim[,1])){

  temp <- taxa_names[seqtab.nochim[sample,] != 0]
    
  if(length(temp) == 0) {temp = NA}
  
  write.csv(
    temp,
    paste0(
      "./assigned_taxa_", 
      strsplit(rownames(seqtab)[sample], "_")[[1]][1],
      ".csv")
  )
    
}  

# REVIEW REQUIRED: We preserved the unique nucleotide sequence identifiers
# generated by DADA2 as column names of the ASV table.If you use taxonomy-derived IDs assigned as above, they may
# be duplicated and may not correspond to nucleotide sequence identifiers.You should verify before moving on.
#dada2 back to qiime2; see https://forum.qiime2.org/t/dada2-back-to-qiime-2/1316/11

write.table(t(seqtab.nochim), "/Users/macstudio/Desktop/merged_fastqs/seqtab_nochim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE) #save table with all denoised, dereplicated, and removed chimeric seqs; it contains NAs
uniquesToFasta(seqtab.nochim, fout='/Users/macstudio/Desktop/merged_fastqs/rep-seqs_nochim.fna', ids=colnames(seqtab.nochim))

# ============================================================
# End of Script 1
# ============================================================