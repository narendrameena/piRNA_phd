# Load necessary libraries with logging
message("[", Sys.time(), "] Loading required libraries...")
library(PICB)
library(Biostrings)
library(GenomicRanges)
library(GenomeInfoDb)
library(seqinr)
message("[", Sys.time(), "] Libraries loaded successfully.")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 3) {
    stop("Usage: Rscript script.R <FASTA file> <BAM file> <Output Excel file>")
}

# Assign arguments to variables with logging
FASTA.NAME <- args[1]
BAM.NAME <- args[2]
OUTPUT.XLSX <- args[3]
message("[", Sys.time(), "] Input FASTA: ", FASTA.NAME)
message("[", Sys.time(), "] Input BAM: ", BAM.NAME)
message("[", Sys.time(), "] Output Excel: ", OUTPUT.XLSX)

# Function to load the FASTA file and return sequence information with logging
PICBloadfasta <- function(FASTA.NAME = NULL) {
    message("[", Sys.time(), "] Loading FASTA file: ", FASTA.NAME)
    if (is.null(FASTA.NAME)) {
        stop("Please provide FASTA.NAME!")
    }
    FAdata <- seqinr::read.fasta(FASTA.NAME)
    FAnames <- names(FAdata)
    FAlengths <- lengths(FAdata)
    seqinfo <- GenomeInfoDb::Seqinfo(seqnames = FAnames, seqlengths = FAlengths)
    message("[", Sys.time(), "] FASTA file loaded successfully. Number of sequences: ", length(seqinfo))
    return(seqinfo)
}

# Load the genome using PICBloadfasta with logging
message("[", Sys.time(), "] Loading genome from FASTA...")
mygenome <- PICBloadfasta(FASTA.NAME)
message("[", Sys.time(), "] Genome loaded successfully.")

# Load the alignments using PICBload with logging
message("[", Sys.time(), "] Loading alignments from BAM...")
myAlignments <- PICBload(BAM.NAME, REFERENCE.GENOME = mygenome, VERBOSE = TRUE, SEQ.LEVELS.STYLE = "NCBI",WHAT=c("flag"))
message("[", Sys.time(), "] Alignments loaded successfully.")

# Set custom sequence levels to UCSC style with logging
message("[", Sys.time(), "] Setting sequence levels style to UCSC...")
seqlevelsStyle(myAlignments$unique) <- "UCSC"
seqlevelsStyle(myAlignments$multi.primary) <- "UCSC"
seqlevelsStyle(myAlignments$multi.secondary) <- "UCSC"
message("[", Sys.time(), "] Sequence levels style set to UCSC.")

# Form clusters using PICBbuild with logging
message("[", Sys.time(), "] Building clusters...")
myClusters <- PICBbuild(myAlignments, REFERENCE.GENOME = mygenome, VERBOSITY = 0, SEQ.LEVELS.STYLE = "UCSC")#,COMPUTE.1U.10A.FRACTIONS = TRUE)
message("[", Sys.time(), "] Clusters built successfully. Number of clusters: ", length(myClusters))

# Export clusters to Excel file with logging
message("[", Sys.time(), "] Exporting clusters to Excel file...")
PICBexporttoexcel(myClusters, OUTPUT.XLSX)
message("[", Sys.time(), "] Clusters exported to ", OUTPUT.XLSX, " successfully.")

# Final message indicating completion
message("[", Sys.time(), "] Script execution completed successfully.")
