#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(openxlsx))
suppressMessages(library(GenomicRanges))
# Load your custom PICBcombine/PICBannotate/PICBload here if needed
# source("/path/to/PICBcombine.R")

option_list = list(
  make_option("--xlsxs", type="character", help="Input Excel files for replicates", nargs = "+"),
  make_option("--bams", type="character", help="Input BAM files for replicates (optional)", nargs = "+"),
  make_option("--ref", type="character", help="Reference FASTA file"),
  make_option("--out", type="character", help="Output Excel file")
)
opt = parse_args(OptionParser(option_list=option_list))

# Helper: convert Excel data.frame to GRanges, adjust as per your actual Excel column names
df_to_gr <- function(df) {
  # Try to automatically find the appropriate columns
  chr_col <- grep("^chr|^seq", names(df), value=TRUE, ignore.case=TRUE)[1]
  start_col <- grep("^start", names(df), value=TRUE, ignore.case=TRUE)[1]
  end_col <- grep("^end", names(df), value=TRUE, ignore.case=TRUE)[1]
  strand_col <- grep("^strand", names(df), value=TRUE, ignore.case=TRUE)[1]
  if (is.na(chr_col) | is.na(start_col) | is.na(end_col)) stop("Could not find chr/start/end columns!")
  GRanges(
    seqnames = df[[chr_col]],
    ranges = IRanges(start = df[[start_col]], end = df[[end_col]]),
    strand = if (!is.na(strand_col)) df[[strand_col]] else "*"
  )
}

# 1. Read Excel files and convert to rangesList
rangesList = lapply(opt$xlsxs, function(f) {
  df = read.xlsx(f)
  list(clusters = df_to_gr(df))   # Change 'clusters' if your region type differs
})

# 2. (Optional) Read BAMs as alignmentsList (if coverage annotation wanted)
# If you want to run without BAMs, just provide empty list, or implement loading as needed.
alignmentsList = list()
if (!is.null(opt$bams) && length(opt$bams) > 0) {
  # Example: Load as Rsamtools::BamFile objects or your own PICBload function
  # alignmentsList = lapply(opt$bams, function(bam) Rsamtools::BamFile(bam))
  # Names should match replicates if required by your pipeline
  # For now, leave empty unless you want mean coverage annotation
  message("Coverage annotation is not implemented in this template.")
}

# 3. Determine REFERENCE.GENOME string from FASTA, or set manually
# For example, if you use BSgenome, set the correct genome name:
REFERENCE.GENOME = "BSgenome.Mmusculus.UCSC.mm10"  # <-- change as appropriate!

# 4. Call PICBcombine
# TYPE.OF.REGION = "clusters" (change if your Excel region key is different)
combined = PICBcombine(
  rangesList = rangesList,
  alignmentsList = alignmentsList,
  REFERENCE.GENOME = REFERENCE.GENOME,
  TYPE.OF.REGION = "clusters",
  SEQ.LEVELS.STYLE = "UCSC"
)

# 5. Write output as Excel
write.xlsx(as.data.frame(combined), file = opt$out)
message("Written: ", opt$out)
