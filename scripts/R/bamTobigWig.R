####### 
#@author narumeena
#@description This script converts a BAM file, which is a binary file format that stores nucleotide 
# sequence alignment data, into a bigWig file, which is a compact, indexed, and binary version of a bedGraph file.
#running the script
#Rscript bam2bw.R --bamfile example.bam --outfile example.bw



# Check if the necessary packages are installed and install them if not




######
library(Rsamtools, quietly = TRUE)
library(rtracklayer, quietly = TRUE)
library(GenomicRanges,  quietly = TRUE)
library(GenomicAlignments,  quietly = TRUE)
library(argparse,  quietly = TRUE)

bam2bw <- function(){
  
  # Create an argument parser
  parser <- ArgumentParser()

  # Add the input file argument
  addArgument(parser, "bamfile", help="input BAM file")

  # Add the output file argument
  addArgument(parser, "outfile", help="output bigWig file")

  # Parse the command line arguments
  args <- parseArgs(parser)

  extLen <- 150
  tryCatch({
    cat("opening:", args$bamfile, sep="\n")
    bd <- readGAlignments(args$bamfile)
  }, fileNotFound = function(e) {
    stop("Error: BAM file not found. Please check the file path and try again.")
  })

  cat("convert to GRanges\n")
  mygr <- as(bd,"GRanges")
  cat("extending reads\n")
  mygr <- resize(mygr, extLen)
  mygr <- trim(mygr)

  totalReads <- length(mygr)

  cat("getting coverage\n")
  # get coverage                                                                                                      
  cov <- coverage(mygr)

  rpm <- lapply(cov, function(x) signif(10^6 * x/totalReads,3))
  rpm <- as(rpm,"SimpleRleList")

  
  # create output file if not exists
  if (!file.exists(args$outfile)) {
    file.create(args$outfile)
  }
  
  # export rpm to bigWig                                                                                              
  cat(paste("exporting to bigwig", args$outfile, "\n", sep="\t"))
  export.bw(rpm, args$outfile)
}

