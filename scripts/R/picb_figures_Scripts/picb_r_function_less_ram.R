# Author: Aleksandr Friman
# Optimized by OpenAI Assistant

# Load necessary libraries
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(GenomicAlignments)
    library(Biostrings)
    library(ggplot2)
    library(parallel)
    library(dplyr)
    library(stringr)
    library(ggpubr)
    library(patchwork)
})

# Function to get names of all reads explained in a batch
namesOfallReadsExplainedBatch <- function(gr, alignments) {
    usedAlignments <- subsetByOverlaps(alignments, gr)
    foundnames <- names(usedAlignments)
    return(unique(foundnames))
}

# Function to get the number of all reads explained in a batch
numberOfallReadsExplainedBatch <- function(gr, alignments) {
    foundnames <- namesOfallReadsExplainedBatch(gr, alignments)
    numbr <- length(foundnames)
    return(numbr)
}

# Optimized forceMap function
forceMap <- function(gr, algn, orderby, numcores, multiplier, maxBatch, insuffix = "") {
    ord1 <- order(mcols(gr)[[orderby]], decreasing = TRUE)
    message("Sorting done")
    sordedGR <- gr[ord1]
    
    if (is.list(algn)) {
        AllAlignments <- c(algn[["unique"]], algn[["multi.primary"]], algn[["multi.secondary"]])
    } else {
        AllAlignments <- algn
    }
    message("AllAlignments ready")
    
    unexplainedAlgn <- AllAlignments
    explainedReads <- numeric(length(sordedGR))
    NumBinsPerBatch <- numcores - 1
    processedBins <- 0
    
    while (processedBins < length(sordedGR)) {
        message("Processing batch")
        
        upperLimit <- min(length(sordedGR), processedBins + NumBinsPerBatch)
        lowerLimit <- processedBins + 1
        message(paste("Limits:", lowerLimit, "-", upperLimit))
        
        BinsToProcess <- sordedGR[lowerLimit:upperLimit]
        message(paste("Number of bins to process:", length(BinsToProcess)))
        message(paste("Alignments remaining:", length(unexplainedAlgn)))
        
        # Find overlaps between the bins and unexplained alignments
        overlaps <- findOverlaps(BinsToProcess, unexplainedAlgn)
        
        # Get the reads explained by these bins
        readsExplained <- unique(names(unexplainedAlgn)[subjectHits(overlaps)])
        
        # Count reads per bin
        binReadCounts <- as.numeric(table(queryHits(overlaps)))
        cumulativeCount <- ifelse(processedBins > 0, explainedReads[processedBins], 0)
        explainedReads[lowerLimit:upperLimit] <- cumulativeCount + cumsum(binReadCounts)
        
        # Remove explained reads from unexplainedAlgn
        unexplainedAlgn <- unexplainedAlgn[!(names(unexplainedAlgn) %in% readsExplained)]
        message("Explained reads removed")
        
        processedBins <- upperLimit
        NumBinsPerBatch <- min(round(multiplier * NumBinsPerBatch), maxBatch)
        gc()
        message("Memory freed")
    }
    
    input.ranges.ordered <- gr[ord1]
    
    if (insuffix != "") {
        insuffix <- paste0("_", insuffix)
    }
    
    readsExplainedDiff <- c(explainedReads[1], diff(explainedReads))
    totalReads <- sum(length(algn$unique), length(algn$multi.primary))
    mcols(input.ranges.ordered)[[paste0("reads_explained", insuffix)]] <- readsExplainedDiff
    mcols(input.ranges.ordered)[[paste0("reads_FPKM", insuffix)]] <- readsExplainedDiff / 
        ((width(input.ranges.ordered) / 1e3) * (totalReads / 1e6))
    
    return(input.ranges.ordered)
}

# Optimized plotDistributionOfAllClusters function
plotDistributionOfAllClusters <- function(inRanges) {
    df2Bnew <- as.data.frame(inRanges)
    maxVal <- max(df2Bnew$width) * 1.1
    DefaultStepVal <- 5000

    # Efficient binning function
    makeBins <- function(inDF, stepVal, maxVal) {
        binEdges <- seq(0, maxVal, by = stepVal)
        inDF$bin <- cut(inDF$width, breaks = binEdges, include.lowest = TRUE, labels = FALSE)
        binCounts <- as.data.frame(table(inDF$bin))
        binCounts$bin <- as.numeric(as.character(binCounts$Var1))
        binCounts$width <- binEdges[-length(binEdges)] + stepVal / 2
        binCounts$normcount <- binCounts$Freq / nrow(inDF)
        return(binCounts)
    }

    df2BsummaryOrd <- makeBins(df2Bnew, DefaultStepVal, maxVal)
    StepVal1 <- 500
    df2BsummaryOrd500 <- makeBins(df2Bnew, StepVal1, maxVal)
    
    # Plotting code (adjust as needed for your specific plot)
    # ...

    # Function to calculate the maximum value for zoomed histograms
    maxValZoom <- function(inDF) {
        maxIndex <- which.max(inDF$Freq)
        maxIndex <- min(maxIndex + 5, nrow(inDF))
        maxVal <- inDF$width[maxIndex]
        return(maxVal)
    }

    Max.Val.AllCl.Zoom <- maxValZoom(df2BsummaryOrd500)
    Max.Val.AllCl.Zoom <- ceiling(Max.Val.AllCl.Zoom / 1000)

    # Create the plots (adjust according to your plotting needs)
    # ...

    return((histAllClustersZoom / histAllClusters / scatterAllClusters) + 
           plot_layout(heights = c(0.477, 0.6, 1.79)))
}

# Optimized clusterCompositionViolinsAndCorrelations function
clusterCompositionViolinsAndCorrelations <- function(inRanges, inGenomeName, maxKmerLen, numCores) {
    # Check if using BSgenome or fasta file
    useBSGenome <- startsWith(inGenomeName, "BSgenome")
    
    if (useBSGenome) {
        message("Using BSgenome")
        inBSgenome <- get(inGenomeName, envir = asNamespace("BSgenome"))
        inRanges$seq <- getSeq(inBSgenome, inRanges, as.character = TRUE)
    } else {
        message("Using fasta file")
        myFAfile <- FaFile(inGenomeName)
        open(myFAfile)
        inRanges$seq <- getSeq(myFAfile, inRanges, as.character = TRUE)
        close(myFAfile)
    }
    
    message("Done loading sequences")
    
    compALL <- addCompositionToGRanges(inRanges, maxKmerLen, numCores)
    visDFkmersALL <- compALL$visDFkmers
    kmersDFALL <- compALL$kmersDF
    
    UpperComposition <- ggplot(visDFkmersALL, aes(
        x = factor(kmer, level = c("A/T", "C/G", 'A', "T", "C", "G")),
        y = score * 100,
        fill = factor(kmer, level = c("A/T", "C/G", 'A', "T", "C", "G"))
    )) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.1, outlier.size = 0.25) +
        theme_classic() +
        ylab("Piwi-piRNA CL (%)") +
        xlab("") +
        scale_fill_manual(values = c("grey", "grey", "green", "red", "blue", "#fff200")) +
        theme(
            legend.position = "none", 
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            axis.text.x = element_text(size = 8), 
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 8),
            axis.title = element_text(size = 8)
        )
    
    corcoefs.Visualization.df <- makeViskmerDF(kmersDFALL, inRanges, "All Clusters")
    
    for (t in sort(unique(inRanges$type))) {
        inRangesType <- inRanges[inRanges$type == t]
        compType <- addCompositionToGRanges(inRangesType, maxKmerLen, numCores)
        kmersDFType <- compType$kmersDF
        corcoefs.Visualization.df.tmp <- makeViskmerDF(kmersDFType, inRangesType, t)
        corcoefs.Visualization.df <- rbind(corcoefs.Visualization.df, corcoefs.Visualization.df.tmp)
    }
    
    LowerComposition <- ggplot(corcoefs.Visualization.df, aes(x = ordered(kmer, level = c("A/T", "C/G", 'A', "T", "C", "G")))) +
        geom_text(aes(y = ordered(type, level = c("SingleCore", 'ExtendedCore', "MultiCore", "All Clusters")), label = corcoef), size = 8 * 5 / 14) +
        geom_tile(aes(y = ordered(type, level = c("SingleCore", 'ExtendedCore', "MultiCore", "All Clusters"))), fill = "white", alpha = .01, color = "black") +
        labs(y = "", x = NULL) +
        scale_y_discrete(labels = c("SCCL", 'ECCL', "MCCL", "All CL")) +
        theme_minimal() +
        theme(
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            strip.text = element_blank(),
            legend.position = "none",
            plot.margin = margin(t = -0.35, r = 0, b = 0, l = 0, unit = "cm"),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 8)
        )
    
    return(list(UpperComposition, LowerComposition))
}

# Helper functions for clusterCompositionViolinsAndCorrelations
FindKmerInRegions <- function(inkmer, inseqvector, inlenvector) {
    tmpcounts <- vapply(seq_along(inseqvector), function(i) {
        str_count(inseqvector[i], fixed(inkmer))
    }, integer(1))
    normcounts <- tmpcounts / inlenvector
    return(normcounts)
}

addCompositionToGRanges <- function(inRanges, maxKmerLen, numCores) {
    bases <- c('A', 'T', 'G', 'C')
    kmersDF <- data.frame()
    
    for (kmerlen in 1:maxKmerLen) {
        tmpkmers <- apply(expand.grid(rep(list(bases), kmerlen)), 1, paste0, collapse = "")
        tmpdf <- mclapply(tmpkmers, FindKmerInRegions, inseqvector = inRanges$seq, inlenvector = width(inRanges), mc.cores = numCores)
        tmpdf <- do.call(cbind, tmpdf)
        colnames(tmpdf) <- tmpkmers
        kmersDF <- cbind(kmersDF, tmpdf)
    }
    
    kmerRegionsDF <- cbind(as.data.frame(inRanges), kmersDF)
    
    # Prepare data for visualization
    visDFkmers <- data.frame(
        kmer = rep(c('A', 'T', 'C', 'G', 'A/T', 'C/G'), each = nrow(kmerRegionsDF)),
        score = c(kmerRegionsDF$A, kmerRegionsDF$T, kmerRegionsDF$C, kmerRegionsDF$G, 
                  kmerRegionsDF$A + kmerRegionsDF$T, kmerRegionsDF$C + kmerRegionsDF$G)
    )
    
    outputList <- list(kmersDF = kmersDF, visDFkmers = visDFkmers, kmerRegionsDF = kmerRegionsDF)
    return(outputList)
}

makeViskmerDF <- function(kmersDF, inRanges, inType) {
    kmersDF.Visualization <- data.frame(
        `A/T` = kmersDF$A + kmersDF$T,
        `C/G` = kmersDF$C + kmersDF$G,
        A = kmersDF$A,
        T = kmersDF$T,
        C = kmersDF$C,
        G = kmersDF$G
    )
    
    corcoefs.Visualization <- sapply(kmersDF.Visualization, cor, y = inRanges$reads_FPKM)
    corcoefs.Visualization.df <- data.frame(
        kmer = names(corcoefs.Visualization),
        corcoef = round(corcoefs.Visualization, digits = 3),
        type = inType
    )
    
    return(corcoefs.Visualization.df)
}

# Note: Include any additional functions you need, optimized similarly.

# Remember to adjust the plotting code according to your specific requirements.

