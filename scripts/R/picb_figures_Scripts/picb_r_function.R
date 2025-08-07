
#author Aleksandr Friman
# Function to get names of all reads explained in a batch
namesOfallReadsExplainedBatch <- function(gr, alignments) {
    usedAlignments <- subsetByOverlaps(alignments, gr)
    foundnames <- names(usedAlignments)
    return(unlist(unique(foundnames)))
}

# Function to get the number of all reads explained in a batch
numberOfallReadsExplainedBatch <- function(gr, alignments) {
    foundnames <- namesOfallReadsExplainedBatch(gr, alignments)
    numbr <- length(foundnames)
    return(numbr)
}

# Function to force map alignments to genomic regions
forceMap <- function(gr, algn, orderby, numcores, multiplier, maxBatch, insuffix = "") {
    ord1 <- order(mcols(gr)[[orderby]], decreasing = TRUE)
    message("Sorting done")
    sordedGR <- gr[ord1]

    if (typeof(algn) == "list") {
        AllAlignments <- c(algn[["unique"]], algn[["multi.primary"]], algn[["multi.secondary"]])
    } else {
        AllAlignments <- algn
    }
    message("AllAlignments ready")

    unexplainedAlgn <- AllAlignments[names(AllAlignments) %in% namesOfallReadsExplainedBatch(gr, AllAlignments)]

    numberOfReads <- function(gr) {
        return(numberOfallReadsExplainedBatch(gr, unexplainedAlgn))
    }

    explainedReads <- c()
    NumBinsPerBatch <- numcores - 1
    processedBins <- 0
    addedCounts <- 0

    while (processedBins < length(sordedGR)) {
        message("Doing batch")
        
        upperLimit <- min(length(sordedGR), (processedBins + 1 + NumBinsPerBatch))
        lowerLimit <- (processedBins + 1)
        message(paste("Limits:", lowerLimit, " - ", upperLimit))

        BinsToProcess <- lapply(lowerLimit:upperLimit, function(idx) sordedGR[1:idx])
        message(paste("Len:", length(BinsToProcess)))
        message(paste("Alignments to do:", length(unexplainedAlgn)))

        explainedReadsTMP <- mclapply(BinsToProcess, numberOfReads, mc.cores = numcores, mc.preschedule = TRUE) # Change mc.cores to the value suitable for your machine. Too high value leads to running out of memory
        explainedReads <- c(explainedReads, unlist(explainedReadsTMP) + addedCounts)
        addedCounts <- explainedReads[length(explainedReads)]
        message("Batch ready")

        # Identifying the explained read names
        theLastBin <- BinsToProcess[[length(BinsToProcess)]]
        readsToRemove <- namesOfallReadsExplainedBatch(theLastBin, unexplainedAlgn)
        unexplainedAlgn <- unexplainedAlgn[!names(unexplainedAlgn) %in% readsToRemove]
        message("Explained reads removed")

        processedBins <- upperLimit
        NumBinsPerBatch <- min(round(multiplier * NumBinsPerBatch), maxBatch)
        rm(BinsToProcess)
        gc()
        message("Memory freed")
    }

    input.ranges.ordered <- gr[ord1]

    if (insuffix != "") {
        insuffix <- paste0("_", insuffix)
    }

    mcols(input.ranges.ordered)[[paste0("reads_explained", insuffix)]] <- c(explainedReads[1], diff(explainedReads))
    mcols(input.ranges.ordered)[[paste0("reads_FPKM", insuffix)]] <- mcols(input.ranges.ordered)[[paste0("reads_explained", insuffix)]] / 
        ((width(input.ranges.ordered) / 1e3) * ((sum(length(algn$unique)) + sum(length(algn$multi.primary))) / 1e6))
    
    return(input.ranges.ordered)
}




#author Aleksandr Friman
plotBarplotsWithErrorBarsFig1B <- function(inLociList, inAlignmentsList, chosenGenome) {
    # Sizes for annotations, titles, and ticks
    annotationSize <- 8
    titleSize <- 8
    ticksSize <- 8

    # Y-axis scales and flips
    YScales <- list(
        Explained = scale_y_continuous(
            breaks = c(0, 0.25, 0.50, 0.75, 1.00),
            labels = c(0, "", 50, "", 100),
            limits = c(0, 1)
        ),
        Number = NULL,
        Width = NULL
    )

    YFlips <- list(
        Explained = coord_flip(),
        Number = scale_y_continuous(n.breaks = 4),
        Width = scale_y_continuous(n.breaks = 3)
    )

    # Y-axis labels
    YLabs <- list(
        Explained = "% of the library",
        Number = "Number (thousands)",
        Width = "Genome (%)"
    )

    # Processing input data
    masterDF <- data.frame(
        matrix(ncol = 5, nrow = 0,
            dimnames = list(NULL, c("type", "replicate", "Explained", "Number", "Width"))
        )
    )

    for (t in c("seeds", "cores", "clusters")) {
        for (replN in 1:length(inLociList)) {
            curAlignments <- inAlignmentsList[[replN]]
            curLoci <- inLociList[[replN]][[t]]
            curSumReads <- sum(curLoci$reads_explained)
            curTotalReads <- length(curAlignments$unique) + length(curAlignments$multi.primary)

            masterDF[nrow(masterDF) + 1, ] <- c(
                t, replN,
                curSumReads / curTotalReads,
                length(curLoci) / 1000,
                sum(width(curLoci)) * 100 / (2 * sum(seqlengths(chosenGenome))) # 2* to account for the minus strand
            )
        }
    }

    # Convert columns to numeric
    for (f in c("Explained", "Number", "Width")) {
        masterDF[[f]] <- as.numeric(masterDF[[f]])
    }

    # Define Y-axis scales for Number and Width
    YScales$Number <- coord_flip(
        ylim = c(0.8 * min(masterDF$Number), ceiling(11 * max(masterDF$Number)) / 10),
        expand = TRUE
    )
    
    YScales$Width <- coord_flip(
        ylim = c(0, 1.1 * max(masterDF$Width)),
        expand = FALSE
    )

    # Prepare the data for plotting
    PlotDF <- data.frame(
        matrix(ncol = 7, nrow = 0,
            dimnames = list(NULL, c("type", "Explained", "Number", "Width", "Explainedsd", "Numbersd", "Widthsd"))
        )
    )

    for (t in c('seeds', 'cores', 'clusters')) {
        PlotDF[nrow(PlotDF) + 1, ] <- c(t, NA, NA, NA, NA, NA, NA)
        
        for (f in c("Explained", "Number", "Width")) {
            myval <- mean(masterDF[[f]][masterDF$type == t])
            mysd <- sd(masterDF[[f]][masterDF$type == t])
            
            PlotDF[[f]][PlotDF[["type"]] == t] <- as.numeric(myval)
            PlotDF[[paste0(f, "sd")]][PlotDF[["type"]] == t] <- as.numeric(mysd)
        }
    }

    # Create a list to store the figures
    FigList <- list()

    # Generate bar plots with error bars
    for (f in c("Explained", "Width", "Number")) {
        fsd <- paste0(f, "sd")
        
        PlotDF[[f]] <- as.numeric(PlotDF[[f]])
        PlotDF[[fsd]] <- as.numeric(PlotDF[[fsd]])
        
        FigList[[f]] <- ggplot(PlotDF, aes(x = factor(type, level = c('clusters', 'cores', 'seeds')), y = !!sym(f), fill = type)) +
            geom_bar(stat = "identity", color = "black", position = position_dodge(), show.legend = FALSE) +
            geom_errorbar(aes(ymin = !!sym(f) - !!sym(fsd), ymax = !!sym(f) + !!sym(fsd)), width = .8, position = position_dodge(.9), linewidth = 0.5) +
            scale_fill_manual(values = c(color.clusters, color.cores, color.seeds)) +
            theme(
                axis.ticks.length = unit(0.055, "in"),
                axis.ticks = element_line(colour = "black", linewidth = 0.5),
                axis.line.x = element_line(colour = "black", linewidth = 0.5)
            ) +
            YScales[[f]] +
            theme(legend.position = "none") + 
            xlab("") + 
            ggtitle(YLabs[[f]]) + 
            ylab("") +
            theme(
                plot.title.position = 'plot',
                plot.title = element_text(hjust = 0.5, size = titleSize),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text.y = element_text(size = 0),
                axis.text.x = element_text(size = ticksSize),
                axis.ticks.y = element_blank()
            ) +
            YFlips[[f]] +
            theme(
                panel.background = element_rect(fill = 'transparent'),
                plot.background = element_rect(fill = 'transparent', color = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.margin = unit(c(0, 0, 0, 0), "cm")
            )
    }

    # Arrange the plots in a grid and return
    return(
        ggarrange(
            NULL, FigList[["Number"]], NULL, FigList[["Width"]], NULL,
            FigList[["Explained"]], NULL,
            nrow = 1,
            widths = c(0.2, 1, 0.1, 0.5, 0.1, 1, 0.2)
        )
    )
}





#Author Aleksandr Friman
barPlotOfClustersByTypehoriz <- function(inRanges) {
  df2Anew <- as.data.frame(inRanges)
  df2Apie <- data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("clustertype", "value", "valuetype"))))
  
  for (t in sort(unique(inRanges$type))) {
    mycount <- nrow(df2Anew[df2Anew$type == t,])
    mysumreads <- sum(unlist(df2Anew$reads_explained[df2Anew$type == t]))
    mysumlen <- sum(unlist(df2Anew$width_in_nt[df2Anew$type == t]))
    
    df2Apie[nrow(df2Apie) + 1, ] <- c(t, as.numeric(mycount), "Number of\nclusters")
    df2Apie[nrow(df2Apie) + 1, ] <- c(t, as.numeric(mysumreads), "Reads\nexplained")
    df2Apie[nrow(df2Apie) + 1, ] <- c(t, as.numeric(mysumlen), "Length")
  }
  
  df2Apie$value <- as.numeric(df2Apie$value)
  
  for (myvaltype in unique(df2Apie$valuetype)) {
    message(df2Apie$value[df2Apie$valuetype == myvaltype])
    normCoef <- sum(df2Apie$value[df2Apie$valuetype == myvaltype])
    message(normCoef)
    df2Apie$value[df2Apie$valuetype == myvaltype] <- df2Apie$value[df2Apie$valuetype == myvaltype] / normCoef
  }
  
  message(df2Apie)
  
  colorLevels <- c(color.clusters, color.eC, color.cores) # Ensure these color variables are defined: color.clusters, color.eC, color.cores
  
  pHoriz <- ggplot(df2Apie, aes(x = factor(valuetype, levels = rev(c("Number of\nclusters", "Length", "Reads\nexplained")), ordered = TRUE),
                                y = value, fill = factor(clustertype, level = c("MultiCore", "ExtendedCore", "SingleCore"), ordered = TRUE))) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(mult = c(0, .0))) +
    theme(axis.line.y = element_blank(), 
          axis.line.x = element_line(),
          panel.background = element_blank(), 
          axis.ticks.y = element_blank()) +
    scale_fill_manual(values = colorLevels) + 
    theme(legend.position = "none", plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
    theme(axis.text.y = element_text(size = 8), 
          axis.text.x = element_text(size = 8),
          axis.title = element_blank()) +
    coord_flip()
  
  return(pHoriz)
}




#Author Aleksandr Friman
perChrClustersAndProductivity <- function(inRanges, inAlignments, inGenome) {
  # Convert inGenome to GRanges
  GenomeGR <- GRanges(inGenome)
  
  # Calculate total reads
  TotalReads <- length(inAlignments$unique) + length(inAlignments$multi.primary)
  
  # Initialize dataframe for Genome Space
  GSdf <- data.frame(matrix(ncol = 5, nrow = 0, dimnames = list(NULL, c("chr", "strand", "widthFraction", "LibraryFraction", "type"))))
  
  # Loop over strands, chromosomes, and types
  for (s in c("+", "-")) {
    for (mychr in as.vector(seqnames(GenomeGR))) {
      for (mytype in c("SingleCore", "ExtendedCore", "MultiCore")) {
        chrWidth <- width(GenomeGR[seqnames(GenomeGR) == mychr])
        subsetGR <- inRanges[(seqnames(inRanges) == mychr) & (strand(inRanges) == s) & (inRanges$type == mytype)]
        GSdf[nrow(GSdf) + 1, ] <- c(gsub("chr", "", mychr), s, sum(width(subsetGR)) / chrWidth, sum(subsetGR$reads_explained) / TotalReads, mytype)
      }
    }
  }
  
  # Convert columns to numeric
  GSdf$widthFraction <- as.numeric(GSdf$widthFraction)
  GSdf$LibraryFraction <- as.numeric(GSdf$LibraryFraction)
  
  # Define chromosome order
  uniqCHR <- sort(unique(GSdf$chr))
  
  if ("X" %in% GSdf$chr) {
    orderofChr <- c("X", uniqCHR[!uniqCHR %in% "X"])
  } else {
    orderofChr <- uniqCHR
  }
  
  NumbericCHRs <- uniqCHR[!is.na(as.numeric(uniqCHR))]
  NumbericCHRs <- as.character(sort(as.numeric(NumbericCHRs)))
  LetterCHRs <- uniqCHR[is.na(as.numeric(uniqCHR))]
  orderofChr <- c(LetterCHRs, NumbericCHRs)
  
  # Define color levels
  colorLevels <- c(color.cores, color.eC, color.clusters)  # Ensure these color variables are defined: color.cores, color.eC, color.clusters
  
  # Create Genome Space plot
  WidthPlot <- ggplot(GSdf, aes(x = factor(chr, levels = rev(orderofChr)), fill = factor(type, levels = c("SingleCore", "ExtendedCore", "MultiCore")))) +
    geom_bar(data = subset(GSdf, strand == "+"), aes(y = widthFraction * 100, fill = factor(type, levels = c("SingleCore", "ExtendedCore", "MultiCore"))), stat = "identity") +
    geom_bar(data = subset(GSdf, strand == "-"), aes(y = -widthFraction * 100, fill = factor(type, levels = c("SingleCore", "ExtendedCore", "MultiCore"))), stat = "identity") +
    scale_y_continuous(labels = abs) +
    geom_hline(yintercept = 0, colour = "grey90") +
    scale_fill_manual(values = colorLevels) +
    theme_classic() +
    xlab("") +
    ylab("Genome space (%)") +
    coord_flip() +
    theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, -0.5), "cm")) +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8))
  
  # Create Library Fraction plot
  LibraryPlot <- ggplot(GSdf, aes(x = factor(chr, levels = rev(orderofChr)))) +
    geom_bar(data = subset(GSdf, strand == "+"), aes(y = LibraryFraction * 100, fill = strand), stat = "identity", position = "dodge") +
    geom_bar(data = subset(GSdf, strand == "-"), aes(y = -LibraryFraction * 100, fill = strand), stat = "identity", position = "dodge") +
    scale_y_continuous(labels = abs) +
    geom_hline(yintercept = 0, colour = "grey90") +
    scale_fill_manual(values = c("blue", "red")) +
    theme_classic() +
    xlab("") +
    ylab("piRNAs (%)") +
    coord_flip() +
    theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, -0.2), "cm")) +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8))
  
  # Arrange and return the plots
  return(ggarrange(WidthPlot, LibraryPlot, nrow = 1))
}



FPKMCumsumOfReadsWithExlReadsPercentilesClustersHL <- function(inRanges, inAlignemnts, cuts, inFPKMfactor, inHLRanges, perc = c(0.9, 0.95)) {
  totalReads <- length(inAlignemnts$unique) + length(inAlignemnts$multi.primary)
  explainedReads <- sum(inRanges$clusters$reads_explained)
  GRsorted <- inRanges$clusters[order(inRanges$clusters$reads_explained, decreasing = TRUE)]
  GRsorted$rank <- 1:length(GRsorted)
  
  # Overlaps with HL ranges
  colorMap <- list(productivity = color.clusters, FPKM = color.FPKM)
  GRsorted$varname <- "productivity"
  
  if (length(inHLRanges) > 0) {
    for (scind in 1:length(inHLRanges)) {
      sc <- inHLRanges[scind]
      scName <- sc[1]$names
      scColor <- sc[1]$color
      GRsorted$specialOverlap <- countOverlaps(GRsorted, sc)
      GRsorted$varname[GRsorted$specialOverlap > 0] <- scName
      GRsorted$specialOverlap <- NULL
      colorMap[[scName]] <- scColor
    }
  }
  
  # End overlaps with HL
  dfplotmelted1 <- data.frame(
    value = cumsum(GRsorted$reads_explained) * 100 / totalReads,
    variable = GRsorted$varname,
    number = GRsorted$rank
  )
  
  scalingCoefficientFPKM <- max(GRsorted$reads_FPKM) / 100 # To make sure both axes are on the same scale
  
  if (inFPKMfactor > 0) {
    dfplotmelted2 <- data.frame(
      value = GRsorted$reads_FPKM / scalingCoefficientFPKM,
      variable = "FPKM",
      number = 1:length(GRsorted)
    )
    dfplotmelted <- rbind(dfplotmelted2, dfplotmelted1)
  } else {
    dfplotmelted <- dfplotmelted1
  }
  
  p1 <- ggplot(dfplotmelted, aes(x = number, y = value)) +
    geom_point(aes(color = factor(variable)), size = 0.25) +
    xlab(paste0("Piwi-piC (1-", length(GRsorted), ") (ranked)")) +
    theme_classic() +
    scale_colour_manual(values = colorMap, name = "") +
    scale_y_continuous(
      name = "Reads, % of library", 
      breaks = seq(0, 100, 10), 
      limits = c(0, 100),
      sec.axis = sec_axis(trans = ~ . * scalingCoefficientFPKM / inFPKMfactor, 
                          name = paste0("piRNAs (x ", inFPKMfactor, ") (FPKM)"))
    ) + 
    theme(
      legend.position = c(0.1, 0.1),
      axis.text.x.bottom = element_text(size = 8), 
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      legend.text = element_text(size = 8), 
      axis.title.x.top = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x.top = element_blank(), 
      legend.title = element_blank()
    )
  
  # Adding explained reads
  maxCut <- 1
  if (length(cuts) > 0) {
    maxCut <- max(unlist(cuts))
  }
  
  if (length(perc) > 0) {
    perc <- perc[order(perc, decreasing = TRUE)]
    expl100 <- 100 * explainedReads / totalReads
    textEndX <- length(GRsorted) - 0.1 * (length(GRsorted) - maxCut)
    dataToFindClusterOfPerc <- cumsum(GRsorted$reads_explained) * 100 / totalReads
    textEndX <- textEndX - 0.2 * (length(GRsorted) - maxCut)
    
    for (targetperc in perc) {
      textEndX <- textEndX - 0.1 * (length(GRsorted) - maxCut)
      targetCluster <- min(which(dataToFindClusterOfPerc >= targetperc * expl100))
      message(paste(targetperc, targetCluster))
      
      p1 <- p1 +
        annotate("text", x = targetCluster, vjust = -0.3, y = 0.50 * targetperc * expl100, color = "red",
                 label = TeX(paste("$", targetperc * 100, "^{th}$ percentile")), angle = 90, size = 8 / .pt) +
        annotate("text", x = targetCluster, vjust = 1, y = 0, color = "red", label = targetCluster, size = 8 / .pt) +
        annotate("segment", x = targetCluster, xend = targetCluster, y = targetperc * expl100, yend = 0, color = "red", linetype = 2)
    }
  }
  
  for (myCut in cuts) {
    cutScale <- abs(diff(myCut)) / length(GRsorted)
    p1 <- p1 +
      scale_x_break(unlist(myCut), scales = "fixed", expand = expansion(add = c(1, 0)), space = 0.178 * 2.54) +
      geom_vline(xintercept = myCut, linetype = "dotted")
  }
  
  mybreaks <- unique(sort(c(1, length(GRsorted), unlist(cuts))))
  
  if (length(cuts) > 0) {
    p1 <- p1 + scale_x_continuous(breaks = mybreaks, expand = c(0.01, 0.01))
  }
  
  return(p1)
}





# Authors: Zuzana Loubalova, Aleksandr Friman

strandFlip <- function(chosenAlignments) {
  alignmentsFLIP <- chosenAlignments
  
  for (t in names(chosenAlignments)) {
    GenomicRanges::strand(alignmentsFLIP[[t]]) <- ifelse(
      GenomicRanges::strand(alignmentsFLIP[[t]]) == "+", 
      "-", 
      "+"
    )
  }
  
  return(alignmentsFLIP)
}

getStrandeness <- function(chosenLoci, chosenAlignments, chosenGenome) {
  chosenAlignmentsFlip <- strandFlip(chosenAlignments)
  
  chosenLoci$clusters <- PICBannotate(
    INPUT.GRANGES = chosenLoci$clusters,
    ALIGNMENTS = chosenAlignmentsFlip,
    REPLICATE.NAME = "AntiSense",
    REFERENCE.GENOME = chosenGenome,
    PROVIDE.NON.NORMALIZED = TRUE
  )
  
  chosenLoci$clusters$Unique.S_AS_ratio_log10 <- log10(
    (chosenLoci$clusters$uniq_reads + 1) / 
    (chosenLoci$clusters$uniq_reads.AntiSense + 1)
  )
  
  return(chosenLoci)
}


# Author: Aleksandr Friman

getStrandnessPlot <- function(inRanges, sampleName, coln = "S_AS_ratio_log10", weightsCol = "reads_explained") {
  
  # Sort inRanges by the specified weights column in descending order
  inRangesSorted <- inRanges[order(mcols(inRanges)[[weightsCol]], decreasing = TRUE)]
  
  # Set up the weights for ggplot
  inRangesSorted$weightsForGGplot <- mcols(inRangesSorted)[[weightsCol]]
  
  # Define Y-axis ticks and labels for the plot
  Yticks <- floor(min(mcols(inRangesSorted)[[coln]])):ceiling(max(mcols(inRangesSorted)[[coln]]))
  Ylabels <- TeX(paste0("$10^{", Yticks, "}$"))
  
  # Set up the visual column for plotting
  inRangesSorted$visual <- mcols(inRangesSorted)[[coln]]
  
  # Create the ggplot object
  p2 <- ggplot(inRangesSorted, aes(x = sampleName, y = visual)) +
    geom_violin(fill = color.raw, color = "#58595B") +
    geom_violin(aes(weights = weightsForGGplot), fill = color.weigthed, color = color.weigthed, alpha = 0.3) +
    theme_classic() +
    annotation_logticks(sides = "l") +
    scale_y_continuous(breaks = Yticks, labels = Ylabels) +
    ylab("Sense/Antisense ratio") +
    xlab("") +
    annotate("rect", ymin = -1, ymax = 1, xmax = Inf, xmin = -Inf, alpha = 0.15, fill = color.doublestrandness) +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
    theme(axis.title = element_text(size = 8))
  
  return(p2@ggplot)
}



# Author: Aleksandr Friman

plotDistributionOfAllClusters <- function(inRanges) {
  df2Bnew <- as.data.frame(inRanges)
  maxVal <- max(df2Bnew$width) * 1.1
  DefaultStepVal <- 5000
  
  df2BsummaryOrd <- makeBins(df2Bnew, DefaultStepVal, maxVal)
  
  StepVal1 <- 500
  df2BsummaryOrd500 <- makeBins(df2Bnew, StepVal1, maxVal)
  
  StepVal2 <- 1000
  df2BsummaryOrd1000 <- makeBins(df2Bnew, StepVal2, maxVal)
  
  # Function to calculate the maximum value for zoomed histograms
  maxValZoom <- function(inDF) {
    maxIndex <- which(inDF$count == max(inDF$count))
    maxIndex <- min(maxIndex + 5, length(inDF$count))
    maxVal <- inDF$width[maxIndex]
    return(maxVal)
  }
  
  myY.Max <- ceiling(log10(max(inRanges$reads_FPKM)))
  myY.Min <- floor(log10(min(inRanges$reads_FPKM)))
  myY.Breaks <- myY.Min:myY.Max
  
  scatterAllClusters <- ggplot(df2Bnew, aes(x = width / 1000, y = log10(reads_FPKM), color = color.clusters)) +
    geom_point(size = 0.25) +
    xlab("Length, kb")
  
  if (max(df2Bnew$width) > 30000) {
    scatterAllClusters <- scatterAllClusters + scale_x_continuous(breaks = breaks_width(10 * ceiling(max(df2Bnew$width / 1000) / 30)))
  }
  
  scatterAllClusters <- scatterAllClusters +
    ylab("piRNA production (FPKM)") +
    scale_color_manual(values = c(color.clusters), labels = c("Cluster")) +
    theme_classic() +
    themeAddon +
    scale_y_continuous(breaks = myY.Breaks, labels = 10^myY.Breaks, limits = c(myY.Min, myY.Max)) +
    annotation_logticks(sides = "l", short = unit(0.016, "in"), mid = unit(0.032, "in"), long = unit(0.05, "in")) +
    theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
    xlim(0, NA)
  
  scatterAllClusters.built <- ggplot_build(scatterAllClusters)
  scatterAllClusters.raw.Xticks <- scatterAllClusters.built$layout$panel_params[[1]]$x$breaks # nolint
  scatterAllClusters.Xticks <- scatterAllClusters.raw.Xticks[!is.na(scatterAllClusters.raw.Xticks)]
  
  scatterAllClusters <- scatterAllClusters + 
    annotate("text", x = 0.9 * max(scatterAllClusters.Xticks), y = 0.95 * max(myY.Breaks), 
             label = paste0("r=", round(cor(df2Bnew$width, df2Bnew$reads_FPKM), digits = 3)), 
             size = 8 / .pt, color = "black", family = "Helvetica")
  
  # Histograms
  histAllClusters <- ggplot(data = df2BsummaryOrd[df2BsummaryOrd$count > 0,], aes(x = width / 1000, y = normcount, fill = color.cores)) +
    geom_bar(stat = "identity") + 
    scale_x_continuous(breaks = scatterAllClusters.Xticks) +
    xlab("Length, kb") + 
    ylab("Frequency") +
    scale_fill_manual(name = "Type", values = c(color.clusters)) +
    theme_classic() +
    themeAddon +
    theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  
  Max.Val.AllCl.Zoom <- maxValZoom(df2BsummaryOrd)
  Max.Val.AllCl.Zoom <- ceiling(Max.Val.AllCl.Zoom / 1000)
  
  histAllClustersZoom <- ggplot(data = df2BsummaryOrd500, aes(x = width / 1000, y = normcount, fill = color.cores)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, Max.Val.AllCl.Zoom)) +
    scale_fill_manual(name = "Type", values = c(color.clusters)) +
    xlab("Length, kb") + 
    ylab("Frequency") +
    theme_classic() +
    themeAddon +
    theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_blank(), axis.title.x = element_text(size = 8)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  
  return((histAllClustersZoom / histAllClusters / scatterAllClusters) + plot_layout(heights = c(0.477, 0.6, 1.79)))
}

# Function to create bins for histograms
makeBins <- function(inDF, stepVal, maxVal) {
  clusterTypes <- sort(unique(inDF$type))
  vals2B <- seq(0, maxVal + stepVal / 2, stepVal)
  
  df2Bsummary <- data.frame(
    matrix(
      ncol = 5, 
      nrow = 0,
      dimnames = list(NULL, c("startW", "stopW", "type", "count", "normcount"))
    )
  )
  
  for (i in 1:(length(vals2B) - 1)) {
    startW <- vals2B[i]
    stopW <- vals2B[i + 1]
    
    for (t in clusterTypes) {
      mycount <- nrow(inDF[(inDF$width > startW) & (inDF$width <= stopW) & (inDF$type == t),])
      mynormcount <- mycount / nrow(inDF[(inDF$type == t),])
      df2Bsummary[nrow(df2Bsummary) + 1, ] <- c(startW, stopW, t, mycount, mynormcount)
    }
  }
  
  df2Bsummary$width <- (as.numeric(df2Bsummary$stopW) + as.numeric(df2Bsummary$startW)) / 2
  df2Bsummary$normcount <- as.numeric(df2Bsummary$normcount)
  df2Bsummary$startW <- as.numeric(df2Bsummary$startW)
  df2Bsummary$stopW <- as.numeric(df2Bsummary$stopW)
  df2Bsummary$color <- NA
  df2Bsummary$color[df2Bsummary$type == clusterTypes[1]] <- color.eC
  df2Bsummary$color[df2Bsummary$type == clusterTypes[2]] <- color.clusters
  df2Bsummary$color[df2Bsummary$type == clusterTypes[3]] <- color.cores
  
  order2Bnew <- order(df2Bsummary$type, decreasing = TRUE)
  df2BsummaryOrd <- df2Bsummary[order2Bnew,]
  df2BsummaryOrd$count <- as.numeric(df2BsummaryOrd$count)
  
  return(df2BsummaryOrd)
}


# Author: Aleksandr Friman

clusterCompositionViolinsAndCorrelations <- function(inRanges, inGenomeName, maxKmerLen, numCores) {
  BSgenomeTestSTR <- "BSgenome"
  useBSGenome <- FALSE
  
  if (nchar(inGenomeName) > nchar(BSgenomeTestSTR)) {
    message("DEBUG1")
    message(substr(inGenomeName, 1, nchar(BSgenomeTestSTR)))
    
    if (substr(inGenomeName, 1, nchar(BSgenomeTestSTR)) == BSgenomeTestSTR) {
      message("DEBUG2")
      useBSGenome <- TRUE
    }
  }
  
  if (useBSGenome) {
    message("Using BSgenome")
    inBSgenome <- getBSgenome(inGenomeName)
    inRanges$seq <- getSeq(inBSgenome, inRanges, as.character = TRUE)
  } else {
    message("Using fasta file")
    myFAfile <- Rsamtools::FaFile(inGenomeName)
    mySeq <- as.character(getSeq(myFAfile, inRanges, as.character = TRUE))
    inRanges$seq <- as.character(mySeq)
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
    ) +
    themeAddon
  
  corcoefs.Visualization.df <- makeViskmerDF(kmersDFALL, inRanges, "All Clusters")
  
  for (t in sort(unique(inRanges$type))) {
    inRangesType <- inRanges[inRanges$type == t]
    compType <- addCompositionToGRanges(inRangesType, maxKmerLen, numCores)
    visDFkmersType <- compType$visDFkmers
    kmersDFType <- compType$kmersDF
    corcoefs.Visualization.df.tmp <- makeViskmerDF(kmersDFType, inRangesType, t)
    corcoefs.Visualization.df <- rbind(corcoefs.Visualization.df, corcoefs.Visualization.df.tmp)
  }
  
  LowerComposition <- corcoefs.Visualization.df %>%
    ggplot(aes(x = ordered(kmer, level = c("A/T", "C/G", 'A', "T", "C", "G")))) +
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
      axis.title = element_text(size = 8),
      text = element_text(size = 8)
    ) +
    themeAddon
  
  return(list(UpperComposition, LowerComposition))
}

FindKmerInRegions <- function(inkmer, inseqvector, inlenvector) {
  tmpcounts <- unlist(lapply(str_locate_all(inseqvector, inkmer), nrow))
  normcounts <- tmpcounts / inlenvector
  return(normcounts)
}

mypaste <- function(inputV) {
  return(paste0(inputV, collapse = ""))
}

addCompositionToGRanges <- function(inRanges, maxKmerLen, numCores) {
  bases <- c('A', 'T', 'G', 'C')
  kmersDF <- data.frame()
  
  for (kmerlen in 1:maxKmerLen) {
    tmpperm <- permutations(n = length(bases), r = kmerlen, v = bases, repeats.allowed = TRUE)
    tmpkmers <- apply(tmpperm, 1, mypaste)
    tmpdf <- as.data.frame(mclapply(tmpkmers, FindKmerInRegions, inseqvector = inRanges$seq, inlenvector = width(inRanges), mc.cores = numCores))
    colnames(tmpdf) <- tmpkmers
    
    if (nrow(kmersDF) == 0) {
      kmersDF <- tmpdf
    } else {
      kmersDF <- cbind(kmersDF, tmpdf)
    }
  }
  
  kmerRegionsDF <- as.data.frame(inRanges)
  kmerRegionsDF <- cbind(kmerRegionsDF, kmersDF)
  
  # Filtering the main bases
  visDFkmers <- data.frame(kmer = NULL, score = NULL)
  
  # Add A
  visDFkmersTMP <- data.frame(kmer = rep("A", length(kmerRegionsDF$A)), score = kmerRegionsDF$A)
  visDFkmers <- rbind(visDFkmers, visDFkmersTMP)
  
  # Add T
  visDFkmersTMP <- data.frame(kmer = rep("T", length(kmerRegionsDF$T)), score = kmerRegionsDF$T)
  visDFkmers <- rbind(visDFkmers, visDFkmersTMP)
  
  # Add C
  visDFkmersTMP <- data.frame(kmer = rep("C", length(kmerRegionsDF$C)), score = kmerRegionsDF$C)
  visDFkmers <- rbind(visDFkmers, visDFkmersTMP)
  
  # Add G
  visDFkmersTMP <- data.frame(kmer = rep("G", length(kmerRegionsDF$G)), score = kmerRegionsDF$G)
  visDFkmers <- rbind(visDFkmers, visDFkmersTMP)
  
  # Add A or T
  visDFkmersTMP <- data.frame(kmer = rep("A/T", length(kmerRegionsDF$A + kmerRegionsDF$T)), score = (kmerRegionsDF$A + kmerRegionsDF$T))
  visDFkmers <- rbind(visDFkmers, visDFkmersTMP)
  
  # Add C or G
  visDFkmersTMP <- data.frame(kmer = rep("C/G", length(kmerRegionsDF$G + kmerRegionsDF$C)), score = (kmerRegionsDF$G + kmerRegionsDF$C))
  visDFkmers <- rbind(visDFkmers, visDFkmersTMP)
  
  outputList <- list()
  outputList$kmersDF <- kmersDF
  outputList$visDFkmers <- visDFkmers
  outputList$kmerRegionsDF <- kmerRegionsDF
  
  return(outputList)
}

makeViskmerDF <- function(kmersDF, inRanges, inType) {
  kmersDF.Visualization <- kmersDF['A'] + kmersDF['T']
  colnames(kmersDF.Visualization) <- "A/T"
  
  kmersDF.Visualization["C/G"] <- kmersDF['C'] + kmersDF['G']
  kmersDF.Visualization[c('A', "T", "C", "G")] <- kmersDF[c('A', "T", "C", "G")]
  
  corcoefs.Visualization <- apply(kmersDF.Visualization, MARGIN = 2, FUN = cor, y = inRanges$reads_FPKM)
  corcoefs.Visualization.df <- data.frame(
    kmer = names(corcoefs.Visualization),
    corcoef = round(corcoefs.Visualization, digits = abs(ceiling(abs(log10(abs(corcoefs.Visualization)))))),
    type = inType
  )
  
  return(corcoefs.Visualization.df)
}


filterBam_TK <- function(
  ## INPUTS
  BAMFILE = NULL,
  BSSPECIES = NULL,
  EXTENTION = ".Aligned.sortedByCoord.out.bam",
  ## OPTIONS
  SIMPLECIGAR = TRUE,
  INCLUDE.SECONDARY.ALIGNEMNT = FALSE, # FALSE: only primary alignments are imported, TRUE: only secondary alignments are imported, NA: all alignments are imported
  GET.ORIGINAL.SEQUENCE = FALSE,
  STANDARD.CONTIGS.ONLY = TRUE,
  PERFECT.MATCH.ONLY = TRUE,
  FILTER.BY.FLAG = TRUE,
  SELECTFLAG = c(0, 16),
  USE.SIZE.FILTER = TRUE,
  READ.SIZE.RANGE = c(18, 50),
  TAGS = c("NH", "NM", "MD"),
  WHAT = c("flag"),
  ## name specific
  SPLIT.NAME.BY = "-"
) {
  ## AUTHOR: Pavol Genzor
  ## Use: Load BAM file into R
  ## 02.08.22; Version 8; chromosome filter, including multi-mappers
  ## 06.25.21; Version 7; refined
  ## NOTE ON BAM INFO
  ## Flag: 256 = not primary alignment; 272 = reverse strand not primary alignment;
  ## Flag: 0 = forward unpaired unique alignment; 16 = reverse unpaired unique alignment;
  ## Tags: NH:i:1 = unique alignment
  ## Tags: NM = edit distance to the reference
  ##
  ## NOTE: This program requires prepareFastq.R to be able to extract multiplicity information
  ## Libraries
  suppressPackageStartupMessages({
    library("data.table")
    library("dplyr")
    library("Rsamtools")
    library("GenomicAlignments")
    library("BSgenome.Hsapiens.UCSC.hg38")
    library("BSgenome.Dmelanogaster.UCSC.dm6")
    library("BSgenome.Mmusculus.UCSC.mm10")
  })

  ## Check input
  if (is.null(BAMFILE)) stop("Please provide full path to a .bam file !!!")
  if (is.null(BSSPECIES)) stop("Please provide BSSPECIES name !!!")
  if (isTRUE(GET.ORIGINAL.SEQUENCE)) { WHAT <- c("flag", "seq") }

  ## For report
  PROGRESS.L <- list()
  FILE.NAME <- gsub(EXTENTION, "", basename(BAMFILE))

  message("Processing ...")
  message("***\n")
  message("NOTE: To load all available reads (unique and multimapping first and other positions), set the following tag:")
  message("\tINCLUDE.SECONDARY.ALIGNEMNT = NA")
  message("***\n")

  ## PARAMETERS FOR LOADING BAM FILE
  PARAM <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = INCLUDE.SECONDARY.ALIGNEMNT),
    tag = TAGS, 
    simpleCigar = SIMPLECIGAR, 
    what = WHAT
  )
  
  message(" prepared loading parameters")
  message(paste0("\tTAGS:\t", paste(TAGS, collapse = ", ")))
  message(paste0("\tCIGAR:\t", ifelse(isTRUE(SIMPLECIGAR), "simple cigar", "all cigar")))
  message(paste0("\tWHAT:\t", paste0(WHAT, collapse = ", ")))
  message("")
  
  message(" loading .bam file into GAlignments")
  GA <- GenomicAlignments::readGAlignments(file = BAMFILE, use.names = TRUE, param = PARAM)

  ## ***
  GA.IN <- length(GA)
  PROGRESS.L[["INPUT"]] <- GA.IN
  message(paste0("\tIMPORTED: ", GA.IN))

  if (isTRUE(USE.SIZE.FILTER)) {
    message(" filtering by read size")
    message(paste0("\tRANGE:\t", paste0(READ.SIZE.RANGE, collapse = "-")))
    GA <- GA[width(GA) %in% seq(READ.SIZE.RANGE[1], READ.SIZE.RANGE[2], by = 1)]
    ## ***
    REMAINDER <- (length(GA) / GA.IN) * 100
    PROGRESS.L[["SIZE_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2)))
  }

  if (isTRUE(STANDARD.CONTIGS.ONLY)) {
    message(" removing non-standard contigs")
    GAR <- keepStandardChromosomes(x = GA, pruning.mode = "coarse")
    ## ***
    REMAINDER <- (length(GAR) / GA.IN) * 100
    PROGRESS.L[["CONTIG_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2)))
  } else { 
    GAR <- GA 
  }

  if (isTRUE(FILTER.BY.FLAG)) {
    if (is.na(INCLUDE.SECONDARY.ALIGNEMNT)) {
      SELECTFLAG <- c(0, 16, 256, 272)
      message(" keeping all alignments")
      GARP <- GAR
    } else {
      message(" selecting ONLY primary alignments")
      GARP <- GAR[mcols(GAR)[["flag"]] %in% SELECTFLAG]
    }
    ## ***
    REMAINDER <- (length(GARP) / GA.IN) * 100
    PROGRESS.L[["FLAG_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2)))
  } else { 
    GARP <- GAR 
  }

  if (isTRUE(PERFECT.MATCH.ONLY)) {
    message(" removing reads with mismatches")
    GARP <- GARP[mcols(GARP)[["NM"]] %in% c(0)]
    mcols(GARP)[["NM"]] <- NULL
    mcols(GARP)[["MD"]] <- NULL
    ## ***
    REMAINDER <- (length(GARP) / GA.IN) * 100
    PROGRESS.L[["MISMATCH_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2)))
  }

  if (isTRUE(GET.ORIGINAL.SEQUENCE)) {
    message(" retrieving original read sequences")
    BAMSEQ <- mcols(GARP)[["seq"]]
    ISONMINUS <- as.logical(GenomicAlignments::strand(GARP) == "-")
    BAMSEQ[ISONMINUS] <- Biostrings::reverseComplement(BAMSEQ[ISONMINUS])
    mcols(GARP)[["seq"]] <- BAMSEQ
  }

  message(" converting to GRanges")
  GARP.GR <- GenomicRanges::granges(GARP, use.names = TRUE, use.mcols = TRUE)
  PROGRESS.L[["FINAL"]] <- length(GARP.GR)

  message(" adding multiplicity column [MULT]")
  mcols(GARP.GR)[["MULT"]] <- as.integer(dplyr::nth(
    data.table::tstrsplit(dplyr::nth(data.table::tstrsplit(names(GARP.GR), split = SPLIT.NAME.BY), -1), split = "M"),
    -1
  ))

  if (length(GARP.GR[is.na(GARP.GR$MULT), ]) > 0) {
    GARP.GR[is.na(GARP.GR$MULT), ]$MULT <- 1
  }

  ## RETURN
  message("Done!")
  message("")
  print(as.data.table(PROGRESS.L))
  return(GARP.GR)
}



Bidirectional_clusters_prediction_TK <- function(
  clusters.GR = NULL,
  piRNAs.GR = NULL,
  Five_prime_percentage_overlap = 40,
  Three_prime_percentage_for_piRNAs_strand = 70,
  sense_antisense_ratio_threshold = 100
) {
  # Usage example:
  # Bidirectional_clusters_prediction_TK(clusters.GR = clusters.GR, piRNAs.GR = piRNAs.GR)
  
  # Libraries
  library(GenomicRanges)
  library(data.table)
  
  # Input checking
  if (is.null(clusters.GR))
    stop("please provide clusters coordinates in GRanges format")
  if (is.null(piRNAs.GR))
    stop("please provide mapped piRNAs in GRanges format")
  if (!is.numeric(Five_prime_percentage_overlap))
    stop("please provide valid percentage of cluster length that you want to use for investigating the diverging 5' overlaps - it has to be numeric eg: if you want to focus on the 5' most 20% of each cluster you need to type 20")
  if (!is.numeric(Three_prime_percentage_for_piRNAs_strand))
    stop("please provide valid percentage of cluster length that you want to use for checking the strand ratio of piRNAs in the 3'parts - it has to be numeric eg: if you want to focus on the 3' most 70% of each cluster you need to type 70, if you want to check piRNAs strand ratio throughout the cluster just type 100")
  if (!is.numeric(sense_antisense_ratio_threshold))
    stop("please provide numeric value for the threshold of the sense over antisense piRNAs ratio for bi-directional clusters to be made from 2 unistrand clusters and not include any dual-strand. Default is 100")
  if (!"MULT" %in% colnames(mcols(piRNAs.GR))) {
    piRNAs.GR$MULT <- 1
  } else if (length(piRNAs.GR[is.na(piRNAs.GR$MULT)]) > 0) {
    piRNAs.GR[is.na(piRNAs.GR$MULT)]$MULT <- 1
  }
  
  threshold <- sense_antisense_ratio_threshold
  Ovrlp_percent <- Five_prime_percentage_overlap / 100
  Three_prime <- Three_prime_percentage_for_piRNAs_strand / 100
  clusters.GR <- clusters.GR[order(clusters.GR)]
  clusters.GR$bidirectional <- FALSE
  
  # Keep the 5' most 40% of each cluster (to later check if they overlap with a beginning of a cluster on the opposite strand)
  clusters.GR_5ends <- GenomicRanges::resize(x = clusters.GR, width = width(clusters.GR) * Ovrlp_percent, fix = "start")
  clusters.GR_5endsFLIP <- GenomicRanges::invertStrand(clusters.GR_5ends)
  DT3 <- as.data.table(GenomicRanges::findOverlaps(query = clusters.GR_5ends, subject = clusters.GR_5endsFLIP, type = "any", ignore.strand = FALSE))
  clusters.GR[DT3$queryHits, ]$bidirectional <- TRUE
  
  # Remove the bidirectional flag from clusters that fall completely within another cluster in the opposite orientation
  clusters.GR_FLIP <- GenomicRanges::invertStrand(clusters.GR)
  DT4 <- as.data.table(GenomicRanges::findOverlaps(query = clusters.GR, subject = clusters.GR_FLIP, type = "within", ignore.strand = FALSE))
  if (length(clusters.GR[DT4$queryHits, ]$bidirectional) > 0) {
    clusters.GR[DT4$queryHits, ]$bidirectional <- FALSE
  }
  
  ### Incorporate piRNAs to distinguish bidirectional from dual-strand clusters ###
  ## get ratio of sense vs antisense piRNAs per cluster
  clusters2.GR <- clusters.GR
  clusters2.GR$resize_by <- width(clusters2.GR)
  clusters2.GR[clusters2.GR$bidirectional %in% TRUE]$resize_by <- Three_prime * width(clusters2.GR[clusters2.GR$bidirectional %in% TRUE])
  clusters2 <- GenomicRanges::resize(x = clusters2.GR, width = clusters2.GR$resize_by, fix = "end")
  
  # Use uniquely mapping piRNAs to figure out directionality
  piRNAs <- piRNAs.GR[piRNAs.GR$NH == 1]
  
  # Get plus-strand piRNAs and minus-strand piRNAs within each cluster
  piRNAs_plus <- piRNAs[strand(piRNAs) %in% "+"]
  piRNAs_minus <- piRNAs[strand(piRNAs) %in% "-"]
  
  # Find overlaps of clusters with + piRNAs
  DT_plus <- as.data.table(findOverlaps(clusters2, piRNAs_plus, ignore.strand = TRUE))
  # Retrieve multiplicity info
  DT_plus$MULT <- piRNAs_plus[DT_plus$subjectHits, ]$MULT
  DT_plus.2 <- as.data.table(DT_plus[, sum(MULT), by = "queryHits"])
  # Add numbers of plus-strand piRNAs to clusters
  clusters2$plus_piRNAs <- 0
  clusters2[DT_plus.2$queryHits, ]$plus_piRNAs <- DT_plus.2$V1
  
  ## Repeat for minus-strand piRNAs
  # Find overlaps of clusters with - piRNAs
  DT_minus <- as.data.table(findOverlaps(clusters2, piRNAs_minus, ignore.strand = TRUE))
  # Retrieve multiplicity info
  DT_minus$MULT <- piRNAs_minus[DT_minus$subjectHits, ]$MULT
  DT_minus.2 <- as.data.table(DT_minus[, sum(MULT), by = "queryHits"])
  # Add numbers of minus-strand piRNAs to clusters
  clusters2$minus_piRNAs <- 0
  clusters2[DT_minus.2$queryHits, ]$minus_piRNAs <- DT_minus.2$V1
  
  # Add ratio of sense/antisense piRNAs per cluster
  clusters3 <- unlist(GRangesList(lapply(as.list(seq(1:length(clusters2))), function(i) {
    GR <- clusters2[i, ]
    GR$sense_antisense_ratio <- 1
    if (as.character(strand(GR)) == "+") {
      GR$sense_antisense_ratio <- (GR$plus_piRNAs / GR$minus_piRNAs)
    } else {
      GR$sense_antisense_ratio <- (GR$minus_piRNAs / GR$plus_piRNAs)
    }
    return(GR)
  })))
  
  # Add sense_antisense_ratio of piRNAs to the initial GRanges file (which contain unaltered ranges - no resizing)
  clusters.GR$sense_antisense_ratio <- NA
  clusters.GR[order(clusters.GR)]$sense_antisense_ratio <- clusters3[order(clusters3)]$sense_antisense_ratio
  
  test <- clusters.GR
  test$partner_sense_antisense_ratio <- NA
  testb <- test[test$bidirectional %in% TRUE]
  test2 <- GenomicRanges::resize(x = testb, width = (testb$width_in_nt) * Ovrlp_percent, fix = "start")
  test2_FLIP <- GenomicRanges::invertStrand(test2)
  DT3 <- as.data.table(GenomicRanges::findOverlaps(query = test2, subject = test2_FLIP, type = "any", ignore.strand = FALSE))
  test[test$bidirectional %in% TRUE][DT3$queryHits, ]$partner_sense_antisense_ratio <- test2_FLIP[DT3$subjectHits, ]$sense_antisense_ratio
  
  # Add partner_sense_antisense_ratio to the initial GRanges file (which contain unaltered ranges - no resizing)
  clusters.GR$partner_sense_antisense_ratio <- NA
  clusters.GR[order(clusters.GR)]$partner_sense_antisense_ratio <- test[order(test)]$partner_sense_antisense_ratio
  
  clusters.GR_w_strandedness <- clusters.GR
  clusters.GR_w_strandedness$veryfied_bidirectional <- FALSE
  clusters.GR_w_strandedness[!is.na(clusters.GR_w_strandedness$partner_sense_antisense_ratio) &
                               clusters.GR_w_strandedness$bidirectional %in% TRUE &
                               clusters.GR_w_strandedness$sense_antisense_ratio > threshold &
                               clusters.GR_w_strandedness$partner_sense_antisense_ratio > threshold]$veryfied_bidirectional <- TRUE
  
  return(clusters.GR_w_strandedness)
}


lengthPerCluster <- function(
  bam_file, 
  chosenReplicate, 
  sampleName, 
  lengthDistrColor = NA
) {
  # Initialize dataframe to store entries
  column_names <- paste0(seq(20, 40), "M")
  lengths_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(lengths_df) <- column_names
  
  # Get reads per cluster
  reads <- scanBam(bam_file, param = ScanBamParam(
    which = chosenReplicate,
    what = c("cigar"),
    flag = scanBamFlag(isSecondaryAlignment = FALSE)
  ))
  
  for (cl_num in 1:length(reads)) {
    # Check if alignment gap, deletion, and insertion is a large fraction (more than 2%)
    totalAlgn <- length(reads[[cl_num]]$cigar)
    subsetAlgn <- reads[[cl_num]]$cigar[!grepl("N|D|I", reads[[cl_num]]$cigar, ignore.case = TRUE)]
    
    if (length(subsetAlgn) / totalAlgn < 0.98) {
      message(paste0(
        names(reads)[cl_num], 
        " has more than 2% alignments with alignment gap (N), Deletion (D), Insertion (I).\nExact percentage is: ", 
        length(subsetAlgn) / totalAlgn
      ))
    }
    
    # Create table with frequency of each length and rename columns
    lengthFrequency <- as.data.frame(table(subsetAlgn))
    names(lengthFrequency) <- c("Length", "Count")
    lengthFrequency$Count <- lengthFrequency$Count / totalAlgn
    
    # Initialize all columns with zero
    lengths_df[cl_num, ] <- 0
    
    # Go through each M length in the frequency table
    for (i in 1:nrow(lengthFrequency)) {
      col_name <- as.character(lengthFrequency$Length[i])
      if (col_name %in% names(lengths_df)) {
        lengths_df[cl_num, col_name] <- lengthFrequency$Count[i]
      }
    }
    
    # Include information about cluster, rank, total read counts, and used read counts (those without N/D/I)
    lengths_df$cluster[cl_num] <- names(reads[cl_num])
    lengths_df$strand[cl_num] <- as.character(strand(chosenReplicate[cl_num]))
    lengths_df$rank[cl_num] <- chosenReplicate$rank_readsExplained[cl_num]
    lengths_df$readCount_total[cl_num] <- totalAlgn
    lengths_df$readCount_Monly[cl_num] <- length(subsetAlgn)
  }
  
  # Column: peak of piRNA lengths
  lengths_df$globalMax <- apply(
    lengths_df[, !names(lengths_df) %in% c('cluster', 'strand', 'rank', 'readCount_total', 'readCount_Monly')], 
    1, 
    function(row) {
      max_col <- which.max(row)
      names(lengths_df)[max_col]
    }
  )
  
  # Column: sum of all considered reads without deletions/insertions/alignment gaps
  lengths_df$percentageTotal <- rowSums(lengths_df[, !names(lengths_df) %in% c('cluster', 'strand', 'rank', 'readCount_total', 'readCount_Monly', 'globalMax')])
  
  # Sort by rank and save as csv
  lengths_df <- lengths_df[order(lengths_df$rank), ]
  
  # Start creating Global Max freq table and stacked bar graph
  row.names(lengths_df) <- NULL
  
  # Frequency table of global maxima
  freq_GlobalMax <- as.data.frame(table(lengths_df$globalMax))
  
  # Make stacked bar graph with ratio of clusters
  freq_GlobalMax$ratio <- freq_GlobalMax$Freq / sum(freq_GlobalMax$Freq)
  freq_GlobalMax$group <- 'GlobalMax'
  
  # Set this order in the dataframe
  final_order <- sort(unique(freq_GlobalMax$Var1), decreasing = TRUE)
  freq_GlobalMax$Var1 <- factor(freq_GlobalMax$Var1, levels = final_order)
  
  # Filter by only read match lengths between 20 and 35
  col20to35 <- grep("^(2[0-9]M|30M|31M|32M|33M|34M|35M)$", names(lengths_df), value = TRUE)
  colWanted <- gsub("M", "", col20to35)
  lengthPerCluster_filtered <- lengths_df[, col20to35]
  
  # Calculate mean, standard deviation, and standard error for each piRNA length
  stats_df <- data.frame(column = colWanted)
  stats_df$mean <- sapply(lengthPerCluster_filtered, mean, na.rm = TRUE)
  stats_df$sd <- sapply(lengthPerCluster_filtered, sd, na.rm = TRUE)
  stats_df$se <- stats_df$sd / sqrt(nrow(lengthPerCluster_filtered))
  
  plotLengthDistr <- ggplot(stats_df, aes(x = column, y = mean)) +
    geom_bar(stat = "identity", position = position_dodge(), fill = lengthDistrColor, color = "black") +
    labs(x = "piRNA length (nt)", y = "Percentage of piC", title = paste0("Length distr/cluster — ", gsub("_", " ", sampleName))) +
    theme_pubr(base_size = 12, base_family = "Helvetica")
  
  violinPlotLengthDistr <- NULL
  
  # Only create error bars and violin plot if more than 2 clusters present
  if (nrow(lengths_df) > 2) {
    plotLengthDistr <- plotLengthDistr + 
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = "black", alpha = 0.9)
    
    # Create box for under length distribution
    lengthFreq_df <- as.data.frame(select(as.data.frame(freq_GlobalMax), Var1, Freq))
    
    # Remove 'M' and convert Var1 to integer
    lengthFreq_df$Var1 <- as.integer(sub("M", "", freq_GlobalMax$Var1))
    long_data <- rep(lengthFreq_df$Var1, times = lengthFreq_df$Freq)
    
    # Creating a density violin plot
    violinPlotLengthDistr <- ggplot(data = data.frame(Var1 = long_data), aes(x = factor(1), y = Var1)) +
      geom_violin(trim = TRUE, adjust = 3, fill = lengthDistrColor, alpha = 0.3) +
      coord_flip() +
      scale_y_continuous(limits = c(20, 35)) +
      labs(title = "", x = "Sign. length", y = "piRNA length (nt)") +
      theme_pubr(base_size = 12, base_family = "Helvetica")
  }
  
  barAndViolinLengthDistr <- plotLengthDistr / violinPlotLengthDistr + plot_layout(heights = c(7, 0.8))
  
  return(list(
    df_allInfo = lengths_df, 
    freq_GlobalMax = t(freq_GlobalMax),
    barLengthDistr = barAndViolinLengthDistr
  ))
}
