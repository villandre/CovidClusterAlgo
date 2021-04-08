#' Clusters SARS-Cov-2 sequencing data
#'
#' The function uses MCMC to produce cluster estimates based on sequencing data and covariate information. The function automatically derives suitable starting values for all parameters.
#'
#' @param DNAbinData DNAbin object, the sequencing data,
#' @param rootSequenceName Name of the sequence used to root the tree, has to match a name in DNAbinData
#' @param clusterRegion optional string indicating region for which clusters should be computed; if left NULL, only chain results are returned
#' @param covariateFrame (ignored for now) data.frame containing covariate values in the order of elements in DNAbinData, used for scoring sample partitions
#' @param seqsTimestampsPOSIXct named vector of sequence collection times, with names matching sequence names in DNAbinData
#' @param seqsRegionStamps named vector of region labels for sequences, with names matching sequence names in DNAbinData
#' @param perSiteClockRate fixed mutation rate in nucleotide substitution per site per year
#' @param startingValuePhylo (optional) phylo object, starting value for the phylogeny
#' @param evoParsList (optional) list of values for the evolutionary parameters used in the likelihood function specified in control$logLikFun, see for example ?ape::pml
#' @param clusterScoringFun (ignored for now) function with four arguments (in order): cluster labels (numeric), phylogeny (phylo) genotyping data (DNAbin), covariate information (data.frame).
#' @param control List of control parameters, cf. findBayesianClusters.control
#'
#' @details By default, the function uses `phangorn::pml` to compute log-likelihoods, and `phangorn::optim.pml` to obtain a starting value for the phylogeny. The function then fixes the evolutionary parameters at their ML value. The user can input their own phylogenetic likelihood function with `control$logLikFun`. In that case however, values will also have to be provided for `startingValuesPhylo` and `evoParsList`.
#'
#' We strongly recommend setting `control$MCMC.control$folderToSaveIntermediateResults`. Running a reasonably long chain on a sample of a moderate size, e.g. 2,000, will take at least a number of days: that option allows the user to interrupt and resume a run. Each chain is assigned a glyph: its value is specified in the function's printout, e.g.
#' "Chain is called Ph9t7z. Specify this string in MCMC.control if you want to resume simulations."
#'
#' You can find cluster estimates from any tree in `chain` and for any region by using the `getRegionClusters` function.
#'
#' @return A list with two elements
#' \itemize{
#' \item{`chain`}{list giving the thinned chain results,}
#' \item{`MAPclusters`}{vector giving the cluster membership indices for sequences in region clusterRegion based on the maximum posterior probability tree; a 0 indicates that a sequence is not in that region.}
#' }
#'
#' @examples
#' \dontrun{
#' # See example in vignette.
#' }
#' @export

covidCluster <- function(
  DNAbinData,
  MLtree = FALSE,
  outgroup,
  clusterRegion = NULL,
  rootRegion = NULL,
  covariateFrame = NULL,
  seqsTimestampsPOSIXct = NULL,
  folderForTreeEstFiles = NULL,
  seqsRegionStamps = rep(1, ifelse(is.matrix(DNAbinData), nrow(DNAbinData), length(DNAbinData))),
  perSiteClockRate = NULL,
  clusteringCriterion = c("mrca", "cophenetic", "consecutive"),
  distLimit = NULL,
  raxmlExecWithPath = "/home/luc/standard-RAxML/raxmlHPC-PTHREADS-AVX2",
  control = list(),
  MrBayes.control = gen.MrBayes.control(),
  RAxML.control = gen.RAxML.control(),
  priors.control = gen.priors.control()) {
  clusteringCriterion <- clusteringCriterion[[1]]
  outgroupPos <- match(outgroup, rownames(DNAbinData))
  rownames(DNAbinData) <- .editTaxonNames(rownames(DNAbinData))
  names(seqsTimestampsPOSIXct) <- .editTaxonNames(names(seqsTimestampsPOSIXct))
  names(seqsRegionStamps) <- .editTaxonNames(names(seqsRegionStamps))
  outgroup <- rownames(DNAbinData)[[outgroupPos]] # We fix the outgroup name since it may contain characters the MrBayes does not understand.
  perSiteClockRate <- perSiteClockRate/365 # Time is expressed in days in the code, whereas perSiteClockRate is expressed in substitutions per site per *year*.

  control <- do.call("gen.covidCluster.control", control)
  control$clusteringCriterion <- clusteringCriterion
  control$distLimit <- distLimit
  control$clusterRegion <- clusterRegion
  MrBayes.control <- do.call("gen.MrBayes.control", MrBayes.control)
  RAxML.control <- do.call("gen.RAxML.control", RAxML.control)
  priors.control <- do.call("gen.priors.control", priors.control)
  treeSampleList <- NULL
  if (!MLtree) {
    # We call MrBayes
    nexusFilename <- "mrbayesData.nex"
    nexusFilenameWithFolder <- paste(folderForTreeEstFiles, nexusFilename, sep = "/")
    if (!control$skipTreeEstimation) {
      # We write necessary files to disk
      scriptFilenameWithFolder <- .writeMrBayesFiles(DNAbinData, nexusFilename, folderForTreeEstFiles, outgroup, MrBayes.control)

      commandForSystem2 <- MrBayes.control$MrBayesShellCommand
      argsForSystem2 <- scriptFilenameWithFolder
      if (MrBayes.control$MPIenabled) {
        commandForSystem2 <- "mpirun"
        argsForSystem2 <- c(paste("-n", control$numThreads), paste(MrBayes.control$MrBayesShellCommand, scriptFilenameWithFolder))
      }
      system2(command = commandForSystem2, args = argsForSystem2)
    }
    MrBayesOutputFilenamePrefix <- paste(nexusFilenameWithFolder, ".run1", sep = "")
    # We read in the tree samples and parameter files
    treeSampleFile <- paste(MrBayesOutputFilenamePrefix, ".t", sep = "")
    treeSample <- ape::read.nexus(treeSampleFile)
    numItersToDrop <- ceiling(MrBayes.control$burninfrac * length(treeSample))
    itersToKeep <- seq(from = numItersToDrop + 1, to = length(treeSample), by = ceiling(1/control$MrBayesOutputThinningRate))
    treeSample <- treeSample[itersToKeep]
    treeSampleList <- lapply(seq_along(treeSample), function(index) treeSample[[index]])
    rm(treeSample)
    gc() # Tree sample can be huge. Better to call gc explicitly
    # We also read in associated parameter values...
  } else {
    if (!control$skipTreeEstimation) {
    treeSampleList <- list(obtainMLtreeRAxML(alignment = DNAbinData, raxmlCommand = raxmlExecWithPath, directory = folderForTreeEstFiles, outgroupName = outgroup, control = RAxML.control))
    } else {
      treeSampleFile <- list.files(path = folderForTreeEstFiles, pattern = "result.optimizedBranches", full.names = TRUE)
      treeSampleList <- list(ape::read.tree(treeSampleFile))
    }
  }

  # parameterValuesFile <- paste(MrBayesOutputFilenamePrefix, ".p", sep = "")
  # parameterValues <- .formatParameterFiles(parameterValuesFile, itersToKeep = itersToKeep, control = MrBayes.control)

  clusMembershipVecList <- .simulateClusMembership(phyloList = treeSampleList, targetRegion = clusterRegion, rootRegion = rootRegion, timestamps = seqsTimestampsPOSIXct, regionStamps = seqsRegionStamps, clockRate = perSiteClockRate, covidCluster.control = control, priors.control = priors.control)
  output <- produceClusters(clusMembershipVecList, control = control)
  if (control$saveArgs) {
    output$args <- mget(names(formals()))
  }
  output$introSizeFreqTables <- clusMembershipVecList$introSizeFreqTables
  output$clusSizeDist <- clusMembershipVecList$clusSizeDist
  output$numTrees <- clusMembershipVecList$numTrees
  output
}

gen.covidCluster.control <- function(lengthForNullExtBranchesInPhylo = 1e-8, numReplicatesForNodeTimes = 500, numReplicatesForCoalRates = 500,  numThreads = 1, clusterCriterion = c("mrca", "cophenetic"), skipTreeEstimation = FALSE, MrBayesOutputThinningRate = 1, hclustMethod = "single", linkageRequirement = 0.90, saveArgs = FALSE, logLambdaPriorSD = log(2)/2, RNGseed = 10L) {
  list(lengthForNullExtBranchesInPhylo = lengthForNullExtBranchesInPhylo, numReplicatesForNodeTimes = numReplicatesForNodeTimes,  numReplicatesForCoalRates = numReplicatesForCoalRates, numThreads = numThreads, clusterCriterion = clusterCriterion[[1]], skipTreeEstimation = skipTreeEstimation, MrBayesOutputThinningRate = MrBayesOutputThinningRate, hclustMethod = hclustMethod, linkageRequirement = linkageRequirement, saveArgs = saveArgs, logLambdaPriorSD = logLambdaPriorSD, RNGseed = RNGseed)
}

clusterFromMrBayesOutput <- function(seqsTimestampsPOSIXct, seqsRegionStamps, MrBayesTreesFilename, MrBayesParametersFilename = NULL, clusterRegion, rootRegion, clusterCriterion, burninFraction = 0.5, linkageRequirement = 0.5, distLimit, perSiteClockRate, control = gen.covidCluster.control()) {
  perSiteClockRate <- perSiteClockRate/365 # Time is expressed in days in the code, whereas perSiteClockRate is expressed in substitutions per site per *year*.

  control <- do.call("gen.covidCluster.control", control)
  control$linkageRequirement <- linkageRequirement
  control$clusteringCriterion <- clusterCriterion
  control$distLimit <- distLimit
  control$clusterRegion <- clusterRegion

  # if (is.null(MrBayesParametersFilename)) {
  #   fileRoot <- stringr::str_extract(MrBayesTreesFilename, pattern = ".*(?=\\.t$)")
  #   MrBayesParametersFilename <- paste(fileRoot, ".p", sep = "")
  #   if (!file.exists(MrBayesParametersFilename)) stop("Could not infer parameter values filename. Please enter value for argument MrBayesParametersFilename.")
  # }
  treeSample <- ape::read.nexus(MrBayesTreesFilename)
  numItersToDrop <- ceiling(burninFraction * length(treeSample))
  itersToKeep <- seq(from = numItersToDrop + 1, to = length(treeSample), by = ceiling(1/control$MrBayesOutputThinningRate))
  treeSample <- treeSample[itersToKeep]
  treeSampleList <- lapply(seq_along(treeSample), function(index) treeSample[[index]])
  rm(treeSample)
  gc() # Tree sample can be huge. Better to call gc explicitly

  # parameterValues <- .formatParameterFiles(MrBayesParametersFilename, itersToKeep = itersToKeep)

  clusMembershipVecList <- .simulateClusMembership(phyloList = treeSampleList, targetRegion = clusterRegion, rootRegion = rootRegion, timestamps = seqsTimestampsPOSIXct, regionStamps = seqsRegionStamps, clockRate = perSiteClockRate, covidCluster.control = control)
  output <- produceClusters(clusMembershipVecList, control = control)
  if (control$saveArgs) {
    output$args <- mget(names(formals()))
  }
  output$introSizeFreqTables <- clusMembershipVecList$introSizeFreqTables
  output$clusSizeDist <- clusMembershipVecList$clusSizeDist
  output$numTrees <- clusMembershipVecList$numTrees
  output
}

.convertClusterListToVecOfIndices <- function(clusterList, seqNames) {
  clusterNums <- seq_along(clusterList)
  clusterVec <- rep(0, length(seqNames))
  names(clusterVec) <- seqNames
  for (clusNumber in clusterNums) {
    clusterVec[clusterList[[clusNumber]]] <- clusNumber
  }
  clusterVec
}

# The default run uses the following parameters:
# GTR model (nst = 6)
# Discrete gamma rate variation with a proportion of invariable sites (rates = "invgamma")
# Number of rate variation categories: 4 (ngammacat = 4)
# Two runs are executed simultaneously for diagnostic purposes (nruns = 2)
# Four heated chains are produced simultaneously, 1 cold, 3 hot (nchains = 4)
# Total number of iterations: 1e6 (ngen = 1000000)
# Thinning ratio: 1 out of 1000 iteration (samplefreq = 1000)
# Convergence diagnostics are computed once every 5000 iterations (diagnfreq = 5000)
# Chain characteristics are printed once every 500 iterations (printfreq = 500)
# The first 20% of iterations are discarded as a burn in (burninfrac = 0.2)

gen.MrBayes.control <- function(MrBayesShellCommand = "mb", nst = 6L, rates = "invgamma", ngammacat = 4L, nruns = 2L, nchains = 4L, ngen = 1000000L, samplefreq = 1000L, diagnfreq = 5000L, printfreq = 500L, burninfrac = 0.2, beaglesse = "yes", usebeagle = "no", beaglescaling = "dynamic", seed = 42L, swapseed = 1L, temp = 0.1, MPIenabled = FALSE, append = "no") {
  list(MrBayesShellCommand = MrBayesShellCommand, nst = nst, rates = rates, ngammacat = ngammacat, nruns = nruns, nchains = nchains, ngen = ngen, samplefreq = samplefreq, diagnfreq = diagnfreq, printfreq = printfreq, burninfrac = burninfrac, beaglesse = beaglesse, usebeagle = usebeagle, beaglescaling = beaglescaling, seed = seed, swapseed = swapseed, MPIenabled = MPIenabled, append = append, temp = temp)
}

gen.RAxML.control <- function(numRateCats = 4, numThreads = 4, seedForBS = 42, seedForParsemony = 24) {
  list(numRateCats = numRateCats, numThreads = numThreads, seedForBS = seedForBS, seedForParsemony = seedForParsemony)
}

gen.priors.control <- function() {
 list()
}

.simulateClusMembership <- function(phyloList, targetRegion, rootRegion, timestamps, regionStamps, clockRate, covidCluster.control, priors.control) {
  timestampsInDays <- as.numeric(timestamps)/86400

  names(timestampsInDays) <- names(timestamps)
  cl <- clusMembershipCondOnPhylo <- NULL

  if (covidCluster.control$numThreads > 1) {
    cat("Simulating transmission trees... ")
    cl <- parallel::makePSOCKcluster(covidCluster.control$numThreads)
    parallel::clusterSetRNGStream(cl = cl, iseed = covidCluster.control$RNGseed)
    clusMembershipCondOnPhylo <- parallel::parLapply(cl = cl, X = phyloList, fun = .simulateClustersFromStartingTree, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = log(clockRate), targetRegion = targetRegion, rootRegion = rootRegion, control = covidCluster.control)
    cat("Done \n")
    parallel::stopCluster(cl)
  } else {
    cat("Simulating transmission trees... ")
  # clusMembershipCondOnPhylo <- lapply(phyloList, FUN = function(listElement) tryCatch(expr = .simulateClustersFromStartingTree(phylogeny = listElement, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = log(clockRate), targetRegion = targetRegion, control = covidCluster.control), error = function(e) e))
    set.seed(covidCluster.control$RNGseed)
    clusMembershipCondOnPhylo <- lapply(phyloList, FUN = .simulateClustersFromStartingTree, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = log(clockRate), targetRegion = targetRegion, rootRegion = rootRegion, control = covidCluster.control)
  cat("Done \n")
  }
  indicesToKeep <- sapply(clusMembershipCondOnPhylo, FUN = function(listElement) !("error" %in% class(listElement)))

  combinedVecs <- do.call("c", lapply(clusMembershipCondOnPhylo[indicesToKeep], "[[", "clusVectors"))
  seqsInOrderNames <- names(combinedVecs[[1]])
  clusSizeDist <- summariseClusSizeDistsRcpp(combinedVecs)
  clusSizeDist <- clusSizeDist[order(as.numeric(names(clusSizeDist)))]
  list(
    clusMemVecList = lapply(combinedVecs, function(clusMemVec) as.integer(factor(clusMemVec[seqsInOrderNames]))),
    seqNames = seqsInOrderNames,
    introSizeFreqTables = lapply(clusMembershipCondOnPhylo, "[[", "introSizesDist"), # Call to factor removes sequence names... For memory reasons, better to keep names separate.
    clusSizeDist = clusSizeDist,
    numTrees = length(combinedVecs))
}

.writeMrBayesFiles <- function(DNAbinData, nexusFilename, folderForTreeEstFiles, outgroup, control) {
  nexusFilenameWithFolder <- paste(folderForTreeEstFiles, nexusFilename, sep = "/")
  scriptForMrBayes <- .produceMrBayesScript(outgroup = outgroup, nexusDataFilename = nexusFilenameWithFolder, control = control)

  ape::write.nexus.data(x = DNAbinData, file = nexusFilenameWithFolder)
  scriptFilenameWithFolder <- paste(folderForTreeEstFiles, "MrBayesScript.nex", sep = "/")
  fileConn <- file(scriptFilenameWithFolder)
  writeLines(scriptForMrBayes, fileConn)
  close(fileConn)
  scriptFilenameWithFolder
}

.produceMrBayesScript <- function(outgroup, nexusDataFilename, control = gen.MrBayes.control()) {
  paste("begin mrbayes;\n set autoclose=yes nowarn=yes seed=", as.integer(control$seed), " swapseed=", as.integer(control$swapseed), ";\n execute ", nexusDataFilename, ";\n lset nst=", as.integer(control$nst), " rates=", control$rates, ";\n outgroup ", outgroup, ";\n set usebeagle=", control$usebeagle, " beaglescaling=", control$beaglescaling," beaglesse=", control$beaglesse, ";\n mcmc nruns=", as.integer(control$nruns), " nchains=", as.integer(control$nchains), " ngen=", as.integer(control$ngen), " samplefreq=", as.integer(control$samplefreq), " diagnfreq=", as.integer(control$diagnfreq), " printfreq=", as.integer(control$printfreq), " append=", control$append, " temp=",  control$temp, ";\n sump relburnin=yes burninfrac=", control$burninfrac, ";\n end;", sep = "")
}

.formatParameterFiles <- function(filenames, itersToKeep, control = NULL) {
  parameterTables <- lapply(filenames, function(filename) {
    paraTable <- read.table(file = filename, header = TRUE, sep = "\t", skip = 1)
    paraTable[itersToKeep, ]
  })
  mergedTables <- do.call("rbind", parameterTables)
  list(logLik = mergedTables$LnL, logPrior = mergedTables$LnPr, treeLength = mergedTables$TL, evoPars = mergedTables[ , -(1:4)])
}

# subtreeClusterFun is a function that takes a phyloAndTransTree object and a region index and produces a vector of cluster membership indices.

.simulateClustersFromStartingTree <- function(phylogeny,  timestampsInDays, regionStamps, logClockRatePriorMean, targetRegion, rootRegion, control) {

  phyloAndTransTree <- .genStartPhyloAndTransTree(phylogeny = phylogeny, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = logClockRatePriorMean, targetRegion = targetRegion, rootRegion = rootRegion, control = control)

  funToReplicate <- function(phyloAndTransTree, control) {
    medianRates <- sapply(phyloAndTransTree$LambdaList, "[[", "Lambda")
    coalRates <- exp(rnorm(n = length(phyloAndTransTree$LambdaList), mean = log(medianRates), sd = control$logLambdaPriorSD))
    # coalRates <- rlnorm(n = length(medianRates), meanlog = log(medianRates), sdlog = abs(log(medianRates)))
    for (i in seq_along(phyloAndTransTree$LambdaList)) {
      phyloAndTransTree$LambdaList[[i]]$Lambda <- coalRates[[i]]
    }
    mode(phyloAndTransTree$edge) <- "integer"
    funToReplicateInner <- function(phyloAndTransTree, control) {
      numTips <- length(phyloAndTransTree$tip.label)
      orderedVertices <- rev(phyloAndTransTree$vertexOrderByDepth)
      tipTimes <- sapply(phyloAndTransTree$tip.label, "[[", "time")

      nodeTimesAndEdgeLengths <- simulateNodeTimesRcpp(
        numTips = as.integer(numTips),
        baseRatePerIntroduction = coalRates,
        orderedVertices = as.integer(orderedVertices),
        subtreeIndexVec = as.integer(phyloAndTransTree$subtreeIndexVec),
        tipTimes = tipTimes,
        edgeMatrix = phyloAndTransTree$edge,
        childrenNumList = phyloAndTransTree$childrenNumList,
        branchMatchIndexVec = phyloAndTransTree$branchMatchIndex,
        parentNumVec = phyloAndTransTree$parentNumVec)

      subtreeClusterFun <- function(phyloAndTransTree, subtreeIndex, distLimit, clusterRegion, clusteringCriterion, distTipsAncestorsMatrix) {
        numTips <- length(phyloAndTransTree$tip.label)
        tipsInSubtreeAndRegion <- phyloAndTransTree$tipsInRegionBySubtree[[subtreeIndex]]
        seqNames <- sapply(tipsInSubtreeAndRegion, function(tipNum) phyloAndTransTree$tip.label[[tipNum]]$name)
        output <- seq_along(seqNames)
        names(output) <- seqNames
        if (identical(phyloAndTransTree$node.label[[phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum - numTips]]$region, clusterRegion)) {
          clusterFunction <- getMRCAclustersRcpp
          if (control$clusterCriterion == "cophenetic") {
            clusterFunction <- getCopheneticClustersRcpp
          }

          output <- clusterFunction(
            parentNumVec = phyloAndTransTree$parentNumVec,
            childrenNumList = phyloAndTransTree$childrenNumList,
            descendedTipsList = phyloAndTransTree$descendedTips,
            subtreeIndexVec = phyloAndTransTree$subtreeIndexVec,
            vertexRegionVec = phyloAndTransTree$vertexRegionVec,
            tipNamesVec = phyloAndTransTree$tipNamesVec,
            subtreeRootNum = phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum,
            distTipsAncestorsMatrix = nodeTimesAndEdgeLengths$distTipsAncestorsMatrix,
            subtreeIndex = subtreeIndex,
            numTips = numTips,
            regionLabel = clusterRegion,
            distLimit = distLimit)
          output <- as.integer(output[seqNames])
          names(output) <- seqNames # as.integer removes the names
        }
        output
      }
      clustersBySubtree <- sapply(seq_along(phyloAndTransTree$LambdaList), subtreeClusterFun, phyloAndTransTree = phyloAndTransTree, distLimit = control$distLimit, clusteringCriterion = control$clusteringCriterion, clusterRegion = control$clusterRegion, distTipsAncestorsMatrix = nodeTimesAndEdgeLengths$distTipsAncestorsMatrix)
      shiftValue <- 0
      for (i in 2:length(clustersBySubtree)) {
        if (length(clustersBySubtree[[i - 1]]) > 0) {
          shiftValue <- max(clustersBySubtree[[i - 1]])
        }
        clustersBySubtree[[i]] <- clustersBySubtree[[i]] + shiftValue
      }
      do.call("c", clustersBySubtree)
    }
    replicate(n = control$numReplicatesForCoalRates, expr = funToReplicateInner(phyloAndTransTree = phyloAndTransTree, control = control), simplify = FALSE)
  }
  funToGetIntroSize <- function(subtreeIndex) {
    subtreeRootNum <- phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum
    introSize <- 0
    if (phyloAndTransTree$node.label[[subtreeRootNum - length(phyloAndTransTree$tip.label)]]$region == targetRegion) {
      introSize <- sum(sapply(phyloAndTransTree$tip.label, function(listElement) (listElement$region == targetRegion) & (listElement$subtreeIndex) == subtreeIndex))
    }
    introSize
  }
  introSizes <- sapply(seq_along(phyloAndTransTree$LambdaList), funToGetIntroSize)
 list(clusVectors = do.call("c", replicate(n = control$numReplicatesForNodeTimes, funToReplicate(phyloAndTransTree, control = control), simplify = FALSE)), introSizesDist = table(introSizes[introSizes != 0]))
}

.clusterIndicesFromAdjMat <- function(adjacencyMatrix) {
  clusterMembershipVector <- rep(0, nrow(adjacencyMatrix))
  clusterIndex <- 1
  repeat {
    lineNum <- match(0, clusterMembershipVector)
    updateIndices <- which(adjacencyMatrix[lineNum, ] != 0)
    clusterMembershipVector[updateIndices] <- clusterIndex
    clusterIndex <- clusterIndex + 1
    if (!(0 %in% clusterMembershipVector)) break
  }
  clusterMembershipVector
}

# The following function computes log(sum(vector)) when only the logarithms of each element are available.
# This is used to avoid computational zeros.

computeLogSum <- function(logValues) {
  maximumValue <- max(logValues)
  maximumValue + log(sum(exp(logValues - maximumValue)))
}

# MrBayes does not like "/" in taxon names. We will replace them with underscores.

.editTaxonNames <- function(taxonNames) {
  newTaxonNames <- gsub(x = taxonNames, pattern = "/", replacement = "_")
  newTaxonNames <- gsub(x = newTaxonNames, pattern = "\\|", replacement = "__")
  newTaxonNames <- gsub(x = newTaxonNames, pattern = " ", replacement = "")
  newTaxonNames
}

.genStartPhyloAndTransTree <- function(phylogeny, timestampsInDays, regionStamps, logClockRatePriorMean, targetRegion, rootRegion, control) {
  phyloAndTransTree <- phylogeny
  phyloAndTransTree$edge.length <- replace(phyloAndTransTree$edge.length, which(phyloAndTransTree$edge.length == 0), control$lengthForNullExtBranchesInPhylo)
  phyloAndTransTree$tip.label <- lapply(phyloAndTransTree$tip.label, function(x) list(name = x, time = timestampsInDays[[x]], region = regionStamps[[x]]))
  phyloAndTransTree$node.label <- lapply(1:ape::Nnode(phyloAndTransTree), function(node) {
    list(time = NA, region = NA)
  })
  phyloAndTransTree$edge.length <- lapply(seq_along(phyloAndTransTree$edge.length), function(edgeIndex) {
    list(phylogeny = phyloAndTransTree$edge.length[[edgeIndex]], transmissionTree = NA, logXi = NA)
  })
  logXi <- logClockRatePriorMean
  if (length(logXi) == 1) logXi <- rep(logClockRatePriorMean, length(phyloAndTransTree$edge.length))
  for (i in seq_along(phyloAndTransTree$edge.length)) {
    phyloAndTransTree$edge.length[[i]]$logXi <- logXi[[i]]
  }
  phyloAndTransTree <- .incrementPhyloStructOnly(phyloAndTransTree)
  phyloAndTransTree <- .identifyNodeRegions(phyloAndTransTree, rootRegion = rootRegion)
  # Initialises transmission tree branch lengths conditional on phylogenetic branch lengths and an equal clock rate across branches
  subtreeRootNodes <- .computeSubtreeRootNums(phyloAndTransTree)
  subtreeRootNodeDepths <- sapply(subtreeRootNodes, function(nodeNum) length(phangorn::Ancestors(phyloAndTransTree, nodeNum, "all")) + 1)
  updateOrder <- order(subtreeRootNodeDepths)
  subtreeRootNodes <- subtreeRootNodes[updateOrder]

  for (subtreeIndex in seq_along(subtreeRootNodes)) {
    nodeNum <- subtreeRootNodes[[subtreeIndex]]
    phyloAndTransTree$node.label[[nodeNum - ape::Ntip(phyloAndTransTree)]]$subtreeIndex <- subtreeIndex
    for (descNodeNum in phangorn::Descendants(phyloAndTransTree, nodeNum, "all")) {
      if (descNodeNum <= ape::Ntip(phyloAndTransTree)) {
        phyloAndTransTree$tip.label[[descNodeNum]]$subtreeIndex <- subtreeIndex
      } else {
        phyloAndTransTree$node.label[[descNodeNum - ape::Ntip(phyloAndTransTree)]]$subtreeIndex <- subtreeIndex
      }
    }
  }
  phyloAndTransTree <- .incrementPhyloRegionInfo(phyloAndTransTree, clusterRegion = targetRegion)
  phyloAndTransTree$LambdaList <- lapply(seq_along(subtreeRootNodes), function(i) list(Lambda = NULL, rootNodeNum = subtreeRootNodes[[i]]))
  # coalRates <- sapply(seq_along(phyloAndTransTree$LambdaList), .computeCoalRateSubtree, phyloAndTransTree = phyloAndTransTree)
  coalRates <- sapply(seq_along(phyloAndTransTree$LambdaList), .computeCoalRateSubtreeAlt, phyloAndTransTree = phyloAndTransTree)
  for (i in seq_along(coalRates)) {
    phyloAndTransTree$LambdaList[[i]]$Lambda <- coalRates[[i]]
  }

  phyloAndTransTree
}

.simulateNodeTimes <- function(phyloAndTransTree) {
  numTips <- length(phyloAndTransTree$tip.label)
  baseRatePerIntroduction <- sapply(phyloAndTransTree$LambdaList, "[[", "Lambda")
  # orderedVertices <- .getVertexOrderByDepth(phyloAndTransTree, reverse = TRUE) # Represents the topology
  orderedVertices <- rev(phyloAndTransTree$vertexOrderByDepth)
  orderedNodes <- orderedVertices[orderedVertices > numTips]
  for (nodeNum in orderedNodes) {
    childrenNums <- NULL
    childrenNums <- phyloAndTransTree$childrenNumList[[nodeNum]]

    introForMerge <- phyloAndTransTree$node.label[[nodeNum - length(phyloAndTransTree$tip.label)]]$subtreeIndex
    minChildrenTimes <- min(sapply(childrenNums, function(childNum) {
      returnTime <- NULL
      if (childNum <= numTips) {
        returnTime <- phyloAndTransTree$tip.label[[childNum]]$time
      } else {
        returnTime <- phyloAndTransTree$node.label[[childNum - numTips]]$time
      }
      returnTime
    }))
    phyloAndTransTree$node.label[[nodeNum - numTips]]$time <- minChildrenTimes - rexp(n = 1, rate = baseRatePerIntroduction[[introForMerge]])
  }
  phyloAndTransTree <- .addTransTreeEdgeLengths(phyloAndTransTree)
  phyloAndTransTree
}

.computeCoalRateSubtreeAlt <- function(phyloAndTransTree, subtreeIndex, control) {
  numTips <- length(phyloAndTransTree$tip.label)
  subtreeRootNum <- phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum
  nodeIndicesInSubtree <- which(sapply(phyloAndTransTree$node.label, "[[", "subtreeIndex") == subtreeIndex) + numTips

  vertexTimes <- numeric(length(phyloAndTransTree$edge.length) + 1)
  vertexTimes[[subtreeRootNum]] <- 1e-200 # Just to distinguish it from unspecified values.
  vertexInDepthOrder <- intersect(phyloAndTransTree$vertexOrderByDepth, nodeIndicesInSubtree)

  for (vertexNum in vertexInDepthOrder) {
    if (vertexNum <= numTips) next
    childrenNums <- phyloAndTransTree$childrenNumList[[vertexNum]]
    currentTime <- vertexTimes[[vertexNum]]
    childrenTimes <- currentTime + sapply(phyloAndTransTree$edge.length[phyloAndTransTree$branchMatchIndex[childrenNums]], "[[", "phylogeny")/exp(sapply(phyloAndTransTree$edge.length[phyloAndTransTree$branchMatchIndex[childrenNums]], "[[", "logXi"))
    vertexTimes[childrenNums] <- childrenTimes
  }

  timeRanges <- sapply(nodeIndicesInSubtree, function(nodeNum) {
    childrenNums <- phyloAndTransTree$childrenNumList[[nodeNum]]
    currentTime <- vertexTimes[[nodeNum]]
    childrenTimes <- vertexTimes[childrenNums]
    min(abs(childrenTimes - currentTime))
  })

  length(timeRanges)/sum(timeRanges)
}

.computeCoalRateSubtree <- function(phyloAndTransTree, subtreeIndex) {
  subtreeRootNum <- phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum
  subtreeRegion <- phyloAndTransTree$node.label[[subtreeRootNum - length(phyloAndTransTree$tip.label)]]$region
  timesToConsider <- 0
  nodeOrTip <- "node"
  verticesToReview <- subtreeRootNum
  index <- 1
  phyloAndTransTreeCopy <- phyloAndTransTree
  phyloAndTransTreeCopy$node.label[[subtreeRootNum - length(phyloAndTransTreeCopy$tip.label)]]$time <- 0
  repeat {
    vertexNum <- verticesToReview[[index]]
    currentTime <- phyloAndTransTreeCopy$node.label[[vertexNum - length(phyloAndTransTreeCopy$tip.label)]]$time
    childrenNums <- phyloAndTransTree$childrenNumList[[vertexNum]]
    matchingBranchNums <- phyloAndTransTree$branchMatchIndex[childrenNums]
    if (is.null(childrenNums)) { # We're initialising, so both childrenNums and matchingBranchNums will be NULL...
      childrenNums <- phyloAndTransTree$childrenNumList[[vertexNum]]
      matchingBranchNums <- phyloAndTransTree$branchMatchIndex[childrenNums]
    }

    transTreeBranchLengths <- sapply(matchingBranchNums, function(edgeNum) phyloAndTransTreeCopy$edge.length[[edgeNum]]$phylogeny/exp(phyloAndTransTreeCopy$edge.length[[edgeNum]]$logXi))
    childrenTimes <- currentTime + transTreeBranchLengths
    timesToConsider <- c(timesToConsider, childrenTimes)
    for (i in seq_along(childrenTimes)) {
      if (childrenNums[[i]] > length(phyloAndTransTreeCopy$tip.label)) {
        phyloAndTransTreeCopy$node.label[[childrenNums[[i]] - length(phyloAndTransTreeCopy$tip.label)]]$time <- childrenTimes[[i]]
      } # No need to update tip times...
    }
    childrenRegions <- sapply(childrenNums, function(childNum) .getVertexLabel(phyloAndTransTreeCopy, childNum)$region)
    addChildToReviewVec <- (childrenRegions == subtreeRegion) & (childrenNums > ape::Ntip(phyloAndTransTreeCopy))
    childrenStatus <- ifelse(addChildToReviewVec, "node", "tip") # A child in region y whose parent is in region x is the root of another subtree, and will be classified as such. At the same time, that child corresponds to a tip in the supporting tree. A tip is a tip, notwithstanding the region it belongs to.
    nodeOrTip <- c(nodeOrTip, childrenStatus)
    verticesToReview <- c(verticesToReview, childrenNums[addChildToReviewVec])
    index <- index + 1
    if (index > length(verticesToReview)) break
  }
  LambdaStart <- phyloAndTransTree$LambdaList[[subtreeIndex]]$priorMean
  if (is.null(LambdaStart)) { # We're initialising LambdaList...
    LambdaStart <- sum(nodeOrTip == "node")/(max(timesToConsider) - min(timesToConsider))
  }
  optimResult <- nloptr::lbfgs(x0 = log(LambdaStart), fn = function(x) -.funForLambdaOptim(logLambda = x,  times = timesToConsider, nodeOrTip = nodeOrTip), lower = log(LambdaStart) - 100)
  if (optimResult$convergence < 1) {
    warning("Convergence issue with the estimation of the mean of the prior for coalescence rates! \n")
  }
  exp(optimResult$par)
}

.produceDistTipsAncestorsMatrix <- function(phylogeny) {
  if (length(phylogeny$tip.label[[1]]) > 1) phylogeny <- .convertToTransTree(phylogeny)
  tipNums <- seq_along(phylogeny$tip.label)
  containerMatrix <- matrix(0, length(phylogeny$tip.label), length(phylogeny$node.label))
  for (tipNum in tipNums) {
    currentPos <- tipNum
    totalDist <- 0
    repeat {
      branchMatchIndex <- phylogeny$branchMatchIndex[[currentPos]]
      parentNum <- phylogeny$parentNumVec[[currentPos]]
      totalDist <- totalDist + phylogeny$edge.length[[branchMatchIndex]]
      containerMatrix[tipNum, parentNum - length(tipNums)] <- totalDist
      currentPos <- parentNum
      if (parentNum == (length(phylogeny$tip.label) + 1)) break
    }
  }
  containerMatrix
}

.incrementPhyloStructOnly <- function(phylogeny) {
  if ("error" %in% class(phylogeny)) return(list(error = phylogeny)) # Returning a list will allow it to be incremented
  phyloCopy <- phylogeny
  branchMatchParentMatrix <- sapply(1:(nrow(phyloCopy$edge) + 1), function(vertexNum) {
    branchMatch <- match(vertexNum, phyloCopy$edge[ , 2])
    parentNum <- phyloCopy$edge[branchMatch, 1]
    c(branchMatch, parentNum)
  })
  nodeChildrenList <- lapply(length(phyloCopy$tip.label):length(phyloCopy$edge.length) + 1, function(nodeNum) {
    as.integer(phyloCopy$edge[which(phyloCopy$edge[ , 1] == nodeNum), 2])
  })
  childrenList <- vector(mode = "list", length = length(phyloCopy$edge.length) + 1)
  childrenList[length(phyloCopy$tip.label):length(phyloCopy$edge.length) + 1] <- nodeChildrenList
  phyloCopy$tipNamesVec <- sapply(phyloCopy$tip.label, "[[", "name")
  phyloCopy$descendedTips <- phangorn::Descendants(x = phyloCopy, node = 1:(length(phyloCopy$edge.length) + 1), type = "tips")
  phyloCopy$descendedTips <- lapply(phyloCopy$descendedTips, "as.integer")
  phyloCopy$vertexOrderByDepth <- as.integer(.getVertexOrderByDepth(phyloCopy))
  phyloCopy$branchMatchIndex <- as.integer(branchMatchParentMatrix[1, ])
  phyloCopy$parentNumVec <- as.integer(branchMatchParentMatrix[2, ])
  phyloCopy$childrenNumList <- childrenList
  phyloCopy
}

.incrementPhyloRegionInfo <- function(phylogeny, clusterRegion) {
  phyloCopy <- phylogeny
  numTips <- length(phylogeny$tip.label)
  funForDescTips <- function(nodeNumber) {
    allDesc <- phyloCopy$descendedTips[[nodeNumber]]
    if (nodeNumber > length(phyloCopy$tip.label)) {
      nodeSubtree <- phyloCopy$node.label[[nodeNumber - numTips]]$subtreeIndex
      tipRegions <- sapply(allDesc, function(tipNumber) phyloCopy$tip.label[[tipNumber]]$region)
      tipSubtrees <- sapply(allDesc, function(tipNumber) phyloCopy$tip.label[[tipNumber]]$subtreeIndex)
      output <- allDesc[(tipSubtrees == nodeSubtree) & (tipRegions == clusterRegion)]
    } else {
      output <- nodeNumber
    }
    output
  }
  phyloCopy$descendedTips <- lapply(1:(length(phyloCopy$edge.length) + 1), FUN = funForDescTips)
  # phyloCopy$tipNumsInSubtree <- lapply(seq_along(phyloCopy$LambdaList), function(subtreeIndex) {
  #   which(sapply(phyloCopy$tip.label, "[[", "subtreeIndex") == subtreeIndex)
  # })
  phyloCopy$subtreeIndexVec <- c(sapply(phyloCopy$tip.label, "[[", "subtreeIndex"), sapply(phyloCopy$node.label, "[[", "subtreeIndex"))
  phyloCopy$vertexRegionVec <- c(sapply(phyloCopy$tip.label, "[[", "region"), sapply(phyloCopy$node.label, "[[", "region"))

  tipRegions <- sapply(phyloCopy$tip.label, "[[", "region")
  tipRegionsCheck <- tipRegions == clusterRegion
  tipSubtrees <- sapply(phyloCopy$tip.label, "[[", "subtreeIndex")

  phyloCopy$tipsInRegionBySubtree <- lapply(1:max(phyloCopy$subtreeIndexVec), function(subtreeIndex) {
    which((tipSubtrees == subtreeIndex) & tipRegionsCheck)
  })

  # This identifies the "tip" numbers for each subtree. An internal node may be a tip for a given subtree.
  numTips <- length(phyloCopy$tip.label)
  # getTipNumbersBySubtree <- function(subtreeIndex) {
  #   startNode <- phyloCopy$LambdaList[[subtreeIndex]]$rootNodeNum
  #   nodeNumsForSubtree <- numeric(0)
  #   nodesToCheck <- startNode
  #   repeat {
  #     nodeNum <- nodesToCheck[[1]]
  #     nodesToCheck <- nodesToCheck[-1]
  #     childrenNum <- phyloCopy$childrenNumList[[nodeNum]]
  #     nodeNumsToAddTest <- sapply(childrenNum, function(childNum) {
  #       childSubtree <- NULL
  #       if (childNum <= numTips) {
  #         childSubtree <- phyloCopy$tip.label[[childNum]]$subtreeIndex
  #       } else {
  #         childSubtree <- phyloCopy$node.label[[childNum - numTips]]$subtreeIndex
  #       }
  #       (childNum > numTips) & (childSubtree == subtreeIndex)
  #     })
  #     nodesToCheck <- c(nodesToCheck, childrenNum[nodeNumsToAddTest])
  #     nodeNumsForSubtree <- c(nodeNumsForSubtree, childrenNum[!nodeNumsToAddTest])
  #     if (length(nodesToCheck) == 0) break
  #   }
  #   nodeNumsForSubtree
  # }
  # phyloCopy$tipNumbersBySubtree <- lapply(seq_along(phyloCopy$LambdaList), getTipNumbersBySubtree) # The function identifies tips of the *subtree* without having to break up the complete tree into separate components. In this case, a tip either corresponds to a tip in the complete tree, or to an internal node belonging to another subtree supported by a parent that belongs to the subtree numbered subtreeIndex, which represents an introduction of the virus into a new region.
  phyloCopy$nodesInSubtreeNums <- lapply(seq_along(phyloCopy$LambdaList), function(subtreeIndex) {
    which(sapply(phyloCopy$node.label, "[[", "subtreeIndex") == subtreeIndex) + numTips
  })

  phyloCopy
}

produceClusters <- function(clusMembershipList, control) {
  cat("Summarising cluster configurations... \n")
  summaryMat <- getSumMatRcpp(clusMembershipList$clusMemVecList)/length(clusMembershipList$clusMemVecList)
  hashClusInd <- sapply(clusMembershipList$clusMemVecList, digest::digest)
  # names(clusMembershipList$clusMemVecList) <- hashClusInd
  frequencies <- table(hashClusInd)
  elementToExtract <- match(names(frequencies)[[which.max(frequencies)]], hashClusInd)
  MAPclusters <- clusMembershipList$clusMemVecList[[elementToExtract]]
  names(MAPclusters) <- clusMembershipList$seqNames
  cat("The MAP configuration was produced in", max(frequencies), "trees out of", length(clusMembershipList$clusMemVecList), "\n")
  rownames(summaryMat) <- colnames(summaryMat) <- clusMembershipList$seqNames
  reorderedSummaryMatAndHclustObj <- reorder_cormat(summaryMat, method = control$hclustMethod) # Involves a temporary switch to a dense matrix. Should work if number of sequences to cluster is under 5,000.

  hierCluster <- cutree(reorderedSummaryMatAndHclustObj$hclustObject, h = 1 - control$linkageRequirement)
  cat("Done \n")
  list(
    adjMatrix = summaryMat,
    hclustObject = reorderedSummaryMatAndHclustObj$hclustObject,
    MAPclusters = MAPclusters,
    hierarchicalClusters = hierCluster)
}

.getCoclusterMat <- function(clusMembershipVector) {
  clusterIndices <- sort(unique(clusMembershipVector))
  getAllCombns <- function(clusIndex) {
    indicesInClus <- which(clusMembershipVector == clusIndex)
    combinations <- matrix(rep(indicesInClus, 2), ncol = 2)
    if (length(indicesInClus) > 1) {
      combinations <- rbind(combinations, t(combn(indicesInClus, 2)))
    }
    combinations
  }
  nonZeroIndicesList <- lapply(clusterIndices, getAllCombns)
  nonZeroIndices <- do.call("rbind", nonZeroIndicesList)
  resultMatrix <- Matrix::sparseMatrix(i = nonZeroIndices[ , 1], j = nonZeroIndices[ , 2], dims = rep(length(clusMembershipVector), 2), symmetric = TRUE)
  colnames(resultMatrix) <- rownames(resultMatrix) <- names(clusMembershipVector)
  resultMatrix
}

reorder_cormat <- function(sparseCormat, method) {
  # Use correlation between variables as distance
  denseMatrix <- as.matrix(sparseCormat)
  dd <- as.dist(1 - denseMatrix)
  hc <- hclust(dd, method = method)
  denseMatrixReordered <- denseMatrix[hc$order, hc$order]
  list(sparseMatrix = as(denseMatrixReordered, "sparseMatrix"), hclustObject = hc)
}

# Not very fast, but only needs to be run once, which should take at most a couple of minutes
.summariseClusSizeDists <- function(clusMemVecList) {
  maxFromList <- max(sapply(clusMemVecList, function(clusVec) length(unique(clusVec))))
  envirForMatch <- new.env(size = 300) # I don't think the size argument really matters...
  sizeNameVec <- as.character(1:maxFromList)

  for (listElement in clusMemVecList) {
    freqTable <- table(table(listElement))
    freqTable <- freqTable/(sum(freqTable) * length(clusMemVecList))
    for (sizeName in names(freqTable)) {
      if (sizeName %in% ls(envirForMatch)) {
        envirForMatch[[sizeName]] <- envirForMatch[[sizeName]] + freqTable[[sizeName]]
      } else {
        assign(x = sizeName, value = freqTable[[sizeName]], envir = envirForMatch)
      }
    }
  }
  envirAsVec <- mget(x = ls(envirForMatch), envir = envirForMatch)
  objectToReturn <- do.call("c", envirAsVec)
  objectToReturn[order(as.numeric(names(objectToReturn)))]
}

obtainMLtreeRAxML <- function(alignment, raxmlCommand, directory, outgroupName, control)
{
  alignmentDigest <- digest::digest(alignment, algo = "md5", serialize = TRUE, file = FALSE, length = Inf, skip = "auto", ascii = FALSE, raw = FALSE)
  subDirectory <- paste(directory, as.character(alignmentDigest), sep = "/")
  system2(command = "mkdir", args = subDirectory)
  filename <- paste(subDirectory, "/RAxML_alignedDNA.inter", sep = "")
  ape::write.dna(x = alignment, file = filename, format = "interleaved")

  raxml1args <- paste("-s", filename, "-f d -k -m GTRCAT -c", control$numRateCats, "-n covidTree -T", control$numThreads, "-U -p", control$seedForParsemony, "-o", outgroupName, "-w", subDirectory, sep = " ")
  system2(command = raxmlCommand, args = raxml1args)

  bestTreeFile <- list.files(subDirectory, pattern = "bestTree", full.names = TRUE)

  raxml2args <- paste("-s", filename, "-f e -k -m GTRCAT -c", control$numRateCats, "-n optimizedBranches -T", control$numThreads, "-U -p", control$seedForParsemony, "-o", outgroupName, "-w", subDirectory, "-t", bestTreeFile, sep = " ")
  system2(command = raxmlCommand, args = raxml2args)

  RAxMLtree <- ape::read.tree(list.files(path = subDirectory, pattern = "result.optimizedBranches", full.names = TRUE))
  # raxml3args <- paste("-s", filename, "-f d -k -m GTRCAT -c", numRateCats, "-n treeFromSim_Boot -T", numThreads, "-U -p", seedForParsemony, "-o", outgroupName, "-x", seedForBS, "-#", numBoot,"-w", subDirectory, sep = " ")
  # system2(command = raxmlCommand, args = raxml3args)
  RAxMLtree
}

