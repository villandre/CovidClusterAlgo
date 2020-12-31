#' Clusters SARS-Cov-2 sequencing data
#'
#' The function uses MCMC to produce cluster estimates based on sequencing data and covariate information. The function automatically derives suitable starting values for all parameters.
#'
#' @param DNAbinData DNAbin object, the sequencing data,
#' @param rootSequenceName Name of the sequence used to root the tree, has to match a name in DNAbinData
#' @param clusterRegion optional string indicating region for which clusters should be computed; if left NULL, only chain results are returned
#' @param covariateFrame (ignored for now) data.frame containing covariate values in the order of elements in DNAbinData, used for scoring sample partitions
#' @param seqsTimestampsPOSIXct named vector of sequence collection times, with names matching sequence names in DNAbinData
#' @param epidemicRootTimePOSIXct POSIXct value, date on which epidemic started; affects the transmission tree topology prior; keep NULL to ignore,
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
  outgroup,
  clusterRegion = NULL,
  covariateFrame = NULL,
  seqsTimestampsPOSIXct = NULL,
  epidemicRootTimePOSIXct = NULL,
  folderForMrBayesFiles = NULL,
  seqsRegionStamps = rep(1, ifelse(is.matrix(DNAbinData), nrow(DNAbinData), length(DNAbinData))),
  perSiteClockRate = NULL,
  clusteringCriterion = c("mrca", "cophenetic", "consecutive"),
  distLimit = NULL,
  control = list(),
  MrBayes.control = gen.MrBayes.control(),
  priors.control = gen.priors.control()) {
  clusteringCriterion <- clusteringCriterion[[1]]
  outgroupPos <- match(outgroup, rownames(DNAbinData))
  rownames(DNAbinData) <- .editTaxonNames(rownames(DNAbinData))
  names(seqsTimestampsPOSIXct) <- .editTaxonNames(names(seqsTimestampsPOSIXct))
  names(seqsRegionStamps) <- .editTaxonNames(names(seqsRegionStamps))
  outgroup <- rownames(DNAbinData)[[outgroupPos]] # We fix the outgroup name since it may contain characters the MrBayes does not understand.
  perSiteClockRate <- perSiteClockRate/365 # Time is expressed in days in the code, whereas perSiteClockRate is expressed in substitutions per site per *year*.
  estRootTime <- NULL
  if (!is.null(epidemicRootTimePOSIXct)) {
    estRootTime <- as.numeric(epidemicRootTimePOSIXct)/86400
  }
  control <- do.call("gen.covidCluster.control", control)
  control$clusteringCriterion <- clusteringCriterion
  control$distLimit <- distLimit
  control$clusterRegion <- clusterRegion
  MrBayes.control <- do.call("gen.MrBayes.control", MrBayes.control)
  priors.control <- do.call("gen.priors.control", priors.control)
  nexusFilename <- "mrbayesData.nex"
  nexusFilenameWithFolder <- paste(folderForMrBayesFiles, nexusFilename, sep = "/")

  # We write necessary files to disk
  scriptFilenameWithFolder <- .writeMrBayesFiles(DNAbinData, nexusFilename, folderForMrBayesFiles, outgroup, MrBayes.control)

  # We call MrBayes
  if (!control$skipMrBayes) {
    commandForSystem2 <- MrBayes.control$MrBayesShellCommand
    argsForSystem2 <- scriptFilenameWithFolder
    if (MrBayes.control$MPIenabled) {
      commandForSystem2 <- "mpirun"
      argsForSystem2 <- c(paste("-n", control$numThreads), paste(MrBayes.control$MrBayesShellCommand, scriptFilenameWithFolder))
    }

    system2(command = commandForSystem2, args = argsForSystem2)
  }
  MrBayesOutputFilenamePrefix <- paste(nexusFilenameWithFolder, ".run1", sep = "")
  numItersToDrop <- floor(MrBayes.control$burninfrac * MrBayes.control$ngen / MrBayes.control$samplefreq)
  # We read in the tree samples and parameter files
  treeSampleFile <- paste(MrBayesOutputFilenamePrefix, ".t", sep = "")

  treeSample <- ape::read.nexus(treeSampleFile)
  itersToKeep <- seq(from = numItersToDrop + 1, to = length(treeSample), by = ceiling(1/control$MrBayesOutputThinningRate))
  treeSample <- treeSample[itersToKeep]
  # We also read in associated parameter values...

  parameterValuesFile <- paste(MrBayesOutputFilenamePrefix, ".p", sep = "")
  parameterValues <- .formatParameterFiles(parameterValuesFile, MrBayes.control)

  unstandardisedLogWeights <- (parameterValues$logLik + parameterValues$logPrior)
  standardisedLogWeights <- unstandardisedLogWeights - computeLogSum(unstandardisedLogWeights)

  # .computeClusMembershipDistribution(phyloList = mergedChains, targetRegion = clusterRegion, logWeights = standardisedLogWeights, timestamps = seqsTimestampsPOSIXct, regionStamps = seqsRegionStamps, clockRate = perSiteClockRate, rootTime = estRootTime, covidCluster.control = control, priors.control = priors.control)
  clusMembershipAndWeights <- .computeClusMembershipDistribution(phyloList = treeSample, targetRegion = clusterRegion, logWeights = standardisedLogWeights, timestamps = seqsTimestampsPOSIXct, regionStamps = seqsRegionStamps, clockRate = perSiteClockRate, rootTime = estRootTime, covidCluster.control = control, priors.control = priors.control)
  clusterAndReturnPlotObject(clusMembershipAndWeights, control = control)

}

gen.covidCluster.control <- function(lengthForNullExtBranchesInPhylo = 1e-8, numReplicatesForClusMemScoring = 5000, clusIndReportDomainSize = 10, numThreads = 1, clusterCriterion = c("mrca", "cophenetic"), skipMrBayes = FALSE, MrBayesOutputThinningRate = 1, hclustMethod = "ward.D", linkageRequirement = 0.90) {
  list(lengthForNullExtBranchesInPhylo = lengthForNullExtBranchesInPhylo, numReplicatesForClusMemScoring = numReplicatesForClusMemScoring, clusIndReportDomainSize = clusIndReportDomainSize, numThreads = numThreads, clusterCriterion = clusterCriterion[[1]], skipMrBayes = skipMrBayes, MrBayesOutputThinningRate = MrBayesOutputThinningRate, hclustMethod = hclustMethod, linkageRequirement = linkageRequirement)
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

gen.priors.control <- function() {
 list()
}

.computeClusMembershipDistribution <- function(phyloList, targetRegion, logWeights, timestamps, regionStamps, clockRate, rootTime, covidCluster.control, priors.control) {
  timestampsInDays <- as.numeric(timestamps)/86400
  names(timestampsInDays) <- names(timestamps)
  phyloAndTransTreeList <- cl <- NULL
  if (covidCluster.control$numThreads > 1) {
    cl <- parallel::makeForkCluster(covidCluster.control$numThreads)
    phyloAndTransTreeList <- parallel::parLapply(X = phyloList, cl = cl, fun = .genStartPhyloAndTransTree, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = log(clockRate), estRootTime = rootTime, control = covidCluster.control)
    phyloAndTransTreeList <- parallel::parLapply(X = phyloAndTransTreeList, cl = cl, fun = .incrementPhylo, clusterRegion = targetRegion)
  } else {
    phyloAndTransTreeList <- lapply(phyloList, .genStartPhyloAndTransTree, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = log(clockRate), estRootTime = rootTime, control = covidCluster.control)
    phyloAndTransTreeList <- lapply(phyloAndTransTreeList, .incrementPhylo, clusterRegion = targetRegion)
  }

  for (i in seq_along(logWeights)) {
    phyloAndTransTreeList[[i]]$logWeight <- logWeights[[i]]
  }
  funToComputeDistribCondOnPhylo <- function(phyloAndTransTree, estRootTime, control) {
    logScalingFactor <- phyloAndTransTree$logWeight
    distrib <- .computeCondClusterScore(phyloAndTransTree = phyloAndTransTree, estRootTime = estRootTime, control = control)
    for (i in seq_along(distrib)) {
      distrib[[i]]$logScore <- distrib[[i]]$logScore + logScalingFactor
    }
    distrib
  }
  if (covidCluster.control$numThreads > 1) {
    distribCondOnPhyloList <- parallel::parLapply(cl = cl, X = phyloAndTransTreeList, fun = funToComputeDistribCondOnPhylo, estRootTime = rootTime, control = covidCluster.control)
    parallel::stopCluster(cl)
  } else {
  distribCondOnPhyloList <- lapply(phyloAndTransTreeList, funToComputeDistribCondOnPhylo, estRootTime = rootTime, control = covidCluster.control)
  }
  combinedDistribs <- do.call("c", distribCondOnPhyloList)
  namesForOrder <- names(combinedDistribs[[1]]$config)
  for (i in seq_along(combinedDistribs)) {
    combinedDistribs[[i]]$config <- .standardiseClusterIndices(combinedDistribs[[i]]$config[namesForOrder])
  }
  # Now we have each element of the sum to obtain p(c | y). We must now identify all identical partitions and sum their scores to produce the final distribution, which also will be a long list whose elements are lists with two elements, config, the cluster membership vector, and logScore, the logarithm of its posterior probability.

  digestCodes <- sapply(combinedDistribs, function(listElement) digest::digest(listElement$config))
  logScores <- sapply(combinedDistribs, "[[", "logScore")
  combinedScores <- tapply(logScores, INDEX = factor(digestCodes), FUN = "computeLogSum")
  clusterMembershipVecs <- lapply(combinedDistribs, "[[", "config")
  uniqueClusMembershipVecs <- lapply(split(clusterMembershipVecs, f = factor(digestCodes)), "[[", 1)
  lapply(seq_along(uniqueClusMembershipVecs), function(index) list(config = uniqueClusMembershipVecs[[index]], logScore = combinedScores[[index]]))
}

.writeMrBayesFiles <- function(DNAbinData, nexusFilename, folderForMrBayesFiles, outgroup, control) {
  nexusFilenameWithFolder <- paste(folderForMrBayesFiles, nexusFilename, sep = "/")
  scriptForMrBayes <- .produceMrBayesScript(outgroup = outgroup, nexusDataFilename = nexusFilenameWithFolder, control = control)

  ape::write.nexus.data(x = DNAbinData, file = nexusFilenameWithFolder)
  scriptFilenameWithFolder <- paste(folderForMrBayesFiles, "MrBayesScript.nex", sep = "/")
  fileConn <- file(scriptFilenameWithFolder)
  writeLines(scriptForMrBayes, fileConn)
  close(fileConn)
  scriptFilenameWithFolder
}

.produceMrBayesScript <- function(outgroup, nexusDataFilename, control = gen.MrBayes.control()) {
  paste("begin mrbayes;\n set autoclose=yes nowarn=yes seed=", as.integer(control$seed), " swapseed=", as.integer(control$swapseed), ";\n execute ", nexusDataFilename, ";\n lset nst=", as.integer(control$nst), " rates=", control$rates, ";\n outgroup ", outgroup, ";\n set usebeagle=", control$usebeagle, " beaglescaling=", control$beaglescaling," beaglesse=", control$beaglesse, ";\n mcmc nruns=", as.integer(control$nruns), " nchains=", as.integer(control$nchains), " ngen=", as.integer(control$ngen), " samplefreq=", as.integer(control$samplefreq), " diagnfreq=", as.integer(control$diagnfreq), " printfreq=", as.integer(control$printfreq), " append=", control$append, " temp=",  control$temp, ";\n sump relburnin=yes burninfrac=", control$burninfrac, ";\n end;", sep = "")
}

.formatParameterFiles <- function(filenames, control) {
  numItersToDrop <- floor(control$burninfrac * control$ngen / control$samplefreq)
  parameterTables <- lapply(filenames, function(filename) {
    paraTable <- read.table(file = filename, header = TRUE, sep = "\t", skip = 1)
    paraTable[-(1:numItersToDrop), ]
  })
  mergedTables <- do.call("rbind", parameterTables)
  list(logLik = mergedTables$LnL, logPrior = mergedTables$LnPr, treeLength = mergedTables$TL, evoPars = mergedTables[ , -(1:4)])
}

# subtreeClusterFun is a function that takes a phyloAndTransTree object and a region index and produces a vector of cluster membership indices.

.computeCondClusterScoreBySubtree <- function(phyloAndTransTree, n = 5000, estRootTime, control) {

  funToReplicate <- function(phyloAndTransTree, estRootTime, control) {
    numTips <- length(phyloAndTransTree$tip.label)
    coalRates <- sapply(phyloAndTransTree$LambdaList, "[[", "Lambda")
    orderedVertices <- rev(phyloAndTransTree$vertexOrderByDepth)
    tipTimes <- sapply(phyloAndTransTree$tip.label, "[[", "time")
   # newTree <- .simulateNodeTimes(phyloAndTransTree)
    nodeTimesAndEdgeLengths <- simulateNodeTimesRcpp(
      numTips = numTips,
      baseRatePerIntroduction = coalRates,
      orderedVertices = orderedVertices,
      subtreeIndexVec = phyloAndTransTree$subtreeIndexVec,
      tipTimes = tipTimes,
      edgeMatrix = phyloAndTransTree$edge,
      childrenNumList = phyloAndTransTree$childrenNumList)
    newTree <- phyloAndTransTree
    for (i in seq_along(newTree$node.label)) {
      newTree$node.label[[i]]$time <- nodeTimesAndEdgeLengths$vertexTimes[[i]]
    }
    for (i in seq_along(newTree$edge.length)) {
      newTree$edge.length[[i]]$transmissionTree <- nodeTimesAndEdgeLengths$edgeLengths[[i]]
    }
    # newTree <- .addTransTreeEdgeLengths(newTree)
    # edgeLengths <- sapply(newTree$edge.length, "[[", "transmissionTree")
    newTree$distTipsAncestorsMatrix <- produceDistTipsAncestorsMatrixRcpp(
     numTips = length(phyloAndTransTree$tip.label),
     numNodes = length(phyloAndTransTree$node.label),
     branchMatchIndexVec = phyloAndTransTree$branchMatchIndex - 1,
     edgeLengthsVec = nodeTimesAndEdgeLengths$edgeLengths,
     parentNumVec = phyloAndTransTree$parentNumVec - 1) # The -1 is there because this is a C++ function, and indexing starts at 0 instead of 1.
   subtreeClusterFun <- function(phyloAndTransTree, subtreeIndex, distLimit, clusterRegion, clusteringCriterion) {
     numTips <- length(phyloAndTransTree$tip.label)
     # tipsInSubtreeAndRegion <- which((sapply(phyloAndTransTree$tip.label, "[[", "region") == clusterRegion) & sapply(phyloAndTransTree$tip.label, "[[", "subtreeIndex") == subtreeIndex)
     tipsInSubtreeAndRegion <- phyloAndTransTree$tipsInRegionBySubtree[[subtreeIndex]]
     seqNames <- sapply(tipsInSubtreeAndRegion, function(tipNum) phyloAndTransTree$tip.label[[tipNum]]$name)
     output <- seq_along(seqNames)
     names(output) <- seqNames
     if (identical(phyloAndTransTree$node.label[[phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum - numTips]]$region, clusterRegion)) {
       clusterFunction <- getMRCAclustersRcpp
       if (control$clusterCriterion == "cophenetic") {
         clusterFunction <- getCopheneticClustersRcpp
       }
       # clusterList <- getMRCAclustersRcpp(
       clusterList <- clusterFunction(
         parentNumVec = phyloAndTransTree$parentNumVec,
         childrenNumList = phyloAndTransTree$childrenNumList,
         descendedTipsList = phyloAndTransTree$descendedTips,
         subtreeIndexVec = phyloAndTransTree$subtreeIndexVec,
         vertexRegionVec = phyloAndTransTree$vertexRegionVec,
         tipNamesVec = phyloAndTransTree$tipNamesVec,
         subtreeRootNum = phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum,
         distTipsAncestorsMatrix = phyloAndTransTree$distTipsAncestorsMatrix,
         subtreeIndex = subtreeIndex,
         numTips = numTips,
         regionLabel = clusterRegion,
         distLimit = distLimit)
       # clusterList <- getDistanceBasedClusters(phyloAndTransTree = phyloAndTransTree, subtreeIndex = subtreeIndex, distLimit = distLimit, regionLabel = clusterRegion, criterion = clusteringCriterion)
       output <- integer(0)
       if (length(clusterList) > 0) {
         output <- .convertClusterListToVecOfIndices(clusterList, seqNames)
       }
     }
     output
   }
   clustersBySubtree <- sapply(seq_along(newTree$LambdaList), subtreeClusterFun, phyloAndTransTree = newTree, distLimit = control$distLimit, clusteringCriterion = control$clusteringCriterion, clusterRegion = control$clusterRegion)
   subtreeScoresBySubtree <- sapply(seq_along(newTree$LambdaList), function(subtreeIndex) {
     output <- NA
     if (length(clustersBySubtree[[subtreeIndex]]) > 0) {
      output <- .nodeTimesSubtreeLogPriorFun(phyloAndTransTree = newTree, estRootTime = estRootTime, subtreeIndex = subtreeIndex, control = control)
     }
     output
   })
   lapply(seq_along(clustersBySubtree), function(index) list(clusters = clustersBySubtree[[index]], logScore = subtreeScoresBySubtree[[index]]))
 }
 sampledClustersAndScores <- replicate(n = control$numReplicatesForClusMemScoring, funToReplicate(phyloAndTransTree, estRootTime = estRootTime, control = control), simplify = FALSE)

 sampledClustersAndScoresBySubtree <- lapply(seq_along(sampledClustersAndScores[[1]]), function(subtreeIndex) {
   lapply(seq_along(sampledClustersAndScores), function(replicateIndex) list(clusters = sampledClustersAndScores[[replicateIndex]][[subtreeIndex]]$clusters, logScore = sampledClustersAndScores[[replicateIndex]][[subtreeIndex]]$logScore))
 })
 lapply(sampledClustersAndScoresBySubtree, .obtainClusMembershipDistrib)
}

.obtainClusMembershipDistrib <- function(clustersAndScoresList) {
  clusterMembershipVecs <- lapply(clustersAndScoresList, function(listElement) {
    .standardiseClusterIndices(listElement$clusters) # As elements of the list are replicates, re-ordering listElement$clusters should not be required.
  })
  clusterCodes <- factor(sapply(clusterMembershipVecs, digest::digest))
  associatedScores <- sapply(clustersAndScoresList, "[[", "logScore")
  clusterLogScores <- tapply(associatedScores, list(clusterCodes), FUN = "computeLogSum")
  clusterMembershipVecs <- lapply(split(clusterMembershipVecs, f = clusterCodes), "[[", 1)
  lapply(seq_along(clusterMembershipVecs), function(index) list(config = clusterMembershipVecs[[index]], logScore = clusterLogScores[[index]]))
}

.computeCondClusterScore <- function(phyloAndTransTree, estRootTime, control) {
  n <- control$numReplicatesForClusMemScoring
  clustersAndScoresBySubtree <- .computeCondClusterScoreBySubtree(phyloAndTransTree = phyloAndTransTree, n = n, estRootTime = estRootTime, control = control)
  .combineRegionalEstimatesAndScores(clustersAndScoresBySubtree, control = control)
}

.combineRegionalEstimatesAndScores <- function(clusterDistribsBySubtree, control) {
  subtreesToKeep <- which(sapply(clusterDistribsBySubtree, function(listElement) {
    length(listElement[[1]]$config) > 0
  }))
  logScoresToBySubtree <- lapply(clusterDistribsBySubtree[subtreesToKeep], function(listElement) sapply(listElement, "[[" ,"logScore"))
  listWithIndicesAndLogScores <- .sortNumberListDecreasing(logScoresToBySubtree, control$clusIndReportDomainSize) # Each element of the list is itself a list, with two elements: a vector of indices (giving the vector of indices to select for each subtree to form the complete cluster membership indices vector) and the associated log-score.
  concatenateClusInd <- function(listElement) {
    clusIndList <- lapply(seq_along(subtreesToKeep), function(i) clusterDistribsBySubtree[[subtreesToKeep[[i]]]][[listElement$configIndexVec[[i]]]]$config)
    list(config = .combineClusMembershipIndices(clusIndList), logScore = listElement$logScore)
  }
  lapply(listWithIndicesAndLogScores, concatenateClusInd)
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

.genStartPhyloAndTransTree <- function(phylogeny, timestampsInDays, regionStamps, logClockRatePriorMean, estRootTime, control) {
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
  phyloAndTransTree <- .identifyNodeRegions(phyloAndTransTree)
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
  phyloAndTransTree$LambdaList <- lapply(seq_along(subtreeRootNodes), function(i) list(Lambda = NULL, rootNodeNum = subtreeRootNodes[[i]]))
  coalRates <- sapply(seq_along(phyloAndTransTree$LambdaList), .computeCoalRateSubtree, phyloAndTransTree = phyloAndTransTree)

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
    # childrenNums <- phyloAndTransTree$edge[which(phyloAndTransTree$edge[ , 1] == nodeNum), 2]
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
      childrenNums <- phangorn::Children(phyloAndTransTree, vertexNum)
      matchingBranchNums <- match(childrenNums, phyloAndTransTreeCopy$edge[ , 2])
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
    childrenStatus <- ifelse(addChildToReviewVec, "node", "tip") # A child in region y whose parent is in region x is the root of another subtree, and will be classified as such. At the same time, that child corresponds to a tip in the supporting tree. A tip is a tip, notwithstanding the region it belongs to. This will affect the computation of the prior mean.
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

.combineClusMembershipIndices <- function(clusIndList) {
  modClusIndList <- lapply(clusIndList, function(x) {
    output <- as.numeric(as.factor(x))
    names(output) <- names(x)
    output
  })
  offsetValue <- cumsum(sapply(modClusIndList[1:(length(modClusIndList) - 1)], max))
  for (i in 2:length(modClusIndList)) {
    modClusIndList[[i]] <- modClusIndList[[i]] + offsetValue[[i - 1]]
  }
  do.call("c", modClusIndList)
}

# We assume that the vector elements of listToSort do not contain ties.

.sortNumberListDecreasing <- function(listToSort, numConfigs) {
  if (all(sapply(listToSort, length) == 1)) {
    return(list(list(configIndexVec = rep(1, length(listToSort)), logScore = Reduce("+", listToSort))))
  }
  numConfigs <- min(numConfigs, prod(sapply(listToSort, "length")))

  configList <- accessibleValues <- vector(mode = "list", length = numConfigs)
  configList[[1]] <- sapply(listToSort, which.max)
  scoresVec <- numeric(numConfigs)
  scoresVec[[1]] <- sum(sapply(seq_along(configList[[1]]), function(i) listToSort[[i]][[configList[[1]][[i]]]]))

  for (i in 2:numConfigs) {
    funForApplyFun <- function(subtreeIndex) {
      originalValue <- listToSort[[subtreeIndex]][[configList[[i - 1]][[subtreeIndex]]]]
      potentialValues <- listToSort[[subtreeIndex]]
      scoreIncrements <- potentialValues - originalValue
      scoreIncrements <- replace(scoreIncrements, which(scoreIncrements >= 0), -Inf) # To ensure we don't return those states anymore.
      scoresVec[[i - 1]] + scoreIncrements
    }
    accessibleValues[[i - 1]] <- lapply(seq_along(configList[[i - 1]]), funForApplyFun)
    maxValueFromEachConfig <- sapply(accessibleValues[1:(i - 1)], function(listElement) max(sapply(listElement, max)))
    configFromWhichToMoveIndex <- which.max(maxValueFromEachConfig)
    subtreeToChangeIndex <- which.max(sapply(accessibleValues[[configFromWhichToMoveIndex]], max))
    elementToMoveToIndex <- which.max(accessibleValues[[configFromWhichToMoveIndex]][[subtreeToChangeIndex]])
    configList[[i]] <- configList[[configFromWhichToMoveIndex]]
    configList[[i]][[subtreeToChangeIndex]] <- elementToMoveToIndex
    scoresVec[[i]] <- accessibleValues[[configFromWhichToMoveIndex]][[subtreeToChangeIndex]][[elementToMoveToIndex]]
    accessibleValues[[configFromWhichToMoveIndex]][[subtreeToChangeIndex]][[elementToMoveToIndex]] <- -Inf
  }
  lapply(seq_along(configList), function(i) {
    list(configIndexVec = configList[[i]], logScore = scoresVec[[i]])
  })
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

.incrementPhylo <- function(phylogeny, clusterRegion) {
  phyloCopy <- phylogeny
  branchMatchParentMatrix <- sapply(1:(nrow(phyloCopy$edge) + 1), function(vertexNum) {
    branchMatch <- match(vertexNum, phyloCopy$edge[ , 2])
    parentNum <- phyloCopy$edge[branchMatch, 1]
    c(branchMatch, parentNum)
  })
  nodeChildrenList <- lapply(length(phyloCopy$tip.label):length(phyloCopy$edge.length) + 1, function(nodeNum) {
    phyloCopy$edge[which(phyloCopy$edge[ , 1] == nodeNum), 2]
  })
  childrenList <- vector(mode = "list", length = length(phyloCopy$edge.length) + 1)
  childrenList[length(phyloCopy$tip.label):length(phyloCopy$edge.length) + 1] <- nodeChildrenList
  phyloCopy$tipNumsInSubtree <- lapply(seq_along(phyloCopy$LambdaList), function(subtreeIndex) {
    which(sapply(phyloCopy$tip.label, "[[", "subtreeIndex") == subtreeIndex)
  })
  phyloCopy$subtreeIndexVec <- c(sapply(phyloCopy$tip.label, "[[", "subtreeIndex"), sapply(phyloCopy$node.label, "[[", "subtreeIndex"))
  phyloCopy$vertexRegionVec <- c(sapply(phyloCopy$tip.label, "[[", "region"), sapply(phyloCopy$node.label, "[[", "region"))
  phyloCopy$tipNamesVec <- sapply(phyloCopy$tip.label, "[[", "name")
  phyloCopy$descendedTips <- phangorn::Descendants(x = phyloCopy, node = 1:(length(phyloCopy$edge.length) + 1), type = "tips")
  phyloCopy$vertexOrderByDepth <- .getVertexOrderByDepth(phyloCopy)
  phyloCopy$branchMatchIndex <- branchMatchParentMatrix[1, ]
  phyloCopy$parentNumVec <- branchMatchParentMatrix[2, ]
  phyloCopy$childrenNumList <- childrenList

  tipRegions <- sapply(phyloCopy$tip.label, "[[", "region")
  tipRegionsCheck <- tipRegions == clusterRegion
  tipSubtrees <- sapply(phyloCopy$tip.label, "[[", "subtreeIndex")

  phyloCopy$tipsInRegionBySubtree <- lapply(1:max(phyloCopy$subtreeIndexVec), function(subtreeIndex) {
    which((tipSubtrees == subtreeIndex) & tipRegionsCheck)
  })

  # This identifies the "tip" numbers for each subtree. An internal node may be a tip for a given subtree.
  numTips <- length(phyloCopy$tip.label)
  getTipNumbersBySubtree <- function(subtreeIndex) {
    startNode <- phyloCopy$LambdaList[[subtreeIndex]]$rootNodeNum
    nodeNumsForSubtree <- numeric(0)
    nodesToCheck <- startNode
    repeat {
      nodeNum <- nodesToCheck[[1]]
      nodesToCheck <- nodesToCheck[-1]
      childrenNum <- phyloCopy$childrenNumList[[nodeNum]]
      nodeNumsToAddTest <- sapply(childrenNum, function(childNum) {
        childSubtree <- NULL
        if (childNum <= numTips) {
          childSubtree <- phyloCopy$tip.label[[childNum]]$subtreeIndex
        } else {
          childSubtree <- phyloCopy$node.label[[childNum - numTips]]$subtreeIndex
        }
        (childNum > numTips) & (childSubtree == subtreeIndex)
      })
      nodesToCheck <- c(nodesToCheck, childrenNum[nodeNumsToAddTest])
      nodeNumsForSubtree <- c(nodeNumsForSubtree, childrenNum[!nodeNumsToAddTest])
      if (length(nodesToCheck) == 0) break
    }
    nodeNumsForSubtree
  }
  phyloCopy$tipNumbersBySubtree <- lapply(seq_along(phyloCopy$LambdaList), getTipNumbersBySubtree) # The function identifies tips of the *subtree* without having to break up the complete tree into separate components. In this case, a tip either corresponds to a tip in the complete tree, or to an internal node belonging to another subtree supported by a parent that belongs to the subtree numbered subtreeIndex, which represents an introduction of the virus into a new region.
  phyloCopy$nodesInSubtreeNums <- lapply(seq_along(phyloCopy$LambdaList), function(subtreeIndex) {
    which(sapply(phyloCopy$node.label, "[[", "subtreeIndex") == subtreeIndex) + numTips
  })

  phyloCopy
}

clusterAndReturnPlotObject <- function(clusMembershipListAndWeights, control) {
  adjMats <- lapply(clusMembershipListAndWeights, function(listElement) .getCoclusterMat(listElement$config))
  logScoresVec <- sapply(clusMembershipListAndWeights, "[[", "logScore")
  logStandardConstant <- computeLogSum(logScoresVec)
  logScoresVecStandard <- logScoresVec - logStandardConstant
  summaryMat <- Reduce("+", mapply(logScore = logScoresVecStandard, adjMat = adjMats, FUN = function(logScore, adjMat) exp(logScore) * adjMat))
  reorderedSummaryMatAndHclustObj <- reorder_cormat(summaryMat, method = control$hclustMethod) # Involves a temporary switch to a dense matrix. Should work if number of sequences to cluster is under 5,000.
  seqNames <- names(clusMembershipListAndWeights[[1]]$config)
  matrixToPlot <- as.data.frame(Matrix::mat2triplet(reorderedSummaryMatAndHclustObj$sparseMatrix))
  colnames(matrixToPlot) <- c("Sequence_x", "Sequence_y", "CoclusteringRate")
  matrixToPlot$Sequence_x <- factor(matrixToPlot$Sequence_x, levels = seq_along(seqNames), labels = seqNames)
  matrixToPlot$Sequence_y <- factor(matrixToPlot$Sequence_y, levels = seq_along(seqNames), labels = seqNames)

  plotObject <-
    ggplot2::ggplot(data = matrixToPlot, ggplot2::aes(x = Sequence_x, y = Sequence_y, fill = CoclusteringRate)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "white", high = "red", limit = c(0,1), space = "Lab", name = "Coclustering\nrate") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 4, hjust = 1)) +
    ggplot2::coord_fixed()
  hierCluster <- cutree(reorderedSummaryMatAndHclustObj$hclustObject, h = 1 - control$linkageRequirement)
  list(
    adjMatrix = reorderedSummaryMatAndHclustObj$sparseMatrix,
    MAPclusters = clusMembershipListAndWeights[[which.max(logScoresVecStandard)]]$config,
    hierarchicalClusters = , objectToPlot = plotObject)
}

.getCoclusterMat <- function(clusMembershipVector) {
  clusterIndices <- sort(unique(clusMembershipVector))
  getAllCombns <- function(clusIndex) {
    indicesInClus <- which(clusMembershipVector == clusIndex)
    combinations <- rep(indicesInClus[[1]], 2)
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
