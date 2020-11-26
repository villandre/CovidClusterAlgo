covidCluster <- function(
  DNAbinData,
  outgroup,
  clusterRegion = NULL,
  covariateFrame = NULL,
  seqsTimestampsPOSIXct = NULL,
  epidemicRootTimePOSIXct = NULL,
  folderForMrBayesFiles = NULL,
  seqsRegionStamps = rep(1, ifelse(is.matrix(DNAbinData), nrow(DNAbinData), length(DNAbinData))),
  perSiteClockRate,
  control = list(),
  MrBayes.control = gen.MrBayes.control(),
  priors.control = gen.priors.control()) {
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
  MrBayes.control <- do.call("gen.MrBayes.control", MrBayes.control)
  priors.control <- do.call("gen.priors.control", priors.control)
  nexusFilename <- "mrbayesData.nex"
  nexusFilenameWithFolder <- paste(folderForMrBayesFiles, nexusFilename, sep = "/")

  # We write necessary files to disk
  scriptFilenameWithFolder <- .writeMrBayesFiles(DNAbinData, nexusFilename, folderForMrBayesFiles, outgroup, MrBayes.control)

  # We call MrBayes

  system2(command = MrBayes.control$MrBayesShellCommand, args = scriptFilenameWithFolder)
  MrBayesOutputFilenamePrefix <- paste(nexusFilenameWithFolder, ".run", 1:MrBayes.control$nruns, sep = "")

  # We read in the tree samples and parameter files
  treeSampleFiles <- paste(MrBayesOutputFilenamePrefix, ".t", sep = "")
  treeSamples <- lapply(treeSampleFiles, ape::read.nexus)
  # We also read in associated parameter values...
  parameterValuesFiles <- paste(MrBayesOutputFilenamePrefix, ".p", sep = "")
  parameterValues <- .formatParameterFiles(parameterValuesFiles)
  mergedChains <- do.call("c", treeSamples)
  unstandardisedWeights <- (parameterValues$logLik + parameterValues$logPrior)
  standardisedWeights <- unstandardisedWeights/sum(unstandardisedWeights)
  .computeClusMembershipDistribution(phyloList = mergedChains, logWeights =  standardisedWeights, timestamps = seqsTimestampsPOSIXct, regionStamps = seqsRegionStamps, clockRate = perSiteClockRate, rootTime = estRootTime, covidCluster.control = control, priors.control = priors.control)
}

gen.covidCluster.control <- function(lengthForNullExtBranchesInPhylo = 1e-8) {
  list(lengthForNullExtBranchesInPhylo = lengthForNullExtBranchesInPhylo)
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

gen.MrBayes.control <- function(MrBayesShellCommand = "mb", nst = 6, rates = "invgamma", ngammacat = 4, nruns = 2, nchains = 4, ngen = 1000000, samplefreq = 1000, diagnfreq = 5000, printfreq = 500, burninfrac = 0.2, beaglesse = "yes", usebeagle = "no", beaglescaling = "dynamic") {
  list(MrBayesShellCommand = MrBayesShellCommand, nst = nst, rates = rates, ngammacat = ngammacat, nruns = nruns, nchains = nchains, ngen = ngen, samplefreq = samplefreq, diagnfreq = diagnfreq, printfreq = printfreq, burninfrac = burninfrac, beaglesse = beaglesse, usebeagle = usebeagle, beaglescaling = beaglescaling)
}

gen.priors.control <- function() {
 list()
}

.computeClusMembershipDistribution <- function(phyloList, logWeights, timestamps, regionStamps, clockRate, rootTime, covidCluster.control, priors.control) {
  timestampsInDays <- as.numeric(timestamps)/86400
  names(timestampsInDays) <- names(timestamps)
  phyloAndTransTreeList <- lapply(phyloList, .genStartPhyloAndTransTree, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = clockRate, estRootTime = rootTime, control = covidCluster.control)

  regionClusterFun <- function(x, regionIndex) {
    # TO_DO
  }
  funToComputeDistribCondOnPhylo <- function(index) {
    logScalingFactor <- logWeights[[index]]
    distrib <- .computeCondClusterScore(phyloAndTransTree = phyloAndTransTreeList[[index]], regionClusterFun = regionClusterFun, estRootTime = rootTime)
    for (i in seq_along(distrib)) {
      distrib[[i]]$logScore <- distrib[[i]]$logScore + logScalingFactor
    }
    distrib
  }
  distribCondOnPhyloList <- lapply(seq_along(phyloAndTransTreeList), funToComputeDistribCondOnPhylo)
  combinedDistribs <- do.call("c", distribCondOnPhyloList)

  # Now we have each element of the sum to obtain p(c | y). We must now identify all identical partitions and sum their scores to produce the final distribution, which also will be a long list whose elements are lists with two elements, config, the cluster membership vector, and logScore, the logarithm of its posterior probability.

  digestCodes <- sapply(combinedDistribs, function(listElement) digest::digest(listElement$config))
  logScores <- sapply(combinedDistribs, "[[", "logScore")
  combinedScores <- tapply(logScores, INDEX = factor(digestCodes), FUN = "logSum")
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
  paste("begin mrbayes;\n set autoclose=yes nowarn=yes;\n execute ", nexusDataFilename, ";\n lset nst=", control$nst, " rates=", control$rates, ";\n outgroup ", outgroup, ";\n set usebeagle=", control$usebeagle, " beaglescaling=", control$beaglescaling," beaglesse=", control$beaglesse, ";\n mcmc nruns=", control$nruns, " nchains=", control$nchains, " ngen=", control$nchains, " samplefreq=", control$samplefreq, " diagnfreq=", control$diagnfreq, " printfreq=", control$printfreq, " append=no;\n sump relburnin=yes burninfrac=", control$burninfrac, ";\n end;", sep = "")
}

.formatParameterFiles <- function(filenames) {
  parameterTables <- lapply(filenames, read.table, header = TRUE, sep = "\t", skip = 1)
  mergedTables <- do.call("rbind", parameterTables)
  list(logLik = mergedTables$LnL, logPrior = mergedTables$LnPr, treeLength = mergedTables$TL, evoPars = mergedTables[ , -(1:4)])
}

# regionClusterFun is a function that takes a phyloAndTransTree object and a region index and produces a vector of cluster membership indices.

.computeCondRegionClusterScore <- function(phyloAndTransTree, regionIndex, regionClusterFun, n = 5000, estRootTime, control) {
 funToReplicate <- function(phyloAndTransTree) {
   newTree <- .nodeTimesSubtreeTransFun(phyloAndTransTree, regionIndex, control)
   subtreeScore <- .nodeTimesSubtreeLogPriorFun(newTree, regionIndex, estRootTime, control)
   clustering <- regionClusterFun(newTree, regionIndex)
   list(clusters = clustering, logScore = subtreeScore)
 }
 sampledClustersAndScores <- replicate(n = n, funToReplicate(phyloAndTransTree), simplify = FALSE)
 .obtainClusMembershipDistrib(sampledClustersAndScores)
}

.obtainClusMembershipDistrib <- function(clustersAndScoresList) {
  clusterMembershipVecs <- lapply(clustersAndScoresList, function(listElement) {
    adjMat <- .getAdjacencyFromVec(listElement$clusters)
    .clusterIndicesFromAdjMat(adjMat) # This is to ensure that cluster membership vectors for identical adjacency matrices are also identical.
  })
  clusterCodes <- factor(sapply(clusterMembershipVecs, digest::digest))
  associatedScores <- sapply(clustersAndScoresList, "[[", "score")
  clusterLogScores <- tapply(associatedScores, list(clusterCodes), FUN = "computeLogSum")
  clusterMembershipVecs <- lapply(split(clusterMembershipVecs, f = clusterCodes), "[[", 1)
  lapply(seq_along(clusterMembershipVecs), function(index) list(config = clusterMembershipVecs[[index]], logScore = clusterLogScores[[index]]))
}

.computeCondClusterScore <- function(phyloAndTransTree, regionClusterFun, estRootTime, control) {
  n <- control$numReplicatesForClusMemScoring
  clustersAndScoresByRegion <- lapply(seq_along(phyloAndTransTree$LambdaList), .computeCondRegionClusterScore, phyloAndTransTree = phyloAndTransTree, regionClusterFun = regionClusterFun, n = n, estRootTime = estRootTime, control = control)
  .combineRegionalEstimatesAndScores(clustersAndScoresByRegion)
}

.combineRegionalEstimatesAndScores <- function(clusterDistribsByRegion) {
  numConfigsByRegion <- sapply(clusterDistribsByRegion, "length")
  listForExpandGrid <- lapply(numConfigsByRegion, function(x) 1:x)
  indicesToCombine <- expand.grid(listForExpandGrid)
  lapply(1:nrow(indicesToCombine), function(rowIndex) {
    elementsToSelect <- indicesToCombine[rowIndex, ]
    logScore <- sum(sapply(seq_along(elementsToSelect), function(elementIndex) {
      clusterDistribsByRegion[[elementIndex]][[elementsToSelect[[elementIndex]]]]$logScore
    }))
    combinedClusMembership <- do.call("c", lapply(seq_along(elementsToSelect), function(elementIndex) {
      clusterDistribsByRegion[[elementIndex]][[elementsToSelect[[elementIndex]]]]$config
    }))
    list(config = combinedClusMembership, logScore = logScore)
  })
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

  phyloAndTransTree <- .computeCoalRatesAddSubtreeIndex(phyloAndTransTree, estRootTime = estRootTime, control = control)
  phyloAndTransTree
}

.simulateNodeTimes <- function(phyloAndTransTree) {
  numTips <- length(phyloAndTransTree$tip.label)
  baseRatePerIntroduction <- sapply(phyloAndTransTree$LambdaList, "[[", "Lambda")
  orderedVertices <- .getVertexOrderByDepth(phyloAndTransTree, reverse = TRUE) # Represents the topology
  orderedNodes <- orderedVertices[orderedVertices > numTips]
  for (nodeNum in orderedNodes) {
    childrenNums <- NULL
    if (!is.null(phyloAndTransTree$childrenList)) {
      childrenNums <- phyloAndTransTree$childrenList[[nodeNum - numTips]]
    } else {
      childrenNums <- phyloAndTransTree$edge[which(phyloAndTransTree$edge[ , 1] == nodeNum), 2]
    }
    introForMerge <- phyloAndTransTree$node.label[[nodeNum - length(phyloAndTransTree$tip.label)]]$region
    minChildrenTimes <- min(sapply(childrenNums, function(childNum) {
      returnTime <- NULL
      if (childNum <= numTips) {
        returnTime <- phyloAndTransTree$tip.label[[childNum]]$time
      } else {
        returnTime <- phyloAndTransTree$node.label[[childNum - numTips]]$time
      }
      returnTime
    }))
    phyloAndTransTree$node.label[[nodeNum]]$time <- minChildrenTimes - rexp(n = 1, rate = baseRatePerIntroduction[[introForMerge]])
  }
  phyloAndTransTree
}
