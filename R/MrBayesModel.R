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

  # We call MrBayes
  if (!control$skipMrBayes) {
    # We write necessary files to disk
    scriptFilenameWithFolder <- .writeMrBayesFiles(DNAbinData, nexusFilename, folderForMrBayesFiles, outgroup, MrBayes.control)

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
  numItersToDrop <- floor(MrBayes.control$burninfrac * length(treeSample))
  itersToKeep <- seq(from = numItersToDrop + 1, to = length(treeSample), by = ceiling(1/control$MrBayesOutputThinningRate))
  treeSample <- treeSample[itersToKeep]
  # We also read in associated parameter values...

  parameterValuesFile <- paste(MrBayesOutputFilenamePrefix, ".p", sep = "")
  parameterValues <- .formatParameterFiles(parameterValuesFile, itersToKeep = itersToKeep, control = MrBayes.control)

  clusMembershipVecList <- .simulateClusMembership(phyloList = treeSample, targetRegion = clusterRegion, timestamps = seqsTimestampsPOSIXct, regionStamps = seqsRegionStamps, clockRate = perSiteClockRate, rootTime = estRootTime, covidCluster.control = control, priors.control = priors.control)
  output <- produceClusters(clusMembershipVecList, control = control)
  if (control$saveArgs) {
    output$args <- mget(names(formals()))
  }
  output
}

gen.covidCluster.control <- function(lengthForNullExtBranchesInPhylo = 1e-8, numReplicatesForNodeTimes = 500, numReplicatesForCoalRates = 500,  numThreads = 1, clusterCriterion = c("mrca", "cophenetic"), skipMrBayes = FALSE, MrBayesOutputThinningRate = 1, hclustMethod = "single", linkageRequirement = 0.90, saveArgs = FALSE, lambdaPriorCV = 3) {
  list(lengthForNullExtBranchesInPhylo = lengthForNullExtBranchesInPhylo, numReplicatesForNodeTimes = numReplicatesForNodeTimes,  numReplicatesForCoalRates = numReplicatesForCoalRates, numThreads = numThreads, clusterCriterion = clusterCriterion[[1]], skipMrBayes = skipMrBayes, MrBayesOutputThinningRate = MrBayesOutputThinningRate, hclustMethod = hclustMethod, linkageRequirement = linkageRequirement, saveArgs = saveArgs, lambdaPriorCV = lambdaPriorCV)
}

clusterFromMrBayesOutput <- function(seqsTimestampsPOSIXct, seqsRegionStamps, MrBayesTreesFilename, MrBayesParametersFilename = NULL, clusterRegion, clusterCriterion, burninFraction = 0.5, linkageRequirement, distLimit, epidemicRootTimePOSIXct, perSiteClockRate, control = gen.covidCluster.control()) {
  perSiteClockRate <- perSiteClockRate/365 # Time is expressed in days in the code, whereas perSiteClockRate is expressed in substitutions per site per *year*.
  estRootTime <- as.numeric(epidemicRootTimePOSIXct)/86400
  control <- do.call("gen.covidCluster.control", control)
  control$clusteringCriterion <- clusterCriterion
  control$distLimit <- distLimit
  control$clusterRegion <- clusterRegion

  if (is.null(MrBayesParametersFilename)) {
    fileRoot <- stringr::str_extract(MrBayesTreesFilename, pattern = ".*(?=\\.t$)")
    MrBayesParametersFilename <- paste(fileRoot, ".p", sep = "")
    if (!file.exists(MrBayesParametersFilename)) stop("Could not infer parameter values filename. Please enter value for argument MrBayesParametersFilename.")
  }
  treeSample <- ape::read.nexus(MrBayesTreesFilename)
  numItersToDrop <- floor(burninFraction * length(treeSample))
  itersToKeep <- seq(from = numItersToDrop + 1, to = length(treeSample), by = ceiling(1/control$MrBayesOutputThinningRate))
  treeSample <- treeSample[itersToKeep]

  parameterValues <- .formatParameterFiles(MrBayesParametersFilename, itersToKeep = itersToKeep)

  clusMembershipVecList <- .simulateClusMembership(phyloList = treeSample, targetRegion = clusterRegion, timestamps = seqsTimestampsPOSIXct, regionStamps = seqsRegionStamps, clockRate = perSiteClockRate, rootTime = estRootTime, covidCluster.control = control)
  output <- produceClusters(clusMembershipVecList, control = control)
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

gen.priors.control <- function() {
 list()
}

.simulateClusMembership <- function(phyloList, targetRegion, timestamps, regionStamps, clockRate, rootTime, covidCluster.control, priors.control) {
  timestampsInDays <- as.numeric(timestamps)/86400

  names(timestampsInDays) <- names(timestamps)
  phyloAndTransTreeList <- cl <- NULL
  if (covidCluster.control$numThreads > 1) {
    cl <- parallel::makeForkCluster(covidCluster.control$numThreads)
    phyloAndTransTreeList <- parallel::parLapply(X = phyloList, cl = cl, fun = function(listElement) tryCatch(expr = .genStartPhyloAndTransTree(phylogeny = listElement, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = log(clockRate), estRootTime = rootTime, control = covidCluster.control), error = function(e) e))
    phyloAndTransTreeList <- parallel::parLapply(X = phyloAndTransTreeList, cl = cl, fun = .incrementPhylo, clusterRegion = targetRegion)
  } else {
    phyloAndTransTreeList <- lapply(phyloList, function(listElement) tryCatch(expr = .genStartPhyloAndTransTree(phylogeny = listElement, timestampsInDays = timestampsInDays, regionStamps = regionStamps, logClockRatePriorMean = log(clockRate), estRootTime = rootTime, control = covidCluster.control), error = function(e) e))
    phyloAndTransTreeList <- lapply(phyloAndTransTreeList, .incrementPhylo, clusterRegion = targetRegion)
  }

  clusMembershipCondOnPhylo <- NULL
  if (covidCluster.control$numThreads > 1) {
    clusMembershipCondOnPhylo <- parallel::parLapply(cl = cl, X = phyloAndTransTreeList, fun = function(listElement) tryCatch(expr = .simulateClustersFromStartingTree(phyloAndTransTree = listElement, estRootTime = rootTime, control = covidCluster.control), error = function(e) e))
    parallel::stopCluster(cl)
  } else {
  clusMembershipCondOnPhylo <- lapply(phyloAndTransTreeList, FUN = function(listElement) tryCatch(expr = .simulateClustersFromStartingTree(phyloAndTransTree = listElement, estRootTime = rootTime, control = covidCluster.control), error = function(e) e))
  }
  indicesToKeep <- sapply(clusMembershipCondOnPhylo, FUN = function(listElement) !("error" %in% class(listElement)))
  combinedVecs <- do.call("c", clusMembershipCondOnPhylo[indicesToKeep])
  seqsInOrderNames <- names(combinedVecs[[1]])
  lapply(combinedVecs, function(clusMemVec) .standardiseClusterIndices(clusMemVec[seqsInOrderNames]))
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

.formatParameterFiles <- function(filenames, itersToKeep, control = NULL) {
  parameterTables <- lapply(filenames, function(filename) {
    paraTable <- read.table(file = filename, header = TRUE, sep = "\t", skip = 1)
    paraTable[itersToKeep, ]
  })
  mergedTables <- do.call("rbind", parameterTables)
  list(logLik = mergedTables$LnL, logPrior = mergedTables$LnPr, treeLength = mergedTables$TL, evoPars = mergedTables[ , -(1:4)])
}

# subtreeClusterFun is a function that takes a phyloAndTransTree object and a region index and produces a vector of cluster membership indices.

.simulateClustersFromStartingTree <- function(phyloAndTransTree,  estRootTime, control) {
  funToReplicate <- function(phyloAndTransTree, estRootTime, control) {
    medianRates <- sapply(phyloAndTransTree$LambdaList, "[[", "Lambda")
    # coalRates <- rnorm(n = length(phyloAndTransTree$LambdaList), mean = meanCoalRates, sd = meanCoalRates * control$lambdaPriorCV)
    coalRates <- rlnorm(n = length(medianRates), meanlog = log(medianRates), sdlog = abs(log(medianRates)))
    for (i in seq_along(phyloAndTransTree$LambdaList)) {
      phyloAndTransTree$LambdaList[[i]]$Lambda <- coalRates[[i]]
    }
    funToReplicateInner <- function(phyloAndTransTree, estRootTime, control) {
      numTips <- length(phyloAndTransTree$tip.label)
      orderedVertices <- rev(phyloAndTransTree$vertexOrderByDepth)
      tipTimes <- sapply(phyloAndTransTree$tip.label, "[[", "time")
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
      newTree$distTipsAncestorsMatrix <- produceDistTipsAncestorsMatrixRcpp(
        numTips = length(phyloAndTransTree$tip.label),
        numNodes = length(phyloAndTransTree$node.label),
        branchMatchIndexVec = phyloAndTransTree$branchMatchIndex - 1,
        edgeLengthsVec = nodeTimesAndEdgeLengths$edgeLengths,
        parentNumVec = phyloAndTransTree$parentNumVec - 1) # The -1 is there because this is a C++ function, and indexing starts at 0 instead of 1.
      subtreeClusterFun <- function(phyloAndTransTree, subtreeIndex, distLimit, clusterRegion, clusteringCriterion) {
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
          output <- integer(0)
          if (length(clusterList) > 0) {
            output <- .convertClusterListToVecOfIndices(clusterList, seqNames)
          }
        }
        output
      }
      clustersBySubtree <- sapply(seq_along(newTree$LambdaList), subtreeClusterFun, phyloAndTransTree = newTree, distLimit = control$distLimit, clusteringCriterion = control$clusteringCriterion, clusterRegion = control$clusterRegion)
      shiftValue <- 0
      for (i in 2:length(clustersBySubtree)) {
        if (length(clustersBySubtree[[i - 1]]) > 0) {
          shiftValue <- max(clustersBySubtree[[i - 1]])
        }
        clustersBySubtree[[i]] <- clustersBySubtree[[i]] + shiftValue
      }
      do.call("c", clustersBySubtree)
    }
    replicate(n = control$numReplicatesForCoalRates, expr = funToReplicateInner(phyloAndTransTree = phyloAndTransTree, estRootTime = estRootTime, control = control), simplify = FALSE)
  }
 do.call("c", replicate(n = control$numReplicatesForNodeTimes, funToReplicate(phyloAndTransTree, estRootTime = estRootTime, control = control), simplify = FALSE))
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
  if ("error" %in% class(phylogeny)) return(list(error = phylogeny)) # Returning a list will allow it to be incremented
  phyloCopy <- phylogeny
  numTips <- length(phylogeny$tip.label)
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
  allDescList <- phangorn::Descendants(x = phyloCopy, node = 1:(length(phyloCopy$edge.length) + 1), type = "tips")
  funForDescTips <- function(nodeNumber) {
    allDesc <- allDescList[[nodeNumber]]
    if (nodeNumber > length(phyloCopy$tip.label)) {
      nodeRegion <- phyloCopy$node.label[[nodeNumber - numTips]]$region
      nodeSubtree <- phyloCopy$node.label[[nodeNumber - numTips]]$subtreeIndex
      tipRegions <- sapply(allDesc, function(tipNumber) phyloCopy$tip.label[[tipNumber]]$region)
      tipSubtrees <- sapply(allDesc, function(tipNumber) phyloCopy$tip.label[[tipNumber]]$subtreeIndex)
      output <- allDesc[(tipSubtrees == nodeSubtree) & (tipRegions == nodeRegion)]
    } else {
      output <- nodeNumber
    }
    output
  }
  phyloCopy$descendedTips <- lapply(1:(length(phyloCopy$edge.length) + 1), FUN = funForDescTips)
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

produceClusters <- function(clusMembershipList, control) {
  adjMats <- vector("list", 0)
  if (control$numThreads > 1) {
    cl <- parallel::makeForkCluster(control$numThreads)
    adjMats <- parallel::parLapply(X = clusMembershipList, cl = cl, fun = .getCoclusterMat)
    parallel::stopCluster(cl)
  } else {
    adjMats <- lapply(clusMembershipList, .getCoclusterMat)
  }
  clusMembershipCategs <- unique(clusMembershipList)

  combinFrequencies <- table(match(clusMembershipList, clusMembershipCategs))
  MAPclusters <- clusMembershipCategs[[as.numeric(names(combinFrequencies)[[which.max(combinFrequencies)]])]]
  cat("The MAP configuration was produced in", max(combinFrequencies), "trees out of", length(clusMembershipList), ". \n")
  summaryMat <- Reduce("+", adjMats)/length(adjMats)
  reorderedSummaryMatAndHclustObj <- reorder_cormat(summaryMat, method = control$hclustMethod) # Involves a temporary switch to a dense matrix. Should work if number of sequences to cluster is under 5,000.

  hierCluster <- cutree(reorderedSummaryMatAndHclustObj$hclustObject, h = 1 - control$linkageRequirement)
  list(
    adjMatrix = reorderedSummaryMatAndHclustObj$sparseMatrix,
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
