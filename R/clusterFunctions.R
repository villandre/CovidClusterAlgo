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

findBayesianClusters <- function(
  DNAbinData,
  rootSequenceName,
  clusterRegion = NULL,
  covariateFrame = NULL,
  seqsTimestampsPOSIXct,
  epidemicRootTimePOSIXct = NULL,
  seqsRegionStamps = rep(1, ifelse(is.matrix(DNAbinData), nrow(DNAbinData), length(DNAbinData))),
  perSiteClockRate,
  startingValuesFromChainIter = NULL,
  evoParsList = NULL,
  clusterScoringFun = NULL,
  control = list()) {
  # .performBasicChecks()
  DNAbinData <- as.matrix(DNAbinData)
  perSiteClockRate <- perSiteClockRate/365 # Time is expressed in days in the code, whereas perSiteClockRate is expressed in substitutions per site per *year*.
  # Converting to value in days.
  if (is.null(control$lengthForNullExtBranchesInPhylo)) control$lengthForNullExtBranchesInPhylo <- perSiteClockRate/100
  control <- do.call(findBayesianClusters.control, args = control)

  # For backwards compatibility
  if ("clockRatesLogKernelSD" %in% names(control$MCMC.control)) {
    names(control$MCMC.control) <- replace(names(control$MCMC.control), match("clockRatesLogKernelSD", names(control$MCMC.control)), "logClockRatesKernelSD")
  }
  control$MCMC.control <- do.call(MCMC.control, args = control$MCMC.control)
  estRootTime <- NULL
  if (!is.null(epidemicRootTimePOSIXct)) {
    estRootTime <- as.numeric(epidemicRootTimePOSIXct)/86400
  }

  if (is.null(control$numMigrationsPoissonPriorMean)) {
    control$numMigrationsPoissonPriorMean <- length(unique(seqsRegionStamps)) - 1 # We prefer transmission trees with only one introduction per region; the -1 is for the outgroup, which starts off with a region without an introduction.
  }
  startingValues <- list()
  # The "edge.length" component of "phylogeny" is a list containing branch lengths for both the transmission tree and phylogeny.
  chainId <- chainIdCopy <- control$MCMC.control$chainId # Copy exists to distinguish a new id (which will be generated later on) and an id used previously.
  if (is.null(control$MCMC.control$folderToSaveIntermediateResults) | is.null(chainId)) {
    chainId <- stringi::stri_rand_strings(1, 6)
    if (is.null(startingValuesFromChainIter)) {
      MLphyloAndEvoPars <- .genMLphyloAndEvoPars(DNAbinData, rootSequenceName)
      startingValuePhylo <- MLphyloAndEvoPars$phylogeny
      if (is.null(evoParsList)) evoParsList <- MLphyloAndEvoPars$evoParsList
    } else {
      evoParsList <- startingValuesFromChainIter$evoParsList
    }
    if (!is.null(control$MCMC.control$folderToSaveIntermediateResults)) {
      filename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_controlParameters.Rdata", sep = "")
      save(control, file = filename)
      evoParsFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_evoParameters.Rdata", sep = "")
      save(evoParsList, file = evoParsFilename)
    }
    if (is.null(startingValuesFromChainIter)) {
      # startingValues$phyloAndTransTree <- .genStartphyloAndTransTree(phylogeny = startingValuePhylo, seqsTimestampsPOSIXct = seqsTimestampsPOSIXct, seqsRegionStamps = seqsRegionStamps, perSiteClockRate = perSiteClockRate, control = control$controlForGenStartTransmissionTree)
      timestampsInDays <- as.numeric(seqsTimestampsPOSIXct)/86400
      names(timestampsInDays) <- names(seqsTimestampsPOSIXct)
      startingValues$phyloAndTransTree <- .genStartPhyloAndTransTreeAlt(phylogeny = startingValuePhylo, timestampsInDays = timestampsInDays, regionStamps = seqsRegionStamps, logClockRatePriorMean = log(perSiteClockRate), estRootTime = estRootTime, control = control)

      startingValues$Lambda <- .genStartCoalescenceRates(startingValues$phyloAndTransTree, estRootTime = estRootTime, control = control)
    } else {
      startingValues$phyloAndTransTree <- startingValuesFromChainIter$paraValues$phyloAndTransTree
      startingValues$Lambda <- startingValuesFromChainIter$paraValues$Lambda
    }
    # logXi is now being added to edge.length within .genStartPhyloAndTransTreeAlt. The following should be deleted once it is confirmed that the function works.
    # if (!control$fixedClockRate & is.null(startingValuesFromChainIter)) {
    #   startXi <- rep(perSiteClockRate, length(startingValues$phyloAndTransTree$edge.length))
    #   startingValues$phyloAndTransTree$edge.length <- lapply(seq_along(startingValues$phyloAndTransTree$edge.length), function(branchIndex) {
    #     updatedBranch <- startingValues$phyloAndTransTree$edge.length[[branchIndex]]
    #     updatedBranch$logXi <- startXi[[branchIndex]]
    #     updatedBranch
    #   })
    # }
    if (!is.null(control$MCMC.control$folderToSaveIntermediateResults)) {
      startingValuesFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_startingValues.Rdata", sep = "")
      save(startingValues, file = startingValuesFilename)
    }
  } else {
    cat("Restoring chain", chainId, "control parameters... \n", sep = " ")
    burnin <- control$MCMC.control$burnin
    n <- control$MCMC.control$n
    stepSize <- control$MCMC.control$stepSize
    nIterPerSweep <- control$MCMC.control$nIterPerSweep

    filename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_controlParameters.Rdata", sep = "")
    evoParsFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_evoParameters.Rdata", sep = "")
    startingValuesFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_startingValues.Rdata", sep = "")
    load(filename) # This will restore "control" to its original values.
    control$MCMC.control$chainId <- chainIdCopy
    # The user should be allowed to change these chain parameters, as they do not affect the transitions.
    if (!is.null(n)) control$MCMC.control$n <- n
    if (!is.null(burnin)) control$MCMC.control$burnin <- burnin
    if (!is.null(stepSize)) control$MCMC.control$stepSize <- stepSize
    if (!is.null(nIterPerSweep)) control$MCMC.control$nIterPerSweep <- nIterPerSweep

    load(evoParsFilename) # This will restore evoParsList.
    load(startingValuesFilename) # This will restore startingValues.
  }

  sampledTreesWithPP <- .optimTreeMCMC(
    startingValues = startingValues,
    chainId = chainId,
    logLikFun = control$logLikFun,
    DNAbinData = DNAbinData,
    evoParsList = evoParsList,
    covariateFrame = covariateFrame,
    perSiteClockRate = perSiteClockRate,
    estRootTime = estRootTime,
    control = control)

  clusterValues <- NULL
  if (!is.null(clusterRegion)) {
    logPPvalues <- sapply(sampledTreesWithPP, FUN = "[[", "logPP")
    clusterValues <- getRegionClusters(phyloAndTransTree = sampledTreesWithPP[[which.max(logPPvalues)]]$paraValues$phyloAndTransTree, clusterRegionCode = clusterRegion)
  }
  outputValue <- list(chain = sampledTreesWithPP, MAPclusters = clusterValues)
  if (control$saveData) {
    outputValue$data <- DNAbinData
    outputValue$timestampsPOSIXct <- seqsTimestampsPOSIXct
    outputValue$regionStamps <- seqsRegionStamps
    outputValue$rootSequenceName <- rootSequenceName
  }
  outputValue
}

#' Control parameters for findBayesianClusters
#'
#' @param logLikFun function producing the phylogenetic log-likelihood, cf. presetPML for default values
#' @param numMigrationsPoissonPriorMean Mean of the Poisson distribution used to penalise for the number of clusters in the sample
#' @param MCMC.control control parameters for the MCMC run, cf. ?MCMC.control
#' @param transTreeCondOnPhylo Determines the factorisation of the posterior distribution; when TRUE, the transmission tree branch length prior is conditional on phylogenetic branch lengths,
#' @param controlForGenStartTransmissionTree (not used for now) control parameters for the function producing a starting value for the transmission tree.
#'
#' @details You should not need to call this function directly: it is merely a convenient and formal way to let users know how the function behaviour can be customised
#' @return A list of control parameters for findBayesianClusters.
#'
#' }
#'
#' @export
findBayesianClusters.control <- function(
  logLikFun = presetPML,
  numMigrationsPoissonPriorMean = 1,
  MCMC.control = MCMC.control(),
  transTreeCondOnPhylo = TRUE,
  controlForGenStartTransmissionTree = .controlForGenStartTransmissionTree(),
  fixedClockRate = TRUE,
  strictClockModel = TRUE,
  saveData = TRUE,
  resolvePolytomies = FALSE,
  epidemicRootTimePriorVar = 1/86400,
  lengthForNullExtBranchesInPhylo = 0) {
  list(logLikFun = logLikFun, numMigrationsPoissonPriorMean = numMigrationsPoissonPriorMean, MCMC.control = do.call("MCMC.control", MCMC.control), controlForGenStartTransmissionTree = do.call(".controlForGenStartTransmissionTree", controlForGenStartTransmissionTree), transTreeCondOnPhylo = transTreeCondOnPhylo, fixedClockRate = fixedClockRate, strictClockModel = strictClockModel, saveData = saveData, resolvePolytomies = resolvePolytomies, epidemicRootTimePriorVar = epidemicRootTimePriorVar, lengthForNullExtBranchesInPhylo = lengthForNullExtBranchesInPhylo)
}

.performBasicChecks <- function() {

}

# evoParsList is a list of varying parameters for the log-likelihood model. Control parameters are hard-coded.
presetPML <- function(phyloObj, phyDatObj, evoParsList) {
  phangorn::pml(tree = phyloObj,
                data = phyDatObj,
                bf = evoParsList$bf,
                Q = evoParsList$Q,
                inv = evoParsList$propInv,
                k = 4,
                shape = evoParsList$gammaShape,
                model = "GTR")$logLik
}


.genMLphyloAndEvoPars <- function(DNAbinData, rootSequenceName) {
  distMatrix <- ape::dist.dna(x = DNAbinData, model = "raw", pairwise.deletion = TRUE)
  startingPhylo <- ape::bionj(distMatrix)
  startPML <- phangorn::pml(tree = startingPhylo, data = phangorn::as.phyDat(DNAbinData), k = 4, model = "GTR")
  optimisedPhylo <- phangorn::optim.pml(object = startPML, optNni = TRUE, optBf = TRUE, optQ = TRUE, optInv = TRUE, optGamma = TRUE, optEdge = TRUE, optRate = TRUE)
  optimisedPhylo$tree <- ape::root(optimisedPhylo$tree, outgroup = rootSequenceName, resolve.root = TRUE)

  list(phylogeny = optimisedPhylo$tree, evoParsList = list(Q = optimisedPhylo$Q, bf = optimisedPhylo$bf, gammaShape = optimisedPhylo$shape, propInv = optimisedPhylo$inv, ncat = optimisedPhylo$k))
}

.controlForGenStartTransmissionTree <- function() {
  list()
}

# .genStartphyloAndTransTree <- function(phylogeny, seqsTimestampsPOSIXct, seqsRegionStamps, perSiteClockRate, estRootTime = NULL, control) {
#   if (length(perSiteClockRate) == 0) {
#     timestamps <- seqsTimestampsPOSIXct[phylogeny$tip.label]
#     distInDaysToRoot <- as.numeric(timestamps - estRootTime) # Difference already in days
#     distsToRoot <- dist.tips.root.phylo(phylogeny, 1:ape::Ntip(phylogeny))
#     perSiteClockRate <- mean(distsToRoot/distInDaysToRoot)
#   }
#   numTips <- ape::Ntip(phylogeny)
#   numNodes <- ape::Nnode(phylogeny)
#   transmissionTree <- phylogeny
#   transmissionTree$tip.label <- lapply(phylogeny$tip.label, function(seqName) {
#     list(name = seqName, time = as.numeric(seqsTimestampsPOSIXct[[seqName]])/86400, region = seqsRegionStamps[[seqName]]) # Time is re-expressed in days.
#   })
#
#   transmissionTree$node.label <- lapply(1:numNodes, function(x) list(time = NULL, region = NA))
#   recursiveFunction <- function(nodeIndex) {
#     childrenIndices <- phangorn::Children(transmissionTree, node = nodeIndex)
#     resolvedChildrenTime <- sapply(childrenIndices, function(childIndex) {
#       resolved <- TRUE
#       if (childIndex > numTips) resolved <- !is.null(transmissionTree$node.label[[childIndex - numTips]]$time) > 0
#       resolved
#     })
#     if (!all(resolvedChildrenTime)) {
#       unresolvedChildrenIndices <- childrenIndices[!resolvedChildrenTime]
#       sapply(unresolvedChildrenIndices, recursiveFunction)
#     }
#     childNodeTimes <- sapply(childrenIndices, FUN = function(x) .getVertexLabel(phylogeny = transmissionTree, vertexNum = x)$time)
#     transBranchLength <- mean(sapply(childrenIndices, function(childIndex) phylogeny$edge.length[match(childIndex, phylogeny$edge[ , 2])]))/perSiteClockRate
#     if (transBranchLength == 0) transBranchLength <- 1
#     topNodeTime <- min(childNodeTimes) - transBranchLength + rnorm(1, mean = 0, sd = 1e-10) # The jittering ensures that coalescence events do not happen simultaneously.
#     transmissionTree$node.label[[nodeIndex - numTips]]$time <<- topNodeTime
#     NULL
#   }
#   startIndex <- numTips + 1
#   recursiveFunction(startIndex) # Function works with side-effects.
#
#   getEdgeListElements <- function(edgeNumber) {
#     nodeNumber <- transmissionTree$edge[edgeNumber, 2]
#     nodeTime <- .getVertexLabel(transmissionTree, vertexNum = nodeNumber)$time
#     # parentNumber <- phangorn::Ancestors(transmissionTree, node = nodeNumber, type = "parent")
#     parentNumber <- transmissionTree$edge[match(nodeNumber, transmissionTree$edge[ , 2]), 1]
#     parentTime <- transmissionTree$node.label[[parentNumber - numTips]]$time
#     list(transmissionTree = nodeTime - parentTime, phylogeny = transmissionTree$edge.length[[edgeNumber]])
#   }
#   transmissionTree$edge.length <- lapply(1:nrow(transmissionTree$edge), FUN = getEdgeListElements)
#   transmissionTree <- .identifyNodeRegions(transmissionTree)
#   transmissionTree
# }

.genStartCoalescenceRates <- function(phyloAndTransTree, estRootTime, control) {
  transmissionTree <- .convertToTransTree(phyloAndTransTree)
  regionNames <- unique(sapply(transmissionTree$tip.label, '[[', "region"))
  getCoalVec <- function(nodeIndex) {
    basicCoalVec <- rep(0, length(regionNames))
    names(basicCoalVec) <- regionNames
    childNodes <- phangorn::Children(transmissionTree, nodeIndex)
    childRegions <- sapply(childNodes, function(x) .getVertexLabel(transmissionTree, x)$region)
    if (length(unique(childRegions)) == 1) basicCoalVec[[childRegions[[1]]]] <- 1
    basicCoalVec
  }
  numWithinCoalescenceEvents <- rowSums(sapply(seq_along(transmissionTree$node.label) + length(phyloAndTransTree$tip.label), getCoalVec))

  startWithin <- numWithinCoalescenceEvents/(sum(transmissionTree$edge.length) * length(regionNames))
  startWithin <- replace(startWithin, which(startWithin == 0), min(startWithin[startWithin > 0])) # Coalescence rates of 0 could create issues down the line.
  names(startWithin) <- regionNames

  funForOptim <- function(x) {
    names(x) <- regionNames # For some reason, the name attribute of the argument x is stripped when the function is called within lbfgs...
    x <- exp(x)
    returnValue <- .topologyLogPriorFunAndGradient(phyloAndTransTree = phyloAndTransTree, Lambda = x, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean, estRootTime = estRootTime)
     # return(list("objective" = -returnValue$objective,
     #            "gradient" = -returnValue$gradient[regionNames]) ) # lbfgs minimises, hence the "-". Mimicking the syntax in the nloptr vignette. FIX THE GRADIENT-BASED ESTIMATION.
    -returnValue$objective
  }

  optimResult <- nloptr::lbfgs(x0 = log(startWithin), fn = funForOptim)
  if (optimResult$convergence < 0) {
    stop("Could not generate starting values for region-specific coalescence rates. Stopping...")
  }
  rescaledResults <- exp(optimResult$par)
  names(rescaledResults) <- regionNames
  rescaledResults
}

# This is a constructor-like function for phylo objects.
phylo <- function(edge, edge.length, tip.label, node.label = NULL) {
  phyloObject <- list(edge = edge, tip.label = tip.label, edge.length = edge.length, Nnode = nrow(edge) - length(tip.label) + 1, node.label = node.label)
  class(phyloObject) <- "phylo"
  phyloObject
}

.optimTreeMCMC <- function(
  startingValues,
  chainId,
  logLikFun,
  DNAbinData,
  evoParsList,
  perSiteClockRate,
  estRootTime,
  covariateFrame,
  control) {
  MCMCcontrol <- control$MCMC.control
  phyDatData <- phangorn::as.phyDat(DNAbinData)
  # The following assumes that the order of tip labels matches the order of sequences in DNAbinData.
  startIterNum <- 1
  MCMCcontainer <- vector(mode = "list", length = floor((MCMCcontrol$n + MCMCcontrol$burnin)/MCMCcontrol$nIterPerSweep))
  filesToRestore <- NULL
  if (!is.null(MCMCcontrol$folderToSaveIntermediateResults) & !is.null(MCMCcontrol$chainId)) {
    filesToRestore <- list.files(path = MCMCcontrol$folderToSaveIntermediateResults, pattern = paste(MCMCcontrol$chainId, "_atIter", sep = ""), full.names = TRUE)
  }
  if (length(filesToRestore) > 0) {
    iterNums <- as.numeric(stringr::str_extract(filesToRestore, "[:digit:]+(?=.Rdata)"))
    filesToRestoreOrder <- order(iterNums)
    MCMCcontainer[1:length(filesToRestore)] <- lapply(filesToRestore[filesToRestoreOrder], function(filename) {
      loadName <- load(filename)
      get(loadName)
    })
    startIterNum <- max(iterNums) + 1
    currentState <- MCMCcontainer[[length(filesToRestore)]]
    cat("Resuming MCMC at iteration ", startIterNum, ". \n", sep = "")
  } else {
    phyloObj <- .convertToPhylo(startingValues$phyloAndTransTree)
    startLogLikValue <- logLikFun(phyloObj, phyDatData, evoParsList)

    currentState <- lapply(1:MCMCcontrol$nChains, function(x) list(paraValues = startingValues, logLik = startLogLikValue)) # All chains start in the same state.
    cat("Launching MCMC... \n")
  }

  logPriorAndTransFunList <- list(
   topology = list(
      # logPriorFun = function(x) {
      #   .topologyLogPriorFun(phyloAndTransTree = x, Lambda = currentState$paraValues$Lambda, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean, estRootTime = estRootTime)
      #   },
      logPriorFun = function(x, Lambda = currentState[[1]]$paraValues$Lambda) {
          .topologyLogPriorFun(phyloAndTransTree = x, Lambda = Lambda, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean, estRootTime = estRootTime)
        },
      transFun = .topologyTransFun),
    b = list( # The mean here corresponds to 7 days for the time between the time the sequence is sampled and the time at which the first transmission linking it to another sequence in the sample occurred.
      logPriorFun =  function(x) .phyloBranchLengthsLogPriorFun(phyloAndTransTree = x, extBranchLengthPriorMean = perSiteClockRate * 7, extBranchLengthPriorSD = perSiteClockRate * 100),
      transFun = function(x) .phyloBranchLengthsTransFun(x, logTransKernSD = MCMCcontrol$phyloBranchLengthsLogTransKernSD, propToModify = MCMCcontrol$propPhyloBranchesToModify)),
    l = list(
      logPriorFun = function(x) .transTreeBranchLengthsConditionalLogPrior(phyloAndTransTree = x, MVratio = MCMCcontrol$transTreeBranchLengthsPriorMVratio),
      transFun = function(x) .transTreeBranchLengthsTransFun(phyloAndTransTree = x, tuningPara = MCMCcontrol$transTreeTuningPara, propToModify = MCMCcontrol$propTransTreeBranchesToModify, rootTransitionPar = MCMCcontrol$rootTransitionPar)),
    Lambda = list(
      logPriorFun = .coalescenceRatesLogPriorFun,
      transFun = function(x) .coalescenceRatesTransFun(x, sd = MCMCcontrol$coalRateKernelSD)))
  if (!control$fixedClockRate) {
    logPriorAndTransFunList$l$logPriorFun <- function(x) .transTreeBranchLengthsConditionalLogPrior(phyloAndTransTree = x, MVratio = MCMCcontrol$transTreeBranchLengthsPriorMVratio)

  logPriorAndTransFunList$logXi <- list(
    logPriorFun = function(x) .logClockRatesLogPrior(phyloAndTransTree = x, logClockRatePriorMean = log(perSiteClockRate), stdDev = control$MCMC.control$logClockRatesPriorSD, strict = control$strictClockModel),
    transFun = function(x) {
      .logClockRatesTransFun(phyloAndTransTree = x, propToModify = MCMCcontrol$propClockRatesToModify, logKernelSD = MCMCcontrol$logClockRatesKernelSD, strict = control$strictClockModel)
    })
  }
  if (length(filesToRestore) == 0) { # New chains...
    logPriorValue <- sapply(names(logPriorAndTransFunList), FUN = function(paraName) logPriorAndTransFunList[[paraName]]$logPriorFun(.getCurrentState(currentStateVector = currentState[[1]], paraName = paraName)))
    names(logPriorValue) <- names(logPriorAndTransFunList)
    currentState <- lapply(1:MCMCcontrol$nChains, function(chainNum) {
      currentStateMember <- currentState[[chainNum]]
      currentStateMember$logPrior <- logPriorValue
      currentStateMember$logPP <- (startLogLikValue + sum(currentStateMember$logPrior)) * (1/MCMCcontrol$temperatureParFun(chainNum))
      currentStateMember
    })
  }
  cat("Chain is called ", chainId, ". Specify this string in MCMC.control if you want to resume simulations.\n", sep = "")
  cat("Starting values for log-priors: \n")
  print(currentState[[1]]$logPrior) # "1" is the cold chain
  cat("Starting value for log-lik.: \n")
  print(currentState[[1]]$logLik) # "1" is the cold chain
  cat("Starting value for log-PP: \n")
  print(currentState[[1]]$logLik + sum(currentState[[1]]$logPrior))
  if (MCMCcontrol$nChains > 1) {
    clusterAddress <- parallel::makeForkCluster(nnodes = MCMCcontrol$nChains)
  }
  totalNumSweeps <- ceiling((MCMCcontrol$n + MCMCcontrol$burnin)/MCMCcontrol$nIterPerSweep)
  sweepFun <- function(sweepNum) {
    chainFun <- function(chainNumber) {
      chainState <- currentState[[chainNumber]]
      # This odd scheme with the default value for Lambda is for the treatment of a putative jump in Lambda in the MCMC. We want to specify Lambda only in that case; else, it corresponds to the value in memory.
      logPriorAndTransFunList$topology <- list(
        logPriorFun = function(x, Lambda = chainState$paraValues$Lambda) {
          .topologyLogPriorFun(phyloAndTransTree = x, Lambda = Lambda, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean, estRootTime = estRootTime)
        },
        transFun = .topologyTransFun)
      # for (MCMCiter in 1:(MCMCcontrol$n + MCMCcontrol$burnin)) {
      for (MCMCiter in 1:MCMCcontrol$nIterPerSweep) {
        # if ((MCMCiter %% MCMCcontrol$print.frequency) == 0) cat("This is MCMC iteration ", MCMCiter, sep = "", ".\n")
        for (paraName in names(logPriorAndTransFunList)) {
          proposalValueAndTransKernRatio <- logPriorAndTransFunList[[paraName]]$transFun(.getCurrentState(chainState, paraName))
          updatedLogPrior <- chainState$logPrior
          updatedLogPrior[[paraName]] <- logPriorAndTransFunList[[paraName]]$logPriorFun(proposalValueAndTransKernRatio$value)
          updatedDualTree <- .getCurrentState(chainState, "topology")
          if (paraName %in% c("topology", "b", "l", "logXi")) {
            updatedDualTree <- proposalValueAndTransKernRatio$value
          }
          updatedLogLik <- chainState$logLik
          # That's where we take care of other conditional priors depending on the parameter whose transition is being processed.
          if (paraName %in% c("topology", "b")) {
            updatedLogLik <- logLikFun(.convertToPhylo(updatedDualTree), phyDatData, evoParsList)
            if (paraName == "b") { # Prior for "l" is conditioned on "b"
              updatedLogPrior[["l"]] <- logPriorAndTransFunList[["l"]]$logPriorFun(updatedDualTree)
            }
          } else if (paraName == "l") {
            updatedLogPrior[["topology"]] <- logPriorAndTransFunList[["topology"]]$logPriorFun(updatedDualTree)
          } else if (paraName == "Lambda") {
            updatedLogPrior[["topology"]] <- logPriorAndTransFunList[["topology"]]$logPriorFun(updatedDualTree, Lambda = proposalValueAndTransKernRatio$value)
          } else if (paraName == "logXi") {
            updatedLogPrior[["l"]] <- logPriorAndTransFunList[["l"]]$logPriorFun(updatedDualTree)
          }
          proposalLogPP <- (updatedLogLik + sum(updatedLogPrior)) * 1/MCMCcontrol$temperatureParFun(chainNumber)
          exponentValue <- proposalLogPP - chainState$logPP
          MHratio <- proposalValueAndTransKernRatio$transKernRatio * exp(exponentValue)
          # cat("Para.:", paraName, "\n", sep = " ")
          # cat("transKernRatio:", proposalValueAndTransKernRatio$transKernRatio, "exponent:", exponentValue, "proposal LogPP:", proposalLogPP, "current LogPP:", chainState$logPP, "\n", sep = " ")
          if (runif(1) <= MHratio) {
            if (paraName %in% c("topology", "l", "b", "logXi")) {
              # cat("Transition accepted. Parameter name:", paraName, "\n New log-PP", proposalLogPP, "\n")
              # cat("Log-priors:\n")
              # print(updatedLogPrior)
              # browser()
              chainState$paraValues$phyloAndTransTree <- updatedDualTree
            } else {
              chainState$paraValues[[paraName]] <- proposalValueAndTransKernRatio$value
            }
            chainState$logLik <- updatedLogLik
            chainState$logPrior <- updatedLogPrior
            chainState$logPP <- proposalLogPP
            chainState$evoParsList <- evoParsList # Needlessly repeated, but does not take too much memory, and allows the chains to be easily resumed from any iteration.
          }
        }
        # cat("End iter. Current log-PP:", chainState$logPP, "\n", sep = " ")
      }
      chainState
    }
    if (MCMCcontrol$nChains > 1) {
      currentState <<- parallel::parLapply(X = 1:MCMCcontrol$nChains, cl = clusterAddress, fun = chainFun)
      # Processing exchanges between chains...
      # Based on https://www.cs.ubc.ca/~nando/540b-2011/projects/8.pdf
      for (k in 1:(MCMCcontrol$nChains - 1)) {
        logExpr <- MCMCcontrol$temperatureParFun(k + 1)/MCMCcontrol$temperatureParFun(k) * currentState[[k + 1]]$logPP + MCMCcontrol$temperatureParFun(k)/MCMCcontrol$temperatureParFun(k + 1) * currentState[[k]]$logPP - (currentState[[k]]$logPP + currentState[[k + 1]]$logPP)
        swapRatio <- exp(logExpr)
        cat("Swap ratio:", swapRatio, "\n")
        if (runif(1) < swapRatio) {
          cat("Swapping chains", k, "and", k + 1, "\n")
          currentState[c(k, k + 1)] <- currentState[c(k + 1, k)]
          currentState[[k]]$logPP <- currentState[[k]]$logPP * MCMCcontrol$temperatureParFun(k + 1)/MCMCcontrol$temperatureParFun(k)
          currentState[[k + 1]]$logPP <- currentState[[k + 1]]$logPP * MCMCcontrol$temperatureParFun(k)/MCMCcontrol$temperatureParFun(k + 1)
        }
      }
    } else {
      currentState <<- lapply(X = 1:MCMCcontrol$nChains, FUN = chainFun)
    }
    cat("Processed MCMC iteration ", sweepNum * MCMCcontrol$nIterPerSweep, ".\n", sep = "")
    cat("\n\n\n Current log-lik.:", currentState[[1]]$logLik)
    cat("Current log-PP:", currentState[[1]]$logPP, "\n\n")
    cat("Current log-priors: \n")
    print(currentState[[1]]$logPrior)
    # if ((MCMCiter %% MCMCcontrol$stepSize) == 0) {
    MCMCcontainer[[sweepNum]] <<- currentState
    if (!is.null(MCMCcontrol$folderToSaveIntermediateResults)) {
      MCMCiter <- sweepNum * MCMCcontrol$nIterPerSweep
      save(currentState, file = paste(MCMCcontrol$folderToSaveIntermediateResults, "/chainID_", chainId, "_atIter", MCMCiter, ".Rdata", sep = ""), compress = TRUE)
    }
    # }
  }
  currentSweep <- (startIterNum - 1)/MCMCcontrol$nIterPerSweep + 1
  lapply(currentSweep:totalNumSweeps, FUN = sweepFun)
  cat("MCMC complete. Finalising... \n")
  lastIterToDrop <- floor(MCMCcontrol$burnin/MCMCcontrol$nIterPerSweep)
  elementsToDrop <- 0
  if (lastIterToDrop > 0) {
    elementsToDrop <- 1:lastIterToDrop
  }

  # elementsToKeep <- seq(
  #   from = floor(MCMCcontrol$burnin/MCMCcontrol$stepSize) + 1,
  #   to = floor((MCMCcontrol$burnin + MCMCcontrol$n)/MCMCcontrol$stepSize))
  lapply(MCMCcontainer[-elementsToDrop], '[[', 1) # We only keep the cold chain.
}

#' Control parameters for the MCMC run
#'
#' Control parameters for the MCMC run.
#'
#' @param n total number of iterations after the burn-in
#' @param stepSize frequency at which iterations are kept in memory, relates to thinning, e.g. stepSize = 100 means that the returned chain will comprise iterations 100, 200, 300, after the burn-in
#' @param burnin number of iterations to discard at the start of the chain
#' @param folderToSaveIntermediateResults folder where intermediate chain results are saved, used to allow chains to be interrupted/resumed
#' @param chainId glyph program previously assigned a chain, used to resume a chain, can be found in the printed output of findBayesianClusters after "Launching MCMC..."
#' @param coalRateKernelSD tuning parameter for the Gaussian coalescence rates transition kernel
#' @param print.frequency number indicating at which interval the program prints information on the chain
#' @param phyloBranchLengthsLogTransKernSD number between 0 and 1, tuning parameter for the phylogenetic branch lengths transition kernel, a lower value translates to larger variations, and a lower acceptance rate in the chain
#' @param propPhyloBranchesToModify number between 0 and 1, proportion of branches that are modified by the transition kernel every iteration
#' @param transTreeTuningPara tuning parameter for the transmission tree transition kernel, determines the boundaries of the uniform density used to propose new branch lengths in the transmission tree; a higher value leads to larger steps, and a lower acceptance rate
#' @param propTransTreeBranchesToModify like propPhyloBranchesToModify, but for the transmission tree
#' @param transTreeBranchLengthsPriorMVratio priors for branch lengths in the transmission tree are gamma-distributed with this shape parameter; at the default value (1), the distribution is exponential
#'
#' @details You should not need to call this function directly: it is merely a convenient and formal way to let users know how the MCMC run can be customised.
#' @return A list of control parameters
#' @examples \dontrun{
#'  MCMC.control(n = 1e6, stepSize = 100, burnin = 1e5)
#' }
#' @export
MCMC.control <- function(n = 2e5, nIterPerSweep = 100, nChains = 1, temperatureParFun = function(x) x, stepSize = 50, burnin = 1e4, folderToSaveIntermediateResults = NULL, chainId = NULL, coalRateKernelSD = 0.5, print.frequency = 10, phyloBranchLengthsLogTransKernSD = log(1.05), propPhyloBranchesToModify = 0.1, transTreeTuningPara = 0.99, propTransTreeBranchesToModify = 0.1, propClockRatesToModify = 0.1, rootTransitionPar = 15, transTreeBranchLengthsPriorMVratio = 1.1, logClockRatesPriorSD = log(1.1), logClockRatesKernelSD = 1) { # transTreeBranchLengthsPriorMVratio must be greater than 1 for the initial optimisation of branch lengths in the transmission tree to work.
  if (temperatureParFun(1) != 1) {
    stop("temperatureParFun(1) must return '1' as it is for the cold chain! Please fix this (or leave the default value) and re-run the code. \n")
  }
  if (nChains == 1) {
    nIterPerSweep <- stepSize
  }
  list(n = n, nIterPerSweep = nIterPerSweep, nChains = nChains, temperatureParFun = temperatureParFun, stepSize = stepSize, burnin = burnin, folderToSaveIntermediateResults = folderToSaveIntermediateResults, chainId = chainId, coalRateKernelSD = coalRateKernelSD, print.frequency = print.frequency, phyloBranchLengthsLogTransKernSD = phyloBranchLengthsLogTransKernSD, propPhyloBranchesToModify = propPhyloBranchesToModify, transTreeTuningPara = transTreeTuningPara, propTransTreeBranchesToModify = propTransTreeBranchesToModify, propClockRatesToModify = propClockRatesToModify, rootTransitionPar = rootTransitionPar, transTreeBranchLengthsPriorMVratio = transTreeBranchLengthsPriorMVratio, logClockRatesKernelSD = logClockRatesKernelSD, logClockRatesPriorSD = logClockRatesPriorSD)
}

.getCurrentState <- function(currentStateVector, paraName = c("topology", "b", "l", "Lambda", "logXi")) {
  paraName <- paraName[[1]]
  if (paraName %in% c("topology", "b", "l", "logXi")) return(currentStateVector$paraValues$phyloAndTransTree)
  currentStateVector$paraValues[[paraName]]
}

.topologyTransFun <- function(phyloAndTransTree) {
  counter <- 0
  repeat {
    counter <- counter + 1
    proposedMove <- phangorn::rNNI(phyloAndTransTree)
    newEdgeLengths <- .deriveTransTreeEdgeLengths(proposedMove)
    if (all(newEdgeLengths >= 0)) break
    if (counter == 50) {
      cat("Could not find a suitable move in the topological space... \n")
      proposedMove <- phyloAndTransTree
      break
    }
  }
  proposedMove$edge.length <- lapply(seq_along(newEdgeLengths), function(edgeIndex) {
    edgeLengthList <- proposedMove$edge.length[[edgeIndex]]
    edgeLengthList$transmissionTree <- newEdgeLengths[[edgeIndex]]
    edgeLengthList
  })

  # A change in the topology can affect node region assignment.
  proposedMove <- .clearNodeRegions(proposedMove)
  proposedMove <- .identifyNodeRegions(proposedMove)

  list(value = proposedMove, transKernRatio = 1)
}

.convertToTransTree <- function(phyloAndTransTree) {
  transmissionTree <- phyloAndTransTree
  if (length(transmissionTree$edge.length[[1]]) > 1) {
    transmissionTree$edge.length <- sapply(phyloAndTransTree$edge.length, '[[', "transmissionTree")
  }
  transmissionTree
}

.convertToPhylo <- function(phyloAndTransTree) {
  phylogeny <- phyloAndTransTree
  phylogeny$edge.length <- sapply(phyloAndTransTree$edge.length, '[[', "phylogeny")
  phylogeny$tip.label <- sapply(phylogeny$tip.label, function(x) x$name) # The other components must be discarded for now, to allow pml to work.
  phylogeny$node.label <- sapply(seq_along(phylogeny$node.label), function(nodeLabelIndex) stringr::str_c(nodeLabelIndex + length(phylogeny$tip.label), phylogeny$node.label[[nodeLabelIndex]], sep = ":", collapse = ":"))
  phylogeny
}

.topologyLogPriorFunAndGradient <- function(phyloAndTransTree, Lambda, numMigrationsPoissonPriorMean = 1, estRootTime = NULL) {
  numTips <- length(phyloAndTransTree$tip.label)
  numNodes <- length(phyloAndTransTree$node.label)
  transmissionTree <- .convertToTransTree(phyloAndTransTree)
  nodeTimes <- sapply(transmissionTree$node.label, function(nodeLabel) nodeLabel$time)
  tipTimes <- sapply(transmissionTree$tip.label, function(tipLabel) tipLabel$time)
  getMigrationTimesAndChildNodeNums <- function(nodeNum) {
    nodeIndicator <- nodeNum > numTips
    indexValue <- nodeNum - nodeIndicator * numTips
    listToIndex <- "tip.label"
    if (nodeIndicator) listToIndex <- "node.label"
    # parentNum <- phangorn::Ancestors(x = transmissionTree, node = nodeNum, type = "parent")
    parentNum <- transmissionTree$edge[match(nodeNum, transmissionTree$edge[ , 2]), 1]
    # parentRegion <- .getVertexLabel(transmissionTree, parentNum)$region
    parentRegion <- transmissionTree$node.label[[parentNum - numTips]]$region
    currentRegion <- transmissionTree[[listToIndex]][[indexValue]]$region
    returnValue <- list(childNode = nodeNum, time = NULL)
    if (parentRegion != currentRegion) returnValue$time <- (transmissionTree$node.label[[parentNum - numTips]]$time + transmissionTree[[listToIndex]][[indexValue]]$time)/2
    returnValue
  }
  migrationTimesAndChildNodeNums <- lapply(seq_along(transmissionTree$node.label), FUN = getMigrationTimesAndChildNodeNums)
  migrationTimes <- do.call("c", lapply(migrationTimesAndChildNodeNums, function(listElement) listElement$time))
  migrationChildNodes <- do.call("c", lapply(migrationTimesAndChildNodeNums, function(listElement) {
    if (!is.null(listElement$time)) return(listElement$childNode)
    NULL
  }))
  verticesInvolved <- c(1:(numTips + numNodes), migrationChildNodes)
  eventType <- rep(c("T", "C", "M"), c(numTips, numNodes, length(migrationTimes)))
  eventTimes <- c(tipTimes, nodeTimes, migrationTimes)
  phyloStructCodeFrame <- data.frame(vertexNum = verticesInvolved, type = eventType, time = eventTimes)
  orderIndices <- order(eventTimes, decreasing = TRUE) # We start at the tips, hence the decreasing time.

  sortedPhyloStructCodeFrame <- phyloStructCodeFrame[orderIndices, ]
  regionNames <- unique(sapply(transmissionTree$tip.label, function(x) x$region))
  numRegions <- length(regionNames)
  numLineagesPerRegion <- rep(0, numRegions)
  names(numLineagesPerRegion) <- regionNames
  cumulativeLogProb <- 0
  # cumulativeGradientWithin <- rep(0, length(Lambda))
  # names(cumulativeGradientWithin) <- names(Lambda)
  timeIndices <- as.numeric(factor(sortedPhyloStructCodeFrame$time))
  for (timeIndex in head(unique(timeIndices), n = -1)) {
    rowsConsidered <- which(timeIndices == timeIndex)
    for (rowIndex in rowsConsidered) {
      vertexNumber <- sortedPhyloStructCodeFrame$vertexNum[rowIndex]
      if (sortedPhyloStructCodeFrame$type[rowIndex] == "T") {
        numLineagesPerRegion[[transmissionTree$tip.label[[vertexNumber]]$region]] <- numLineagesPerRegion[[transmissionTree$tip.label[[vertexNumber]]$region]] + 1
      } else if (sortedPhyloStructCodeFrame$type[rowIndex] == "C") {
        numLineagesPerRegion[[transmissionTree$node.label[[vertexNumber - numTips]]$region]] <- numLineagesPerRegion[[transmissionTree$node.label[[vertexNumber - numTips]]$region]] - 1
      } else {
        childNodeNum <- sortedPhyloStructCodeFrame$vertexNum[rowIndex]
        # parentNodeNum <- phangorn::Ancestors(x = transmissionTree, node = childNodeNum, type = "parent")
        parentNodeNum <- transmissionTree$edge[match(childNodeNum, transmissionTree$edge[ , 2]), 1]
        nodeIndicator <- childNodeNum > numTips
        vertexIndex <- childNodeNum - numTips * nodeIndicator
        listName <- "tip.label"
        if (nodeIndicator) listName <- "node.label"
        childRegion <- transmissionTree[[listName]][[vertexIndex]]$region
        parentRegion <- transmissionTree$node.label[[parentNodeNum - numTips]]$region
        # A migration along a branch in reverse time eliminates one lineage in the child region and adds one in the parent region.
        numLineagesPerRegion[[childRegion]] <- numLineagesPerRegion[[childRegion]] - 1
        numLineagesPerRegion[[parentRegion]] <- numLineagesPerRegion[[parentRegion]] + 1
      }
    }
    intervalDuration <- sortedPhyloStructCodeFrame$time[[rowsConsidered[[1]]]] - sortedPhyloStructCodeFrame$time[[rowsConsidered[[1]] + length(rowsConsidered)]]

    withinRegionCoalescencePossible <- numLineagesPerRegion > 1
    # withinGradient <- rep(0, length(Lambda))
    # names(withinGradient) <- names(Lambda)
    cumulWithinRatePerRegion <- sapply(regionNames[withinRegionCoalescencePossible], function(regionName) {
      choose(numLineagesPerRegion[[regionName]], 2) * Lambda[[regionName]] * intervalDuration
    })
    names(cumulWithinRatePerRegion) <- regionNames[withinRegionCoalescencePossible]
    # withinGradient[names(cumulWithinRatePerRegion)] <- -cumulWithinRatePerRegion/Lambda[names(cumulWithinRatePerRegion)]
    if (length(cumulWithinRatePerRegion) == 0) cumulWithinRatePerRegion <- 0
    totalRate <- sum(cumulWithinRatePerRegion)
    cumulativeLogProb <- cumulativeLogProb - totalRate # That takes into account the exponent of exp(), hence the "-".
    # cumulativeGradientWithin <- cumulativeGradientWithin + withinGradient[names(cumulativeGradientWithin)]

    if (sortedPhyloStructCodeFrame$type[[rowsConsidered[[1]] + length(rowsConsidered)]] == "C") {
      regionCode <- transmissionTree$node.label[[sortedPhyloStructCodeFrame$vertexNum[[rowsConsidered[[1]] + length(rowsConsidered)]] - numTips]]$region
      cumulativeLogProb <- cumulativeLogProb + log(Lambda[[regionCode]] * intervalDuration)
      # cumulativeGradientWithin[[regionCode]] <- cumulativeGradientWithin[[regionCode]] + 1/(Lambda[[regionCode]] * intervalDuration)
    }
  }
  # getMigrationIndicator <- function(vertexIndex) {
  #   if (vertexIndex > numTips) {
  #     currentRegion <- phyloAndTransTree$node.label[[vertexIndex - numTips]]$region
  #   } else {
  #     currentRegion <- phyloAndTransTree$tip.label[[vertexIndex]]$region
  #   }
  #   # parentIndex <- phangorn::Ancestors(phyloAndTransTree, vertexIndex, "parent")
  #   parentIndex <- phyloAndTransTree$edge[match(vertexIndex, phyloAndTransTree$edge[ , 2]), 1]
  #   # parentRegion <- .getVertexLabel(phyloAndTransTree, parentIndex)$region
  #   parentRegion <- phyloAndTransTree$node.label[[parentIndex - numTips]]$region
  #   parentRegion != currentRegion
  # }
  if (!is.null(estRootTime)) {
    cumulativeLogProb <- cumulativeLogProb + .rootTimeLogPriorFun(x = phyloAndTransTree$node.label[[1]]$time, meanValue = estRootTime)
  }

  list(objective = cumulativeLogProb + .numMigrationsLogPriorFunWithMeanPar(x = length(migrationTimes), meanPar = numMigrationsPoissonPriorMean), gradient = NULL)
}

.topologyLogPriorFun <- function(phyloAndTransTree, Lambda, numMigrationsPoissonPriorMean, estRootTime) {
  .topologyLogPriorFunAndGradient(phyloAndTransTree, Lambda, numMigrationsPoissonPriorMean, estRootTime)$objective
}

.getNodeLabel <- function(phylogeny, nodeNum) {
  phylogeny$node.label[[nodeNum - length(phylogeny$tip.label)]]
}

.getTipLabel <- function(phylogeny, tipNum) {
  phylogeny$tip.label[[tipNum]]
}

.numMigrationsLogPriorFunWithMeanPar <- function(x, meanPar = 1) {
  dpois(x = x, lambda = meanPar, log = TRUE)
}

.getVertexLabel <- function(phylogeny, vertexNum) {
  if (vertexNum <= length(phylogeny$tip.label)) {
    return(.getTipLabel(phylogeny, vertexNum))
  }
  return(.getNodeLabel(phylogeny, vertexNum))
}

# .getEvoParsLogPriorsFunList <- function() {
#   K80transversionLogPrior <- function(x) {
#     return(0) # Improper
#   }
#   K80transitionLogPrior <- function(x) {
#     return(0) # Improper
#   }
#   propInvLogPrior <- function(x) {
#     return(0) # Improper
#   }
#   gammaShapeLogPrior <- function(x) {
#     return(0) # Improper
#   }
#   list(K80transitionLogPrior, K80transversionLogPrior, propInvLogPrior, gammaShapeLogPrior)
# }

.getEvoParsTransFunList <- function() {
  K80transversionTransFun <- K80transitionTransFun <- function(currentValue) {
    newPos <- exp(rnorm(n = 1, mean = log(currentValue), sd = 0.01))
    transKernRatio <- dnorm(x = log(currentValue), mean = log(newPos), sd = 0.01)/dnorm(x = log(newPos), mean = log(currentValue), sd = 0.01)
    list(value = newPos,
         transKernRatio = transKernRatio)
  }
  propInvTransFun <- function(x) {

  }
  gammaShapeTransFun <- function(x) {

  }
  list(K80transitionTransFun, K80transversionTransFun, propInvTransFun, gammaShapeTransFun)
}

# .phyloBranchLengthsConditionalLogPriorFun <- function(phyloAndTransTree) {
#   branchLengths <- sapply(phyloAndTransTree$edge.length, '$', "phylogeny")
#   dnorm(x = log(branchLengths), mean = -7, sd = 5, log = TRUE)
# }

.phyloBranchLengthsConditionalLogPriorFun <- function(phyloAndTransTree, perSiteClockRate, numSites) {
  # branchLengths <- sapply(phyloAndTransTree$edge.length, '[[', "phylogeny")
  # transTreeBranchLengths <- sapply(phyloAndTransTree$edge.length, '[[', "transmissionTree")
  # sum(dpois(x = floor(branchLengths * numSites), lambda = perSiteClockRate * transTreeBranchLengths, log = TRUE))
  stop("Prior distribution for transmission tree branch lengths should be conditional on phylogenetic branch lengths, not the other way around. \n")
}

# ML estimate of phylogeny supports external branch lengths of 0, which leads to computational issues with the estimation of the transmission tree.
# A prior for branch lengths which imposes non-zero branch lengths is therefore necessary.
# External branches only have a prior which moves their lengths away from 0.
# Internal branches of length 0 would produce a multifurcation, so should not be found.

.phyloBranchLengthsLogPriorFun <- function(phyloAndTransTree, extBranchLengthPriorMean, extBranchLengthPriorSD) {
    branchLengths <- sapply(phyloAndTransTree$edge.length, "[[", "phylogeny")
    externalEdges <- which(phyloAndTransTree$edge[ , 2] <= ape::Ntip(phyloAndTransTree))
    logNormPars <- getLogNormMMpars(meanValue = extBranchLengthPriorMean, varValue = extBranchLengthPriorSD^2)
    sum(dlnorm(branchLengths[externalEdges], meanlog = logNormPars$mu, sdlog = logNormPars$sigma, log = TRUE))
}

# The mode of the displaced log-normal kernel is set at current branch lengths.
# It is then enough to pick a value for sigma, the variance of the random variable expressed on the log-scale, to get a suitable value for mu.
# A better setup might involve working on the logarithmic scale and simulating normal numbers, making the transition kernel symmetric.

.phyloBranchLengthsTransFun <- function(phyloAndTransTree, logTransKernSD = log(1.05), propToModify = 0.1) {
  branchLengths <- sapply(phyloAndTransTree$edge.length, FUN = '[[', "phylogeny")
  numToModify <- 1
  if (propToModify > 0) {
    numToModify <- ceiling(propToModify * length(branchLengths))
  }
  branchesToModify <- sample.int(n = length(branchLengths), size = numToModify, replace = FALSE)
  previousLengths <- branchLengths[branchesToModify]
  # muPar <- log(previousLengths) + logTransKernSD^2
  # newLengths <- rlnorm(n = length(previousLengths), meanlog = muPar, sdlog = logTransKernSD)
  # logTransKernValue <- dlnorm(x = newLengths, meanlog = muPar, sdlog = logTransKernSD, log = TRUE)
  # muParRev <- log(newLengths) + logTransKernSD^2
  # logTransKernValueRev <- dlnorm(x = previousLengths, meanlog = muParRev, sdlog = logTransKernSD, log = TRUE)
  # transKernRatio <- exp(logTransKernValueRev - logTransKernValue)
  logScaleMean <- log(previousLengths)
  newLengths <- exp(rnorm(n = length(logScaleMean), mean = logScaleMean, sd = logTransKernSD))
  transKernRatio <- 1
  for (i in seq_along(branchesToModify)) {
    phyloAndTransTree$edge.length[[branchesToModify[[i]]]]$phylogeny <- newLengths[[i]]
  }
  list(value = phyloAndTransTree, transKernRatio = transKernRatio)
}

.identifyNodeRegions <- function(transmissionTree) {
  numTips <- ape::Ntip(transmissionTree)
  nodeDepths <- sapply(1:(nrow(transmissionTree$edge) + 1), function(vertexNum) {
    length(phangorn::Ancestors(transmissionTree, vertexNum, "all"))
  })
  # We first perform a bottom-up sweep of the tree.
  treeVertexNums <- seq(1, nrow(transmissionTree$edge) + 1)[order(nodeDepths, decreasing = TRUE)]
  for (vertexNum in head(treeVertexNums, -1)) { # We exclude the root node.
    parentNum <- transmissionTree$edge[match(vertexNum, transmissionTree$edge[ , 2]) , 1]
    currentRegion <- NA
    if (vertexNum <= numTips) {
      currentRegion <- transmissionTree$tip.label[[vertexNum]]$region
    } else {
      currentRegion <- transmissionTree$node.label[[vertexNum - numTips]]$region
    }
    parentRegion <- transmissionTree$node.label[[parentNum - numTips]]$region
    if (is.na(parentRegion[[1]])) {
      transmissionTree$node.label[[parentNum - numTips]]$region <- currentRegion
    } else {
      regionIntersection <- intersect(currentRegion, parentRegion)
      if (length(regionIntersection) > 0) {
        transmissionTree$node.label[[parentNum - numTips]]$region <- regionIntersection
      } else {
        transmissionTree$node.label[[parentNum - numTips]]$region <- c(currentRegion, parentRegion) # The union of regions.
      }
    }
  }

  # We now perform a top-down sweep to resolve ambiguities.
  treeVertexNumsTD <- seq(1, nrow(transmissionTree$edge) + 1)[order(nodeDepths, decreasing = FALSE)]
  for (vertexNum in treeVertexNumsTD) {
    if (vertexNum <= numTips) next # Tips are already resolved
    regionSelected <- currentRegion <- transmissionTree$node.label[[vertexNum - numTips]]$region
    if (length(currentRegion) > 1) {
      regionSelected <- sample(currentRegion, size = 1)
      transmissionTree$node.label[[vertexNum - numTips]]$region <- regionSelected
    }
    nodeChildren <- transmissionTree$edge[match(vertexNum, transmissionTree$edge[ , 1]) , 2]
    nodeChildren <- nodeChildren[nodeChildren > numTips] # Tips are already resolved
    for (child in nodeChildren) { # Is supposed to work even if nodeChildren is empty
      childRegion <- transmissionTree$node.label[[child - numTips]]$region
      if (length(childRegion) > 1) { # Child is unresolved
        regionIntersect <- intersect(regionSelected, childRegion)
        if (length(regionIntersect) > 0) {
          transmissionTree$node.label[[child - numTips]]$region <- regionSelected
        }
      }
    }
  }
  transmissionTree
}

.clearNodeRegions <- function(phyloAndTransTree) {
  phyloAndTransTree$node.label <- lapply(phyloAndTransTree$node.label, function(nodeLabel) {
    nodeLabel$region <- NA
    nodeLabel
  })
  phyloAndTransTree
}

.transTreeBranchLengthsLogPrior <- function(phyloAndTransTree, perSiteClockRate) {
  transTreeSingleBranchLengthLogPrior <- function(transTreeBranchLength) {
    return(0)
  }
  sum(sapply(seq_along(phyloAndTransTree$edge.length), FUN = transTreeSingleBranchLengthLogPrior))
}

.transTreeBranchLengthsConditionalLogPrior <- function(phyloAndTransTree, MVratio = 1) {
  perSiteClockRate <- exp(sapply(phyloAndTransTree$edge.length, "[[", "logXi"))
  phyloBranchLengths <- sapply(phyloAndTransTree$edge.length, "[[", "phylogeny")
  # phyloBranchLengths <- replace(phyloBranchLengths, which(phyloBranchLengths == 0), 1/(2 * numSites)) # Branches of length 0 have already been fixed.
  transTreeBranchLengths <- sapply(phyloAndTransTree$edge.length, "[[", "transmissionTree")
  meanValue <- phyloBranchLengths/perSiteClockRate
  varValue <- meanValue/MVratio
  muPar <- log(meanValue^2/sqrt(varValue + meanValue^2))
  sigmaPar <- sqrt(log(varValue/meanValue^2 + 1))
  sum(dlnorm(x = transTreeBranchLengths, meanlog = muPar, sdlog = sigmaPar, log = TRUE))
  # sum(dnorm(x = log(transTreeBranchLengths), mean = meanValue, sd = log(5), log = TRUE))
}

.transTreeBranchLengthsTransFun <- function(phyloAndTransTree, tuningPara = 0.1, propToModify = 0.1, rootTransitionPar = 15) {
  numTips <- length(phyloAndTransTree$tip.label)
  numNodes <- length(phyloAndTransTree$node.label)
  updateNodeTime <- function(nodeIndex) {
    currentTime <- phyloAndTransTree$node.label[[nodeIndex - numTips]]$time
    childrenIndices <- phangorn::Children(phyloAndTransTree, nodeIndex)
    childrenTimes <- sapply(childrenIndices, function(childIndex) .getVertexLabel(phyloAndTransTree, childIndex)$time)
    branchIndex <- match(nodeIndex, phyloAndTransTree$edge[ , 2])
    if (!is.na(branchIndex)) {
      lowerBound <- currentTime - phyloAndTransTree$edge.length[[branchIndex]]$transmissionTree * tuningPara
    } else {
      lowerBound <- currentTime - runif(1, 0, rootTransitionPar) # Time is in days.
    }
      upperBound <- currentTime + abs(min(childrenTimes) - currentTime) * tuningPara
    updatedTime <- runif(1, lowerBound, upperBound)
    updatedTime
  }
  numNodesToModify <- 1
  if (propToModify > 0) {
    numNodesToModify <- ceiling(propToModify * numNodes)
  }
  nodesToUpdate <- sample(seq_along(phyloAndTransTree$node.label) + numTips, size = numNodesToModify, replace = FALSE)
  updatedNodeTimes <- sapply(nodesToUpdate, updateNodeTime)

  phyloAndTransTree$node.label[nodesToUpdate - numTips] <- lapply(seq_along(nodesToUpdate), function(alongNodeIndex) {
    updatedNode <- phyloAndTransTree$node.label[[nodesToUpdate[alongNodeIndex] - numTips]]
    updatedNode$time <- updatedNodeTimes[[alongNodeIndex]]
    updatedNode
  })
  updatedTransTreeEdges <- .deriveTransTreeEdgeLengths(phyloAndTransTree)
  phyloAndTransTree$edge.length <- lapply(seq_along(phyloAndTransTree$edge.length), function(edgeIndex) {
    updatedEdgeLength <- phyloAndTransTree$edge.length[[edgeIndex]]
    updatedEdgeLength$transmissionTree <- updatedTransTreeEdges[[edgeIndex]]
    updatedEdgeLength
  })

  list(value = phyloAndTransTree, transKernRatio = 1) # The uniform transition kernel is symmetrical.
}

.deriveTransTreeEdgeLengths <- function(phyloAndTransTree) {
  childNodeNumbers <- phyloAndTransTree$edge[ , 2]
  parentNodeNumbers <- phyloAndTransTree$edge[ , 1]
  numTips <- length(phyloAndTransTree$tip.label)
  sapply(childNodeNumbers, function(x) .getVertexLabel(phylogeny = phyloAndTransTree, vertexNum = x)$time) - sapply(parentNodeNumbers, function(x) phyloAndTransTree$node.label[[x - numTips]]$time)
}

.getVertexTimes <- function(phyloAndTransTree, type = c("all", "nodes", "tips")) {
  type <- type[[1]]
  labelList <- phyloAndTransTree$tip.label
  if (type == "all") {
    labelList <- c(phyloAndTransTree$tip.label, phyloAndTransTree$node.label)
  } else if (type == "nodes") {
    labelList <- labelList$node.label
  }
  sapply(labelList, "[[", "time")
}

.coalescenceRatesLogPriorFun <- function(rates) {
  return(0)
}

.coalescenceRatesTransFun <- function(currentState, sd) {
  logNewState <- rnorm(n = length(currentState), mean = log(currentState), sd = sd)
  names(logNewState) <- names(currentState)
  # logProbWeightNewState <- sum(dnorm(logNewState, mean = logCurrentState, sd = sd, log = TRUE))
  # inverseLogProbWeightNewState <- sum(dnorm(logCurrentState, mean = logNewState, sd = sd, log = TRUE))
  list(value = exp(logNewState), transKernRatio = 1) # The Gaussian transition kernel is symmetrical.
}

plotTransmissionTree <- function(dualPhyloAndTransmissionTree, targetRegion = NULL, plotTipLabels = FALSE, plotTitle = NULL, showLabelTime = FALSE, showLabelRegion = FALSE, device = jpeg, filename, argsForDevice = list(), argsForPlotPhylo = NULL, showLegend = TRUE) {
  if (!requireNamespace("ggtree")) {
    stop("Library ggtree is required to plot the transmission tree, but has not been found on this system. Please install it and try again. \n")
  }
  timestampsAsNums <- sapply(dualPhyloAndTransmissionTree$tip.label, "[[", "time")
  timestamps <- as.POSIXct(timestampsAsNums * 86400, origin = "1970-01-01") # From days to seconds, 1970-01-01 is the default origin for POSIXct.
  transmissionTree <- .convertToTransTree(dualPhyloAndTransmissionTree)

  transmissionTree$node.label <- sapply(transmissionTree$node.label, function(x) x$time)
  transmissionTree$tip.label <- sapply(transmissionTree$tip.label, function(x) {
    seqLabels <- x$name
    if (showLabelTime) {
      timeLabels <- stringr::str_c("Time = ", x$time)
      seqLabels <- stringr::str_c(seqLabels, timeLabels, sep = ", ")
    }
    if (showLabelRegion) {
      regionLabels <- stringr::str_c("Region = ", x$region)
      stringr::str_c(seqLabels, regionLabels, sep = ", ")
    }
    seqLabels
  })

  transmissionTree <- treeio::as.treedata(transmissionTree)
  getNodeEdgeColourLinetypesTibble <- function(dualPhyloObj) {
    nodeEdgeColourLabels <- sapply(seq_along(c(dualPhyloObj$tip.label, dualPhyloObj$node.label)), function(nodeNum) {
      if (nodeNum == ape::Ntip(dualPhyloObj) + 1) {
        colourIndic <- dualPhyloObj$node.label[[1]]$region
        if (!is.null(targetRegion) & !identical(targetRegion, colourIndic)) {
          colourIndic <- "Other"
        }
        return(rep(colourIndic, 2))
      }
      parentNode <- phangorn::Ancestors(x = dualPhyloObj, node = nodeNum, type = "parent")
      parentRegion <- .getVertexLabel(dualPhyloObj, parentNode)$region
      currentRegion <- .getVertexLabel(dualPhyloObj, nodeNum)$region
      colourLabel <- nodeColourLabel <- currentRegion
      if (!is.null(targetRegion) & !identical(targetRegion, currentRegion)) {
        colourLabel <- nodeColourLabel <- "Other"
      }
      if (!identical(parentRegion, currentRegion)) {
        colourLabel <- "Transition"
      }
      c(edgeColourLabel = colourLabel, nodeColourLabel = nodeColourLabel)
    })
    linetypes <- rep(0, ncol(nodeEdgeColourLabels))
    linetypes <- replace(linetypes, which(nodeEdgeColourLabels["edgeColourLabel", ] == "Transition"), 1)
    tibble::tibble(edgeRegion = nodeEdgeColourLabels["edgeColourLabel", ], nodeRegion = nodeEdgeColourLabels["nodeColourLabel", ], transition = as.factor(linetypes), node = as.character(seq_along(linetypes)))
  }
  maxSamplingDate <- max(timestamps)
  transmissionTree@phylo$edge.length <- transmissionTree@phylo$edge.length/365 # Time is in days. We need it in years for time-scaled phylogeny.
  transmissionTree@data <- getNodeEdgeColourLinetypesTibble(dualPhyloAndTransmissionTree)
  plottedTree <- do.call(ggtree::ggtree, args = c(list(tr = transmissionTree, mrsd = maxSamplingDate, as.Date = TRUE, mapping = ggplot2::aes(colour = edgeRegion, linetype = transition)), argsForPlotPhylo)) + ggtree::geom_nodepoint(mapping = ggplot2::aes(colour = nodeRegion), shape = 16, size = 2) + ggtree::geom_tippoint(mapping = ggplot2::aes(colour = nodeRegion), shape = 15, size = 2) + ggplot2::scale_color_discrete(name = "Region") + ggplot2::scale_linetype_discrete(name = "Transition", breaks = c(0, 1), labels = c("No", "Yes")) + ggtree::theme_tree2()
  if (plotTipLabels) {
    plottedTree <- plottedTree + ggtree::geom_tiplab()
  }
  if (!is.null(plotTitle)) {
    plottedTree <- plottedTree + ggplot2::ggtitle(plotTitle)
  }
  if (!showLegend) {
    plottedTree <- plottedTree + ggplot2::theme(legend.position = "none")
  }
  do.call(ggtree::ggsave, c(list(filename = filename, plot = plottedTree), argsForDevice))
  NULL
}

#' Returns region-specific clusters
#'
#' Returns region-specific transmission clusters from a tree, with a cluster defined as the largest set of cases registered in a certain region descended from the same introduction of the pathogen.
#'
#' @param phyloAndTransTree `phylo`-like object found in `paraValues$phyloAndTransTree` in the `chain` element of the `findBayesianClusters` output
#' @param clusterRegionCode string indicating for which region cluster estimates should be produced
#'
#' @details Note that clusters are not clades, as they exclude tips descended from the a cluster's MRCA that do not have the right regional label.
#' @return A numeric vector, with a $0$ indicating that the sequence does not belong to the identified region.
#' @examples \dontrun{
#' TO_DO
#' }
#'
#'@export

getRegionClusters <- function(phyloAndTransTree, clusterRegionCode) {
  clusterIndices <- rep(0, ape::Ntip(phyloAndTransTree))
  vertexProcessed <- rep(FALSE, ape::Nedge(phyloAndTransTree) + 1)
  clusterNumber <- 0
  outerFindDescendants <- function(vertexIndex) {
    innerFindDescendants <- function(vertexIndex) {
      vertexProcessed[[vertexIndex]] <<- TRUE
      vertexCountry <- .getVertexLabel(phyloAndTransTree, vertexIndex)$region
      if (vertexCountry == clusterRegionCode) {
        if (vertexIndex <= ape::Ntip(phyloAndTransTree)) {
          clusterIndices[[vertexIndex]] <<- clusterNumber
        } else {
          childrenVertices <- phangorn::Children(phyloAndTransTree, vertexIndex)
          sapply(childrenVertices, innerFindDescendants)
        }
      }
      NULL
    }
    if (!vertexProcessed[[vertexIndex]] & identical(.getVertexLabel(phyloAndTransTree, vertexIndex)$region, clusterRegionCode)) {
      clusterNumber <<- clusterNumber + 1
      innerFindDescendants(vertexIndex)
    }

    if (vertexIndex > ape::Ntip(phyloAndTransTree)) {
      childrenVertices <- phangorn::Children(phyloAndTransTree, vertexIndex)
      sapply(childrenVertices, outerFindDescendants)
    }
    NULL
  }
  outerFindDescendants(ape::Ntip(phyloAndTransTree) + 1)
  clusterIndices
}

.plotTransmissionTree.control <- function(plotTipLabels = TRUE, device = "pdf", filename = filename, argsForDevice = list(), argsForPlotPhylo = NULL) {
  list(plotTipLabels = plotTipLabels, device = device, argsForDevice = argsForDevice, argsForPlotPhylo = argsForPlotPhylo)
}

plotClusters <- function(clusterList, phyloAndTransTree, chainID = NULL, clusterRegionCode, minClusterSize = 5, folderName, plotClusterParent = FALSE, controlListForPlotTransmissionTree = list()) {
  tipNames <- sapply(phyloAndTransTree$tip.label, "[[", "name")
  timestampsAsNums <- sapply(phyloAndTransTree$tip.label, "[[", "time")
  timestamps <- as.POSIXct(timestampsAsNums * 86400, origin = "1970-01-01") # From days to seconds
  clusterIndicesList <- lapply(seq_along(clusterList), function(x) {
    clusterIndexVec <- rep(x, length(clusterList[[x]]))
    names(clusterIndexVec) <- clusterList[[x]]
    clusterIndexVec
  })
  clusterIndices <- do.call("c", clusterIndicesList)
  clusterIndices <- clusterIndices[tipNames]
  clusterSizes <- table(clusterIndices)
  clusterSizes <- clusterSizes[-1] # First element in the table is for unclustered sequences
  keepClusterNumbers <- as.numeric(names(clusterSizes)[clusterSizes >= minClusterSize])
  pruneFunction <- function(treeObj) {
    newTreeObj <- treeObj
    repeat {
      vertexVec <- ape::Ntip(newTreeObj) + 1
      repeat {
        currentVertex <- vertexVec[[1]]
        descendantTipsRegions <- sapply(phangorn::Descendants(newTreeObj, currentVertex, "tips")[[1]], function(tipNumber) newTreeObj$tip.label[[tipNumber]]$region)
        descendantTipsTest <- any(descendantTipsRegions == clusterRegionCode)
        if (!descendantTipsTest) break
        childrenVertices <- NULL
        if (currentVertex > ape::Ntip(newTreeObj)) {
          childrenVertices <- phangorn::Children(newTreeObj, currentVertex)
          childrenVertices <- childrenVertices[childrenVertices > ape::Ntip(newTreeObj)]
        }
        vertexVec <- c(vertexVec[-1], childrenVertices)
        if (length(vertexVec) == 0) break
      }
      if (length(vertexVec) == 0) break
      newTreeObj <- prune.tree(newTreeObj, currentVertex)
    }
    newTreeObj
  }

  getClusterTree <- function(clusterNumber) {
    tipNums <- which(clusterIndices == clusterNumber)
    mrcaNum <- ape::getMRCA(phyloAndTransTree, tipNums)
    if (plotClusterParent) {
      parentNum <- phangorn::Ancestors(phyloAndTransTree, mrcaNum, "parent")
      siblingNum <- phangorn::Siblings(phyloAndTransTree, mrcaNum)
      subTree <- ape::extract.clade(phyloAndTransTree, parentNum)
      subTree <- prune.tree(subTree, siblingNum)
    } else {
      subTree <- ape::extract.clade(phyloAndTransTree, mrcaNum)
    }
    pruneFunction(subTree)
  }
  controlArgsForPlotTransmissionTree <- do.call(.plotTransmissionTree.control, controlListForPlotTransmissionTree)
  clusterTrees <- lapply(keepClusterNumbers, FUN = getClusterTree)
  filenames <- paste(folderName, "/cluster", keepClusterNumbers, "_chainID", chainID, ".", controlArgsForPlotTransmissionTree$device, sep = "")
  plotNames <- paste(clusterRegionCode, ": cluster ", keepClusterNumbers, sep = "")

  plotFunction <- function(treeObj, filename, plotTitle) {
    do.call(plotTransmissionTree, args = c(list(dualPhyloAndTransmissionTree = treeObj, timestamps = timestamps, filename = filename, plotTitle = plotTitle), controlArgsForPlotTransmissionTree, showLegend = FALSE))
    NULL
  }
  lapply(seq_along(clusterTrees), FUN = function(clusterIndex) {
    plotFunction(clusterTrees[[clusterIndex]], filename = filenames[[clusterIndex]], plotTitle = plotNames[[clusterIndex]])
    invisible(NULL)
  })
  invisible(NULL)
}

#' Prunes a phylogeny
#'
#' Prunes a phylo object at the node whose number is given by node
#'
#' @param phylogeny phylo object
#' @param node node number at which the tree should be pruned, must be greater than the number of tips in the tree plus one,
#'
#' @details The function transforms the selected node into a tip whose label will be NA, unless node.label is specified. Note that entering node = ape::Ntip(phylogeny) + 1 will produce an error, as it involves pruning the tree at the root.
#' @return A phylo object
#' @examples \dontrun{
#' set.seed(10)
#' aPhylo <- rtree(10)
#' prunedPhylo <- prune.tree(aPhylo, 12)
#' plot(prunedPhylo)
#' }
#'
#'@export

prune.tree <- function(phylogeny, node) {
  if (node <= ape::Ntip(phylogeny)) return(phylogeny)
  if (node == ape::Ntip(phylogeny) + 1) stop("Can't prune the tree at the root level! \n")

  verticesToRemove <- phangorn::Descendants(x = phylogeny, node = node, type = "all")
  tipsToRemove <- verticesToRemove[verticesToRemove <= ape::Ntip(phylogeny)]
  newNumTips <- ape::Ntip(phylogeny) - length(tipsToRemove) + 1
  nodesToRemove <- c(verticesToRemove[verticesToRemove > ape::Ntip(phylogeny)], node) # node has to be removed too, as it becomes a tip.
  edgesToRemove <- match(verticesToRemove, phylogeny$edge[ , 2])
  nodeLabelToTransform <- .getVertexLabel(phylogeny =  phylogeny, vertexNum = node)
  if (is.null(nodeLabelToTransform)) {
    nodeLabelToTransform <- NA
  }
  untailoredEdgeMat <- phylogeny$edge[-edgesToRemove, ]
  matchVec <- 1:length(unique(as.vector(untailoredEdgeMat)))
  names(matchVec) <- as.character(sort(unique(as.vector(untailoredEdgeMat)))) # The node at which the pruning is done should be a tip. It follows that it should map to newNumTips, and all values greater than that should be incremented by one.
  matchVec[(matchVec >= newNumTips) & (matchVec < matchVec[[as.character(node)]])] <- matchVec[(matchVec >= newNumTips) & (matchVec < matchVec[[as.character(node)]])] + 1
  matchVec[[as.character(node)]] <- newNumTips
  edgeMat <- matrix(matchVec[as.character(untailoredEdgeMat)], ncol = 2)
  tipLabels <- c(phylogeny$tip.label[-tipsToRemove], list(nodeLabelToTransform))
  nodeLabels <- phylogeny$node.label[-(nodesToRemove - ape::Ntip(phylogeny))]
  edgeLengths <- phylogeny$edge.length[-edgesToRemove]
  phylo(edge = edgeMat, tip.label = tipLabels, node.label = nodeLabels, edge.length = edgeLengths)
}

# Time values are not restricted to positive. POSIXct sets January 1st 1970 as time 0.
.rootTimeLogPriorFun <- function(x, meanValue, variance = 1) {
  dnorm(x = x, mean = meanValue, sd = sqrt(variance), log = TRUE)
}

# We pick a log-normal prior on the original scale because we would like the prior on the log-scale to be normal.
.logClockRatesLogPrior <- function(phyloAndTransTree, logClockRatePriorMean, stdDev, strict = FALSE) {
  x <- phyloAndTransTree$edge.length[[1]]$logXi
  if (!strict) {
    x <- sapply(phyloAndTransTree$edge.length, "[[", "logXi")
  }
  sum(dnorm(x = x, mean = logClockRatePriorMean, sd = stdDev, log = TRUE))
}

.logClockRatesTransFun <- function(phyloAndTransTree, propToModify = 0.01, logKernelSD = 1, strict = FALSE) {
  if (!strict) {
    numToModify <- ceiling(propToModify * length(phyloAndTransTree$edge.length))
    posToModify <- sort(sample.int(n = length(phyloAndTransTree$edge.length), size = numToModify, replace = FALSE))

    newEdgeLengths <- lapply(seq_along(phyloAndTransTree$edge.length), function(posIndex) {
      returnValue <- phyloAndTransTree$edge.length[[posIndex]]
      if (posIndex %in% posToModify) {
        returnValue <- phyloAndTransTree$edge.length[[posIndex]]
        currentLogXi <- returnValue$logXi
        returnValue$logXi <- rnorm(n = 1, mean = currentLogXi, sd = logKernelSD)
      }
      returnValue
    })
  } else {
    currentLogXi <- phyloAndTransTree$edge.length[[1]]$logXi
    newLogXi <-  rnorm(n = 1, mean = currentLogXi, sd = logKernelSD)
    newEdgeLengths <- lapply(seq_along(phyloAndTransTree$edge.length), function(posIndex) {
      returnValue <- phyloAndTransTree$edge.length[[posIndex]]
      returnValue$logXi <- newLogXi
      returnValue
    })
  }
  phyloAndTransTree$edge.length <- newEdgeLengths
  list(value = phyloAndTransTree, transKernRatio = 1)
}

getTransTreeEdgeLengths <- function(phyloAndTransTree) {
  sapply(phyloAndTransTree$edge.length, "[[", "transmissionTree")
}

getLogClockRates <- function(phyloAndTransTree) {
  sapply(phyloAndTransTree$edge.length, "[[", "logXi")
}

.getAdjacencyFromVec <- function(clusterIndices) {
  adjMatrixBase <- diag(length(clusterIndices))
  rowIndices <- row(adjMatrixBase)[lower.tri(adjMatrixBase)]
  colIndices <- col(adjMatrixBase)[lower.tri(adjMatrixBase)]

  lowerTriAdj <- lapply(1:(length(clusterIndices) - 1), function(index) {
    currentValue <- clusterIndices[[index]]
    theVec <- clusterIndices[(index + 1):length(clusterIndices)] == currentValue
    theVec
  })
  adjMatrix <- Matrix::sparseMatrix(i = rowIndices,
                            j = colIndices,
                            x = do.call("c",lowerTriAdj),
                            symmetric = TRUE)
  rownames(adjMatrix) <- colnames(adjMatrix) <- names(clusterIndices)
  adjMatrix
}

getClustersFromChains <- function(findBayesianClusterResultsList, linkageThreshold = 0.5, regionLabel, nThreads = 1) {
  if ("chain" %in% names(findBayesianClusterResultsList)) {
    findBayesianClusterResultsList <- list(findBayesianClusterResultsList)
  }
  numElements <- ape::Ntip(findBayesianClusterResultsList[[1]]$chain[[1]]$paraValues$phyloAndTransTree)

  seqNames <- sapply(findBayesianClusterResultsList[[1]]$chain[[1]]$paraValues$phyloAndTransTree$tip.label, "[[", "name")
  regionNames <- sapply(findBayesianClusterResultsList[[1]]$chain[[1]]$paraValues$phyloAndTransTree$tip.label, "[[", "region")
  regionSeqNames <- seqNames[regionNames == regionLabel]

  getCombinedAdjacency <- function(chain) {
    getClusterIndicesInIter <-  function(iterIndex) {
      seqNames <- sapply(chain[[iterIndex]]$paraValues$phyloAndTransTree$tip.label, "[[", "name")
      regionNames <- sapply(chain[[iterIndex]]$paraValues$phyloAndTransTree$tip.label, "[[", "region")
      clusterVec <- getRegionClusters(chain[[iterIndex]]$paraValues$phyloAndTransTree, regionLabel)
      names(clusterVec) <- seqNames
      vecIndicesSubset <- clusterVec[regionSeqNames]
      sortedUniqueValues <- sort(unique(vecIndicesSubset))
      vecIndicesSubsetRenumbered <- match(vecIndicesSubset, sortedUniqueValues)
      names(vecIndicesSubsetRenumbered) <- names(vecIndicesSubset)
      vecIndicesSubsetRenumbered
    }
    clusterIndicesPerIter <- lapply(seq_along(chain), FUN = getClusterIndicesInIter)

    combinedAdjacency <- .getAdjacencyFromVec(clusterIndicesPerIter[[1]])
    if (length(clusterIndicesPerIter) > 1) {
      for (i in 2:length(clusterIndicesPerIter)) {
        combinedAdjacency <- combinedAdjacency + .getAdjacencyFromVec(clusterIndicesPerIter[[i]])
      }
    }
    combinedAdjacency/length(clusterIndicesPerIter)
  }

  if (nThreads == 1) {
    adjacencyFromEachChain <- lapply(findBayesianClusterResultsList, function(result) getCombinedAdjacency(result$chain))
  } else {
    cl <- parallel::makePSOCKcluster(nThreads)
    parallel::clusterEvalQ(cl = cl, expr = library(CovidCluster))
    adjacencyFromEachChain <- parallel::parLapply(cl = cl, X = findBayesianClusterResultsList, fun = function(result) getCombinedAdjacency(result$chain))
    parallel::stopCluster(cl)
  }
  adjMatrixAcrossChains <- Reduce("+", adjacencyFromEachChain)/length(adjacencyFromEachChain)
  adjMatrixAcrossChains@x <- as.numeric(adjMatrixAcrossChains@x >= linkageThreshold)
  roundedAdjMatGraph <- igraph::graph_from_adjacency_matrix(adjMatrixAcrossChains, weighted = NULL, "undirected")
  modules <- igraph::cluster_walktrap(graph = roundedAdjMatGraph)
  modules$membership
}

getDistanceBasedClusters <- function(phyloAndTransTree, distLimit, regionLabel, criterion = c("mrca", "cophenetic", "consecutive")) {
  criterion <- criterion[[1]]
  if (!criterion %in% c("cophenetic", "mrca", "consecutive")) {
    stop("Clustering criterion must be either 'cophenetic', 'mrca', or 'consecutive'.")
  }
  transmissionTree <- .convertToTransTree(phyloAndTransTree)
  if (identical(criterion, "consecutive")) {
    updateTransmissionTree <- function(phylogeny) { # Seems inefficient, but it shouldn't take too much time, even in large samples...
      transmissionTree <- ape::multi2di(phylogeny)
      nodesToFix <- which(sapply(transmissionTree$node.label, identical, "")) + ape::Ntip(transmissionTree)
      while (length(nodesToFix) > 0) {
        currentNode <- phangorn::Ancestors(transmissionTree, nodesToFix[[1]], "parent")
        parentSequence <- nodesToFix[[1]]
        while (identical(.getVertexLabel(phylogeny = transmissionTree, vertexNum = currentNode), "")) {
          parentSequence <- c(parentSequence, currentNode)
          currentNode <- phangorn::Ancestors(transmissionTree, currentNode, "parent")
        }
        vertexToCopy <- .getVertexLabel(transmissionTree, currentNode)

        for (nodeIndex in parentSequence) { # parentSequence cannot include a tip
          transmissionTree$node.label[[nodeIndex - ape::Ntip(transmissionTree)]] <- vertexToCopy
        }
        nodesToFix <- setdiff(nodesToFix, parentSequence)
      }
      transmissionTree
    }
    transmissionTree <- updateTransmissionTree(transmissionTree)
  }

  clusterList <- list()
  nodesToCheck <- ape::Ntip(transmissionTree) + 1
  repeat {
    incrementNodesToCheckFlag <- TRUE
    nodeNumber <- nodesToCheck[[1]]
    # Not very memory efficient, but easier to handle than a list of lists with varying depth.
    if (nodeNumber <= ape::Ntip(transmissionTree)) {
      if (transmissionTree$tip.label[[nodeNumber]]$region == regionLabel) {
        clusterList[[length(clusterList) + 1]] <- transmissionTree$tip.label[[nodeNumber]]$name
      }
      incrementNodesToCheckFlag <- FALSE
    } else {
      if (transmissionTree$node[[nodeNumber - ape::Ntip(transmissionTree)]]$region == regionLabel) {
        allDescendantTips <- phangorn::Descendants(transmissionTree, node = nodeNumber, type = "tips")[[1]]
        allDescendantTipsRegions <- sapply(allDescendantTips, FUN = function(tipNum) transmissionTree$tip.label[[tipNum]]$region)
        descendantTips <- allDescendantTips[allDescendantTipsRegions == regionLabel]
        # descendantTips <- .regionalDescendants(transmissionTree, nodeNumber, FALSE)
        distances <- NULL
        if (identical(criterion, "cophenetic")) {
          distances <- dist.tipPairs.mrca(phylogeny = transmissionTree, tipNumbers = descendantTips)$distance
        } else if (identical(criterion, "mrca")) {
          # distances <- dist.tips.mrca(phylogeny = transmissionTree, tipNumbers = allDescendantTips) # The MRCA should be computed based on all tips, including those not belonging to the target region.
          # distances <- distances[allDescendantTipsRegions == regionLabel]
          if (length(descendantTips) > 1) {
            distances <- dist.tips.mrca(phylogeny = transmissionTree, tipNumbers = descendantTips)
          } else {
            distances <- 0 # We're dealing with a singleton. Should be noted as such.
          }
        } else if (identical(criterion, "consecutive")) {
          # nodeDescendants <- phangorn::Descendants(transmissionTree, nodeNumber, "all")
          # nodeDescendants <- nodeDescendants[nodeDescendants > ape::Ntip(transmissionTree)]
          # nodesToConsiderTest <- sapply(nodeDescendants, function(nodeNum) {
          #   currentRegion <- transmissionTree$node.label[[nodeNum - ape::Ntip(transmissionTree)]]$region
          #   parentNode <- phangorn::Ancestors(transmissionTree, nodeNum, "parent")
          #   parentRegion <- transmissionTree$node.label[[parentNode - ape::Ntip(transmissionTree)]]$region
          #   (currentRegion == regionLabel) & (parentRegion == regionLabel)
          # })
          # if (length(nodesToConsiderTest) > 0) {
          #   distances <- transmissionTree$edge.length[match(nodeDescendants[nodesToConsiderTest], transmissionTree$edge[ , 2])]
          # } else {
          #   distances <- Inf # We've reached a transmission pair. The "consecutive" criterion excludes that these structures could be clusters.
          # }
          cladeSubTree <- ape::extract.clade(transmissionTree, node = nodeNumber)
          tipRegions <- sapply(cladeSubTree$tip.label, "[[", "region")
          tipsOutsideRegion <- which(tipRegions != regionLabel)
          if ((length(tipsOutsideRegion) > 0) & (length(tipsOutsideRegion) < ape::Ntip(cladeSubTree))) {
            cladeSubTree <- ape::drop.tip(cladeSubTree, tip = tipsOutsideRegion)
          } else {
            cladeSubTree <- ape::rtree(1) # This is a placeholder.
          }
          edgesToConsiderIndex <- cladeSubTree$edge[ , 2] > ape::Ntip(cladeSubTree) # Internal branches only
          distances <- cladeSubTree$edge.length[edgesToConsiderIndex]

          if (length(distances) == 0) distances <- Inf # This is a transmission pair: can't be a cluster under the "consecutive" criterion.
        }
        if (all(distances < distLimit)) {
          clusterList[[length(clusterList) + 1]] <- sapply(descendantTips, function(x) transmissionTree$tip.label[[x]]$name)
          incrementNodesToCheckFlag <- FALSE
        }
      }
    }
    if (incrementNodesToCheckFlag) {
      verticesToAdd <- phangorn::Children(transmissionTree, node = nodeNumber)
      nodesToCheck <- c(nodesToCheck, verticesToAdd)
    }

    nodesToCheck <- nodesToCheck[-1]
    if (length(nodesToCheck) == 0) break
  }
  clusterList
}

# Identifies descendants of "node" with a path through the same region.
.regionalDescendants <- function(phyloAndTransTree, node, returnTipNames = FALSE) {
  nodesToCheck <- node
  mrcaRegion <- .getVertexLabel(phyloAndTransTree, node)$region
  descendants <- NULL
  repeat {
    currentNode <- nodesToCheck[[1]]
    cat("Checking node", currentNode, "\n", sep = " ")
    currentRegion <- .getVertexLabel(phyloAndTransTree, currentNode)$region
    if (currentNode <= ape::Ntip(phyloAndTransTree)) {
      if (identical(currentRegion, mrcaRegion)) {
        descendants <- c(descendants, currentNode)
      }
    } else {
      if (identical(currentRegion, mrcaRegion)) {
        nodesToCheck <- c(nodesToCheck, phangorn::Children(phyloAndTransTree, currentNode))
      }
    }
    nodesToCheck <- nodesToCheck[-1]
    if (length(nodesToCheck) == 0) break
  }
  descendants <- sort(descendants)

  if (returnTipNames) {
    descendants <- sapply(descendants, function(tipNumber) phyloAndTransTree$tip.label[[tipNumber]]$name)
  }
  # if (length(descendants) == 0) stop(cat("Node", node, "has no descendant in region", mrcaRegion,"\n Please check the provided tree. \n", sep = " "))
  descendants
}

dist.node.ancestor <- function(phylogeny, node1, node2) {
  ancestorsNode1 <- phangorn::Ancestors(phylogeny, node1)
  ancestorsNode2 <- phangorn::Ancestors(phylogeny, node2)
  node1ancestorNode2 <- node1 %in% ancestorsNode2
  node2ancestorNode1 <- node2 %in% ancestorsNode1
  if (!(node1ancestorNode2 | node2ancestorNode1)) {
    stop("node1 and node2 are not on the same lineage (or are invalid numbers). \n")
  }
  if (node1ancestorNode2) {
    node2copy <- node2
    node2 <- node1
    node1 <- node2copy
    ancestorsNode1 <- ancestorsNode2
  }
  ancestorPos <- match(node2, ancestorsNode1)
  nodesToProcess <- c(node1, ancestorsNode1[ancestorPos - 1])
  branchLengths <- phylogeny$edge.length[match(nodesToProcess, phylogeny$edge[ , 2])]
  sum(branchLengths)
}

dist.tipPairs.mrca <- function(phylogeny, tipNumbers) {
  allCombn <- combn(seq_along(tipNumbers), 2)
  distanceByPair <- apply(allCombn, 2, function(indexPair) {
    sum(dist.tips.mrca(phylogeny, indexPair))
  })
  data.frame(node1 = allCombn[1, ], node2 = allCombn[2, ], distance = distanceByPair)
}

dist.tips.mrca <- function(phylogeny, tipNumbers) {
  tipsMRCA <- ape::getMRCA(phylogeny, tipNumbers)
  sapply(tipNumbers, function(tipNum) dist.node.ancestor(phylogeny, tipNum, tipsMRCA))
}

dist.tips.root.phylo <- function(phylogeny, tips) {
  if (is.character(tips)) {
    tips <- match(tips, phylogeny$tip.label)
  }
  ancestorsList <- phangorn::Ancestors(phylogeny, tips)
  if (length(tips) == 1) {
    ancestorsList <- list(ancestorsList)
  }
  ancestorsList <- lapply(seq_along(ancestorsList), function(index) {
    if (index == ape::Ntip(phylogeny) + 1) return(NA)
    c(index, ancestorsList[[index]])
  })
  edgeMatchVec <- match(1:(nrow(phylogeny$edge) + 1), phylogeny$edge[ , 2])

  sapply(ancestorsList, function(nodeAncestors) {
    sum(phylogeny$edge.length[edgeMatchVec[head(nodeAncestors, -1)]]) # Last element is the root, which does not have a supporting edge
  })
}

.genStartPhyloAndTransTreeAlt <- function(phylogeny, timestampsInDays, regionStamps, logClockRatePriorMean, estRootTime, control) {
  MCMCcontrol <- do.call(MCMC.control, control$MCMC.control)
  MVratio <- MCMCcontrol$transTreeBranchLengthsPriorMVratio
  phyloAndTransTree <- phylogeny
  phyloAndTransTree$edge.length <- replace(phyloAndTransTree$edge.length, which(phyloAndTransTree$edge.length == 0), control$lengthForNullExtBranchesInPhylo)
  phyloAndTransTree$tip.label <- lapply(phyloAndTransTree$tip.label, function(x) list(name = x, time = timestampsInDays[[x]], region = regionStamps[[x]]))
  phyloAndTransTree$node.label <- lapply(1:ape::Nnode(phyloAndTransTree), function(node) {
    list(time = NA, region = NA)
  })
  phyloAndTransTree <- .identifyNodeRegions(phyloAndTransTree)
  # Initialises transmission tree branch lengths conditional on phylogenetic branch lengths and an equal clock rate across branches
  phyloAndTransTree <- .initialiseTree(phyloAndTransTree, timestampsInDays, regionStamps, startLogClockRate = logClockRatePriorMean, estRootTime, control)

  # Improves transmission tree branch lengths conditional on phylogenetic branch lengths and an equal clock rate across branches
  currentScore <- (.transTreeBranchLengthsConditionalLogPrior(phyloAndTransTree, MVratio = MCMCcontrol$transTreeBranchLengthsPriorMVratio) + .rootTimeLogPriorFun(phyloAndTransTree$node.label[[1]]$time, meanValue = estRootTime) + .logClockRatesLogPrior(phyloAndTransTree, logClockRatePriorMean = logClockRatePriorMean, stdDev = control$MCMC.control$logClockRatesPriorSD, strict = control$strictClockModel))
  phyloAndTransTree <- .refineTimeEstimates(phyloAndTransTree = phyloAndTransTree, estRootTime = estRootTime, logClockRatePriorMean = logClockRatePriorMean, control = control)
  phyloAndTransTree
}

.initialiseTree <- function(phyloAndTransTree, timestampsInDays, regionStamps, startLogClockRate, estRootTime, control) {
  MVratio <- do.call(MCMC.control, control$MCMC.control)$transTreeBranchLengthsPriorMVratio

  resolvedNodeIndices <- sapply((ape::Ntip(phyloAndTransTree) + 1):(nrow(phyloAndTransTree$edge) + 1), function(nodeNumber) {
    childrenNodes <- phangorn::Children(phyloAndTransTree, nodeNumber)
    all(childrenNodes <= ape::Ntip(phyloAndTransTree))
  })
  nodesToCheck <- which(resolvedNodeIndices) + ape::Ntip(phyloAndTransTree)
  repeat {
    nodeNumber <- nodesToCheck[[1]]
    childrenNums <- phangorn::Children(phyloAndTransTree, nodeNumber)
    matchPhyloBranches <- match(childrenNums, phyloAndTransTree$edge[ , 2])
    childrenLengths <- phyloAndTransTree$edge.length[matchPhyloBranches]
    childrenMeans <- childrenLengths/exp(startLogClockRate)
    varValue <- childrenMeans/MVratio
    muPar <- log(childrenMeans^2/sqrt(varValue + childrenMeans^2))
    sigmaPar <- sqrt(log(varValue/childrenMeans^2 + 1))

    childrenTimes <- sapply(childrenNums, function(childNum) .getVertexLabel(phyloAndTransTree, childNum)$time)
    maximumTime <- min(childrenTimes)
    funToOptimise <- function(rootTime) {
      -sum(dlnorm(x = childrenTimes - rootTime, meanlog = muPar, sdlog = sigmaPar, log = TRUE))
    }
    gradFct <- function(rootTime) {
      -sum(1/(childrenTimes - rootTime) * (1 + 1/sigmaPar^2 * (log(childrenTimes - rootTime) - muPar)))
    }
    optRootTime <- nloptr::lbfgs(x0 = maximumTime - 1, fn = funToOptimise, gr = gradFct, upper = maximumTime, control = list(xtol_rel = 1e-4))
    # optRootTime <- nloptr::lbfgs(x0 = maximumTime - 1, fn = funToOptimise, upper = maximumTime)
    phyloAndTransTree$node.label[[nodeNumber - ape::Ntip(phyloAndTransTree)]]$time <- optRootTime$par + rnorm(1, sd = 1e-15) # A bit of jittering to ensure that we don't get equal branching times...

    # Updating nodesToCheck (if possible)
    # The condition is necessary since the root node has no sibling.
    if (nodeNumber > ape::Ntip(phyloAndTransTree) + 1) {
      nodeSiblings <- phangorn::Siblings(phyloAndTransTree, nodeNumber, include.self = FALSE)
      siblingTimes <- sapply(nodeSiblings, function(nodeNum) .getVertexLabel(phyloAndTransTree, nodeNum)$time)
      if (all(!is.na(siblingTimes))) {
        vertexParent <- phangorn::Ancestors(phyloAndTransTree, nodeNumber, type = "parent")
        nodesToCheck <- c(nodesToCheck, vertexParent)
      }
    }
    nodesToCheck <- nodesToCheck[-1]
    if (length(nodesToCheck) == 0) break
  }
  phyloAndTransTree$edge.length <- lapply(seq_along(phyloAndTransTree$edge.length), function(x) {
    childNum <- phyloAndTransTree$edge[x, 2]
    parentNum <- phyloAndTransTree$edge[x, 1]
    transTreeBranchLength <- .getVertexLabel(phyloAndTransTree, childNum)$time - .getVertexLabel(phyloAndTransTree, parentNum)$time
    startLogXi <- log(phyloAndTransTree$edge.length[[x]])  - log(transTreeBranchLength)
    list(phylogeny = phyloAndTransTree$edge.length[[x]], transmissionTree = transTreeBranchLength, logXi = startLogXi)
  })
  phyloAndTransTree
}

.refineTimeEstimates <- function(phyloAndTransTree, estRootTime, logClockRatePriorMean, control) {
  cat("Refining starting value for transmission tree...\n")
  MCMCcontrol <- do.call("MCMC.control", args = control$MCMC.control)
  MVratio <- MCMCcontrol$transTreeBranchLengthsPriorMVratio
  nodeChildrenParentParsList <- lapply(seq_along(phyloAndTransTree$node.label), function(nodeIndex) {
    nodeLabel <- nodeIndex + ape::Ntip(phyloAndTransTree)
    parentNum <- NA
    supportBranchLength <- NA
    supportEdgeNum <- NA
    if (nodeLabel > ape::Ntip(phyloAndTransTree) + 1) {
      parentNum <- phangorn::Ancestors(phyloAndTransTree, nodeLabel, "parent")
      supportEdgeNum <- match(nodeLabel, phyloAndTransTree$edge[ , 2])
      supportBranchLength <- phyloAndTransTree$edge.length[[supportEdgeNum]]$phylogeny
    }

    childrenNums <- phangorn::Children(phyloAndTransTree, nodeLabel)
    childrenBranchLengths <- sapply(childrenNums, function(vertexNum) {
      matchValue <- match(vertexNum,  phyloAndTransTree$edge[ , 2])
      phyloAndTransTree$edge.length[[matchValue]]$phylogeny
    })
    childrenEdgeNums <- match(childrenNums, phyloAndTransTree$edge[ , 2])
    list(parent = parentNum, supportEdgeNum = supportEdgeNum, children = childrenNums, supportBranchLength = supportBranchLength, childrenBranchLengths = childrenBranchLengths, childrenEdgeNums = childrenEdgeNums)
  })

  for (nodeIndex in seq_along(phyloAndTransTree$node.label)) {
    nodeInfoList <- nodeChildrenParentParsList[[nodeIndex]]
    childrenNodes <- nodeInfoList$children
    childrenTimes <- sapply(childrenNodes, function(childIndex) .getVertexLabel(phylogeny = phyloAndTransTree, vertexNum = childIndex)$time)
    parentTime <- -Inf
    if (nodeIndex != 1) {
      parentTime <- .getVertexLabel(phylogeny = phyloAndTransTree, vertexNum = nodeInfoList$parent)$time
    }
    objectiveFunForOptim <- function(timePointAndLogXiValues, phyloAndTransTree) {
      timePoint <- timePointAndLogXiValues[[1]]
      phyloAndTransTree$node.label[[nodeIndex]]$time <- timePoint
      logXiValues <- tail(timePointAndLogXiValues, n = -1)
      childrenLogXiValues <- head(logXiValues, n = length(childrenTimes))
      for (index in seq_along(childrenLogXiValues)) {
        childEdgeNum <- nodeInfoList$childrenEdgeNums[[index]]
        phyloAndTransTree$edge.length[[childEdgeNum]]$logXi <- childrenLogXiValues[[index]]
      }

      if (nodeIndex > 1) {
        supportLogXi <- tail(logXiValues, n = 1)
        supportEdgeNum <- nodeInfoList$supportEdgeNum
        phyloAndTransTree$edge.length[[supportEdgeNum]]$logXi <- supportLogXi
      }
      objective <- -1 * (.transTreeBranchLengthsConditionalLogPrior(phyloAndTransTree = phyloAndTransTree, MVratio = MVratio) + .logClockRatesLogPrior(phyloAndTransTree = phyloAndTransTree, logClockRatePriorMean = logClockRatePriorMean, stdDev = control$MCMC.control$logClockRatesPriorSD, strict = control$strictClockModel) + .rootTimeLogPriorFun(phyloAndTransTree$node.label[[1]]$time, meanValue = estRootTime))
      objective
    }
    currentTime <- .getVertexLabel(phyloAndTransTree, nodeIndex + ape::Ntip(phyloAndTransTree))$time
    currentLogXis <- sapply(nodeInfoList$children, function(childNum) {
      branchNum <- match(childNum, phyloAndTransTree$edge[ , 2])
      phyloAndTransTree$edge.length[[branchNum]]$logXi
    })
    if (nodeIndex > 1) {
      parentNum <- nodeInfoList$parent
      parentBranchNum <- match(parentNum, phyloAndTransTree$edge[ , 2])
      currentLogXis <- c(currentLogXis, phyloAndTransTree$edge.length[[parentBranchNum]]$logXi)
    }
    startVec <- c(currentTime, currentLogXis)

    optTimepointAndLogXi <- nloptr::lbfgs(x0 = startVec, fn = objectiveFunForOptim, lower = c(parentTime, rep(-23, length(currentLogXis))), upper = c(min(childrenTimes), rep(0, length(currentLogXis))), phyloAndTransTree = phyloAndTransTree, control = list(xtol_rel = 1e-4))

    phyloAndTransTree$node.label[[nodeIndex]]$time <- optTimepointAndLogXi$par[[1]]
    for (index in seq_along(nodeInfoList$children)) {
      phyloAndTransTree$edge.length[[nodeInfoList$childrenEdgeNums[[index]]]]$logXi <- optTimepointAndLogXi$par[[index + 1]] # First element of par is a time point, hence the + 1.
    }
    if (nodeIndex > 1) {
      phyloAndTransTree$edge.length[[nodeInfoList$supportEdgeNum]]$logXi <- tail(optTimepointAndLogXi$par, n = 1)
    }
  }
  cat("Done!")
  phyloAndTransTree
}

getLogNormMMpars <- function(meanValue, varValue) {
  meanValueSq <- meanValue^2
  list(mu = log(meanValueSq/sqrt(varValue + meanValueSq)), sigma = sqrt(log(varValue/meanValueSq + 1)))
}

.optimStrictOrRelaxedXi <- function(phyloAndTransTree, logClockRatePriorMean, control) {
  funForOptimiser <- function(logXiVec) {
    for (i in seq_along(phyloAndTransTree$edge.length)) {
      phyloAndTransTree$edge.length[[i]]$logXi <- logXiVec[[i]]
    }
    -1 * (.transTreeBranchLengthsConditionalLogPrior(phyloAndTransTree, MVratio = control$MCMC.control$transTreeBranchLengthsPriorMVratio) + .logClockRatesLogPrior(phyloAndTransTree, logClockRatePriorMean = logClockRatePriorMean, stdDev = control$MCMC.control$logClockRatesPriorSD, strict = control$strictClockModel))
  }
  startValue <- sapply(phyloAndTransTree$edge.length, "[[", "logXi")
  optimResults <- nloptr::lbfgs(x0 = startValue, fn = funForOptimiser)
  for (i in seq_along(phyloAndTransTree$edge.length)) {
    phyloAndTransTree$edge.length[[i]]$logXi <- optimResults$par[[i]]
  }
  phyloAndTransTree
}


