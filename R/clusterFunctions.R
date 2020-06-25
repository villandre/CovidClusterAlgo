#' Clusters SARS-Cov-2 sequencing data
#'
#' The function uses MCMC to produce cluster estimates based on sequencing data and covariate information. The function automatically derives suitable starting values for all parameters.
#'
#' @param DNAbinData DNAbin object, the sequencing data,
#' @param covariateFrame A data.frame containing covariate values in the order of elements in DNAbinData.
#' @param branchLengthPriorsFun A function with one argument, a scalar indicating the actual branch length, specifying the prior distribution for the branch lengths.
#' @param evoParsPriorsFun A function with one argument, a numeric vector with named elements if priors vary among *log-evolutionary parameters*, can be left blank if fixEvoPars is set to TRUE in 'control'.
#' @param clusterScoringFun A function with four arguments (in order): cluster labels (numeric), phylogeny (phylo) genotyping data (DNAbin), covariate information (data.frame).
#' @param control List of control parameters, cf. findBayesianClusters.control
#'
#' @details Spatial coordinates must use the longitude/latitude ("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") or sinusoidal ("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") projection. In the latter case, they will be automatically converted to the lon./lat. projection.
#'
#' Some of the control parameters should be tuned to ensure better computational or predictive performance, or to make it possible to stop and resume model fitting. See INLAMRA.control.
#'
#'
#' @return A list with three components:
#' \itemize{
#'  \item{hyperMarginalMoments} {A data.frame giving the mean, and standard deviation of the marginal log-hyperparameter posteriors, as well as their 95\% credible intervals. Note that the scaling of the time hyperparameters depends on the provided time values. If they are inputted as POSIX* objects, time hyperparameters will relate to time measured in days. If they are inputted as numeric, the original scale is used instead.}
#'  \item{FEmarginalMoments} {A data.frame giving the mean and standard deviation of the marginal fixed effects posteriors, as well as their 95\% credibility intervals.}
#'  \item{predMoments} {A data.frame with two columns, Mean and SD. The order of the predictions matches the one in predCovariateFrame.}
#' }
#'
#' @examples
#' \dontrun{
#' # See example in vignette.
#' }
#' @export

findBayesianClusters <- function(
  DNAbinData,
  covariateFrame = NULL,
  seqsTimestampsPOSIXct,
  seqsRegionStamps = rep(1, ifelse(is.matrix(DNAbinData), nrow(DNAbinData), length(DNAbinData))),
  mutationRate,
  startingValuePhylo = NULL,
  rootTimeLogPriorFunction = NULL,
  evoParsList = NULL,
  clusterScoringFun = NULL,
  control = list()) {
  # .performBasicChecks()

  control <- do.call('.defineBayesianClustersControl', control)

  if (is.null(startingValuePhylo)) {
    MLphyloAndEvoPars <- .genMLphyloAndEvoPars(DNAbinData)
    startingValuePhylo <- MLphyloAndEvoPars$phylogeny
    if (is.null(evoParsList)) evoParsList <- MLphyloAndEvoPars$evoParsList
  }

  startingValues <- list()

  # The "edge.length" component of "phylogeny" is a list containing branch lengths for both the transmission tree and phylogeny.
  #
  startingValues$phyloAndTransTree <- .genStartDualPhyloAndTransTree(phylogeny = startingValuePhylo, seqsTimestampsPOSIXct = seqsTimestampsPOSIXct, seqsRegionStamps = seqsRegionStamps, mutationRate = mutationRate, control = control$controlForGenStartTransmissionTree)

  startingValues$Lambda <- .genStartCoalescenceRates(startingValues$transmissionTree)

  sampledTreesWithPP <- .optimTreeMCMC(
    startingValues = startingValues,
    logLikFun = control$logLikFun,
    DNAbinData = DNAbinData,
    evoParsList = evoParsList,
    covariateFrame = covariateFrame,
    control = control$MCMC.control)
  .clusterFun(MCMCoutput = sampledTreesWithPP, clusterScoringFun = control$clusterScoringFun)
  sampledTreesWithPP
}

.performBasicChecks <- function() {

}

.defineBayesianClustersControl <- function(
  logLikFun = presetPML,
  MCMC.control = MCMC.control(),
  controlForGenStartTransmissionTree = .controlForGenStartTransmissionTree()) {
  list(logLikFun = logLikFun, MCMC.control = do.call("MCMC.control", MCMC.control), controlForGenStartTransmissionTree = do.call(".controlForGenStartTransmissionTree", controlForGenStartTransmissionTree))
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

.genMLphyloAndEvoPars <- function(DNAbinData) {
  distMatrix <- ape::dist.dna(x = DNAbinData, gamma = TRUE, pairwise.deletion = TRUE)
  startingPhylo <- ape::bionj(distMatrix)
  startPML <- phangorn::pml(tree = startingPhylo, data = phangorn::as.phyDat(DNAbinData), k = 4, model = "GTR")
  optimisedPhylo <- phangorn::optim.pml(object = startPML, optNni = TRUE, optBf = TRUE, optQ = TRUE, optInv = TRUE, optGamma = TRUE, optEdge = TRUE, optRate = TRUE)
  list(phylogeny = optimisedPhylo$tree, evoParsList = list(Q = optimisedPhylo$Q, bf = optimisedPhylo$bf, gammaShape = optimisedPhylo$gammaShape, propInv = optimisedPhylo$propInv, ncat = 4))
}

.controlForGenStartTransmissionTree <- function(startTransTreeBranchLength = 1/24) {
  list(startTransTreeBranchLength = startTransTreeBranchLength)
}

.genStartDualPhyloAndTransTree <- function(phylogeny, seqsTimestampsPOSIXct, seqsRegionStamps, mutationRate, control) {
  numTips <- ape::Ntip(phylogeny)
  numNodes <- ape::Nnode(phylogeny)
  transmissionTree <- phylogeny
  transmissionTree$tip.label <- lapply(phylogeny$tip.label, function(seqName) {
    list(name = seqName, time = as.numeric(seqsTimestampsPOSIXct[[seqName]])/86400, region = seqsRegionStamps[[seqName]]) # Time is re-expressed in days.
  })

  transmissionTree$node.label <- lapply(1:numNodes, function(x) list(time = NULL, region = NULL))
  recursiveFunction <- function(nodeIndex) {
    childrenIndices <- phangorn::Children(transmissionTree, node = nodeIndex)
    resolvedChildrenTime <- sapply(childrenIndices, function(childIndex) {
      resolved <- TRUE
      if (childIndex > ape::Ntip(transmissionTree)) resolved <- !is.null(transmissionTree$node.label[[childIndex - numTips]]$time) > 0
      resolved
    })
    if (!all(resolvedChildrenTime)) {
      unresolvedChildrenIndices <- childrenIndices[!resolvedChildrenTime]
      sapply(unresolvedChildrenIndices, recursiveFunction)
    }
    childNodeTimes <- sapply(childrenIndices, function(x) {
      if (x <= ape::Ntip(transmissionTree)) {
        return(transmissionTree$tip.label[[x]]$time)
      }
      return(transmissionTree$node.label[[x - ape::Ntip(transmissionTree)]]$time)
    })
    topNodeTime <- min(childNodeTimes) - control$startTransTreeBranchLength # One hour, as time is expressed in days.
    transmissionTree$node.label[[nodeIndex - numTips]]$time <<- topNodeTime
    NULL
  }
  startIndex <- numTips + 1
  recursiveFunction(startIndex) # Function works with side-effects.

  getEdgeListElements <- function(edgeNumber) {
    nodeNumber <- transmissionTree$edge[edgeNumber, 2]
    nodeTime <- ifelse(nodeNumber <= ape::Ntip(transmissionTree), transmissionTree$tip.label[[nodeNumber]]$time, transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$time)
    parentNumber <- phangorn::Ancestors(transmissionTree, node = nodeNumber, type = "parent")
    parentTime <- transmissionTree$node.label[[parentNumber - ape::Ntip(transmissionTree)]]$time
    list(transmissionTree = nodeTime - parentTime, phylogeny = transmissionTree$edge.length[[edgeNumber]])
  }
  transmissionTree$edge.length <- lapply(1:ape::Nedge(transmissionTree), FUN = getEdgeListElements)
  transmissionTree <- .identifyNodeRegions(transmissionTree)
  transmissionTree
}

.genStartCoalescenceRates <- function(dualPhyloAndTransTree) {
  transmissionTree <- .convertToTransTree(dualPhyloAndTransTree)
  tipTimes <- sapply(transmissionTree$tip.label, function(x) x$time)
  rootTime <- transmissionTree$node.label[[1]]$time
  maxTimeRange <- max(tipTimes - rootTime)/86400 # Rate is in days.
  numWithinCoalescenceEvents <- sum(sapply(transmissionTree$node.label, function(x) length(stringr::str_split(x, pattern = ",")) == 1))
  numBetweenCoalescenceEvents <- length(transmissionTree$node.label) - numWithinCoalescenceEvents
  regionNames <- unique(sapply(transmissionTree$tip.label, '$', "region"))
  numRegions <- length(regionNames)
  startWithin <- rep(maxTimeRange/numWithinCoalescenceEvents, numRegions)
  names(startWithin) <- regionNames
  startBetween <- c(between = maxTimeRange/numBetweenCoalescenceEvents)
  startValueForOptim <- c(startWithin, startBetween)
  funForOptim <- function(x) {
    LambdaList <- list(within = head(x, n = -1), between = tail(x, n = 1))
    .topologyLogPriorFun(dualPhyloAndTransTree = dualPhyloAndTransTree, Lambda = LambdaList)
  }
  optimResult <- nloptr::lbfgs(x0 = startValueForOptim, fn = funForOptim)
  list(within = head(optimResult$par, n = -1), between = tail(optimResult$par, n = 1))
}

# This is a constructor-like function for phylo objects.
phylo <- function(edge, edge.length, tip.label, node.label = NULL) {
  phyloObject <- list(edge = edge, tip.label = tip.label, edge.length = edge.length, Nnode = nrow(edge) - length(tip.label) + 1, node.label = node.label)
  class(phyloObject) <- "phylo"
  phyloObject
}

.optimTreeMCMC <- function(
  startingValues,
  logLikFun,
  DNAbinData,
  evoParsList,
  mutationRate,
  covariateFrame,
  control) {
  control <- do.call(MCMC.control, control)
  phyDatData <- phangorn::as.phyDat(DNAbinData)
  # The following assumes that the order of tip labels matches the order of sequences in DNAbinData.
  startLogLikValue <- logLikFun(startingValues$phylogeny, phyDatData, evoParsList)

  MCMCcontainer <- vector(mode = list(), length = control$n + control$burnin)

  currentState <- list(paraValues = startingValues, logLik = startLogLikValue)

  # Functions are defined here with an external dependence on currentState, but this is voluntary.
  logPriorAndTransFunList <- list(
    topology = list(
      logPriorFun = function(x) .topologyLogPriorFun(dualPhyloAndTransTree = x, Lambda = currentState$paraValues$Lambda),
      transFun = .topologyTransFun),
    b = list(
      logPriorFun = function(x) .logPhyloBranchLengthsLogPriorFun(dualPhyloAndTransTree = x, mutationRate = mutationRate, numSites = ncol(DNAbinData)),
      transFun = .logPhyloBranchLengthsTransFun),
    l = list(
      logPriorFun = function(x) .transTreeBranchLengthsLogPrior(dualPhyloAndTransTree = x, mutationRate = mutationRate),
      transFun = .transTreeBranchLengthsTransFun),
    LambdaWithin = list(
      logPriorFun = .coalescenceRatesWithinLogPriorFun,
      transFun = function(x) .coalescenceRatesWithinTransFun(x, sd = control$withinCoalRateKernelSD)),
    LambdaBetween = list(
      logPriorFun = .coalescenceRatesWithinLogPriorFun,
      transFun = function(x) .coalescenceRatesBetweenTransFun(x, sd = control$betweenCoalRateKernelSD)))
  currentState$logPrior <- sum(sapply(names(logPriorAndTransFunList), FUN = function(paraName) logPriorAndTransFunList[[paraName]]$logPriorFun(.getCurrentState(currentStateVector = currentState, paraName = paraName))))

  currentState$logPP <- startLogLikValue + currentState$logPrior

  MCMCchainRandomId <- stringi::stri_rand_strings(1, 6) # This should really be different every time, hence its position before set.seed.
  set.seed(control$seed)
  startIterNum <- 1

  if (!is.null(control$chainId)) {
    if (!is.null(control$folderToSaveIntermediateResults)) {
      filesToRestore <- list.files(path = control$folderToSaveIntermediateResults, pattern = control$chainId, full.names = TRUE)
      if (length(filesToRestore) == 0) {
        stop("Specified chain ID, but couldn't find associated files in folderToSaveIntermediateResults. Make sure chainID is correctly specified, or remove chain ID if you want to start the chain from scratch.")
      }
      orderedFilesToRestore <- filesToRestore[order(as.numeric(stringr::str_extract(filesToRestore, pattern = "[[:digit:]]+(?=.Rdata)")))]
      partialChain <- lapply(orderedFilesToRestore, function(fileToRestore) {
        loadName <- load(fileToRestore)
        get(loadName)
      })
      startIterNum <- length(orderedFilesToRestore) + 1
      MCMCchainRandomId[seq_along(orderedFilesToRestore)] <- partialChain
      cat("Resuming MCMC at iteration ", startIterNum, ". \n", sep = "")
    } else {
      stop("Specified a chain ID, but no folder where to find the intermediate results. Remove the chain ID if you want to start the chain from scratch, or specify a value for folderToSaveIntermediateResults in control. \n")
    }
  }
  cat("Chain is called ", MCMCchainRandomId, ". Specify this string in control if you want to resume simulations. \n", sep = "")

  for (MCMCiter in startIterNum:(control$n + control$burnin)) {
    for (paraName in names(startingValues)) {
      if (paraName == "topology") {
        valueForCallToTransFun <- currentState$transmissionTree
      } else if (paraName == "b") {
        valueForCallToTransFun <- log(currentState$phylogeny$edge.length)
      } else if (paraName == "l") {
        valueForCallToTransFun <- log(currentState$transmissionTree$edge.length)
      } else {
        valueForCallToTransFun <- log(currentState$paraValues[[paraName]])
      }

      proposalValueAndTransKernRatio <- logPriorAndTransFunList[[paraName]]$transFun(valueForCallToTransFun)
      updatedLogPrior <- logPriorAndTransFunList[[paraName]]$logPriorFun(proposalValueAndTransKernRatio$value)
      updatedPhylogeny <- currentState$phylogeny
      updatedLogLik <- currentState$logLik

      if (identical(paraName, "topology")) {
        updatedPhylogeny <- proposalValueAndTransKernRatio$value
        updatedLogLik <- logLikFun(updatedPhylogeny, phyDatData, currentState$theta)
      } else if (paraName == "b") {
        updatedPhylogeny$edge.length <- proposalValueAndTransKernRatio$value
        updatedLogLik <- logLikFun(updatedPhylogeny, phyDatData, currentState$theta)
      } else if (paraName == "l") {
        updatedPhylogeny$node.label <- proposalValueAndTransKernRatio$value
      }

      proposalLogPP <- updatedLogLik + updatedLogPrior
      MHratio <- proposalValueAndTransKernRatio$transKernRatio * exp(proposalLogPP - currentState$logPP)
      if (runif(1) <= MHratio) {
        if (paraName %in% c("topology", "l", "b")) {
          currentState$paraValues$phylogeny <- updatedPhylogeny
          ## How to merge the topologies of the transmission tree and the phylogeny?
        } else {
          currentState$paraValues[[paraName]] <- proposalValueAndTransKernRatio$value
        }
        currentState$logLik <- updatedLogLik
        currentState$logPrior <- updatedLogPrior
        currentState$logPP <- updatedLogLik + updatedLogPrior
      }
    }
    MCMCcontainer[[MCMCiter]] <- currentState
  }
  stepSize <- ceiling(control$n * control$thinning)
  itersToKeep <- seq(from = control$burnin + stepSize + 1, to = control$burnin + control$n, by = stepSize)
  MCMCcontainer[itersToKeep]
}

MCMC.control <- function(n = 1e6, thinning = 0.1, burnin = 1e4, seed = 24, folderToSaveIntermediateResults = NULL, chainId = NULL, withinCoalRateKernelSD = 0.5, betweenCoalRateKernelSD = 0.5) {
  list(n = n, thinning = thinning, burnin = burnin, seed = seed, folderForIntermediateResults = folderToSaveIntermediateResults, chainId = chainId, withinCoalRateKernelSD = withinCoalRateKernelSD, betweenCoalRateKernelSD = betweenCoalRateKernelSD)
}

.getCurrentState <- function(currentStateVector, paraName = c("topology", "b", "l", "LambdaWithin", "LambdaBetween")) {
  if (paraName %in% c("topology", "b", "l")) return(currentStateVector$phyloAndTransTree)
  currentStateVector[[paraName]]
}

.clusterFun <- function(MCMCoutput, clusterScoringFun) {

}

.defaultClusterScoringFun <- function(x, genotypingData, phylogeny, covariateFrame) {

}

.topologyTransFun <- function(phylogeny) {
  list(position = phangorn::rNNI(phylogeny), transKernRatio = 1)
}

.convertToTransTree <- function(dualPhyloAndTransTree) {
  transmissionTree <- dualPhyloAndTransTree
  transmissionTree$edge.length <- sapply(dualPhyloAndTransTree$edge.length, '$', "transmissionTree")
  transmissionTree
}

.convertToPhylo <- function(dualPhyloAndTransTree) {
  phylogeny <- dualPhyloAndTransTree
  phylogeny$edge.length <- sapply(dualPhyloAndTransTree$edge.length, '$', "phylogeny")
  phylogeny
}

.topologyLogPriorFun <- function(dualPhyloAndTransTree, Lambda) {
  transmissionTree <- .convertToTransTree(dualPhyloAndTransTree)
  nodeTimes <- sapply(transmissionTree$node.label, function(nodeLabel) nodeLabel$time)
  tipTimes <- sapply(transmissionTree$tip.label, function(tipLabel) tipLabel$time)
  orderIndices <- order(c(nodeTimes, tipTimes), decreasing = TRUE) # We start at the tips, hence the decreasing time.
  nodesIncremented <- lapply(transmissionTree$node.label, function(nodeLabel) {
    nodeLabel$type <- "C"
    nodeLabel
  })

  tipsIncremented <- lapply(transmissionTree$tip.label, function(tipLabel) {
    tipLabel$type <- "T"
    tipLabel
  })
  reorderedTipsAndNodes <- c(nodesIncremented, tipsIncremented)[orderIndices]
  regionNames <- unique(sapply(tipsIncremented, function(x) x$region))
  numRegions <- length(regionNames)
  numLineagesPerRegion <- rep(0, numRegions)
  names(numLineagesPerRegion) <- regionNames
  cumulativeLogProb <- 0

  for (timeIndex in 1:(length(reorderedTipsAndNodes) - 1)) {

    if (reorderedTipsAndNodes[[timeIndex]]$type == "T") {
      numLineagesPerRegion[[reorderedTipsAndNodes[[timeIndex]]$region]] <- numLineagesPerRegion[[reorderedTipsAndNodes[[timeIndex]]$region]] + 1
    } else {
      numLineagesPerRegion[[reorderedTipsAndNodes[[timeIndex]]$region]] <- numLineagesPerRegion[[reorderedTipsAndNodes[[timeIndex]]$region]] - 1
    }
    intervalDuration <- reorderedTipsAndNodes[[timeIndex + 1]]$time - reorderedTipsAndNodes[[timeIndex]]$time
    coalescenceEventIndicator <- reorderedTipsAndNodes[[timeIndex + 1]]$type == "C"

    withinRegionCoalescencePossible <- numLineagesPerRegion > 1

    probsWithinRegionsComponents <- sapply(regionNames[withinRegionCoalescencePossible], function(regionName) {
      -choose(numLineagesPerRegion[[regionName]], 2) * Lambda$within[[regionName]] * intervalDuration
    })
    if (length(probsWithinRegionsComponents) == 0) probsWithinRegionsComponents <- 0
    regionsToCombine <- regionNames[numLineagesPerRegion > 1]
    probsBetweenRegionsComponents <- 0
    if (length(regionsToCombine) > 1) {
      probsBetweenRegionsComponents <- sapply(combn(regionsToCombine, m = 2, simplify = FALSE), FUN = function(regionPair) {
        -prod(numLineagesPerRegion[regionPair]) * Lambda$between * intervalDuration
      })
    }
    cumulativeLogProb <- cumulativeLogProb + sum(probsWithinRegionsComponents) + sum(probsWithinRegionsComponents)
    if (coalescenceEventIndicator) {
      if (reorderedTipsAndNodes[[timeIndex]]$coalescenceType == "within") {
        cumulativeLogProb <- cumulativeLogProb + log(Lambda$within[[reorderedTipsAndNodes[[timeIndex]]$region]])
      } else {
        cumulativeLogProb <- cumulativeLogProb + log(Lambda$between)
      }
    }
  }
  cumulativeLogProb
}

.logPhyloBranchLengthsTransFun <- function(logBranchLengths) {
  transKernSD <- 1
  newPos <- rnorm(n = length(logBranchLengths), mean = logBranchLengths, sd = transKernSD)
  logTransKernRatio <- sum(dnorm(x = logBranchLengths, mean = log(newPos), sd = transKernSD, log = TRUE) - dnorm(x = newPos, mean = logBranchLengths, sd = transKernSD, log = TRUE))
  list(position = newPos, transKernRatio = exp(logTransKernRatio))
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
    list(position = newPos,
         transKernRatio = transKernRatio)
  }
  propInvTransFun <- function(x) {

  }
  gammaShapeTransFun <- function(x) {

  }
  list(K80transitionTransFun, K80transversionTransFun, propInvTransFun, gammaShapeTransFun)
}

# .logPhyloBranchLengthsLogPriorFun <- function(dualPhyloAndTransTree) {
#   branchLengths <- sapply(dualPhyloAndTransTree$edge.length, '$', "phylogeny")
#   dnorm(x = log(branchLengths), mean = -7, sd = 5, log = TRUE)
# }

.logPhyloBranchLengthsLogPriorFun <- function(dualPhyloAndTransTree, mutationRate, numSites) {
  branchLengths <- sapply(dualPhyloAndTransTree$edge.length, '$', "phylogeny")
  transTreeBranchLengths <- sapply(dualPhyloAndTransTree$edge.length, '$', "transmissionTree")
  sum(dpois(x = floor(branchLengths * numSites), lambda = mutationRate * transTreeBranchLengths, log = TRUE))
}

.identifyNodeRegions <- function(transmissionTree) {
  resolveLowerLevels <- function(nodeNumber, resolveValue) {
    currentNodeRegion <- transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region
    splitRegionName <- stringr::str_split(currentNodeRegion, pattern = ",")
    if (resolveValue %in% splitRegionName) {
      transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region <<- resolveValue
      childrenNodes <- phangorn::Children(transmissionTree, nodeNumber)
      childrenForNextStep <- childrenNodes[childrenNodes > ape::Ntip(transmissionTree)]
      if (length(childrenForNextStep) > 0) {
        lapply(childrenForNextStep, resolveLowerLevels, resolveValue = resolveValue)
      }
    }
    NULL
  }

  checkRecursive <- function(nodeNumber) { # nodeNumber will never be for a tip, as tip regions are all known.
    childrenNodes <- phangorn::Children(transmissionTree, nodeNumber)
    getChildrenRegions <- function(childIndex) {
      if (childIndex <= ape::Ntip(transmissionTree)) return(transmissionTree$tip.label[[childIndex]]$region)

      if (is.null(transmissionTree$node.label[[childIndex - ape::Ntip(transmissionTree)]]$region)) return(NA)

      return(transmissionTree$node.label[[childIndex - ape::Ntip(transmissionTree)]]$region)
    }
    childrenRegions <- sapply(childrenNodes, getChildrenRegions)

    if (any(is.na(childrenRegions))) {
      unresolvedChildren <- childrenNodes[is.na(childrenRegions)]
      lapply(unresolvedChildren, checkRecursive)
      childrenRegions <- sapply(childrenNodes, getChildrenRegions)
    }

    splitChildrenRegions <- stringr::str_split(childrenRegions, pattern = ",")
    regionInCommon <- Reduce(intersect, splitChildrenRegions)
    if (length(regionInCommon) > 0) {
      transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region <<- stringr::str_c(sort(regionInCommon), collapse = ",")

      resolveLowerLevels(nodeNumber, resolveValue =  transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region)
    } else {
      transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region <<- stringr::str_c(sort(unique(do.call("c", splitChildrenRegions))), collapse = ",")
    }

    NULL
  }

  # Function works with side-effects
  checkRecursive(ape::Ntip(transmissionTree) + 1)
  resolveNodeRegionsRandomly <- function(nodeNumber) {
    nodeRegions <- stringr::str_split(transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region, pattern = ",")
    nodeChildren <- phangorn::Children(transmissionTree, nodeNumber)
    if (length(nodeRegions) > 1) {
      selectedRegion <- sample(x = nodeRegions, size = 1)
      transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region <<- selectedRegion
      lapply(nodeChildren, FUN = resolveLowerLevels, resolveValue = selectedRegion)
    }
    childrenToConsider <- nodeChildren[nodeChildren > ape::Ntip(transmissionTree)]
    lapply(childrenToConsider, FUN = resolveNodeRegionsRandomly)
    NULL
  }
  # Also works with side-effects
  resolveNodeRegionsRandomly(ape::Ntip(transmissionTree) + 1)
  # identifyWithinBetweenCoalescenceEvents <- function(nodeNumber) {
  #   childrenNodes <- phangorn::Children(transmissionTree, nodeNumber)
  #   childrenRegions <- sapply(childrenNodes, FUN = function(nodeNumber) {
  #     if (nodeNumber > ape::Ntip(transmissionTree)) {
  #       return(transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region)
  #     } else {
  #       return(transmissionTree$tip.label[[nodeNumber]]$region)
  #     }
  #   })
  #   if (length(unique(childrenRegions)) == 1) {
  #     transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$coalescenceType <<- "within"
  #   } else {
  #     transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$coalescenceType <<- "between"
  #   }
  #   childrenToIterateOver <- childrenNodes[childrenNodes > ape::Ntip(transmissionTree)]
  #   lapply(childrenToIterateOver, identifyWithinBetweenCoalescenceEvents) # Returns an empty list if childrenToIterateOver has length 0.
  #   NULL
  # }
  # identifyWithinBetweenCoalescenceEvents(ape::Ntip(transmissionTree) + 1)
  transmissionTree$node.label <- lapply(transmissionTree$node.label, FUN = function(nodeLabelElement) {
    if (stringr::str_detect(nodeLabelElement$region, pattern = ",")) {
      nodeLabelElement$coalescenceType <- "between"
    } else {
      nodeLabelElement$coalescenceType <- "within"
    }
    nodeLabelElement
  })
  transmissionTree
}

.transTreeBranchLengthsLogPrior <- function(dualPhyloAndTransTree, mutationRate) {
  transTreeEdgeLengths <- sapply(dualPhyloAndTransTree$edge.length, '$', "transmissionTree")
  sum(sapply(seq_along(transTreeEdgeLengths), FUN = .transTreeSingleBranchLengthLogPrior))
}

.transTreeSingleBranchLengthLogPrior <- function(transTreeBranchLength) {
  return(1)
}

.transTreeBranchLengthsTransFun <- function(currentState) {
  logNewState <- rnorm(n = length(currentState), mean = log(currentState), sd = 0.5)
  logProbWeightNewState <- sum(dnorm(logNewState, mean = logNewState, sd = 0.5, log = TRUE))
  inverseLogProbWeightNewState <- sum(dnorm(log(currentState), mean = logNewState, sd = 0.5, log = TRUE))
  list(position = exp(logNewState), transKernRatio = exp(sum(inverseLogProbWeightNewState - logProbWeightNewState)))
}

.coalescenceRatesWithinLogPriorFun <- function(withinRates) {
  return(0)
}

.coalescenceRatesWithinTransFun <- function(currentState, sd) {
  logNewState <- rnorm(n = length(currentState), mean = log(currentState), sd = sd)
  logProbWeightNewState <- sum(dnorm(logNewState, mean = logNewState, sd = sd, log = TRUE))
  inverseLogProbWeightNewState <- sum(dnorm(log(currentState), mean = logNewState, sd = sd, log = TRUE))
  list(position = exp(logNewState), transKernRatio = exp(sum(inverseLogProbWeightNewState - logProbWeightNewState)))
}

.coalescenceRatesBetweenLogPriorFun <- function(betweenRate) {
  return(0)
}

.coalescenceRatesBetweenTransFun <- function(currentState, sd) {
  logNewState <- rnorm(n = length(currentState), mean = log(currentState), sd = sd)
  logProbWeightNewState <- sum(dnorm(logNewState, mean = logNewState, sd = sd, log = TRUE))
  inverseLogProbWeightNewState <- sum(dnorm(log(currentState), mean = logNewState, sd = sd, log = TRUE))
  list(position = exp(logNewState), transKernRatio = exp(sum(inverseLogProbWeightNewState - logProbWeightNewState)))
}

plotTransmissionTree <- function(dualPhyloAndTransmissionTree, device = jpeg, filename, argsForDevice = list(), argsForPlotPhylo = list(show.tip.label = TRUE, edge.width = 2)) {
  transmissionTree <- dualPhyloAndTransmissionTree
  transmissionTree$edge.length <- sapply(transmissionTree$edge.length, function(x) x$transmissionTree)
  vertexRegions <- sapply(c(transmissionTree$tip.label, transmissionTree$node.label), function(x) x$region)
  branchRegions <- sapply(seq_along(vertexRegions), FUN = function(vertexIndex) {
    updatedRegion <- vertexRegions[[vertexIndex]]
    if (vertexIndex > ape::Ntip(transmissionTree)) {
      if (stringr::str_detect(transmissionTree$node.label[[vertexIndex - ape::Ntip(transmissionTree)]]$region, pattern = ",")) {
        vertexSiblings <- phangorn::Siblings(transmissionTree, vertexIndex, include.self = TRUE)
        siblingRegions <- sapply(vertexSiblings, function(x) ifelse(x > ape::Ntip(transmissionTree), transmissionTree$node.label[[x - ape::Ntip(transmissionTree)]]$region, transmissionTree$tip.label[[x]]$region))
        splitRegions <- stringr::str_split(siblingRegions, pattern = ",")
        commonRegion <- Reduce(intersect, splitRegions)
        if (length(commonRegion) == 1) {
          updatedRegion <- commonRegion
        } else {
          updatedRegion <- vertexRegions[[vertexIndex]]
        }
      }
    }
    updatedRegion
  })

  transmissionTree$node.label <- sapply(transmissionTree$node.label, function(x) x$time)
  transmissionTree$tip.label <- sapply(transmissionTree$tip.label, function(x) {
    seqLabels <- stringr::str_c("Name = ", x$name)
    timeLabels <- stringr::str_c("Time = ", x$time)
    regionLabels <- stringr::str_c("Region = ", x$region)
    stringr::str_c(seqLabels, timeLabels, regionLabels, sep = ", ")
  })
  reorderedEdgeRegions <- branchRegions[transmissionTree$edge[ , 2]]
  colourScale <- viridis::viridis(n = length(unique(branchRegions)))
  edgeColours <- colourScale[as.numeric(factor(reorderedEdgeRegions))]
  do.call(device, args = c(list(filename = filename), argsForDevice))
  do.call(plot, args = c(list(x = transmissionTree, edge.color = edgeColours), argsForPlotPhylo))
  dev.off()
  NULL
}
