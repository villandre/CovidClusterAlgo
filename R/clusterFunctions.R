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
  rootSequenceName,
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

  if (is.null(control$numMigrationsPoissonPriorMean)) {
    control$numMigrationsPoissonPriorMean <- length(unique(seqsRegionStamps)) - 1 # We prefer transmission trees with only one introduction per region; the -1 is for the outgroup, which starts off with a region without an introduction.
  }

  control <- do.call('findBayesianClusters.control', control)

  if (is.null(startingValuePhylo)) {
    MLphyloAndEvoPars <- .genMLphyloAndEvoPars(DNAbinData, rootSequenceName)
    startingValuePhylo <- MLphyloAndEvoPars$phylogeny
    if (is.null(evoParsList)) evoParsList <- MLphyloAndEvoPars$evoParsList
  }

  startingValues <- list()

  # The "edge.length" component of "phylogeny" is a list containing branch lengths for both the transmission tree and phylogeny.

  startingValues$phyloAndTransTree <- .genStartDualPhyloAndTransTree(phylogeny = startingValuePhylo, seqsTimestampsPOSIXct = seqsTimestampsPOSIXct, seqsRegionStamps = seqsRegionStamps, mutationRate = mutationRate, control = control$controlForGenStartTransmissionTree)

  startingValues$Lambda <- .genStartCoalescenceRates(startingValues$phyloAndTransTree, control = control)

  sampledTreesWithPP <- .optimTreeMCMC(
    startingValues = startingValues,
    logLikFun = control$logLikFun,
    DNAbinData = DNAbinData,
    evoParsList = evoParsList,
    covariateFrame = covariateFrame,
    mutationRate = mutationRate,
    control = control)
  # .clusterFun(MCMCoutput = sampledTreesWithPP, clusterScoringFun = control$clusterScoringFun)
  sampledTreesWithPP
}

.performBasicChecks <- function() {

}

findBayesianClusters.control <- function(
  logLikFun = presetPML,
  numMigrationsPoissonPriorMean = 1,
  MCMC.control = MCMC.control(),
  transTreeCondOnPhylo = TRUE,
  controlForGenStartTransmissionTree = .controlForGenStartTransmissionTree()) {
  list(logLikFun = logLikFun, numMigrationsPoissonPriorMean = numMigrationsPoissonPriorMean, MCMC.control = do.call("MCMC.control", MCMC.control), controlForGenStartTransmissionTree = do.call(".controlForGenStartTransmissionTree", controlForGenStartTransmissionTree), transTreeCondOnPhylo = transTreeCondOnPhylo)
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
      if (childIndex > numTips) resolved <- !is.null(transmissionTree$node.label[[childIndex - numTips]]$time) > 0
      resolved
    })
    if (!all(resolvedChildrenTime)) {
      unresolvedChildrenIndices <- childrenIndices[!resolvedChildrenTime]
      sapply(unresolvedChildrenIndices, recursiveFunction)
    }
    childNodeTimes <- sapply(childrenIndices, FUN = function(x) .getVertexLabel(phylogeny = transmissionTree, vertexNum = x)$time)
    topNodeTime <- min(childNodeTimes) - control$startTransTreeBranchLength + rnorm(1, mean = 0, sd = 1e-10) # The jittering ensure that coalescence events do not happen simultaneously.
    transmissionTree$node.label[[nodeIndex - numTips]]$time <<- topNodeTime
    NULL
  }
  startIndex <- numTips + 1
  recursiveFunction(startIndex) # Function works with side-effects.

  getEdgeListElements <- function(edgeNumber) {
    nodeNumber <- transmissionTree$edge[edgeNumber, 2]
    nodeTime <- .getVertexLabel(transmissionTree, vertexNum = nodeNumber)$time
    # parentNumber <- phangorn::Ancestors(transmissionTree, node = nodeNumber, type = "parent")
    parentNumber <- transmissionTree$edge[match(nodeNumber, transmissionTree$edge[ , 2]), 1]
    parentTime <- transmissionTree$node.label[[parentNumber - numTips]]$time
    list(transmissionTree = nodeTime - parentTime, phylogeny = transmissionTree$edge.length[[edgeNumber]])
  }
  transmissionTree$edge.length <- lapply(1:nrow(transmissionTree$edge), FUN = getEdgeListElements)
  transmissionTree <- .identifyNodeRegions(transmissionTree)
  transmissionTree
}

.genStartCoalescenceRates <- function(dualPhyloAndTransTree, control) {
  transmissionTree <- .convertToTransTree(dualPhyloAndTransTree)
  regionNames <- unique(sapply(transmissionTree$tip.label, '[[', "region"))
  getCoalVec <- function(nodeIndex) {
    basicCoalVec <- rep(0, length(regionNames))
    names(basicCoalVec) <- regionNames
    childNodes <- phangorn::Children(transmissionTree, nodeIndex)
    childRegions <- sapply(childNodes, function(x) .getVertexLabel(transmissionTree, x)$region)
    if (length(unique(childRegions)) == 1) basicCoalVec[[childRegions[[1]]]] <- 1
    basicCoalVec
  }
  numWithinCoalescenceEvents <- rowSums(sapply(seq_along(transmissionTree$node.label) + length(dualPhyloAndTransTree$tip.label), getCoalVec))

  startWithin <- numWithinCoalescenceEvents/(sum(transmissionTree$edge.length) * length(regionNames))
  names(startWithin) <- regionNames

  funForOptim <- function(x) {
    names(x) <- regionNames # For some reason, the name attribute of the argument x is stripped when the function is called within lbfgs...
    x <- exp(x)
    returnValue <- .topologyLogPriorFunAndGradient(dualPhyloAndTransTree = dualPhyloAndTransTree, Lambda = x, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean)
     # return(list("objective" = -returnValue$objective,
     #            "gradient" = -returnValue$gradient[regionNames]) ) # lbfgs minimises, hence the "-". Mimicking the syntax in the nloptr vignette. FIX THE GRADIENT-BASED ESTIMATION.
    -returnValue$objective
  }

  optimResult <- nloptr::lbfgs(x0 = log(startWithin), fn = funForOptim)
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
  logLikFun,
  DNAbinData,
  evoParsList,
  mutationRate,
  covariateFrame,
  control) {
  MCMCcontrol <- do.call(MCMC.control, control$MCMC.control)
  phyDatData <- phangorn::as.phyDat(DNAbinData)
  # The following assumes that the order of tip labels matches the order of sequences in DNAbinData.
  phyloObj <- .convertToPhylo(startingValues$phyloAndTransTree)
  startLogLikValue <- logLikFun(phyloObj, phyDatData, evoParsList)

  MCMCcontainer <- vector(mode = "list", length = floor((MCMCcontrol$n + MCMCcontrol$burnin)/MCMCcontrol$stepSize))

  currentState <- list(paraValues = startingValues, logLik = startLogLikValue)
  bLogPrior <- lLogPrior <- function(x) 0
  if (control$transTreeCondOnPhylo) {
    bLogPrior <- function(x) .phyloBranchLengthsLogPriorFun(dualPhyloAndTransTree = x)
    lLogPrior <- function(x) .transTreeBranchLengthsConditionalLogPrior(dualPhyloAndTransTree = x, mutationRate = mutationRate, numSites = ncol(DNAbinData))
  } else {
    bLogPrior <- function(x) .phyloBranchLengthsConditionalLogPriorFun(dualPhyloAndTransTree = x, mutationRate = mutationRate, numSites = ncol(DNAbinData))
    lLogPrior <- function(x) .transTreeBranchLengthsLogPrior(dualPhyloAndTransTree = x, mutationRate = mutationRate)
  }

  # Functions are defined here with an external dependence on currentState, but this is voluntary.
  logPriorAndTransFunList <- list(
    topology = list(
      logPriorFun = function(x, Lambda = currentState$paraValues$Lambda) .topologyLogPriorFun(dualPhyloAndTransTree = x, Lambda = Lambda, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean),
      transFun = .topologyTransFun),
    b = list(
      logPriorFun = bLogPrior,
      transFun = function(x) .phyloBranchLengthsTransFun(x, deflateCoef = MCMCcontrol$phyloBranchLengthsDeflateCoef, propToModify = MCMCcontrol$propPhyloBranchesToModify)),
    l = list(
      logPriorFun = lLogPrior,
      transFun = function(x) .transTreeBranchLengthsTransFun(dualPhyloAndTransTree = x, tuningPara = MCMCcontrol$transTreeTuningPara, propToModify = MCMCcontrol$propTransTreeBranchesToModify)),
    Lambda = list(
      logPriorFun = .coalescenceRatesLogPriorFun,
      transFun = function(x) .coalescenceRatesTransFun(x, sd = MCMCcontrol$coalRateKernelSD)))
  currentState$logPrior <- sapply(names(logPriorAndTransFunList), FUN = function(paraName) logPriorAndTransFunList[[paraName]]$logPriorFun(.getCurrentState(currentStateVector = currentState, paraName = paraName)))
  names(currentState$logPrior) <- names(logPriorAndTransFunList)

  currentState$logPP <- startLogLikValue + sum(currentState$logPrior)

  MCMCchainRandomId <- stringi::stri_rand_strings(1, 6) # This should really be different every time, hence its position before set.seed.
  set.seed(MCMCcontrol$seed)
  startIterNum <- 1

  if (!is.null(MCMCcontrol$chainId)) {
    if (!is.null(MCMCcontrol$folderToSaveIntermediateResults)) {
      filesToRestore <- list.files(path = MCMCcontrol$folderToSaveIntermediateResults, pattern = MCMCcontrol$chainId, full.names = TRUE)
      if (length(filesToRestore) == 0) {
        stop("Specified chain ID, but couldn't find associated files in folderToSaveIntermediateResults. Make sure chainID is correctly specified, or remove chain ID if you want to start the chain from scratch.")
      }
      iterNums <- as.numeric(stringr::str_extract(filesToRestore, "[:digit:]+(?=.Rdata)"))
      filesToRestoreOrder <- order(iterNums)
      MCMCcontainer[1:length(filesToRestore)] <- do.call("c", lapply(filesToRestore[filesToRestoreOrder], function(filename) {
        loadName <- load(filename)
        get(loadName)
      }))
      startIterNum <- max(iterNums) + 1
      MCMCchainRandomId <- MCMCcontrol$chainId
      cat("Resuming MCMC at iteration ", startIterNum, ". \n", sep = "")
    } else {
      stop("Specified a chain ID, but no folder where to find the intermediate results. Remove the chain ID if you want to start the chain from scratch, or specify a value for folderToSaveIntermediateResults in MCMCcontrol. \n")
    }
  } else {
    cat("Launching MCMC... \n")
  }
  cat("Chain is called ", MCMCchainRandomId, ". Specify this string in MCMCcontrol if you want to resume simulations.\n", sep = "")
  cat("Starting values for log-priors: \n")
  print(currentState$logPrior)
  cat("Starting value for log-lik.: \n")
  print(currentState$logLik)
  for (MCMCiter in startIterNum:(MCMCcontrol$n + MCMCcontrol$burnin)) {
    if ((MCMCiter %% MCMCcontrol$print.frequency) == 0) cat("This is MCMC iteration ", MCMCiter, sep = "", ".\n")
    for (paraName in names(logPriorAndTransFunList)) {
      proposalValueAndTransKernRatio <- logPriorAndTransFunList[[paraName]]$transFun(.getCurrentState(currentState, paraName))
      updatedLogPrior <- currentState$logPrior
      updatedLogPrior[[paraName]] <- logPriorAndTransFunList[[paraName]]$logPriorFun(proposalValueAndTransKernRatio$value)
      updatedDualTree <- .getCurrentState(currentState, "topology")
      updatedLogLik <- currentState$logLik

      if (paraName %in% c("topology", "b")) {
        updatedDualTree <- proposalValueAndTransKernRatio$value
        updatedLogLik <- logLikFun(.convertToPhylo(updatedDualTree), phyDatData, evoParsList)
      } else if (paraName == "l") { # The prior for 'b' is conditioned on 'l'.
        updatedDualTree <- proposalValueAndTransKernRatio$value
        updatedLogPrior[["topology"]] <- logPriorAndTransFunList[["topology"]]$logPriorFun(updatedDualTree)
        updatedLogPrior[["b"]] <- logPriorAndTransFunList[["b"]]$logPriorFun(updatedDualTree)
      } else if (paraName == "Lambda") {
        updatedLogPrior[["topology"]] <- logPriorAndTransFunList[["topology"]]$logPriorFun(updatedDualTree, Lambda = proposalValueAndTransKernRatio$value)
      }
      proposalLogPP <- updatedLogLik + sum(updatedLogPrior)
      MHratio <- proposalValueAndTransKernRatio$transKernRatio * exp(proposalLogPP - currentState$logPP)
      if (runif(1) <= MHratio) {
        # if (paraName == "b") cat("Modified phylo. branch lengths! \n")
        # if (paraName == "l") cat("Modified trans. tree branch lengths! \n")
        # cat("Move accepted! Parameter name: ", paraName, "\n", sep = "")
        # cat("Updated log-priors: \n")
        # print(updatedLogPrior)
        # cat("Updated log-lik.: \n")
        # cat(updatedLogLik)
        if (paraName %in% c("topology", "l", "b")) {
          currentState$paraValues$phyloAndTransTree <- updatedDualTree
        } else {
          currentState$paraValues[[paraName]] <- proposalValueAndTransKernRatio$value
        }
        currentState$logLik <- updatedLogLik
        currentState$logPrior <- updatedLogPrior
        currentState$logPP <- proposalLogPP
      }
    }
    if ((MCMCiter %% MCMCcontrol$stepSize) == 0) {
      MCMCcontainer[[MCMCiter/MCMCcontrol$stepSize]] <- currentState
      save(currentState, file = paste(MCMCcontrol$folderToSaveIntermediateResults, "/chainID_", MCMCchainRandomId, "_atIter", MCMCiter, ".Rdata", sep = ""), compress = TRUE)
    }
  }
  cat("MCMC complete. Finalising... \n")
  elementsToKeep <- seq(
    from = floor(MCMCcontrol$burnin/MCMCcontrol$stepSize) + 1,
    to = floor((MCMCcontrol$burnin + MCMCcontrol$n)/MCMCcontrol$stepSize))
  MCMCcontainer[elementsToKeep]
}

MCMC.control <- function(n = 1e6, stepSize = 50, burnin = 1e4, seed = 24, folderToSaveIntermediateResults = NULL, save.frequency = 100, chainId = NULL, coalRateKernelSD = 0.5, print.frequency = 10, phyloBranchLengthsDeflateCoef = 0.98, propPhyloBranchesToModify = 0.1, transTreeTuningPara = 0.1, propTransTreeBranchesToModify = 0.2) {
  list(n = n, stepSize = stepSize, burnin = burnin, seed = seed, folderToSaveIntermediateResults = folderToSaveIntermediateResults, chainId = chainId, coalRateKernelSD = coalRateKernelSD, print.frequency = print.frequency, save.frequency = save.frequency, phyloBranchLengthsDeflateCoef = phyloBranchLengthsDeflateCoef, propPhyloBranchesToModify = propPhyloBranchesToModify, transTreeTuningPara = transTreeTuningPara, propTransTreeBranchesToModify = propTransTreeBranchesToModify)
}

.getCurrentState <- function(currentStateVector, paraName = c("topology", "b", "l", "Lambda")) {
  if (paraName %in% c("topology", "b", "l")) return(currentStateVector$paraValues$phyloAndTransTree)
  currentStateVector$paraValues[[paraName]]
}

.clusterFun <- function(MCMCoutput, clusterScoringFun) {

}

.defaultClusterScoringFun <- function(x, genotypingData, phylogeny, covariateFrame) {

}

.topologyTransFun <- function(dualPhyloAndTransTree) {
  counter <- 0
  repeat {
    counter <- counter + 1
    proposedMove <- phangorn::rNNI(dualPhyloAndTransTree)
    newEdgeLengths <- .deriveTransTreeEdgeLengths(proposedMove)
    if (all(newEdgeLengths >= 0)) break
    if (counter == 50) {
      cat("Could not find a suitable move in the topological space... \n")
      proposedMove <- dualPhyloAndTransTree
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

.convertToTransTree <- function(dualPhyloAndTransTree) {
  transmissionTree <- dualPhyloAndTransTree
  transmissionTree$edge.length <- sapply(dualPhyloAndTransTree$edge.length, '[[', "transmissionTree")
  transmissionTree
}

.convertToPhylo <- function(dualPhyloAndTransTree) {
  phylogeny <- dualPhyloAndTransTree
  phylogeny$edge.length <- sapply(dualPhyloAndTransTree$edge.length, '[[', "phylogeny")
  phylogeny$tip.label <- sapply(phylogeny$tip.label, function(x) x$name) # The other components must be discarded for now, to allow pml to work.
  phylogeny$node.label <- sapply(seq_along(phylogeny$node.label), function(nodeLabelIndex) stringr::str_c(nodeLabelIndex + length(phylogeny$tip.label), phylogeny$node.label[[nodeLabelIndex]], sep = ":", collapse = ":"))
  phylogeny
}

.topologyLogPriorFunAndGradient <- function(dualPhyloAndTransTree, Lambda, numMigrationsPoissonPriorMean = 1) {
  numTips <- length(dualPhyloAndTransTree$tip.label)
  numNodes <- length(dualPhyloAndTransTree$node.label)
  transmissionTree <- .convertToTransTree(dualPhyloAndTransTree)
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
  cumulativeGradientWithin <- rep(0, length(Lambda))
  names(cumulativeGradientWithin) <- names(Lambda)
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
        numLineagesPerRegion[[childRegion]] <- numLineagesPerRegion[[childRegion]] + 1
        numLineagesPerRegion[[parentRegion]] <- numLineagesPerRegion[[parentRegion]] + 1
      }
    }
    intervalDuration <- sortedPhyloStructCodeFrame$time[[rowsConsidered[[1]]]] - sortedPhyloStructCodeFrame$time[[rowsConsidered[[1]] + length(rowsConsidered)]]

    withinRegionCoalescencePossible <- numLineagesPerRegion > 1
    withinGradient <- rep(0, length(Lambda))
    names(withinGradient) <- names(Lambda)
    cumulWithinRatePerRegion <- sapply(regionNames[withinRegionCoalescencePossible], function(regionName) {
      choose(numLineagesPerRegion[[regionName]], 2) * Lambda[[regionName]] * intervalDuration
    })
    names(cumulWithinRatePerRegion) <- regionNames[withinRegionCoalescencePossible]
    withinGradient[names(cumulWithinRatePerRegion)] <- -cumulWithinRatePerRegion/Lambda[names(cumulWithinRatePerRegion)]
    if (length(cumulWithinRatePerRegion) == 0) cumulWithinRatePerRegion <- 0
    totalRate <- sum(cumulWithinRatePerRegion)
    cumulativeLogProb <- cumulativeLogProb - totalRate # That takes into account the exponent of exp(), hence the "-".
    cumulativeGradientWithin <- cumulativeGradientWithin + withinGradient[names(cumulativeGradientWithin)]

    if (sortedPhyloStructCodeFrame$type[[rowsConsidered[[1]] + length(rowsConsidered)]] == "C") {
      regionCode <- transmissionTree$node.label[[sortedPhyloStructCodeFrame$vertexNum[[rowsConsidered[[1]] + length(rowsConsidered)]] - numTips]]$region
      cumulativeLogProb <- cumulativeLogProb + log(Lambda[[regionCode]] * intervalDuration)
      cumulativeGradientWithin[[regionCode]] <- cumulativeGradientWithin[[regionCode]] + 1/(Lambda[[regionCode]] * intervalDuration)
    }
  }
  getMigrationIndicator <- function(vertexIndex) {
    if (vertexIndex > numTips) {
      currentRegion <- dualPhyloAndTransTree$node.label[[vertexIndex - numTips]]$region
    } else {
      currentRegion <- dualPhyloAndTransTree$tip.label[[vertexIndex]]$region
    }
    # parentIndex <- phangorn::Ancestors(dualPhyloAndTransTree, vertexIndex, "parent")
    parentIndex <- dualPhyloAndTransTree$edge[match(vertexIndex, dualPhyloAndTransTree$edge[ , 2]), 1]
    # parentRegion <- .getVertexLabel(dualPhyloAndTransTree, parentIndex)$region
    parentRegion <- dualPhyloAndTransTree$node.label[[parentIndex - numTips]]$region
    parentRegion != currentRegion
  }
  numMigrations <- sum(sapply(c(1:numTips, (numTips + 2):(numTips + numNodes)), FUN = getMigrationIndicator))
  list(objective = cumulativeLogProb + .numMigrationsLogPriorFunWithMeanPar(x = numMigrations, meanPar = numMigrationsPoissonPriorMean), gradient = cumulativeGradientWithin)
}

.topologyLogPriorFun <- function(dualPhyloAndTransTree, Lambda, numMigrationsPoissonPriorMean) {
  .topologyLogPriorFunAndGradient(dualPhyloAndTransTree, Lambda, numMigrationsPoissonPriorMean)$objective
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

# .phyloBranchLengthsConditionalLogPriorFun <- function(dualPhyloAndTransTree) {
#   branchLengths <- sapply(dualPhyloAndTransTree$edge.length, '$', "phylogeny")
#   dnorm(x = log(branchLengths), mean = -7, sd = 5, log = TRUE)
# }

.phyloBranchLengthsConditionalLogPriorFun <- function(dualPhyloAndTransTree, mutationRate, numSites) {
  branchLengths <- sapply(dualPhyloAndTransTree$edge.length, '[[', "phylogeny")
  transTreeBranchLengths <- sapply(dualPhyloAndTransTree$edge.length, '[[', "transmissionTree")
  sum(dpois(x = floor(branchLengths * numSites), lambda = mutationRate * transTreeBranchLengths, log = TRUE))
}

.phyloBranchLengthsLogPriorFun <- function(dualPhyloAndTransTree) {
  sum(sapply(dualPhyloAndTransTree$edge.length, function(edgeLengthElement) return(0))) # We're assuming a uniform prior. The densities are non-standardised but it doesn't matter for MCMC
}

.phyloBranchLengthsTransFun <- function(dualPhyloAndTransTree, deflateCoef = 0.98, propToModify = 0.1) {
  logBranchLengths <- log(sapply(dualPhyloAndTransTree$edge.length, FUN = '[[', "phylogeny"))
  nonZeroLengthBranchPos <- which(!is.infinite(logBranchLengths))
  numToModify <- ceiling(propToModify * length(nonZeroLengthBranchPos))
  branchesToModify <- sample(nonZeroLengthBranchPos, size = numToModify, replace = FALSE)
  modCoef <- c(deflateCoef, 1/deflateCoef)
  newPos <- logBranchLengths
  #newPos[nonZeroLengthBranchPos] <- rnorm(n = length(nonZeroLengthBranchPos), mean = logBranchLengths[nonZeroLengthBranchPos], sd = transKernSD)
  newPos[branchesToModify] <- sample(x = modCoef, size = length(branchesToModify), replace = TRUE) * newPos[branchesToModify]
  # logTransKernRatio <- sum(dnorm(x = logBranchLengths[nonZeroLengthBranchPos], mean = newPos[nonZeroLengthBranchPos], sd = transKernSD, log = TRUE) - dnorm(x = newPos[nonZeroLengthBranchPos], mean = logBranchLengths[nonZeroLengthBranchPos], sd = transKernSD, log = TRUE))
  dualPhyloAndTransTree$edge.length <- lapply(seq_along(dualPhyloAndTransTree$edge.length), function(edgeIndex) {
    edgeElement <- dualPhyloAndTransTree$edge.length[[edgeIndex]]
    edgeElement$phylogeny <- exp(newPos[[edgeIndex]])
    edgeElement
  })

  list(value = dualPhyloAndTransTree, transKernRatio = 1)
}

.identifyNodeRegions <- function(transmissionTree) {
  numTips <- length(transmissionTree$tip.label)
  resolveLowerLevels <- function(nodeNumber, resolveValue) {
    currentNodeRegion <- transmissionTree$node.label[[nodeNumber - numTips]]$region
    splitRegionName <- stringr::str_split(currentNodeRegion, pattern = ",")
    if (resolveValue %in% splitRegionName) {
      transmissionTree$node.label[[nodeNumber - numTips]]$region <<- resolveValue
      childrenNodes <- phangorn::Children(transmissionTree, nodeNumber)
      childrenForNextStep <- childrenNodes[childrenNodes > numTips]
      if (length(childrenForNextStep) > 0) {
        lapply(childrenForNextStep, resolveLowerLevels, resolveValue = resolveValue)
      }
    }
    NULL
  }

  checkRecursive <- function(nodeNumber) { # nodeNumber will never be for a tip, as tip regions are all known.
    childrenNodes <- phangorn::Children(transmissionTree, nodeNumber)
    getChildrenRegions <- function(childIndex) {
      if (childIndex <= numTips) return(transmissionTree$tip.label[[childIndex]]$region)

      if (is.null(transmissionTree$node.label[[childIndex - numTips]]$region)) return(NA)

      return(transmissionTree$node.label[[childIndex - numTips]]$region)
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
      transmissionTree$node.label[[nodeNumber - numTips]]$region <<- stringr::str_c(sort(regionInCommon), collapse = ",")

      resolveLowerLevels(nodeNumber, resolveValue = transmissionTree$node.label[[nodeNumber - numTips]]$region)
    } else {
      transmissionTree$node.label[[nodeNumber - numTips]]$region <<- stringr::str_c(sort(unique(do.call("c", splitChildrenRegions))), collapse = ",")
    }
    NULL
  }
  # Function works with side-effects
  checkRecursive(ape::Ntip(transmissionTree) + 1)

  resolveNodeRegionsRandomly <- function(nodeNumber) {
    if (nodeNumber <= numTips) return(NULL)
    nodeRegions <- stringr::str_split(transmissionTree$node.label[[nodeNumber - numTips]]$region, pattern = ",")[[1]]
    nodeChildren <- phangorn::Children(transmissionTree, nodeNumber)
    if (length(nodeRegions) > 1) {
      selectedRegion <- sample(x = nodeRegions, size = 1)
      transmissionTree$node.label[[nodeNumber - numTips]]$region <<- selectedRegion
      lapply(nodeChildren[nodeChildren > numTips], FUN = resolveLowerLevels, resolveValue = selectedRegion)
    }
    childrenToConsider <- nodeChildren[nodeChildren > numTips]
    lapply(childrenToConsider, FUN = resolveNodeRegionsRandomly)
    NULL
  }
  transmissionTree$node.label <- lapply(transmissionTree$node.label, function(nodeLabel) {
    nodeLabel$originalRegion <- nodeLabel$region
    nodeLabel
  }) # I copy the information before resolving the tree randomly, just in case.

  # Also works with side-effects
  resolveNodeRegionsRandomly(numTips + 1)

  transmissionTree
}

.clearNodeRegions <- function(dualPhyloAndTransTree) {
  dualPhyloAndTransTree$node.label <- lapply(dualPhyloAndTransTree$node.label, function(nodeLabel) {
    nodeLabel$region <- NA
    nodeLabel
  })
  dualPhyloAndTransTree
}

.transTreeBranchLengthsLogPrior <- function(dualPhyloAndTransTree, mutationRate) {
  transTreeSingleBranchLengthLogPrior <- function(transTreeBranchLength) {
    return(0)
  }
  sum(sapply(seq_along(dualPhyloAndTransTree$edge.length), FUN = transTreeSingleBranchLengthLogPrior))
}

.transTreeBranchLengthsConditionalLogPrior <- function(dualPhyloAndTransTree, mutationRate, numSites) {

  transTreeSingleBranchLengthConditionalLogPrior <- function(branchIndex) {
    phyloBranchLength <- dualPhyloAndTransTree$edge.length[[branchIndex]]$phylogeny
    transTreeBranchLength <- dualPhyloAndTransTree$edge.length[[branchIndex]]$transmissionTree

    if (identical(phyloBranchLength, 0)) {
      logProb <- dexp(transTreeBranchLength, rate = mutationRate/2, log = TRUE)
    } else {
      logProb <- dexp(transTreeBranchLength, rate = mutationRate * numSites, log = TRUE)
    }
    logProb
  }
  sum(sapply(seq_along(dualPhyloAndTransTree$edge.length), FUN = transTreeSingleBranchLengthConditionalLogPrior))
}

.transTreeBranchLengthsTransFun <- function(dualPhyloAndTransTree, tuningPara = 0.1, propToModify = 0.1) {
  numTips <- length(dualPhyloAndTransTree$tip.label)
  numNodes <- length(dualPhyloAndTransTree$node.label)
  updateNodeTime <- function(nodeIndex) {
    currentTime <- dualPhyloAndTransTree$node.label[[nodeIndex - numTips]]$time
    if (nodeIndex != (numTips + 1)) {
      # parentIndex <- phangorn::Ancestors(dualPhyloAndTransTree, nodeIndex, "parent")
      parentIndex <- dualPhyloAndTransTree$edge[match(nodeIndex, dualPhyloAndTransTree$edge[ , 2]), 1]
      parentTime <- dualPhyloAndTransTree$node.label[[parentIndex - numTips]]$time
    } else {
      parentTime <- currentTime
    }
    childrenIndices <- phangorn::Children(dualPhyloAndTransTree, nodeIndex)
    childrenTimes <- sapply(childrenIndices, function(childIndex) .getVertexLabel(dualPhyloAndTransTree, childIndex)$time)
    lowerBound <- currentTime - (currentTime - parentTime) * tuningPara
    upperBound <- currentTime + (min(childrenTimes) - currentTime) * tuningPara
    updatedTime <- runif(1, lowerBound, upperBound)
    updatedTime
  }
  numNodesToModify <- ceiling(propToModify * numNodes)
  nodesToUpdate <- sample(seq_along(dualPhyloAndTransTree$node.label) + numTips, size = numNodesToModify, replace = FALSE)
  updatedNodeTimes <- sapply(nodesToUpdate, updateNodeTime)

  dualPhyloAndTransTree$node.label[nodesToUpdate - numTips] <- lapply(seq_along(nodesToUpdate), function(alongNodeIndex) {
    updatedNode <- dualPhyloAndTransTree$node.label[[nodesToUpdate[alongNodeIndex] - numTips]]
    updatedNode$time <- updatedNodeTimes[[alongNodeIndex]]
    updatedNode
  })
  updatedTransTreeEdges <- .deriveTransTreeEdgeLengths(dualPhyloAndTransTree)
  dualPhyloAndTransTree$edge.length <- lapply(seq_along(dualPhyloAndTransTree$edge.length), function(edgeIndex) {
    updatedEdgeLength <- dualPhyloAndTransTree$edge.length[[edgeIndex]]
    updatedEdgeLength$transmissionTree <- updatedTransTreeEdges[[edgeIndex]]
    updatedEdgeLength
  })

  list(value = dualPhyloAndTransTree, transKernRatio = 1) # The uniform transition kernel is symmetrical.
}

.deriveTransTreeEdgeLengths <- function(dualPhyloAndTransTree) {
  childNodeNumbers <- dualPhyloAndTransTree$edge[ , 2]
  parentNodeNumbers <- dualPhyloAndTransTree$edge[ , 1]
  numTips <- length(dualPhyloAndTransTree$tip.label)
  sapply(childNodeNumbers, function(x) .getVertexLabel(phylogeny = dualPhyloAndTransTree, vertexNum = x)$time) - sapply(parentNodeNumbers, function(x) dualPhyloAndTransTree$node.label[[x - numTips]]$time)
}

.getVertexTimes <- function(dualPhyloAndTransTree, type = c("all", "nodes", "tips")) {
  type <- type[[1]]
  labelList <- dualPhyloAndTransTree$tip.label
  if (type == "all") {
    labelList <- c(dualPhyloAndTransTree$tip.label, dualPhyloAndTransTree$node.label)
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

plotTransmissionTree <- function(dualPhyloAndTransmissionTree, device = jpeg, filename, argsForDevice = list(), argsForPlotPhylo = list(show.tip.label = TRUE, edge.width = 2)) {
  numTips <- length(dualPhyloAndTransmissionTree$tip.label)
  transmissionTree <- dualPhyloAndTransmissionTree
  transmissionTree$edge.length <- sapply(transmissionTree$edge.length, function(x) x$transmissionTree)
  vertexRegions <- sapply(c(transmissionTree$tip.label, transmissionTree$node.label), function(x) x$region)
  branchRegions <- sapply(seq_along(vertexRegions), FUN = function(vertexIndex) {
    updatedRegion <- vertexRegions[[vertexIndex]]
    if (vertexIndex > numTips) {
      if (stringr::str_detect(.getVertexLabel(transmissionTree, vertexIndex)$region, pattern = ",")) {
        vertexSiblings <- phangorn::Siblings(transmissionTree, vertexIndex, include.self = TRUE)
        siblingRegions <- sapply(vertexSiblings, function(x) .getVertexLabel(phylogeny = transmissionTree, vertexNum = x)$region)
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
  do.call(device, args = c(list(filename), argsForDevice))
  do.call(plot, args = c(list(x = transmissionTree, edge.color = edgeColours), argsForPlotPhylo))
  dev.off()
  NULL
}
