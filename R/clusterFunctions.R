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
    topNodeTime <- min(childNodeTimes) - control$startTransTreeBranchLength + rnorm(1, mean = 0, sd = 1e-10) # The jittering ensure that coalescence events do not happen simultaneously.
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
  numWithinCoalescenceEvents <- rowSums(sapply(seq_along(transmissionTree$node.label) + ape::Ntip(transmissionTree), getCoalVec))

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
      filesToRestoreOrder <- order(as.numeric(stringr::str_extract(filesToRestore, "[:digit:]+(?=.Rdata)")))
      MCMCcontainer[1:length(filesToRestore)] <- do.call("c", lapply(filesToRestore[filesToRestoreOrder], function(filename) {
        loadName <- load(filename)
        get(loadName)
      }))
      startIterNum <- length(MCMCcontainer) * MCMCcontrol$stepSize + 1
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
  phylogeny$node.label <- sapply(seq_along(phylogeny$node.label), function(nodeLabelIndex) stringr::str_c(nodeLabelIndex + ape::Ntip(phylogeny), phylogeny$node.label[[nodeLabelIndex]], sep = ":", collapse = ":"))
  phylogeny
}

.topologyLogPriorFunAndGradient <- function(dualPhyloAndTransTree, Lambda, numMigrationsPoissonPriorMean = 1) {
  transmissionTree <- .convertToTransTree(dualPhyloAndTransTree)
  nodeTimes <- sapply(transmissionTree$node.label, function(nodeLabel) nodeLabel$time)
  tipTimes <- sapply(transmissionTree$tip.label, function(tipLabel) tipLabel$time)
  getMigrationTimesAndChildNodeNums <- function(nodeNum) {
    parentNum <- phangorn::Ancestors(x = transmissionTree, node = nodeNum, type = "parent")
    parentRegion <- .getVertexLabel(transmissionTree, parentNum)$region
    currentRegion <- .getVertexLabel(transmissionTree, nodeNum)$region
    returnValue <- list(childNode = nodeNum, time = NULL)
    if (parentRegion != currentRegion) returnValue$time <- (.getVertexLabel(transmissionTree, parentNum)$time + .getVertexLabel(transmissionTree, nodeNum)$time)/2
    returnValue
  }
  migrationTimesAndChildNodeNums <- lapply(seq_along(transmissionTree$node.label), FUN = getMigrationTimesAndChildNodeNums)
  migrationTimes <- do.call("c", lapply(migrationTimesAndChildNodeNums, function(listElement) listElement$time))
  migrationChildNodes <- do.call("c", lapply(migrationTimesAndChildNodeNums, function(listElement) {
    if (!is.null(listElement$time)) return(listElement$childNode)
    NULL
  }))
  verticesInvolved <- c(1:(ape::Ntip(transmissionTree) + ape::Nnode(transmissionTree)), migrationChildNodes)
  eventType <- rep(c("T", "C", "M"), c(ape::Ntip(transmissionTree), ape::Nnode(transmissionTree), length(migrationTimes)))
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
        numLineagesPerRegion[[.getVertexLabel(transmissionTree, vertexNumber)$region]] <- numLineagesPerRegion[[.getVertexLabel(transmissionTree, vertexNumber)$region]] + 1
      } else if (sortedPhyloStructCodeFrame$type[rowIndex] == "C") {
        numLineagesPerRegion[[.getVertexLabel(transmissionTree, vertexNumber)$region]] <- numLineagesPerRegion[[.getVertexLabel(transmissionTree, vertexNumber)$region]] - 1
      } else {
        childNodeNum <- sortedPhyloStructCodeFrame$vertexNum[rowIndex]
        parentNodeNum <- phangorn::Ancestors(x = transmissionTree, node = childNodeNum, type = "parent")
        childRegion <- .getVertexLabel(transmissionTree, childNodeNum)$region
        parentRegion <- .getVertexLabel(transmissionTree, parentNodeNum)$region
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
      regionCode <- .getVertexLabel(transmissionTree, sortedPhyloStructCodeFrame$vertexNum[[rowsConsidered[[1]] + length(rowsConsidered)]])$region
      cumulativeLogProb <- cumulativeLogProb + log(Lambda[[regionCode]] * intervalDuration)
      cumulativeGradientWithin[[regionCode]] <- cumulativeGradientWithin[[regionCode]] + 1/(Lambda[[regionCode]] * intervalDuration)
    }
  }
  getMigrationIndicator <- function(vertexIndex) {
    vertexLabel <- .getVertexLabel(dualPhyloAndTransTree, vertexIndex)
    currentRegion <- vertexLabel$region
    parentIndex <- phangorn::Ancestors(dualPhyloAndTransTree, vertexIndex, "parent")
    parentRegion <- .getVertexLabel(dualPhyloAndTransTree, parentIndex)$region
    parentRegion != currentRegion
  }
  numMigrations <- sum(sapply(c(1:ape::Ntip(dualPhyloAndTransTree), (ape::Ntip(dualPhyloAndTransTree) + 2):(ape::Ntip(dualPhyloAndTransTree) + ape::Nnode(dualPhyloAndTransTree))), FUN = getMigrationIndicator))
  list(objective = cumulativeLogProb + .numMigrationsLogPriorFunWithMeanPar(x = numMigrations, meanPar = numMigrationsPoissonPriorMean), gradient = cumulativeGradientWithin)
}

.topologyLogPriorFun <- function(dualPhyloAndTransTree, Lambda, numMigrationsPoissonPriorMean) {
  .topologyLogPriorFunAndGradient(dualPhyloAndTransTree, Lambda, numMigrationsPoissonPriorMean)$objective
}

.getNodeLabel <- function(phylogeny, nodeNum) {
  phylogeny$node.label[[nodeNum - ape::Ntip(phylogeny)]]
}

.getTipLabel <- function(phylogeny, tipNum) {
  phylogeny$tip.label[[tipNum]]
}

.numMigrationsLogPriorFunWithMeanPar <- function(x, meanPar = 1) {
  dpois(x = x, lambda = meanPar, log = TRUE)
}

.getVertexLabel <- function(phylogeny, vertexNum) {
  if (vertexNum <= ape::Ntip(phylogeny)) {
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
    nodeRegions <- stringr::str_split(transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region, pattern = ",")[[1]]
    nodeChildren <- phangorn::Children(transmissionTree, nodeNumber)
    if (length(nodeRegions) > 1) {
      selectedRegion <- sample(x = nodeRegions, size = 1)
      transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region <<- selectedRegion
      lapply(nodeChildren[nodeChildren > ape::Ntip(transmissionTree)], FUN = resolveLowerLevels, resolveValue = selectedRegion)
    }
    childrenToConsider <- nodeChildren[nodeChildren > ape::Ntip(transmissionTree)]
    lapply(childrenToConsider, FUN = resolveNodeRegionsRandomly)
    NULL
  }
  transmissionTree$node.label <- lapply(transmissionTree$node.label, function(nodeLabel) {
    nodeLabel$originalRegion <- nodeLabel$region
    nodeLabel
  }) # I copy the information before resolving the tree randomly, just in case.

  # Also works with side-effects
  resolveNodeRegionsRandomly(ape::Ntip(transmissionTree) + 1)
  addCoalescenceType <- function(nodeIndex) {
    nodeCopy <- transmissionTree$node.label[[nodeIndex]]
    childrenVertices <- phangorn::Children(transmissionTree, nodeIndex + ape::Ntip(transmissionTree))
    childrenRegions <- sapply(childrenVertices, FUN = function(childIndex) {
      if (childIndex <= ape::Ntip(transmissionTree)) {
        return(transmissionTree$tip.label[[childIndex]]$region)
      }
      transmissionTree$node.label[[childIndex - ape::Ntip(transmissionTree)]]$region
    })
    if (length(unique(childrenRegions)) > 1) {
      nodeCopy$coalescenceType <- "between"
    } else {
      nodeCopy$coalescenceType <- "within"
    }
    nodeCopy
  }

  transmissionTree$node.label <- lapply(seq_along(transmissionTree$node.label), FUN = addCoalescenceType)
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
  updateNodeTime <- function(nodeIndex) {
    currentTime <- .getVertexLabel(dualPhyloAndTransTree, nodeIndex)$time
    if (nodeIndex != (ape::Ntip(dualPhyloAndTransTree) + 1)) {
      parentIndex <- phangorn::Ancestors(dualPhyloAndTransTree, nodeIndex, "parent")
      parentTime <- .getVertexLabel(dualPhyloAndTransTree, parentIndex)$time
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
  numNodesToModify <- ceiling(propToModify * ape::Nnode(dualPhyloAndTransTree))
  nodesToUpdate <- sample(seq_along(dualPhyloAndTransTree$node.label) + ape::Ntip(dualPhyloAndTransTree), size = numNodesToModify, replace = FALSE)
  updatedNodeTimes <- sapply(nodesToUpdate, updateNodeTime)

  dualPhyloAndTransTree$node.label[nodesToUpdate - ape::Ntip(dualPhyloAndTransTree)] <- lapply(seq_along(nodesToUpdate), function(alongNodeIndex) {
    updatedNode <- .getVertexLabel(dualPhyloAndTransTree, nodesToUpdate[alongNodeIndex])
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
  sapply(childNodeNumbers, function(x) .getVertexLabel(phylogeny = dualPhyloAndTransTree, vertexNum = x)$time) - sapply(parentNodeNumbers, function(x) .getVertexLabel(phylogeny = dualPhyloAndTransTree, vertexNum = x)$time)
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
