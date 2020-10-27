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
  startingTreeFromChainIter = NULL,
  evoParsList = NULL, # Must be specified if startingTreeFromChainIter is specified. It should be found in folderToSaveIntermediateResults if it was specified.
  clusterScoringFun = NULL,
  control = list()) {
  #.performBasicChecks()
  if (!is.null(startingTreeFromChainIter) & is.null(evoParsList)) {
    stop("You have specified a starting value for the tree (startingTreeFromChainIter) without the associated evolutionary parameters (evoParsList). \n")
  }
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
  startingTree <- startingTreeFromChainIter
  if (is.null(control$numMigrationsPoissonPriorMean)) {
    control$numMigrationsPoissonPriorMean <- length(unique(seqsRegionStamps)) - 1 # We prefer transmission trees with only one introduction per region; the -1 is for the outgroup, which starts off with a region without an introduction.
  }
  chainId <- chainIdCopy <- control$MCMC.control$chainId # Copy exists to distinguish a new id (which will be generated later on) and an id used previously.
  if (is.null(control$MCMC.control$folderToSaveIntermediateResults) | is.null(chainId)) {
    chainId <- stringi::stri_rand_strings(1, 6)
    if (is.null(startingTree)) {
      MLphyloAndEvoPars <- .genMLphyloAndEvoPars(DNAbinData, rootSequenceName)
      startingValuePhylo <- MLphyloAndEvoPars$phylogeny
      timestampsInDays <- as.numeric(seqsTimestampsPOSIXct)/86400
      names(timestampsInDays) <- names(seqsTimestampsPOSIXct)
      startingTree <- .genStartPhyloAndTransTreeAlt(phylogeny = startingValuePhylo, timestampsInDays = timestampsInDays, regionStamps = seqsRegionStamps, logClockRatePriorMean = log(perSiteClockRate), estRootTime = estRootTime, control = control)
      evoParsList <- MLphyloAndEvoPars$evoParsList
    }

    if (!is.null(control$MCMC.control$folderToSaveIntermediateResults)) {
      filename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_controlParameters.Rdata", sep = "")
      save(control, file = filename)
      evoParsFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_evoParameters.Rdata", sep = "")
      save(evoParsList, file = evoParsFilename)
      startingTreeFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_startingTree.Rdata", sep = "")
      save(startingTree, file = startingTreeFilename)
    }
  } else {
    cat("Restoring chain", chainId, "control parameters... \n", sep = " ")
    burnin <- control$MCMC.control$burnin
    n <- control$MCMC.control$n
    stepSize <- control$MCMC.control$stepSize
    nIterPerSweep <- control$MCMC.control$nIterPerSweep

    controlFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_controlParameters.Rdata", sep = "")
    evoParsFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_evoParameters.Rdata", sep = "")
    startingTreeFilename <- paste(control$MCMC.control$folderToSaveIntermediateResults, "/chain_", chainId, "_startingTree.Rdata", sep = "")
    load(controlFilename) # This will restore "control" to its original values.
    control$MCMC.control$chainId <- chainIdCopy
    # The user should be allowed to change these chain parameters, as they do not affect the transitions.
    if (!is.null(n)) control$MCMC.control$n <- n
    if (!is.null(burnin)) control$MCMC.control$burnin <- burnin
    if (!is.null(stepSize)) control$MCMC.control$stepSize <- stepSize
    if (!is.null(nIterPerSweep)) control$MCMC.control$nIterPerSweep <- nIterPerSweep

    load(evoParsFilename) # This will restore evoParsList.
    load(startingTreeFilename) # This will restore startingTree.
  }

  sampledTreesWithPP <- .optimTreeMCMC(
    startingTree = startingTree,
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
    phyloAndTransTree$Lambda <- x
    returnValue <- .topologyLogPriorFun(phyloAndTransTree = phyloAndTransTree, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean, estRootTime = estRootTime, control = control)
    -returnValue
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
  startingTree,
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
  phyloObj <- .convertToPhylo(startingTree)
  startLogLikValue <- logLikFun(phyloObj, phyDatData, evoParsList)
  nodeOrder <- .getOrderedNodeNumsForTimeUpdate(startingTree)
  currentState <- lapply(1:MCMCcontrol$nChains, function(x) list(phyloAndTransTree = startingTree, logLik = startLogLikValue, nodeUpdateOrder = nodeOrder - ape::Ntip(startingTree))) # All chains start in the same state.
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
    cat("Launching MCMC... \n")
  }

  logPriorAndTransFunList <- list(
    b = list( # The mean here corresponds to 14 days for the time between the time the sequence is sampled and the time at which the first transmission linking it to another sequence in the sample occurred.
      logPriorFun =  function(x) .phyloBranchLengthsLogPriorFun(phyloAndTransTree = x, extBranchLengthPriorMean = perSiteClockRate * 14, extBranchLengthPriorSD = perSiteClockRate * 50),
      transFun = function(x) .phyloBranchLengthsTransFun(x, logTransKernSD = MCMCcontrol$phyloBranchLengthsLogTransKernSD, propToModify = MCMCcontrol$propPhyloBranchesToModify)),
   Lambda = list(
     logPriorFun = .coalescenceRatesLogPriorFun,
     transFun = function(x) .coalescenceRatesTransFun(x, sd = MCMCcontrol$coalRateKernelSD)))
  nodeNumbers <- (ape::Ntip(startingTree) + 1):(length(startingTree$edge.length) + 1)
  nodeTimesLogPriorAndTransFunList <- lapply(nodeNumbers, function(nodeNumber) {
    list(
      logPriorFun = function(x) .nodeTimeConditionalLogPrior(phyloAndTransTree = x, nodeNumber = nodeNumber, meanLogClockRate = log(perSiteClockRate), control = control),
      transFun = function(x) .nodeTimeTransFun(phyloAndTransTree = x, tuningPara = MCMCcontrol$transTreeTuningPara, rootTransitionPar = MCMCcontrol$rootTransitionPar, nodeNumber = nodeNumber))
  })
  timeNames <- paste("t", nodeNumbers, sep = "")
  names(nodeTimesLogPriorAndTransFunList) <- timeNames
  logPriorAndTransFunList <- c(logPriorAndTransFunList, nodeTimesLogPriorAndTransFunList[nodeOrder - ape::Ntip(startingTree)])
  logPriorAndTransFunTopology <- list(
    topology = list(
    logPriorFun = function(x) {
      .topologyLogPriorFun(phyloAndTransTree = x, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean, estRootTime = estRootTime, control = control)
    },
    transFun = .topologyTransFun))
  logPriorAndTransFunList <- c(logPriorAndTransFunList, logPriorAndTransFunTopology) # Topology is processed last since it might change the order in which parameters must be updated (conditioning for node time priors is with respect to children node times).

  if (length(filesToRestore) == 0) { # New chains...
    logPriorValue <- sapply(names(logPriorAndTransFunList), FUN = function(paraName) logPriorAndTransFunList[[paraName]]$logPriorFun(currentState[[1]]$phyloAndTransTree))
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
      # for (MCMCiter in 1:(MCMCcontrol$n + MCMCcontrol$burnin)) {
      for (MCMCiter in 1:MCMCcontrol$nIterPerSweep) {
        # Re-ordering priors for node times
        logPriorAndTransFunList[tail(seq_along(logPriorAndTransFunList), n = ape::Nnode(chainState$phyloAndTransTree)) - 1] <- nodeTimesLogPriorAndTransFunList[chainState$nodeUpdateOrder] # Last element of logPriorAndTransFunList is for the topology, everything before that is for the node times, hence the re-shuffling of those elements.
        for (paraName in names(logPriorAndTransFunList)) {
          proposalValueAndTransKernRatio <- logPriorAndTransFunList[[paraName]]$transFun(chainState$phyloAndTransTree)
          updatedLogPrior <- chainState$logPrior
          updatedDualTree <- proposalValueAndTransKernRatio$value
          updatedLogLik <- chainState$logLik
          if (paraName %in% c("topology", "b")) {
            updatedLogLik <- logLikFun(.convertToPhylo(updatedDualTree), phyDatData, evoParsList)
          }
          updatedLogPrior <- .updateLogPrior(updatedDualTree, chainState$logPrior, paraName = paraName, timeNames = timeNames, logPriorAndTransFunList = logPriorAndTransFunList)

          proposalLogPP <- (updatedLogLik + sum(updatedLogPrior)) * 1/MCMCcontrol$temperatureParFun(chainNumber)
          exponentValue <- proposalLogPP - chainState$logPP
          MHratio <- proposalValueAndTransKernRatio$transKernRatio * exp(exponentValue)
          # cat("Para.:", paraName, "\n", sep = " ")
          # cat("transKernRatio:", proposalValueAndTransKernRatio$transKernRatio, "exponent:", exponentValue, "proposal LogPP:", proposalLogPP, "current LogPP:", chainState$logPP, "\n", sep = " ")
          if (runif(1) <= MHratio) {
            chainState$phyloAndTransTree <- updatedDualTree
            chainState$logLik <- updatedLogLik
            chainState$logPrior <- updatedLogPrior
            chainState$logPP <- proposalLogPP
            chainState$evoParsList <- evoParsList # Needlessly repeated, but does not take too much memory, and allows the chains to be easily resumed from any iteration.
            if (paraName == "topology") {
              chainState$nodeUpdateOrder <- .getOrderedNodeNumsForTimeUpdate(updatedDualTree) - ape::Ntip(updatedDualTree)
            }
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

.updateLogPrior <- function(updatedDualTree, previousLogPriorValues, paraName, timeNames, logPriorAndTransFunList) {
  nodePrefix <- stringr::str_extract(timeNames[[1]], pattern = "[:alpha:]+(?=[:digit:])")
  timeNameRoot <- stringr::str_c(nodePrefix, ape::Ntip(updatedDualTree) + 1, sep = "")
  updatedLogPrior <- previousLogPriorValues
  updatedLogPrior[[paraName]] <- logPriorAndTransFunList[[paraName]]$logPriorFun(updatedDualTree)
  if (paraName == "b") { # Not very efficient: a modification of a given value of b only affects one branch...
    updatedLogPrior[timeNames] <- sapply(timeNames, function(timeName) {
      logPriorAndTransFunList[[timeName]]$logPriorFun(updatedDualTree)
    })
  } else if (paraName %in% setdiff(timeNames, timeNameRoot)) {
    nodeNumber <- as.numeric(stringr::str_extract(paraName, pattern = "(?<=[:alpha:])[:digit:]+"))
    parentNum <- phangorn::Ancestors(updatedDualTree, nodeNumber, "parent")
    parentLabel <- stringr::str_c(nodePrefix, parentNum, sep = "")
    updatedLogPrior[[parentLabel]] <- logPriorAndTransFunList[[parentLabel]]$logPriorFun(updatedDualTree)
  } else if (paraName == "Lambda") {
    updatedLogPrior[["topology"]] <- logPriorAndTransFunList[["topology"]]$logPriorFun(updatedDualTree)
  }
  updatedLogPrior
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
#' @param transTreeBranchLengthsPriorMVratio priors for branch lengths in the transmission tree are gamma-distributed with this shape parameter; at the default value (1), the distribution is exponential
#'
#' @details You should not need to call this function directly: it is merely a convenient and formal way to let users know how the MCMC run can be customised.
#' @return A list of control parameters
#' @examples \dontrun{
#'  MCMC.control(n = 1e6, stepSize = 100, burnin = 1e5)
#' }
#' @export
MCMC.control <- function(n = 2e5, nIterPerSweep = 100, nChains = 1, temperatureParFun = function(x) x, stepSize = 50, burnin = 1e4, folderToSaveIntermediateResults = NULL, chainId = NULL, coalRateKernelSD = 0.5, print.frequency = 10, phyloBranchLengthsLogTransKernSD = log(1.05), propPhyloBranchesToModify = 0.1, transTreeTuningPara = 0.99, rootTransitionPar = 15, transTreeBranchLengthsPriorMVratio = 1.1, logClockRatesPriorSD = log(1.1), logClockRatesKernelSD = 1, rootTimePriorSD = 7) { # transTreeBranchLengthsPriorMVratio must be greater than 1 for the initial optimisation of branch lengths in the transmission tree to work.
  if (temperatureParFun(1) != 1) {
    stop("temperatureParFun(1) must return '1' as it is for the cold chain! Please fix this (or leave the default value) and re-run the code. \n")
  }
  if (nChains == 1) {
    nIterPerSweep <- stepSize
  }
  list(n = n, nIterPerSweep = nIterPerSweep, nChains = nChains, temperatureParFun = temperatureParFun, stepSize = stepSize, burnin = burnin, folderToSaveIntermediateResults = folderToSaveIntermediateResults, chainId = chainId, coalRateKernelSD = coalRateKernelSD, print.frequency = print.frequency, phyloBranchLengthsLogTransKernSD = phyloBranchLengthsLogTransKernSD, propPhyloBranchesToModify = propPhyloBranchesToModify, transTreeTuningPara = transTreeTuningPara, rootTransitionPar = rootTransitionPar, transTreeBranchLengthsPriorMVratio = transTreeBranchLengthsPriorMVratio, logClockRatesKernelSD = logClockRatesKernelSD, logClockRatesPriorSD = logClockRatesPriorSD, rootTimePriorSD = rootTimePriorSD)
}

.topologyTransFun <- function(phyloAndTransTree) {
  counter <- 0
  # Help file says the tree must be bifurcating for rNNI to work.
  phyloAndTransTree <- .resolveMulti(phyloAndTransTree)
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
  proposedMove <- .collapseIntoMulti(proposedMove)

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

.topologyLogPriorFun <- function(phyloAndTransTree, numMigrationsPoissonPriorMean = 1, estRootTime = NULL, control) {
  Lambda <- phyloAndTransTree$Lambda
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
    cumulWithinRatePerRegion <- sapply(regionNames[withinRegionCoalescencePossible], function(regionName) {
      choose(numLineagesPerRegion[[regionName]], 2) * Lambda[[regionName]] * intervalDuration
    })
    names(cumulWithinRatePerRegion) <- regionNames[withinRegionCoalescencePossible]
    if (length(cumulWithinRatePerRegion) == 0) cumulWithinRatePerRegion <- 0
    totalRate <- sum(cumulWithinRatePerRegion)
    cumulativeLogProb <- cumulativeLogProb - totalRate # That takes into account the exponent of exp(), hence the "-".
    if (sortedPhyloStructCodeFrame$type[[rowsConsidered[[1]] + length(rowsConsidered)]] == "C") {
      regionCode <- transmissionTree$node.label[[sortedPhyloStructCodeFrame$vertexNum[[rowsConsidered[[1]] + length(rowsConsidered)]] - numTips]]$region
      cumulativeLogProb <- cumulativeLogProb + log(Lambda[[regionCode]] * intervalDuration)
    }
  }
  if (!is.null(estRootTime)) {
    cumulativeLogProb <- cumulativeLogProb + .rootTimeLogPriorFun(x = phyloAndTransTree$node.label[[1]]$time, meanValue = estRootTime, sd = control$MCMC.control$rootTimePriorSD)
  }
  cumulativeLogProb + .numMigrationsLogPriorFunWithMeanPar(x = length(migrationTimes), meanPar = numMigrationsPoissonPriorMean)
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


.phyloBranchLengthsConditionalLogPriorFun <- function(phyloAndTransTree, perSiteClockRate, numSites) {
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

.coalescenceRatesLogPriorFun <- function(phyloAndTransTree) {
  return(0)
}

.coalescenceRatesTransFun <- function(phyloAndTransTree, sd) {
  currentState <- phyloAndTransTree$Lambda
  logNewState <- rnorm(n = length(currentState), mean = log(currentState), sd = sd)
  names(logNewState) <- names(currentState)
  phyloAndTransTree$Lambda <- exp(logNewState)
  list(value = phyloAndTransTree, transKernRatio = 1) # The Gaussian transition kernel is symmetrical.
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
.rootTimeLogPriorFun <- function(x, meanValue, sd = 1) {
  dnorm(x = x, mean = meanValue, sd = sd, log = TRUE)
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
  numElements <- ape::Ntip(findBayesianClusterResultsList[[1]]$chain[[1]]$phyloAndTransTree)

  seqNames <- sapply(findBayesianClusterResultsList[[1]]$chain[[1]]$phyloAndTransTree$tip.label, "[[", "name")
  regionNames <- sapply(findBayesianClusterResultsList[[1]]$chain[[1]]$phyloAndTransTree$tip.label, "[[", "region")
  regionSeqNames <- seqNames[regionNames == regionLabel]

  getCombinedAdjacency <- function(chain) {
    getClusterIndicesInIter <-  function(iterIndex) {
      seqNames <- sapply(chain[[iterIndex]]$phyloAndTransTree$tip.label, "[[", "name")
      regionNames <- sapply(chain[[iterIndex]]$phyloAndTransTree$tip.label, "[[", "region")
      clusterVec <- getRegionClusters(chain[[iterIndex]]$phyloAndTransTree, regionLabel)
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
  # phyloAndTransTree <- .refineTimeEstimates(phyloAndTransTree = phyloAndTransTree, estRootTime = estRootTime, logClockRatePriorMean = logClockRatePriorMean, control = control)
  phyloAndTransTree$Lambda <- .genStartCoalescenceRates(phyloAndTransTree, estRootTime = estRootTime, control = control)
  phyloAndTransTree
}

.initialiseTree <- function(phyloAndTransTree, timestampsInDays, regionStamps, startLogClockRate, estRootTime, control) {

  nodesToCheck <- .getOrderedNodeNumsForTimeUpdate(phyloAndTransTree)
  phyloAndTransTree$edge.length <- lapply(seq_along(phyloAndTransTree$edge.length), function(edgeIndex) {
    list(phylogeny = phyloAndTransTree$edge.length[[edgeIndex]], transmissionTree = NA, logXi = NA)
  })

  for (nodeNumber in nodesToCheck) {
    childrenNums <- phangorn::Children(phyloAndTransTree, nodeNumber)
    childrenEdgeNums <- match(childrenNums, phyloAndTransTree$edge[ , 2])
    childrenTimes <- sapply(childrenNums, function(childNum) .getVertexLabel(phyloAndTransTree, childNum)$time)
    maximumTime <- min(childrenTimes)
    funToOptimise <- function(rootTime) {
      phyloAndTransTree$node.label[[nodeNumber - ape::Ntip(phyloAndTransTree)]]$time <- rootTime
      phyloAndTransTree$edge.length[childrenEdgeNums] <- lapply(seq_along(childrenEdgeNums), function(childIndex) {
        newEdge <- phyloAndTransTree$edge.length[[childrenEdgeNums[[childIndex]]]]
        newEdge$transmissionTree <- childrenTimes[[childIndex]] - rootTime
        newEdge$logXi <- log(newEdge$phylogeny) - log(newEdge$transmissionTree)
        newEdge
      })
      nodeLogPrior <- .nodeTimeConditionalLogPrior(phyloAndTransTree, nodeNumber, startLogClockRate, control = control)
      if (nodeNumber == ape::Ntip(phyloAndTransTree) + 1) nodeLogPrior <- nodeLogPrior + .rootTimeLogPriorFun(rootTime, meanValue = estRootTime, sd = control$MCMC.control$rootTimePriorSD)
      -nodeLogPrior
    }
    optRootTime <- nloptr::lbfgs(x0 = maximumTime - 1, fn = funToOptimise, upper = maximumTime - 1/86400, control = list(xtol_rel = 1e-4))

    phyloAndTransTree$node.label[[nodeNumber - ape::Ntip(phyloAndTransTree)]]$time <- optRootTime$par + rnorm(1, sd = 1e-15) # A bit of jittering to ensure that we don't get equal branching times...
  }
  .addTransTreeEdgeLengths(phyloAndTransTree)
}

.addTransTreeEdgeLengths <- function(phyloAndTransTree) {
  phyloAndTransTree$edge.length <- lapply(seq_along(phyloAndTransTree$edge.length), function(x) {
    childNum <- phyloAndTransTree$edge[x, 2]
    parentNum <- phyloAndTransTree$edge[x, 1]
    transTreeBranchLength <- .getVertexLabel(phyloAndTransTree, childNum)$time - .getVertexLabel(phyloAndTransTree, parentNum)$time
    logXi <- log(phyloAndTransTree$edge.length[[x]]$phylogeny)  - log(transTreeBranchLength)
    list(phylogeny = phyloAndTransTree$edge.length[[x]]$phylogeny, transmissionTree = transTreeBranchLength, logXi = logXi)
  })
  phyloAndTransTree
}

.refineTimeEstimates <- function(phyloAndTransTree, estRootTime, logClockRatePriorMean, control) {
  cat("Refining starting value for transmission tree...\n")
  nodeChildrenParentParsList <- lapply(seq_along(phyloAndTransTree$node.label), function(nodeIndex) {
    nodeLabel <- nodeIndex + ape::Ntip(phyloAndTransTree)
    parentNum <- NA
    supportBranchLength <- NA
    supportEdgeNum <- NA
    if (nodeLabel > ape::Ntip(phyloAndTransTree) + 1) {
      parentNum <- phangorn::Ancestors(phyloAndTransTree, nodeLabel, "parent")
      supportEdgeNum <- match(nodeLabel, phyloAndTransTree$edge[ , 2])
      supportBranchLogXi <- phyloAndTransTree$edge.length[[supportEdgeNum]]$logXi
    }

    childrenNums <- phangorn::Children(phyloAndTransTree, nodeLabel)
    childrenBranchesLogXis <- sapply(childrenNums, function(vertexNum) {
      matchValue <- match(vertexNum,  phyloAndTransTree$edge[ , 2])
      phyloAndTransTree$edge.length[[matchValue]]$logXi
    })
    childrenEdgeNums <- match(childrenNums, phyloAndTransTree$edge[ , 2])
    list(parent = parentNum, supportEdgeNum = supportEdgeNum, children = childrenNums, supportBranchLogXi = supportBranchLogXi, childrenBranchesLogXis = childrenBranchesLogXis, childrenEdgeNums = childrenEdgeNums)
  })
  envForPhyloAndTransTree <- new.env()
  assign("phyloAndTransTree", value = phyloAndTransTree, envir = envForPhyloAndTransTree)
  for (nodeIndex in seq_along(phyloAndTransTree$node.label)) {
    nodeInfoList <- nodeChildrenParentParsList[[nodeIndex]]
    childrenNodes <- nodeInfoList$children
    childrenTimes <- sapply(childrenNodes, function(childIndex) .getVertexLabel(phylogeny = phyloAndTransTree, vertexNum = childIndex)$time)
    parentTime <- -Inf
    if (nodeIndex != 1) {
      parentTime <- .getVertexLabel(phylogeny = phyloAndTransTree, vertexNum = nodeInfoList$parent)$time
    }
    objectiveFunForOptim <- function(timePoint, envForPhyloAndTransTree) {
      phyloAndTransTree <- envForPhyloAndTransTree$phyloAndTransTree
      phyloAndTransTree$node.label[[nodeIndex]]$time <- timePoint
      for (childNodeIndex in seq_along(nodeInfoList$childrenEdgeNums)) {
        numerator <- phyloAndTransTree$edge.length[[nodeInfoList$childrenEdgeNums[[childNodeIndex]]]]$phylogeny
        denominator <- .getVertexLabel(phyloAndTransTree, nodeInfoList$children[[childNodeIndex]])$time - timePoint
        newLogXi <- log(numerator) - log(denominator)
        phyloAndTransTree$edge.length[[nodeInfoList$childrenEdgeNums[[childNodeIndex]]]]$transmissionTree <- denominator
        phyloAndTransTree$edge.length[[nodeInfoList$childrenEdgeNums[[childNodeIndex]]]]$logXi <- newLogXi
      }

      if (nodeIndex > 1) {
        numerator <- phyloAndTransTree$edge.length[[nodeInfoList$supportEdgeNum]]$phylogeny
        denominator <- timePoint - .getVertexLabel(phyloAndTransTree, nodeInfoList$parent)$time
        newLogXi <- log(numerator) - log(denominator)
        phyloAndTransTree$edge.length[[nodeInfoList$supportEdgeNum]]$logXi <- newLogXi
        phyloAndTransTree$edge.length[[nodeInfoList$supportEdgeNum]]$transmissionTree <- denominator
      }
      envForPhyloAndTransTree$phyloAndTransTree <- phyloAndTransTree
      objective <- -1 * (sapply(ape::Ntip(phyloAndTransTree):length(phyloAndTransTree$edge.length) + 1, FUN = .nodeTimeConditionalLogPrior, phyloAndTransTree = phyloAndTransTree, meanLogClockRate = logClockRatePriorMean, control = control) + .rootTimeLogPriorFun(phyloAndTransTree$node.label[[1]]$time, meanValue = estRootTime, sd = control$MCMC.control$rootTimePriorSD))
      objective
    }
    currentTime <- .getVertexLabel(phyloAndTransTree, nodeIndex + ape::Ntip(phyloAndTransTree))$time
    nloptr::lbfgs(x0 = currentTime, fn = objectiveFunForOptim, lower = parentTime, upper = min(childrenTimes), envForPhyloAndTransTree = envForPhyloAndTransTree, control = list(xtol_rel = 1e-4))
  }
  cat("Done!")
  envForPhyloAndTransTree$phyloAndTransTree
}

getLogNormMMpars <- function(meanValue, varValue) {
  meanValueSq <- meanValue^2
  list(mu = log(meanValueSq/sqrt(varValue + meanValueSq)), sigma = sqrt(log(varValue/meanValueSq + 1)))
}

.nodeTimeConditionalLogPrior <- function(phyloAndTransTree, nodeNumber, meanLogClockRate, control) {
  childrenNums <- phangorn::Children(phyloAndTransTree, nodeNumber)
  childrenLogXis <- sapply(childrenNums, function(childNum) {
    childEdgeNum <- match(childNum, phyloAndTransTree$edge[, 2])
    phyloAndTransTree$edge.length[[childEdgeNum]]$logXi
  })

  sum(sapply(childrenLogXis, .logXiLogPriorFun, meanValue = meanLogClockRate, stdDev = control$MCMC.control$logClockRatesPriorSD))
}

.logXiLogPriorFun <- function(logXi, meanValue, stdDev) {
  dnorm(logXi, mean = meanValue, sd = stdDev, log = TRUE)
}

.nodeTimeTransFun <- function(phyloAndTransTree, nodeNumber, tuningPara = 0.5, rootTransitionPar = 15) {
  numTips <- length(phyloAndTransTree$tip.label)
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
  updatedNodeTime <- updateNodeTime(nodeNumber)

  phyloAndTransTree$node.label[[nodeNumber - numTips]]$time <- updatedNodeTime
  updatedTransTreeEdges <- .deriveTransTreeEdgeLengths(phyloAndTransTree)
  phyloAndTransTree$edge.length <- lapply(seq_along(phyloAndTransTree$edge.length), function(edgeIndex) {
    updatedEdgeLength <- phyloAndTransTree$edge.length[[edgeIndex]]
    updatedEdgeLength$transmissionTree <- updatedTransTreeEdges[[edgeIndex]]
    updatedEdgeLength$logXi <- log(updatedEdgeLength$phylogeny) - log(updatedEdgeLength$transmissionTree)
    updatedEdgeLength
  })

  list(value = phyloAndTransTree, transKernRatio = 1) # The uniform transition kernel is symmetrical.
}

.getOrderedNodeNumsForTimeUpdate <- function(phyloAndTransTree) {
  nodeCounts <- rep(0, length(phyloAndTransTree$edge.length) + 1)
  nodesToCheck <- rep(0, length(phyloAndTransTree$edge.length) + 1)
  nodeOrderIndex <- 1
  nodesToCheck[1:ape::Ntip(phyloAndTransTree)] <- 1:ape::Ntip(phyloAndTransTree)
  for (i in 1:(length(phyloAndTransTree$edge.length))) { # The root node will always be the last one to be checked, and will be added after the loop. It should not be processed, as it doesn't have a parent. The number of nodes to check therefore corresponds to the number of edges.
    nodeToProcess <- nodesToCheck[[i]]
    parentNum <- phangorn::Ancestors(phyloAndTransTree, nodeToProcess, "parent")
    nodeCounts[[parentNum]] <- nodeCounts[[parentNum]] + 1
    numChildren <- length(phangorn::Children(phyloAndTransTree, parentNum))
    if (nodeCounts[[parentNum]] == numChildren) {
      nodesToCheck[[nodeOrderIndex + ape::Ntip(phyloAndTransTree)]] <- parentNum
      nodeOrderIndex <- nodeOrderIndex + 1
    }
  }
  nodesToCheck[[length(nodesToCheck)]] <- ape::Ntip(phyloAndTransTree) + 1
  tail(nodesToCheck, n = ape::Nnode(phyloAndTransTree))
}

.resolveMulti <- function(phyloAndTransTree) {
  resolvedTree <- ape::multi2di(phyloAndTransTree)
  for (edgeIndex in 1:length(resolvedTree$edge.length)) {
    if (length(resolvedTree$edge.length[[edgeIndex]]) == 1) {
      resolvedTree$edge.length[[edgeIndex]] <- list(phylogeny = 0, transmissionTree = 0, logXi = NA)
    }
  }
  resolvedNodesBool <- sapply(resolvedTree$node.label, is.list)
  repeat {
    undefinedNodeSeq <- currentNode <- match(FALSE, resolvedNodesBool) + ape::Ntip(resolvedTree)
    if (is.na(undefinedNodeSeq)) break
    repeat {
      parentNum <- phangorn::Ancestors(resolvedTree, currentNode, "parent")
      parentNode <- resolvedTree$node.label[[parentNum - ape::Ntip(resolvedTree)]]
      if (is.list(parentNode)) {
        for (nodeNum in undefinedNodeSeq) {
          resolvedTree$node.label[[nodeNum - ape::Ntip(resolvedTree)]] <- list(time = parentNode$time, region = parentNode$region)
        }
        resolvedNodesBool[undefinedNodeSeq - ape::Ntip(resolvedTree)] <- TRUE
        break
      } else {
        undefinedNodeSeq <- c(undefinedNodeSeq, parentNum)
        currentNode <- parentNum
      }
    }
  }
  resolvedTree
}

.collapseIntoMulti <- function(phyloAndTransTree, threshold = 1e-8) {
  phylogeny <- phyloAndTransTree
  phylogeny$edge.length <- sapply(phylogeny$edge.length, "[[", "phylogeny")
  collapsedTree <- ape::di2multi(phylogeny, threshold = threshold)
  collapsedTree$edge.length <- lapply(collapsedTree$edge.length, function(edgeLength) list(phylogeny = edgeLength))
  .addTransTreeEdgeLengths(collapsedTree)
}
