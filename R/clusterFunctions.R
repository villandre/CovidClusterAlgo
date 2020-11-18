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
    cat("Generated new chain ID: ", chainId, "\n")
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
  saveData = TRUE,
  resolvePolytomies = FALSE,
  epidemicRootTimePriorVar = 1/86400,
  lengthForNullExtBranchesInPhylo = 0) {
  list(logLikFun = logLikFun, numMigrationsPoissonPriorMean = numMigrationsPoissonPriorMean, MCMC.control = do.call("MCMC.control", MCMC.control), controlForGenStartTransmissionTree = do.call(".controlForGenStartTransmissionTree", controlForGenStartTransmissionTree), transTreeCondOnPhylo = transTreeCondOnPhylo, saveData = saveData, resolvePolytomies = resolvePolytomies, epidemicRootTimePriorVar = epidemicRootTimePriorVar, lengthForNullExtBranchesInPhylo = lengthForNullExtBranchesInPhylo)
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

.addCoalescenceRates <- function(phyloAndTransTree, estRootTime, control) {
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
  phyloAndTransTree$LambdaList <- lapply(seq_along(subtreeRootNodes), function(i) list(rootNodeNum = subtreeRootNodes[[i]])) # Needed for .computeMeanCoalRateSubtree.
  rateBySubtree <- sapply(seq_along(subtreeRootNodes), .computeMeanCoalRateSubtree, phyloAndTransTree = phyloAndTransTree)
  phyloAndTransTree$LambdaList <- lapply(seq_along(rateBySubtree), function(i) list(Lambda = rateBySubtree[[i]], rootNodeNum = subtreeRootNodes[[i]], priorMean = rateBySubtree[[i]]))
  funToOptimise <- function(logLambdaValue, subtreeIndex) {
    phyloAndTransTree$LambdaList[[subtreeIndex]]$Lambda <- exp(logLambdaValue)
    -(.nodeTimesSubtreeLogPriorFun(phyloAndTransTree, subtreeIndex = subtreeIndex, estRootTime = estRootTime, control = control) + .coalescenceRateLogPriorFun(phyloAndTransTree = phyloAndTransTree, index = subtreeIndex, control = control))
  }
  optimLambdas <- sapply(seq_along(rateBySubtree), function(index) {
    optimValue <- nloptr::lbfgs(x0 = log(rateBySubtree[[index]]), fn = funToOptimise, lower = -20, subtreeIndex = index)
    exp(optimValue$par)
  })
  # LogNormPars <- lapply(rateBySubtree, function(rateValue) computeLogNormMMpars(meanValue = rateValue, varValue = (control$MCMC.control$coalescentPriorCoefVar * rateValue)^2))
  phyloAndTransTree$LambdaList <- lapply(seq_along(rateBySubtree), function(i) list(Lambda = optimLambdas[[i]], rootNodeNum = subtreeRootNodes[[i]], priorMean = rateBySubtree[[i]]))
  phyloAndTransTree
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
  if (MCMCcontrol$nChains == 1) {
    MCMCcontrol$numIterPerSweep <- MCMCcontrol$stepSize # This is to ensure that intermediate results are saved each time a step is completed. By default, they are saved after a swap is attempted.
  }
  startIterNum <- 1
  filesToRestore <- NULL
  if (!is.null(MCMCcontrol$folderToSaveIntermediateResults) & !is.null(MCMCcontrol$chainId)) {
    filesToRestore <- list.files(path = MCMCcontrol$folderToSaveIntermediateResults, pattern = paste(MCMCcontrol$chainId, "_atIter", sep = ""), full.names = TRUE)
  }
  phyloObj <- .convertToPhylo(startingTree)
  startLogLikValue <- logLikFun(phyloObj, phyDatData, evoParsList)
  nodeOrder <- .getOrderedNodeNumsForTimeUpdate(startingTree)
  startingState <- lapply(1:MCMCcontrol$nChains, function(x) list(phyloAndTransTree = startingTree, logLik = startLogLikValue, nodeUpdateOrder = nodeOrder - ape::Ntip(startingTree))) # All chains start in the same state.
  MCMCcontainer <- vector(mode = "list", length = floor((MCMCcontrol$n + MCMCcontrol$burnin)/MCMCcontrol$stepSize))
  if (length(filesToRestore) > 0) {
    iterNums <- as.numeric(stringr::str_extract(filesToRestore, "[:digit:]+(?=.Rdata)"))
    filesToRestoreOrder <- order(iterNums)
    MCMCcontainer[1:length(filesToRestore)] <- lapply(filesToRestore[filesToRestoreOrder], function(filename) {
      loadName <- load(filename)
      get(loadName)
    })
    startIterNum <- max(iterNums) + 1
    cat("Resuming MCMC at iteration ", startIterNum, ". \n", sep = "")
    startingState <- MCMCcontainer[[length(filesToRestore)]]
  } else {
    cat("Launching MCMC... \n")
  }
  totalBranchLengthsPhylo <- sum(sapply(startingTree$edge.length, "[[", "phylogeny"))
  branchLengthsLogPriorAndTransFunList <- list(
   b = list( # The mean here corresponds to 14 days for the time between the time the sequence is sampled and the time at which the first transmission linking it to another sequence in the sample occurred.
      # logPriorFun =  function(x) .phyloBranchLengthsLogPriorFun(phyloAndTransTree = x, extBranchLengthPriorMean = perSiteClockRate, extBranchLengthPriorSD = perSiteClockRate * 50), # SHOULD NOT BE HARD-CODED
     logPriorFun = function(x) .phyloBranchLengthsLogPriorFun(phyloAndTransTree = x, extBranchLengthPriorMean = totalBranchLengthsPhylo, extBranchLengthPriorSD = totalBranchLengthsPhylo/100),
      transFun = function(x) .phyloBranchLengthsTransFun(x, phyloBranchLengthsTransFunTuningPara = MCMCcontrol$phyloBranchLengthsTransFunTuningPara, propToModify = MCMCcontrol$propPhyloBranchesToModify)))

  coalescenceRatesLogPriorAndTransFunList <-  lapply(seq_along(startingTree$LambdaList), function(i) {
    list(
      logPriorFun = function(x) .coalescenceRateLogPriorFun(phyloAndTransTree = x, index = i, control = control),
      transFun = function(x) .coalescenceRateTransFun(x, index = i, sd = MCMCcontrol$coalRateKernelSD))
  })
  names(coalescenceRatesLogPriorAndTransFunList) <- paste("Lambda", seq_along(coalescenceRatesLogPriorAndTransFunList), sep = "")

  nodeTimesLogPriorAndTransFunList <- lapply(seq_along(coalescenceRatesLogPriorAndTransFunList), function(subtreeIndex) list(
    logPriorFun = function(x) .nodeTimesSubtreeLogPriorFun(phyloAndTransTree = x, subtreeIndex = subtreeIndex, estRootTime = estRootTime, control = control),
    transFun = function(x) .nodeTimesSubtreeTransFun(phyloAndTransTree = x, subtreeIndex = subtreeIndex, control = control)
  ))
  names(nodeTimesLogPriorAndTransFunList) <- paste("nodeTimes", seq_along(coalescenceRatesLogPriorAndTransFunList), sep = "")
  logPriorAndTransFunList <- c(branchLengthsLogPriorAndTransFunList, coalescenceRatesLogPriorAndTransFunList, nodeTimesLogPriorAndTransFunList)
  if (MCMCcontrol$topologyTransition) {
    stop("Topology transitions not implemented yet. Set topologyTransition to FALSE in MCMC.control. \n")
    # logPriorAndTransFunTopology <- list(
    #   topology = list(
    #     logPriorFun = function(x) {
    #       .topologyLogPriorFun(phyloAndTransTree = x, numMigrationsPoissonPriorMean = control$numMigrationsPoissonPriorMean, estRootTime = estRootTime, control = control)
    #     },
    #     transFun = .topologyTransFun))
    #
    # logPriorAndTransFunList <- c(logPriorAndTransFunList, logPriorAndTransFunTopology) # Topology is processed last since it might change the order in which parameters must be updated (conditioning for node time priors is with respect to children node times).
  }
  if (length(filesToRestore) == 0) { # New chains...
    logPriorValue <- sapply(names(logPriorAndTransFunList), FUN = function(paraName) logPriorAndTransFunList[[paraName]]$logPriorFun(startingState[[1]]$phyloAndTransTree))
    names(logPriorValue) <- names(logPriorAndTransFunList)
    startingState <- lapply(1:MCMCcontrol$nChains, function(chainNum) {
      startingStateMember <- startingState[[chainNum]]
      startingStateMember$logPrior <- logPriorValue
      startingStateMember$logPP <- (startLogLikValue + sum(startingStateMember$logPrior)) * (1/MCMCcontrol$temperatureParFun(chainNum))
      startingStateMember
    })
  }
  cat("Chain is called ", chainId, ". Specify this string in MCMC.control if you want to resume simulations.\n", sep = "")
  cat("Starting values for log-priors: \n")
  print(startingState[[1]]$logPrior) # "1" is the cold chain
  cat("Starting value for log-lik.: \n")
  print(startingState[[1]]$logLik) # "1" is the cold chain
  cat("Starting value for log-PP: \n")
  print(startingState[[1]]$logLik + sum(startingState[[1]]$logPrior))
  clusterAddress <- NULL
  if (MCMCcontrol$nChains > 1) {
    clusterAddress <- parallel::makeForkCluster(nnodes = MCMCcontrol$nChains)
  }
  totalNumSweeps <- ceiling((MCMCcontrol$n + MCMCcontrol$burnin)/MCMCcontrol$nIterPerSweep)
  sweepFun <- function(sweepNum) {
    chainFun <- function(chainNumber = 1) {
      stateContainer <- vector(mode = "list", length = ceiling(MCMCcontrol$nIterPerSweep/MCMCcontrol$stepSize))
      chainState <- startingState[[chainNumber]]
      if (length(MCMCcontainer[[1]]) > 0) chainState <- MCMCcontainer[[tail(which(sapply(MCMCcontainer, length) > 0), n = 1)]][[chainNumber]]
      for (MCMCiter in 1:MCMCcontrol$nIterPerSweep) {
        for (paraName in names(logPriorAndTransFunList)) {
          proposalValueAndTransKernRatio <- logPriorAndTransFunList[[paraName]]$transFun(chainState$phyloAndTransTree)
          updatedLogPrior <- chainState$logPrior
          updatedDualTree <- proposalValueAndTransKernRatio$value
          updatedLogLik <- chainState$logLik
          if (paraName %in% c("topology", "b")) {
            updatedLogLik <- logLikFun(.convertToPhylo(updatedDualTree), phyDatData, evoParsList)
            if (paraName == "b") {
              diffBranchNums <- which(getPhyloEdgeLengths(updatedDualTree) != getPhyloEdgeLengths(chainState$phyloAndTransTree))
              parentNums <- updatedDualTree$edge[diffBranchNums, 1]
              subtreeNums <- unique(sapply(updatedDualTree$node.label[parentNums - ape::Ntip(updatedDualTree)], "[[", "subtreeIndex"))
              for (subtreeNum in subtreeNums) updatedDualTree$LambdaList[[subtreeNum]]$priorMean <- .computeMeanCoalRateSubtree(updatedDualTree, subtreeNum)
            } else {
              stop("Topology transitions not implemented! \n")
            }
          }
          updatedLogPrior <- .updateLogPrior(updatedDualTree, chainState$phyloAndTransTree, chainState$logPrior, paraName = paraName, logPriorAndTransFunList = logPriorAndTransFunList, control = control)

          proposalLogPP <- (updatedLogLik + sum(updatedLogPrior)) * 1/MCMCcontrol$temperatureParFun(chainNumber)
          exponentValue <- proposalLogPP - chainState$logPP
          MHratio <- proposalValueAndTransKernRatio$transKernRatio * exp(exponentValue)
          # cat("Para.:", paraName, "\n", sep = " ")
          # cat("transKernRatio:", proposalValueAndTransKernRatio$transKernRatio, "exponent:", exponentValue, "proposal LogPP:", proposalLogPP, "current LogPP:", chainState$logPP, "\n", sep = " ")
          if (runif(1) <= MHratio) {
            # cat("Jump accepted:", paraName, "\n")
            # cat("Ratio:", MHratio, "\n")
            # cat("Original PP:", chainState$logPP, "\n")
            # cat("New log-PP:", proposalLogPP, "\n")
            # browser()
            chainState$phyloAndTransTree <- updatedDualTree
            chainState$logLik <- updatedLogLik
            chainState$logPrior <- updatedLogPrior
            chainState$logPP <- proposalLogPP
            chainState$evoParsList <- evoParsList # Needlessly repeated, but does not take too much memory, and allows the chains to be easily resumed from any iteration.
            if (paraName == "topology") {
              #NOT IMPLEMENTED
            }
          }
        }
        if ((MCMCiter %% MCMCcontrol$stepSize) == 0) {
          stateContainer[[MCMCiter/MCMCcontrol$stepSize]] <- chainState
        }
      }
      stateContainer
    }
    if (MCMCcontrol$nChains > 1) {
      sweepResultPerChain <- parallel::parLapply(X = 1:MCMCcontrol$nChains, cl = clusterAddress, fun = chainFun)
      # sweepResultPerChain <- lapply(1:MCMCcontrol$nChains, FUN = chainFun) # For testing...
    } else {
      sweepResultPerChain <- list(chainFun())
    }

    if (MCMCcontrol$nChains > 1) {
      # Processing exchanges between chains...
      # Based on https://www.cs.ubc.ca/~nando/540b-2011/projects/8.pdf
      for (k in 1:(MCMCcontrol$nChains - 1)) {
        logExpr <- MCMCcontrol$temperatureParFun(k + 1)/MCMCcontrol$temperatureParFun(k) * tail(sweepResultPerChain[[k + 1]], n = 1)[[1]]$logPP + MCMCcontrol$temperatureParFun(k)/MCMCcontrol$temperatureParFun(k + 1) * tail(sweepResultPerChain[[k]], n = 1)[[1]]$logPP - (tail(sweepResultPerChain[[k]], n = 1)[[1]]$logPP + tail(sweepResultPerChain[[k + 1]], n = 1)[[1]]$logPP)
        swapRatio <- exp(logExpr)
        cat("Swap ratio:", swapRatio, "\n")
        if (runif(1) < swapRatio) {
          cat("Swapping chains", k, "and", k + 1, "\n")
          sweepResultPerChain[c(k, k + 1)] <- sweepResultPerChain[c(k + 1, k)]
          sweepResultPerChain[[k]]$logPP <- sweepResultPerChain[[k]]$logPP * MCMCcontrol$temperatureParFun(k + 1)/MCMCcontrol$temperatureParFun(k)
          sweepResultPerChain[[k]] <- lapply(seq_along(sweepResultPerChain[[k]]), function(index) {
            sweepResultPerChain[[k]][[index]]$logPP <- sweepResultPerChain[[k]][[index]]$logPP * MCMCcontrol$temperatureParFun(k + 1)/MCMCcontrol$temperatureParFun(k)
            sweepResultPerChain[[k]][[index]]
          })

          sweepResultPerChain[[k + 1]] <- lapply(seq_along(sweepResultPerChain[[k + 1]]), function(index) {
            sweepResultPerChain[[k + 1]][[index]]$logPP <- sweepResultPerChain[[k + 1]][[index]]$logPP * MCMCcontrol$temperatureParFun(k)/MCMCcontrol$temperatureParFun(k + 1)
            sweepResultPerChain[[k + 1]][[index]]
          })
        }
      }
    }

    cat("Processed MCMC iteration ", sweepNum * MCMCcontrol$nIterPerSweep, ".\n", sep = "")
    cat("\n\n\n Current log-lik.:", tail(sweepResultPerChain[[1]], n = 1)[[1]]$logLik, "\n")
    cat("Current log-PP:", tail(sweepResultPerChain[[1]], n = 1)[[1]]$logPP, "\n")

    reformattedSweepResult <- lapply(1:ceiling(MCMCcontrol$nIterPerSweep/MCMCcontrol$stepSize), function(iterNumber) {
      lapply(1:length(sweepResultPerChain), function(chainNumber) {
        sweepResultPerChain[[chainNumber]][[iterNumber]]
      })
    })
    elementNums <- (sweepNum - 1) * ceiling(MCMCcontrol$nIterPerSweep/MCMCcontrol$stepSize) + 1:ceiling(MCMCcontrol$nIterPerSweep/MCMCcontrol$stepSize)
    MCMCcontainer[elementNums] <<- reformattedSweepResult
    if (!is.null(MCMCcontrol$folderToSaveIntermediateResults)) {
      lapply(elementNums, function(elementNum) {
        currentState <- MCMCcontainer[[elementNum]]
        save(currentState, file = paste(MCMCcontrol$folderToSaveIntermediateResults, "/chainID_", chainId, "_atIter", elementNum * MCMCcontrol$stepSize, ".Rdata", sep = ""), compress = TRUE)
      })
    }
  }
  currentSweep <- (startIterNum - 1)/MCMCcontrol$nIterPerSweep + 1
  lapply(currentSweep:totalNumSweeps, FUN = sweepFun)
  cat("MCMC complete. Finalising... \n")
  lastIterToDrop <- floor(MCMCcontrol$burnin/MCMCcontrol$stepSize)
  elementsToDrop <- 0
  if (lastIterToDrop > 0) {
    elementsToDrop <- 1:lastIterToDrop
  }
  outputResults <- lapply(MCMCcontainer[-elementsToDrop], '[[', 1) # We only keep the cold chain.
  if (MCMCcontrol$nChains > 1) {
    parallel::stopCluster(clusterAddress)
  }
  outputResults
}

.updateLogPrior <- function(updatedDualTree, originalDualTree, previousLogPriorValues, paraName, logPriorAndTransFunList, control) {
  updatedLogPrior <- previousLogPriorValues
  updatedLogPrior[[paraName]] <- logPriorAndTransFunList[[paraName]]$logPriorFun(updatedDualTree)
  if (paraName == "b") { # Lambda prior means are already modified
    diffBranchNums <- which(getPhyloEdgeLengths(updatedDualTree) != getPhyloEdgeLengths(originalDualTree))
    if (length(diffBranchNums) > 0) {
      parentNums <- updatedDualTree$edge[diffBranchNums, 1]
      subtreeNums <- unique(sapply(updatedDualTree$node.label[parentNums - ape::Ntip(updatedDualTree)], "[[", "subtreeIndex"))
      LambdaNames <- paste("Lambda", subtreeNums, sep = "")
      updatedLogPrior[LambdaNames] <- sapply(LambdaNames, function(LambdaName) logPriorAndTransFunList[[LambdaName]]$logPriorFun(updatedDualTree))
    }
  } else if (paraName %in% grep(pattern = "Lambda", x = names(logPriorAndTransFunList), value = TRUE)) {
    subtreeIndex <- as.numeric(gsub(x = paraName, pattern = "[a-zA-Z]+", replacement = ""))
    priorName <- paste("nodeTimes", subtreeIndex, sep = "")
    updatedLogPrior[[priorName]] <- logPriorAndTransFunList[[priorName]]$logPriorFun(updatedDualTree)
  }
  updatedLogPrior
}

#' Control parameters for the MCMC run.
#'
#' @param n total number of iterations after the burn-in
#' @param stepSize frequency at which iterations are kept in memory, relates to thinning, e.g. stepSize = 100 means that the returned chain will comprise iterations 100, 200, 300, after the burn-in
#' @param burnin number of iterations to discard at the start of the chain
#' @param folderToSaveIntermediateResults folder where intermediate chain results are saved, used to allow chains to be interrupted/resumed
#' @param chainId glyph program previously assigned a chain, used to resume a chain, can be found in the printed output of findBayesianClusters after "Launching MCMC..."
#' @param coalRateKernelSD tuning parameter for the Gaussian coalescence rates transition kernel
#' @param print.frequency number indicating at which interval the program prints information on the chain
#' @param phyloBranchLengthsTransFunTuningPara number between 0 and 1, tuning parameter for the phylogenetic branch lengths transition kernel, a lower value translates to larger variations, and a lower acceptance rate in the chain
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
MCMC.control <- function(n = 2e5, nIterPerSweep = 100, nChains = 1, temperatureParFun = function(x) x, stepSize = 50, burnin = 1e4, folderToSaveIntermediateResults = NULL, chainId = NULL, coalRateKernelSD = 0.5, print.frequency = 10, phyloBranchLengthsTransFunTuningPara = 0.5, propPhyloBranchesToModify = 0.1, transTreeTuningPara = 0.99, rootTransitionPar = 15, logClockRatesPriorSD = log(5), rootTimePriorSD = 7, topologyTransition = FALSE, coalescentPriorCoefVar = 10, propNodeTimesToModify = 0) {
  if (temperatureParFun(1) != 1) {
    stop("temperatureParFun(1) must return '1' as it is for the cold chain! Please fix this (or leave the default value) and re-run the code. \n")
  }
  if (nChains == 1) {
    nIterPerSweep <- stepSize
  }
  list(n = n, nIterPerSweep = nIterPerSweep, nChains = nChains, temperatureParFun = temperatureParFun, stepSize = stepSize, burnin = burnin, folderToSaveIntermediateResults = folderToSaveIntermediateResults, chainId = chainId, coalRateKernelSD = coalRateKernelSD, print.frequency = print.frequency, phyloBranchLengthsTransFunTuningPara = phyloBranchLengthsTransFunTuningPara, propPhyloBranchesToModify = propPhyloBranchesToModify, transTreeTuningPara = transTreeTuningPara, rootTransitionPar = rootTransitionPar, rootTimePriorSD = rootTimePriorSD, logClockRatesPriorSD = logClockRatesPriorSD, topologyTransition = topologyTransition, coalescentPriorCoefVar = coalescentPriorCoefVar, propNodeTimesToModify = propNodeTimesToModify)
}

.topologyTransFun <- function(phyloAndTransTree) {
  transKernRatio <- 1
  proposedMove <- .rNNItransTree(phyloAndTransTree)
  # I don't want a move to be processed if the proposal is identical to the original state.
  if (identical(proposedMove, phyloAndTransTree)) transKernRatio <- 1e-300
  list(value = proposedMove, transKernRatio = transKernRatio)
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

# .topologyLogPriorFun <- function(phyloAndTransTree, numMigrationsPoissonPriorMean = 1, estRootTime = NULL, control) {
#   Lambda <- phyloAndTransTree$Lambda
#   numTips <- length(phyloAndTransTree$tip.label)
#   numNodes <- length(phyloAndTransTree$node.label)
#   transmissionTree <- .convertToTransTree(phyloAndTransTree)
#   nodeTimes <- sapply(transmissionTree$node.label, function(nodeLabel) nodeLabel$time)
#   tipTimes <- sapply(transmissionTree$tip.label, function(tipLabel) tipLabel$time)
#   getMigrationTimesAndChildNodeNums <- function(nodeNum) {
#     nodeIndicator <- nodeNum > numTips
#     indexValue <- nodeNum - nodeIndicator * numTips
#     listToIndex <- "tip.label"
#     if (nodeIndicator) listToIndex <- "node.label"
#     # parentNum <- phangorn::Ancestors(x = transmissionTree, node = nodeNum, type = "parent")
#     parentNum <- transmissionTree$edge[match(nodeNum, transmissionTree$edge[ , 2]), 1]
#     # parentRegion <- .getVertexLabel(transmissionTree, parentNum)$region
#     parentRegion <- transmissionTree$node.label[[parentNum - numTips]]$region
#     currentRegion <- transmissionTree[[listToIndex]][[indexValue]]$region
#     returnValue <- list(childNode = nodeNum, time = NULL)
#     if (parentRegion != currentRegion) returnValue$time <- (transmissionTree$node.label[[parentNum - numTips]]$time + transmissionTree[[listToIndex]][[indexValue]]$time)/2
#     returnValue
#   }
#   migrationTimesAndChildNodeNums <- lapply(seq_along(transmissionTree$node.label), FUN = getMigrationTimesAndChildNodeNums)
#   migrationTimes <- do.call("c", lapply(migrationTimesAndChildNodeNums, function(listElement) listElement$time))
#   migrationChildNodes <- do.call("c", lapply(migrationTimesAndChildNodeNums, function(listElement) {
#     if (!is.null(listElement$time)) return(listElement$childNode)
#     NULL
#   }))
#   verticesInvolved <- c(1:(numTips + numNodes), migrationChildNodes)
#   eventType <- rep(c("T", "C", "M"), c(numTips, numNodes, length(migrationTimes)))
#   eventTimes <- c(tipTimes, nodeTimes, migrationTimes)
#   phyloStructCodeFrame <- data.frame(vertexNum = verticesInvolved, type = eventType, time = eventTimes)
#   orderIndices <- order(eventTimes, decreasing = TRUE) # We start at the tips, hence the decreasing time.
#
#   sortedPhyloStructCodeFrame <- phyloStructCodeFrame[orderIndices, ]
#   regionNames <- unique(sapply(transmissionTree$tip.label, function(x) x$region))
#   numRegions <- length(regionNames)
#   numLineagesPerRegion <- rep(0, numRegions)
#   names(numLineagesPerRegion) <- regionNames
#   cumulativeLogProb <- 0
#   timeIndices <- as.numeric(factor(sortedPhyloStructCodeFrame$time))
#   for (timeIndex in head(unique(timeIndices), n = -1)) {
#     rowsConsidered <- which(timeIndices == timeIndex)
#     for (rowIndex in rowsConsidered) {
#       vertexNumber <- sortedPhyloStructCodeFrame$vertexNum[rowIndex]
#       if (sortedPhyloStructCodeFrame$type[rowIndex] == "T") {
#         numLineagesPerRegion[[transmissionTree$tip.label[[vertexNumber]]$region]] <- numLineagesPerRegion[[transmissionTree$tip.label[[vertexNumber]]$region]] + 1
#       } else if (sortedPhyloStructCodeFrame$type[rowIndex] == "C") {
#         numLineagesPerRegion[[transmissionTree$node.label[[vertexNumber - numTips]]$region]] <- numLineagesPerRegion[[transmissionTree$node.label[[vertexNumber - numTips]]$region]] - 1
#       } else {
#         childNodeNum <- sortedPhyloStructCodeFrame$vertexNum[rowIndex]
#         # parentNodeNum <- phangorn::Ancestors(x = transmissionTree, node = childNodeNum, type = "parent")
#         parentNodeNum <- transmissionTree$edge[match(childNodeNum, transmissionTree$edge[ , 2]), 1]
#         nodeIndicator <- childNodeNum > numTips
#         vertexIndex <- childNodeNum - numTips * nodeIndicator
#         listName <- "tip.label"
#         if (nodeIndicator) listName <- "node.label"
#         childRegion <- transmissionTree[[listName]][[vertexIndex]]$region
#         parentRegion <- transmissionTree$node.label[[parentNodeNum - numTips]]$region
#         # A migration along a branch in reverse time eliminates one lineage in the child region and adds one in the parent region.
#         numLineagesPerRegion[[childRegion]] <- numLineagesPerRegion[[childRegion]] - 1
#         numLineagesPerRegion[[parentRegion]] <- numLineagesPerRegion[[parentRegion]] + 1
#       }
#     }
#     intervalDuration <- sortedPhyloStructCodeFrame$time[[rowsConsidered[[1]]]] - sortedPhyloStructCodeFrame$time[[rowsConsidered[[1]] + length(rowsConsidered)]]
#
#     withinRegionCoalescencePossible <- numLineagesPerRegion > 1
#     cumulWithinRatePerRegion <- sapply(regionNames[withinRegionCoalescencePossible], function(regionName) {
#       choose(numLineagesPerRegion[[regionName]], 2) * Lambda[[regionName]] * intervalDuration
#     })
#     names(cumulWithinRatePerRegion) <- regionNames[withinRegionCoalescencePossible]
#     if (length(cumulWithinRatePerRegion) == 0) cumulWithinRatePerRegion <- 0
#     totalRate <- sum(cumulWithinRatePerRegion)
#     cumulativeLogProb <- cumulativeLogProb - totalRate # That takes into account the exponent of exp(), hence the "-".
#     if (sortedPhyloStructCodeFrame$type[[rowsConsidered[[1]] + length(rowsConsidered)]] == "C") {
#       regionCode <- transmissionTree$node.label[[sortedPhyloStructCodeFrame$vertexNum[[rowsConsidered[[1]] + length(rowsConsidered)]] - numTips]]$region
#       cumulativeLogProb <- cumulativeLogProb + log(Lambda[[regionCode]] * intervalDuration)
#     }
#   }
#   if (!is.null(estRootTime)) {
#     cumulativeLogProb <- cumulativeLogProb + .rootTimeLogPriorFun(x = phyloAndTransTree$node.label[[1]]$time, meanValue = estRootTime, sd = control$MCMC.control$rootTimePriorSD)
#   }
#   cumulativeLogProb + .numMigrationsLogPriorFunWithMeanPar(x = length(migrationTimes), meanPar = numMigrationsPoissonPriorMean)
# }

.topologyLogPriorFun <- function(phyloAndTransTree) {
  return(0)
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
  sum(dnorm(sum(branchLengths), mean = extBranchLengthPriorMean, sd = extBranchLengthPriorSD, log = TRUE))
    # externalEdges <- which(phyloAndTransTree$edge[ , 2] <= ape::Ntip(phyloAndTransTree))
    # paretoPars <- computeParetoMMpars(meanValue = extBranchLengthPriorMean, minValue = 1e-15)
    # sum(dpareto(x = branchLengths[externalEdges], alpha = paretoPars$alpha, xm = paretoPars$xm, log = TRUE))
    # logNormPars <- computeLogNormMMpars(meanValue = extBranchLengthPriorMean, varValue = extBranchLengthPriorSD^2)
    # sum(dlnorm(branchLengths[externalEdges], meanlog = logNormPars$mu, sdlog = logNormPars$sigma, log = TRUE))
    # logMean <- log(extBranchLengthPriorMean)
    # sum(dnorm(log(branchLengths[externalEdges]), mean = logMean, sd = log(5), log = TRUE))
}

# The mode of the displaced log-normal kernel is set at current branch lengths.
# It is then enough to pick a value for sigma, the variance of the random variable expressed on the log-scale, to get a suitable value for mu.
# A better setup might involve working on the logarithmic scale and simulating normal numbers, making the transition kernel symmetric.

.phyloBranchLengthsTransFun <- function(phyloAndTransTree, phyloBranchLengthsTransFunTuningPara = 0.5, propToModify = 0.1) {
  # subtreeIndices <- seq_along(phyloAndTransTree$LambdaList)
  # for (subtreeIndex in subtreeIndices) {
  #   nodesInSubtree <- which(sapply(phyloAndTransTree$node.label, "[[", "subtreeIndex") == subtreeIndex)
  #   nodeChildren <- do.call("c", phyloAndTransTree$childrenList[nodesInSubtree])
  #   matchingBranchIndices <- match(nodeChildren, phyloAndTransTree$edge[ , 2])
  #   branchLengths <- sapply(matchingBranchIndices, function(branchIndex) phyloAndTransTree$edge.length[[branchIndex]]$phylogeny)
  #   # scalingFactor <- runif(1, min = 1/(1 + phyloBranchLengthsTransFunTuningPara), max = 1 + phyloBranchLengthsTransFunTuningPara)
  #   scalingFactor <- sample(c(1/(1 + phyloBranchLengthsTransFunTuningPara), 1 + phyloBranchLengthsTransFunTuningPara), size = 1)
  #   newBranchLengths <- scalingFactor * branchLengths
  #   for (j in seq_along(matchingBranchIndices)) {
  #     phyloAndTransTree$edge.length[[matchingBranchIndices[[j]]]]$phylogeny <- newBranchLengths[[j]]
  #   }
  # }
  # numToModify <- 1
  # if (propToModify > 0) {
  #   numToModify <- ceiling(propToModify * length(branchLengths))
  # }
  # branchesToModify <- sample.int(n = length(branchLengths), size = numToModify, replace = FALSE)
  # previousLengths <- branchLengths[branchesToModify]
  # # logScaleMean <- log(previousLengths)
  # # newLengths <- exp(rnorm(n = length(logScaleMean), mean = logScaleMean, sd = phyloBranchLengthsTransFunTuningPara))
  # multipliers <- runif(length(previousLengths), min = 1 - phyloBranchLengthsTransFunTuningPara, max = 1 + phyloBranchLengthsTransFunTuningPara)
  # newLengths <- multipliers * previousLengths
  # for (i in seq_along(branchesToModify)) {
  #   phyloAndTransTree$edge.length[[branchesToModify[[i]]]]$phylogeny <- newLengths[[i]]
  #   phyloAndTransTree$edge.length[[branchesToModify[[i]]]]$logXi <- log(newLengths[[i]]) - log(phyloAndTransTree$edge.length[[branchesToModify[[i]]]]$transmissionTree)
  # }
  #
  list(value = phyloAndTransTree, transKernRatio = 1)
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

.updateTransTreeEdgeLengths <- function(phyloAndTransTree) {
  updatedPhyloAndTransTree <- phyloAndTransTree
  childNodeNumbers <- updatedPhyloAndTransTree$edge[ , 2]
  parentNodeNumbers <- updatedPhyloAndTransTree$edge[ , 1]
  numTips <- length(updatedPhyloAndTransTree$tip.label)
  newLengths <- sapply(childNodeNumbers, function(x) .getVertexLabel(phylogeny = updatedPhyloAndTransTree, vertexNum = x)$time) - sapply(parentNodeNumbers, function(x) updatedPhyloAndTransTree$node.label[[x - numTips]]$time)

  for (i in seq_along(newLengths)) {
    updatedPhyloAndTransTree$edge.length[[i]]$transmissionTree <- newLengths[[i]]
    updatedPhyloAndTransTree$edge.length[[i]]$logXi <- log(updatedPhyloAndTransTree$edge.length[[i]]$phylogeny) - log(updatedPhyloAndTransTree$edge.length[[i]]$transmissionTree)
  }
  updatedPhyloAndTransTree
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

# Assumes that migration event is immediately above children node along a transition branch.

.computeSubtreeRootNums <- function(phyloAndTransTree) {
  orderedVertexNums <- .getVertexOrderByDepth(phyloAndTransTree)
  orderedVertexNums <- orderedVertexNums[orderedVertexNums > ape::Ntip(phyloAndTransTree)]
  snipPoints <- vector(mode = "numeric")
  for (vertexNum in tail(orderedVertexNums, n = -1)) {
    parentNum <- phangorn::Ancestors(phyloAndTransTree, vertexNum, "parent")
    parentRegion <- phyloAndTransTree$node.label[[parentNum - ape::Ntip(phyloAndTransTree)]]$region
    currentRegion <- phyloAndTransTree$node.label[[vertexNum - ape::Ntip(phyloAndTransTree)]]$region
    if (!identical(parentRegion, currentRegion)) {
      if (parentNum %in% snipPoints) {
        parentPos <- match(parentNum, snipPoints)
        snipPoints[[parentPos]] <- vertexNum
      } else {
        snipPoints <- c(snipPoints, vertexNum)
      }
    }
  }
  c(snipPoints, ape::Ntip(phyloAndTransTree) + 1)
}

# splitPhyloAndTransTreeByRegion <- function(phyloAndTransTree) {
#   orderedVertexNums <- .getVertexOrderByDepth(phyloAndTransTree)
#   snipPoints <- vector(mode = "numeric")
#   for (vertexNum in tail(orderedVertexNums, n = -1)) {
#     parentNum <- phangorn::Ancestors(phyloAndTransTree, vertexNum, "parent")
#     parentRegion <- phyloAndTransTree$node.label[[parentNum]]$region
#     currentRegion <- phyloAndTransTree$node.label[[vertexNum]]$region
#     if (!identical(parentRegion, currentRegion)) {
#       snipPoints <- c(snipPoints, vertexNum)
#     }
#   }
#   snipPointsDepths <- sapply(snipPoints, function(snipPoint) length(phangorn::Ancestors(phyloAndTransTree, snipPoint, "all") + 1))
#   snipPointsOrder <- order(snipPointsDepths, decreasing = TRUE)
#   treeList <- vector('list', length(snipPoints) + 1)
#   snipPointsReordered <- snipPoints[snipPointsOrder]
#   prunedPhyloAndTransTree <- phyloAndTransTree
#
#   for (nodeIndex in seq_along(snipPointsReordered)) {
#     nodeNum <- snipPointsReordered[[nodeIndex]]
#     treeList[[nodeIndex]] <- ape::extract.clade(prunedPhyloAndTransTree, nodeNum)
#     prunedPhyloAndTransTree <- prune.tree(prunedPhyloAndTransTree, nodeNum)
#   }
#   treeList[[length(treeList) + 1]] <- prunedPhyloAndTransTree
#   keepIndices <- sapply(treeList, function(tree) nrow(tree$edge) > 1)
#   # rootNodes identifies the nodes supporting each introduction of the virus into a different region
#   list(treeList = treeList[keepIndices], rootNodeNums = c(snipPointsReordered[keepIndices], ape::Ntip(phyloAndTransTree) + 1))
# }

.computeMeanCoalRateSubtree <- function(phyloAndTransTree, subtreeIndex) {
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
    childrenNums <- phangorn::Children(phyloAndTransTree, vertexNum)
    matchingBranchNums <- match(childrenNums, phyloAndTransTreeCopy$edge[ , 2])
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

.funForLambdaOptim <- function(logLambda, times, nodeOrTip) {
  timesOrder <- order(times, decreasing = TRUE)
  times <- times[timesOrder]
  nodeOrTip <- nodeOrTip[timesOrder]
  Lambda <- exp(logLambda)
  associationVec <- c(tip = 1, node = -1)
  numLineages <- cumsum(associationVec[nodeOrTip])
  logContributions <- sapply(1:(length(nodeOrTip) - 1), function(index) {
    if (numLineages[[index]] == 1) return(0)
    mergeEvent <- nodeOrTip[[index + 1]] == "node"
    timeElapsed <- times[[index]] - times[[index + 1]]
    numPairs <- choose(numLineages[[index]], 2)
    logContribution <- -numPairs * Lambda * timeElapsed
    if (mergeEvent) {
      logContribution <- logContribution + log(Lambda * numPairs)
    }
    logContribution
  })
  sum(logContributions)
}

.computeMeanCoalRate <- function(phyloAndTransTree) {
  phyloAndTransTreeCopy <- phyloAndTransTree
  orderedVertexNums <- .getVertexOrderByDepth(phyloAndTransTree)
  for (vertexNum in tail(orderedVertexNums, n = -1)) {
    parentNodeNum <- phangorn::Ancestors(phyloAndTransTree, vertexNum, "parent")
    parentEstTime <- phyloAndTransTreeCopy$node.label[[parentNodeNum - ape::Ntip(phyloAndTransTree)]]$time
    edgeNum <- match(vertexNum, phyloAndTransTree$edge[ , 2])
    timeIncrement <- phyloAndTransTree$edge.length[[edgeNum]]$phylogeny/exp(phyloAndTransTree$edge.length[[edgeNum]]$logXi)
    if (vertexNum > ape::Ntip(phyloAndTransTreeCopy)) {
      phyloAndTransTreeCopy$node.label[[vertexNum - ape::Ntip(phyloAndTransTree)]]$time <- parentEstTime + timeIncrement
    } else {
      phyloAndTransTreeCopy$tip.label[[vertexNum]]$time <- parentEstTime + timeIncrement
    }
  }
  funToOptimise <- function(Lambda) {
    .regionNodeTimesLogPrior(phyloAndTransTreeCopy, Lambda)
  }
  startLambda <- (max(sapply(phyloAndTransTreeCopy$tip.label, "[[", "time")) - phyloAndTransTreeCopy$node.label[[ape::Ntip(phyloAndTransTreeCopy) + 1]]$time)/ape::Nnode(phyloAndTransTreeCopy)
  optimResult <- nloptr::lbfgs(x0 = startLambda, fn = funToOptimise)
  optimResult$par
}

.getVertexOrderByDepth <- function(phyloAndTransTree, reverse = FALSE) {
  nodeOrder <- rep(NA, length(phyloAndTransTree$edge.length) + 1)
  nodeOrder[[1]] <- currentNode <- ape::Ntip(phyloAndTransTree) + 1
  readIndex <- 1
  writeIndex <- 2

  while (NA %in% nodeOrder) {
    currentNode <- nodeOrder[[readIndex]]
    if (currentNode > ape::Ntip(phyloAndTransTree)) {
      childrenNums <- phyloAndTransTree$edge[which(phyloAndTransTree$edge[ , 1] == currentNode) , 2]
      nodeOrder[seq_along(childrenNums) + writeIndex - 1] <- childrenNums
      writeIndex <- writeIndex + length(childrenNums)
    }
    readIndex <- readIndex + 1
  }
  if (reverse) nodeOrder <- rev(nodeOrder)
  nodeOrder
}

.regionNodeTimesLogPrior <- function(phyloAndTransTree, Lambda) {
  tipTimes <- getTipTimes(phyloAndTransTree)
  nodeTimes <- getNodeTimes(phyloAndTransTree)
  orderVec <- order(c(tipTimes, nodeTimes), decreasing = TRUE)
  orderedTimes <- c(tipTimes, nodeTimes)[orderVec]
  tipOrNodeVec <- rep(c("tip", "node"), c(length(tipTimes), length(nodeTimes)))[orderVec]
  associationVec <- c(tip = 1, node = -1)
  numLineages <- cumsum(tipOrNodeVec[associationVec])
  logContributions <- sapply(1:(length(tipOrNodeVec) - 1), function(index) {
    if (numLineages == 1) return(0)
    mergeEvent <- tipOrNodeVec[[index + 1]] == "node"
    timeElapsed <- orderedTimes[[index + 1]] - orderedTimes[[index]]
    numPairs <- choose(numLineages[[index]], 2)
    logContribution <- -numPairs * Lambda * timeElapsed
    if (mergeEvent) {
      logContribution <- logContribution + log(Lambda * numPairs)
    }
    logContribution
  })
  sum(logContributions)
}

.nodeTimesSubtreeLogPriorFun <- function(phyloAndTransTree, subtreeIndex, estRootTime, control) {
  Lambda <- phyloAndTransTree$LambdaList[[subtreeIndex]]$Lambda
  subtreeRootNum <- phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum
  subtreeIndices <- sapply(phyloAndTransTree$node.label, "[[", "subtreeIndex")
  nodesInSubtreeNums <- which(subtreeIndices == subtreeIndex) + length(phyloAndTransTree$tip.label)
  getTipTimesList <- function(nodeNum) {
    childrenNums <- phyloAndTransTree$childrenList[[nodeNum - length(phyloAndTransTree$tip.label)]]
    # childrenNums <- phyloAndTransTree$edge[phyloAndTransTree$edge[ , 1] == nodeNum , 2]
    childrenSubtreeIndices <- sapply(childrenNums, function(childNum) {
      if (childNum <= length(phyloAndTransTree$tip.label)) {
        return(phyloAndTransTree$tip.label[[childNum]]$subtreeIndex)
      } else {
        return(phyloAndTransTree$node.label[[childNum - length(phyloAndTransTree$tip.label)]]$subtreeIndex)
      }
    })
    childrenTimes <- sapply(childrenNums, function(childNum) {
      if (childNum <= length(phyloAndTransTree$tip.label)) {
        return(phyloAndTransTree$tip.label[[childNum]]$time)
      } else {
        return(phyloAndTransTree$node.label[[childNum - length(phyloAndTransTree$tip.label)]]$time)
      }
    })
    childrenTimes[((childrenNums > ape::Ntip(phyloAndTransTree)) & (childrenSubtreeIndices != subtreeIndex)) | (childrenNums <= length(phyloAndTransTree$tip.label))]
  }
  tipTimesList <- lapply(nodesInSubtreeNums, getTipTimesList) # The function identifies tips of the *subtree* without having to break up the complete tree into separate components. In this case, a tip either corresponds to a tip in the complete tree, or to an internal node belonging to another subtree supported by a parent that belongs to the subtree numbered subtreeIndex, which represents an introduction of the virus into a new region.
  tipTimes <- do.call("c", tipTimesList)
  nodeTimes <- sapply(nodesInSubtreeNums, function(x) phyloAndTransTree$node.label[[x - length(phyloAndTransTree$tip.label)]]$time)
  nodeOrTip <- rep(c("node", "tip"), c(length(nodeTimes), length(tipTimes)))
  timesToConsider <- c(nodeTimes, tipTimes)
  .funForLambdaOptim(logLambda = log(Lambda), times = timesToConsider, nodeOrTip = nodeOrTip)
}

.coalescenceRateLogPriorFun <- function(phyloAndTransTree, index, control) {
  lambdaValue <- phyloAndTransTree$LambdaList[[index]]$Lambda
  priorMean <-  phyloAndTransTree$LambdaList[[index]]$priorMean
  priorSD <- priorMean * control$MCMC.control$coalescentPriorCoefVar
  logNormPars <- computeLogNormMMpars(meanValue = priorMean, varValue = priorSD^2)
  dlnorm(lambdaValue, meanlog = logNormPars$mu, sdlog = logNormPars$sigma, log = TRUE)
}

.coalescenceRateTransFun <- function(phyloAndTransTree, index, sd) {
  currentState <- phyloAndTransTree$LambdaList[[index]]$Lambda
  logNewState <- rnorm(n = 1, mean = log(currentState), sd = sd)
  names(logNewState) <- names(currentState)
  phyloAndTransTree$LambdaList[[index]]$Lambda <- exp(logNewState)
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
  if (!control$MCMC.control$topologyTransition) {
    phyloAndTransTree$childrenList <- lapply(seq_along(phyloAndTransTree$node.label) + ape::Ntip(phyloAndTransTree), function(nodeNum) phangorn::Children(phyloAndTransTree, nodeNum)) # Useful for speeding up the code when the topology does not change.
  }
  phyloAndTransTree <- .identifyNodeRegions(phyloAndTransTree)
  # Initialises transmission tree branch lengths conditional on phylogenetic branch lengths and an equal clock rate across branches
  phyloAndTransTree <- .initialiseTree(phyloAndTransTree, timestampsInDays, regionStamps, startLogClockRate = logClockRatePriorMean, estRootTime, control)
  # phyloAndTransTree <- .refineTimeEstimates(phyloAndTransTree = phyloAndTransTree, estRootTime = estRootTime, logClockRatePriorMean = logClockRatePriorMean, control = control)
  phyloAndTransTree <- .addCoalescenceRates(phyloAndTransTree, estRootTime = estRootTime, control = control)
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
      nodeLogPrior <- .childrenBranchLengthsWeights(phyloAndTransTree, nodeNumber, startLogClockRate, control = control)
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
    supportPhyloBranchLength <- NA
    supportEdgeNum <- NA
    if (nodeLabel > ape::Ntip(phyloAndTransTree) + 1) {
      parentNum <- phangorn::Ancestors(phyloAndTransTree, nodeLabel, "parent")
      supportEdgeNum <- match(nodeLabel, phyloAndTransTree$edge[ , 2])
      supportPhyloBranchLength <- phyloAndTransTree$edge.length[[supportEdgeNum]]$phylogeny
    }

    childrenNums <- phangorn::Children(phyloAndTransTree, nodeLabel)
    childrenPhyloBranchLengths <- sapply(childrenNums, function(vertexNum) {
      matchValue <- match(vertexNum,  phyloAndTransTree$edge[ , 2])
      phyloAndTransTree$edge.length[[matchValue]]$phylogeny
    })
    childrenEdgeNums <- match(childrenNums, phyloAndTransTree$edge[ , 2])
    list(parent = parentNum, supportEdgeNum = supportEdgeNum, children = childrenNums, supportPhyloBranchLength = supportPhyloBranchLength, childrenPhyloBranchLengths = childrenPhyloBranchLengths, childrenEdgeNums = childrenEdgeNums)
  })

  currentLogPrior <- sum(sapply(1:ape::Nnode(phyloAndTransTree) + ape::Ntip(phyloAndTransTree), .childrenBranchLengthsWeights, phyloAndTransTree = phyloAndTransTree, meanLogClockRate = logClockRatePriorMean, control = control))
  repeat {
    for (nodeIndex in 2:length(phyloAndTransTree$node.label)) {
      nodeInfoList <- nodeChildrenParentParsList[[nodeIndex]]
      childrenNodes <- nodeInfoList$children
      childrenTimes <- sapply(childrenNodes, function(childIndex) .getVertexLabel(phylogeny = phyloAndTransTree, vertexNum = childIndex)$time)
      childrenPhyloLengths <- phyloLengths <- nodeChildrenParentParsList[[nodeIndex]]$childrenPhyloBranchLength
      parentTime <- estRootTime - 365 # One year before the estimated root time.
      parentPhyloLength <- NA
      if (nodeIndex > 1) {
        parentTime <- .getVertexLabel(phylogeny = phyloAndTransTree, vertexNum = nodeInfoList$parent)$time
        parentPhyloLength <- nodeChildrenParentParsList[[nodeIndex]]$supportPhyloBranchLength
        phyloLengths <- c(phyloLengths, parentPhyloLength)
      }
      objectiveFunForOptim <- function(timePoint) {
        newLengths <- childrenTimes - timePoint
        if (nodeIndex > 1) {
          parentLength <- timePoint - parentTime
          newLengths <- c(newLengths, parentLength)
        }
        newLogXis <- log(phyloLengths) - log(newLengths)
        objective <- sum(sapply(newLogXis, .logXiLogPriorFun, mean = logClockRatePriorMean, std = control$MCMC.control$logClockRatesPriorSD))
        if (nodeIndex == 1) {
          objective <- objective + .rootTimeLogPriorFun(timePoint, meanValue = estRootTime, sd = control$MCMC.control$rootTimePriorSD)
        }
        -objective
      }
      currentTime <- .getVertexLabel(phyloAndTransTree, nodeIndex + ape::Ntip(phyloAndTransTree))$time
      optimResult <- nloptr::lbfgs(x0 = currentTime, fn = objectiveFunForOptim, lower = parentTime + 1/86400, upper = min(childrenTimes) - 1/86400, control = list(xtol_rel = 1e-6))
      phyloAndTransTree$node.label[[nodeIndex]]$time <- optimResult$par
      for (childEdgeIndex in seq_along(nodeChildrenParentParsList[[nodeIndex]]$childrenEdgeNums)) {
        childEdgeNum <- nodeChildrenParentParsList[[nodeIndex]]$childrenEdgeNums[[childEdgeIndex]]
        phyloAndTransTree$edge.length[[childEdgeNum]]$transmissionTree <- childrenTimes[[childEdgeIndex]] - optimResult$par
        phyloAndTransTree$edge.length[[childEdgeNum]]$logXi <-  log(childrenPhyloLengths[[childEdgeIndex]]) - log(phyloAndTransTree$edge.length[[childEdgeNum]]$transmissionTree)
      }
      if (nodeIndex > 1) {
        parentEdgeNum <- nodeChildrenParentParsList[[nodeIndex]]$supportEdgeNum
        phyloAndTransTree$edge.length[[parentEdgeNum]]$transmissionTree <- optimResult$par - parentTime
        phyloAndTransTree$edge.length[[parentEdgeNum]]$logXi <- log(parentPhyloLength) - log(phyloAndTransTree$edge.length[[parentEdgeNum]]$transmissionTree)
      }
    }
    newLogPrior <- sum(sapply(1:ape::Nnode(phyloAndTransTree) + ape::Ntip(phyloAndTransTree), .childrenBranchLengthsWeights, phyloAndTransTree = phyloAndTransTree, meanLogClockRate = logClockRatePriorMean, control = control))
    if ((improvementRatio <- abs((newLogPrior - currentLogPrior)/currentLogPrior)) < 1e-8) break
    cat("Improvement ratio:", improvementRatio, "\n")
    currentLogPrior <- newLogPrior
  }
  cat("Done!")
  phyloAndTransTree
}

computeLogNormMMpars <- function(meanValue, varValue) {
  meanValueSq <- meanValue^2
  list(mu = log(meanValueSq/sqrt(varValue + meanValueSq)), sigma = sqrt(log(varValue/meanValueSq + 1)))
}

computeParetoMMpars <- function(meanValue, minValue) {
  list(alpha = -meanValue/(minValue - meanValue), xm = minValue)
}

dpareto <- function(x, alpha, xm, log = TRUE) {
  alpha * xm/(x^(alpha + 1))
}

.childrenBranchLengthsWeights <- function(phyloAndTransTree, nodeNumber, meanLogClockRate, control) {
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

.logXiTransFun <- function(phyloAndTransTree, control) {
  if (control$fixedClockRate) {

  }
}

.nodeTimesSubtreeTransFun <- function(phyloAndTransTree, subtreeIndex, control) {
  if (control$MCMC.control$propNodeTimesToModify == 0) {
    control$MCMC.control$propNodeTimesToModify <- 1e-200
  } # This will ensure that only one time is modified at a time
  nodeNumsToSampleFrom <- which(sapply(phyloAndTransTree$node.label, "[[", "subtreeIndex") == subtreeIndex) + ape::Ntip(phyloAndTransTree)
  nodesToModNums <- nodeNumsToSampleFrom[[1]]
  if (length(nodeNumsToSampleFrom) > 1) {
    nodesToModNums <- sample(nodeNumsToSampleFrom, size = ceiling(control$MCMC.control$propNodeTimesToModify * length(nodeNumsToSampleFrom)))
  }
  tuningPara <- control$MCMC.control$transTreeTuningPara
  rootTransitionPar <- control$MCMC.control$rootTransitionPar

  for (nodeNumber in nodesToModNums) {
    numTips <- length(phyloAndTransTree$tip.label)
    parentNum <- phyloAndTransTree$edge[match(nodeNumber, phyloAndTransTree$edge[ , 2]), 1]
    parentTime <- phyloAndTransTree$node.label[[parentNum - numTips]]$time
    currentTime <- phyloAndTransTree$node.label[[nodeNumber - numTips]]$time
    childrenNums <- phyloAndTransTree$childrenList[[nodeNumber - numTips]]
    childrenTimes <- sapply(childrenNums, function(childIndex) .getVertexLabel(phyloAndTransTree, childIndex)$time)
    branchIndex <- match(nodeNumber, phyloAndTransTree$edge[ , 2])
    if (!is.na(branchIndex)) {
      lowerBound <- currentTime - (currentTime - parentTime) * tuningPara
    } else {
      lowerBound <- currentTime - runif(1, 0, rootTransitionPar) # Time is in days.
    }
    upperBound <- currentTime + abs(min(childrenTimes) - currentTime) * tuningPara
    updatedNodeTime <- runif(1, lowerBound, upperBound)

    phyloAndTransTree$node.label[[nodeNumber - numTips]]$time <- updatedNodeTime
    if (!is.na(branchIndex)) { # Root node has no parent, hence the check.
      phyloAndTransTree$edge.length[[branchIndex]]$transmissionTree <- phyloAndTransTree$edge.length[[branchIndex]]$transmissionTree + (updatedNodeTime - currentTime)
      phyloAndTransTree$edge.length[[branchIndex]]$logXi <- log(phyloAndTransTree$edge.length[[branchIndex]]$phylogeny) - log(phyloAndTransTree$edge.length[[branchIndex]]$transmissionTree)
    }

    childrenIndices <- match(childrenNums, phyloAndTransTree$edge[ , 2])
    for (index in seq_along(childrenIndices)) {
      phyloAndTransTree$edge.length[[childrenIndices[[index]]]]$transmissionTree <- childrenTimes[[index]] - updatedNodeTime
      phyloAndTransTree$edge.length[[childrenIndices[[index]]]]$logXi <- log(phyloAndTransTree$edge.length[[childrenIndices[[index]]]]$phylogeny) - log(phyloAndTransTree$edge.length[[childrenIndices[[index]]]]$transmissionTree)
    }
  }

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
  bifurcatingTree <- ape::multi2di(phyloAndTransTree)
  newEdgeNums <- which(!sapply(bifurcatingTree$edge.length, is.list))
  bifurcatingTree$edge.length[newEdgeNums] <- lapply(seq_along(newEdgeNums), function(x) list(phylogeny = 0, transmissionTree = 0, logXi = NA))
  unresolvedNodesInOrder <- .getNodeUpdateOrderForMulti(bifurcatingTree)

  for (nodeNum in unresolvedNodesInOrder) {
    parentNum <- phangorn::Ancestors(bifurcatingTree, nodeNum, "parent")
    bifurcatingTree$node.label[[nodeNum - ape::Ntip(bifurcatingTree)]] <- bifurcatingTree$node.label[[parentNum - ape::Ntip(bifurcatingTree)]]
  }

  bifurcatingTree
}

.getNodeUpdateOrderForMulti <- function(tree) {
  unresolvedNodes <- which(!sapply(tree$node.label, is.list)) + ape::Ntip(tree)
  resolvedInfoMat <- sapply(unresolvedNodes, function(nodeNum) {
    parentNum <- phangorn::Ancestors(tree, nodeNum, "parent")
    c(parentNum = parentNum, resolved = is.list(tree$node.label[[parentNum - ape::Ntip(tree)]]), processed = FALSE)
  })
  nodesToResolve <- rep(0, length(unresolvedNodes))
  for (i in seq_along(unresolvedNodes)) {
    colIndex <- match(TRUE, resolvedInfoMat["resolved", ] & !resolvedInfoMat["processed", ])
    nodeNum <- unresolvedNodes[[colIndex]]
    nodesToResolve[[i]] <- nodeNum
    resolvedInfoMat["processed", colIndex] <- TRUE
    childrenNodes <- phangorn::Children(tree, nodeNum)
    unresolvedChildrenNums <- intersect(childrenNodes, unresolvedNodes)
    if (length(unresolvedChildrenNums) > 0) {
      matchingColNums <- match(unresolvedChildrenNums, unresolvedNodes)
      resolvedInfoMat["resolved", matchingColNums] <- TRUE
    }
  }
  nodesToResolve
}

.collapseIntoMulti <- function(phyloAndTransTree, threshold = 1e-300) {
  phylogeny <- phyloAndTransTree
  phylogeny$edge.length <- sapply(phylogeny$edge.length, "[[", "phylogeny")
  collapsedTree <- ape::di2multi(phylogeny, tol = threshold)
  collapsedTree$edge.length <- lapply(collapsedTree$edge.length, function(edgeLength) list(phylogeny = edgeLength))
  .addTransTreeEdgeLengths(collapsedTree)
}

# For a NNI move to be possible with respect to the transmission tree, the sibling of the child node of the internal branch selected for the move be assigned a time point which is higher (lower in the tree).
# Function should work with multifurcations.
.rNNItransTree <- function(phyloAndTransTree, moves = 1) {
  childNodesToTry <- sample(seq(from = ape::Ntip(phyloAndTransTree) + 2, to = length(phyloAndTransTree$edge.length) + 1), size = ape::Nnode(phyloAndTransTree) - 1)
  index <- 1
  while (index <= length(childNodesToTry)) {
    childNum <- childNodesToTry[[index]]
    nodeTime <- phyloAndTransTree$node.label[[childNum - ape::Ntip(phyloAndTransTree)]]$time
    nodeSiblings <- phangorn::Siblings(phyloAndTransTree, childNum)
    nodeSiblingsTryOrder <- nodeSiblings
    if (length(nodeSiblingsTryOrder) > 1) {
      nodeSiblingsTryOrder <- sample(nodeSiblings, size = length(nodeSiblings)) # 'sample' is called in case we have a multifurcation
    }

    for (siblingIndex in seq_along(nodeSiblingsTryOrder)) {
      nodeSibling <- nodeSiblingsTryOrder[[siblingIndex]]
      siblingTime <- .getVertexLabel(phyloAndTransTree, nodeSibling)$time
      if ((siblingTime > nodeTime)) break
    }
    if ((siblingTime > nodeTime)) break
    index <- index + 1
  }
  treeToReturn <- phyloAndTransTree
  if (index > length(childNodesToTry)) {
    warning("Could not find a suitable NNI move for the transmission tree! \n")
  } else {
    siblingEdgeNum <- match(nodeSibling, phyloAndTransTree$edge[ , 2])
    grandChildrenEdgeNums <- which(phyloAndTransTree$edge[ , 1] == childNum)
    edgeNumForInterchange <- sample(grandChildrenEdgeNums, size = 1)
    treeToReturn$edge[c(siblingEdgeNum, edgeNumForInterchange), 2] <- treeToReturn$edge[c(edgeNumForInterchange, siblingEdgeNum), 2]
    # edge.length gives supporting branch lengths for nodes listed in the second column of edge. It follows that exchanging two elements in that column must result in a similar exchange of the elements of edge.length.
    treeToReturn$edge.length[c(siblingEdgeNum, edgeNumForInterchange)] <- treeToReturn$edge.length[c(edgeNumForInterchange, siblingEdgeNum)]
  }
  treeToReturn <- .updateTransTreeEdgeLengths(treeToReturn)
  treeToReturn <- .clearNodeRegions(treeToReturn)
  treeToReturn <- .identifyNodeRegions(treeToReturn)
  treeToReturn
}

logPPfromChain <- function(result) {
  sapply(result$chain, "[[", "logPP")
}

logLikfromChain <- function(result) {
  sapply(result$chain, "[[", "logLik")
}

logPriorsFromChain <- function(result) {
  sapply(result$chain, function(chainIter) {
    nodeTimePriors <- stringr::str_detect(names(chainIter$logPrior), pattern = "(?<=[:alpha:])[:digit:]+")
    nodeTimesPriorsSum <- sum(chainIter$logPrior[nodeTimePriors])
    c(chainIter$logPrior[!nodeTimePriors], nodeTimes = nodeTimesPriorsSum)
  })
}

getTipTimes <- function(phyloAndTransTree) {
  sapply(phyloAndTransTree$tip.label, "[[", "time")
}

getLambdas <- function(phyloAndTransTree) {
  sapply(phyloAndTransTree$LambdaList, "[[", "Lambda")
}

.getSubtreeRootNums <- function(phyloAndTransTree) {
  sapply(phyloAndTransTree$LambdaList, "[[", "rootNodeNum")
}

getPhyloEdgeLengths <- function(phyloAndTransTree) {
  sapply(phyloAndTransTree$edge.length, "[[", "phylogeny")
}

getNodeTimes <- function(phyloAndTransTree) {
  sapply(phyloAndTransTree$node.label, "[[", "time")
}

.optimiseExtBranchLengths <- function(phyloAndTransTree, logPriorAndTransFunList) {

}
