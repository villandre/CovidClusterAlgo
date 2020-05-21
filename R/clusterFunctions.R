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


findBayesianClusters <- function(DNAbinData, covariateFrame, branchLengthPriorsFun = NULL, evoParsPriorsFun = NULL, clusterScoringFun = NULL, control = list()) {
  control <- do.call('findBayesianClustersControl', control)

  startingValues <- control$startingValuesFun(DNAbinData)

  sampledTreesWithPP <- .optimTreeMCMC(
    startPhylo = startingValues$phylogeny,
    logLikModelParsList = startingValues$logLikModelParsList,
    DNAbinData = DNAbinData,
    covariateFrame = covariateFrame,
    treeTopologyLogPriorFun = control$treeTopologyLogPriorFun,
    branchLengthsLogPriorsFun = control$branchLengthLogPriorsFun,
    evoParsLogPriorsList = control$evoParsLogPriorsList,
    topologyTransFun = control$topologyTransFun,
    branchLengthsTransFun = control$branchLengthsTransFun,
    fixEvoPars = control$fixEvoPars,
    evoParsTransFunList = control$evoParsTransFunList,
    control = control$MCMC.control)
}

findBayesianClustersControl <- function(
  logLikFun = presetPML, logLikFun.control = list(),
  startingValuesFun = NJandOptim.pml, startingValuesFun.control = list(),
  timedTreeFun = presetTreedater, timedTreeFun.control = list(), fixEvoPars = FALSE) {

  if (identical(startingValuesFun, phangorn::optim.pml) & identical(startingValuesFun.control, list())) {
    startingValuesFun.control <- list(optNni = TRUE, optBf = TRUE, optQ = TRUE, optInv = TRUE, optGamma = TRUE, optEdge = TRUE, optRate = TRUE)
  }

  startingValuesFun <- function(DNAbinData) {
    do.call('startingValuesFun', c(list(DNAbinData), startingValuesFun.control))
  }

  timedTreeFun <- function(phylogeny) {
    do.call('timedTreeFun', c(list(phylogeny), timedTreeFun.control))
  }

  list(logLikFun = logLikFun, startingValuesFun = startingValuesFun, timedTreeFun = timedTreeFun, fixEvoPars = fixEvoPars)
}

# logLikModelParsList is a list of varying parameters for the log-likelihood model. Control parameters are hard-coded.
presetPML <- function(phyloObj, DNAbinObj, logLikModelParsList) {
  phangorn::pml(tree = phyloObj, data = DNAbinObj,
                bf = logLikModelParsList$bf,
                Q = logLikModelParsList$Q,
                inv = logLikModelParsList$propInv,
                k = 4,
                shape = logLikModelParsList$gammaShape,
                model = "GTR")$logLik
}

presetTreedater <- function(phyloObj, tipDates, DNAbinData) {
  treedater::dater(tre = phyloObj, sts = tipDates, s = ncol(DNAbinData), clock = "uncorrelated", ncpu = parallel::detectCores(), temporalConstraints = FALSE, searchRoot = 20)
}

NJandOptim.pml <- function(DNAbinData) {
  distMatrix <- ape::dist.dna(x = DNAbinData, gamma = TRUE, pairwise.deletion = TRUE)
  startingPhylo <- ape::bionj(distMatrix)
  startPML <- phangorn::pml(tree = startingPhylo, k = 4, model = "GTR")
  optimisedPhylo <- phangorn::optim.pml(object = startPML, optNni = TRUE, optBf = TRUE, optQ = TRUE, optInv = TRUE, optGamma = TRUE, optEdge = TRUE, optRate = TRUE)
  list(phylogeny = optimisedPhylo$tree, logLikModelParsList = list(Q = optimisedPhylo$Q, bf = optimisedPhylo$bf, gammaShape = optimisedPhylo$gammaShape, propInv = optimisedPhylo$propInv)) # Fix this.
}

.optimTreeMCMC <- function(
  startPhylo,
  logLikFun,
  logLikModelParsStartList,
  DNAbinData,
  covariateFrame,
  treeTopologyLogPriorFun,
  branchLengthsLogPriorsFun,
  evoParsLogPriorsFun,
  topologyTransFun,
  branchLengthsTransFun,
  fixEvoPars,
  evoParsTransFunList,
  control) {

  logPPfun <- function(phylogeny, evoParsList) {
    logLikArgs <- logLikFun(phylogeny, DNAbinData, evoParsList)
    logPrior <- treeTopologyLogPriorFun(phylogeny, covariateFrame) + branchLengthsLogPriorsFun(log(phylogeny$edge.length))
    if (!fixEvoPars) {
      logPrior <- logPrior + evoParsLogPriorsFun(evoParsList)
    }
    logLik + logPrior
  }
  startLogPP <- logPPfun(phylogeny = startPhylo, evoParsList = logLikModelParsStartList)
  MCMCcontainer <- vector(mode = list(), length = control$n + control$burnin)

  stateList <- c(phylogeny = startPhylo, logPP = startLogPP, evoParsList = logLikModelParsStartList)

  set.seed(control$seed)
  for (MCMCiter in 1:(control$n + control$burnin)) {
    # Handling branch lengths
    for (branchIndex in seq_along(stateList$phylogeny$edge.length)) {
      proposalPhylo <- stateList$phylogeny
      proposalPhylo$edge.length[[branchIndex]] <- branchLengthsTransFun(stateList$phylogeny$edge.length[[branchIndex]])
      newLogPP <- logPPfun(phylogeny = proposalPhylo, evoParsList = stateList$evoParsList)
      acceptanceRatio <- newLogPP/stateList$logPP # Assuming an invertible kernel.
      if (runif(1) < acceptanceRatio) {
        stateList$logPP <- newLogPP
        stateList$phylogeny <- proposalPhylo
      }
    }

    # Handling the topology
    proposalPhylo <- topologyTransFun(stateList$phylogeny)
    newLogPP <- logPPfun(phylogeny = proposalPhylo, evoParsList = stateList$evoParsList)
    if (runif(1) < acceptanceRatio) {
      stateList$logPP <- newLogPP
      stateList$phylogeny <- proposalPhylo
    }

    # Handling the evolutionary parameters
    if (!fixEvoPars) {
      for (evoParIndex in seq_along(stateList$evoParsList)) {
        proposalEvoPars <- stateList$evoParsList
        proposalEvoPars[[evoParIndex]] <- evoParsTransFunList[[evoParIndex]](stateList$evoParsList[[evoParIndex]])
        newLogPP <- logPPfun(phylogeny = stateList$phylogeny, evoParsList = proposalEvoPars)
        if (runif(1) < acceptanceRatio) {
          stateList$logPP <- newLogPP
          stateList$evoParsList <- proposalEvoPars
        }
      }
    }

    if (!is.null(control$folderToSaveIntermediateResults)) {
      filename <- paste(control$folderToSaveIntermediateResults, "/MCMCrunIter", MCMCiter, ".Rdata", sep = "")
      save(stateList, file = filename, compress = TRUE)
    }
    MCMCcontainer[[MCMCiter]] <- stateList
  }
  stepSize <- ceiling(control$n * control$thinning)
  itersToKeep <- seq(from = control$burnin + stepSize + 1, to = control$burnin + control$n, by = stepSize)
  MCMCcontainer[itersToKeep]
}

MCMC.control <- function(n = 1e6, thinning = 0.1, burnin = 1e4, seed = 24, folderToSaveIntermediateResults = NULL) {
  list(n = n, thinning = thinning, burnin = burnin, seed = seed, folderForIntermediateResults = folderForIntermediateResults)
}

defaultClusterScoringFun <- function(x, genotypingData, phylogeny, covariateFrame) {

}
