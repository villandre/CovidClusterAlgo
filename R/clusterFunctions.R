findBayesianClusters <- function(DNAbinData, covariateMatrix, control = list()) {
  control <- do.call('findBayesianClustersControl', control)

  originalDistEstimate <- control$distFun(DNAbinData)
  startPhylo <- control$treeConstructionFun(originalDistEstimate)
  objForOptim <- control$pmlFun(startPhylo)
  if (control$optimStartingTree) {
    startPhylo <- control$optimPhyloFun(objForOptim)
  }
  sampledTreesWithPP <- .optimTreeMCMC(
    startPhylo = startPhylo,
    DNAbinData = DNAbinData,
    treeMoveFun = control$treeMoveFun,
    treeTopologyPriorFun = control$treeTopologyPriorFun,
    branchLengthPriorsFun = control$branchLengthPriorsFun,
    evoParsPriorsFun = control$evoParsPriorFun,
    fixEvoPars = control$fixEvoPars)

}

pmlOpts <- function(k = 4, model = "K80", ...) {
  list(k = k, model = model, ...)
}

optim.pmlOpts <- function(optNni = TRUE, optBf = TRUE, optQ = TRUE, optInv = TRUE, optGamma = TRUE, optEdge = TRUE, optRate = TRUE, ...) {
  list(optNni = optNni, optBf = optBf, optQ = optQ, optInv = optInv, optGamma = optGamma, optEdge = optEdge, optRate = optRate, ...)
}

findBayesianClustersControl <- function(
  treeConstructionFun = ape::bionj, treeConstructionFunControl = list(),
  distFun = ape::dist.dna, distFunControl = list(),
  pmlFun = phangorn::pml, pml.control = list(),
  optim.pmlFun = phangorn::optim.pml, optim.pml.control = list(),
  clusterFun = NULL, clusterFunControl = list(),
  timedTreeFun = treedater::dater, timedTreeFunControl = list(),
  branchLengthPriorsFun = .branchLengthPriors,
  evoParsPriorsFun = .evoParsPriors,
  optimStartingTree = TRUE, fixEvoPars = FALSE) {
  pml.control <- do.call('pmlOpts', pml.control)
  optim.pml.control <- do.call('optim.pmlOpts', optim.pml.control)
  distFun <- function(DNAbinObject) {
    listForDistFunCall <- list(DNAbinObject)
    names(listForDistFunCall) <- names(formals(distFun))[[1]] # Assumes that the DNAbin data argument is the first, which is normally the case, but not always... Should stress it in the help file.
    distFunArguments <- c(listForDistFunCall, distFunControl)
    do.call('distFun', distFunArguments)
  }

  treeConstructionFun <- function(distMatrix) {
    listForTreeConstructionCall <- list(x = distMatrix)
    names(listForTreeConstructionCall) <- names(formals(treeConstructionFun))[[1]] # Assumes that the DNAbin data argument is the first, which is normally the case, but not always...
    treeConstructionArguments <- c(listForTreeConstructionCall, treeConstructionFunControl)
    do.call('treeConstructionFun', treeConstructionArguments)
  }

  pmlFun <- function(phyloObj) {
    do.call('pmlFun', c(list(phyloObj), pml.control))
  }

  optim.pmlFun < function(objForOptim) {
    do.call('optim.pmlFun', c(list(objForOptim), optim.pml.control))
  }

  clusterFun <- function(phyloObject, covariateMatrix) {
    do.call('clusterFun', c(list(phyloObject, covariateMatrix), clusterFunControl))
  }
  list(treeConstructionFun = treeConstructionFun, distFun = distFun, pmlFun = pmlFun, optim.pmlFun = optim.pmlFun, clusterFun = clusterFun, optimiseStartingTree = optimStartingTree)
}

.optimTreeMCMC <- function(startPhylo, DNAbinData,  treeMoveFun, treeTopologyPriorFun, branchLengthPriorsFun, evoParsPriorsFun, fixEvoPars) {

}
