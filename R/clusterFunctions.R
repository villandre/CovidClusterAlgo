# This is a constructor-like function for phylo objects.
phylo <- function(edge, edge.length, tip.label, node.label = NULL) {
  phyloObject <- list(edge = edge, tip.label = tip.label, edge.length = edge.length, Nnode = nrow(edge) - length(tip.label) + 1, node.label = node.label)
  class(phyloObject) <- "phylo"
  phyloObject
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
  # orderedVertexNums <- phyloAndTransTree$vertexOrderByDepth
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
  numTips <- length(phyloAndTransTree$tip.label)
  Lambda <- phyloAndTransTree$LambdaList[[subtreeIndex]]$Lambda

  # getTipTimesList <- function(nodeNum) {
  #   # childrenNums <- phyloAndTransTree$edge[which(nodeNum == phyloAndTransTree$edge[ , 1]), 2]
  #   childrenNums <- phyloAndTransTree$childrenNumList[[nodeNum]]
  #   childrenSubtreeIndices <- sapply(childrenNums, function(childNum) {
  #     if (childNum <= length(phyloAndTransTree$tip.label)) {
  #       return(phyloAndTransTree$tip.label[[childNum]]$subtreeIndex)
  #     } else {
  #       return(phyloAndTransTree$node.label[[childNum - length(phyloAndTransTree$tip.label)]]$subtreeIndex)
  #     }
  #   })
  #   childrenTimes <- sapply(childrenNums, function(childNum) {
  #     if (childNum <= length(phyloAndTransTree$tip.label)) {
  #       return(phyloAndTransTree$tip.label[[childNum]]$time)
  #     } else {
  #       return(phyloAndTransTree$node.label[[childNum - length(phyloAndTransTree$tip.label)]]$time)
  #     }
  #   })
  #   childrenTimes[((childrenNums > ape::Ntip(phyloAndTransTree)) & (childrenSubtreeIndices != subtreeIndex)) | (childrenNums <= length(phyloAndTransTree$tip.label))]
  # }
  # tipTimesList <- lapply(nodesInSubtreeNums, getTipTimesList) # The function identifies tips of the *subtree* without having to break up the complete tree into separate components. In this case, a tip either corresponds to a tip in the complete tree, or to an internal node belonging to another subtree supported by a parent that belongs to the subtree numbered subtreeIndex, which represents an introduction of the virus into a new region.
  # tipTimes <- do.call("c", tipTimesList)
  tipTimes <- sapply(phyloAndTransTree$tipNumbersBySubtree[[subtreeIndex]], function(vertexNum) {
    returnValue <- NULL
    if (vertexNum <= numTips) {
      returnValue <- phyloAndTransTree$tip.label[[vertexNum]]$time
    } else {
      returnValue <- phyloAndTransTree$node.label[[vertexNum - numTips]]$time
    }
    returnValue
  })
  nodeTimes <- sapply(phyloAndTransTree$nodesInSubtreeNums[[subtreeIndex]], function(x) phyloAndTransTree$node.label[[x - numTips]]$time)
  nodeOrTip <- rep(c("node", "tip"), c(length(nodeTimes), length(tipTimes)))
  timesToConsider <- c(nodeTimes, tipTimes)
  .funForLambdaOptim(logLambda = log(Lambda), times = timesToConsider, nodeOrTip = nodeOrTip)
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

# getDistanceBasedClusters <- function(phyloAndTransTree, subtreeIndex = NULL, distLimit, regionLabel, criterion = c("mrca", "cophenetic", "consecutive")) {
#   numTips <- length(phyloAndTransTree$tip.label)
#   if (!is.null(subtreeIndex)) {
#     subtreeRootNodeNum <- phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum
#     subtreeRootCondition <- phyloAndTransTree$node.label[[subtreeRootNodeNum - numTips]]$region != regionLabel
#     subtreeTipsCondition <- !any((sapply(phyloAndTransTree$tip.label, "[[", "region") == regionLabel) & (sapply(phyloAndTransTree$tip.label, "[[", "subtreeIndex") == subtreeIndex)) # That condition should not be necessary...
#     if (subtreeRootCondition &  subtreeTipsCondition) return(list()) # A cluster must be contained within an introduction into the target region.
#   }
#   criterion <- criterion[[1]]
#   if (!criterion %in% c("cophenetic", "mrca", "consecutive")) {
#     stop("Clustering criterion must be either 'cophenetic', 'mrca', or 'consecutive'.")
#   }
#   transmissionTree <- .convertToTransTree(phyloAndTransTree)
#   # transmissionTree$branchMatchIndex <- sapply(1:(nrow(transmissionTree$edge) + 1), function(vertexNum) {
#   #   match(vertexNum, transmissionTree$edge[ , 2])
#   # })
#   # transmissionTree$distTipsAncestorsMatrix <- .produceDistTipsAncestorsMatrix(transmissionTree)
#   if (identical(criterion, "consecutive")) {
#     updateTransmissionTree <- function(phylogeny) { # Seems inefficient, but it shouldn't take too much time, even in large samples...
#       transmissionTree <- ape::multi2di(phylogeny)
#       nodesToFix <- which(sapply(transmissionTree$node.label, identical, "")) + length(transmissionTree$tip.label)
#       while (length(nodesToFix) > 0) {
#         currentNode <- phangorn::Ancestors(transmissionTree, nodesToFix[[1]], "parent")
#         parentSequence <- nodesToFix[[1]]
#         while (identical(.getVertexLabel(phylogeny = transmissionTree, vertexNum = currentNode), "")) {
#           parentSequence <- c(parentSequence, currentNode)
#           currentNode <- phangorn::Ancestors(transmissionTree, currentNode, "parent")
#         }
#         vertexToCopy <- .getVertexLabel(transmissionTree, currentNode)
#
#         for (nodeIndex in parentSequence) { # parentSequence cannot include a tip
#           transmissionTree$node.label[[nodeIndex - ape::Ntip(transmissionTree)]] <- vertexToCopy
#         }
#         nodesToFix <- setdiff(nodesToFix, parentSequence)
#       }
#       transmissionTree
#     }
#     transmissionTree <- updateTransmissionTree(transmissionTree)
#   }
#
#   clusterList <- list()
#
#   nodesToCheck <- length(phyloAndTransTree$tip.label) + 1
#   if (!is.null(subtreeIndex)) {
#     nodesToCheck <- phyloAndTransTree$LambdaList[[subtreeIndex]]$rootNodeNum
#   }
#
#   repeat {
#     incrementNodesToCheckFlag <- TRUE
#     nodeNumber <- nodesToCheck[[1]]
#     # Not very memory efficient, but easier to handle than a list of lists with varying depth.
#     if (nodeNumber <= ape::Ntip(transmissionTree)) {
#       if ((transmissionTree$tip.label[[nodeNumber]]$region == regionLabel) & (is.null(subtreeIndex) | identical(subtreeIndex, transmissionTree$tip.label[[nodeNumber]]$subtreeIndex))) {
#         clusterList[[length(clusterList) + 1]] <- transmissionTree$tip.label[[nodeNumber]]$name
#       }
#       incrementNodesToCheckFlag <- FALSE
#     } else {
#       if ((transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$region == regionLabel) & (is.null(subtreeIndex) | identical(subtreeIndex, transmissionTree$node.label[[nodeNumber - ape::Ntip(transmissionTree)]]$subtreeIndex))) {
#         allDescendantTips <- NULL
#         if (is.null(subtreeIndex)) {
#           allDescendantTips <- phangorn::Descendants(transmissionTree, node = nodeNumber, type = "tips")[[1]]
#         } else {
#           # allDescendantTips <- which(sapply(transmissionTree$tip.label, "[[", "subtreeIndex") == subtreeIndex)
#           allDescendantTips <- intersect(transmissionTree$tipNumsInSubtree[[subtreeIndex]], transmissionTree$descendedTips[[nodeNumber]])
#         }
#         allDescendantTipsRegions <- sapply(transmissionTree$tip.label[allDescendantTips], "[[", "region")
#         descendantTips <- allDescendantTips[allDescendantTipsRegions == regionLabel]
#         # if (!is.null(subtreeIndex)) {
#         #   descendantTipsSubtrees <- sapply(descendantTips, FUN = function(tipNum) transmissionTree$tip.label[[tipNum]]$subtreeIndex)
#         #   descendantTips <- descendantTips[descendantTipsSubtrees == subtreeIndex]
#         # }
#
#         distances <- NULL
#         if (identical(criterion, "cophenetic")) {
#           distances <- dist.tipPairs.mrca(phylogeny = transmissionTree, tipNumbers = descendantTips)$distance
#         } else if (identical(criterion, "mrca")) {
#           # distances <- dist.tips.mrca(phylogeny = transmissionTree, tipNumbers = allDescendantTips) # The MRCA should be computed based on all tips, including those not belonging to the target region.
#           # distances <- distances[allDescendantTipsRegions == regionLabel]
#           if (length(descendantTips) > 1) {
#             distances <- dist.tips.mrca(phylogeny = transmissionTree, tipNumbers = descendantTips)
#           } else if (length(descendantTips) == 1) {
#             distances <- 0 # We're dealing with a singleton. Should be noted as such.
#           } else {
#             distances <- Inf # This is an internal node whose descended tips are in the wrong subtree/region.
#           }
#         } else if (identical(criterion, "consecutive")) {
#           cladeSubTree <- ape::extract.clade(transmissionTree, node = nodeNumber)
#           tipRegions <- sapply(cladeSubTree$tip.label, "[[", "region")
#           tipsOutsideRegion <- which(tipRegions != regionLabel)
#           if ((length(tipsOutsideRegion) > 0) & (length(tipsOutsideRegion) < ape::Ntip(cladeSubTree))) {
#             cladeSubTree <- ape::drop.tip(cladeSubTree, tip = tipsOutsideRegion)
#           } else {
#             cladeSubTree <- ape::rtree(1) # This is a placeholder.
#           }
#           edgesToConsiderIndex <- cladeSubTree$edge[ , 2] > ape::Ntip(cladeSubTree) # Internal branches only
#           distances <- cladeSubTree$edge.length[edgesToConsiderIndex]
#
#           if (length(distances) == 0) distances <- Inf # This is a transmission pair: can't be a cluster under the "consecutive" criterion.
#         }
#         if (all(distances < distLimit)) {
#           clusterList[[length(clusterList) + 1]] <- sapply(descendantTips, function(x) transmissionTree$tip.label[[x]]$name)
#           incrementNodesToCheckFlag <- FALSE # Cluster has been found: stop exploring that section of the tree.
#         }
#         if (length(distances) == 0) incrementNodesToCheckFlag <- FALSE
#       }
#     }
#     if (incrementNodesToCheckFlag) {
#       # nodeChildren <- phangorn::Children(transmissionTree, node = nodeNumber)
#       nodeChildren <- transmissionTree$childrenNumList[[nodeNumber]]
#       childrenRegions <- sapply(nodeChildren, function(vertexNum) {
#         if (vertexNum <= numTips) {
#           return(transmissionTree$tip.label[[vertexNum]]$region)
#         } else {
#           return(transmissionTree$node.label[[vertexNum - numTips]]$region)
#         }
#       })
#       keepChild <- childrenRegions == regionLabel
#       if (!is.null(subtreeIndex)) {
#         childrenSubtrees <- sapply(nodeChildren, function(vertexNum) {
#           if (vertexNum <= numTips) {
#             return(transmissionTree$tip.label[[vertexNum]]$subtreeIndex)
#           } else {
#             return(transmissionTree$node.label[[vertexNum - numTips]]$subtreeIndex)
#           }
#         })
#         keepChild <- keepChild & (childrenSubtrees == subtreeIndex)
#       }
#       nodesToCheck <- c(nodesToCheck, nodeChildren[keepChild])
#     }
#
#     nodesToCheck <- nodesToCheck[-1]
#     if (length(nodesToCheck) == 0) break
#   }
#   clusterList
# }

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

# dist.tipPairs.mrca <- function(phylogeny, tipNumbers) {
#   allCombn <- combn(seq_along(tipNumbers), 2)
#   distanceByPair <- apply(allCombn, 2, function(indexPair) {
#     sum(dist.tips.mrca(phylogeny, indexPair))
#   })
#   data.frame(node1 = allCombn[1, ], node2 = allCombn[2, ], distance = distanceByPair)
# }

# dist.tips.mrca <- function(phylogeny, tipNumbers) {
#   # tipsMRCA <- ape::getMRCA(phylogeny, tipNumbers)
#   tipsMRCA <- getMRCA_Rcpp(phylogeny$parentNumVec, tipNumbers, length(phylogeny$tip.label))
#   phylogeny$distTipsAncestorsMatrix[tipNumbers, tipsMRCA - length(phylogeny$tip.label)]
# }



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

.addTransTreeEdgeLengths <- function(phyloAndTransTree) {
  numTips <- length(phyloAndTransTree$tip.label)
  edgeLengths <- sapply(seq_along(phyloAndTransTree$edge.length), function(edgeIndex) {
    childNum <- phyloAndTransTree$edge[edgeIndex, 2]
    parentNum <- phyloAndTransTree$edge[edgeIndex, 1]
    if (childNum <= numTips) {
      return(phyloAndTransTree$tip.label[[childNum]]$time - phyloAndTransTree$node.label[[parentNum - numTips]]$time)
    } else {
      return(phyloAndTransTree$node.label[[childNum - numTips]]$time - phyloAndTransTree$node.label[[parentNum - numTips]]$time)
    }
  })
  for (edgeIndex in seq_along(edgeLengths)) {
    phyloAndTransTree$edge.length[[edgeIndex]]$transmissionTree <- edgeLengths[[edgeIndex]]
  }
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
# Function does not allow a move that would allow a multifurcating split to have both tips and nodes as children.
# Function only allows transitions within subtrees.

.rNNItransTree <- function(phyloAndTransTree, moves = 1, subtreeIndex = NULL, control) {
  if (is.null(subtreeIndex)) {
    nodesInSubtrees <- table(sapply(phyloAndTransTree$node.label, "[[", "subtreeIndex"))
    sampleProbs <- (nodesInSubtrees - 1)/sum(nodesInSubtrees - 1)
    subtreeIndex <- sample(seq_along(nodesInSubtrees), prob = sampleProbs, size = 1)
  }
  nodesInSubtreeNums <- which(sapply(phyloAndTransTree$node.label, "[[", "subtreeIndex") == subtreeIndex) + length(phyloAndTransTree$tip.label)
  if (length(nodesInSubtreeNums) < 2) return(phyloAndTransTree) # No NNI move within the subtree if it's a transmission pair or a singleton.

  parentNums <- phyloAndTransTree$edge[match(nodesInSubtreeNums, phyloAndTransTree$edge[ , 2]), 1]  # The subtree root numbers will produce a FALSE.
  parentRegions <- sapply(parentNums, function(nodeNum) phyloAndTransTree$node.label[[nodeNum - length(phyloAndTransTree$tip.label)]]$region)
  currentRegions <- sapply(nodesInSubtreeNums, function(nodeNum) phyloAndTransTree$node.label[[nodeNum - length(phyloAndTransTree$tip.label)]]$region)
  parentInSubtree <- parentRegions == currentRegions
  nodesToSampleFrom <- nodesInSubtreeNums[parentInSubtree]
  childNodesToTry <- nodesToSampleFrom
  if (length(childNodesToTry) > 1)
    childNodesToTry <- sample(nodesToSampleFrom, size = length(nodesToSampleFrom))
  # It is possible that several subtrees will not allow permit a NNI move, due to tip time constraints.
  index <- 1
  while (index <= length(childNodesToTry)) {
    childNum <- childNodesToTry[[index]]
    if (control$MCMC.control$topologyTransition) {
      grandchildrenNums <- phyloAndTransTree$edge[which(childNum == phyloAndTransTree$edge[ , 1]), 2]
    } else {
      grandchildrenNums <- phyloAndTransTree$childrenList[[childNum - length(phyloAndTransTree$tip.label)]]
    }
    if (length(grandchildrenNums) > 2) { # We can't allow a transition that would create a multifurcation that would include both tips and an internal node. childNum must be for an internal node.
      index <- index + 1
      next
    }
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
    return(phyloAndTransTree)
  } else {
    siblingEdgeNum <- match(nodeSibling, phyloAndTransTree$edge[ , 2])
    grandChildrenEdgeNums <- which(phyloAndTransTree$edge[ , 1] == childNum)
    edgeNumForInterchange <- sample(grandChildrenEdgeNums, size = 1)
    treeToReturn$edge[c(siblingEdgeNum, edgeNumForInterchange), 2] <- treeToReturn$edge[c(edgeNumForInterchange, siblingEdgeNum), 2]
    # edge.length gives supporting branch lengths for nodes listed in the second column of edge. It follows that exchanging two elements in that column must result in a similar exchange of the elements of edge.length.
    treeToReturn$edge.length[c(siblingEdgeNum, edgeNumForInterchange)] <- treeToReturn$edge.length[c(edgeNumForInterchange, siblingEdgeNum)]
  }
  treeToReturn <- .updateTransTreeEdgeLengths(treeToReturn)

  inconsistencyTests <- sapply(1:(ape::Nnode(treeToReturn) - 1) + ape::Ntip(treeToReturn) + 1, function(x) {
    currentTime <- .getVertexLabel(treeToReturn, x)$time
    parentTime <- .getVertexLabel(treeToReturn, phangorn::Ancestors(treeToReturn, x, "parent"))$time
    parentTime >= currentTime
  })
  if (any(inconsistencyTests)) {
    stop("Inconsistent time values in .rNNItransTree! \n")
  }
  treeToReturn
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
