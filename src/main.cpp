// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include <stdio.h>
#include <RcppArmadillo.h>

// #include <gperftools/profiler.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

arma::mat produceDistTipsAncestorsMatrixRcpp(int numTips,
                                       int numNodes,
                                       const IntegerVector & branchMatchIndexVec,
                                       const std::vector<double> & edgeLengthsVec,
                                       const IntegerVector & parentNumVec) {
  int branchMatchIndex = 0 ;
  arma::mat containerMatrix(numTips, numNodes, arma::fill::zeros) ;
  int parentNum = 0 ;

  for (uint tipNum = 0; tipNum < numTips; tipNum++) {
    int currentPos = tipNum ;
    double totalDist = 0 ;
    do {
      parentNum = parentNumVec.at(currentPos) - 1;
      branchMatchIndex = branchMatchIndexVec.at(currentPos) - 1 ;
      totalDist = totalDist + edgeLengthsVec[branchMatchIndex] ;
      containerMatrix.at(tipNum, parentNum - numTips) = totalDist ;
      currentPos = parentNum ;
    } while (currentPos > numTips) ; // numTips is the root index, rather than numTips + 1
  }
  return containerMatrix ;
}

// The getMRCA_* functions are broken.
// On some occasions, they can return a tip index as the MRCA when more than one tip index is provided for tipsNumVec.

// uint getMRCA_Rcpp(IntegerVector parentNumVec, IntegerVector tipsNumVec, uint numTips) {
//   arma::uvec parentVecArma = as<arma::uvec>(parentNumVec) ;
//   arma::uvec tipsNumVecArma = as<arma::uvec>(tipsNumVec) ;
//   arma::uvec tipsDepths(tipsNumVecArma.size()) ;
//   parentVecArma = parentVecArma - 1 ;
//   tipsNumVecArma = tipsNumVecArma - 1 ;
//   uint returnValue = numTips + 1;
//   bool allEqual = false ;
//   uint depth, currentPos ;
//
//   for (uint i = 0; i < tipsDepths.size(); i++) {
//     depth = 0 ;
//     currentPos = tipsNumVecArma(i) ;
//     do {
//       currentPos = parentVecArma(currentPos) ;
//       depth += 1 ;
//     } while (currentPos != numTips) ;
//     tipsDepths(i) = depth ;
//   }
//   uint maxDepth = max(tipsDepths) ;
//   arma::umat ancestorsMat(tipsNumVecArma.size(), maxDepth, arma::fill::zeros) ;
//
//   for (uint i = 0; i < ancestorsMat.n_rows; i++) {
//     uint offset = maxDepth - tipsDepths(i) ;
//     currentPos = tipsNumVecArma(i) ;
//     for (uint currentDepth = offset; currentDepth < maxDepth; currentDepth++) {
//       ancestorsMat(i, currentDepth) = currentPos ;
//       currentPos = parentVecArma(currentPos) ;
//     }
//   }
//
//   // We check if all values in a column are equal.
//   for (uint i = 0; i < ancestorsMat.n_cols; i++) {
//     for (uint j = 0; j < ancestorsMat.n_rows - 1; j++) {
//       allEqual = (ancestorsMat(j, i) == ancestorsMat(j + 1, i)) ;
//       if (!allEqual) break ;
//     }
//     if (allEqual) {
//       returnValue = ancestorsMat(0, i) + 1 ;
//       break ;
//     }
//   }
//
//   return returnValue ;
// }

// Function is bugged
// On some occasions, they can return a tip index as the MRCA when more than one tip index is provided for tipsNumVec.
// Fix before uncommenting and using.
// uint getMRCA(std::vector<uint> & parentNumVec, std::vector<uint> & tipsNumVec, int & numTips) {
//   arma::uvec tipsDepths(tipsNumVec.size()) ;
//   uint returnValue = numTips ;
//   bool allEqual = false ;
//   uint depth, currentPos ;
//
//   for (uint i = 0; i < tipsDepths.size(); i++) {
//     depth = 0 ;
//     currentPos = tipsNumVec[i] ;
//     do {
//       currentPos = parentNumVec[currentPos] ;
//       depth += 1 ;
//     } while (currentPos != numTips) ;
//     tipsDepths(i) = depth ;
//   }
//   uint maxDepth = max(tipsDepths) ;
//   arma::umat ancestorsMat(tipsNumVec.size(), maxDepth, arma::fill::zeros) ;
//
//   for (uint i = 0; i < ancestorsMat.n_rows; i++) {
//     uint offset = maxDepth - tipsDepths(i) ;
//     currentPos = tipsNumVec[i] ;
//     for (uint currentDepth = offset; currentDepth < maxDepth; currentDepth++) {
//       ancestorsMat.at(i, currentDepth) = currentPos ;
//       currentPos = parentNumVec[currentPos] ;
//     }
//   }
//
//   // We check if all values in a column are equal.
//   for (uint i = 0; i < ancestorsMat.n_cols; i++) {
//     for (uint j = 0; j < ancestorsMat.n_rows - 1; j++) {
//       allEqual = (ancestorsMat.at(j, i) == ancestorsMat.at(j + 1, i)) ;
//       if (!allEqual) break ;
//     }
//     if (allEqual) {
//       returnValue = ancestorsMat.at(0, i) ;
//       break ;
//     }
//   }
//
//   return returnValue ;
// }

// [[Rcpp::export]]

std::unordered_map<std::string, int> getMRCAclustersRcpp(
    const IntegerVector & parentNumVec,
    const List & childrenNumList,
    const List & descendedTipsList,
    const IntegerVector & subtreeIndexVec,
    const StringVector & vertexRegionVec,
    const StringVector & tipNamesVec,
    int subtreeRootNum,
    const NumericMatrix & distTipsAncestorsMatrix,
    int subtreeIndex,
    int numTips,
    std::string regionLabel,
    int distLimit,
    bool printDetails = false) {
  // ProfilerStart("/home/luc/temp/profile.log") ;
  // std::vector<std::vector<std::string>> clusterList ;
  // for (uint i = 0; i < 45000; i++) {
  bool incrementNodesToCheckFlag = TRUE ;
  int nodeNumber = 0 ;
  std::vector<int> nodesToCheck ;
  nodesToCheck.push_back(subtreeRootNum - 1) ;
  std::unordered_map<std::string, int> clusMemIndices ;
  for (uint i = 0; i < tipNamesVec.size(); i++) {
    if ((subtreeIndexVec.at(i) == subtreeIndex) & (vertexRegionVec.at(i) == regionLabel)) {
      std::string mapKey = Rcpp::as<std::string>(tipNamesVec.at(i)) ;
      clusMemIndices[mapKey] = int(0) ;
    }
  }
  int clusterNumber = 1 ;
  do {
    incrementNodesToCheckFlag = TRUE ;
    nodeNumber = nodesToCheck.back() ;
    nodesToCheck.pop_back() ; // We pre-emptively remove the element we are checking ;
    if (nodeNumber < numTips) {
      if ((vertexRegionVec(nodeNumber) == regionLabel) & (subtreeIndexVec(nodeNumber) == subtreeIndex)) {
        clusMemIndices[Rcpp::as<std::string>(tipNamesVec(nodeNumber))] = clusterNumber ;
        clusterNumber += 1 ;
      }
      incrementNodesToCheckFlag = false ;
    } else {
      bool testValue ;
      if (subtreeIndexVec(nodeNumber) == subtreeIndex) {
        if (Rcpp::as<NumericVector>(descendedTipsList[nodeNumber]).size() > 1) {
          bool checkValue = TRUE;
          for (uint i = 0; i < Rcpp::as<NumericVector>(descendedTipsList[nodeNumber]).size(); i++) {
            checkValue = distTipsAncestorsMatrix(Rcpp::as<IntegerVector>(descendedTipsList[nodeNumber])[i] - 1, nodeNumber - numTips) < distLimit ;
            if (!checkValue) break ;
          }
          if (checkValue) {
            for (uint i = 0; i < Rcpp::as<NumericVector>(descendedTipsList[nodeNumber]).size(); i++) {
              std::string mapKey = Rcpp::as<std::string>(tipNamesVec(Rcpp::as<IntegerVector>(descendedTipsList[nodeNumber])[i] - 1)) ;
              clusMemIndices[mapKey] = clusterNumber ;
            }
            clusterNumber += 1 ;
            incrementNodesToCheckFlag = false ; // Cluster has been found: stop exploring that section of the tree.
          }
        } else if (Rcpp::as<NumericVector>(descendedTipsList[nodeNumber]).size() == 1) {
          // It is automatically a singleton that should be added...
          std::string mapKey = Rcpp::as<std::string>(tipNamesVec.at(Rcpp::as<IntegerVector>(descendedTipsList[nodeNumber])[0] - 1)) ;
          clusMemIndices[mapKey] = clusterNumber;
          clusterNumber += 1 ;
          incrementNodesToCheckFlag = false ;
        } // descendedTipsListStd[nodeNumber] can have size 0
      }
    }
    if (incrementNodesToCheckFlag) {
      bool keepChildTest ;

      for (uint i = 0; i < Rcpp::as<IntegerVector>(childrenNumList[nodeNumber]).size(); i++) {
        int childSubtree = subtreeIndexVec(Rcpp::as<IntegerVector>(childrenNumList[nodeNumber])[i] - 1) ;
        keepChildTest = (childSubtree == subtreeIndex) ;
        if (keepChildTest) {
          nodesToCheck.push_back(Rcpp::as<IntegerVector>(childrenNumList[nodeNumber])[i] - 1) ;
        }
      }
    }
  } while (nodesToCheck.size() > 0) ;
  // }
  // ProfilerStop() ;
  // Rcpp::stop("Stop for profiling.") ;
  return clusMemIndices ;
}

// TO DO: Remove call to getMRCA. Not needed anymore.
// Rcpp::List getCopheneticClustersRcpp(
//     IntegerVector & parentNumVec,
//     List & childrenNumList,
//     List & descendedTipsList,
//     IntegerVector & subtreeIndexVec,
//     CharacterVector & vertexRegionVec,
//     CharacterVector & tipNamesVec,
//     uint & subtreeRootNum,
//     NumericMatrix & distTipsAncestorsMatrix,
//     int subtreeIndex,
//     int numTips,
//     std::string regionLabel,
//     int distLimit) {
//   std::vector<std::vector<std::string>> clusterList ;
//   bool incrementNodesToCheckFlag = TRUE ;
//   uint nodeNumber = 0 ;
//   std::vector<uint> parentNumVecStd = Rcpp::as<std::vector<uint>>(parentNumVec) ;
//   std::transform(parentNumVecStd.begin(), parentNumVecStd.end(), parentNumVecStd.begin(),
//                  [] (uint & i) { return i - 1 ;}) ;
//   std::vector<std::vector<uint>> childrenNumListTypecast ;
//
//   for (uint i = 0; i < childrenNumList.size(); i++) {
//     if (Rf_isNull(childrenNumList[i])) {
//       std::vector<uint> placeholder ; // A 0-length vector, since tips don't have children.
//       childrenNumListTypecast.push_back(placeholder) ;
//     } else {
//       std::vector<uint> vectorToInclude = Rcpp::as<std::vector<uint>>(Rcpp::as<IntegerVector>(childrenNumList[i])) ;
//       std::transform(vectorToInclude.begin(), vectorToInclude.end(), vectorToInclude.begin(),
//                      [] (uint & i) { return i - 1 ;}) ;
//       childrenNumListTypecast.push_back(vectorToInclude) ;
//     }
//   }
//   std::vector<std::vector<uint>> descendedTipsListStd ;
//   for (uint i = 0 ; i < descendedTipsList.size() ; i++) {
//     descendedTipsListStd.push_back(Rcpp::as<std::vector<uint>>(descendedTipsList(i))) ;
//     std::transform(descendedTipsListStd.back().begin(), descendedTipsListStd.back().end(), descendedTipsListStd.back().begin(),
//                    [] (uint & i) { return i - 1 ;}) ;
//   }
//
//   std::vector<uint> tipsSubtreeIndexVec  ;
//   std::vector<std::string> tipsRegionVec ;
//   std::string currentRegion ;
//   for (uint i = 0; i < numTips; i++) {
//     tipsSubtreeIndexVec.push_back(subtreeIndexVec(i)) ;
//     currentRegion = vertexRegionVec(i) ;
//     tipsRegionVec.push_back(currentRegion) ;
//   }
//   std::string tipName ;
//   std::vector<uint> nodesToCheck ;
//   nodesToCheck.push_back(subtreeRootNum - 1) ;
//   do {
//     incrementNodesToCheckFlag = TRUE ;
//     nodeNumber = nodesToCheck.back() ;
//     nodesToCheck.pop_back() ; // We pre-emptively remove the element we are checking ;
//     if (nodeNumber < numTips) {
//       currentRegion = vertexRegionVec(nodeNumber) ;
//       if ((currentRegion == regionLabel) & (subtreeIndex == subtreeIndexVec(nodeNumber))) {
//         tipName = tipNamesVec(nodeNumber) ;
//         std::vector<std::string> outputVec ;
//         outputVec.push_back(tipName) ;
//         clusterList.push_back(outputVec) ;
//       }
//       incrementNodesToCheckFlag = false ;
//     } else {
//       bool testValue ;
//
//       if ((vertexRegionVec(nodeNumber) == regionLabel) & (subtreeIndexVec(nodeNumber) == subtreeIndex)) {
//         std::vector<uint> descendantTips ;
//         for (uint i = 0; i < descendedTipsListStd.at(nodeNumber).size(); i++) {
//           testValue = (tipsSubtreeIndexVec.at(descendedTipsListStd.at(nodeNumber).at(i)) == subtreeIndex) & (tipsRegionVec.at(descendedTipsListStd.at(nodeNumber).at(i)) == regionLabel) ;
//           if (testValue) descendantTips.push_back(descendedTipsListStd.at(nodeNumber).at(i)) ;
//         }
//
//         if (descendantTips.size() > 1) {
//           std::vector<std::vector<uint>> tipPairsVec ;
//           for (uint i = 0; i < descendantTips.size() - 1; i++) {
//             for (uint j = i + 1; j < descendantTips.size(); j++) {
//               std::vector<uint> tipPair ;
//               tipPair.push_back(descendantTips.at(i)) ;
//               tipPair.push_back(descendantTips.at(j)) ;
//               tipPairsVec.push_back(tipPair) ;
//             }
//           }
//           arma::vec distances(tipPairsVec.size()) ;
//           // distances = dist.tips.mrca(phylogeny = transmissionTree, tipNumbers = descendantTips) ;
//
//           for (uint i = 0 ; i < distances.size(); i++) {
//             uint mrca = getMRCA(parentNumVecStd, tipPairsVec.at(i), numTips) ;
//             distances(i) =  distTipsAncestorsMatrix(tipPairsVec.at(i).at(0), mrca - numTips) +
//               distTipsAncestorsMatrix(tipPairsVec.at(i).at(1), mrca - numTips) ;
//           }
//           bool checkValue = TRUE;
//           for (uint i = 0; i < distances.size(); i++) {
//             // for (auto & distValue : distances) {
//             checkValue = distances(i) < distLimit ;
//             if (!checkValue) break ;
//           }
//           if (checkValue) {
//             std::vector<std::string> tipsToKeep(distances.size()) ;
//             for (uint i = 0; i < tipsToKeep.size(); i++) {
//               tipsToKeep.at(i) = tipNamesVec(descendantTips.at(i)) ;
//             }
//             clusterList.push_back(tipsToKeep) ;
//             incrementNodesToCheckFlag = false ; // Cluster has been found: stop exploring that section of the tree.
//           }
//         } else if (descendantTips.size() == 1) {
//           // It is automatically a singleton that should be added...
//           std::vector<std::string> tipsToKeep(1) ;
//           tipsToKeep.at(0) = tipNamesVec.at(descendantTips.at(0)) ;
//           clusterList.push_back(tipsToKeep) ;
//           incrementNodesToCheckFlag = false ;
//         }
//       }
//     }
//     if (incrementNodesToCheckFlag) {
//       std::vector<uint> nodeChildren = childrenNumListTypecast.at(nodeNumber) ;
//       bool keepChildTest ;
//
//       for (uint i = 0; i < nodeChildren.size(); i++) {
//         std::string childRegion = Rcpp::as<std::string>(vertexRegionVec(nodeChildren.at(i))) ;
//         uint childSubtree = subtreeIndexVec(nodeChildren.at(i)) ;
//         keepChildTest = (childRegion == regionLabel) & (childSubtree == subtreeIndex) ;
//         if (keepChildTest) {
//           nodesToCheck.push_back(nodeChildren.at(i)) ;
//         }
//       }
//     }
//   } while (nodesToCheck.size() > 0) ;
//   return Rcpp::wrap(clusterList) ;
// }

// [[Rcpp::export]]

List simulateNodeTimesRcpp(
    int numTips,
    const NumericVector & baseRatePerIntroduction,
    const IntegerVector & orderedVertices,
    const IntegerVector & subtreeIndexVec,
    const NumericVector & tipTimes,
    const IntegerMatrix & edgeMatrix,
    const List & childrenNumList,
    const IntegerVector & branchMatchIndexVec,
    const IntegerVector & parentNumVec) {
  std::vector<double> vertexTimes(subtreeIndexVec.size()) ;
  for (uint i = 0; i < tipTimes.size(); i++) {
    vertexTimes[i] = tipTimes(i) ;
  }

  for (auto & nodeNum : orderedVertices) {
    if (nodeNum <= numTips) continue ;
    int subtreeForMerge = subtreeIndexVec(nodeNum - 1) ; // -1 will have to be applied for indexing
    int elementSize = Rcpp::as<IntegerVector>(childrenNumList[nodeNum - 1]).size() ;
    arma::vec childrenTimes(elementSize) ;
    for (uint j = 0; j < childrenTimes.size(); j++) {
      childrenTimes.at(j) = vertexTimes[Rcpp::as<IntegerVector>(childrenNumList[nodeNum - 1])[j] - 1];
    }

    double minChildrenTimes = min(childrenTimes) ;
    double coalRate = baseRatePerIntroduction(subtreeForMerge - 1) ;
    double expValue = R::rexp(1 / coalRate) ;
    vertexTimes[nodeNum - 1] = minChildrenTimes - expValue ;
  }

  std::vector<double> edgeLengths(subtreeIndexVec.size() - 1) ;
  for (uint i = 0 ; i < edgeMatrix.rows(); i++) {
    edgeLengths[i] = vertexTimes[edgeMatrix(i, 1) - 1] - vertexTimes[edgeMatrix(i, 0) - 1] ;
  }

  std::vector<double> nodeTimes ;

  for (uint i = numTips; i < vertexTimes.size(); i++) {
    nodeTimes.push_back(vertexTimes[i]) ;
  }

  arma::mat distTipsAncestors = produceDistTipsAncestorsMatrixRcpp(numTips,
                                     subtreeIndexVec.size() - numTips,
                                     branchMatchIndexVec,
                                     edgeLengths,
                                     parentNumVec) ;

  return List::create(Named("vertexTimes") = Rcpp::wrap(nodeTimes), Named("edgeLengths") = Rcpp::wrap(edgeLengths), Named("distTipsAncestorsMatrix") = distTipsAncestors) ;
}

arma::sp_imat getCoclusterMat(const IntegerVector & clusMemVec) {
  arma::ivec uniqueValues = arma::unique(Rcpp::as<arma::ivec>(clusMemVec)) ;
  arma::umat locations = arma::umat(2, clusMemVec.size()) ;
  for (uint colIndex = 0; colIndex < locations.n_cols; colIndex++) {
    locations.at(0, colIndex) = locations.at(1, colIndex) = colIndex ;
  }
  arma::uvec indices ;
  for (auto & clusNum : uniqueValues) {
    indices = arma::find(Rcpp::as<arma::ivec>(clusMemVec) == clusNum) ;
    int numCols = int(indices.size() * (indices.size() - 1) / 2) ;
    arma::umat matToMerge(2, numCols) ;
    int matToMergeCol = 0 ;
    for (uint i = 0; i < indices.size() - 1; i++) {
      for (uint j = i + 1; j < indices.size(); j++) {
        matToMerge.at(0, matToMergeCol) = indices.at(j) ;
        matToMerge.at(1, matToMergeCol) = indices.at(i) ; // We want a lower-triangular matrix, hence j then i.
        matToMergeCol++ ;
      }
    }
    locations = arma::join_rows(locations, matToMerge) ;
  }
  return arma::sp_imat(locations, arma::ivec(locations.n_cols, arma::fill::ones)) ;
}

// [[Rcpp::export]]

arma::sp_imat getSumMatRcpp(const List & clusMemVecList) {
  IntegerVector vecToTransform = Rcpp::as<IntegerVector>(clusMemVecList.at(0)) ;
  arma::sp_imat resultMatrix = getCoclusterMat(vecToTransform) ;
  for (uint i = 1; i < clusMemVecList.size(); i++) {
    vecToTransform = Rcpp::as<IntegerVector>(clusMemVecList.at(i)) ;
    resultMatrix += getCoclusterMat(vecToTransform) ;
  }
  return resultMatrix ;
}

// [[Rcpp::export]]

std::unordered_map<int, double> summariseClusSizeDistsRcpp(const List & clusMemVecList) {
  arma::ivec sampleVec = Rcpp::as<arma::ivec>(clusMemVecList.at(0)) ;
  std::unordered_map<int, double> freqVec ;

  for (uint i = 0; i < clusMemVecList.size(); i++) {
    sampleVec = Rcpp::as<arma::ivec>(clusMemVecList.at(i)) ;
    arma::uvec clusSizes = arma::hist(sampleVec, unique(sampleVec)) ;
    for (auto & clusSize : clusSizes) {
      double incrementValue = 1/double(clusSizes.size() * clusMemVecList.size()) ;
      if (freqVec.find(clusSize) == freqVec.end()) {
        freqVec[clusSize] = incrementValue ;
      } else {
        freqVec[clusSize] += incrementValue ;
      }
    }
  }
  return freqVec ;
}
