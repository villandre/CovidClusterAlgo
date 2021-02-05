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

// [[Rcpp::export]]
arma::mat produceDistTipsAncestorsMatrixRcpp(uint numTips,
                                       uint numNodes,
                                       IntegerVector branchMatchIndexVec,
                                       NumericVector edgeLengthsVec,
                                       IntegerVector parentNumVec) {
  uint branchMatchIndex = 0 ;
  arma::mat containerMatrix(numTips, numNodes, arma::fill::zeros) ;
  uint parentNum = 0 ;

  for (uint tipNum = 0; tipNum < numTips; tipNum++) {
    uint currentPos = tipNum ;
    double totalDist = 0 ;
    do {
      parentNum = parentNumVec.at(currentPos) ;
      branchMatchIndex = branchMatchIndexVec.at(currentPos) ;
      totalDist = totalDist + edgeLengthsVec.at(branchMatchIndex) ;
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

Rcpp::List getMRCAclustersRcpp(
    IntegerVector & parentNumVec,
    List & childrenNumList,
    List & descendedTipsList,
    IntegerVector & subtreeIndexVec,
    CharacterVector & vertexRegionVec,
    CharacterVector & tipNamesVec,
    uint & subtreeRootNum,
    NumericMatrix & distTipsAncestorsMatrix,
    int subtreeIndex,
    int numTips,
    std::string regionLabel,
    int distLimit) {
  // ProfilerStart("/home/luc/temp/profile.log") ;
  // std::vector<std::vector<std::string>> clusterList ;
  // for (uint i = 0; i < 20000; i++) {
  std::vector<std::vector<std::string>> clusterList ;
  bool incrementNodesToCheckFlag = TRUE ;
  uint nodeNumber = 0 ;
  std::vector<uint> parentNumVecStd = Rcpp::as<std::vector<uint>>(parentNumVec) ;
  std::transform(parentNumVecStd.begin(), parentNumVecStd.end(), parentNumVecStd.begin(),
                 [] (uint & i) { return i - 1 ;}) ;
  std::vector<std::vector<uint>> childrenNumListTypecast ;

  for (uint i = 0; i < childrenNumList.size(); i++) {
    if (Rf_isNull(childrenNumList[i])) {
      std::vector<uint> placeholder ; // A 0-length vector, since tips don't have children.
      childrenNumListTypecast.push_back(placeholder) ;
    } else {
      std::vector<uint> vectorToInclude = Rcpp::as<std::vector<uint>>(Rcpp::as<IntegerVector>(childrenNumList[i])) ;
      std::transform(vectorToInclude.begin(), vectorToInclude.end(), vectorToInclude.begin(),
                     [] (uint & i) { return i - 1 ;}) ;
      childrenNumListTypecast.push_back(vectorToInclude) ;
    }
  }
  std::vector<std::vector<uint>> descendedTipsListStd ;
  for (uint i = 0 ; i < descendedTipsList.size() ; i++) {
    descendedTipsListStd.push_back(Rcpp::as<std::vector<uint>>(descendedTipsList(i))) ;
    std::transform(descendedTipsListStd.back().begin(), descendedTipsListStd.back().end(), descendedTipsListStd.back().begin(),
                   [] (uint & i) { return i - 1 ;}) ;
  }

  std::vector<uint> tipsSubtreeIndexVec  ;
  std::vector<std::string> tipsRegionVec ;
  std::string currentRegion ;
  for (uint i = 0; i < numTips; i++) {
    tipsSubtreeIndexVec.push_back(subtreeIndexVec(i)) ;
    currentRegion = vertexRegionVec(i) ;
    tipsRegionVec.push_back(currentRegion) ;
  }
  std::string tipName ;
  std::vector<uint> nodesToCheck ;
  nodesToCheck.push_back(subtreeRootNum - 1) ;
  do {
    incrementNodesToCheckFlag = TRUE ;
    nodeNumber = nodesToCheck.back() ;
    nodesToCheck.pop_back() ; // We pre-emptively remove the element we are checking ;
    if (nodeNumber < numTips) {
      currentRegion = vertexRegionVec(nodeNumber) ;
      if ((currentRegion == regionLabel) & (subtreeIndex == subtreeIndexVec(nodeNumber))) {
        tipName = tipNamesVec(nodeNumber) ;
        std::vector<std::string> outputVec ;
        outputVec.push_back(tipName) ;
        clusterList.push_back(outputVec) ;
      }
      incrementNodesToCheckFlag = false ;
    } else {
      bool testValue ;

      if ((vertexRegionVec(nodeNumber) == regionLabel) & (subtreeIndexVec(nodeNumber) == subtreeIndex)) {

        arma::vec distances(descendedTipsListStd[nodeNumber].size())  ;

        if (descendedTipsListStd[nodeNumber].size() > 1) {
          // distances = dist.tips.mrca(phylogeny = transmissionTree, tipNumbers = descendantTips) ;
          for (uint i = 0 ; i < distances.size(); i++) {
            distances(i) =  distTipsAncestorsMatrix(descendedTipsListStd[nodeNumber][i], nodeNumber - numTips) ;
          }
          bool checkValue = TRUE;
          for (uint i = 0; i < distances.size(); i++) {
            // for (auto & distValue : distances) {
            checkValue = distances(i) < distLimit ;
            if (!checkValue) break ;
          }
          if (checkValue) {
            std::vector<std::string> tipsToKeep(distances.size()) ;
            for (uint i = 0; i < tipsToKeep.size(); i++) {
              tipsToKeep[i] = tipNamesVec(descendedTipsListStd[nodeNumber][i]) ;
            }
            clusterList.push_back(tipsToKeep) ;
            // distances.print("Distances in cluster:") ;
            // for (auto & i : tipsToKeep) Rcout << "Tip in cluster" << i << std::endl ;
            incrementNodesToCheckFlag = false ; // Cluster has been found: stop exploring that section of the tree.
          }
        } else if (descendedTipsListStd[nodeNumber].size() == 1) {
          // It is automatically a singleton that should be added...
          std::vector<std::string> tipsToKeep(1) ;
          tipsToKeep[0] = tipNamesVec.at(descendedTipsListStd[nodeNumber][0]) ;
          clusterList.push_back(tipsToKeep) ;
          incrementNodesToCheckFlag = false ;
        }
      }
    }
    if (incrementNodesToCheckFlag) {
      std::vector<uint> nodeChildren = childrenNumListTypecast[nodeNumber] ;
      bool keepChildTest ;

      for (uint i = 0; i < nodeChildren.size(); i++) {
        std::string childRegion = Rcpp::as<std::string>(vertexRegionVec(nodeChildren[i])) ;
        uint childSubtree = subtreeIndexVec(nodeChildren[i]) ;
        keepChildTest = (childRegion == regionLabel) & (childSubtree == subtreeIndex) ;
        if (keepChildTest) {
          nodesToCheck.push_back(nodeChildren[i]) ;
        }
      }
    }
  } while (nodesToCheck.size() > 0) ;
  // }
  // ProfilerStop() ;
  // Rcpp::stop("Stop for profiling.") ;
  return Rcpp::wrap(clusterList) ;
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
    uint numTips,
    NumericVector & baseRatePerIntroduction,
    IntegerVector & orderedVertices,
    IntegerVector & subtreeIndexVec,
    NumericVector & tipTimes,
    IntegerMatrix & edgeMatrix,
    List & childrenNumList) {
  std::vector<std::vector<uint>> childrenNumListTypecast ;
  std::vector<double> vertexTimes(subtreeIndexVec.size()) ;
  for (uint i = 0; i < tipTimes.size(); i++) {
    vertexTimes[i] = tipTimes(i) ;
  }

  for (uint i = 0; i < childrenNumList.size(); i++) {
    if (Rf_isNull(childrenNumList[i])) {
      std::vector<uint> placeholder ; // A 0-length vector, since tips don't have children.
      childrenNumListTypecast.push_back(placeholder) ;
    } else {
      std::vector<uint> vectorToInclude = Rcpp::as<std::vector<uint>>(Rcpp::as<IntegerVector>(childrenNumList[i])) ;
      std::transform(vectorToInclude.begin(), vectorToInclude.end(), vectorToInclude.begin(),
                     [] (uint & i) { return i - 1 ;}) ;
      childrenNumListTypecast.push_back(vectorToInclude) ;
    }
  }
  std::vector<uint> orderedNodes ;
  for (auto & vertexNum : orderedVertices) {
    if (vertexNum > numTips) { // R indexing, so >, not >=.
      orderedNodes.push_back(vertexNum - 1) ; // Bringing it back to C++ indexing.
    }
  }

  for (auto & nodeNum : orderedNodes) {
    std::vector<uint> childrenNums = childrenNumListTypecast[nodeNum] ;
    uint subtreeForMerge = subtreeIndexVec(nodeNum) ; // -1 will have to be applied for indexing
    arma::vec childrenTimes(childrenNums.size()) ;
    for (uint j = 0; j < childrenTimes.size(); j++) {
      childrenTimes.at(j) = vertexTimes[childrenNums.at(j)];
    }

    arma::vec minChildrenTimes(1) ; // The odd setup here is due to arma::randg returning a vector.
    minChildrenTimes(0) = min(childrenTimes) ;
    arma::vec vecWithTimeValue = minChildrenTimes - arma::randg(1, arma::distr_param(double(1), baseRatePerIntroduction(subtreeForMerge - 1))) ;
    vertexTimes[nodeNum] = vecWithTimeValue(0) ;
  }

  std::vector<double> edgeLengths(subtreeIndexVec.size() - 1) ;
  for (uint i = 0 ; i < edgeMatrix.rows(); i++) {
    edgeLengths[i] = vertexTimes[edgeMatrix(i, 1) - 1] - vertexTimes[edgeMatrix(i, 0) - 1] ;
  }

  std::vector<double> nodeTimes ;

  for (uint i = numTips; i < vertexTimes.size(); i++) {
    nodeTimes.push_back(vertexTimes[i]) ;
  }

  return List::create(Named("vertexTimes") = Rcpp::wrap(nodeTimes), Named("edgeLengths") = Rcpp::wrap(edgeLengths)) ;
}


