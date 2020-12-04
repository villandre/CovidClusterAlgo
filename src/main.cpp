// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include <stdio.h>
#include <RcppArmadillo.h>

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

// [[Rcpp::export]]

uint getMRCA_Rcpp(IntegerVector parentNumVec, IntegerVector tipsNumVec, uint numTips) {
  arma::uvec parentVecArma = as<arma::uvec>(parentNumVec) ;
  arma::uvec tipsNumVecArma = as<arma::uvec>(tipsNumVec) ;
  arma::uvec tipsDepths(tipsNumVecArma.size()) ;
  parentVecArma = parentVecArma - 1 ;
  tipsNumVecArma = tipsNumVecArma - 1 ;
  uint returnValue = numTips + 1;
  bool allEqual = false ;
  uint depth, currentPos ;

  for (uint i = 0; i < tipsDepths.size(); i++) {
    depth = 0 ;
    currentPos = tipsNumVecArma(i) ;
    do {
      currentPos = parentVecArma(currentPos) ;
      depth += 1 ;
    } while (currentPos != numTips) ;
    tipsDepths(i) = depth ;
  }
  uint maxDepth = max(tipsDepths) ;
  arma::umat ancestorsMat(tipsNumVecArma.size(), maxDepth, arma::fill::zeros) ;

  for (uint i = 0; i < ancestorsMat.n_rows; i++) {
    uint offset = maxDepth - tipsDepths(i) ;
    currentPos = tipsNumVecArma(i) ;
    for (uint currentDepth = offset; currentDepth < maxDepth; currentDepth++) {
      ancestorsMat(i, currentDepth) = currentPos ;
      currentPos = parentVecArma(currentPos) ;
    }
  }

  // We check if all values in a column are equal.
  for (uint i = 0; i < ancestorsMat.n_cols; i++) {
    for (uint j = 0; j < ancestorsMat.n_rows - 1; j++) {
      allEqual = (ancestorsMat(j, i) == ancestorsMat(j + 1, i)) ;
      if (!allEqual) break ;
    }
    if (allEqual) {
      returnValue = ancestorsMat(0, i) + 1 ;
      break ;
    }
  }

  return returnValue ;
}

uint getMRCA(std::vector<uint> & parentNumVec, std::vector<uint> & tipsNumVec, int & numTips) {
  arma::uvec tipsDepths(tipsNumVec.size()) ;
  for (auto & vecElement : parentNumVec) vecElement = vecElement - 1 ;
  for (auto & vecElement : tipsNumVec) vecElement = vecElement - 1 ;
  uint returnValue = numTips ;
  bool allEqual = false ;
  uint depth, currentPos ;

  for (uint i = 0; i < tipsDepths.size(); i++) {
    depth = 0 ;
    currentPos = tipsNumVec.at(i) ;
    do {
      currentPos = parentNumVec.at(currentPos) ;
      depth += 1 ;
    } while (currentPos != numTips) ;
    tipsDepths(i) = depth ;
  }
  uint maxDepth = max(tipsDepths) ;
  arma::umat ancestorsMat(tipsNumVec.size(), maxDepth, arma::fill::zeros) ;

  for (uint i = 0; i < ancestorsMat.n_rows; i++) {
    uint offset = maxDepth - tipsDepths(i) ;
    currentPos = tipsNumVec.at(i) ;
    for (uint currentDepth = offset; currentDepth < maxDepth; currentDepth++) {
      ancestorsMat(i, currentDepth) = currentPos ;
      currentPos = parentNumVec.at(currentPos) ;
    }
  }

  // We check if all values in a column are equal.
  for (uint i = 0; i < ancestorsMat.n_cols; i++) {
    for (uint j = 0; j < ancestorsMat.n_rows - 1; j++) {
      allEqual = (ancestorsMat(j, i) == ancestorsMat(j + 1, i)) ;
      if (!allEqual) break ;
    }
    if (allEqual) {
      returnValue = ancestorsMat(0, i) ;
      break ;
    }
  }

  return returnValue ;
}

// [[Rcpp::export]]

Rcpp::List getMRCAclustersRcpp(
    IntegerVector & parentNumVec,
    List & childrenNumList,
    IntegerVector & branchMatchIndexVec,
    NumericVector & branchLengthsVec,
    IntegerVector & subtreeIndexVec,
    CharacterVector & vertexRegionVec,
    CharacterVector & tipNamesVec,
    IntegerVector & subtreeRootNum,
    NumericMatrix & distTipsAncestorsMatrix,
    int subtreeIndex,
    int numTips,
    std::string regionLabel,
    int distLimit) {
  std::vector<std::vector<std::string>> clusterList ;
  std::vector<uint> nodesToCheck ;
  nodesToCheck.push_back(subtreeRootNum(subtreeIndex) - 1) ;
  bool incrementNodesToCheckFlag = TRUE ;
  uint nodeNumber = 0 ;
  std::string currentRegion ;
  std::string tipName ;
  std::vector<uint> tipsSubtreeIndexVec  ;
  std::vector<std::string> tipsRegionVec ;
  std::vector<uint> parentNumVecStd = Rcpp::as<std::vector<uint>>(parentNumVec) ;

  std::vector<arma::uvec> childrenNumListTypecast = Rcpp::as<std::vector<arma::uvec>>(childrenNumList) ;

  for (uint i = 0; i < numTips; i++) {
    tipsSubtreeIndexVec.push_back(subtreeIndexVec(i)) ;
    currentRegion = vertexRegionVec(i) ;
    tipsRegionVec.push_back(currentRegion) ;
  }

  do {
    incrementNodesToCheckFlag = TRUE ;
    nodeNumber = nodesToCheck.back() ;
    nodesToCheck.pop_back() ; // We pre-emptively remove the element we are checking ;
    if (nodeNumber <= numTips) {
      currentRegion = vertexRegionVec(nodeNumber) ;
      if ((currentRegion == regionLabel) & (subtreeIndex == subtreeIndexVec(nodeNumber))) {
        tipName = tipNamesVec(nodeNumber) ;
        std::vector<std::string> outputVec ;
        outputVec.push_back(tipName) ;
        clusterList.push_back(outputVec) ;
      }
      incrementNodesToCheckFlag = FALSE ;
    } else {
      bool testValue ;
      std::vector<uint> descendantTips ;
      if ((vertexRegionVec(nodeNumber) == regionLabel) & (subtreeIndexVec(nodeNumber) == subtreeIndex)) {
        for (uint i = 0; i < tipsSubtreeIndexVec.size(); i++) {
          testValue = (tipsSubtreeIndexVec.at(i) == subtreeIndex) & (tipsRegionVec.at(i) == regionLabel) ;
          if (testValue) descendantTips.push_back(i) ;
        }

        arma::vec distances(descendantTips.size())  ;

        if (descendantTips.size() > 1) {
          // distances = dist.tips.mrca(phylogeny = transmissionTree, tipNumbers = descendantTips) ;
          uint mrca = getMRCA(parentNumVecStd, descendantTips, numTips) ;
          for (uint i = 0 ; i < distances.size(); i++) {
            distances(i) =  distTipsAncestorsMatrix(descendantTips.at(i), mrca) ;
          }
        } else if (descendantTips.size() == 1) {
          distances(0) = 0 ; // We're dealing with a singleton. Should be noted as such.
        } else {
          distances(0) = arma::datum::inf ; // This is an internal node whose descended tips are in the wrong subtree/region.
        }
        bool checkValue = TRUE;
        for (auto & distValue : distances) {
          checkValue = distValue < distLimit ;
          if (!checkValue) break ;
        }
        if (checkValue) {
          std::vector<std::string> tipsToKeep(distances.size()) ;
          for (uint i = 0; i < tipsToKeep.size(); i++) {
            tipsToKeep.at(i) = tipNamesVec(i) ;
          }
          clusterList.push_back(tipsToKeep) ;
          incrementNodesToCheckFlag = false ; // Cluster has been found: stop exploring that section of the tree.
        }
        if (distances.size() == 0) incrementNodesToCheckFlag = false ;
      }
    }
    if (incrementNodesToCheckFlag) {
      arma::uvec nodeChildren = childrenNumListTypecast.at(nodeNumber) ;
      std::vector<std::string> childrenRegions ;
      for (auto & nodeChild : nodeChildren) {
        childrenRegions.push_back(Rcpp::as<std::string>(vertexRegionVec(nodeChild))) ;
      }
      arma::uvec childrenSubtrees = Rcpp::as<arma::uvec>(subtreeIndexVec).elem(nodeChildren) ;

      bool keepChildTest ;

      for (uint i = 0; i < childrenRegions.size(); i++) {
        keepChildTest = (childrenRegions.at(i) == regionLabel) & (childrenSubtrees(i) == subtreeIndex) ;
        if (keepChildTest) {
          nodesToCheck.push_back(nodeChildren(i)) ;
        }
      }
    }
  } while (nodesToCheck.size() > 0) ;
  return Rcpp::wrap(clusterList) ;
}
