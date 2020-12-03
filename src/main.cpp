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




