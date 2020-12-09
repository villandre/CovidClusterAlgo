// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// produceDistTipsAncestorsMatrixRcpp
arma::mat produceDistTipsAncestorsMatrixRcpp(uint numTips, uint numNodes, IntegerVector branchMatchIndexVec, NumericVector edgeLengthsVec, IntegerVector parentNumVec);
RcppExport SEXP _CovidCluster_produceDistTipsAncestorsMatrixRcpp(SEXP numTipsSEXP, SEXP numNodesSEXP, SEXP branchMatchIndexVecSEXP, SEXP edgeLengthsVecSEXP, SEXP parentNumVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint >::type numTips(numTipsSEXP);
    Rcpp::traits::input_parameter< uint >::type numNodes(numNodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type branchMatchIndexVec(branchMatchIndexVecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type edgeLengthsVec(edgeLengthsVecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type parentNumVec(parentNumVecSEXP);
    rcpp_result_gen = Rcpp::wrap(produceDistTipsAncestorsMatrixRcpp(numTips, numNodes, branchMatchIndexVec, edgeLengthsVec, parentNumVec));
    return rcpp_result_gen;
END_RCPP
}
// getMRCA_Rcpp
uint getMRCA_Rcpp(IntegerVector parentNumVec, IntegerVector tipsNumVec, uint numTips);
RcppExport SEXP _CovidCluster_getMRCA_Rcpp(SEXP parentNumVecSEXP, SEXP tipsNumVecSEXP, SEXP numTipsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type parentNumVec(parentNumVecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tipsNumVec(tipsNumVecSEXP);
    Rcpp::traits::input_parameter< uint >::type numTips(numTipsSEXP);
    rcpp_result_gen = Rcpp::wrap(getMRCA_Rcpp(parentNumVec, tipsNumVec, numTips));
    return rcpp_result_gen;
END_RCPP
}
// getMRCAclustersRcpp
Rcpp::List getMRCAclustersRcpp(IntegerVector& parentNumVec, List& childrenNumList, List& descendedTipsList, IntegerVector& subtreeIndexVec, CharacterVector& vertexRegionVec, CharacterVector& tipNamesVec, uint& subtreeRootNum, NumericMatrix& distTipsAncestorsMatrix, int subtreeIndex, int numTips, std::string regionLabel, int distLimit);
RcppExport SEXP _CovidCluster_getMRCAclustersRcpp(SEXP parentNumVecSEXP, SEXP childrenNumListSEXP, SEXP descendedTipsListSEXP, SEXP subtreeIndexVecSEXP, SEXP vertexRegionVecSEXP, SEXP tipNamesVecSEXP, SEXP subtreeRootNumSEXP, SEXP distTipsAncestorsMatrixSEXP, SEXP subtreeIndexSEXP, SEXP numTipsSEXP, SEXP regionLabelSEXP, SEXP distLimitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type parentNumVec(parentNumVecSEXP);
    Rcpp::traits::input_parameter< List& >::type childrenNumList(childrenNumListSEXP);
    Rcpp::traits::input_parameter< List& >::type descendedTipsList(descendedTipsListSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type subtreeIndexVec(subtreeIndexVecSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type vertexRegionVec(vertexRegionVecSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type tipNamesVec(tipNamesVecSEXP);
    Rcpp::traits::input_parameter< uint& >::type subtreeRootNum(subtreeRootNumSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distTipsAncestorsMatrix(distTipsAncestorsMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type subtreeIndex(subtreeIndexSEXP);
    Rcpp::traits::input_parameter< int >::type numTips(numTipsSEXP);
    Rcpp::traits::input_parameter< std::string >::type regionLabel(regionLabelSEXP);
    Rcpp::traits::input_parameter< int >::type distLimit(distLimitSEXP);
    rcpp_result_gen = Rcpp::wrap(getMRCAclustersRcpp(parentNumVec, childrenNumList, descendedTipsList, subtreeIndexVec, vertexRegionVec, tipNamesVec, subtreeRootNum, distTipsAncestorsMatrix, subtreeIndex, numTips, regionLabel, distLimit));
    return rcpp_result_gen;
END_RCPP
}
// getCopheneticClustersRcpp
Rcpp::List getCopheneticClustersRcpp(IntegerVector& parentNumVec, List& childrenNumList, List& descendedTipsList, IntegerVector& subtreeIndexVec, CharacterVector& vertexRegionVec, CharacterVector& tipNamesVec, uint& subtreeRootNum, NumericMatrix& distTipsAncestorsMatrix, int subtreeIndex, int numTips, std::string regionLabel, int distLimit);
RcppExport SEXP _CovidCluster_getCopheneticClustersRcpp(SEXP parentNumVecSEXP, SEXP childrenNumListSEXP, SEXP descendedTipsListSEXP, SEXP subtreeIndexVecSEXP, SEXP vertexRegionVecSEXP, SEXP tipNamesVecSEXP, SEXP subtreeRootNumSEXP, SEXP distTipsAncestorsMatrixSEXP, SEXP subtreeIndexSEXP, SEXP numTipsSEXP, SEXP regionLabelSEXP, SEXP distLimitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type parentNumVec(parentNumVecSEXP);
    Rcpp::traits::input_parameter< List& >::type childrenNumList(childrenNumListSEXP);
    Rcpp::traits::input_parameter< List& >::type descendedTipsList(descendedTipsListSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type subtreeIndexVec(subtreeIndexVecSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type vertexRegionVec(vertexRegionVecSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type tipNamesVec(tipNamesVecSEXP);
    Rcpp::traits::input_parameter< uint& >::type subtreeRootNum(subtreeRootNumSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distTipsAncestorsMatrix(distTipsAncestorsMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type subtreeIndex(subtreeIndexSEXP);
    Rcpp::traits::input_parameter< int >::type numTips(numTipsSEXP);
    Rcpp::traits::input_parameter< std::string >::type regionLabel(regionLabelSEXP);
    Rcpp::traits::input_parameter< int >::type distLimit(distLimitSEXP);
    rcpp_result_gen = Rcpp::wrap(getCopheneticClustersRcpp(parentNumVec, childrenNumList, descendedTipsList, subtreeIndexVec, vertexRegionVec, tipNamesVec, subtreeRootNum, distTipsAncestorsMatrix, subtreeIndex, numTips, regionLabel, distLimit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CovidCluster_produceDistTipsAncestorsMatrixRcpp", (DL_FUNC) &_CovidCluster_produceDistTipsAncestorsMatrixRcpp, 5},
    {"_CovidCluster_getMRCA_Rcpp", (DL_FUNC) &_CovidCluster_getMRCA_Rcpp, 3},
    {"_CovidCluster_getMRCAclustersRcpp", (DL_FUNC) &_CovidCluster_getMRCAclustersRcpp, 12},
    {"_CovidCluster_getCopheneticClustersRcpp", (DL_FUNC) &_CovidCluster_getCopheneticClustersRcpp, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_CovidCluster(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
