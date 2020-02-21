#' @title Power curve comparison
#'
#' @param data List of data sets to be compared
#' @param xcol Vector of covariates column number
#' @param xcol.circ vector of circular covariates column number
#' @param ycol A value representing column number of target
#' @param var1col Column number of first covariate to used in generating test set
#' @param var2col Column number of second covariate to be used in generating test set
#' @param thrs A single value or vector represnting threshold for each covariates
#' @param conflevel A single value as a Bound to be used in hypothesis testing while function comparison
#' @param gridSize A single value to used in constructing test set size
#'
#' @export

ComputePCurve = function(data, xcol, xcol.circ = NULL, ycol, var1col, var2col, thrs = 0.2, conflevel = 0.95, gridSize = 50){

  ResultMatching = CovMatch(data, xcol, xcol.circ, thrs)

  TestSet = GenerateTestset(ResultMatching, var1col, var2col, gridSize )

  ResultGP = funGP(ResultMatching, c(var1col, var2col), ycol, conflevel, TestSet)

  ResultSMetric = ComputeSMetric(ResultGP$mu1, ResultGP$mu2, ResultGP$band)

  ResultWMetric = ComputeWMetric(data, ResultMatching, ResultGP$mu1, ResultGP$mu2, TestSet, var1col, var2col)

  returnList = list(muDiff = ResultGP$muDiff, mu2 = ResultGP$diffCov$mu2, mu1 = ResultGP$diffCov$mu1, band = ResultGP$band, confLevel = confLevel, testset = testset, estimatedParams = ResultGP$params, weightedDiff = ResultWMetric, statisticalDiff = ResultSMetric)

  return(returnList)
}
