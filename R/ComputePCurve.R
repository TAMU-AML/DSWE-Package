#' @title Power curve comparison
#'
#' @param data a List of data sets to be compared
#' @param xcol a numeric or vector stating column number of covariates
#' @param xcol.circ a numeric or vector stating column number of circular covariates
#' @param ycol A numeric value stating the column number of target
#' @param var1col a numeric value stating column number of first covariate to used in generating test set
#' @param var2col a numeric value stating column number of second covariate to be used in generating test set
#' @param thrs A single value or vector represnting threshold for each covariates
#' @param conflevel A single value as a Bound to be used in hypothesis testing while function comparison
#' @param gridSize A single value to used in constructing test set size
#'
#' @return a list containing :
#'  \itemize{
#'   \item muDiff - The testset prediction for each of the data set
#'   \item mu1 - The test prediction for first data set
#'   \item mu2 - The test prediction for second data set
#'   \item band - The allowed statistical difference between functions
#'   \item confLevel - The statistical boundation
#'   \item testset - The test prediction for first data set
#'   \item estimatedParams - The function parameter values
#'   \item statisticalDiff - The \% statistical difference between functions in comaprison
#'   \item weightedDiff - The \% wighted difference between functions in comparison
#'   \item ratioVarcol1 - The shrinkage ratio of varcol1
#'   \item ratioVarcol1 - The shrinkage ratio of varcol2
#'   \item extrapolatedPower - The matched power difference scaled on entire domain of a covariate
#' }
#'
#' @export

ComputePCurve = function(data, xcol, xcol.circ = NULL, ycol, var1col, var2col, thrs = 0.2, conflevel = 0.95, gridSize = 50){

  ResultMatching = CovMatch(data, xcol, xcol.circ, thrs)

  TestSet = GenerateTestset(ResultMatching, var1col, var2col, gridSize )

  ResultGP = funGP(ResultMatching, c(var1col, var2col), ycol, conflevel, TestSet)

  ResultSMetric = ComputeSMetric(ResultGP$mu1, ResultGP$mu2, ResultGP$band)

  ResultWMetric = ComputeWMetric(data, ResultGP$mu1, ResultGP$mu2, TestSet, var1col, var2col)

  ReductionRatio = ComputeRatio(data, ResultMatching$matchedData, var1col, var2col)

  ExtrapolatedPower = ComputeExtrapolation(data, ycol, ResultGP$mu1, ResultGP$mu2)

  returnList = list(muDiff = ResultGP$muDiff, mu2 = ResultGP$diffCov$mu2, mu1 = ResultGP$diffCov$mu1, band = ResultGP$band, confLevel = confLevel, testset = testset, estimatedParams = ResultGP$params, weightedDiff = ResultWMetric, statisticalDiff = ResultSMetric, ratioVarcol1 = ReductionRatio$ratioVar1, ratioVarcol2 = ReductionRatio$ratioVar2, extrapolatedPower = ExtrapolatedPower)

  return(returnList)
}
