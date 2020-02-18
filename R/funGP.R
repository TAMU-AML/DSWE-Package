#' @title Function estimator using Gaussian Process
#'
#' @param datalist List of data sets
#' @param covCols Vector of column number of covariates
#' @param yCol A numeric value for target
#' @param confLevel Statistical confidence value, such as 0.95
#' @param testset Test points at which the functions will be compared
#'
#' @export
funGP = function(datalist, covCols, yCol, confLevel, testset ){

  params = estimateParameters(datalist, covCols, yCol)$estimatedParams

  diffCov = computeDiffCov(datalist, covCols, yCol, params, testset)

  muDiff = diffCov$mu2 - diffCov$mu1

  diffCovMat = diffCov$doffCovMat

  band = computeConfBand(diffCovMat, confLevel)

  returnList = list(muDiff = muDiff,mu2 = diffCov$mu2, mu1= diffCov$mu1,band = band, confLevel = confLevel, testset = testset, estimatedParams = params)

  return(returnList)
}
