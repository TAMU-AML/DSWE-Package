#' @title Power curve comparison
#'
#' @param data a List of data sets to be compared
#' @param xCol a numeric or vector stating column number of covariates
#' @param xCol.circ a numeric or vector stating column number of circular covariates
#' @param yCol A numeric value stating the column number of target
#' @param testCol a numeric/vector stating column number of covariates to used in generating test set
#' @param testSet a matrix or dataframe consisting of test points
#' @param thrs A single value or vector represnting threshold for each covariates
#' @param conflevel A single value as a Bound to be used in hypothesis testing while function comparison
#' @param gridSize A single value to used in constructing test set size
#'
#' @return a list containing :
#'  \itemize{
#'   \item muDiff - The testSet prediction for each of the data set
#'   \item mu1 - The test prediction for first data set
#'   \item mu2 - The test prediction for second data set
#'   \item band - The allowed statistical difference between functions
#'   \item confLevel - The statistical boundation
#'   \item testSet - The test prediction for first data set
#'   \item estimatedParams - The function parameter values
#'   \item statisticalDiff - The \% statistical difference between functions in comaprison
#'   \item weightedDiff - The \% wighted difference between functions in comparison
#'   \item ratioVarcol1 - The shrinkage ratio of varcol1
#'   \item ratioVarcol1 - The shrinkage ratio of varcol2
#'   \item extrapolatedPower - The matched power difference scaled on entire domain of a covariate
#' }
#'
#' @import dplyr
#' @export

ComparePCurve = function(data, xCol, xCol.circ = NULL, yCol, testCol, testSet = NULL, thrs = 0.2, conflevel = 0.95, gridSize = c(50, 50)){

  if(!is.list(data)){

    stop('The data must be provided as a list containing data sets')

  }

  if(length(data) != 2){


    stop('The data length should be equal to two')

  }

  if(!is.vector(xCol)){

      stop('The xcol.circ must be provided as a numeric/vector')

    }

  if(length(xCol.circ) > 0){

    if(!is.vector(xCol.circ)){

      stop('The xCol must be provided as a numeric/vector')

    }
  }

  if(length(thrs) > 1){

    if(!(length(thrs) == length(xCol))){

      stop('The thrs must be provided as a single value or vector with weight for each covariate')

    }
  }

  if(!is.vector(testCol)){

    stop('The testCol must be provided as a numeric/vector')

  }else if(length(testCol) != 2){

    stop('The length of testcol vector should be of size two')
  }

  if(!is.vector(gridSize)){

    stop('The gridsize must be provided as a vector')

  }else if(length(gridSize) != 2){

    stop('The length of gridSize vector should be of size two')

  }else if((gridSize[1] * gridSize[2] > 2500)){

    stop('The product of gridSize should not be more than 2500')
  }

  resultMatching = CovMatch(data, xCol, xCol.circ, thrs)

  if(is.null(testSet)){

  testSet = GeneratetestSet(resultMatching, testCol, gridSize )

  }

  resultGP = funGP(resultMatching, testCol, yCol, conflevel, testSet)

  resultSMetric = ComputeSMetric(resultGP$mu1, resultGP$mu2, resultGP$band)

  resultWMetric = ComputeWMetric(data, resultGP$mu1, resultGP$mu2, testSet, testCol)

  reductionRatio = ComputeRatio(data, resultMatching$matchedData, testCol)

  extrapolatedPower = ComputeExtrapolation(data, yCol, resultGP$mu1, resultGP$mu2)

  returnList = list(weightedDiff = resultWMetric, statisticalDiff = resultSMetric, ratioVarcol1 = reductionRatio$ratioVar1, ratioVarcol2 = reductionRatio$ratioVar2, extrapolatedPower = extrapolatedPower, muDiff = resultGP$muDiff, mu2 = resultGP$diffCov$mu2, mu1 = resultGP$diffCov$mu1, band = resultGP$band, confLevel = confLevel, testSet = testSet, estimatedParams = resultGP$params)

  return(returnList)
}
