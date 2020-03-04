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

ComparePCurve = function(data, xcol, xcol.circ = NULL, ycol, testCol, testSet = NULL, thrs = 0.2, conflevel = 0.95, gridSize = c(50, 50)){

  if(!is.list(data)){

    stop('Data set provided should be a list containing data sets')

  }

  if(length(data) != 2){


    stop('The number of data sets to match should be equal to two')

  }

  if(!is.vector(xCol)){

      stop('Non circular covariates column number should be provided as a vector')

    }

  if(length(xCol.circ) > 0){

    if(!is.vector(xCol.circ)){

      stop('Circular covariates column number should be provided as a vector')

    }
  }

  if(length(thrs) > 1){

    if(!(length(thrs) == length(xCol))){

      stop('The weight provided should be a single value or vector with weight for each covariate')

    }
  }

  if(!is.vector(testCol)){

    stop('The testCol should be provided as a vector')

  }else if(length(testCol) != 2){

    stop('The length of testcol vector should be of size two')
  }

  if(!is.vector(gridSize)){

    stop('The gridsize should be provided as a vector')

  }else if(length(gridSize) != 2){

    stop('The length of gridSize vector should be of size two')

  }else if((gridSize[1] * gridSize[2] > 2500)){

    stop('The product of gridSize should not be more than 2500')
  }

  resultMatching = CovMatch(data, xcol, xcol.circ, thrs)

  if(is.null(testSet)){

  testSet = GeneratetestSet(resultMatching, testCol, gridSize )

  }

  resultGP = funGP(resultMatching, testCol, ycol, conflevel, testSet)

  resultSMetric = ComputeSMetric(resultGP$mu1, resultGP$mu2, resultGP$band)

  resultWMetric = ComputeWMetric(data, resultGP$mu1, resultGP$mu2, testSet, var1col, var2col)

  reductionRatio = ComputeRatio(data, resultMatching$matchedData, var1col, var2col)

  extrapolatedPower = ComputeExtrapolation(data, ycol, resultGP$mu1, resultGP$mu2)

  returnList = list(weightedDiff = resultWMetric, statisticalDiff = resultSMetric, ratioVarcol1 = reductionRatio$ratioVar1, ratioVarcol2 = reductionRatio$ratioVar2, extrapolatedPower = extrapolatedPower, muDiff = resultGP$muDiff, mu2 = resultGP$diffCov$mu2, mu1 = resultGP$diffCov$mu1, band = resultGP$band, confLevel = confLevel, testSet = testSet, estimatedParams = resultGP$params)

  return(returnList)
}
