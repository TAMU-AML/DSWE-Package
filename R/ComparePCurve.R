#' @title Power curve comparison
#'
#' @param data A list of data sets to be compared, the difference in the mean function is always computed as (f(data2) - f(data1)) 
#' @param xCol A numeric or vector stating column number of covariates
#' @param xCol.circ A numeric or vector stating column number of circular covariates
#' @param yCol A numeric value stating the column number of the response
#' @param testCol A numeric/vector stating column number of covariates to used in generating test set
#' @param testSet A matrix or dataframe consisting of test points, default value NULL, if NULL computes test points internally using testCol variables
#' @param thrs A numeric or vector representing threshold for each covariates
#' @param conflevel A numeric between (0,1) representing the statistical significance level for constructing the band
#' @param gridSize A numeric / vector to be used in constructing test set, should be provided when testSet is NuLL, else it is ignored
#' @param powerbins A numeric stating the number of power bins for computing the scaled difference, default is 15.
#' @param baseline An integer between 0 to 2, where 1 indicates to use power curve of first dataset as the base for metric calculation, 2 indicates to use the power curve of second dataset as the base, and 0 indicates to use the average of both power curves as the base. Default is set to 1.
#' @param limitMemory A boolean (True/False) indicating whether to limit the memory use or not. Default is true. If set to true, 5000 datapoints are randomly sampled from each dataset under comparison for inference
#'
#' @return a list containing :
#'  \itemize{
#'   \item weightedDiff - a numeric,  \% difference between the functions weighted using the density of the covariates
#'   \item weightedStatDiff - a numeric, \% statistically significant difference between the functions weighted using the density of the covariates
#'   \item scaledDiff -  a numeric, \% difference between the functions scaled to the orginal data
#'   \item scaledStatDiff - a numeric, \% statistically significant difference between the functions scaled to the orginal data
#'   \item unweightedDiff - a numeric,  \% difference between the functions unweighted
#'   \item unweightedStatDiff - a numeric,  \% statistically significant difference between the functions unweighted
#'   \item reductionRatio -  a list consisting of shrinkage ratio of features used in testSet
#'   \item mu1 - a vector of prediction on testset using the first data set
#'   \item mu2 - a vector of prediction on testset using the second data set
#'    \item muDiff - a vector of the difference in prediction (mu2 - mu1) for each test point
#'   \item band - a vector for the confidence band at all the testpoints for the two functions to be the same at a given cofidence level.
#'   \item confLevel - a numeric representing the statistical significance level for constructing the band
#'   \item testSet - a vector/matrix of the test points either provided by user, or generated internally
#'   \item estimatedParams - a list of estimated hyperaparameters for the Gaussian process model
#' }
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @export

ComparePCurve = function(data, xCol, xCol.circ = NULL, yCol, testCol, testSet = NULL, thrs = 0.2, conflevel = 0.95, gridSize = c(50, 50), powerbins = 15, baseline = 1, limitMemory = T ){

  if (class(limitMemory)!="logical"){
    stop('limitMemory should either be TRUE or FALSE')
  }

  if(!is.list(data)){

    stop('The data must be provided as a list containing data sets')

  }

  if(length(data) != 2){


    stop('The data length should be equal to two')

  }

  if(!is.vector(xCol)){

    stop('The xcol.circ must be provided as a numeric/vector')

  }else if(!all(xCol %in% 1:ncol(data[[1]]))){

    stop('The xCol values should be the column number of data set')
  }


  if(!is.null(xCol.circ)){

    if(!is.vector(xCol.circ)){

      stop('The xCol must be provided as a numeric/vector')

    }else if(!all(xCol.circ %in% xCol)){

      stop('xCol.circ should be present in xCol')
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

  if(!is.vector(conflevel)){

    stop('conflevel must be provided as a numeric/vector')

  }else{

    if(length(conflevel) != 1){

      stop('conflevel must be provided as a single numeric value')

    }else if(!(conflevel > 0 & conflevel < 1)){

      stop('conflevel must be between 0 to 1')
    }
  }

  if ((baseline %in% c(0:2)) == FALSE){
    stop('baseline must be an integer between 0 to 2')
  }

  resultMatching = CovMatch(data, xCol, xCol.circ, thrs)

  if(is.null(testSet)){

    if(!is.vector(gridSize)){

      stop('The gridsize must be provided as a vector')

    }else if(length(gridSize) != 2){

      stop('The length of gridSize vector should be of size two')

    }else if((gridSize[1] * gridSize[2] > 2500)){

      stop('The product of gridSize should not be more than 2500')

    }

    testSet = GenerateTestset(resultMatching$matchedData, testCol, gridSize )

  }else if(!is.matrix(testSet) & !is.data.frame(testSet)){

    stop('The test set provided should be a matrix or a dataframe')

  }else if (ncol(testCol) != ncol(testSet)){

    stop('The length of testCol should be equal to the number of columns in testSet')
  }

  resultGP = funGP(resultMatching$matchedData, testCol, yCol, conflevel, testSet, limitMemory)

  weightedDiff = ComputeWeightedDiff(data, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline)

  weightedStatDiff = ComputeWeightedStatDiff(data, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline)

  scaledDiff = ComputeScaledDiff(data, yCol, resultGP$mu1, resultGP$mu2, powerbins, baseline)

  scaledStatDiff = ComputeScaledStatDiff(data, yCol, resultGP$mu1, resultGP$mu2, resultGP$band, powerbins, baseline)

  unweightedDiff = ComputeDiff(resultGP$mu1, resultGP$mu2, baseline)

  unweightedStatDiff = ComputeStatDiff(resultGP$mu1, resultGP$mu2, resultGP$band, baseline)

  reductionRatio = ComputeRatio(data, resultMatching$matchedData, testCol)

  returnList = list(weightedDiff = weightedDiff, weightedStatDiff = weightedStatDiff, scaledDiff = scaledDiff, scaledStatDiff = scaledStatDiff, unweightedDiff = unweightedDiff, unweightedStatDiff = unweightedStatDiff, reductionRatio = reductionRatio, muDiff = resultGP$muDiff, mu2 = resultGP$mu2, mu1 = resultGP$mu1, band = resultGP$band, confLevel = conflevel, testSet = testSet, estimatedParams = resultGP$estimatedParams)

  return(returnList)
}
