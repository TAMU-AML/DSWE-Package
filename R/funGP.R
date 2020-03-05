#' @title Function estimator using Gaussian Process
#'
#' @param datalist a list of data sets to compute a function for each of them
#' @param covCols a numeric or vector stating the column number of covariates
#' @param yCol A numeric value stating the column number of target
#' @param confLevel a statistical value, such as 0.95, to be used for hypothesis testing by computing the band
#' @param testset Test points at which the functions will be compared
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
#' }
#' @export
funGP = function(datalist, covCols, yCol, confLevel = 0.95, testset ){

  if(!is.list(datalist)){

    stop('datalist must be a list containing data sets')

  }

  if(length(datalist) != 2){


    stop('The number of data sets to match should be equal to two')

  }

  if(!is.vector(covCols)){

    stop('covCols must be provided as a numeric/vector')

  }

  if(!is.vector(yCol)){

    stop('covCols must be provided as a numeric/vector')

  }else{

    if(length(yCol) != 1){

      stop('yCol must be provided as a single numeric value')
    }
  }

  if(!is.vector(confLevel)){

    stop('confLevel must be provided as a numeric/vector')

  }else{

    if(length(confLevel) != 1){

      stop('confLevel must be provided as a single numeric value')

    }else if(!(confLevel > 0 & confLevel < 1)){

        stop('confLevel must be between 0 to 1')
      }
    }



  params = estimateParameters(datalist, covCols, yCol)$estimatedParams

  diffCov = computeDiffCov(datalist, covCols, yCol, params, testset)

  muDiff = diffCov$mu2 - diffCov$mu1

  diffCovMat = diffCov$diffCovMat

  band = computeConfBand(diffCovMat, confLevel)

  returnList = list(muDiff = muDiff,mu2 = diffCov$mu2, mu1= diffCov$mu1,band = band, confLevel = confLevel, testset = testset, estimatedParams = params)

  return(returnList)
}
