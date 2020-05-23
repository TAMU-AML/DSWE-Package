# MIT License
# 
# Copyright (c) 2020 Nitesh Kumar, Abhinav Prakash, and Yu Ding
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' @title Function comparison using Gaussian Process and Hypothesis testing
#'
#' @param datalist a list of data sets to compute a function for each of them
#' @param xCol a numeric or vector stating the column number of covariates
#' @param yCol A numeric value stating the column number of target
#' @param confLevel a single value representing the statistical significance level for constructing the band
#' @param testset Test points at which the functions will be compared
#' @param limitMemory A boolean (True/False) indicating whether to limit the memory use or not. Default is true. If set to true, 5000 datapoints are randomly sampled from each dataset under comparison for inference.
#'
#' @return a list containing :
#'  \itemize{
#'   \item muDiff - A vector of pointwise difference between the predictions from the two datasets (mu1- mu2)
#'   \item mu1 - A vector of test prediction for first data set
#'   \item mu2 - A vector of test prediction for second data set
#'   \item band - A vector of the allowed statistical difference between functions at testpoints in testset
#'   \item confLevel - A numeric representing the statistical significance level for constructing the band
#'   \item testset - A matrix of test points to compare the functions
#'   \item estimatedParams - A list of estimated hyperparameters for GP
#' }
#'@export
funGP = function(datalist, xCol, yCol, confLevel = 0.95, testset, limitMemory = T){

  if (class(limitMemory)!="logical"){
    stop('limitMemory should either be TRUE or FALSE')
  }

  if(!is.list(datalist)){

    stop('datalist must be a list containing data sets')

  }

  if(length(datalist) != 2){


    stop('The number of data sets to match should be equal to two')

  }

  if(!is.vector(xCol)){

    stop('xCol must be provided as a numeric/vector')

  }

  if(!is.vector(yCol)){

    stop('xCol must be provided as a numeric/vector')

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



  params = estimateParameters(datalist, xCol, yCol)$estimatedParams

  diffCov = computeDiffCov(datalist, xCol, yCol, params, testset, limitMemory)

  muDiff = diffCov$mu2 - diffCov$mu1

  band = computeConfBand(diffCov$diffCovMat, confLevel)

  returnList = list(muDiff = muDiff,mu2 = diffCov$mu2, mu1= diffCov$mu1,band = band, confLevel = confLevel, testset = testset, estimatedParams = params)

  return(returnList)
}
