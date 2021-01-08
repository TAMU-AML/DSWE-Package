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
#' @param datalist A list of data sets to compute a function for each of them
#' @param xCol A numeric or vector stating the column number of covariates
#' @param yCol A numeric value stating the column number of target
#' @param confLevel A single value representing the statistical significance level for constructing the band
#' @param testset Test points at which the functions will be compared
#' @param limitMemory A boolean (True/False) indicating whether to limit the memory use or not. Default is true. If set to true, 5000 datapoints are randomly sampled from each dataset under comparison for inference.
#' @param opt_method A string specifying the optimization method to be used for hyperparameter estimation. Current options are: \code{'L-BFGS-B'}, \code{'BFGS'}, and \code{'nlminb'}. Default is set to \code{'nlminb'}.
#' @param sampleSize A named list of two integer items: \code{optimSize} and \code{bandSize}, denoting the sample size for each dataset for hyperparameter optimization and confidence band computation, respectively, when \code{limitMemory = TRUE}. Default value is \code{list(optimSize = 500, bandSize = 5000)}. 
#' @param rngSeed Random seed for sampling data when \code{limitMemory = TRUE}. Default is 1.
#' 
#' @return a list containing :
#'  \itemize{
#'   \item muDiff - A vector of pointwise difference between the predictions from the two datasets \code{(mu2- mu1)}
#'   \item mu1 - A vector of test prediction for first data set
#'   \item mu2 - A vector of test prediction for second data set
#'   \item band - A vector of the allowed statistical difference between functions at testpoints in testset
#'   \item confLevel - A numeric representing the statistical significance level for constructing the band
#'   \item testset - A matrix of test points to compare the functions
#'   \item estimatedParams - A list of estimated hyperparameters for GP
#' }
#' @examples 
#' 
#' datalist = list(data1[1:100,], data2[1:100, ])
#' xCol = 2
#' yCol = 7
#' confLevel = 0.95
#' testset = seq(4,10,length.out = 20)

#' function_diff = funGP(datalist, xCol, yCol, confLevel, testset)
#' 
#' @references Prakash, A., Tuo, R., & Ding, Y. (2020). "Gaussian process aided function comparison using noisy scattered data." arXiv preprint arXiv:2003.07899.  <\url{https://arxiv.org/abs/2003.07899}>.
#'@export
funGP = function(datalist, xCol, yCol, confLevel = 0.95, testset, limitMemory = TRUE, opt_method = 'nlminb', sampleSize = list(optimSize = 500, bandSize = 5000), rngSeed = 1){

  if (class(limitMemory)!="logical"){
    stop('limitMemory should either be TRUE or FALSE')
  }
  
  if (limitMemory){
    if(!is.list(sampleSize)){
      stop('If limitMemory is TRUE, sampleSize must be a list with two named items: optimSize and bandSize.')
    }
    if(length(sampleSize) != 2){
      stop('If limitMemory is TRUE, sampleSize must be a list with two named items: optimSize and bandSize.')
    }
    if(!all(names(sampleSize)%in%c("optimSize","bandSize"))){
      stop('If limitMemory is TRUE, sampleSize must be a list with two named items: optimSize and bandSize.')
    }
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

  if (opt_method != "L-BFGS-B" && opt_method != "BFGS" && opt_method != "nlminb"){
    stop("opt_method must be 'L-BFGS-B', 'BFGS', or 'nlminb'.")
  }

  params = estimateParameters(datalist, xCol, yCol, opt_method, limitMemory, sampleSize$optimSize, rngSeed)$estimatedParams

  diffCov = computeDiffCov(datalist, xCol, yCol, params, testset, limitMemory, sampleSize$bandSize, rngSeed)

  muDiff = diffCov$mu2 - diffCov$mu1

  band = computeConfBand(diffCov$diffCovMat, confLevel)

  returnList = list(muDiff = muDiff,mu2 = diffCov$mu2, mu1= diffCov$mu1,band = band, confLevel = confLevel, testset = testset, estimatedParams = params)

  return(returnList)
}
