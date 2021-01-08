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

#' @title Power curve comparison
#'
#' @param data A list of data sets to be compared, the difference in the mean function is always computed as (f(data2) - f(data1)) 
#' @param xCol A numeric or vector stating column number of covariates
#' @param xCol.circ A numeric or vector stating column number of circular covariates
#' @param yCol A numeric value stating the column number of the response
#' @param testCol A numeric/vector stating column number of covariates to used in generating test set. Maximum of two columns to be used.
#' @param testSet A matrix or dataframe consisting of test points, default value NULL, if NULL computes test points internally using testCol variables. If not NULL, total number of test points must be less than or equal to 2500.
#' @param thrs A numeric or vector representing threshold for each covariates
#' @param conflevel A numeric between (0,1) representing the statistical significance level for constructing the band
#' @param gridSize A numeric / vector to be used in constructing test set, should be provided when testSet is NuLL, else it is ignored. Default is \code{c(50,50)} for 2-dim input which is converted internally to a default of \code{c(1000)} for 1-dim input. Total number of test points (product of gridSize vector components) must be less than or equal to 2500. 
#' @param powerbins A numeric stating the number of power bins for computing the scaled difference, default is 15.
#' @param baseline An integer between 0 to 2, where 1 indicates to use power curve of first dataset as the base for metric calculation, 2 indicates to use the power curve of second dataset as the base, and 0 indicates to use the average of both power curves as the base. Default is set to 1.
#' @param limitMemory A boolean (True/False) indicating whether to limit the memory use or not. Default is true. If set to true, 5000 datapoints are randomly sampled from each dataset under comparison for inference
#' @param opt_method A string specifying the optimization method to be used for hyperparameter estimation. Current options are: \code{'L-BFGS-B'}, \code{'BFGS'}, and \code{'nlminb'}. Default is set to \code{'nlminb'}.
#' @param sampleSize A named list of two integer items: \code{optimSize} and \code{bandSize}, denoting the sample size for each dataset for hyperparameter optimization and confidence band computation, respectively, when \code{limitMemory = TRUE}. Default value is \code{list(optimSize = 500, bandSize = 5000)}. 
#' @param rngSeed Random seed for sampling data when \code{limitMemory = TRUE}. Default is 1.
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
#'   \item matchedData - a list of two matched datasets as generated by covariate matching
#' }
#' @examples 
#' 
#' data1 = data1[1:100, ]
#' data2 = data2[1:100, ]
#' data = list(data1, data2)
#' xCol = 2
#' xCol.circ = NULL
#' yCol = 7
#' testCol = 2
#' testSet = NULL
#' thrs = 0.2
#' confLevel = 0.95
#' gridSize = 20

#' function_comparison = ComparePCurve(data, xCol, xCol.circ, yCol,
#' testCol, testSet, thrs, confLevel, gridSize)
#' 
#' @references For details, see Ding et al. (2020) available on \code{arxiv} at <\url{https://arxiv.org/abs/2005.08652}>.
#' @export

ComparePCurve = function(data, xCol, xCol.circ = NULL, yCol, testCol, testSet = NULL, thrs = 0.2, conflevel = 0.95, gridSize = c(50, 50), powerbins = 15, baseline = 1, limitMemory = TRUE, opt_method = 'nlminb',sampleSize = list(optimSize = 500, bandSize = 5000), rngSeed = 1 ){
  
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
  
  if (opt_method != "L-BFGS-B" && opt_method != "BFGS" && opt_method != "nlminb"){
    stop("opt_method must be 'L-BFGS-B', 'BFGS', or 'nlminb'.")
  }
  
  if(!is.list(data)){
    
    stop('The data must be provided as a list containing data sets')
    
  }
  
  if(length(data) != 2){
    
    
    stop('The data length should be equal to two')
    
  }
  
  if(!is.numeric(xCol)){
    
    stop('The xcol.circ must be provided as a numeric/vector')
    
  }else if(!all(xCol %in% 1:ncol(data[[1]]))){
    
    stop('The xCol values should be the column number of data set')
  }
  
  
  if(!is.null(xCol.circ)){
    
    if(!is.numeric(xCol.circ)){
      
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
  
  if(!is.numeric(testCol)){
    
    stop('The testCol must be provided as a numeric')
    
  }else if(length(testCol) > 2){
    
    stop('The length of testcol vector should be less than or equal to two')
    
  }
  
  if(!is.numeric(conflevel)){
    
    stop('conflevel must be provided as a numeric')
    
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
  
  if(is.null(testSet)){
    
    if(!is.numeric(gridSize)){
      
      stop('The gridsize must be provided as a numeric/vector')
      
    } else if(length(gridSize) != 2 && length(testCol) == 2){
      
      stop('The length of gridSize vector should be equal to two when length of testCol is equal to two')
      
    }
    
    if(length(testCol) == 1 && length(gridSize) == 2 && all(gridSize == c(50,50))){
      #Convert the 2-dim default gridSize to 1-dim default internally when length(testCol) == 1.
      gridSize = 1000
      
    } else if (length(testCol) == 1 && length(gridSize) != 1){
      
      stop('The length of gridSize vector should be equal to one when length of testCol is equal to one, or use the default gridSize option')
    }
    
    if (prod(gridSize) > 2500){
      
      stop('The number of test points should be less than or equal to 2500; reduce gridsize')
      
    }
    
  }else if(!is.matrix(testSet) & !is.data.frame(testSet)){
    
    stop('The test set provided should be a matrix or a data frame')
    
  }else if (length(testCol) != ncol(testSet)){
    
    stop('The length of testCol should be equal to the number of columns in testSet')
  }else if(nrow(testSet) > 2500){
    
       stop('The number of test points should be less than 2500')
  }
  
  if(is.null(testSet)){
    
    resultMatching = CovMatch(data, xCol, xCol.circ, thrs)
    
    testSet = GenerateTestset(resultMatching$matchedData, testCol, gridSize )
    
  }else{
    
    resultMatching = CovMatch(data, xCol, xCol.circ, thrs)
  }
  
  resultGP = funGP(resultMatching$matchedData, testCol, yCol, conflevel, testSet, limitMemory, opt_method)
  
  weightedDiff = ComputeWeightedDiff(data, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline)
  
  weightedStatDiff = ComputeWeightedStatDiff(data, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline)
  
  scaledDiff = ComputeScaledDiff(data, yCol, resultGP$mu1, resultGP$mu2, powerbins, baseline)
  
  scaledStatDiff = ComputeScaledStatDiff(data, yCol, resultGP$mu1, resultGP$mu2, resultGP$band, powerbins, baseline)
  
  unweightedDiff = ComputeDiff(resultGP$mu1, resultGP$mu2, baseline)
  
  unweightedStatDiff = ComputeStatDiff(resultGP$mu1, resultGP$mu2, resultGP$band, baseline)
  
  reductionRatio = ComputeRatio(data, resultMatching$matchedData, testCol)
  
  returnList = list(weightedDiff = weightedDiff, weightedStatDiff = weightedStatDiff, scaledDiff = scaledDiff, scaledStatDiff = scaledStatDiff, unweightedDiff = unweightedDiff, unweightedStatDiff = unweightedStatDiff, reductionRatio = reductionRatio, muDiff = resultGP$muDiff, mu2 = resultGP$mu2, mu1 = resultGP$mu1, band = resultGP$band, confLevel = conflevel, testSet = testSet, estimatedParams = resultGP$estimatedParams, matchedData = resultMatching$matchedData)
  
  return(returnList)
}
