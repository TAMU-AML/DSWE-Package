# MIT License
# 
# Copyright (c) 2020 Abhinav Prakash, Rui Tuo, and Yu Ding
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

#' @title temporal Gaussian process
#' @description A Gaussian process based power curve model which explicitly models the temporal aspect of the power curve.
#' @param trainX A matrix with each column corresponding to one input variable. 
#' @param trainY A vector with each element corresponding to the output at the corresponding row of \code{trainX}.
#'
#' @return An object of class \code{tempGP} with the following attributes:
#' \itemize{
#' \item trainX
#' \item trainY
#' \item thinningNumber
#' \item estimatedParams
#' \item llval
#' \item gradval
#' }
#' 
#' @importFrom stats pacf
#' @export
#' 

tempGP = function(trainX, trainY){
  
  trainX = as.matrix((trainX)) #trying to coerce trainX into a matrix, if not already.
  
  if (!is.matrix(trainX)){
    stop('trainX must be a matrix.')
  }
  
  if (!is.numeric(trainX)){
    stop('trainX must be a numeric matrix.')
  }
  
  if (any(!is.finite(trainX))){
    stop('trainX must have finite value values only.')
  }
  
  if (!is.numeric(trainY)){
    stop('trainY must be numeric.')
  }
  
  if(any(!is.finite(trainY))){
    stop('trainY must have finite values only.')
  }
  
  if(nrow(trainX) != length(trainY)){
    stop('Number of rows in trainX and the length of trainY must match.')
  }
  
  cat("All test passed.\n")
  thinningNumber = computeThinningNumber(trainX)
  cat("thinning number =",thinningNumber,'\n')
  if (thinningNumber > 0){
    thinnedBins = createThinnedBins(trainX,trainY,thinningNumber)  
  } else{
    thinnedBins = list(list(x = trainX, y = trainY))
  }
  cat("Estimating hyperparameters...\n")
  optimResult = estimateBinnedParams(thinnedBins)
  output = list(trainX = trainX, trainY = trainY, thinningNumber = thinningNumber, estimatedParams = optimResult$estimatedParams, llval = -optimResult$objVal, gradval = -optimResult$gradVal)
  class(output) = "tempGP"
  return(output)
}


#' @title predict from temporal Gaussian process
#' @description predict function for tempGP objects
#' @param tempGPObj An object of class tempGP
#' @param testX A matrix with each column corresponding to one input variable.
#'
#' @return a vector of predictions at the testpoints in \code{testX}.
#' @export
#' 
predict.tempGP = function(tempGPObj, testX){
  
  if (class(tempGPObj) != "tempGP"){
    stop("First argument (tempGPObj) must be an object of class tempGP.")
  }
  
  testX = as.matrix((testX)) #trying to coerce testX into a matrix, if not already.
  
  if (!is.matrix(testX)){
    stop('testX must be a matrix.')
  }
  
  if (!is.numeric(testX)){
    stop('testX must be a numeric matrix.')
  }
  
  if (any(!is.finite(testX))){
    stop('testX must have finite value values only.')
  }
  
  if (ncol(testX) != ncol(trainX)){
    stop("Number of columns in testX must be the same as that of the tempGP's trainX.")
  }
  
  pred = predictGP(tempGPObj$trainX, tempGPObj$trainY, testX, tempGPObj$estimatedParams)
  
  return(pred)
}




