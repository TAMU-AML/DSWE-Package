# MIT License
# 
# Copyright (c) 2024 Ahmadreza Chokhachian, and Yu Ding
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

#' @title xgboost based power curve modelling
#'
#' @param trainX a matrix or dataframe to be used in modelling
#' @param trainY a numeric or vector as a target
#' @param testX a matrix or dataframe, to be used in computing the predictions
#' @param max.depth  maximum depth of a tree
#' @param eta learning rate
#' @param nthread This parameter specifies the number of CPU threads that XGBoost
#' @param nrounds  number of boosting rounds or trees to build
#' @return a vector or numeric predictions on user provided test data
#'
#' @importFrom xgboost xgboost
#' @examples 
#' 
#' data = data1
#' trainX = as.matrix(data[c(1:100),2])
#' trainY = data[c(1:100),7]
#' testX = as.matrix(data[c(101:110),2])
#' 
#' Xgb_prediction = XgbPCFit(trainX, trainY, testX)
#' 
#' @references Chen, T., & Guestrin, C. (2016). "XGBoost: A Scalable Tree Boosting System." Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, 785-794. \doi{10.1145/2939672.2939785}.
#' @export

XgbPCFit = function(trainX, trainY, testX, max.depth = 8, eta = 0.25, nthread = 2, nrounds = 5){
  
  if (!is.matrix(trainX) && !is.data.frame(trainX)) {
    
    stop("trainX must be a matrix or a dataframe.")
    
  }
  
  if (!is.numeric(trainY)){
    
    stop("trainY must be numeric/vector.")
    
  }
  
  if (length(trainY) != nrow(trainX)){
    
    stop("number of datapoints in trainX and trainY must be the same.")
    
  }
  
  if (!is.matrix(testX) && !is.data.frame(testX)) {
    
    stop("testX must be a matrix or a dataframe.")
    
  }
  
  if(ncol(trainX) != ncol(testX)){
    
    stop("testX and trainX must have same number of columns")
  }

  if (!is.numeric(max.depth) || max.depth <= 0) {
    stop("max.depth must be a positive numeric value.")
  }
  
  if (!is.numeric(eta) || eta <= 0 || eta > 1) {
    stop("eta must be a numeric value between 0 and 1.")
  }
  
  if (!is.numeric(nthread) || nthread <= 0 || nthread != round(nthread)) {
    stop("nthread must be a positive integer.")
  }
  
  if (!is.numeric(nrounds) || nrounds <= 0 || nrounds != round(nrounds)) {
    stop("nrounds must be a positive integer.")
  }

  modelFit = xgboost::xgboost(data = trainX, label = trainY, max.depth = max.depth, eta = eta, nthread = nthread, nrounds = nrounds, verbose=FALSE)

  testPred = predict(modelFit, testX)
  
  return(testPred)
}
