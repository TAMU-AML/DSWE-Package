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

#' @title Tree based power curve estimate
#'
#' @param trainX a matrix or dataframe to be used in modelling
#' @param trainY a numeric or vector as a target
#' @param testX a matrix or dataframe, to be used in computing the predictions
#' @param nTree a numeric value specifying number of trees to be constructed in model 
#'
#' @return a vector or numeric predictions on user provided test data
#'
#' @importFrom BayesTree bart
#' @importFrom stats predict
#' @examples 
#' 
#' data = data1
#' trainX = as.matrix(data[c(1:100),2])
#' trainY = data[c(1:100),7]
#' testX = as.matrix(data[c(100:110),2])
#' 
#' Bart_prediction = BayesTreePCFit(trainX, trainY, testX)
#' 
#' @export


BayesTreePCFit = function(trainX, trainY, testX, nTree = 50){
  
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
  
  modelFit = bart(x.train = trainX, y.train = trainY, x.test = testX, ntree = nTree, verbose = FALSE)
  
  return(modelFit$yhat.test.mean)
}
