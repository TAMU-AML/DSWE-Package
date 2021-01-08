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

#' @title SVM based power curve modelling
#'
#' @param trainX a matrix or dataframe to be used in modelling
#' @param trainY a numeric or vector as a target
#' @param testX a matrix or dataframe, to be used in computing the predictions
#' @param kernel default is 'radial' else can be 'linear', 'polynomial' and 'sigmoid'
#' @return a vector or numeric predictions on user provided test data
#'
#' @importFrom e1071 svm
#' @examples 
#' 
#' data = data1
#' trainX = as.matrix(data[c(1:100),2])
#' trainY = data[c(1:100),7]
#' testX = as.matrix(data[c(101:110),2])
#' 
#' Svm_prediction = SvmPCFit(trainX, trainY, testX)
#' 
#' @export

SvmPCFit = function(trainX, trainY, testX, kernel = 'radial'){
  
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
  
  if(!is.character(kernel)){
    
    stop("kernel can only take string input")
    
  }else if(kernel != 'radial' && kernel != 'linear' && kernel != 'polynomial' && kernel != 'sigmoid'){
    
    stop("kernel can only take followings as input: linear, radial, polynomial and sigmoid")
    
  }
  
  modelFit = svm(x = trainX, y = trainY, kernel = kernel)
  
  testPred = predict(modelFit, testX)
  
  return(testPred)
}
