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

#' @title Smoothing spline Anova method
#'
#' @param data a matrix or dataframe to be used in modelling
#' @param xCol a numeric or vector stating the column number of feature covariates
#' @param yCol a numeric value stating the column number of target
#' @param testX a matrix or dataframe, to be used in computing the predictions
#' @param modelFormula default is NULL else a model formula specifying target and features.Please refer 'gss' package documentation for more details
#' 
#' @return a vector or numeric predictions on user provided test data
#'
#' @importFrom gss ssanova
#' @examples 
#' 
#' data = data1[c(1:100),]
#' xCol = 2
#' yCol = 7
#' testX = data1[c(101:110), ]

#' Spline_prediction = SplinePCFit(data, xCol, yCol, testX)
#' 
#' @export

SplinePCFit = function(data, xCol, yCol, testX, modelFormula = NULL){
  
  if(!is.matrix(data) && !is.data.frame(data)) {
    
    stop("data must be a matrix or a dataframe.")
  }
  
  if(!is.vector(xCol) & !is.numeric(xCol)){
    
    stop('xCol must be provided as numeric/vector')
    
  }else if(!(all(xCol %in% 1:ncol(data)))){
    
    stop('xCol value should not be more than the column in data')
  }
  
  
  if(!is.vector(yCol) & !is.numeric(yCol)){
    
    stop('yCol must be provided as a numeric/vector')
    
  }else if(!(all(yCol %in% 1:ncol(data)))){
    
    stop('yCol value should not be more than the column in data')
    
  }else if(length(yCol) != 1){
    
    stop('yCol must be provided as a single numeric value')
  }
  
  if (!is.matrix(testX) && !is.data.frame(testX)) {
    
    stop("testX must be a matrix or a dataframe.")
  }else if(ncol(data) != ncol(testX)){
    
    stop("testX and data should have same number of columns")
  }
  
  # adjusting test points which are outliers
  for(col in xCol){  
    
    upperRange = max(data[, col]) * 1.04
    lowerRange = min(data[, col]) * 0.96
    testX[testX[, col] > upperRange, col] = upperRange
    testX[testX[, col] < lowerRange, col] = lowerRange
    
  }
  
  if(is.null(modelFormula)){
    
    # manipulating data and test set column names
    colnames(data) = paste('col', 1:ncol(data), sep = '')
    colnames(testX) = paste('col', 1:ncol(testX), sep = '')
    
    # preparing x and y for the formula
    xCombined = (paste('col', xCol, sep = '', collapse = '+'))
    yTarget = paste('col', yCol, sep = '')
    
    # model formula
    modelFormula = paste(yTarget, xCombined, sep = '~')
    
    # modified xCol and yCol
    xCol = (paste('col', xCol, sep = ''))
    yCol = (paste('col', yCol, sep = ''))
  }
  
  #model fitting
  modelFit = ssanova(data = data, formula = modelFormula, skip.iter = FALSE)
  
  #predcition on test points
  testPred = predict(modelFit, testX[, xCol, drop = FALSE])
  return(testPred)
}
