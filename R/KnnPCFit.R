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

#' @title KNN : Fit
#' @description The function models the powercurve using KNN, against supplied arguments
#' @param data a dataframe or a matrix, to be used in modelling
#' @param xCol a vector or numeric values stating the column number of features
#' @param yCol a numerical or a vector value stating the column number of target
#' @param subsetSelection a boolean, default value is FALSE, if TRUE returns the best feature column number as xCol
#'
#' @return a list containing :
#'  \itemize{
#'   \item data - The data set provided by user
#'   \item xCol - The column number of features provided by user or the best subset column number
#'   \item yCol - The column number of target provided by user
#'   \item bestK - The best k nearest neighbor calculated using the function
#'   \item RMSE - The RMSE calculated using the function for provided data using user defined features and best obtained K
#'   \item MAE - The MAE calculated using the function for provided data using user defined features and best obtained K
#'}
#' @examples 
#' 
#' data = data1[c(1:100),]
#' xCol = 2
#' yCol = 7
#' subsetSelection = FALSE
#' 
#' knn_model = KnnPCFit(data, xCol, yCol, subsetSelection)
#' 
#' @export
#' @importFrom FNN knn.reg knnx.index
KnnPCFit = function(data, xCol, yCol, subsetSelection = FALSE){

  if(!is.matrix(data) & !is.data.frame(data)){

    stop('data provided should either be a matrix or data frame')
  }

  if(!is.numeric(xCol)  & !is.vector(xCol)){

    stop('column number of features should be provided as a numeric or vector')

  }else{

    if(sum(xCol %in% 1:ncol(data)) < length(xCol)){

      stop('column numbers provided are not in the data')
    }
  }

  if(!is.numeric(yCol)  & !is.vector(yCol)){

    stop('column number of target should be provided as a numeric or vector')

  }else{

    if(length(yCol) > 1){

      stop('a signle numeric or vector input should be provided for target')
    }
  }

  normalizedData = data

  for (feature in xCol) {

    normalizedData[, feature] = (data[, feature] - min(data[, feature])) / (max(data[, feature]) - min(data[, feature]))

  }

  rangeK = seq(5,50,5)

  if(subsetSelection == FALSE){

    result = computeBestK(normalizedData[, xCol, drop = FALSE], normalizedData[, yCol], rangeK)
    pred = FNN::knn.reg(normalizedData[, xCol, drop = FALSE], normalizedData[, xCol, drop = FALSE], normalizedData[, yCol], result$bestK)
    mae = mean(abs(normalizedData[, yCol] - pred$pred))
    returnList = list(bestK = result$bestK, RMSE = result$bestRMSE, MAE = mae, data = data, xCol = xCol, yCol = yCol)

  }else{

    result = computeBestSubset(normalizedData, xCol, yCol, rangeK)
    pred = FNN::knn.reg(normalizedData[, xCol, drop = FALSE], normalizedData[, xCol, drop = FALSE], normalizedData[, yCol], result$bestK)
    mae = mean(abs((normalizedData[, yCol] - pred$pred)/(1 - (1 / result$bestK))))
    returnList = list(bestK = result$bestK, RMSE = result$bestRMSE, MAE = mae, data = data, xCol = result$bestSubset, yCol = yCol )
  }

  return(returnList)
}
