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

#' @title KNN : Update
#' @description The function can be used to update KNN model when new data is provided
#' @param knnMdl a list containing:
#'  \itemize{
#'   \item knnMdl$data - The data set provided by user
#'   \item knnMdl$xCol - The column number of features provided by user or the best subset column number
#'   \item knnMdl$yCol - The column number of target provided by user
#'   \item knn$bestK - The best k nearest neighbor calculated using the function KnnFit
#'}
#' @param newData a dataframe or a matrix, to be used for updating the model
#'
#' @return a list containing :
#'  \itemize{
#'   \item data - The updated data using old data set and new data
#'   \item xCol - The column number of features provided by user or the best subset column number
#'   \item yCol - The column number of target provided by user
#'   \item bestK - The best k nearest neighbor calculated for the new data using user specified features and target
#'}
#' @examples 
#' 
#' data = data1[c(1:100),]
#' xCol = 2
#' yCol = 7
#' subsetSelection = FALSE
#' 
#' knn_model = KnnPCFit(data, xCol, yCol, subsetSelection)
#' newData = data1[c(101:110), ]
#' 
#' knn_newmodel = KnnUpdate(knn_model, newData)
#' 
#' @export
#'
KnnUpdate = function(knnMdl, newData){
  data = knnMdl$data
  bestK = knnMdl$bestK
  yCol = knnMdl$yCol
  xCol = knnMdl$xCol
  data = data[-c(1:nrow(newData)),]
  data = rbind(data, newData)
  normalizedData = data
  for (i in 1:length(xCol)){
    normalizedData[,i] = (data[,i] - min(data[,i])) / (max(data[,i]) - min(data[,i]))
  }
  ubK = 1.2
  lbK = 0.8
  intervalK = 5
  maxK = ceiling(ubK*bestK)
  maxK = maxK + (intervalK - (maxK %% intervalK))
  minK = floor(lbK*bestK)
  minK = minK - (minK %% intervalK)
  rangeK = seq(minK,maxK,intervalK)
  dataX = normalizedData[,xCol]
  dataY = normalizedData[,yCol]
  bestK = computeBestK(dataX, dataY, rangeK)
  knnMdl = list(data = data, bestK = bestK$bestK, xCol = xCol, yCol = yCol)
  return(knnMdl)
}
