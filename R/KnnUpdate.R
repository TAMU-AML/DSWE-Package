#' @title KNN : Update
#' @description The function can be used to update KNN model when new data is provided
#' @param KnnMdl a list containing:
#'               1) data
#'               2) column number of features
#'               3) column number of target
#'               4) bestK from trained model
#' @param newData a dataframe or a matrix, to be used for updating the model
#'
#' @return a list containing :
#'         1) data
#'         2) column number of features
#'         3) column number of target
#'         4) best K
#' @export
updateKNN = function(knnMdl, newData){
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
  knnMdl = list(data = data, bestK = bestK, xCol = xCol, yCol = yCol)
  return(knnMdl)
}
