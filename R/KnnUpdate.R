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
