#' @title KNN : Predict
#' @description The function can be used to make prediction on test data using trained model
#' @param KnnMdl a list containing:
#'               1) data
#'               2) column number of features
#'               3) column number of target
#'               4) bestK from trained model
#'
#' @param testData a data frame or matrix, to compute the predictions
#'
#' @return a list containing :
#'         1) data
#'         2) column number of features
#'         3) column number of target
#'         4) prediction on test data
#' @export
#'
KnnPredict = function(knnMdl, testData){

  trainData = knnMdl$data
  xCol = knnMdl$xCol
  yCol = knnMdl$yCol
  bestK = knnMdl$bestK

  prediction = knn.reg(trainData[, xCol], testData[, xCol], trainData[, yCol], k = bestK )

  return(list(data = data, xCol = result$bestSubset, yCol = yCol, prediction = prediction))
}
