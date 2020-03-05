# Computes best k using generalized cross validation
computeBestK = function(dataX, dataY, rangeK ){

  bestK = NULL
  dataX = as.matrix(dataX)
  maxK = max(rangeK)
  nnIdx = knnx.index(dataX, query = dataX, k = maxK)
  gcv = rep(0,length(rangeK))
  gcv_mae = rep(0, length(rangeK))
  for (i in 1:length(rangeK)){
    predY = rowMeans(matrix(dataY[nnIdx[, 1:rangeK[i]]], ncol = ncol(nnIdx[, 1:rangeK[i]])))
    gcv[i] = sqrt(mean(((dataY - predY) / (1 - (1 / rangeK[i])))^2))
    gcv_mae[i] = mean(dataY - predY)
  }
  bestK = rangeK[which.min(gcv)]
  bestRMSE = min(gcv)
  bestMAE = min(gcv_mae)
  returnList = list(bestK = bestK, bestRMSE = bestRMSE, bestMAE = bestMAE)
  if (bestK == maxK){
    rangeK = maxK + seq(5,50,5)
    returnList = computeBestK(dataX, dataY, rangeK)
  }
  return(returnList)

}


# Computes best Subset from given features
computeBestSubset = function(data, xCol, yCol,rangeK){

  bestSubset = NULL
  bestRMSE = Inf
  bestK = NULL

  .computeBestSubset = function(data, xCol, yCol, rangeK, bestSubset, bestRMSE, bestK){
    nCov = length(xCol)
    bestCol = NULL
    for (i in 1:nCov){
      result = computeBestK(data[, c(bestSubset, xCol[i])], data[, yCol], rangeK)
      RMSE = result$bestRMSE
      if (RMSE < bestRMSE){
        bestRMSE = RMSE
        bestK = result$bestK
        bestCol = xCol[i]
      }
    }

    returnList = list(bestSubset = bestSubset, bestK = bestK, bestRMSE = bestRMSE )

    if (length(bestCol)>0){
      bestSubset = c(bestSubset, bestCol)
      xColDiff = setdiff(xCol, bestSubset)
      returnList = .computeBestSubset(data, xColDiff, yCol, rangeK, bestSubset, bestRMSE, bestK)
    }

    return(returnList)
  }

  returnList = .computeBestSubset(data, xCol, yCol, rangeK, bestSubset, bestRMSE)
  return(returnList)

}
