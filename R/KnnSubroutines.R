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

# Computes best k using generalized cross validation
computeBestK = function(dataX, dataY, rangeK ){
  bestK = NULL
  dataX = as.matrix(dataX)
  maxK = max(rangeK)
  nnIdx = knnx.index(dataX, query = dataX, k = maxK)
  gcv = rep(0,length(rangeK))
  for (i in 1:length(rangeK)){
    predY = rowMeans(matrix(dataY[nnIdx[, 1:rangeK[i]]], ncol = ncol(nnIdx[, 1:rangeK[i]])))
    gcv[i] = sqrt(mean(((dataY - predY) / (1 - (1 / rangeK[i])))^2))
  }
  bestK = rangeK[which.min(gcv)]
  bestRMSE = min(gcv)
  returnList = list(bestK = bestK, bestRMSE = bestRMSE)
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

    returnList = list(bestSubset = bestSubset, bestK = bestK, bestRMSE = bestRMSE)

    if (length(bestCol)>0){
      bestSubset = c(bestSubset, bestCol)
      xColDiff = setdiff(xCol, bestSubset)
      if (length(xColDiff)>0){
        returnList = .computeBestSubset(data, xColDiff, yCol, rangeK, bestSubset, bestRMSE, bestK)
      }else {
        returnList = list(bestSubset = bestSubset, bestK = bestK, bestRMSE = bestRMSE)
      }
    }

    return(returnList)
  }

  returnList = .computeBestSubset(data, xCol, yCol, rangeK, bestSubset, bestRMSE)
  return(returnList)

}
