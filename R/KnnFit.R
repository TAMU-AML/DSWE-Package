#' @title KNN : Fit
#' @description The function models the powercurve using KNN, against supplied arguments
#' @param data a dataframe or a matrix, to be used in modelling
#' @param xCol a vector or numeric values stating the column number of features
#' @param yCol a numerical or a vector value stating the column number of target
#' @param subsetSelection a boolean, default value is FALSE, if TRUE returns the best beature combination
#'
#' @return a list containing :
#'         1) data
#'         2) column number of features
#'         3) column number of target
#'         4) best K
#' @export
#' @import FNN
KnnFit = function(data, xCol, yCol, subsetSelection = FALSE){

  if(!is.matrix(data) || !is.data.frame(data)){

    stop('data provided should either be a matrix or data frame')
  }

  if(!is.numeric(xCol)  || !is.vector(xCol)){

    stop('column number of features should be provided as a numeric or vector')

  }else{

    if(sum(xCol %in% 1:ncol(data)) < length(xCol)){

      stop('column numbers provided are not in the data')
    }
  }

  if(!is.numeric(yCol)  || !is.vector(yCol)){

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

    result = computeBestK(normalizedData[, xCol], normalizedData[, yCol], rangeK)
    returnList = list(bestK = result$bestK, bestRMSE = result$bestRMSE, data = data, xCol = xCol, yCol = yCol)

  }else{

    result = computeBestSubset(normalizedData, xCol, yCol, rangeK)
    returnList = list(bestK = result$bestK, bestRMSE = result$bestRMSE, data = data, xCol = result$bestSubset, yCol = yCol )
  }

  return(returnList)
}
