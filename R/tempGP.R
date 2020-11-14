# MIT License
# 
# Copyright (c) 2020 Abhinav Prakash, Rui Tuo, and Yu Ding
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

#' @title temporal Gaussian process
#' @description A Gaussian process based power curve model which explicitly models the temporal aspect of the power curve.
#' @param trainX A matrix with each column corresponding to one input variable. 
#' @param trainY A vector with each element corresponding to the output at the corresponding row of \code{trainX}.
#'
#' @return An object of class \code{tempGP} with the following attributes:
#' \itemize{
#' \item trainX
#' \item trainY
#' \item thinningNumber
#' \item estimatedParams
#' \item llval
#' \item gradval
#' }
#' 
#' @seealso [predict.tempGP()] for computing predictions and [updateData.tempGP()] for updating data in a tempGP object.
#' @importFrom stats pacf
#' @export
#' 

tempGP = function(trainX, trainY, trainT = NULL){
  
  trainX = as.matrix((trainX)) #trying to coerce trainX into a matrix, if not already.
  
  if (!is.matrix(trainX)){
    stop('trainX must be a matrix.')
  }
  
  if (!is.numeric(trainX)){
    stop('trainX must be a numeric matrix.')
  }
  
  if (any(!is.finite(trainX))){
    stop('trainX must have finite value values only.')
  }
  
  if (!is.numeric(trainY)){
    stop('trainY must be numeric.')
  }
  
  if(any(!is.finite(trainY))){
    stop('trainY must have finite values only.')
  }
  
  if(nrow(trainX) != length(trainY)){
    stop('Number of rows in trainX and the length of trainY must match.')
  }
  
  if (!is.null(trainT)){
    if(!is.numeric(trainT)){
      stop('trainT must be numeric or integer vector.')
    }
    
    if(any(!is.finite(trainT))){
      stop('trainT must have finite values only.')
    }
    
    if(length(trainT) != length(trainY)){
      stop('trainX, trainY, and trainT must have the same number of data points.')
    }
    
    trainX = trainX[order(trainT),,drop = F]
    trainY = trainY[order(trainT)]
    trainT = trainT[order(trainT)]
    
  } else {
    trainT = c(1:length(trainY))
  }
  
  #cat("All test passed.\n")
  thinningNumber = computeThinningNumber(trainX)
  cat("thinning number =",thinningNumber,'\n')
  if (thinningNumber > 0){
    thinnedBins = createThinnedBins(trainX,trainY,thinningNumber)  
  } else{
    thinnedBins = list(list(x = trainX, y = trainY))
  }
  cat("Estimating hyperparameters...\n")
  optimResult = estimateBinnedParams(thinnedBins)
  weightedY = computeWeightedY(trainX, trainY, optimResult$estimatedParams)
  modelF = list(X = trainX, y = trainY, weightedY = weightedY)
  trainResiduals = trainY - predictGP(trainX, weightedY, trainX, optimResult$estimatedParams)
  modelG = list(residuals = trainResiduals, time_index = trainT)
  output = list(trainX = trainX, trainY = trainY, trainT = trainT, 
                thinningNumber = thinningNumber, modelF = modelF, 
                modelG = modelG, estimatedParams = optimResult$estimatedParams, 
                llval = -optimResult$objVal, gradval = -optimResult$gradVal)
  class(output) = "tempGP"
  return(output)
}


#' @title predict from temporal Gaussian process
#' @description predict function for tempGP objects
#' @param tempGPObj An object of class tempGP
#' @param testX A matrix with each column corresponding to one input variable.
#'
#' @return a vector of predictions at the testpoints in \code{testX}.
#' @export
#' 
predict.tempGP = function(tempGPObj, testX, testT = NULL, trainT = NULL){
  
  testX = as.matrix((testX)) #trying to coerce testX into a matrix, if not already.
  
  if (!is.matrix(testX)){
    
    stop('testX must be a matrix.')
  }
  
  if (!is.numeric(testX)){
    
    stop('testX must be a numeric matrix.')
  }
  
  if (any(!is.finite(testX))){
    
    stop('testX must have finite value values only.')
  }
  
  if (ncol(testX) != ncol(trainX)){
    
    stop("Number of columns in testX must be the same as that of the tempGP's trainX.")
  }
  
  if (!is.null(trainT)){
    
    if(!is.numeric(trainT)){
      stop('trainT must be numeric or integer vector.')
    }
    
    if(any(!is.finite(trainT))){
      stop('trainT must have finite values only.')
    }
    
    if(length(trainT) != length(trainY)){
      stop('trainT must have the same number of data points as the tempGP object.')
    }
    
    cat("Replacing old trainT with the new one.")
    tempGPObj$trainT = trainT
    tempGPObj$modelG$time_index = trainT
    
  }
  
  if (!is.null(testT)){
    
    if(!is.numeric(testT)){
    
        stop('trainT must be numeric or integer vector.')
    }
    
    if(any(!is.finite(testT))){
      
      stop('trainT must have finite values only.')
    }
    if (nrow(testX) != length(testT)){
      
      stop('Number of rows in testX and length of testT must be the same.')
    }
    
  }
  
  #cat("All test passed.\n")
  
  predF = predictGP(tempGPObj$modelF$X, tempGPObj$modelF$weightedY, testX, tempGPObj$estimatedParams)
  
  if (is.null(testT)){
      
    return(predF)
    
  } else {
    
    predG = computeLocalFunction(tempGPObj$modelG$residuals, tempGPObj$modelG$time_index, testT, tempGPObj$thinningNumber)
    
    return(predF + predG)
    
  }
  
}

#' @title Update the data in a tempGP object
#' @description This function updates the trainX, trainY and trainT in a tempGP object. By default, if the new data has \code{m} data points, the function replaces top \code{m} data points from the tempGP object and appends the new data at the bottom, thus keeping the total number of data points same. This can be overwritten by setting \code{replace = FALSE}.
#' @param tempGPObj An object of class tempGP
#' @param newX A matrix with each column corresponding to one input variable.
#' @param newY A vector with each element corresponding to the output at the corresponding row of \code{newX}.
#' @param replace A boolean to specify whether to replace the old data with the new one, or to append the new data. Default is TRUE, which replaces the old data.
#'
#' @return updated object of class tempGP.
#' @export
#' 
updateData.tempGP = function(tempGPObj,newX, newY, newT = NULL, replace = TRUE, updateModelF = FALSE){
  
  newX = as.matrix((newX)) #trying to coerce newX into a matrix, if not already.
  
  if (!is.matrix(newX)){
    stop('newX must be a matrix.')
  }
  
  if (!is.numeric(newX)){
    stop('newX must be a numeric matrix.')
  }
  
  if (any(!is.finite(newX))){
    stop('newX must have finite value values only.')
  }
  
  if (!is.numeric(newY)){
    stop('newY must be numeric.')
  }
  
  if(any(!is.finite(newY))){
    stop('newY must have finite values only.')
  }
  
  if(nrow(newX) != length(newY)){
    stop('Number of rows in newX and the length of newY must match.')
  }
  
  if (ncol(newX) != ncol(tempGPObj$trainX)){
    
    stop('Number of columns in newX must be the same as that trainX in tempGPObj.')
  }
  
  if (!is.null(newT)){
    if(!is.numeric(newT)){
      stop('newT must be numeric or integer vector.')
    }
    
    if(any(!is.finite(newT))){
      stop('newT must have finite values only.')
    }
    
    if(length(newT) != length(newY)){
      stop('newX, newY, and newT must have the same number of data points.')
    }
    newX = newX[order(newT),,drop=F]
    newY = newY[order(newT)]
    newT = newT[order(newT)]
    
  } else {
    
    newT = c((tempGPObj$trainT[length(tempGPObj$trainY)]+1):(tempGPObj$trainT[length(tempGPObj$trainY)]+length(newY)))
  }
  
  if (replace){
    
    if (length(newY) < length(tempGPObj$trainY)){
      
      tempGPObj$trainX = rbind(tempGPObj$trainX[-c(1:length(newY)),,drop = F],newX)
      tempGPObj$trainY = c(tempGPObj$trainY[-c(1:length(newY))],newY)
      tempGPObj$trainT = c(tempGPObj$trainT[-c(1:length(newY))],newT)
      
    } else {
      
      tempGPObj$trainX = newX
      tempGPObj$trainY = newY
      tempGPObj$trainT = newT
      
    }
  } else {
    
    tempGPObj$trainX = rbind(tempGPObj$trainX, newX)
    tempGPObj$trainY = c(tempGPObj$trainY, newY)
    tempGPObj$trainT = c(tempGPObj$trainT, newT)
    
  }
  
  if (updateModelF){
    
    tempGPObj$modelF$X = tempGPObj$trainX
    tempGPObj$modelF$y = tempGPObj$trainY
    weightedY = computeWeightedY(tempGPObj$modelF$X, tempGPObj$modelF$y, tempGPObj$estimatedParams)
    tempGPObj$modelF$weightedY = weightedY
    residuals = tempGPObj$trainY - predictGP(tempGPObj$modelF$X, tempGPObj$modelF$weightedY, tempGPObj$trainX, tempGPObj$estimatedParams)
    tempGPObj$modelG$residuals = residuals
    
    
  } else {
    
    newResiduals = newY - predictGP(tempGPObj$modelF$X, tempGPObj$modelF$weightedY, newX, tempGPObj$estimatedParams)
    if (replace){
      if (length(newY) < length(tempGPObj$trainY)){
        tempGPObj$modelG$residuals = c( tempGPObj$modelG$residuals[-c(1:length(newY))],newResiduals)
        
      } else {
      
        tempGPObj$modelG$residuals = newResiduals
      }
    } else {
      
      tempGPObj$modelG$residuals = c( tempGPObj$modelG$residuals,newResiduals)
      
    }
    
  }
  
  tempGPObj$modelG$time_index = tempGPObj$trainT

  return(tempGPObj)
}

#' @export
updateData = function(object, ...){
  
  UseMethod("updateData")
  
}
