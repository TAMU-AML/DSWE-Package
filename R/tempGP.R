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
#' @description A Gaussian process based power curve model which explicitly models the temporal aspect of the power curve. The model consists of two parts: \code{f(x)} and \code{g(t)}.
#' @param trainX A matrix with each column corresponding to one input variable. 
#' @param trainY A vector with each element corresponding to the output at the corresponding row of \code{trainX}.
#' @param trainT A vector for time indices of the data points. By default, the function assigns natural numbers starting from 1 as the time indices. 
#'
#' @return An object of class \code{tempGP} with the following attributes:
#' \itemize{
#' \item trainX - same as the input matrix \code{trainX}.
#' \item trainY - same as the input vector \code{trainY}.
#' \item thinningNumber - the thinning number computed by the algorithm.
#' \item modelF - A list containing the details of the model for predicting function \code{f(x)}:
#' \itemize{
#' \item X - The input variable matrix for computing the cross-covariance for predictions, same as \code{trainX} unless the model is updated. See \code{\link{updateData.tempGP}} method for details on updating the model.
#' \item y - The response vector, again same as \code{trainY} unless the model is updated.
#' \item weightedY - The weighted response, that is, the response left multiplied by the inverse of the covariance matrix.
#' } 
#' \item modelG - A list containing the details of the model for predicting function \code{g(t)}:
#' \itemize{
#' \item residuals - The residuals after subtracting function \code{f(x)} from the response. Used to predict \code{g(t)}. See \code{\link{updateData.tempGP}} method for updating the residuals.
#' \item time_index - The time indices of the residuals, same as \code{trainT}.
#' }
#' \item estimatedParams - Estimated hyperparameters for function \code{f(x)}.
#' \item llval - log-likelihood value of the hyperparameter optimization for \code{f(x)}.
#' \item gradval - gradient vector at the optimal log-likelihood value.
#' }
#' 
#' @seealso \code{\link{predict.tempGP}} for computing predictions and \code{\link{updateData.tempGP}} for updating data in a tempGP object.
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
#' @description predict function for tempGP objects. This function computes the prediction \code{f(x)} or \code{f(x) + g(t)} depending on the temporal distance between training and test points and whether the time indices for the test points are provided. 
#' @param tempGPObj An object of class tempGP.
#' @param testX A matrix with each column corresponding to one input variable.
#' @param testT A vector of time indices of the test points. When \code{NULL}, only function \code{f(x)} is used for prediction. Default is \code{NULL}.
#' @param trainT Optional argument to override the existing trainT indices of the \code{tempGP} object.
#'
#' @return A vector of predictions at the testpoints in \code{testX}.
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
#' @description This function updates \code{trainX}, \code{trainY}, and \code{trainT} in a \code{tempGP} object. By default, if the new data has \code{m} data points, the function removes top \code{m} data points from the tempGP object and appends the new data at the bottom, thus keeping the total number of data points the same. This can be overwritten by setting \code{replace = FALSE} to keep all the data points (old and new). The method also updates \code{modelG} by computing and updating residuals at the new data points. \code{modelF} can be also be updated by setting the argument \code{updateModelF} to \code{TRUE}, though not required generally (see comments in the \code{Arguments}.)
#' @param tempGPObj An object of class tempGP.
#' @param newX A matrix with each column corresponding to one input variable.
#' @param newY A vector with each element corresponding to the output at the corresponding row of \code{newX}.
#' @param newT A vector with time indices of the new datapoints. If \code{NULL}, the function assigns natural numbers starting with one larger than the existing time indices in \code{trainT}.
#' @param replace A boolean to specify whether to replace the old data with the new one, or to add the new data while still keeping all the old data. Default is TRUE, which replaces the top \code{m} rows from the old data, where \code{m} is the number of data points in the new data.
#' @param updateModelF A boolean to specify whether to update \code{modelF} as well. If the original \code{tempGP} model is trained on a sufficiently large dataset (say one year), updating \code{modelF} regularly may not result in any significant improvement, but can be computationally expensive.
#'
#' @return An updated object of class \code{tempGP}.
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