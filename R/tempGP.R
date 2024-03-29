# MIT License
# 
# Copyright (c) 2020-2022 Abhinav Prakash and Yu Ding
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
#' @param max_thinning_number An integer specifying the max lag to compute the thinning number. If the PACF does not become insignificant till \code{max_thinning_number}, then \code{max_thinning_number} is used for thinning.
#' @param fast_computation A Boolean that specifies whether to do exact inference or fast approximation. Default is \code{TRUE}.
#' @param vecchia A Boolean that specifies whether to do exact inference or vecchia approximation. Default is \code{TRUE}.
#' @param limit_memory An integer or \code{NULL}. The integer is used sample training points during prediction to limit the total memory requirement. Setting the value to \code{NULL} would result in no sampling, that is, full training data is used for prediction. Default value is \code{5000}.
#' @param optim_control A list parameters passed to the Adam optimizer when \code{fast_computation} is set to \code{TRUE}. The default values have been tested rigorously and tend to strike a balance between accuracy and speed. \itemize{
#' \item \code{batch_size}: Number of training points sampled at each iteration of Adam.
#' \item \code{learn_rate}: The step size for the Adam optimizer.
#' \item \code{max_iter}: The maximum number of iterations to be performed by Adam.
#' \item \code{tol}: Gradient tolerance.
#' \item \code{beta1}: Decay rate for the first moment of the gradient.
#' \item \code{beta2}: Decay rate for the second moment of the gradient.
#' \item \code{epsilon}: A small number to avoid division by zero.
#' \item \code{logfile}: A string specifying a file name to store hyperparameters value for each iteration.
#' }
#' 
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
#' @importFrom stats pacf sd predict
#' @importFrom utils write.table
#' @examples
#' 
#'     data = DSWE::data1
#'     trainindex = 1:50 #using the first 50 data points to train the model
#'     traindata = data[trainindex,]
#'     xCol = 2 #input variable columns
#'     yCol = 7 #response column
#'     trainX = as.matrix(traindata[,xCol])
#'     trainY = as.numeric(traindata[,yCol])
#'     tempGPObject = tempGP(trainX, trainY)
#' 
#' 
#' @references Prakash, A., Tuo, R., & Ding, Y. (2022). "The temporal overfitting problem with applications in wind power curve modeling." Technometrics. \doi{10.1080/00401706.2022.2069158}.
#' @references Katzfuss, M., & Guinness, J. (2017). A general framework for Vecchia approximations of Gaussian processes. \doi{1708.06302}.
#' @export
#' 

tempGP = function(trainX, trainY, trainT = NULL, 
                  fast_computation = TRUE,
                  limit_memory = 5000L,
                  max_thinning_number = 20L,
                  vecchia=TRUE,
                  optim_control = list(batch_size = 100L, 
                                       learn_rate = 0.05, 
                                       max_iter = 5000L, 
                                       tol = 1e-6,
                                       beta1 = 0.9, 
                                       beta2 = 0.999, 
                                       epsilon = 1e-8,
                                       logfile = NULL
                                       )){
  
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
    
    if (!inherits(max_thinning_number, "integer")){
      stop('max_thinning_number must be an integer.') 
    }
    
    trainX = trainX[order(trainT),,drop = FALSE]
    trainY = trainY[order(trainT)]
    trainT = trainT[order(trainT)]
    
  } else {
    trainT = c(1:length(trainY))
  }
  
  if (!(class(limit_memory) %in% c("integer", "numeric", "NULL"))){
    stop("limit memory must be an integer or NULL")
  } else if (inherits(limit_memory, "numeric")){
    limit_memory = as.integer(limit_memory)
  }
  
  thinningNumber = computeThinningNumber(trainX, max_thinning_number)

  if (thinningNumber > 0){
    thinnedBins = createThinnedBins(trainX,trainY,thinningNumber)  
  } else{
    thinnedBins = list(list(x = trainX, y = trainY))
  }
  
  if (!is.logical(vecchia)) {
    stop('vecchia must be boolean.')
  }
  
  if (vecchia == FALSE){
  optimResult = estimateBinnedParams(thinnedBins, fast_computation, optim_control)
  if (inherits(limit_memory, "integer")){
    ntrain = nrow(trainX)
    if (limit_memory < ntrain) {
      pred_index = sample(ntrain, limit_memory)
      activeX = trainX[pred_index,, drop=FALSE]
      activeY = trainY[pred_index]
    } else {
      activeX = trainX
      activeY = trainY
    }
  } else {
    activeX = trainX
    activeY = trainY
  } 
  weightedY = computeWeightedY(activeX, activeY, optimResult$estimatedParams)
  modelF = list(X = activeX, y = activeY, weightedY = weightedY)
  trainResiduals = trainY - predictGP(modelF$X, modelF$weightedY, trainX, optimResult$estimatedParams)
  modelG = list(residuals = trainResiduals, time_index = trainT)
  
  }else {
    optimResult=fit_scaled_thinned(y=trainY,inputs=trainX,thinnedBins=thinnedBins,T=thinningNumber)
    modelF=optimResult
    trainResiduals = trainY - predictions_scaled_thinned(optimResult,trainX,m=100,joint=TRUE,nsims=0,
                                                         predvar=FALSE,scale='parms')
    modelG = list(residuals = trainResiduals, time_index = trainT)
  }
    
  output = list(trainX = trainX, trainY = trainY, trainT = trainT, 
                thinningNumber = thinningNumber, modelF = modelF, 
                modelG = modelG, estimatedParams = optimResult$estimatedParams, 
                llval = optimResult$objVal, gradval = optimResult$gradVal,vecchia = vecchia)
  class(output) = "tempGP"
  return(output)
}


#' @title predict from temporal Gaussian process
#' @description predict function for tempGP objects. This function computes the prediction \code{f(x)} or \code{f(x) + g(t)} depending on the temporal distance between training and test points and whether the time indices for the test points are provided. 
#' @param object An object of class tempGP.
#' @param testX A matrix with each column corresponding to one input variable.
#' @param testT A vector of time indices of the test points. When \code{NULL}, only function \code{f(x)} is used for prediction. Default is \code{NULL}.
#' @param trainT Optional argument to override the existing trainT indices of the \code{tempGP} object.
#' @param ... additional arguments for future development
#' @return A vector of predictions at the testpoints in \code{testX}.
#' @examples 
#'    data = DSWE::data1
#'    trainindex = 1:50 #using the first 50 data points to train the model
#'    traindata = data[trainindex,]
#'    xCol = 2 #input variable columns
#'    yCol = 7 #response column
#'    trainX = as.matrix(traindata[,xCol])
#'    trainY = as.numeric(traindata[,yCol])
#'    tempGPObject = tempGP(trainX, trainY)
#'    testdata = DSWE::data1[101:110,] # defining test data 
#'    testX = as.matrix(testdata[,xCol, drop = FALSE])
#'    predF = predict(tempGPObject, testX)
#' 
#' @export
#' 
predict.tempGP = function(object, testX, testT = NULL, trainT = NULL,...){
  
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
  
  if (ncol(testX) != ncol(object$trainX)){
    
    stop("Number of columns in testX must be the same as that of the tempGP's trainX.")
  }
  
  if (!is.null(trainT)){
    
    if(!is.numeric(trainT)){
      stop('trainT must be numeric or integer vector.')
    }
    
    if(any(!is.finite(trainT))){
      stop('trainT must have finite values only.')
    }
    
    if(length(trainT) != length(object$trainY)){
      stop('trainT must have the same number of data points as the tempGP object.')
    }
    
    cat("Replacing old trainT with the new one.")
    object$trainT = trainT
    object$modelG$time_index = trainT
    
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
  if (object$vecchia == TRUE){
    predF = predictions_scaled_thinned(object$modelF,locs_pred = testX,m=400,joint=TRUE,nsims=0,
                                       predvar=FALSE,scale='parms')
    
  } else {
  predF = predictGP(object$modelF$X, object$modelF$weightedY, testX, object$estimatedParams)
  }
  if (is.null(testT)){
      
    return(predF)
    
  } else {
    
    predG = computeLocalFunction(object$modelG$residuals, object$modelG$time_index, testT, object$thinningNumber)
    
    return(predF + predG)
    
  }
  
}

#' @title Update the data in a tempGP object
#' @description This function updates \code{trainX}, \code{trainY}, and \code{trainT} in a \code{tempGP} object. By default, if the new data has \code{m} data points, the function removes top \code{m} data points from the tempGP object and appends the new data at the bottom, thus keeping the total number of data points the same. This can be overwritten by setting \code{replace = FALSE} to keep all the data points (old and new). The method also updates \code{modelG} by computing and updating residuals at the new data points. \code{modelF} can be also be updated by setting the argument \code{updateModelF} to \code{TRUE}, though not required generally (see comments in the \code{Arguments}.)
#' @param object An object of class tempGP.
#' @param newX A matrix with each column corresponding to one input variable.
#' @param newY A vector with each element corresponding to the output at the corresponding row of \code{newX}.
#' @param newT A vector with time indices of the new datapoints. If \code{NULL}, the function assigns natural numbers starting with one larger than the existing time indices in \code{trainT}.
#' @param replace A boolean to specify whether to replace the old data with the new one, or to add the new data while still keeping all the old data. Default is TRUE, which replaces the top \code{m} rows from the old data, where \code{m} is the number of data points in the new data.
#' @param updateModelF A boolean to specify whether to update \code{modelF} as well. If the original \code{tempGP} model is trained on a sufficiently large dataset (say one year), updating \code{modelF} regularly may not result in any significant improvement, but can be computationally expensive.
#' @param ... additional arguments for future development
#'
#' @return An updated object of class \code{tempGP}.
#' @examples 
#'    data = DSWE::data1
#'    trainindex = 1:50 #using the first 50 data points to train the model
#'    traindata = data[trainindex,]
#'    xCol = 2 #input variable columns
#'    yCol = 7 #response column
#'    trainX = as.matrix(traindata[,xCol])
#'    trainY = as.numeric(traindata[,yCol])
#'    tempGPObject = tempGP(trainX, trainY)
#'    newdata = DSWE::data1[101:110,] # defining new data  
#'    newX = as.matrix(newdata[,xCol, drop = FALSE])
#'    newY = as.numeric(newdata[,yCol])
#'    tempGPupdated = updateData(tempGPObject, newX, newY)
#' @export
#' 
updateData.tempGP = function(object,newX, newY, newT = NULL, replace = TRUE ,updateModelF = FALSE, ...){
  
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
  
  if (ncol(newX) != ncol(object$trainX)){
    
    stop('Number of columns in newX must be the same as that trainX in object.')
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
    newX = newX[order(newT),,drop = FALSE]
    newY = newY[order(newT)]
    newT = newT[order(newT)]
    
  } else {
    
    newT = c((object$trainT[length(object$trainY)]+1):(object$trainT[length(object$trainY)]+length(newY)))
  }
  
  if (replace){
    
    if (length(newY) < length(object$trainY)){
      
      object$trainX = rbind(object$trainX[-c(1:length(newY)),,drop = FALSE],newX)
      object$trainY = c(object$trainY[-c(1:length(newY))],newY)
      object$trainT = c(object$trainT[-c(1:length(newY))],newT)
      
    } else {
      
      object$trainX = newX
      object$trainY = newY
      object$trainT = newT
      
    }
  } else {
    
    object$trainX = rbind(object$trainX, newX)
    object$trainY = c(object$trainY, newY)
    object$trainT = c(object$trainT, newT)
    
  }
  
  if (updateModelF){
    
    if(!object$vecchia){
      
    object$modelF$X = object$trainX
    object$modelF$y = object$trainY
    weightedY = computeWeightedY(object$modelF$X, object$modelF$y, object$estimatedParams)
    object$modelF$weightedY = weightedY
    residuals = object$trainY - predictGP(object$modelF$X, object$modelF$weightedY, object$trainX, object$estimatedParams)
    object$modelG$residuals = residuals
    
    } else {
      
      object$modelF$X = object$trainX
      object$modelF$y = object$trainY
      residuals = object$trainY - predictions_scaled_thinned(object$modelF,object$modelF$X,m=100,joint=TRUE,nsims=0,predvar=FALSE,scale='parms')
      object$modelG$residuals = residuals
    }
    
  } else {
    if(!object$vecchia){
    
    newResiduals = newY - predictGP(object$modelF$X, object$modelF$weightedY, newX, object$estimatedParams)
    
    }else{
    
        newResiduals = newY - predictions_scaled_thinned(object$modelF,object$modelF$X,m=100,joint=TRUE,nsims=0,predvar=FALSE,scale='parms')
    }
    
    if (replace){
      if (length(newY) < length(object$trainY)){
        object$modelG$residuals = c( object$modelG$residuals[-c(1:length(newY))],newResiduals)
        
      } else {
      
        object$modelG$residuals = newResiduals
      }
    } else {
      
      object$modelG$residuals = c( object$modelG$residuals,newResiduals)
      
    }
    
  
  }
  
  object$modelG$time_index = object$trainT

  return(object)
}

#' @title Updating data in a model
#' @description \code{updateData} is a generic function to update data in a model.
#' @param object A model object
#' @param ... additional arguments for passing to specific methods
#'
#' @return The returned value would depend on the class of its argument \code{object}.
#'
#' @seealso \code{\link{updateData.tempGP}}
#' @export
updateData = function(object, ...){

  UseMethod("updateData")

}

