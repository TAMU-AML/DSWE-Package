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

computeThinningNumber = function(trainX){
  thinningNumber = max(apply(trainX,2,function(col) 
    min(which(c(1,abs(stats::pacf(col, plot = FALSE)$acf[,1,1])) <= (2/sqrt(nrow(trainX)))))))
  return(thinningNumber)
}

createThinnedBins = function(dataX, dataY, thinningNumber){
  nData = nrow(dataX)
  thinnedBins = list()
  if (thinningNumber < 2) {
    thinnedBins[[1]] = list(X = dataX, y = dataY)
  } else {
    for (i in 1:thinningNumber){
      nPoints = floor((nData - i)/thinningNumber)
      lastIdx = i + (nPoints*thinningNumber)
      idx = seq(i, lastIdx, length.out = (nPoints+1))
      thinnedBins[[i]] = list(X = dataX[idx,,drop = FALSE], y = dataY[idx])
    }
  }
  
  return(thinnedBins)
}

estimateBinnedParams= function(databins, fast_computation, optim_control){
  nCov = ncol(databins[[1]]$X)
  theta = rep(0,nCov)
  for (i in 1:nCov){
    theta[i] = mean(unlist(lapply(databins, function(bin) stats::sd(bin$X[,i]))))
  }
  beta = mean(unlist(lapply(databins, function(bin) mean(bin$y))))
  sigma_f = mean(unlist(lapply(databins, function(bin) stats::sd(bin$y)/sqrt(2))))
  sigma_n = mean(unlist(lapply(databins, function(bin) stats::sd(bin$y)/sqrt(2))))
  parInit = c(theta,sigma_f,sigma_n,beta)
  objFun = function(par){computeloglikSumTempGP(databins,
                                                params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  objGrad = function(par){computeloglikGradSumTempGP(databins,
                                                     params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  if (fast_computation){
    optimResult = adam_optimizer(databins, parInit, 
                                 optim_control$batch_size, 
                                 optim_control$learn_rate,
                                 optim_control$max_iter,
                                 optim_control$tol, 
                                 optim_control$beta1, 
                                 optim_control$beta2, 
                                 optim_control$epsilon, 
                                 optim_control$logfile )
    estimatedParams = list(theta = abs(optimResult$par[1:nCov]), sigma_f = abs(optimResult$par[nCov+1]), sigma_n = abs(optimResult$par[nCov+2]), beta = optimResult$par[nCov+3])
    objVal = NULL
    gradVal = NULL
    return(list(estimatedParams = estimatedParams,objVal = objVal, gradVal = gradVal))
  } else {
    optimResult = stats::nlminb(start = parInit, objective = objFun, gradient = objGrad, control = list(eval.max = 1000,iter.max = 1000, trace = 0,step.min = 1e-6))
    estimatedParams = list(theta = abs(optimResult$par[1:nCov]), sigma_f = abs(optimResult$par[nCov+1]), sigma_n = abs(optimResult$par[nCov+2]), beta = optimResult$par[nCov+3])
    objVal = -optimResult$objective
    gradVal = -objGrad(optimResult$par)  
    return(list(estimatedParams = estimatedParams,objVal = objVal, gradVal = gradVal))
  }
}

###
computeloglikSumTempGP = function(databins,params){
  logliksum = sum(unlist(lapply(databins, function(bin) computeLogLikGP_(bin$X,bin$y,params))))
  return(logliksum)
}


###
computeloglikGradSumTempGP = function(databins,params){
  
  loglikGradSum = Reduce("+",lapply(databins, function(bin) computeLogLikGradGP_(bin$X,bin$y,params)))
  return(loglikGradSum)
}

estimateLocalFunctionParams= function(trainT, residual){
  theta = sd(trainT)
  sigma_f = sd(residual)/sqrt(2)
  sigma_n = sigma_f
  parInit = c(theta,sigma_f,sigma_n)
  objFun = function(par){computeLogLikGP_(as.matrix(trainT), residual,
                                          params = list(theta=par[1],sigma_f=par[2],sigma_n=par[3],beta=0))}
  objGrad = function(par){computeLogLikGradGPZeroMean_(as.matrix(trainT), residual,
                                                       params = list(theta=par[1],sigma_f=par[2],sigma_n=par[3],beta=0))}
  optimResult = stats::nlminb(start = parInit, objective = objFun, gradient = objGrad )
  estimatedParams = list(theta = abs(optimResult$par[1]), sigma_f = abs(optimResult$par[2]), sigma_n = abs(optimResult$par[3]), beta = 0)
  objVal = optimResult$objective
  gradVal = objGrad(optimResult$par)
  return(list(estimatedParams = estimatedParams,objVal = objVal, gradVal = gradVal))
}


computeLocalFunction = function(residual, traindataT, testdataT, neighbourhood){
  pred = rep(0,length(testdataT))
  for (i in 1:length(testdataT)){
    distance = abs(testdataT[i] - traindataT)
    trainIdx = which(distance < neighbourhood)
    if (length(trainIdx)>0){
      if (var(residual[trainIdx]) < .Machine$double.eps || is.na(var(residual[trainIdx]))){
        warning("While computing g(t), variance of the training residuals is numerically zero for time index: ",testdataT[i], '\nUsing mean of the response as the prediction.\n')
        pred[i] = mean(residual[trainIdx])
      } else {
        params = try (
          expr = estimateLocalFunctionParams(traindataT[trainIdx],residual[trainIdx]), silent = TRUE
          )
        if (class(params) == "try-error"){
          warning("Computing g(t) is numerically unstable for time index: ",testdataT[i], '\nUsing mean of the response as the prediction.\n')
          pred[i] =  mean(residual[trainIdx])
        } else {
          weightedRes = computeWeightedY(as.matrix(traindataT[trainIdx]),residual[trainIdx],params$estimatedParams)
          pred[i] = predictGP(as.matrix(traindataT[trainIdx]),weightedRes,as.matrix(testdataT[i]),params$estimatedParams)    
        }
      }
      
    } else {
      pred[i] = 0
    }
    
  }
  return(pred)
}

adam_optimizer = function(data_bins, par_init, batch_size, learn_rate, max_iter, tol, beta1, beta2, epsilon, logfile){
  t = 0
  nCov = ncol(data_bins[[1]]$X)
  params_t = par_init
  m_t = rep(0, length(par_init))
  v_t = rep(0, length(par_init))
  
  params_mat = matrix(nrow = max_iter, ncol = length(params_t))  
  
  
  while (TRUE){
    t = t+1
    sampled_bin = sample(length(data_bins), 1)
    if (batch_size < length(data_bins[[sampled_bin]]$y)){
      sampled_idx = sample(length(data_bins[[sampled_bin]]$y), batch_size)
    } else {
      sampled_idx = 1:length(data_bins[[sampled_bin]]$y)
    }
    sample_x = data_bins[[sampled_bin]]$X[sampled_idx, , drop=FALSE]
    sample_y = data_bins[[sampled_bin]]$y[sampled_idx]
    
    grad_t = computeLogLikGradGP_(sample_x, sample_y, 
                                  params = list(theta=params_t[1:nCov],sigma_f=params_t[nCov+1],sigma_n=params_t[nCov+2],beta=params_t[nCov+3]))
    m_t = (beta1 * m_t) + ((1 - beta1) * grad_t)
    v_t = (beta2 * v_t) + ((1 - beta2) * (grad_t**2))
    
    m_hat = m_t/(1 - (beta1**t))
    v_hat = v_t/(1 - (beta2**t))
    
    params_prev = params_t
    params_t = params_t - (learn_rate * (m_hat / (sqrt(v_hat) + epsilon)))
    params_mat[t, ] = params_t
    
    if (max(abs(params_t - params_prev)) < tol || t == max_iter){
      if (!is.null(logfile)){
        utils::write.table(params_mat, file = logfile, sep = ",", row.names = FALSE, col.names = FALSE)
      }
      return(list(par=params_t))
    }
  }
}

