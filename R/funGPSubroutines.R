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

#' @useDynLib DSWE
#' @importFrom Rcpp sourceCpp
#' 
estimateParameters= function(datalist, covCols, yCol, opt_method){
  maxDataSample = 500
  for (i in 1:length(datalist)){
    if (nrow(datalist[[i]]) > maxDataSample){
      set.seed(1)
      datalist[[i]] = datalist[[i]][sample(nrow(datalist[[i]]), maxDataSample),]
    }
  }

  nCov = length(covCols)
  theta = rep(0,nCov)
  for (i in 1:length(theta)){
    theta[i] = mean(unlist(lapply(datalist, function(x) stats::sd(x[,covCols[i]]))))
  }
  beta = mean(unlist(lapply(datalist, function(x) mean(x[,yCol]))))
  sigma_f = mean(unlist(lapply(datalist, function(x) stats::sd(x[,yCol])/sqrt(2))))
  sigma_n = mean(unlist(lapply(datalist, function(x) stats::sd(x[,yCol])/sqrt(2))))
  parInit = c(theta,sigma_f,sigma_n,beta)
  objFun = function(par){computeloglikSum(datalist,covCols,yCol,
                                          params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  objGrad = function(par){computeloglikGradSum(datalist,covCols,yCol,
                                               params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  optimResult = stats::optim(par = parInit, fn = objFun, gr = objGrad, method = opt_method, control = list(maxit = 1000, trace = 1, REPORT = 1)) # , lower = c(rep(0.001,nCov+2),-Inf))
  estimatedParams = list(theta = abs(optimResult$par[1:nCov]), sigma_f = abs(optimResult$par[nCov+1]), sigma_n = abs(optimResult$par[nCov+2]), beta = optimResult$par[nCov+3])
  objVal = optimResult$value
  return(list(estimatedParams = estimatedParams,objVal = objVal))
}

###
computeDiffCov = function(datalist, covCols, yCol, params, testset, limitMemory){
  theta = params$theta
  sigma_f = params$sigma_f
  sigma_n = params$sigma_n
  beta = params$beta
  if (limitMemory == T){
    maxDataSample = 5000
    for (i in 1:length(datalist)){
      if (nrow(datalist[[i]]) > maxDataSample){
        set.seed(1)
        datalist[[i]] = datalist[[i]][sample(nrow(datalist[[i]]), maxDataSample),]
      }
    }
  }

  X1 = as.matrix(datalist[[1]][,covCols])
  y1 = as.numeric(datalist[[1]][,yCol])

  X2 = as.matrix(datalist[[2]][,covCols])
  y2 = as.numeric(datalist[[2]][,yCol])
  
  testset = as.matrix(testset)
  
  returnList = computeDiffCov_(X1, y1, X2, y2, testset, theta, sigma_f, sigma_n, beta)
  
  return(returnList)
}

###
computeConfBand = function(diffCovMat, confLevel){
  set.seed(1)
  band = computeConfBand_(diffCovMat, confLevel)
  return(as.numeric(band))

}

###
computeloglikSum = function(datalist,covCols,yCol,params){
  logliksum = sum(unlist(lapply(datalist, function(X) computeLogLikGP_(as.matrix(X[,covCols]),as.numeric(X[,yCol]),params))))
  return(logliksum)
}


###
computeloglikGradSum = function(datalist,covCols,yCol,params){

  loglikGradSum = Reduce("+",lapply(datalist, function(X) computeLogLikGradGP_(as.matrix(X[,covCols]),as.numeric(X[,yCol]),params)))
  return(loglikGradSum)
}

