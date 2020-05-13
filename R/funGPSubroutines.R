
estimateParameters= function(datalist, covCols, yCol){
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
  optimResult = stats::optim(par = parInit, fn = objFun, gr = objGrad, method = 'BFGS', control = list(maxit = 1000, trace = 1, REPORT = 1)) # , lower = c(rep(0.001,nCov+2),-Inf))
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

