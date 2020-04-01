###
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
    theta[i] = mean(unlist(lapply(datalist, function(x) sd(x[,covCols[i]]))))
  }
  beta = mean(unlist(lapply(datalist, function(x) mean(x[,yCol]))))
  sigma_f = mean(unlist(lapply(datalist, function(x) sd(x[,yCol])/sqrt(2))))
  sigma_n = mean(unlist(lapply(datalist, function(x) sd(x[,yCol])/sqrt(2))))
  parInit = c(theta,sigma_f,sigma_n,beta)
  objFun = function(par){computeloglikSum(datalist,covCols,yCol,
                                          params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  objGrad = function(par){computeloglikGradSum(datalist,covCols,yCol,
                                               params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  optimResult = optim(par = parInit, fn = objFun, gr = objGrad, method = 'BFGS', control = list(maxit = 1000, trace = 1, REPORT = 1)) # , lower = c(rep(0.001,nCov+2),-Inf))
  estimatedParams = list(theta = abs(optimResult$par[1:nCov]), sigma_f = abs(optimResult$par[nCov+1]), sigma_n = abs(optimResult$par[nCov+2]), beta = optimResult$par[nCov+3])
  objVal = optimResult$value
  return(list(estimatedParams = estimatedParams,objVal = objVal))
}

###
computeDiffCov = function(datalist, covCols, yCol, params, testset){
  theta = params$theta
  sigma_f = params$sigma_f
  sigma_n = params$sigma_n
  beta = params$beta
  X1 = as.matrix(datalist[[1]][,covCols])
  y1 = as.numeric(datalist[[1]][,yCol])
  correlMat1 = computeCorrelMat(X1,X1,theta)
  KX1X1 = (sigma_f^2)*correlMat1
  diag(KX1X1) = diag(KX1X1) + (sigma_n^2)
  upperCholKX1X1 = chol(KX1X1)
  rm(KX1X1)
  gc()
  KXTX1 = (sigma_f^2)*computeCorrelMat(testset,X1,theta)
  mu1 = beta + (KXTX1%*%backsolve(upperCholKX1X1,forwardsolve(t(upperCholKX1X1),y1-beta)))


  X2 = as.matrix(datalist[[2]][,covCols])
  y2 = as.numeric(datalist[[2]][,yCol])
  correlMat2 = computeCorrelMat(X2,X2,theta)
  KX2X2 = (sigma_f^2)*correlMat2
  diag(KX2X2) = diag(KX2X2) + (sigma_n^2)
  upperCholKX2X2 = chol(KX2X2)
  rm(KX2X2)
  gc()
  KXTX2 = (sigma_f^2)*computeCorrelMat(testset,X2,theta)
  mu2 = beta + (KXTX2%*%backsolve(upperCholKX2X2,forwardsolve(t(upperCholKX2X2),y2-beta)))


  KX2X1 = (sigma_f^2)*computeCorrelMat(X2,X1,theta)


  ls1 = backsolve(upperCholKX2X2,forwardsolve(t(upperCholKX2X2),t(KXTX2)))
  ls24 = backsolve(upperCholKX1X1,forwardsolve(t(upperCholKX1X1),t(KXTX1)))
  ls3 = backsolve(upperCholKX2X2,forwardsolve(t(upperCholKX2X2),KX2X1))

  K2 = KXTX2%*%ls1
  K1 =KXTX1%*%ls24
  K21 =KXTX2%*%ls3%*%ls24

  diffCovMat =  K2 + K1 - (2*K21)
  diffCovMat = (diffCovMat + t(diffCovMat))/2

  returnList = list(diffCovMat = diffCovMat, mu2 = mu2, mu1 = mu1)
  return(returnList)
}

###
computeConfBand = function(diffCovMat, confLevel){
  eigDecomp = eigen(diffCovMat, symmetric = T)
  eigvals = eigDecomp$values
  lastEigIdx = min(which(eigvals < 1e-08 ))
  lambda = diag(eigvals[1:lastEigIdx])
  eigVec = eigDecomp$vectors[,1:lastEigIdx]
  radius = sqrt(qchisq(confLevel,lastEigIdx))
  nSamples = 1000
  Z = matrix(0,nrow = nSamples, ncol = lastEigIdx)
  n = 1
  set.seed(1)
  while (n <= nSamples){
    zSample = rnorm(lastEigIdx)
    if (sqrt(sum(zSample^2))<=radius){
      Z[n,] = zSample
      n = n+1
    }
  }
  G = eigVec%*%sqrt(lambda)%*%t(Z)
  band = cbind(apply(G,1,max),apply(G,1,min))
  band = abs(band)
  band = apply(band,1,max)
  return(band)
}

###
computeloglikSum = function(datalist,covCols,yCol,params){
  logliksum = sum(unlist(lapply(datalist, function(X) computeloglik(as.matrix(X[,covCols]),as.numeric(X[,yCol]),params))))
  return(logliksum)
}


###
computeloglik = function(x,y,params){
  theta = params$theta
  sigma_f = params$sigma_f
  sigma_n = params$sigma_n
  beta = params$beta
  correlMat = computeCorrelMat(x,x,theta)
  covMat = (sigma_f^2)*correlMat
  diag(covMat) = diag(covMat) + (sigma_n^2)
  upperCholMat = chol(covMat)
  rm(covMat)
  gc()
  loglik = ((1/2)*t(y-beta)%*%backsolve(upperCholMat,forwardsolve(t(upperCholMat),y-beta))) + (sum(log(abs(diag(upperCholMat))))) + (nrow(upperCholMat)*log(2*pi)/2)
  rm(upperCholMat)
  gc()
  return(as.numeric(loglik))
}

###
computeloglikGradSum = function(datalist,covCols,yCol,params){

  loglikGradSum = Reduce("+",lapply(datalist, function(X) computeloglikGradient(as.matrix(X[,covCols]),as.numeric(X[,yCol]),params)))
  return(loglikGradSum)
}

###
computeloglikGradient = function(x,y,params){
  theta = params$theta
  sigma_f = params$sigma_f
  sigma_n = params$sigma_n
  beta = params$beta
  correlMat = computeCorrelMat(x,x,theta)
  covMat = (sigma_f^2)*correlMat
  diag(covMat) = diag(covMat) + (sigma_n^2)
  upperCholMat = chol(covMat)

  gradVal = rep(0,length(theta)+3)

  alpha = backsolve(upperCholMat,forwardsolve(t(upperCholMat),y-beta))

  for (i in 1:length(theta)){
    delThetaMat = ((sigma_f^2)*(outer(x[,i],x[,i],"-")^2)/(theta[i]^3))*correlMat
    gradVal[i] = -0.5*sum(diag((alpha%*%t(alpha)%*%delThetaMat)-(backsolve(upperCholMat,forwardsolve(t(upperCholMat),delThetaMat)))))
  }
  rm(delThetaMat)
  gc()
  delSigma_fMat = 2*sigma_f*correlMat
  gradVal[length(theta)+1] = -0.5*sum(diag((alpha%*%t(alpha)%*%delSigma_fMat)-(backsolve(upperCholMat,forwardsolve(t(upperCholMat),delSigma_fMat)))))
  rm(delSigma_fMat)
  gc()
  gradVal[length(theta)+2] = -0.5*sum(diag((alpha%*%t(alpha)%*%(2*sigma_n*diag(nrow(upperCholMat))))-
                                             (backsolve(upperCholMat,forwardsolve(t(upperCholMat),(2*sigma_n*diag(nrow(upperCholMat))))))))
  oneVec = rep(1,nrow(upperCholMat))
  solOneVec =  backsolve(upperCholMat,forwardsolve(t(upperCholMat),oneVec))
  rm(upperCholMat)
  gc()
  gradVal[length(theta)+3] = 0.5*((2*beta*(t(oneVec)%*%solOneVec))-(t(y)%*%solOneVec)-(t(oneVec)%*%(alpha+(beta*solOneVec))))
  return(gradVal)

}

###
computeCorrelMat = function(x1,x2,theta){

    correlMat = matrix(0,nrow = nrow(x1),ncol = nrow(x2))
    for (i in 1:length(theta)){
      correlMat = correlMat + ((outer(x1[,i],x2[,i],"-")/theta[i])^2)
    }
    correlMat = exp(-0.5*correlMat)

  return(correlMat)
}

