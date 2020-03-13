## Function to calculate weights for a testpoint
ComputeGaussianKernel = function(x,y,lambda){
  x = as.numeric(x)
  y = as.numeric(y)
  lambda = as.numeric(lambda)
  kernel = (1/(sqrt(2*pi)*lambda))*exp(-(x-y)**2/(2*lambda**2))
  return(kernel)
}

ComputeVonMisesKernel = function(D,D0, nu){
  D = as.numeric(D)
  D0 = as.numeric(D0)
  nu = as.numeric(nu)
  kernel =  exp(nu*cos(D-D0))/(2*pi*besselI(x=nu,nu=0))
  return(kernel)
}


## Function to calculate weights for a trivariate kernel given the bivariate kernel and the third input variable
calculateWeights = function(trainX,testpoint,bandwidth,nMultiCov,fixedCov,cirCov){
  nCov = ncol(trainX)
  if (nMultiCov == "all"){
    weights = matrix(NA,nrow(trainX),1)
    kernel = rep(1,nrow(trainX))
    for (i in 1:nCov){
      if (i %in% cirCov){
        kernel = kernel * ComputeVonMisesKernel(trainX[,i],testpoint[i],bandwidth[i])
      } else {
        kernel = kernel * ComputeGaussianKernel(trainX[,i],testpoint[i],bandwidth[i])
      }
    }
    weights[,1] = kernel/sum(kernel)
  } else if (nMultiCov == "none"){
    weights = matrix(NA,nrow(trainX),nCov)
    for (i in 1:nCov){
      if (i %in% cirCov){
        kernel = ComputeVonMisesKernel(trainX[,i],testpoint[i],bandwidth[i])
      } else {
        kernel = ComputeGaussianKernel(trainX[,i],testpoint[i],bandwidth[i])
      }
      weights[,i] = kernel/sum(kernel)
    }
  } else {
    nonFixedCov = setdiff(c(1:nCov),fixedCov)
    covCombination = combn(nonFixedCov, (nMultiCov - length(fixedCov)))
    weights = matrix(NA,nrow(trainX),ncol(covCombination))
    for (i in 1:ncol(covCombination)){
      kernel = rep(1,nrow(trainX))
      for (f in fixedCov){
        if (f %in% cirCov){
          kernel = kernel * ComputeVonMisesKernel(trainX[,f],testpoint[f],bandwidth[f])
        } else {
          kernel = kernel * ComputeGaussianKernel(trainX[,f],testpoint[f],bandwidth[f])
        }
      }
      for (j in covCombination[,i]){
        if (j %in% cirCov){
          kernel = kernel * ComputeVonMisesKernel(trainX[,j],testpoint[j],bandwidth[j])
        } else {
          kernel = kernel * ComputeGaussianKernel(trainX[,j],testpoint[j],bandwidth[j])
        }
      }
      weights[,i] = kernel/sum(kernel)
    }
  }
  return(weights)
}

## Function to get bandwidth of x using dpi
computeBandwidth = function(trainY,trainX,cirCov){
  bandwidth = rep(0,ncol(trainX))
  for (i in 1:ncol(trainX)){
    bandwidth[i] = KernSmooth::dpill(trainX[,i],trainY)
  }
  for (i in cirCov){
    bandwidth[i] = bandwidth[i]*pi/180
    bandwidth[i] = 1/((bandwidth[i])^2)
  }
  return(bandwidth)
}



### Function to get mean estimate
kernpred = function(trainX, trainY, testX, bw, nMultiCov, fixedCov, cirCov){

  trainX = as.matrix(trainX)
  trainY = as.numeric(trainY)
  testX = as.matrix(testX)
  if (class(bw)=="character"){
    if (bw == "dpi"){
      bandwidth = computeBandwidth(trainY,trainX,cirCov)
      if (any(!is.finite(bandwidth))){
        message("Bandwidths not finite for some of the covariates. Bandwidths are:")
        for (i in 1:ncol(trainX)){
          message("Covariate ",i," : ",bandwidth[i])
        }
        stop("Estimation stopped. Try with finite bandwidths.")
      }
      pred = computePred(trainX, trainY, testX, bandwidth, nMultiCov, fixedCov, cirCov )
    }else if (bw == "dpi_gap"){
      band = bw.gap(trainY, trainX, id.dir = cirCov)
      if(is.na(band$bw.adp)){

        pred = computePred(trainX, trainY, testX, band$bw.fix, nMultiCov, fixedCov, cirCov )

      }else{

        prediction = rep(NA, nrow(testX))
        for(i in 1:nrow(testX)){
          bandwidth = find.bw(trainY, trainX, testX[i, , drop = F], band)
          prediction[i] = computePredGap(trainX, trainY, testX[i, , drop = F], bandwidth, nMultiCov, fixedCov, cirCov)
        }
        if (any(!is.finite(prediction))){
          warning("some of the testpoints resulted in non-finite predictions.")
        }

        pred = list(pred = prediction)

      }
    }
  }else {
    bandwidth = bw
    pred = computePred(trainX, trainY, testX, bandwidth, nMultiCov, fixedCov, cirCov )
  }
  return(pred)
}

computePred = function(trainX, trainY, testX, bandwidth, nMultiCov, fixedCov, cirCov ){
  if(!is.na(cirCov)){
    for (i in cirCov) {
      trainX[,i] = trainX[,i]*pi/180
      testX[,i] = testX[,i]*pi/180
    }
  }
  pred = rep(NA, nrow(testX))
  for (i in 1:length(pred)){
    weights = calculateWeights(trainX,testX[i,],bandwidth,nMultiCov,fixedCov,cirCov)
    pred[i] = (rowSums(weights)%*%trainY)/ncol(weights)
  }
  if (any(!is.finite(pred))){
    warning("some of the testpoints resulted in non-finite predictions.")
  }
  return(list(bandwidth = bandwidth, pred = pred))
}

