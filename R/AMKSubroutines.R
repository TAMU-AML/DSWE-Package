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
    weights = matrix(1,nrow(trainX),1)/nrow(trainX)
    kernel = rep(1,nrow(trainX))
    for (i in 1:nCov){
      if (i %in% cirCov){
        covKernel = ComputeVonMisesKernel(trainX[,i],testpoint[i],bandwidth[i])
        covKernel[is.nan(covKernel)] = 0
        if (sum(covKernel) != 0){
          kernel = kernel * covKernel
        }
        
      } else {
        covKernel = ComputeGaussianKernel(trainX[,i],testpoint[i],bandwidth[i])
        if (sum(covKernel) != 0){
          kernel = kernel * covKernel
        }
      }
    }
    if (sum(kernel) != 0){
      weights[,1] = kernel/sum(kernel)  
    }
  } else if (nMultiCov == "none"){
    weights = matrix(1,nrow(trainX),nCov)/nrow(trainX)
    kernel = rep(1,nrow(trainX))
    for (i in 1:nCov){
      if (i %in% cirCov){
        covKernel = ComputeVonMisesKernel(trainX[,i],testpoint[i],bandwidth[i])
        covKernel[is.nan(covKernel)] = 0
        if (sum(covKernel) != 0){
          kernel = kernel * covKernel
        }
        
      } else {
        covKernel = ComputeGaussianKernel(trainX[,i],testpoint[i],bandwidth[i])
        if (sum(covKernel) != 0){
          kernel = kernel * covKernel
        }
      }
      if (sum(kernel) != 0){
        weights[,i] = kernel/sum(kernel)   
      }
    }
  } else {
    nonFixedCov = setdiff(c(1:nCov),fixedCov)
    covCombination = utils::combn(nonFixedCov, (nMultiCov - length(fixedCov)))
    weights = matrix(1,nrow(trainX),ncol(covCombination))/nrow(trainX)
    for (i in 1:ncol(covCombination)){
      kernel = rep(1,nrow(trainX))
      for (f in fixedCov){
        if (f %in% cirCov){
          covKernel = ComputeVonMisesKernel(trainX[,f],testpoint[f],bandwidth[f])
          covKernel[is.nan(covKernel)] = 0
          if (sum(covKernel) != 0){
            kernel = kernel * covKernel
          }
        } else {
          covKernel = ComputeGaussianKernel(trainX[,f],testpoint[f],bandwidth[f])
          if (sum(covKernel) != 0){
            kernel = kernel * covKernel
          }
        }
      }
      for (j in covCombination[,i]){
        if (j %in% cirCov){
          covKernel = ComputeVonMisesKernel(trainX[,j],testpoint[j],bandwidth[j])
          covKernel[is.nan(covKernel)] = 0
          if (sum(covKernel) != 0){
            kernel = kernel * covKernel
          }
        } else {
          covKernel = ComputeGaussianKernel(trainX[,j],testpoint[j],bandwidth[j])
          if (sum(covKernel) != 0){
            kernel = kernel * covKernel
          }
        }
      }
      if (sum(kernel) != 0){
        weights[,i] = kernel/sum(kernel) 
      }
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
  if(all(!is.na(cirCov))){
  for (i in cirCov){
    bandwidth[i] = bandwidth[i]*pi/180
    bandwidth[i] = 1/((bandwidth[i])^2)
  }
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
      if(any(!is.finite(bandwidth))){
        for(i in which(!is.finite(bandwidth))){
          bandwidth[i] = sqrt(var(trainX[, i]))
        }
      }
      pred = computePred(trainX, trainY, testX, bandwidth, nMultiCov, fixedCov, cirCov )
    }else if (bw == "dpi_gap"){
      band = bw.gap(trainY, trainX, id.dir = cirCov)
      if(all(!is.na(cirCov))){
      for (i in cirCov){
        band$bw.fix[i] = band$bw.fix[i]* pi/180
        band$bw.fix[i] = 1/((band$bw.fix[i])^2)
      }
      }
      if(any(!is.finite(band$bw.fix))){
      for(i in which(!is.finite(band$bw.fix))){
        band$bw.fix[i] = sqrt(var(trainX[, i]))
      }
      }
      if(!is.finite(band$bw.adp)){
        pred = computePred(trainX, trainY, testX, band$bw.fix, nMultiCov, fixedCov, cirCov )
        
      }else{
        prediction = rep(NA, nrow(testX))
        for(i in 1:nrow(testX)){
          bandwidth = find.bw(trainY, trainX, testX[i, , drop = FALSE], band)
          if(all(!is.na(cirCov))){
          for (i in cirCov){
            bandwidth[i] = bandwidth[i]*pi/180
            bandwidth[i] = 1/((bandwidth[i])^2)
          }
          }
          
          if(any(!is.finite(bandwidth))){
          for(i in which(is.na(bandwidth))){
            bandwidth[i] = sqrt(var(trainX[, i]))
          }
          }
          prediction[i] = computePredGap(trainX, trainY, testX[i, , drop = FALSE], bandwidth, nMultiCov, fixedCov, cirCov)
        }
        
        pred = prediction
        
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
  return(pred)
}
