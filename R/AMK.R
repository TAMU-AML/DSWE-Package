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

#' @title Additive Multiplicative Kernel Regression
#' @description An additive multiplicative kernel regression based on Lee et al. (2015).
#' @param trainX a matrix or dataframe of input variable values in the training dataset.
#' @param trainY a numeric vector for response values in the training dataset.
#' @param testX a matrix or dataframe of test input variable values to compute predictions.
#' @param bw a numeric vector or a character input for bandwidth. If character, bandwidth computed internally; the input should be either \code{'dpi'} or \code{'dpi_gap'}. Default is \code{'dpi_gap'}. See \code{details} for more information.
#' @param nMultiCov an integer or a character input specifying the number of multiplicative covariates in each additive term. Default is 3 (same as Lee et al., 2015). The character inputs can be: \code{'all'} for a completely multiplicative model, or \code{'none'} for a completely additive model. Ignored if the number of covariates is 1.
#' @param fixedCov an integer vector specifying the fixed covariates column number(s), default value is \code{c(1,2)}. Ignored if \code{nMultiCov} is set to \code{'all'} or \code{'none'} or if the number of covariates is less than 3.
#' @param cirCov an integer vector specifying the circular covariates column number(s) in \code{trainX}, default value is \code{NA}.
#'
#' @return a numeric vector for predictions at the data points in \code{testX}.
#' @details This function is based on Lee et al. (2015). Main features  are: 
#' \itemize{
#' \item Flexible number of multiplicative covariates in each additive term, which can be set using \code{nMultiCov}.
#' \item Flexible number and columns for fixed covariates, which can be set using \code{fixedCov}. The default option \code{c(1,2)} sets the first two columns as fixed covariates in each additive term.
#' \item Handling the data with gaps when the direct plug-in estimator used in Lee et al. fails to return a finite bandwidth. This is set using the option \code{bw = 'dpi_gap'} for bandwidth estimation.  
#' }
#' @importFrom KernSmooth dpill
#' @importFrom mixtools normalmixEM
#' @importFrom stats var
#' @examples 
#' data = data1
#' trainX = as.matrix(data[c(1:100),2])
#' trainY = data[c(1:100),7]
#' testX = as.matrix(data[c(101:110),2])
#' AMK_prediction = AMK(trainX, trainY, testX, bw = 'dpi_gap', cirCov = NA)
#' 
#' @references  Lee, Ding, Genton, and Xie, 2015, “Power curve estimation with multivariate environmental factors for inland and offshore wind farms,” Journal of the American Statistical Association, Vol. 110, pp. 56-67, DOI:10.1080/01621459.2014.977385. 
#' @export
AMK = function(trainX, trainY, testX, bw = 'dpi_gap', nMultiCov = 3, fixedCov = c(1, 2), cirCov = NA ){
  
  if (!is.matrix(trainX) && !is.data.frame(trainX)) {
    stop("trainX must be a matrix or a dataframe.")
  }
  
  nCov = ncol(trainX)
  
  if (!is.matrix(testX) && !is.data.frame(testX)) {
    
    stop("testX must be a matrix or a dataframe.")
    
  }else if(ncol(testX) != ncol(trainX)){
    
    stop("number of columns in testX and trainX must be the same")
  }
  
  
  
  if (!is.numeric(trainY)){
    stop("trainY must be numeric/vector.")
  }
  
  if (length(trainY) != nrow(trainX)){
    stop("number of datapoints in trainX and trainY must be the same.")
  }
  
  if (!is.numeric(bw)) {
    if (bw != "dpi" && bw != "dpi_gap"){
      {
        stop("bw must be numeric or set to 'dpi' or 'dpi_gap'.")
      }
    }
  }else if (length(bw)!= nCov){
    stop("length of bw must be same as the number of covariates.")
  }
  
  if (!is.numeric(nMultiCov) ){
    if (nMultiCov != "all" && nMultiCov != "none"){
      stop("nMultiCov must be set to 'all' or 'none' or an integer")
    }
  }
  
  if (nCov == 1){
    nMultiCov = 'all'
  } else if (nCov == 2){
    if (nMultiCov != 'none'){
      nMultiCov = 'all'
    }
  }
  
  if (nMultiCov != "all" && nMultiCov != "none"){
    
    if (!is.numeric(nMultiCov) || nMultiCov%%1 != 0){
      
      stop("if nMultiCov is not set to 'all' or 'none', then it must be set to an integer greater than 1, and less than or equal to the number of covariates.")
      
    }else if(nMultiCov > nCov){
      
      stop('The value of nMultiCov cannot be greater than number of columns in trainX')
      
    }else if(nMultiCov == nCov){
      
      nMultiCov = 'all'
      fixedCov = NULL
      
    }else if(nMultiCov < nCov){
      
      if(!is.null(fixedCov)){
        
        if(!is.numeric(fixedCov)){
          
          stop('fixedCov should either be an integer vector or NULL')
        }else if (sum(fixedCov %in% 1:nCov) != length(fixedCov)){
          
          stop('Any or all the values in fixedCov exceeds the numbr of columns in trainX')
          
        }else if(length(fixedCov) >= nMultiCov){
          
          stop('fixedCov should be less than nMulticov')
          
        }
      }    
    }
  }else if(nMultiCov == 'all' || nMultiCov == 'none'){
    
    fixedCov = NULL
  }
  
  if(!is.na(cirCov) && !is.null(cirCov)){
    
    if(!is.numeric(cirCov) && !is.vector(cirCov)){
      
      stop('cirCov should be provided as NA or numeric/vector')
      
    }else if (sum(cirCov %in% 1:nCov) != length(cirCov)){
      
      stop('Any or all the values in cirCov exceeds the numbr of columns in trainX')
      
    }
    
  }

  returnObj = kernpred(trainX, trainY, testX, bw, nMultiCov, fixedCov, cirCov)
  return(returnObj)
}
