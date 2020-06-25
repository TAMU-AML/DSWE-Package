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
#'
#' @param trainX a matrix or dataframe to be used in modelling
#' @param trainY a numeric or vector as a target
#' @param testX a matrix or dataframe, to be used in computing the predictions
#' @param bw a vector or character input. If character, the input should be 'dpi' or 'dpi_gap'
#' @param nMultiCov a numerical value specifying the number of covariates in multiplicative term
#' @param fixedCov a vector or numeric specifying the fixed covariates column number, default value is c(1,2)
#' @param cirCov a vector or numeric specifying the circular covariates column number, default value is NA
#'
#' @return a vector or numeric predictions on user provided test data
#'
#' @importFrom KernSmooth dpill
#' @importFrom mixtools normalmixEM
#' @export
AMK = function(trainX, trainY, testX, bw = 'dpi_gap', nMultiCov = 3, fixedCov = c(1, 2), cirCov = NA ){
  
  nCov = ncol(trainX)
  
  if (!is.matrix(trainX) && !is.data.frame(trainX)) {
    stop("trainX must be a matrix or a dataframe.")
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
  
  if (nMultiCov != "all" && nMultiCov != "none"){
    
    if (!is.numeric(nMultiCov) || nMultiCov%%1 != 0){
      
      stop("if nMultiCov is not set to 'all' or 'none', then it must be set to an integer greater than 1, and less than or equal to the number of covariates.")
      
    }else if(nMultiCov > nCov){
      
      stop('The value of nMultiCov cannot be greater than number of columns in trainX')
      
    }else if(nMultiCov == nCov){
      
      nMultiCov = 'all'
      fixedCov = NA
      
    }else if(nMultiCov < nCov){
      
      if(!is.null(fixedCov) && !is.na(fixedCov)){
        
        if(!is.numeric(fixedCov) && !is.vector(fixedCov)){
          
          stop('fixedCov should be provided as Null, NA or numeric/vector')
        }else if (sum(fixedCov %in% 1:nCov) != length(fixedCov)){
          
          stop('Any or all the values in fixedCov exceeds the numbr of columns in trainX')
          
        }else if(length(fixedCov) > nMultiCov){
          
          stop('fixedCov should be less than or equal to nMulticov')
          
        }
      }    
    }
  }else if(nMultiCov == 'all' || nMultiCov == 'none'){
    
    fixedCov = NA
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
