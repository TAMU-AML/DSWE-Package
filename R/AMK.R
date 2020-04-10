#' @title Additive Multiplicative Kernel Regression
#'
#' @param trainX a matrix or dataframe to be used in modelling
#' @param trainY a numeric or vector as a target
#' @param testX a matrix or dataframe, to be used in computing the predictions
#' @param bw a vector or character input. If character, the input should be 'dpi' or 'dpi_gap'
#' @param nMultiCov a numerical value specifying the number of covariates in multiplicative term
#' @param fixedCov a vector or numeric specifying the fixed covariates column number, default value is NA
#' @param cirCov a vector or numeric specifying the circular covariates column number, default value is NA
#'
#'' @return a vector or numeric predictions on user provided test data
#'}
#' @import KernSmooth
#' @export
AMK = function(trainX, trainY, testX, bw = 'dpi_gap', nMultiCov = 3, fixedCov = c(1, 2), cirCov = NA ){

  if (!is.matrix(trainX) && !is.data.frame(trainX)) {
    stop("trainX must be a matrix or a dataframe.")
  }

  nCov = ncol(trainX)

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
    }
  }

  if(!is.null(fixedCov) && !is.na(fixedCov)){

    if(!is.numeric(fixedCov) && !is.vector(fixedCov)){

      stop('fixedCov should be provided as Null, NA or numeric/vector')
    }else if (sum(fixedCov %in% 1:nCov) != length(fixedCov)){

      stop('Any or all the values in fixedCov exceeds the numbr of columns in trainX')

    }else if(length(fixedCov) > nMultiCov){

      stop('fixedCov should be less than or equal to nMulticov')

    }else if(length(fixedCov) == nMultiCov){

      nMultiCov = 'all'
    }
  }else{

    nMultiCov = 'all'
  }

  returnObj = kernpred(trainX, trainY, testX, bw, nMultiCov, fixedCov, cirCov)
  return(returnObj)
}
