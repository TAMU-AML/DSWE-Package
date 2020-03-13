#' @title Additive Multiplicative Kernel Regression
#'
#' @param trainX a matrix or dataframe to be used in modelling
#' @param trainY a numeric or vector as a target
#' @param testX a matrix or dataframe, to be used in computing the predictions
#' @param bw a vector of character input. If character, the input should be 'dpi' or 'dpi-gap'
#' @param nMultiCov a numerical value specifying the number of covariates in multiplicative term
#' @param fixedCov a vector or numeric specifying the fixed covariates column number
#' @param cirCov a vector or numeric specifying the circular covariates column number
#'
#' @import KernSmooth
#' @export
AMK = function(trainX, trainY, testX, bw = 'dpi', nMultiCov = ncol(trainX), fixedCov = NA, cirCov = NA ){
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
  if (nCov == 2) {
    if (nMultiCov != "none"){
      nMultiCov = "all"
      message("Setting nMultiCov to 'all', since there are only 2 covariates. It can be set to 'none' for additive kernels.")
    }
  }
  if (nCov == 1) {
    nMultiCov = "all"
  }
  if (nMultiCov != "all" && nMultiCov != "none"){
    if (!is.numeric(nMultiCov) || nMultiCov%%1 != 0){
      stop("if nMultiCov is not set to 'all' or 'none', then it must be set to an integer greater than 1, and less than the number of covariates.")
    }
    if (nMultiCov >= ncol(trainX) || nMultiCov < 2) {
      stop("if nMultiCov is not set to 'all' or 'none', then it must be set to an integer greater than 1, and less than the number of covariates.")
    }
    if (!is.numeric(fixedCov) || any(fixedCov%%1 != 0)){
      stop("fixedCov must be an integer or a vector of integers.")
    }
    if (length(fixedCov) > nMultiCov){
      stop("length of fixedCov must be less than nMultiCov.")
    }
  }

  returnObj = kernpred(trainX, trainY, testX, bw, nMultiCov, fixedCov, cirCov)
  return(returnObj)
}
