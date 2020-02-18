#' @title Additive Multiplicative Kernel Regression
#'
#' @param trainX
#' @param trainY
#' @param testX
#' @param bw
#' @param nMultiCov
#' @param fixedCov
#' @param cirCov
#'
#' @export
AMK = function(trainX, trainY, testX, bw = 'dpi', nMultiCov = 3, fixedCov = c(1,2), cirCov = 2 ){
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
    if (length(fixedCov) >= nMultiCov){
      stop("length of fixedCov must be less than nMultiCov.")
    }
  }
  if (!is.na(cirCov) || !is.null(cirCov) ){
    if (!is.numeric(cirCov) || cirCov%%1 != 0){
      stop("cirCov must either be null or NA or an integer vector.")
    }
    if (length(cirCov)> nCov){
      stop("Number of circular covariates must be less than or equal to the total number of covariates.")
    }
  }

  returnObj = kernpred(trainX, trainY, testX, bw, nMultiCov, fixedCov, cirCov)
  return(returnObj)
}
