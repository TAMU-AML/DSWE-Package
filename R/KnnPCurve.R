#' @title KNN : PowerCurve
#' @description The function models the powercurve using KNN, against supplied arguments.
#' @param dname The data set to be used for modelling, could be a dataframe or matrix.
#' @param xvar a vector stating column number of non circular features.
#' @param xvar.circ a vector stating column number of circular features.
#' @param yvar a numerical value stating column number of target
#' @param kfold split the data set should go through, such as kfold = 10
#' @param forward_stepwise variable selection technique, which gives the best result for 1 feature, 2 features etc
#'
#' @return a dataframe consisting of best features and their rmse for forward stepwise, else a numeric value
#' @export
#' @import FNN
KnnPCurve = function(dname, xvar, xvar.circ = NULL,  yvar, kfold, forward_stepwise = FALSE){

  xcol = c(xvar, xvar.circ)

  if(length(xvar.circ) > 0){

    dname[, xvar.circ] = sin(dname[, xvar.circ])
  }

  for (feature in xcol) {

    dname[, feature] = (dname[, feature] - min(dname[, feature])) / (max(dname[, feature]) - min(dname[, feature]))

  }
  if(forward_stepwise == FALSE){

    result = KNN.Internal1(dname, xcol, yvar, kfold)
  }else{

    result = KNN.Internal2(dname, xcol, yvar, kfold)
  }

  return(result)
}
