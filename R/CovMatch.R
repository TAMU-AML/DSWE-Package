#' @title Covariate Matching
#' @description The function aims to take list of two data sets and returns the after
#' matched data sets using user specified covariates and threshold
#'
#' @param data a list, consisting of data sets to match, also each of the individual data set can be dataframe or a matrix
#' @param wgt a numerical or a vector of threshold values for each covariates, against which matching happens
#' It should be a single value or a vector of values representing threshold for each of the covariate
#' @param xCol a vector stating the column position of covariates used.
#' @param xCol.circ a vector stating the column position of circular variables
#' @usage covMatch(dname, xCol, xCol.circ, wgt, priority)
#' @return a list containing : original data, matched data, MinMax values for covariates in original data and MinMax
#' values for covariates in matched data
#' @export
#' @import foreach

covMatch = function(data, xCol = NULL, xCol.circ = NULL, wgt = 0.2, priority = FALSE){

  # Checks whether the provided data set is a list or not
  if(!is.list(data)){

    stop('Data set provided should be a list containing data sets')

  }

  if(length(data) != 2){


    stop('The number of data sets to match should be equal to two')

  }

  # Checks for non circular covariates
  if(length(xCol) > 0){

    if(!is.vector(xCol)){

      stop('Non circular covariates column number should be provided as a vector')

    }
  }

  # Checks for circular covariates
  if(length(xCol.circ) > 0){

    if(!is.vector(xCol.circ)){

      stop('Circular covariates column number should be provided as a vector')

    }
  }

  # Checks whether any covariate is provided by user or not
  if(!(length(xCol) > 0) && !(length(xCol.circ) > 0)){

    stop('Atleast a single covariate, either circular or non circular should be provided')

  }

  # Checks for dimension compatibility of weight supplied
  if(length(wgt) > 1){

    if(!(length(wgt) == length(c(xCol, xCol.circ)))){

      stop('The weight provided should be a single value or vector with weight for each covariate')

    }
  }

  if(priority == TRUE){

    covDiff = as.numeric(abs(colMeans(data[[1]][, xCol]) - colMeans(data[[2]][, xCol])))

    xCol = xCol[order(covDiff, decreasing = TRUE)]

  }

  ## data set 2 as a baseline
  dname1_ = list(data[[1]], data[[2]])

  ## test set as a baseline
  dname2_ = list(data[[2]], data[[1]])
  filelist_ = list(dname1_, dname2_)

  # sequential computation
  `%do%` = foreach::`%do%`
  matcheddata_ = rep(list(), 2)
  foreach::foreach(i = 1:2) %do% {

    matcheddata_[[i]] = covMatch.Mult(filelist_[[i]], xCol, wgt, xCol.circ)

  }

  ############# Retrieving datasets from 1st matching #########################
  # creating list of matched data set from step 1
  match1_ = matcheddata_[[1]]
  matched1_ = list(match1_[[2]], match1_[[1]])


  ############# Retrieving datasets from 2nd matching #########################
  # creating list of matched data set from step 1
  match2_ = matcheddata_[[2]]
  matched2_ = list(match2_[[1]], match2_[[2]])

  ############ Combining results to generate final matched pairs ###################
  result_ = rep(list(c()), 2)
  result_[[1]] = unique(rbind(matched1_[[1]], matched2_[[1]]))
  result_[[2]] = unique(rbind(matched1_[[2]], matched2_[[2]]))

  MinMaxOriginal = MinMaxData(data, xCol)
  MinMaxMatched = MinMaxData(result_, xCol)

  return(list(originalData = data,  matchedData = result_, MinMaxOriginal = MinMaxOriginal, MinMaxMatched = MinMaxMatched))
}
