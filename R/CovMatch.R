#' @title Covariate Matching
#' @description The function aims to take list of two data sets and returns the after
#' matched data sets using user specified covariates.
#'
#' @param dname This should always be a list, containing data sets to match.
#' @param wgt Vector of threshold values, against which matching happens.
#' It should be a single value (0.2) or a vector of values c(0.2, 0.3).
#' @param cov Vector stating the column position of non circular variables such as wind speed
#' @param cov.circ Vector stating the column position of circular variables such as wind direction,
#'  nacelle position etc.
#' @usage CovMatch(dname, cov, wgt, cov.circ)
#' @return The function returns a list containing after matched data sets.
#' @export
#' @import foreach

CovMatch = function(dname, cov = NULL, wgt = 0.2, cov.circ = NULL ){

  if(length(dname) != 2){


    stop('The number of data sets to match should be equal to two')

  }

  # Checks whether the provided data set is a list or not
  if(!is.list(dname)){

    stop('Data set provided should be a list containing data sets')

  }

  # Checks for non circular covariates
  if(length(cov) > 0){

    if(!is.vector(cov)){

      stop('Non circular covariates column number should be provided as a vector')

    }
  }

  # Checks for circular covariates
  if(length(cov.circ) > 0){

    if(!is.vector(cov.circ)){

      stop('Circular covariates column number should be provided as a vector')

    }
  }

  # Checks whether any covariate is provided by user or not
  if(!(length(cov) > 0) && !(length(cov.circ) > 0)){

    stop('Atleast a single covariate, either circular or non circular should be provided')

  }

  # Checks for dimension compatibility of weight supplied
  if(length(wgt) > 1){

    if(!(length(wgt) == length(c(cov, cov.circ)))){

      stop('The weight provided should be a single value or vector with weight for each covariate')

    }
  }

  ## data set 2 as a baseline
  dname1_ = list(dname[[1]], dname[[2]])

  ## test set as a baseline
  dname2_ = list(dname[[2]], dname[[1]])
  filelist_ = list(dname1_, dname2_)

  # sequential computation
  `%do%` = foreach::`%do%`
  matcheddata_ = rep(list(), 2)
  foreach::foreach(i = 1:2) %do% {

    matcheddata_[[i]] = CovMatch.Mult(dname = file_list[[i]], cov = cov, wgt = wgt, cov.circ = cov.circ)

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

  return(result_)

}
