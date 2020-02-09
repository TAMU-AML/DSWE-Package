#' @title Covariate Matching
#' @description The function aims to take list of two data sets and returns the after
#' matched data sets using user specified covariates.
#'
#' @param dname This should always be a list, containing all the data sets to match.
#' @param cov Vector of column numbers for the covariates to be used in matching.
#' @param weight Vector of threshold values, against which matching happens.
#' It should be a vector such as c(0.2, 0.2, 0.3), if considering three covariates.
#' @param cov_circ Vector stating the column position of circular variables such as wind direction,
#'  nacelle position etc.
#' @usage CovMatch(dname, cov, weight, cov_circ)
#' @return The function returns a list containing after matched data sets.
#' @export
#' @import foreach

CovMatch = function(dname, cov = NULL, weight, cov_circ = NULL ){

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
  if(length(cov_circ) > 0){

    if(!is.vector(cov_circ)){

      stop('Circular covariates column number should be provided as a vector')

    }
  }

  # Checks whether any covariate is provided by user or not
  if(!(length(cov) > 0) && !(length(cov_circ) > 0)){

    stop('Atleast a single covariate, either circular or non circular should be provided')

  }

  # Checks for dimension compatibility of weight supplied
  if(!(length(weight) == length(c(cov, cov_circ)))){

    stop('The weight provided should be a single value or vector with weight for each covariate')

  }

  ## reference set as a baseline
  # file names to be matched
  dname_1 = rep(list(c()),2)
  dname_1[[1]]= dname[[1]]
  dname_1[[2]]= dname[[2]]

  ## test set as a baseline
  # file names to be matched

  dname_2 = rep(list(c()),2)
  dname_2[[1]]= dname[[2]]
  dname_2[[2]]= dname[[1]]


  file_list = list(dname_1,dname_2)

  # Non circular covariates
  cov = cov

  # circular covariates
  cov_circ = cov_circ

  # weight for threshold calculation
  weight = weight

  # sequential computation
  `%do%` = foreach::`%do%`
  matched_data = rep(list(), 2)
  foreach::foreach(i = 1:2) %do% {

  matched_data[[i]] = covmatch.mult(dname = file_list[[i]], cov = cov, weight = weight, cov_circ = cov_circ)

  }

  ############# Retrieving datasets from 1st matching #########################
  # creating list of matched data set from step 1
  match1 = matched_data[[1]]
  matched1 = rep(list(c()),2)

  matched1[[2]]= match1[[2]]
  matched1[[1]]= match1[[1]]

  ############# Retrieving datasets from 2nd matching #########################
  # creating list of matched data set from step 1
  match2 = matched_data[[2]]
  matched2 = rep(list(c()),2)

  matched2[[2]]= match2[[1]]
  matched2[[1]]= match2[[2]]

  ############ Combining results to generate final matched pairs ###################
  matched = rep(list(c()), 2)
  matched[[1]] = unique(rbind(matched1[[1]], matched2[[1]]))
  matched[[2]] = unique(rbind(matched1[[2]], matched2[[2]]))

  return(matched)

}
