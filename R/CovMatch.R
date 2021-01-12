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

#' @title Covariate Matching
#' @description The function aims to take list of two data sets and returns the after
#' matched data sets using user specified covariates and threshold
#'
#' @param data a list, consisting of data sets to match, also each of the individual data set can be dataframe or a matrix
#' @param thrs a numerical or a vector of threshold values for each covariates, against which matching happens
#' It should be a single value or a vector of values representing threshold for each of the covariate
#' @param xCol a vector stating the column position of covariates used
#' @param xCol.circ a vector stating the column position of circular variables
#' @param priority a boolean, default value False, otherwise computes the sequence of matching
#' @usage CovMatch(data, xCol, xCol.circ, thrs, priority)
#' @return a list containing :
#'   \itemize{
#'   \item originalData - The data sets provided for matching
#'   \item matchedData - The data sets after matching
#'   \item MinMaxOriginal - The minimum and maximum value in original data for each covariate used in matching
#'   \item MinMaxMatched - The minimum and maximum value in matched data for each covariates used in matching
#'}
#' @examples 
#' 
#' data1 = data1[1:100, ]
#' data2 = data2[1:100, ]

#' data = list(data1, data2)
#' xCol = 2
#' xCol.circ = NULL
#' thrs = 0.1
#' priority = FALSE

#' matched_data = CovMatch(data, xCol, xCol.circ, thrs, priority)
#' 
#' @references Ding, Y. (2019). Data Science for Wind Energy. Chapman & Hall, Boca Raton, FL.
#' @export

CovMatch = function(data, xCol, xCol.circ = NULL, thrs = 0.2, priority = FALSE){

  # Checks whether the provided data set is a list or not
  if(!is.list(data)){

    stop('data must be a list containing data sets')

  }

  if(length(data) != 2){


    stop('The number of data sets to match should be equal to two')

  }

  # Checks for non circular covariates
  if(!is.vector(xCol)){

    stop('xCol must be provided as numeric/vector')

  }else{

    if(!(all(xCol %in% 1:ncol(data[[1]])))){

      stop('xCol value should not be more than the column in data')
    }
  }

  # Checks for circular covariates
  if(length(xCol.circ) > 0){

    if(!is.vector(xCol.circ)){

      stop('xCol.circ must be provided as a numeric/vector')

    }
  }


  # Checks for dimension compatibility of weight supplied
  if(length(thrs) > 1){

    if(!(length(thrs) == length(xCol))){

      stop('The thrs must be a single value or vector with weight for each covariate')

    }
  }

  if(priority == TRUE){

    covDiff = as.numeric(abs(colMeans(data[[1]][, xCol]) - colMeans(data[[2]][, xCol])))

    xCol = xCol[order(covDiff, decreasing = TRUE)]

  }

  if(!is.null(xCol.circ)) {

    if(!xCol.circ %in% xCol){

      stop('xCol.circ should be present in xCol')
    }else{

      data = lapply(1:length(data), function(x) Circ.Positive(data[[x]], xCol.circ))
    }
  }

  ## data set 2 as a baseline
  dname1_ = list(data[[1]], data[[2]])

  ## test set as a baseline
  dname2_ = list(data[[2]], data[[1]])
  filelist_ = list(dname1_, dname2_)

  # sequential computation
  matcheddata_ = rep(list(), 2)
  for(i in 1:2){

    matcheddata_[[i]] = CovMatch.Mult(filelist_[[i]], xCol, thrs, xCol.circ)
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
  result_[[1]] = unique(rbind(matched1_[[2]], matched2_[[2]]))
  result_[[2]] = unique(rbind(matched1_[[1]], matched2_[[1]]))

  MinMaxOriginal = MinMaxData(data, xCol)
  MinMaxMatched = MinMaxData(result_, xCol)

  return(list(originalData = data,  matchedData = result_, MinMaxOriginal = MinMaxOriginal, MinMaxMatched = MinMaxMatched))
}
