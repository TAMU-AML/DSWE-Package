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

#' @importFrom  matrixStats colMins colMaxs colSds
#' @useDynLib DSWE
#' @importFrom Rcpp sourceCpp

CovMatch.Mult = function(dname, cov, wgt, cov.circ){

  # Store data sets to be compared
  fname_ = dname

  # Weight assigning
  wgt_ = wgt

  # Covariates column number for matching
  covcol_ = cov

  # Circular variable position indicator
  pos = 0

  # Circular variable presence indicator
  flag = 0

  # Ensuring circular variable to be between 0 to 360 degree
  if(length(cov.circ) > 0) {

    # Circular variable after data subsetting position
    pos = which(cov.circ == covcol_)

    # Circular variable indicator
    flag = 1

  }

  # Setting up the reference data set and threshold
  refid_ = length(fname_)
  ref_ = as.matrix(fname_[[refid_]][, covcol_, drop = FALSE])

  # Test files
  testid_ = c(1:length(fname_))[-refid_]

  # Setting up thresholds
  ratio_ = colSds(as.matrix(ref_))
  thres_ = ratio_ * wgt_

  # Matching data sets with ref as reference
  matchID_  = lapply(testid_, function(x) Match.Cov(ref_, fname_[[x]][, covcol_, drop = FALSE], thres_, pos, flag))

  # creating list of matched data set
  matched_ = rep(list(c()), (length(fname_)))

  # selecting indices of matched data sets
  # matched reference set
  ref.id_ = ((matchID_[[1]]) > 0)
  if(length(fname_) < 3){

    ref.id_ = ref.id_

  } else{

    for(i in 2:(length(fname_)-1))
    {
      ref.id_ = ref.id_ & ((matchID_[[i]]) > 0)

    }

  }

  refID_ = which(ref.id_)
  matched_[[refid_]] = fname_[[refid_]][refID_, ]

  # matched test set
  for(j in (1:(length(fname_)-1))){

    matched_[[testid_[j]]] = fname_[[testid_[j]]][matchID_[[j]][refID_], ]

  }

  return(matched_)
}


Match.Cov = function(ref, obj, thres, circ.pos = 0, flag = 0){

  # Ensuring that the data sets are converted to matrix
  ref = as.matrix(ref)
  obj = as.matrix(obj)

  match = matchcov(ref = ref , obj = obj, thres = thres, circ_pos = circ.pos, flag = flag)

  # Returns matched index
  return(match)

}


######### Function to convert circular variable values to positive ##################
# data : data set consiting of features and vlues
# circ : position of circular variable supplied by user
Circ.Positive = function(data, circ){

  # Iterating over circular variabe
  for(i in circ){

    while(sum(data[ , i] < 0) != 0){

      data[, i] = data[, i] + 360

    }

  }

  # Returns manipulated data
  return(data)
}


MinMaxData = function(data, xcol){

  if(is.list(data)){

    comData = rbind(data[[1]], data[[2]])
  }
  filteredData = matrix(nrow = 2, ncol = length(xcol))
  colnames(filteredData) = as.character(xcol)
  row.names(filteredData) = c('Min', 'Max')

  filteredData[1 , ] = colMins(as.matrix(comData[, xcol]))
  filteredData[2 , ] = colMaxs(as.matrix(comData[, xcol]))

  return(filteredData)

}
