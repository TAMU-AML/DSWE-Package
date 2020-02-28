#' @import matrixStats
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

    fname_ = lapply(1:length(fname_), function(x) Circ.Positive(fname_[[x]], cov.circ))

    # Circular variable after data subsetting position
    pos = which(cov.circ == covcol_)

    # Circular variable indicator
    flag = 1

  }

  # Setting up the reference data set and threshold
  refid_ = length(fname_)
  ref_ = as.matrix(fname_[[refid_]][, covcol_, drop = F])

  # Test files
  testid_ = c(1:length(fname_))[-refid_]

  # Setting up thresholds
  ratio_ = matrixStats::colSds(as.matrix(ref_)) / colMeans(ref_)
  thres_ = ratio_ * wgt_

  # Matching data sets with ref as reference
  matchID_  = lapply(testid_, function(x) Match.Cov(ref_, fname_[[x]][, covcol_, drop = F], thres_, pos, flag))

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
