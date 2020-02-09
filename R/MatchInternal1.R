#' @import matrixStats


covmatch.mult = function(dname, cov = NULL, weight, cov_circ = NULL ){

  # Checks for number of data sets provided by user
  if(length(dname) < 2){

    stop('Atleast two data sets must be provided for matching')

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

  # Store data sets to be compared
  f_name = dname

  # Weight assigning
  wgt = weight

  # Covariates column number for matching
  cov_col = c(cov, cov_circ)

  # Circular variable position indicator
  pos = 0

  # Circular variable presence indicator
  flag = 0

  # Ensuring circular variable to be between 0 to 360 degree
  if(length(cov_circ) > 0) {

    f_name = lapply(1:length(f_name), function(x) circ.positive(f_name[[x]], cov_circ))

    # Circular variable after data subsetting position
    pos = (length(cov)+1):length(cov_col)

    # Circular variable indicator
    flag = 1

  }

  # Setting up the reference data set and threshold
  ref_id = length(f_name)
  ref = as.matrix(f_name[[ref_id]][, cov_col, drop = F])

  # Test files
  test_id = c(1:length(f_name))[-ref_id]

  # Setting up thresholds
  ratio = matrixStats::colSds(as.matrix(ref)) / colMeans(ref)
  thres = ratio * wgt

  # Matching data sets with ref as reference
  matchID  = lapply(test_id, function(x) match.cov(ref, f_name[[x]][, cov_col, drop = F], thres, circ_pos = pos, flag = flag))

  # creating list of matched data set
  matched = rep(list(c()), (length(f_name)))

  # selecting indices of matched data sets
  # matched reference set
  refid = ((matchID[[1]]) > 0)
  if(length(f_name) < 3){

    refid = refid

  } else
  {

    for(i in 2:(length(f_name)-1))
    {
      refid = refid & ((matchID[[i]]) > 0)

    }

  }

  refID = which(refid)
  matched[[ref_id]] = f_name[[ref_id]][refID, ]

  # matched test set
  for(j in (1:(length(f_name)-1))){

    matched[[test_id[j]]] = f_name[[test_id[j]]][matchID[[j]][refID], ]

  }

  return(matched)
}
