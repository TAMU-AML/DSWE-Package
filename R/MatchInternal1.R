#' @import matrixStats


covmatch.mult = function(dname, weight, cov_circ = NULL ){

  # Store data sets to be compared
  f_name = dname

  # Weight assigning
  wgt = weight

  # Circular variable presence indicator
  flag = 0

  # Ensuring circular variable to be between 0 to 360 degree
  if(length(cov_circ) > 0) {

    f_name = lapply(1:length(f_name), function(x) circ.positive(f_name[[x]], cov_circ))

    # Circular variable indicator
    flag = 1

  }

  # Setting up the reference data set and threshold
  ref_id = length(f_name)
  ref = as.matrix(f_name[[ref_id]])

  # Test files
  test_id = c(1:length(f_name))[-ref_id]

  # Setting up thresholds
  ratio = matrixStats::colSds(as.matrix(ref)) / colMeans(ref)
  thres = ratio * wgt

  # Matching data sets with ref as reference
  matchID  = lapply(test_id, function(x) match.cov(ref, f_name[[x]], thres, circ_pos = cov_circ, flag = flag))

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
