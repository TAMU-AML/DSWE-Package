#' @useDynLib DSWE
#' @importFrom Rcpp sourceCpp

######### function to match data set##################
# match.cov : main matching function given data sets
# ref : data set considered as a reference
# obj : data set whose each observation searches for a matching observation in ref
# thres : threshold against which similarity is obtained between ref and obj data sets
# circ_pos : position of circular variables (such as wind direction etc), if supplied by user
match.cov = function(ref, obj, thres, circ_pos = 0, flag = 0){

  # Ensuring that the data sets are converted to matrix
  ref = as.matrix(ref)
  obj = as.matrix(obj)

  match = matchcov(ref = ref , obj = obj, thres = thres, circ_pos = circ_pos, flag = flag)

  # Returns matched index
  return(match)

}


######### Function to convert circular variable values to positive ##################
# data : data set consiting of features and vlues
# circ : position of circular variable supplied by user
circ.positive = function(data, circ){

  # Iterating over circular variabe
  for(i in circ){

    while(sum(data[ , i] < 0) != 0){

      data[, i] = data[, i] + 360

    }

  }

  # Returns manipulated data
  return(data)
}



