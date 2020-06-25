ComputeWeightedDifference = function(muDiff, weights, base, statDiff = FALSE, confBand = NULL){
  if (!as.numeric(muDiff)){
    stop('muDiff must be a numeric vector' )
  }
  if (!as.numeric(weights)){
    stop('weights must be a numeric vector')
  }
  if (!as.numeric(base)){
    stop('base must be a numeric vector')
  }
  if (length(muDiff) != length(weights) || length(muDiff) != length(base)){
    stop('length of muDiff, weights and base must be the same')
  }
  if(sum(weights) != 1){
    stop('the weights must sum to 1 for being a valid weights')
  }
  if (!is.logical(statDiff)){
    stop('statDiff must be either TRUE or FALSE')
  } else if (statDiff){
    if (is.null(confBand)){
      stop('confBand must be provided when statDiff is set to TRUE')
    }
    if (!as.numeric(confBand)){
      stop('confBand must be a numeric vector when statDiff is set to TRUE')
    }
    if (length(muDiff) != length(confBand)){
      stop('length of confBand must be the same as muDiff')
    }
  }
  weightedDiff = computeWeightedDiffExtern(muDiff, weights, base, statDiff, confBand)
}
