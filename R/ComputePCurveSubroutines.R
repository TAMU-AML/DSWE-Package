# Generates test set
GenerateTestset = function(data, var1, var2, gridSize){

  Var1Min = max(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, var1], c(0.025,0.975))[1])))
  Var1Max = min(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, var1], c(0.025,0.975))[2])))
  Var1Range = seq(Var1Min, Var1Max, length.out = gridSize )

  Var2Min = max(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, var2], c(0.025,0.975))[1])))
  Var2Max = min(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, var2], c(0.025,0.975))[2])))
  Var2Range = seq(Var1Min, Var1Max, length.out = gridSize )

  return(expand.grid(Var1Range, Var2Range))

}

# Compute statistical metric
ComputeSMetric = function(mu1, mu2, band){

  funcDiff = mu1 - mu2
  funcDiff[funcDiff < band] = 0
  resultP = ((sum(funcDiff) / sum(mu2)) * (100))
  resultA = (sum(funcDiff)) / nrow(mu2)

  return(resultP)
}

# Compute weighted metrics
ComputeWMetric = function(dList1, dList2, mu1, mu2, testdata, var1, var2){

  mixedData = rbind(dList1[[1]], dList[[2]])

  var1Density = density(mixed_data[, var1])
  var1Test = approx(var1Density$x,var1Density$y, xout = testdata[, 1])

  var2Density = density(mixed_data[, var2])
  var2Test = approx(var2Density$x,var2Density$y, xout = testdata[, 2])

  probTest = var1Test$y * var2Test$y / (sum(var1Test$y * var2Test$y))

  resultP = (sum((mu1 - mu2) * (probTest)) / sum((mu2)* (probTest))) * 100
  resultA = sum(((mu1 - mu2)) * (probTest))

  return(resultP)

}
