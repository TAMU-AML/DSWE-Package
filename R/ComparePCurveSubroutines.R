# Generates test set
GenerateTestset = function(data, testCol, gridSize){

  Var1Min = max(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, testCol[1]], c(0.025,0.975))[1])))
  Var1Max = min(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, testCol[1]], c(0.025,0.975))[2])))
  Var1Range = seq(Var1Min, Var1Max, length.out = gridSize[1] )

  Var2Min = max(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, testCol[2]], c(0.025,0.975))[1])))
  Var2Max = min(unlist(lapply(c(1:length(data)), function(x) quantile(data[[x]][, testCol[2]], c(0.025,0.975))[2])))
  Var2Range = seq(Var2Min, Var2Max, length.out = gridSize[2] )

  return(expand.grid(Var1Range, Var2Range))

}

# Compute unweigted difference between the functions
ComputeDiff = function(mu1, mu2){

  diff = sum(mu1 - mu2)/length(mu1)
  avgMu = mean(mu1 + mu2)/2
  percentDiff = diff*100/avgMu
  return(percentDiff)

}

# Compute statistically significant unweigted difference between the functions
ComputeStatDiff = function(mu1, mu2, band){

  diffMu = sum(mu1 - mu2)
  diffMu[which(abs(diffMu)<band)] = 0
  diff = sum(diffMu)/lenght(mu1)
  avgMu = mean(mu1 + mu2)/2
  percentDiff = diff*100/avgMu
  return(percentDiff)

}


# Compute weighted metrics
ComputeWeightedDiff = function(dList, mu1, mu2, testdata, testCol){

  mixedData = rbind(dList[[1]], dList[[2]])

  var1Density = density(mixedData[, testCol[1]])
  var1Test = approx(var1Density$x, var1Density$y, xout = testdata[, 1])

  var2Density = density(mixedData[, testCol[2]])
  var2Test = approx(var2Density$x, var2Density$y, xout = testdata[, 2])

  probTest = var1Test$y * var2Test$y / (sum(var1Test$y * var2Test$y))

  resultP = round((sum((mu1 - mu2) * (probTest)) / (sum((mu1 + mu2)* (probTest)) / 2)) * 100, 2)
  resultA = sum(((mu1 - mu2)) * (probTest))

  return(resultP)

}

# Compute statistical weighted metrics
ComputeWeightedStatDiff = function(dList, mu1, mu2, band, testdata, testCol){

  mixedData = rbind(dList[[1]], dList[[2]])

  var1Density = density(mixedData[, testCol[1]])
  var1Test = approx(var1Density$x, var1Density$y, xout = testdata[, 1])

  var2Density = density(mixedData[, testCol[2]])
  var2Test = approx(var2Density$x, var2Density$y, xout = testdata[, 2])

  probTest = var1Test$y * var2Test$y / (sum(var1Test$y * var2Test$y))

  funcDiff = mu1 - mu2
  funcDiff[abs(funcDiff) <= band] = 0
  resultP = round((sum((funcDiff) * (probTest)) / (sum((mu1 + mu2)* (probTest)) / 2)) * 100, 2)
  resultA = sum(((funcDiff)) * (probTest))

  return(resultP)

}

# Compute difference in the function scaled to the orginal data
ComputeScaledDiff = function(datalist, yCol, mu1, mu2){
  mergedData = rbind(datalist[[1]],datalist[[2]])
  pw = mergedData[,yCol]
  pw[pw < 0] = 0
  binWidth = 100

  #starting bin value
  start = 0

  #ending value based on max power rounded up to the closest multiple of the bin width.
  #This formulation ensures that the last bin of the original data will never be empty.
  end = max(pw) + (binWidth - (max(pw)%%binWidth))
  bins = seq(start,end,binWidth)

  #Merge the empty bins till there are no empty bins in any of the three: pw, mu1, mu2
  cumCount = sapply(c(1:length(bins)), function (x) length(which(pw < bins[x]))) #cumulative count for each bin
  nonEmptyBins = bins[!duplicated(cumCount)] #remove bins with cumulative count same as the previous bin i.e. remove empty bins

  #Repeat the above two lines for mu1 and mu2
  cumCount = sapply(c(1:length(nonEmptyBins)), function (x) length(which(mu1 < nonEmptyBins[x])))
  nonEmptyBins = nonEmptyBins[!duplicated(cumCount)]
  cumCount = sapply(c(1:length(nonEmptyBins)), function (x) length(which(mu2 < nonEmptyBins[x])))
  nonEmptyBins = nonEmptyBins[!duplicated(cumCount)]

  #Check if the last bin is equal to 'end'. If not, replace last bin value with 'end'
  if (max(nonEmptyBins)<end){
    nonEmptyBins[length(nonEmptyBins)] = end
  }

  #Compute probability for each merged bin
  totalCount = length(pw)
  probVector = sapply(c(2:length(nonEmptyBins)), function (x) length(which(pw >= nonEmptyBins[x-1] & pw < nonEmptyBins[x]))/totalCount)

  #pointwise delta and average mu
  delta = mu1 - mu2
  mu = 0.5*(mu1+mu2)

  #bin wise average delta and mu
  avgDelta1 = sapply(c(2:length(nonEmptyBins)), function (x) mean(delta[(which(mu1 >= nonEmptyBins[x-1] & mu1 < nonEmptyBins[x]))]))
  avgMu1 = sapply(c(2:length(nonEmptyBins)), function (x) mean(mu[(which(mu1 >= nonEmptyBins[x-1] & mu1 < nonEmptyBins[x]))]))
  avgDelta2 = sapply(c(2:length(nonEmptyBins)), function (x) mean(delta[(which(mu2 >= nonEmptyBins[x-1] & mu2 < nonEmptyBins[x]))]))
  avgMu2 = sapply(c(2:length(nonEmptyBins)), function (x) mean(mu[(which(mu2 >= nonEmptyBins[x-1] & mu2 < nonEmptyBins[x]))]))

  deltaBin = 0.5*(avgDelta1+avgDelta2)
  muBin = 0.5*(avgMu1+avgMu2)

  scaledDiff = t(probVector)%*%deltaBin
  percentScaledDiff = scaledDiff*100/(t(probVector)%*%muBin)

  return(as.numeric(percentScaledDiff))
}


# Compute statistically significant difference in the function scaled to the orginal data
ComputeScaledStatDiff = function(datalist, yCol, mu1, mu2, band){
  mergedData = rbind(datalist[[1]],datalist[[2]])
  pw = mergedData[,yCol]
  pw[pw < 0] = 0
  binWidth = 100

  #starting bin value
  start = 0

  #ending value based on max power rounded up to the closest multiple of the bin width.
  #This formulation ensures that the last bin of the original data will never be empty.
  end = max(pw) + (binWidth - (max(pw)%%binWidth))
  bins = seq(start,end,binWidth)

  #Merge the empty bins till there are no empty bins in any of the three: pw, mu1, mu2
  cumCount = sapply(c(1:length(bins)), function (x) length(which(pw < bins[x]))) #cumulative count for each bin
  nonEmptyBins = bins[!duplicated(cumCount)] #remove bins with cumulative count same as the previous bin i.e. remove empty bins

  #Repeat the above two lines for mu1 and mu2
  cumCount = sapply(c(1:length(nonEmptyBins)), function (x) length(which(mu1 < nonEmptyBins[x])))
  nonEmptyBins = nonEmptyBins[!duplicated(cumCount)]
  cumCount = sapply(c(1:length(nonEmptyBins)), function (x) length(which(mu2 < nonEmptyBins[x])))
  nonEmptyBins = nonEmptyBins[!duplicated(cumCount)]

  #Check if the last bin is equal to 'end'. If not, replace last bin value with 'end'
  if (max(nonEmptyBins)<end){
    nonEmptyBins[length(nonEmptyBins)] = end
  }

  #Compute probability for each merged bin
  totalCount = length(pw)
  probVector = sapply(c(2:length(nonEmptyBins)), function (x) length(which(pw >= nonEmptyBins[x-1] & pw < nonEmptyBins[x]))/totalCount)

  #pointwise delta and average mu
  delta = mu1 - mu2
  delta[which(abs(delta)<band)] = 0
  mu = 0.5*(mu1+mu2)

  #bin wise average delta and mu
  avgDelta1 = sapply(c(2:length(nonEmptyBins)), function (x) mean(delta[(which(mu1 >= nonEmptyBins[x-1] & mu1 < nonEmptyBins[x]))]))
  avgMu1 = sapply(c(2:length(nonEmptyBins)), function (x) mean(mu[(which(mu1 >= nonEmptyBins[x-1] & mu1 < nonEmptyBins[x]))]))
  avgDelta2 = sapply(c(2:length(nonEmptyBins)), function (x) mean(delta[(which(mu2 >= nonEmptyBins[x-1] & mu2 < nonEmptyBins[x]))]))
  avgMu2 = sapply(c(2:length(nonEmptyBins)), function (x) mean(mu[(which(mu2 >= nonEmptyBins[x-1] & mu2 < nonEmptyBins[x]))]))

  deltaBin = 0.5*(avgDelta1+avgDelta2)
  muBin = 0.5*(avgMu1+avgMu2)

  scaledDiff = t(probVector)%*%deltaBin
  percentScaledDiff = scaledDiff*100/(t(probVector)%*%muBin)

  return(as.numeric(percentScaledDiff))
}


# Compute reduction ratio
ComputeRatio = function(dataList1, dataList2, testCol){

  combList1 = rbind(dataList1[[1]], dataList1[[2]])
  combList2 = rbind(dataList2[[1]], dataList2[[2]])

  ratioCol1 = (max(combList2[, testCol[1]]) - min(combList2[, testCol[1]])) / (max(combList1[, testCol[1]]) - min(combList1[, testCol[1]]))
  ratioCol2 = (max(combList2[, testCol[2]]) - min(combList2[, testCol[2]])) / (max(combList1[, testCol[2]]) - min(combList1[, testCol[2]]))

  return(list(ratioCol1 = ratioCol1, ratioCol2 = ratioCol2))
}
