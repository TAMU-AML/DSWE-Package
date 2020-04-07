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

# Compute statistical metric
ComputeSMetric = function(mu1, mu2, band){

  funcDiff = mu1 - mu2
  funcDiff[abs(funcDiff) < band] = 0
  resultP = ((sum(funcDiff) / (sum(mu1 + mu2) / 2)) * (100))
  resultA = (sum(funcDiff)) / nrow(mu2)

  return(resultP)
}

# Compute unweighted metric
ComputeUWMetric = function(mu1, mu2){

  funcDiff = mu1 - mu2
  resultP = ((sum(funcDiff) / (sum(mu1 + mu2) / 2)) * (100))
  resultA = (sum(funcDiff)) / nrow(mu2)

  return(resultP)
}

# Compute weighted metrics
ComputeWMetric = function(dList, mu1, mu2, testdata, testCol){

  mixedData = rbind(dList[[1]], dList[[2]])

  var1Density = density(mixedData[, testCol[1]])
  var1Test = approx(var1Density$x, var1Density$y, xout = testdata[, 1])

  var2Density = density(mixedData[, testCol[2]])
  var2Test = approx(var2Density$x, var2Density$y, xout = testdata[, 2])

  probTest = var1Test$y * var2Test$y / (sum(var1Test$y * var2Test$y))

  resultP = (sum((mu1 - mu2) * (probTest)) / (sum((mu1 + mu2)* (probTest)) / 2)) * 100
  resultA = sum(((mu1 - mu2)) * (probTest))

  return(resultP)

}

# Compute reduction ratio
ComputeRatio = function(dataList1, dataList2, testCol){

  combList1 = rbind(dataList1[[1]], dataList1[[2]])
  combList2 = rbind(dataList2[[1]], dataList2[[2]])

  ratioCol1 = (max(combList2[, testCol[1]]) - min(combList2[, testCol[1]])) / (max(combList1[, testCol[1]]) - min(combList1[, testCol[1]]))
  ratioCol2 = (max(combList2[, testCol[2]]) - min(combList2[, testCol[2]])) / (max(combList1[, testCol[2]]) - min(combList1[, testCol[2]]))

  return(list(ratioCol1 = ratioCol1, ratioCol2 = ratioCol2))
}

# Compute extrapolatiion
ComputeExtrapolation = function(data, yCol, mu1, mu2){

  # creating bins for combined original data
  combData = do.call(rbind, data)
  combData$bin = cut(combData[, yCol], breaks = seq(0, max(combData[, yCol])+50, 50), labels = FALSE)
  combData$bin = 50 * combData$bin

  # calculating probability in original data set
  combDataBinned = combData[, c('bin'), drop = FALSE] %>%  dplyr::group_by(bin) %>% dplyr::summarise(count = length(bin))
  combDataBinned$prob = (combDataBinned$count) / sum(combDataBinned$count)

  # importing result and binning for each period
  result = data.frame(cbind(mu1, mu2))
  colnames(result) = c('mu1', 'mu2')
  result$mu1[result$mu1 < 0] = 0.1
  result$mu2[result$mu2 < 0] = 0.1
  result$bin1 = cut(result$mu1, breaks = seq(0, max(result$mu1)+50, 50), labels = FALSE) * 50
  result$bin2 = cut(result$mu2, breaks = seq(0, max(result$mu2)+50, 50), labels = FALSE) * 50

  binnedPeriod1 = result[, c('bin1', 'mu1')] %>%  dplyr::group_by(bin1) %>% dplyr::summarise(avg = mean(mu1))
  binnedPeriod2 = result[, c('bin2', 'mu2')] %>%  dplyr::group_by(bin2) %>% dplyr::summarise(avg = mean(mu2))

  num = min(length(combDataBinned$bin),length(binnedPeriod1$bin1), length(binnedPeriod2$bin2))
  extrapolatedPwr = (sum((binnedPeriod1$avg[1:num] - binnedPeriod2$avg[1:num]) * combDataBinned$prob[1:num]) * num)

  return(extrapolatedPwr)
}
