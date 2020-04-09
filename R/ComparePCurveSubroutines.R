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


# Compute statistical weighted metrics
ComputeWSMetric = function(dList, mu1, mu2, band, testdata, testCol){

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

# Compute weighted statistical extrapolatiion
ComputeWSExtrapolation = function(data, yCol, mu1, mu2, band){

  # creating bins for combined original data
  combData = do.call(rbind, data)
  combData$bin = cut(combData[, yCol], breaks = seq(0, max(combData[, yCol])+100, 100), labels = FALSE)
  combData$bin = 100 * combData$bin

  # calculating probability in original data set
  combDataBinned = combData[, c('bin'), drop = FALSE] %>%  dplyr::group_by(bin) %>% dplyr::summarise(count = length(bin))
  combDataBinned$prob = (combDataBinned$count) / sum(combDataBinned$count)

  # importing result and checking bins for each period
  result = data.frame(cbind(mu1, mu2))
  colnames(result) = c('mu1', 'mu2')
  result$mu1[result$mu1 < 0] = 0.1
  result$mu2[result$mu2 < 0] = 0.1
  result$agg = (result$mu1 + result$mu2) / 2
  result$delta = result$mu1 - result$mu2
  result$delta[abs(result$delta) <= band] = 0
  combDataBinned[, c('avg1', 'delta1', 'avg2', 'delta2')] = NA

  for (i in 1:nrow(combDataBinned)){

    if( i == 1){
      combDataBinned$avg1[i] = mean(result$agg[result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$delta1[i] = mean(result$delta[result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$avg2[i] = mean(result$agg[result$mu2 <= combDataBinned$bin[i]])
      combDataBinned$delta2[i] = mean(result$delta[result$mu2 <= combDataBinned$bin[i]])
    }else{

      combDataBinned$avg1[i] = mean(result$agg[result$mu1 > combDataBinned$bin[i-1] & result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$delta1[i] = mean(result$delta[result$mu1 > combDataBinned$bin[i-1] & result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$avg2[i] = mean(result$agg[result$mu2 > combDataBinned$bin[i-1] & result$mu2 <= combDataBinned$bin[i]])
      combDataBinned$delta2[i] = mean(result$delta[result$mu2 > combDataBinned$bin[i-1] & result$mu2 <= combDataBinned$bin[i]])
    }
  }

  while(sum(is.na(combDataBinned)) != 0){

    id = max(unique(which(is.na(combDataBinned), arr.ind=TRUE)[, 1]))
    combDataBinned[is.na(combDataBinned)] = 0
    combDataBinned[min(id)-1, 'prob'] = sum(combDataBinned[c(min(id-1), id), 'prob'])
    combDataBinned[min(id)-1, c('avg1', 'delta1', 'avg2', 'delta2')] = apply(combDataBinned[c(min(id-1), id), c('avg1', 'delta1', 'avg2', 'delta2')], 2, mean)

  }

  combDataBinned = combDataBinned[-id, ]
  combDataBinned$avg = (combDataBinned$avg1 + combDataBinned$avg2) / 2
  combDataBinned$delta = (combDataBinned$delta1 + combDataBinned$delta2) / 2

  extrapolatedDelta = sum(combDataBinned$delta * combDataBinned$prob)
  extrapolatedDeltaP = round((extrapolatedDelta / sum(combDataBinned$avg * combDataBinned$prob)) * 100, 2)
  return(extrapolatedDeltaP)
}


# Compute weighted metrics
ComputeWMetric = function(dList, mu1, mu2, testdata, testCol){

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

# Compute weighted extrapolatiion
ComputeWExtrapolation = function(data, yCol, mu1, mu2){

  # creating bins for combined original data
  combData = do.call(rbind, data)
  combData$bin = cut(combData[, yCol], breaks = seq(0, max(combData[, yCol])+100, 100), labels = FALSE)
  combData$bin = 100 * combData$bin

  # calculating probability in original data set
  combDataBinned = combData[, c('bin'), drop = FALSE] %>%  dplyr::group_by(bin) %>% dplyr::summarise(count = length(bin))
  combDataBinned$prob = (combDataBinned$count) / sum(combDataBinned$count)

  # importing result and checking bins for each period
  result = data.frame(cbind(mu1, mu2))
  colnames(result) = c('mu1', 'mu2')
  result$mu1[result$mu1 < 0] = 0.1
  result$mu2[result$mu2 < 0] = 0.1
  result$agg = (result$mu1 + result$mu2) / 2
  result$delta = result$mu1 - result$mu2
  combDataBinned[, c('avg1', 'delta1', 'avg2', 'delta2')] = NA

  for (i in 1:nrow(combDataBinned)){

    if( i == 1){
      combDataBinned$avg1[i] = mean(result$agg[result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$delta1[i] = mean(result$delta[result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$avg2[i] = mean(result$agg[result$mu2 <= combDataBinned$bin[i]])
      combDataBinned$delta2[i] = mean(result$delta[result$mu2 <= combDataBinned$bin[i]])
    }else{

      combDataBinned$avg1[i] = mean(result$agg[result$mu1 > combDataBinned$bin[i-1] & result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$delta1[i] = mean(result$delta[result$mu1 > combDataBinned$bin[i-1] & result$mu1 <= combDataBinned$bin[i]])
      combDataBinned$avg2[i] = mean(result$agg[result$mu2 > combDataBinned$bin[i-1] & result$mu2 <= combDataBinned$bin[i]])
      combDataBinned$delta2[i] = mean(result$delta[result$mu2 > combDataBinned$bin[i-1] & result$mu2 <= combDataBinned$bin[i]])
    }
  }

  while(sum(is.na(combDataBinned)) != 0){


    id = max(unique(which(is.na(combDataBinned), arr.ind=TRUE)[, 1]))
    combDataBinned[is.na(combDataBinned)] = 0
    combDataBinned[min(id)-1, 'prob'] = sum(combDataBinned[c(min(id-1), id), 'prob'])
    combDataBinned[min(id)-1, c('avg1', 'delta1', 'avg2', 'delta2')] = apply(combDataBinned[c(min(id-1), id), c('avg1', 'delta1', 'avg2', 'delta2')], 2, mean)

  }

  combDataBinned = combDataBinned[-id, ]
  combDataBinned$avg = (combDataBinned$avg1 + combDataBinned$avg2) / 2
  combDataBinned$delta = (combDataBinned$delta1 + combDataBinned$delta2) / 2

  extrapolatedDelta = sum(combDataBinned$delta * combDataBinned$prob)
  extrapolatedDeltaP = round((extrapolatedDelta / sum(combDataBinned$avg * combDataBinned$prob)) * 100, 2)
  return(extrapolatedDeltaP)
}

# Compute reduction ratio
ComputeRatio = function(dataList1, dataList2, testCol){

  combList1 = rbind(dataList1[[1]], dataList1[[2]])
  combList2 = rbind(dataList2[[1]], dataList2[[2]])

  ratioCol1 = (max(combList2[, testCol[1]]) - min(combList2[, testCol[1]])) / (max(combList1[, testCol[1]]) - min(combList1[, testCol[1]]))
  ratioCol2 = (max(combList2[, testCol[2]]) - min(combList2[, testCol[2]])) / (max(combList1[, testCol[2]]) - min(combList1[, testCol[2]]))

  return(list(ratioCol1 = ratioCol1, ratioCol2 = ratioCol2))
}
