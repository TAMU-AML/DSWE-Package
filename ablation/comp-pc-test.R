library(DSWE)
library(rjson)
library("rjson")

data = read.csv("../R/Turbine_Upgrade_Dataset/Turbine Upgrade Dataset(VG Pair).csv")

data1 = data[data$upgrade.status == 0, ]
data2 = data[data$upgrade.status == 1, ]
datalist = list(data1, data2)

xcol = c(4,5,6,7)
yCol = 11
xcol.circ = 5
testCol = c(4,5)
gridSize = c(50,50)
rngSeed = c(1, 10, 150, 343, 769, 1001, 1337, 2222, 3456, 5000)
resultMatching = DSWE::CovMatch(datalist, xCol = xcol, xCol.circ = xcol.circ, thrs = 0.2)
testSet = GenerateTestset(resultMatching$matchedData, testCol, gridSize)

optimIdx = fromJSON(file='optimIdx.json')
bandIdx = fromJSON(file='bandIdx.json')
hyperparameters = fromJSON(file='estimatedParams.json')

weighted_diff = vector("list", length=6)
weighted_stat_diff = vector("list", length=6)

## Experiment 1: optimIdx=False, hyperparameters=False, bandIdx=False (Base case)
wd = vector("list", length=10)
wsd = vector("list", length=10)
for (i in 1:length(rngSeed)) {
  resultGP = funGP(datalist = resultMatching$matchedData, xCol = testCol, yCol = yCol, 
                   confLevel = 0.95, testset = testSet, limitMemory = TRUE, 
                   opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), 
                   rngSeed = rngSeed[[i]], hyperparameters = NULL, optimIdx = NULL, bandIdx = NULL)
  wd[[i]] = ComputeWeightedDiff(datalist, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline=1)$percentDiff
  wsd[[i]] = ComputeWeightedStatDiff(datalist, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline=1)
}
weighted_diff[[1]] = as.numeric(wd)
weighted_stat_diff[[1]] = as.numeric(wsd)

## Experiment 2: optimIdx=True, hyperparameters=False, bandIdx=False
wd = vector("list", length=10)
wsd = vector("list", length=10)
for (i in 1:length(rngSeed)) {
  resultGP = funGP(datalist = resultMatching$matchedData, xCol = testCol, yCol = yCol, 
                   confLevel = 0.95, testset = testSet, limitMemory = TRUE, 
                   opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), 
                   rngSeed = rngSeed[[i]], hyperparameters = NULL, optimIdx = optimIdx, bandIdx = NULL)
  wd[[i]] = ComputeWeightedDiff(datalist, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline=1)$percentDiff
  wsd[[i]] = ComputeWeightedStatDiff(datalist, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline=1)
}
weighted_diff[[2]] = as.numeric(wd)
weighted_stat_diff[[2]] = as.numeric(wsd)

## Experiment 3: optimIdx=True, hyperparameters=False, bandIdx=True
wd = vector("list", length=10)
wsd = vector("list", length=10)
for (i in 1:length(rngSeed)) {
  resultGP = funGP(datalist = resultMatching$matchedData, xCol = testCol, yCol = yCol, 
                   confLevel = 0.95, testset = testSet, limitMemory = TRUE, 
                   opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), 
                   rngSeed = rngSeed[[i]], hyperparameters = NULL, optimIdx = optimIdx, bandIdx = bandIdx)
  wd[[i]] = ComputeWeightedDiff(datalist, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline=1)$percentDiff
  wsd[[i]] = ComputeWeightedStatDiff(datalist, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline=1)
}
weighted_diff[[3]] = as.numeric(wd)
weighted_stat_diff[[3]] = as.numeric(wsd)

## Experiment 4: optimIdx=False, hyperparameters=False, bandIdx=True
wd = vector("list", length=10)
wsd = vector("list", length=10)
for (i in 1:length(rngSeed)) {
  resultGP = funGP(datalist = resultMatching$matchedData, xCol = testCol, yCol = yCol, 
                   confLevel = 0.95, testset = testSet, limitMemory = TRUE, 
                   opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), 
                   rngSeed = rngSeed[[i]], hyperparameters = NULL, optimIdx = NULL, bandIdx = bandIdx)
  wd[[i]] = ComputeWeightedDiff(datalist, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline=1)$percentDiff
  wsd[[i]] = ComputeWeightedStatDiff(datalist, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline=1)
}
weighted_diff[[4]] = as.numeric(wd)
weighted_stat_diff[[4]] = as.numeric(wsd)

## Experiment 5: optimIdx=False, hyperparameters=True, bandIdx=False
wd = vector("list", length=10)
wsd = vector("list", length=10)
for (i in 1:length(rngSeed)) {
  resultGP = funGP(datalist = resultMatching$matchedData, xCol = testCol, yCol = yCol, 
                   confLevel = 0.95, testset = testSet, limitMemory = TRUE, 
                   opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), 
                   rngSeed = rngSeed[[i]], hyperparameters = hyperparameters, optimIdx = NULL, bandIdx = NULL)
  wd[[i]] = ComputeWeightedDiff(datalist, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline=1)$percentDiff
  wsd[[i]] = ComputeWeightedStatDiff(datalist, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline=1)
}
weighted_diff[[5]] = as.numeric(wd)
weighted_stat_diff[[5]] = as.numeric(wsd)


## Experiment 6: optimIdx=False, hyperparameters=True, bandIdx=True
wd = vector("list", length=10)
wsd = vector("list", length=10)
for (i in 1:length(rngSeed)) {
  resultGP = funGP(datalist = resultMatching$matchedData, xCol = testCol, yCol = yCol, 
                   confLevel = 0.95, testset = testSet, limitMemory = TRUE, 
                   opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), 
                   rngSeed = rngSeed[[i]], hyperparameters = hyperparameters, optimIdx = NULL, bandIdx = bandIdx)
  wd[[i]] = ComputeWeightedDiff(datalist, resultGP$mu1, resultGP$mu2, testSet, testCol, baseline=1)$percentDiff
  wsd[[i]] = ComputeWeightedStatDiff(datalist, resultGP$mu1, resultGP$mu2, resultGP$band, testSet, testCol, baseline=1)
}
weighted_diff[[6]] = as.numeric(wd)
weighted_stat_diff[[6]] = as.numeric(wsd)


result = data.frame(
  Experiment = c(1, 2, 3, 4, 5, 6),
  optimIdx = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
  estimatedParams = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  bandIdx = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE),
  weighted_diff_mean = round(unlist(lapply(weighted_diff, mean)), 2),
  weighted_diff_std = round(unlist(lapply(weighted_diff, sd)), 2),
  weighted_stat_diff_mean = round(unlist(lapply(weighted_stat_diff, mean)), 2),
  weighted_stat_diff_std = round(unlist(lapply(weighted_stat_diff, sd)), 2),
  weighted_diff = paste(as.matrix(weighted_diff), ", "),
  weighted_stat_diff = paste(as.matrix(weighted_stat_diff), ", ")
)

write.csv(result, file='outputs.csv', row.names = FALSE)
