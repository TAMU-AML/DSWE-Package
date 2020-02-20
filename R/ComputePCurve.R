ComputePCurve = function(data, xcol, xcol.circ, ycol, var1col, var2col, thrs = 0.2, conflevel = 0.95, gridSize = 50){

  ResultMatching = CovMatch(data, xcol, xcol.circ, thrs)

  TestSet = GenerateTestset(ResultMatching, var1col, var2col, gridSize )

  ResultGP = funGP(ResultMatching, xcol, ycol, conflevel, TestSet)

  ResultSMetric = ComputeSMetric(ResultGP$mu1, ResultGP$mu2, ResultGP$band)

  ResultWMetric = ComputeWMetric(data, ResultMatching, ResultGP$mu1, ResultGP$mu2, TestSet, var1col, var2col)

  returnList = list(muDiff = ResultGP$muDiff, mu2 = ResultGP$diffCov$mu2, mu1 = ResultGP$diffCov$mu1, band = ResultGP$band, confLevel = confLevel, testset = testset, estimatedParams = ResultGP$params, weightedDiff = ResultWMetric, statisticalDiff = ResultSMetric)

  return(returnList)
}
