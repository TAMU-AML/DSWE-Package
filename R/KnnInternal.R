# Computes knn based rmse for given covariates
KNN.Internal1 = function(data, xcol, ycol, kfold){

  set.seed(1)
  data = data[sample(1:nrow(data), nrow(data), replace = FALSE), ]
  folds = cut(1:nrow(data), breaks = kfold, labels = FALSE)
  rmse = 0


  for(i in 1:kfold)
  {

    train_data = data[folds != i, ]
    test_data = data[folds == i, ]
    folds2 = cut(1:nrow(train_data), breaks = kfold, labels = FALSE)

    best_k = computeBestK(train_data[, xcol], train_data[, ycol], seq(5, 50, 5))

    test_pred = knn.reg(as.matrix(train_data[, xcol]), as.matrix(test_data[, xcol]), train_data[, ycol], k = best_k)
    residual = sqrt((sum((test_data[, ycol] - test_pred$pred)^2)) / length(test_data[, ycol]))
    rmse = rmse + residual

  }

  avg_rmse = round(rmse / kfold, 2)

  return(avg_rmse)

}

# Computes knn based forward stepwise : model selection
KNN.Internal2 = function(data, xcol, ycol,kfold){

  dataframe = data.frame(array(0, dim = c(length(xcol),2)))
  best.covariate = NULL

  for(itr in 1:(length(xcol))){
    best.rmse = 1000
    best.cov = NULL
    for( covariate in xcol[!xcol %in% best.covariate])
    {
      x.var = c(best.covariate, covariate)
      rmse = KNN.Internal1(data, x.var, ycol, kfold)

      if(rmse < best.rmse){

        best.rmse = rmse
        best.cov = x.var

      }

    }

    best.covariate = best.cov

    dataframe[itr,1] = paste(best.covariate, collapse="+")
    dataframe[itr,2] = best.rmse

  }
  return(dataframe)
}

# Computes best k using generalized cross validation
computeBestK = function(dataX, dataY, rangeK ){

  bestK = NULL
  dataX = as.matrix(dataX)
  maxK = max(rangeK)
  nnIdx = knnx.index(dataX, query = dataX, k = maxK)
  gcv = rep(0,length(rangeK))
  for (i in 1:length(rangeK)){
    predY = rowMeans(matrix(dataY[nnIdx[,1:rangeK[i]]], ncol = ncol(nnIdx[,1:rangeK[i]])))
    gcv[i] = mean(((dataY-predY)/(1-(1/rangeK[i])))^2)
  }
  bestK = rangeK[which.min(gcv)]
  return(bestK)
}

