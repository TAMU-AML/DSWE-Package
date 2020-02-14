KNN.Internal1 = function(data, xcol, ycol, kfold){

  set.seed(1)
  data = data[sample(1:nrow(data), nrow(data), replace = FALSE), ]
  folds = cut(1:nrow(data), breaks = kfold, labels = FALSE)
  rmse = 0


  for(i in 1:kfold)
  {
    rmse_k = rep(0, length(seq(5,50,5)))
    itr = 0
    train_data = data[folds != i, ]
    test_data = data[folds == i, ]
    folds2 = cut(1:nrow(train_data), breaks = kfold, labels = FALSE)
    for(j in seq(5,50,5))
    {
      rmse2 = 0
      for( h in 1:kfold)
      {
        test_data2 = train_data[folds2 == h, ]
        train_data2 = train_data[folds2 != h, ]
        test_pred2 = knn.reg(as.matrix(train_data2[, xcol]), as.matrix(test_data2[, xcol]), train_data2[,ycol], k = j)
        residual2 = sqrt((sum((test_data2[, ycol] - test_pred2$pred)^2)) / length(test_data2[, ycol]))
        rmse2 = rmse + residual2
      }
      itr = itr + 1
      rmse_k[itr] = rmse2
    }


    test_pred = knn.reg(as.matrix(train_data[, xcol]), as.matrix(test_data[, xcol]), train_data[, ycol], k = seq(5,50,5)[which.min(rmse_k)])
    residual = sqrt((sum((test_data[, ycol] - test_pred$pred)^2)) / length(test_data[, ycol]))
    rmse = rmse + residual

  }

  avg_rmse = round(rmse / kfold, 2)

  return(avg_rmse)

}


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
