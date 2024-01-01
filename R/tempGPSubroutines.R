# MIT License
# 
# Copyright (c) 2020-2022 Abhinav Prakash and Yu Ding
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' @useDynLib DSWE, .registration = TRUE

computeThinningNumber = function(trainX, max_thinning_number){
  thinning_vec = rep(max_thinning_number, ncol(trainX))
  for (col_idx in c(1:length(thinning_vec))){
    col_thinning_vec = which(c(1,abs(stats::pacf(trainX[,col_idx], plot = FALSE, lag.max = max_thinning_number)$acf[,1,1])) <= (2/sqrt(nrow(trainX))))
    if (length(col_thinning_vec) != 0){
      thinning_vec[col_idx] = min(col_thinning_vec)
    }
  }
  thinningNumber = max(thinning_vec)
  return(thinningNumber)
}

createThinnedBins = function(dataX, dataY, thinningNumber){
  nData = nrow(dataX)
  thinnedBins = list()
  if (thinningNumber < 2) {
    thinnedBins[[1]] = list(X = dataX, y = dataY)
  } else {
    for (i in 1:thinningNumber){
      nPoints = floor((nData - i)/thinningNumber)
      lastIdx = i + (nPoints*thinningNumber)
      idx = seq(i, lastIdx, length.out = (nPoints+1))
      thinnedBins[[i]] = list(X = dataX[idx,,drop = FALSE], y = dataY[idx])
    }
  }
  
  return(thinnedBins)
}

estimateBinnedParams= function(databins, fast_computation, optim_control){
  nCov = ncol(databins[[1]]$X)
  theta = rep(0,nCov)
  for (i in 1:nCov){
    theta[i] = mean(unlist(lapply(databins, function(bin) stats::sd(bin$X[,i]))))
  }
  beta = mean(unlist(lapply(databins, function(bin) mean(bin$y))))
  sigma_f = mean(unlist(lapply(databins, function(bin) stats::sd(bin$y)/sqrt(2))))
  sigma_n = mean(unlist(lapply(databins, function(bin) stats::sd(bin$y)/sqrt(2))))
  parInit = c(theta,sigma_f,sigma_n,beta)
  objFun = function(par){computeloglikSumTempGP(databins,
                                                params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  objGrad = function(par){computeloglikGradSumTempGP(databins,
                                                     params = list(theta=par[1:nCov],sigma_f=par[nCov+1],sigma_n=par[nCov+2],beta=par[nCov+3]))}
  if (fast_computation){
    optimResult = adam_optimizer(databins, parInit, 
                                 optim_control$batch_size, 
                                 optim_control$learn_rate,
                                 optim_control$max_iter,
                                 optim_control$tol, 
                                 optim_control$beta1, 
                                 optim_control$beta2, 
                                 optim_control$epsilon, 
                                 optim_control$logfile )
    estimatedParams = list(theta = abs(optimResult$par[1:nCov]), sigma_f = abs(optimResult$par[nCov+1]), sigma_n = abs(optimResult$par[nCov+2]), beta = optimResult$par[nCov+3])
    objVal = NULL
    gradVal = NULL
    return(list(estimatedParams = estimatedParams,objVal = objVal, gradVal = gradVal))
  } else {
    optimResult = stats::nlminb(start = parInit, objective = objFun, gradient = objGrad, control = list(eval.max = 1000,iter.max = 1000, trace = 0,step.min = 1e-6))
    estimatedParams = list(theta = abs(optimResult$par[1:nCov]), sigma_f = abs(optimResult$par[nCov+1]), sigma_n = abs(optimResult$par[nCov+2]), beta = optimResult$par[nCov+3])
    objVal = -optimResult$objective
    gradVal = -objGrad(optimResult$par)  
    return(list(estimatedParams = estimatedParams,objVal = objVal, gradVal = gradVal))
  }
}

###
computeloglikSumTempGP = function(databins,params){
  logliksum = sum(unlist(lapply(databins, function(bin) computeLogLikGP_(bin$X,bin$y,params))))
  return(logliksum)
}


###
computeloglikGradSumTempGP = function(databins,params){
  
  loglikGradSum = Reduce("+",lapply(databins, function(bin) computeLogLikGradGP_(bin$X,bin$y,params)))
  return(loglikGradSum)
}

estimateLocalFunctionParams= function(trainT, residual){
  theta = sd(trainT)
  sigma_f = sd(residual)/sqrt(2)
  sigma_n = sigma_f
  parInit = c(theta,sigma_f,sigma_n)
  objFun = function(par){computeLogLikGP_(as.matrix(trainT), residual,
                                          params = list(theta=par[1],sigma_f=par[2],sigma_n=par[3],beta=0))}
  objGrad = function(par){computeLogLikGradGPZeroMean_(as.matrix(trainT), residual,
                                                       params = list(theta=par[1],sigma_f=par[2],sigma_n=par[3],beta=0))}
  optimResult = stats::nlminb(start = parInit, objective = objFun, gradient = objGrad )
  estimatedParams = list(theta = abs(optimResult$par[1]), sigma_f = abs(optimResult$par[2]), sigma_n = abs(optimResult$par[3]), beta = 0)
  objVal = optimResult$objective
  gradVal = objGrad(optimResult$par)
  return(list(estimatedParams = estimatedParams,objVal = objVal, gradVal = gradVal))
}


computeLocalFunction = function(residual, traindataT, testdataT, neighbourhood){
  pred = rep(0,length(testdataT))
  for (i in 1:length(testdataT)){
    distance = abs(testdataT[i] - traindataT)
    trainIdx = which(distance < neighbourhood)
    if (length(trainIdx)>0){
      if (var(residual[trainIdx]) < .Machine$double.eps || is.na(var(residual[trainIdx]))){
        warning("While computing g(t), variance of the training residuals is numerically zero for time index: ",testdataT[i], '\nUsing mean of the response as the prediction.\n')
        pred[i] = mean(residual[trainIdx])
      } else {
        params = try (
          expr = estimateLocalFunctionParams(traindataT[trainIdx],residual[trainIdx]), silent = TRUE
          )
        if (inherits(params, "try-error")){
          warning("Computing g(t) is numerically unstable for time index: ",testdataT[i], '\nUsing mean of the response as the prediction.\n')
          pred[i] =  mean(residual[trainIdx])
        } else {
          weightedRes = computeWeightedY(as.matrix(traindataT[trainIdx]),residual[trainIdx],params$estimatedParams)
          pred[i] = predictGP(as.matrix(traindataT[trainIdx]),weightedRes,as.matrix(testdataT[i]),params$estimatedParams)    
        }
      }
      
    } else {
      pred[i] = 0
    }
    
  }
  return(pred)
}

adam_optimizer = function(data_bins, par_init, batch_size, learn_rate, max_iter, tol, beta1, beta2, epsilon, logfile){
  t = 0
  nCov = ncol(data_bins[[1]]$X)
  params_t = par_init
  m_t = rep(0, length(par_init))
  v_t = rep(0, length(par_init))
  
  params_mat = matrix(nrow = max_iter, ncol = length(params_t))  
  
  
  while (TRUE){
    t = t+1
    sampled_bin = sample(length(data_bins), 1)
    if (batch_size < length(data_bins[[sampled_bin]]$y)){
      sampled_idx = sample(length(data_bins[[sampled_bin]]$y), batch_size)
    } else {
      sampled_idx = 1:length(data_bins[[sampled_bin]]$y)
    }
    sample_x = data_bins[[sampled_bin]]$X[sampled_idx, , drop=FALSE]
    sample_y = data_bins[[sampled_bin]]$y[sampled_idx]
    
    grad_t = computeLogLikGradGP_(sample_x, sample_y, 
                                  params = list(theta=params_t[1:nCov],sigma_f=params_t[nCov+1],sigma_n=params_t[nCov+2],beta=params_t[nCov+3]))
    m_t = (beta1 * m_t) + ((1 - beta1) * grad_t)
    v_t = (beta2 * v_t) + ((1 - beta2) * (grad_t**2))
    
    m_hat = m_t/(1 - (beta1**t))
    v_hat = v_t/(1 - (beta2**t))
    
    params_prev = params_t
    params_t = params_t - (learn_rate * (m_hat / (sqrt(v_hat) + epsilon)))
    params_mat[t, ] = params_t
    
    if (max(abs(params_t - params_prev)) < tol || t == max_iter){
      if (!is.null(logfile)){
        utils::write.table(params_mat, file = logfile, sep = ",", row.names = FALSE, col.names = FALSE)
      }
      return(list(par=params_t))
    }
  }
}


fit_scaled_thinned <- function(y,inputs,thinnedBins,T,ms=c(10),trend='zero',X,nu=1.5,nug=NULL,scale='none',
                               var.ini,ranges.ini,select=Inf,print.level=0,max.it=1000,tol.dec=5,
                               find.vcf=FALSE,vcf.scorefun=ls,grouped=TRUE) {
  source('https://raw.githubusercontent.com/katzfuss-group/scaledVecchia/master/vecchia_scaled.R')
  
  
  
  loc <- list()
  Y <- list()
  X <- list()
  for (i in 1:T) {
    loc[[i]] <- thinnedBins[[i]]$X
    Y[[i]] <- thinnedBins[[i]]$y
  }
  
  
  n = length(Y[[1]])
  d = ncol(loc[[1]])
  if(trend=='zero'){
    for(i in 1:T){
      X[[i]]=as.matrix(sample(c(-1,1),n,replace=TRUE))
    }
  }else if(trend=='pre'){
    for(i in 1:T){
      X[[i]]=as.matrix(sample(c(-1,1),n,replace=TRUE))
    }
    beta <- mean(y)
    for (i in 1:T) {
      Y[[i]]=Y[[i]]-beta
    }
  }
  
  
  ## default variance parameter
  if(missing(var.ini)) {
    cur.var=summary(stats::lm(Y[[1]]~X[[1]]-1))$sigma^2
  } else cur.var=var.ini
  
  ## default range parameters
  input.ranges=apply(inputs,2,function(x) diff(range(x)))
  if(missing(ranges.ini)) cur.ranges=.2*input.ranges else cur.ranges=ranges.ini
  active=rep(TRUE,d)
  
  ## fixed nugget?
  if(is.null(nug)){
    fix.nug=FALSE; nug=.01*var(y)
  } else fix.nug=TRUE
  
  ## smoothness: fixed? bessel?
  if(is.null(nu)){
    covfun='matern_scaledim'
    cur.oth=c(3.5,nug)
    fix.nu=FALSE
  } else if(nu %in% (.5+(1:4))){
    covfun=paste0("matern",nu*10,"_scaledim")
    cur.oth=nug
    fix.nu=FALSE
  } else {
    covfun='matern_scaledim'
    cur.oth=c(nu,nug)
    fix.nu=TRUE    
  }
  
  
  
  ### for increasing m
  for(i.m in 1:length(ms)){
    
    m=ms[i.m]
    
    if(i.m<length(ms)){ tol=10^(-tol.dec-2) } else {tol=10^(-tol.dec)}
    
    ### increase maxit until convergence
    conv=FALSE
    maxit=2
    while(conv==FALSE & maxit<=max.it){
      
      if(print.level>0) {
        print(paste0('m=',m,', maxit=',maxit)); print(cur.ranges)}
      
      ## check for inactive input dims (large range params)
      #active=(cur.ranges<input.ranges*select)
      if(sum(active,na.rm=TRUE)==0) stop('all loc inactive. increase select?')
      cur.ranges[!active]=Inf
      
      ## specify how to scale input dimensions
      cur.ranges[!active]=Inf
      
      ## order and condition based on current params
      if(scale=='parms'){ scales=1/cur.ranges
      } else if(scale=='ranges'){ scales=1/input.ranges
      } else if(scale=='none'){ scales=1
      } else stop(paste0('invalid argument scale=',scale))
      NNlist<- list()
      NNarray=list()
      loc.ord=list()
      Y.ord=list()
      X.ord=list()
      for(i in 1:T){
        ## order and condition based on current params
        ord=GPvecchia::order_maxmin_exact(t(t(loc[[i]])*scales))
        loc.ord[[i]]=loc[[i]][ord,,drop=FALSE]
        Y.ord[[i]]=Y[[i]][ord]
        X.ord[[i]]=as.matrix(X[[i]][ord,,drop=FALSE])
        output=as.matrix(loc.ord[[i]])
        NNarray[[i]]=GpGp::find_ordered_nn(t(t(as.matrix(output))*scales),m)
      }
      
      ## starting and fixed parameters
      cur.parms=c(cur.var,cur.ranges[active],cur.oth)
      fixed=NULL      
      if(fix.nu) fixed=c(fixed,length(cur.parms)-1)
      if(fix.nug) fixed=c(fixed,length(cur.parms))
      
      
      
      # get link functions
      linkfuns <- GpGp::get_linkfun(covfun)
      link <- linkfuns$link
      dlink <- linkfuns$dlink
      invlink <- linkfuns$invlink
      
      start_parms=cur.parms
      
      invlink_startparms <- invlink(start_parms)
      lonlat <- linkfuns$lonlat
      if(lonlat){
        cat("Assuming columns 1 and 2 of locs are (longitude,latidue) in degrees\n")
      }
      
      penalty=list()
      pen=list()
      dpen=list()
      ddpen=list()
      
      for (i in 1:T){
        NNlist[[i]] <- GpGp::group_obs((NNarray[[i]][,1:(m+1)]))
        penalty[[i]] <- GpGp::get_penalty(Y.ord[[i]],X.ord[[i]],loc.ord[[i]],covfun) # this function takes in y as a vector, returns a list of three functions, # NEED CHANGE
        pen[[i]] <- penalty[[i]]$pen
        dpen[[i]] <- penalty[[i]]$dpen
        ddpen[[i]] <- penalty[[i]]$ddpen
      }
      
      
      
      
      likfun <- function(logparms){
        
        lp <- rep(NA,length(start_parms))
        lp[active] <- logparms
        lp[!active] <- invlink_startparms[!active]
        
        sumLoglik <- 0 # Add
        sumGrad <- 0 # Add
        sumInfo <- 0 # Add
        for (i in 1:T){
          
          if(grouped){
            likobj <-  GpGp::vecchia_grouped_profbeta_loglik_grad_info(link(lp),covfun,y=as.matrix(Y.ord[[i]]),X=as.matrix(X.ord[[i]]),locs=as.matrix(loc.ord[[i]]),(NNlist[[i]]))
          }else{
            likobj <-  GpGp::vecchia_profbeta_loglik_grad_info(link(lp),covfun,y=as.matrix(Y.ord[[i]]),X=as.matrix(X.ord[[i]]),locs=as.matrix(loc.ord[[i]]),(NNarray[[i]]))
            
          }
          likobj$loglik <- -likobj$loglik - pen[[i]](link(lp))
          likobj$grad <- -c(likobj$grad)*dlink(lp) -
            dpen[[i]](link(lp))*dlink(lp)
          likobj$info <- likobj$info*outer(dlink(lp),dlink(lp)) -
            ddpen[[i]](link(lp))*outer(dlink(lp),dlink(lp))
          likobj$grad <- likobj$grad[active]
          likobj$info <- likobj$info[active,active]
          
          sumLoglik <- sumLoglik + likobj$loglik 
          sumGrad <- sumGrad + likobj$grad
          sumInfo <- sumInfo + likobj$info
          
        }
        likobj$loglik <- sumLoglik/T
        likobj$grad <- sumGrad/T
        likobj$info <- sumInfo/T
        
        return(likobj)        
        
      }
      
      
      fit <- GpGp::fisher_scoring( likfun,invlink(start_parms)[active],
                                   link,silent=FALSE, convtol = tol, max_iter = maxit ) 
      invlink_startparms[active] <- fit$logparms
      #start_parms[active] <- fit$covparms
      start_parms <- link(invlink_startparms)
      fit$loglik <- -fit$loglik - pen[[i]](start_parms)
      invlink_startparms <- invlink(start_parms)
      
      
      # return fit and information used for predictions
      fit$covfun_name <- covfun
      #fit$covparms <- start_parms
      lp <- rep(NA,length(start_parms))
      lp[active] <- fit$logparms
      lp[!active] <- invlink_startparms[!active]
      fit$covparms <- link(lp)
      
      
      fit$locs <- inputs
      fit$y <- y
      
      class(fit) <- "GpGp_fit"
      cur.var=fit$covparms[1]
      cur.ranges[active]=fit$covparms[1+(1:sum(active))]
      cur.oth=fit$covparms[-(1:(1+sum(active)))]
      conv=fit$conv
      maxit=maxit*2
      
    }
  }
  
  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var,cur.ranges,cur.oth)
  fit$trend=trend
  
  if(trend=='zero') {
    fit$X=as.matrix(rep(0,nrow(fit$locs)))
  } else if(trend=='pre') {
    fit$betahat=beta
    fit$y=fit$y+beta
    fit$trend='intercept'
    fit$X=as.matrix(rep(1,nrow(fit$locs)))
  }
  
  ### find variance correction factor, if requested
  if(find.vcf){
    fit$vcf=GpGp::fit_vcf(fit,scale=scale,scorefun=vcf.scorefun)
  } else fit$vcf=1
  
  return(fit)
  
}
predictions_scaled_thinned <- function(fit,locs_pred,m=400,joint=TRUE,nsims=0,
                                       predvar=FALSE,X_pred,scale='parms'){
  
  y_obs = fit$y
  locs_obs = fit$locs
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  n_obs <- nrow(locs_obs)
  n_pred <- nrow(locs_pred)
  if(is.null(fit$vcf)) vcf=1 else vcf=fit$vcf
  
  # ## add nugget for numerical stability
  if(covparms[length(covparms)]==0) 
    covparms[length(covparms)]=covparms[1]*1e-12
  
  # specify trend if missing
  if(missing(X_pred)){
    if(fit$trend=='zero'){
      X_pred=as.matrix(rep(0,n_pred))
    } else if(fit$trend=='intercept'){
      X_pred=as.matrix(rep(1,n_pred))
    } else if(fit$trend=='linear'){
      X_pred=cbind(rep(1,n_pred),locs_pred)
    } else stop('X_pred must be specified')
  }
  
  # specify how to scale input dimensions
  if(scale=='parms'){ scales=1/covparms[1+(1:ncol(locs_obs))]
  } else if(scale=='ranges'){ scales=1/apply(locs_obs,2,function(x) diff(range(x)))
  } else stop(paste0('invalid argument scale=',scale))
  
  
  ### 
  if(joint){  # joint predictions
    
    # get orderings
    temp=order_maxmin_pred(t(t(locs_obs)*scales),t(t(locs_pred)*scales))
    ord1=temp$ord
    ord2=temp$ord_pred
    
    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    y  <- y_obs[ord1] - Xord_obs %*% beta
    
    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    
    # get nearest neighbor array
    sm = if (n_pred<1e5) 2 else 1.5
    NNarray_all <- find_ordered_nn_pred(t(t(locs_all)*scales),m,
                                        fix.first=n_obs,searchmult=sm)
    NNarray_pred=NNarray_all[-(1:n_obs),-1]
    
    means=numeric(length=n_pred)
    if(nsims>0) samples=array(dim=c(n_pred,nsims))
    
    # make predictions sequentially
    for(i in 1:n_pred){
      
      # NN conditioning sets
      NN=sort(NNarray_pred[i,])
      NN_obs=NN[NN<=n_obs]
      NN_pred=NN[NN>n_obs]-n_obs
      
      # (co-)variances
      K=get(covfun_name)(covparms,locs_all[c(NN,i+n_obs),])
      cl=t(chol(K))
      
      # prediction
      y.NN=y[NN_obs]
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,means[NN_pred]))
      if(nsims>0){ # conditional simulation
        pred.var=cl[m+1,m+1]^2*vcf
        for(s in 1:nsims){
          pm=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,samples[NN_pred,s]))
          samples[i,s]=stats::rnorm(1,pm,sqrt(pred.var))
        }
      }
    }
    # add (prior) mean and return to original ordering
    means[ord2] = means + c(Xord_pred %*% beta) 
    if(nsims==0){ 
      preds=means 
    } else {
      samples[ord2,] = samples + c(Xord_pred %*% beta)
      preds=list(means=means,samples=samples)
    }
    
  } else {  # separate predictions
    
    if(nsims>0) stop('cannot produce joint samples when joint=FALSE')
    
    y  = y_obs - X_obs %*% beta
    
    # find the NNs 
    m=min(m,nrow(locs_obs))
    
    NNarray=FNN::get.knnx(t(t(locs_obs)*scales),
                          
                          t(t(locs_pred)*scales),m)$nn.index
    
    
    means=vars=numeric(length=n_pred)
    for(i in 1:n_pred){
      
      # NN conditioning sets
      NN=NNarray[i,]
      
      # (co-)variances
      K=get(covfun_name)(covparms,rbind(locs_obs[NN,],locs_pred[i,]))
      cl=t(chol(K))
      
      # prediction
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],y[NN])
      vars[i]=cl[m+1,m+1]^2*vcf
      
    }
    means=means+c(X_pred %*% beta)
    
    if(predvar==FALSE){ 
      preds=means 
    } else {
      preds=list(means=means,vars=vars)
    }
  }
  return(preds)
}


#######   obs.pred maxmin ordering   ########
order_maxmin_pred<-function(locs, locs_pred,refine=FALSE){
  
  ord<-1:nrow(locs) #GPvecchia::order_maxmin_exact(locs)
  ord_pred <-GPvecchia::order_maxmin_exact(locs_pred)
  
  if(refine){
    
    locs_all = rbind(locs, locs_pred)
    
    n <- nrow(locs)
    m <- min( round(sqrt(n)), 200 )
    
    n_pred <- nrow(locs_pred)
    # next is to find 'ord_pred', a maxmin reordering of prediction locations
    NN <- FNN::get.knn( locs_all, k = m )$nn.index
    #NN_pred <- FNN::get.knnx( locs, locs_pred, k = 1 )$nn.dist
    # use ord, then order by NN_pred
    index_in_position <- c( ord, n + ord_pred, rep(NA,n_pred) )
    position_of_index <- order(index_in_position[1:(n+n_pred)])
    
    # move an index to the end if it is a
    # near neighbor of a previous location
    curlen <- n + n_pred
    nmoved <- 0
    for(j in (n+1):(n+2*n_pred) ){
      # nneigh tells us how many neighbors to look at
      # in order to decide whether the current point
      # has a previously ordered neighbor
      nneigh <- round( min(m,1*(n+n_pred)/(j-nmoved+1)) )
      neighbors <- NN[index_in_position[j],1:nneigh]
      if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
        nmoved <- nmoved+1
        curlen <- curlen + 1
        position_of_index[ index_in_position[j] ] <- curlen
        index_in_position[curlen] <- index_in_position[j]
        index_in_position[j] <- NA
      }
    }
    
    ord_pred <- index_in_position[ !is.na( index_in_position ) ][(n+1):(n+n_pred)] - n
    
  }
  
  return(list(ord=ord, ord_pred=ord_pred))
}




#######   find NN for prediction locations   ########
find_ordered_nn_pred <- function(locs,m,fix.first=0,searchmult=2){
  
  # if locs is a vector, convert to matrix
  if( is.null(ncol(locs)) ){
    locs <- as.matrix(locs)
  }
  
  # number of locations
  n <- nrow(locs)
  m <- min(m,n-1)
  mult <- 2
  
  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- min(apply( locs, 2, stats::sd ))
  locs <- locs + matrix( ee*1e-6*stats::rnorm(n*ncol(locs)), n, ncol(locs) )    
  
  # to store the nearest neighbor indices
  NNarray <- matrix(NA,n,m+1)
  
  # to the first mult*m+1 by brute force
  maxval <- min( mult*m + 1, n )
  if(fix.first<=maxval){
    NNarray[1:maxval,] <- GpGp::find_ordered_nn_brute(locs[1:maxval,,drop=FALSE],m)
  } else {
    maxval=fix.first
    NNarray[1:(m+1),] <- GpGp::find_ordered_nn_brute(locs[1:(m+1),,drop=FALSE],m)
    NNarray[1:maxval,1]=1:maxval
    NNarray[(m+1):maxval,1+(1:m)]=matrix(rep(1:m,maxval-m),byrow=TRUE,ncol=m)
  }
  query_inds <- min( maxval+1, n):n
  data_inds <- 1:n
  
  msearch <- m
  
  while( length(query_inds) > 0 ){
    msearch <- min( max(query_inds), round(searchmult*msearch) )
    data_inds <- 1:min( max(query_inds), n )
    NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
    less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
    sum_less_than_k <- apply(less_than_k,1,sum)
    ind_less_than_k <- which(sum_less_than_k >= m+1)
    
    NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))
    
    NNarray[ query_inds[ind_less_than_k], ] <- NN_m
    
    query_inds <- query_inds[-ind_less_than_k]
    
  }
  
  return(NNarray)
}


