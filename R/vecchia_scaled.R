######   scaled Vecchia approximation (Katzfuss, Guinness, Lawrence)  #######


### necessary packages

# install.packages('GpGp')  # need version >= 0.2.2
library(GpGp)
# install.packages('GPvecchia')
library(GPvecchia)



#' fit parameters using scaled Vecchia, assuming matern covariance
#'
#' @param y data vector of length n
#' @param inputs nxd matrix of input coordinates
#' @param ms vector of conditioning-set sizes
#' @param trend options are 'pre' (subtract sample mean as a preprocessing step),
#' 'zero' (no trend), 'intercept', 'linear' (incl intercept)
#' @param X nxp trend matrix (use if more complicated trend is desired)
#' @param nu smoothness parameter. 1.5,2.5,3.5,4.5 avoid bessel (faster). 
#' estimated if nu=NULL.
#' @param nug nugget or noise variance. estimated if nug=NULL.
#' @param scale scaling of inputs for ordering and conditioning. 
#' 'parms': by parameter estimates. 'ranges': to [0,1]. 'none': no scaling
#' @param var.ini initial value for GP variance parameter
#' @param ranges.ini initial values for d range parameters
#' @param select un-select input variables if estimated range parameter is
#' above select (assuming standardized [0,1] inputs)
#' @param print.level 0: no printing. 1: print outer loop. 2: print outer+inner loop
#' @param max.it maximum number of iterations for inner loop
#' @param tol.dec converged if dot product between the step and the gradient is 
#' less than \code{10^(-convtol)}
#' @param n.est subsample size for estimation
#' @param find.vcf find a variance correction factor to be used in prediction?
#' @param vcf.scorefun scoring function to be used for \code{find_vcf()}
#'
#' @return Object containing fit information, including for use in predictions_scaled()
#' @examples
#' inputs=matrix(runif(40),ncol=2)
#' y=sin(rowSums(inputs*5))
#' fit=fit_scaled(y,inputs)
#' summary.GpGp_fit(fit)
#' @export


#####   fitting function   ########

fit_scaled=function(y,inputs,ms=c(30),trend='pre',X,nu=3.5,nug=0,scale='parms',
              var.ini,ranges.ini,select=Inf,print.level=0,max.it=32,tol.dec=4,
              n.est=min(5e3,nrow(inputs)),find.vcf=TRUE,vcf.scorefun=ls) {
  
  ## dimensions
  n=nrow(inputs)
  d=ncol(inputs)

  ## specify trend covariates
  if(missing(X)) {
    if(trend=='zero'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
    } else if(trend=='intercept'){
      X=as.matrix(rep(1,n))
    } else if(trend=='linear'){
      X=cbind(rep(1,n),inputs)
    } else if(trend=='pre'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
      beta=mean(y)
      y=y-beta
    } else stop('invalid trend option specified')
  } else trend='X'
  
  ## default variance parameter
  if(missing(var.ini)) {
    cur.var=summary(stats::lm(y~X-1))$sigma^2
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
  
  ## only use subsample for estimation?
  if(n.est<n){
    ind.est=sample(1:n,n.est)
    y.full=y; inputs.full=inputs; X.full=X
    y=y[ind.est]; inputs=inputs[ind.est,,drop=FALSE]; X=X[ind.est,,drop=FALSE]
  }
  
  ## decrease or remove m values larger than n
  ms=unique(ifelse(ms<length(y),ms,length(y)-1))
  
  
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
      active=(cur.ranges<input.ranges*select)
      if(sum(active,na.rm=TRUE)==0) stop('all inputs inactive. increase select?')
      cur.ranges[!active]=Inf
      
      ## specify how to scale input dimensions
      cur.ranges[!active]=Inf
      
      ## order and condition based on current params
      if(scale=='parms'){ scales=1/cur.ranges
      } else if(scale=='ranges'){ scales=1/input.ranges
      } else if(scale=='none'){ scales=1
      } else stop(paste0('invalid argument scale=',scale))
      
      ## order and condition based on current params
      ord=GPvecchia::order_maxmin_exact(t(t(inputs)*scales))
      inputs.ord=inputs[ord,,drop=FALSE]
      y.ord=y[ord]
      X.ord=X[ord,,drop=FALSE]
      NNarray=GpGp::find_ordered_nn(t(t(inputs.ord)*scales),m)
      
      ## starting and fixed parameters
      cur.parms=c(cur.var,cur.ranges[active],cur.oth)
      fixed=NULL      
      if(fix.nu) fixed=c(fixed,length(cur.parms)-1)
      if(fix.nug) fixed=c(fixed,length(cur.parms))
      
      ## fisher scoring
      fit=GpGp::fit_model(y.ord,inputs.ord[,active,drop=FALSE],X.ord,
                     NNarray=NNarray,m_seq=m,convtol=tol,
                     start_parms=cur.parms,max_iter=maxit,
                     covfun_name=covfun,silent=(print.level<2),
                     reorder=FALSE,fixed_parms=fixed)
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
  if(n.est<n){
    fit$y=y.full
    fit$locs=inputs.full
    fit$X=X.full
  } else {
    fit$locs=inputs.ord
  }
  if(trend=='zero') {
    fit$X=as.matrix(rep(0,n))
  } else if(trend=='pre') {
      fit$betahat=beta
      fit$y=fit$y+beta
      fit$trend='intercept'
      fit$X=as.matrix(rep(1,n))
  }
  
  ### find variance correction factor, if requested
  if(find.vcf){
    fit$vcf=fit_vcf(fit,scale=scale,scorefun=vcf.scorefun)
  } else fit$vcf=1
  
  return(fit)
  
}





#######   prediction   ########

#' prediction using scaled Vecchia, using output from fit_scaled()
#'
#' @param fit object returned from fit_scaled()
#' @param locs_pred n.p x d matrix of test/prediction inputs/locations
#' @param m conditioning-set size (larger is more accurate but slower)
#' @param joint Joint predictions (does not return variances) or 
#' separate/independent predictions (does not produce joint samples)
#' @param n.sims desired number of samples from predictive distributions. 
#' if \code{n.sims=0}, posterior mean is returned.
#' @param predvar return prediction variances? (only if \code{joint=FALSE})
#' @param X_pred n.p x p trend matrix at locs_pred 
#' (if missing, will be generated based on fit object)
#' @param scale scaling of inputs for ordering and conditioning. 
#' 'parms': by parameter estimates. 'ranges': to [0,1]. 'none': no scaling
#'
#' @return Vector of length n.p (\code{n.sims=0}, \code{predvar=FALSE}) or 
#' list with entries \code{means} and/or \code{vars} and/or \code{samples}
#' @examples
#' inputs=matrix(runif(200),ncol=2)
#' y=sin(rowSums(inputs*5))
#' inputs.test=matrix(runif(100),ncol=2)
#' fit=fit_scaled(y,inputs)
#' preds=predictions_scaled(fit,inputs.test)
#' plot(rowSums(inputs.test),preds)
#' @export

predictions_scaled <- function(fit,locs_pred,m=100,joint=TRUE,nsims=0,
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
  # if(covparms[length(covparms)]==0) 
  #   covparms[length(covparms)]=covparms[1]*1e-12
  
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



##########   line search for variance correction factor   ###
fit_vcf=function(fit,m.pred=140,n.test=min(1e3,round(nrow(fit$locs)/5)),
                 scale='parms',scorefun=ls){
  
  # remove test data from fit object
  fitsearch=fit
  inds.test=sample(1:nrow(fit$locs),n.test)
  fitsearch$y=fit$y[-inds.test]
  fitsearch$locs=fit$locs[-inds.test,,drop=FALSE]
  fitsearch$X=fit$X[-inds.test,,drop=FALSE]
  
  # make predictions
  preds=predictions_scaled(fitsearch,locs_pred=fit$locs[inds.test,,drop=FALSE],
                           m=m.pred,joint=FALSE,predvar=TRUE,scale=scale,
                           X_pred=fit$X[inds.test,,drop=FALSE])
  
  # optimize correction factor
  y.test=fit$y[inds.test]
  objfun=function(vcf) scorefun(y.test,preds$means,preds$vars*vcf)
  vcf=optimize(objfun,c(1e-6,1e6))$minimum
  
  return(vcf)
  
}


### log score
ls=function(dat,mu,sig2) -mean(dnorm(dat,mu,sqrt(sig2),log=TRUE))
