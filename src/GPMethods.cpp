/*
 // MIT License
 // 
 //   Copyright (c) 2020 Abhinav Prakash, Rui Tuo, and Yu Ding
 // 
 //   Permission is hereby granted, free of charge, to any person obtaining a copy
 //   of this software and associated documentation files (the "Software"), to deal
 //   in the Software without restriction, including without limitation the rights
 //   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 //   copies of the Software, and to permit persons to whom the Software is
 //   furnished to do so, subject to the following conditions:
 // 
 //     The above copyright notice and this permission notice shall be included in all
 //     copies or substantial portions of the Software.
 // 
 //   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 //     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 //     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 //     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 //     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 //     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 //     SOFTWARE.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat OuterDiff(const arma::vec & x1,const arma::vec & x2){
  arma::mat od = arma::zeros(x1.n_elem,x2.n_elem);
  arma::rowvec x = x2.t();
  od.each_col() = x1;
  od.each_row() -= x;
  return od;
} 

arma::mat computeCorrelMat_(const arma::mat & X1,const arma::mat & X2,const arma::vec & theta){
  arma::mat correlMat = arma::zeros( X1.n_rows, X2.n_rows);
  for (arma::uword i=0; i < theta.n_elem; i++){
    correlMat = correlMat + pow(OuterDiff(X1.col(i),X2.col(i))/theta(i),2);
  }
  correlMat = exp(-0.5*correlMat);
  return correlMat;
} 

// [[Rcpp::export]]
arma::vec computeWeightedY(const arma::mat& X, const arma::vec& y, List params){
  arma::vec theta = params["theta"];
  double sigma_f = params["sigma_f"];
  double sigma_n = params["sigma_n"];
  double beta = params["beta"];
  arma::mat trainMat = pow(sigma_f,2)*computeCorrelMat_(X,X, theta);
  trainMat.diag() +=  pow(sigma_n,2);
  arma::mat upperCholTrainMat = arma::chol(trainMat);
  trainMat = arma::mat();
  arma::vec y_dash = y - beta;
  arma::vec weightedY = arma::solve(arma::trimatu(upperCholTrainMat),arma::solve(arma::trimatl(upperCholTrainMat.t()),y_dash));
  return weightedY;
}

// [[Rcpp::export]]
arma::vec predictGP(const arma::mat& X, const arma::vec& weightedY, const arma::mat& Xnew, List params){
  arma::vec pred = arma::zeros(Xnew.n_rows,1);
  arma::vec theta = params["theta"];
  double sigma_f = params["sigma_f"];
  double beta = params["beta"];
  arma::mat testCovMat = pow(sigma_f,2)*computeCorrelMat_(Xnew,X, theta);
  pred = beta + (testCovMat*weightedY);
  return pred;
}

// [[Rcpp::export]]
double computeLogLikGP_(const arma::mat & X, const arma::vec&  y, const List& params){
  arma::vec theta = params["theta"];
  double sigma_f = params["sigma_f"];
  double sigma_n = params["sigma_n"];
  double beta = params["beta"];
  arma::mat CovMat = pow(sigma_f,2)*computeCorrelMat_(X,X, theta);
  CovMat.diag() +=  pow(sigma_n,2);
  arma::mat UpperCholMat = chol(CovMat);
  arma::vec y_dash = y - beta;
  double t1 = 0.5*as_scalar((y_dash.t()*arma::solve(arma::trimatu(UpperCholMat),arma::solve(arma::trimatl(UpperCholMat.t()),y_dash))));
  double t2 = arma::sum(log(arma::abs(UpperCholMat.diag())));
  double t3 = log(2*arma::datum::pi)*UpperCholMat.n_rows/2;
  return (t1+t2+t3);
}

// [[Rcpp::export]]
arma::vec computeLogLikGradGP_(const arma::mat & X,const arma::vec & y, List params){
  arma::vec theta = params["theta"];
  int n_theta = theta.n_elem;
  double sigma_f = params["sigma_f"];
  double sigma_n = params["sigma_n"];
  double beta = params["beta"];
  arma::mat CorrelMat = computeCorrelMat_(X,X, theta);
  arma::mat CovMat = pow(sigma_f,2)*CorrelMat;
  CovMat.diag() +=  pow(sigma_n,2);
  arma::mat InvMat = arma::inv_sympd(CovMat);
  CovMat = arma::mat();
  arma::vec gradval(n_theta + 3);
  gradval.zeros();
  arma::vec y_dash = y - beta;
  arma::vec alpha = InvMat*y_dash;
  arma::mat diffMat = (alpha*alpha.t()) - InvMat;
  arma::vec OneVec = arma::ones<arma::vec>(y.n_elem);
  arma::vec SolOneVec = InvMat*OneVec;
  InvMat = arma::mat();
  for (int i=0; i<n_theta; i++){
    arma::mat delThetaMat = (pow(sigma_f,2)*pow(OuterDiff(X.col(i),X.col(i)),2)/pow(theta(i),3))%CorrelMat;
    delThetaMat = diffMat*delThetaMat;
    gradval(i) = -0.5*arma::sum(delThetaMat.diag());
  }
  arma::mat delSigma_fMat = 2*sigma_f*CorrelMat;
  delSigma_fMat = diffMat*delSigma_fMat;
  gradval(n_theta) = -0.5*arma::sum(delSigma_fMat.diag());
  delSigma_fMat = arma::mat();
  arma::mat delSigma_nMat = 2*sigma_n*diffMat;
  gradval(n_theta+1) = -0.5*sum(delSigma_nMat.diag());
  delSigma_nMat = arma::mat();
  gradval(n_theta+2) = 0.5*arma::as_scalar((2*beta*OneVec.t()*SolOneVec)-(y.t()*SolOneVec)-(OneVec.t()*(alpha+(beta*SolOneVec))));
  return gradval;
}

// [[Rcpp::export]]
arma::vec computeLogLikGradGPZeroMean_(const arma::mat & X,const arma::vec & y, List params){
  arma::vec theta = params["theta"];
  int n_theta = theta.n_elem;
  double sigma_f = params["sigma_f"];
  double sigma_n = params["sigma_n"];
  double beta = params["beta"];
  arma::mat CorrelMat = computeCorrelMat_(X,X, theta);
  arma::mat CovMat = pow(sigma_f,2)*CorrelMat;
  CovMat.diag() +=  pow(sigma_n,2);
  arma::mat InvMat = arma::inv_sympd(CovMat);
  CovMat = arma::mat();
  arma::vec gradval(n_theta + 2);
  gradval.zeros();
  arma::vec y_dash = y - beta;
  arma::vec alpha = InvMat*y_dash;
  arma::mat diffMat = (alpha*alpha.t()) - InvMat;
  arma::vec OneVec = arma::ones<arma::vec>(y.n_elem);
  arma::vec SolOneVec = InvMat*OneVec;
  InvMat = arma::mat();
  for (int i=0; i<n_theta; i++){
    arma::mat delThetaMat = (pow(sigma_f,2)*pow(OuterDiff(X.col(i),X.col(i)),2)/pow(theta(i),3))%CorrelMat;
    delThetaMat = diffMat*delThetaMat;
    gradval(i) = -0.5*arma::sum(delThetaMat.diag());
  }
  arma::mat delSigma_fMat = 2*sigma_f*CorrelMat;
  delSigma_fMat = diffMat*delSigma_fMat;
  gradval(n_theta) = -0.5*arma::sum(delSigma_fMat.diag());
  delSigma_fMat = arma::mat();
  arma::mat delSigma_nMat = 2*sigma_n*diffMat;
  gradval(n_theta+1) = -0.5*sum(delSigma_nMat.diag());
  //delSigma_nMat = arma::mat();
  //gradval(n_theta+2) = 0.5*arma::as_scalar((2*beta*OneVec.t()*SolOneVec)-(y.t()*SolOneVec)-(OneVec.t()*(alpha+(beta*SolOneVec))));
  return gradval;
}

// [[Rcpp::export]]
List computeDiffCov_(const arma::mat& X1, const arma::vec y1, const arma::mat& X2, const arma::vec y2, const arma::mat& XT, arma::vec theta, double sigma_f, double sigma_n, double beta ){
  arma::mat KX1X1 = pow(sigma_f,2)*computeCorrelMat_(X1,X1,theta);
  KX1X1.diag() +=  pow(sigma_n,2);
  arma::mat invKX1X1 = arma::inv_sympd(KX1X1);
  KX1X1.reset();
  arma::mat KXTX1 = pow(sigma_f,2)*computeCorrelMat_(XT,X1,theta);
  arma::vec mu1 = beta + (KXTX1*invKX1X1*(y1-beta)) ;
  arma::mat K1 = invKX1X1*KXTX1.t() ;
  arma::mat K = KXTX1*K1;
  KXTX1.reset(); 
  invKX1X1.reset();
  arma::mat KX2X2 = pow(sigma_f,2)*computeCorrelMat_(X2,X2,theta);
  KX2X2.diag() +=  pow(sigma_n,2);
  arma::mat invKX2X2 = arma::inv_sympd(KX2X2);
  KX2X2.reset();
  arma::mat KXTX2 = pow(sigma_f,2)*computeCorrelMat_(XT,X2,theta);
  arma::vec mu2 = beta + (KXTX2*invKX2X2*(y2-beta)) ;
  arma::mat K2 = KXTX2*invKX2X2 ;
  invKX2X2.reset();
  K = K + (K2*KXTX2.t());
  KXTX2.reset(); 
  arma::mat KX2X1 = pow(sigma_f,2)*computeCorrelMat_(X2,X1,theta);
  K = K - (2*K2*KX2X1*K1) ;
  K = (K + K.t())/2 ;
  List retList = List::create(_["diffCovMat"]=K, _["mu1"] = mu1, _["mu2"] = mu2);
  return retList;
}

// [[Rcpp::export]]
arma::vec computeConfBand_(const arma::mat& diffCovMat, double confLevel){
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval,eigvec,diffCovMat);
  arma::uword firstIdx = arma::as_scalar(arma::find(eigval < 1E-8,1,"last"));
  arma::mat lambda = arma::diagmat(arma::reverse(eigval.subvec(firstIdx,eigval.n_elem-1)));
  arma::mat eigmat = arma::fliplr(eigvec.cols(firstIdx,eigval.n_elem-1));
  arma::uword n_eig = lambda.n_rows;
  double radius = sqrt(R::qchisq(confLevel,n_eig,1,0));
  arma::uword nSamples = 1000;
  //arma::arma_rng::set_seed(1);
  arma::mat Z = arma::zeros(n_eig, nSamples);
  arma::uword n = 0;
  while (n < nSamples){
    arma::vec zSample = rnorm(n_eig);
    double zSum = sqrt(arma::accu(pow(zSample,2)));
    if (zSum <= radius){
      Z.col(n) = zSample ;
      n += 1;
    }
  }
  arma::mat G = eigmat*arma::sqrt(lambda)*Z;
  G = arma::abs(G);
  arma::vec band = max(G,1);
  return band;
}

