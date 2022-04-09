#include <vector>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppNumerical.h>


using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]



NumericMatrix DFtoNM(DataFrame x) {
  int nRows=x.nrow();  
  NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
    y(_,i)=NumericVector(x[i]);
  }  
  return y;
}


// [[Rcpp::export]]
Rcpp::List SFEc(DataFrame x,DataFrame y,double lambda,NumericVector beta) {
  
  int n=x.nrow(), p=x.ncol();
  arma::rowvec sf(p);
  arma::rowvec s1(p);
  arma::rowvec s2(p);
  arma::rowvec mins(p);
  arma::rowvec  maxs(p);
  arma::rowvec  lams(p);
  
  arma::rowvec ck(p);
  arma::uvec indxset(p);
  //order y
  arma::vec ytime=y[0];
  arma::vec ystatus=y[1];
  arma::uvec timeinx=arma::sort_index(ytime);
  ystatus=ystatus(timeinx);
  ytime=ytime(timeinx);
  // order x
  
  arma::mat x1=as<arma::mat>(DFtoNM(x));
  arma::mat x2(n,p);
  for (int i=0; i<n;i++) {
    x2.row(i)=x1.row(timeinx(i));
  }  
  // status=1 rowinx1
  int event1num=0;
  
  // count event 1 obsverations 
  arma:: uvec event1indx= find(ystatus==1);
  arma:: uvec event2indx= find(ystatus==2);
  arma::mat event1x=x2.rows(event1indx);
  
  
  // sum by column
  ck=sum(event1x,0);
  
  arma:: vec geindx;
  arma:: vec leindx;
  arma::mat tempmat1, tempmat2;
  arma:: uvec lessindx;
  
  for (int i=0; i<event1x.n_rows;i++) {
    //greater than evetn1
    // not the last node
    if (event1indx(i)<n-1){
      tempmat1=x2.rows(event1indx(i)+1,n-1);
    }
    //less than event 1 and with event 2
    lessindx=find(event2indx<event1indx(i));
    tempmat2=x2.rows(event2indx(lessindx));
    auto joined_mat = std::move(arma::join_cols( tempmat1, tempmat2 ));
    // joined_mat.print();
    // find the minmum of each col
    mins=min(joined_mat,0);
    s1=s1+mins;
    // find the maxmum of each col
    maxs=max(joined_mat,0);
    s2=s2+maxs;
  }
  s1=ck-s1;
  
  s2=s2-ck;
  auto joined_mat2 = std::move(arma::join_cols( s1, s2 ));
  // find the maxmum of each col
  sf=max(joined_mat2,0);
  // lambda=lambda*n*log(n)*2;
  lams=abs(lambda/beta);
  arma:: uvec delinx= find(lams>sf)+1;
  // retain at least one
  
  
  Rcpp::List ret;
  
  ret["delinx"] = delinx;
  ret["m"] = delinx.size();
  ret["minsf"] = min(sf);
  return (ret);
  
  
}

