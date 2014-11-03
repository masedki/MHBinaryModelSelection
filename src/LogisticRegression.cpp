//
//  LogisticRegression.cpp
//  
//
//  Created by sedki on 02/11/2014.
//

#include <RcppArmadillo.h>

#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>


using namespace std;
using namespace Rcpp;
using namespace arma;

#include "LogisticRegression.h"

LogisticRegression::LogisticRegression(colvec & yin, mat & Xin, colvec & start)
{
  this->X = Mat<double>(Xin);
  this->n = X.n_rows;
  this->p = X.n_cols;
  X.insert_cols(0,  ones<colvec>(n)); // ajouter la colonne des 1
  this->y=Col<double>(yin);
  this->beta=Col<double>(start);
  if(beta.size()==p)
  {
    beta = zeros<vec>(1);  
    beta.insert_rows(1,Col<double>(start));
  };
  
  this->gradient = zeros<colvec>(p+1);
  this->hessian = zeros<mat>(p+1, p+1);
  this->W = zeros<mat>(n,n);
  this->iter = 0;
  this->maxit = 2000;
  this->tol = 0.0001;
  this->again = true;
  
};

LogisticRegression::LogisticRegression()
{
  this->n=0;  
  this->p=0;  
  this->X=zeros<mat>(n,p+1);
  this->y = zeros<colvec>(n);
  this->beta = zeros<colvec>(p+1);
  this->gradient = zeros<colvec>(p+1);
  this->hessian = zeros<mat>(p+1, p+1);
  this->W = zeros<mat>(n,n);
  this->iter = 0;
  this->maxit = 0;
  this->tol = 0;
  this->again = false; 
};

void LogisticRegression::ComputeGradient(void){
  
  colvec pi = ones<colvec>(n)/(ones<colvec>(n) + exp(-(X*beta))); //compute sigmoid
  gradient = X.t()*(y - pi);
  W.diag() = pi%(ones(n) - pi);
  
};


void LogisticRegression::ComputeHessian(void)
{
  hessian = X.t()*W*X;
};

void LogisticRegression::NewtonRaphsonUpdate(void)
{
  iter++;
  colvec oldbeta(beta);
  beta = oldbeta + inv_sympd(hessian)*gradient;
  if((accu(abs(beta - oldbeta)) < tol) || (iter > maxit))
   again = false;
};

vec LogisticRegression::ShowCoefficient(void)
{
  return(beta);
};


bool LogisticRegression::Continue(void)
{
  return(again);
}

long double LogisticRegression::twiceloglik(void)
{
  colvec xb = X*beta;
  long double l = 2 * sum(y%xb - log(ones<colvec>(n)+exp(xb)));
  return(l);  
};


//arma::vec RcppLogisticRegression(arma::colvec Y, arma::mat X, arma::colvec beta)
//[[Rcpp::export]]
long double RcppLogisticRegression(arma::colvec Y, arma::mat X, arma::colvec beta)
{
  LogisticRegression myReg(Y,X,beta);
  while(myReg.Continue()==true)
  {
    myReg.ComputeGradient();
    myReg.ComputeHessian();
    myReg.NewtonRaphsonUpdate();
  }

  //return(myReg.ShowCoefficient());
  return(myReg.twiceloglik());
};
