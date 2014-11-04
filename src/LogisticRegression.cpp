//
//  LogisticRegression.cpp
//  
//
//  Created by sedki on 02/11/2014.
//

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>


using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace Eigen;
#include "conversion.h"
#include "LogisticRegression.h"

LogisticRegression::LogisticRegression(NumericVector & yin, NumericMatrix & Xin, NumericVector & start)
{
  convertVector<NumericVector,VectorXd>(yin, y);
  convertMatrix<NumericMatrix,MatrixXd>(Xin, X);
  this->n = X.rows();
  this->p = X.cols();
  convertVector<NumericVector,VectorXd>(start, beta);
  this->gradient.setZero(p);
  this->hessian.setZero(p, p);
  this->W.setZero(n,n);
  this->iter = 0;
  this->maxit = 5;
  this->tol = 1e-03;
  this->again = true;
  
};


void LogisticRegression::ComputeGradient(){
  
  
  VectorXd xb = X*beta;
  VectorXd pi = VectorXd::Zero(n);
  for(int i = 0; i < n; ++i)
      pi(i) = exp(xb(i))/(1+exp(xb(i)));

  VectorXd tmp = y - pi;
  gradient = X.transpose()*tmp;
  for(int i = 0; i < n; ++i)
  W(i,i) = pi(i)*(1 - pi(i));
  
};

void LogisticRegression::ComputeHessian()
{
  hessian = (X.transpose()*W)*X;
  
};

void LogisticRegression::NewtonRaphsonUpdate()
{
  iter++;
   VectorXd oldbeta(beta);
   LLT<MatrixXd> llt;
   llt.compute(hessian);
   //beta +=  hessian.ldlt().solve(gradient); //A.ldlt().solve(b))
   beta += llt.solve(gradient);
 if(((beta.array() - oldbeta.array()).array().abs().sum()/p < tol) || (iter > maxit))
   again = false;
};
//
VectorXd LogisticRegression::ShowCoefficient()
{
  return(beta);
};
//
//
bool LogisticRegression::Continue()
{
  return(again);
}
//
long double LogisticRegression::twiceloglik()
{
  VectorXd xb = X*beta;
  
  long double l = 2 * (y.array()*xb.array() - (VectorXd::Ones(n).array() + xb.array().exp()).array().log()).sum();
  return(l);  
};
//
//
void LogisticRegression::Run()
{
  
  while(Continue()==true)
  {
    ComputeGradient();
    ComputeHessian();
    NewtonRaphsonUpdate();
  }

};
//

//[[Rcpp::export]]
List RcppLogisticRegression(NumericVector & Y, NumericMatrix & X, NumericVector & beta)
{
  LogisticRegression myReg(Y,X,beta);
  while(myReg.Continue()==true)
  {
    myReg.ComputeGradient();
    myReg.ComputeHessian();
    myReg.NewtonRaphsonUpdate();
  }


  return List::create(Named("twiceloglik") = myReg.twiceloglik(),
                      Named("coefficients") = myReg.ShowCoefficient());  
  
};
