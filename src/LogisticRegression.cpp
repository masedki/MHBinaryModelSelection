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
  //NumericVector Ry(yin);
  convertVector<NumericVector,VectorXd>(yin, y);
  //NumericMatrix RX(Xin);
  convertMatrix<NumericMatrix,MatrixXd>(Xin, X);
  //Rcout << "X =  \n" << X << endl; 
  
  
  //this->X = Mat<double>(Xin);
  this->n = X.rows();
  this->p = X.cols();
  //Rcout << " n = ... " << n << "  p =..." << p << endl;
  //X.insert_cols(0,  ones<colvec>(n)); // ajouter la colonne des 1
  //this->y=Col<double>(yin);
  //NumericVector Rbeta(start);
  convertVector<NumericVector,VectorXd>(start, beta);
  //this->beta=Col<double>(start);
  //if(beta.size()==p)
  //{
  //  beta = zeros<vec>(1);  
  //  beta.insert_rows(1,Col<double>(start));
  //};
  
  this->gradient.setZero(p);
  this->hessian.setZero(p, p);
  this->W.setZero(n,n);
  this->iter = 0;
  this->maxit = 1;
  this->tol = 1e-07;
  this->again = true;
  
};

//LogisticRegression::LogisticRegression()
//{
//  this->n=0;  
//  this->p=0;  
//  this->X=zeros<mat>(n,p+1);
//  this->y = zeros<colvec>(n);
//  this->beta = zeros<colvec>(p+1);
//  this->gradient = zeros<colvec>(p+1);
//  this->hessian = zeros<mat>(p+1, p+1);
//  this->W = zeros<mat>(n,n);
//  this->iter = 0;
//  this->maxit = 0;
//  this->tol = 0;
//  this->again = false; 
//};
//
void LogisticRegression::ComputeGradient(){
  
  //Rcout << " tmp = \n " << tmp << endl;
  VectorXd pi = VectorXd::Ones(n) + ((X*beta).array().exp()).matrix();
  //Rcout << " pi = \n " << pi << endl;
  pi = VectorXd::Ones(n).array() / pi.array(); //compute sigmoid
  Rcout << " pi = \n " << pi << endl;
  VectorXd tmp = y - pi;
  gradient = X.adjoint()*tmp;//(y.array() - pi.array()).matrix();
  Rcout << " gradient = \n " << gradient << endl;
  W=MatrixXd::Zero(n,n);
  for (int i=0;i<n;i++){
    W(i,i)=pi(i) * (1-pi(i));
  }
//  W.diagonal() = pi.array()*(VectorXd::Ones(n) - pi).array();
  Rcout << " W = \n " << W << endl;
};
//
//
void LogisticRegression::ComputeHessian()
{
  //hessian.setZero(p,p);
  hessian = (X.transpose()*W)*X;
  Rcout << " hessian = \n " << hessian << endl;
};

void LogisticRegression::NewtonRaphsonUpdate()
{
  iter++;
   VectorXd oldbeta(beta);
   Rcout << " oldbeta = \n " << oldbeta << endl;
   //LLT<MatrixXd> llt;
   //llt.compute(hessian);
   //beta = oldbeta.array() + llt.solve(gradient).array();
   //VectorXd gamma(llt.solve(gradient));
   //beta = oldbeta.array() + gamma.array();
   beta +=  hessian.ldlt().solve(gradient); //A.ldlt().solve(b))
   
   Rcout << " beta = \n " << beta << endl;
   //VectorXd tmp = beta.array() - oldbeta.array();
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
