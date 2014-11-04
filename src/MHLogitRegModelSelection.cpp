//
//  MHLogitRegModelSelection.cpp
//  
//
//  Created by sedki on 03/11/2014.
//
//
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppEigen.h>
//#include <Rmath.h>
//#include <Rdefines.h>

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace Eigen;
//#include "LogisticRegression.h"
#include "MHLogitRegModelSelection.h"


// Tool function to sample according U[0,1]
double static runif(){
  double r = ((double) rand())/((double) RAND_MAX);
  return r; 
};


// Tool function to sample according to a multinomial distribution whose parameters are given by P
int rmultinom(const colvec & P){
  double thresh = runif();
  int k  = 0;
  double cdf = P(0);
  while( cdf < thresh ){
    k++;
    cdf += P(k);
  }
  return k;
};

Environment Rstats("package:stats");
Environment base("package:base");
Function glm = Rstats["glm"];
Function asformula = Rstats["as.formula"];
Function dataframe = base["data.frame"];
MHLogitRegModelSelection::MHLogitRegModelSelection(colvec & yin, mat & Xin, colvec & vcin, int  maxitin)
{
  this->y = Col<double>(yin); 
  this->X = Mat<double>(Xin); 
  this->n = X.n_rows;
  this->p = X.n_cols;
  this->VCramer = Col<double>(vcin);
  this->CurrentModel =  Col<int>(p);
  this->BestModel =  Col<int>(p);
  this->CandidateModel =  Col<int>(p);
  this->CurrentFit =  Col<double>(p+1);
  this->BestFit =  Col<double>(p+1);
  this->CandidateFit =  Col<double>(p+1);
  this->CurrentBic = log(0.);
  this->BestBic = log(0.);
  this->CandidateBic = log(0.);
  this->maxit = maxitin;
  this->iter = 0;
};



void MHLogitRegModelSelection::FirstModelFit()
{
  CurrentModel = zeros<ivec>(p);
  for(int j = 0; j < p; ++j)
  if(runif() > 0.5)
  CurrentModel(j) = 1;
  
  while(sum(CurrentModel)==0)
  {
    CurrentModel = zeros<ivec>(p);
    for(int j = 0; j < p; ++j)
    if(runif() > 0.5)
    CurrentModel(j) = 1;
  }
  CandidateModel = Col<int>(CurrentModel);
  uvec idx = find(CurrentModel==1);
  Rcout << " idx  = \n" << idx << endl;
  mat Xc = Mat<double>(X.cols(idx));
  colvec start = zeros<colvec>(sum(CandidateModel)+1);
  List MyReg = glm(Named("formula") = asformula("y~."), 
  Named("data") = dataframe(Named("y") = wrap(y), 
  Named("x") = wrap(Xc)),  
  Named("start") = wrap(start));
  
  //LogisticRegression MyReg(y, Xc, start);
  //MyReg.Run();
  //CurrentBic = MyReg.twiceloglik() - log(n)*(idx.size()+1);
  Rcout << " MyReg[aic]  = " << as<double>(MyReg["aic"]) << endl;
  CurrentBic = -as<double>(MyReg["aic"]) + (2 - log(n))*(idx.size()+1);
  BestBic =  CurrentBic;
  BestModel = Col<int>(CurrentModel);
  //CurrentFit = Col<double>(MyReg.ShowCoefficient());
  CurrentFit = as<colvec>(MyReg["coefficients"]);
  BestFit = CurrentFit;
  
};

double MHLogitRegModelSelection::GenerateCandidate()
{
  CandidateModel = Col<int>(CurrentModel);
  int j=(rand()%CandidateModel.size());
  int add=(rand()%2);
  
  colvec probatire = zeros<colvec>(CandidateModel.size());
  colvec alternative = zeros<colvec>(CandidateModel.size());
  uvec idx1 = find(CandidateModel==1);
  uvec idx0 = find(CandidateModel==0);
  if (add==1){
    probatire.rows(idx0) = VCramer.rows(idx0);
    alternative.rows(idx1) = ones<colvec>(idx1.size()) - VCramer.rows(idx1);
  }
  else{
    probatire.rows(idx1) = ones<colvec>(idx1.size())-VCramer.rows(idx1);
    alternative.rows(idx0) = VCramer.rows(idx0);
  }
  
  probatire = probatire/sum(probatire);
  int drug = rmultinom(probatire);  
  if(add)
  alternative(drug) = VCramer(drug);  
  else
  alternative(drug) = 1 - VCramer(drug);
  
  alternative = alternative/sum(alternative);
  CandidateModel(drug) = 1 - CandidateModel(drug);
  CandidateBic = log(0);
  long double  normalise = alternative(drug)/probatire(drug);
  return normalise;
}

void MHLogitRegModelSelection::UpdateBicCandidate()
{
  uvec idxcurr = find(CurrentModel==1);
  colvec CurrentCoef = zeros<colvec>(p);
  colvec CurrentFit_1 = Col<double>(CurrentFit);
  CurrentFit_1.shed_row(0);
  CurrentCoef.rows(idxcurr) = CurrentFit_1;
  uvec idxcand = find(CandidateModel==1);
  colvec start = zeros<colvec>(idxcand.size()+1);
  start.rows(1,idxcand.size()) = CurrentCoef(idxcand);
  start(0)= CurrentFit(0);
  mat Xc = Mat<double>(X.cols(idxcand));
  
  
  //LogisticRegression MyReg(y, Xc, start);
  //MyReg.Run();
  //CandidateBic = MyReg.twiceloglik() - log(n)*(idxcand.size()+1);
  //CandidateFit = Col<double>(MyReg.ShowCoefficient());  
  List MyReg = glm(Named("formula") = asformula("y~."), 
  Named("data") = dataframe(Named("y") = wrap(y), 
  Named("x") = wrap(Xc)),  
  Named("start") = wrap(start));
  
  //LogisticRegression MyReg(y, Xc, start);
  //MyReg.Run();
  //CurrentBic = MyReg.twiceloglik() - log(n)*(idx.size()+1);
  CandidateBic = - as<double>(MyReg["aic"]) + (2 - log(n))*(idxcand.size()+1);
  CandidateFit = as<colvec>(MyReg["coefficients"]);
}

List MHLogitRegModelSelection::MHRun()
{
  double rho, normalise;
  while(iter < maxit)
  {
    iter++;
    // Generate the proposal 
    normalise = GenerateCandidate();
    UpdateBicCandidate();
    // the acceptance probability
    rho = normalise*exp(CandidateBic - CurrentBic);
    //Rcout << "rho = .... " << rho << endl;
    if(runif() < rho)
    {
      CurrentModel = CandidateModel;
      CurrentFit = CandidateFit;
      CurrentBic = CandidateBic;
      if(CurrentBic > BestBic)
      {
        Rcout <<"I spent ..."<< iter<<" ... iterations to get bic = ..."<<CurrentBic << endl;
        iter = 0;
        BestModel = CurrentModel;
        BestFit = CurrentFit;
        BestBic = CurrentBic;
      }
    }
    
  }
  return List::create(Named("BestBic") = BestBic,
  Named("BestFit") = BestFit,
  Named("BestModel") = BestModel);
}

//[[Rcpp::export]]
List RcppMHModelSelection(arma::colvec Y, arma::mat X, arma::colvec VCin, int MAXIT)
{
  MHLogitRegModelSelection MyMH(Y, X, VCin, MAXIT);
  MyMH.FirstModelFit();
  return(MyMH.MHRun());
}
