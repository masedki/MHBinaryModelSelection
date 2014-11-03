//
//  MHLogitRegModelSelection.cpp
//  
//
//  Created by sedki on 03/11/2014.
//
//
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

#include "LogisticRegression.h"
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




MHLogitRegModelSelection::MHLogitRegModelSelection(colvec & yin, imat & Xin, colvec & vcin, int  maxitin)
{
  this->y = Col<double>(yin); 
  this->X = Mat<int>(Xin); 
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
  Rcout<< " ... j'ai créé mon objet ..."<<endl;
};



void MHLogitRegModelSelection::FirstModelFit()
{
  CurrentModel = zeros<ivec>(p);
  for(int j = 0; j < p; ++j)
     if(runif() > 0.5)
       CurrentModel(j) = 1;
  
  CandidateModel = Col<int>(CurrentModel);
  uvec idx = find(CurrentModel==1);
  mat Xc = Mat<double>(conv_to<mat>::from(X.cols(idx)));
  colvec start = zeros<colvec>(sum(CandidateModel)+1);
  LogisticRegression MyReg(y, Xc, start);
  MyReg.Run();
  CurrentBic = MyReg.twiceloglik() - log(n)*Xc.n_cols;
  BestBic =  CurrentBic;
  BestModel = Col<int>(CurrentModel);
  CurrentFit = Col<double>(MyReg.ShowCoefficient());
  BestFit = CurrentFit;
  Rcout<< " ... j'ai réussi mon first fit ..."<<endl;
};

long double MHLogitRegModelSelection::GenerateCandidate()
{
  CandidateModel = Col<int>(CurrentModel);
  int j=(rand()%CandidateModel.size());
  int add=(rand()%2);
  
  colvec probatire = zeros<colvec>(CandidateModel.size());
  colvec alternative = zeros<colvec>(CandidateModel.size());
  uvec idx1 = find(CandidateModel==1);
  uvec idx0 = find(CandidateModel==0);
  if (add==1){
    //  probatire[candidate.model==0] <- v.cramer[candidate.model==0]  
    //  alternative[candidate.model==1] <- 1 - v.cramer[candidate.model==1]
    
    probatire.rows(idx0) = VCramer.rows(idx0);
    alternative.rows(idx1) = ones<colvec>(idx1.size()) - VCramer.rows(idx1);
  }
  else{
    // probatire[candidate.model==1] <- 1-v.cramer[candidate.model==1]
    //alternative[candidate.model==0] <- v.cramer[candidate.model==0]
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
//    X <- x[,(1:p)[candidate.model==1]]
//    vars <- paste("X[,",1:ncol(X),"]",sep="")
//    fla <- paste("y~", paste(vars, collapse="+"))
//    current.coef <- rep(0,p)
//    current.coef[(1:p)[current.model==1]] <- current.fit$coef[-1]
//    current.coef <- c(current.fit$coef[1], current.coef)
//    coef.start <- c(current.coef[1],current.coef[-1][(1:p)[candidate.model==1]])
//    candidate.fit <- glm(as.formula(fla), family = binomial, data = data.frame(y=y,X=X), start = coef.start)
//    twice.log.lik <- -candidate.fit$aic + 2*length(candidate.fit$coef)
//    candidate.model.bic <-  twice.log.lik - log(n)*length(candidate.fit$coef)
 uvec idxcurr = find(CurrentModel==1);
 colvec CurrentCoef = zeros<colvec>(p);
 colvec CurrentFit_1 = Col<double>(CurrentFit);
 CurrentFit_1.shed_row(0);
 CurrentCoef.cols(idxcurr) = CurrentFit_1;
 
 uvec idxcand = find(CandidateModel==1);
 Rcout << "idxcand.size() = .... " << idxcand.size() << "  sum(CandidateModel) = ..." << sum(CandidateModel) << endl;
 colvec start = CurrentCoef(idxcand);
 Rcout << "start = .... " << start << endl;
 start.insert_rows(0,Col<double>(CurrentFit(0)));
 mat Xc = Mat<double>(conv_to<mat>::from(X.cols(idxcand)));
  
 
 LogisticRegression MyReg(y, Xc, start);
 MyReg.Run();
 CandidateBic = MyReg.twiceloglik() - log(n)*Xc.n_cols;
 //BestBic =  CurrentBic;
 CandidateFit = Col<double>(MyReg.ShowCoefficient());  
  Rcout<< " ... j'ai fitté mon candidat ..."<<endl;
}

List MHLogitRegModelSelection::MHRun()
{
  long double rho, normalise;
  while(iter < maxit)
  {
    iter++;
    // Generate the proposal 
    normalise = GenerateCandidate();
    UpdateBicCandidate();
    // the acceptance probability
    rho = normalise*exp(CandidateBic - CurrentBic);
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
List RcppMHModelSelection(arma::colvec Y, arma::imat X, arma::colvec VCin, int MAXIT)
{
 MHLogitRegModelSelection MyMH(Y, X, VCin, MAXIT);
 MyMH.FirstModelFit();
 return(MyMH.MHRun());
}
