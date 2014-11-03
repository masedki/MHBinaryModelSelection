//
//  myfunctions.cpp
//  
//
//  Created by sedki on 02/11/2014.
//
//

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

#include "myfunctions.h"

colvec sigmoid(colvec beta, mat X)
{
    colvec s(X.n_rows);
    //colvec s = ones<colvec>(X.n_rows)/(ones<colvec>(X.n_rows) + exp(-(X*beta)));
    Rcout <<"X.n_cols = " << X.n_cols << "   beta.size() = " << beta.size() <<endl;
    
    //Rcout << " X * beta = \n" << X*beta << endl;
    //s = X*beta;
    return(s);
};


