//
//  LogisticRegression.h
//  
//
//  Created by sedki on 02/11/2014.
//
//

#ifndef _LogisticRegression_h
#define _LogisticRegression_h
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

class LogisticRegression{
private:
    MatrixXd X;
    VectorXd y;
    VectorXd beta;
    int n;
    int p;
    VectorXd gradient;
    MatrixXd hessian;
    MatrixXd W;
    int iter;
    int maxit;
    long double tol;
    bool again;
public:
    
    
    
    void ComputeGradient();
    void ComputeHessian();
    void NewtonRaphsonUpdate();
    VectorXd ShowCoefficient(); 
    bool Continue(); 
    long double twiceloglik();
    void Run();
    LogisticRegression(NumericVector & yin,  NumericMatrix & Xin, NumericVector & start); // Constructor
    LogisticRegression();
    ~LogisticRegression(){} ; // Destructor
    
};

#endif
