//
//  LogisticRegression.h
//  
//
//  Created by sedki on 02/11/2014.
//
//

#ifndef _LogisticRegression_h
#define _LogisticRegression_h


class LogisticRegression{
private:
    mat X;
    colvec y;
    colvec beta;
    int n;
    int p;
    colvec gradient;
    mat hessian;
    mat W;
    int iter;
    int maxit;
    long double tol;
    bool again;
public:
    
    
    
    void ComputeGradient();
    void ComputeHessian();
    void NewtonRaphsonUpdate();
    colvec ShowCoefficient(); 
    bool Continue(); 
    long double twiceloglik();
    void Run();
    LogisticRegression(colvec & yin,  mat & Xin, colvec & start); // Constructor
    LogisticRegression();
    ~LogisticRegression(){} ; // Destructor
    
};

#endif
