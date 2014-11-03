//
//  MHLogitRegModelSelection.h
//  
//
//  Created by sedki on 03/11/2014.
//
//

#ifndef ____MHLogitRegModelSelection__
#define ____MHLogitRegModelSelection__

#include <stdio.h>
class MHLogitRegModelSelection{
private:
    imat X;
    colvec y;
    int n;
    int p;
    colvec VCramer;
    ivec CurrentModel;
    ivec BestModel;
    ivec CandidateModel;
    colvec CurrentFit;
    colvec BestFit;
    colvec CandidateFit;
    long double CurrentBic;
    long double BestBic;
    long double CandidateBic;
    int maxit;
    int iter;
    

public:
    void FirstModelFit();
    long double GenerateCandidate();
    void UpdateBicCandidate();
    List MHRun();
    MHLogitRegModelSelection(colvec & yin, imat & Xin, colvec & vcin,  int  maxitin);
    ~MHLogitRegModelSelection(){}; //Destructor
    
    
};
#endif /* defined(____MHLogitRegModelSelection__) */
