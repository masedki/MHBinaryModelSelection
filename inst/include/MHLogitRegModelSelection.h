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
    mat X;
    colvec y;
    int n;
    int p;
    vec VCramer;
    vec CurrentModel;
    vec BestModel;
    colvec CurrentFit;
    colvecBestFit;
    long double BestBic;
    long double CurrentBic;
    int maxit;

public:
    double GenerateCandidate();
    void UpdateProbaCandidate();
    MHLogitRegModelSelection(colvec & yin, mat & Xin, int & maxitin);
    ~MHLogitRegModelSelection(){}; //Destructor
    
    
};
#endif /* defined(____MHLogitRegModelSelection__) */
