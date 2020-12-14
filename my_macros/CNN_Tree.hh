
#ifndef __CNN_Tree__
#define __CNN_Tree__

#include <vector>
#include <string>
#include <algorithm>

#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"

struct myCNNTreeVars
{

    /// defining tree variables                              
    
    float  PrimaryParticleEnergy;
    float  PrimaryParticleMomentum[3];
    std::string PrimaryParticleName;
//     Char_t * PrimaryParticleName;
        
    float theta_seed;
    float phi_seed;    
    
    float image_TT[225];
    
    float image_E1[225];
    float image_E2[225];
//     float image_ET[225];
    
    float image_DRT_S[225];
    float image_DRT_C[225];
    
//     float image_E1[30][30];
//     float image_E2[30][30];
//     float image_ET[30][30];
    
          
};

void InitCNNTree(TTree* TreeRun, myCNNTreeVars &treeVars);


#endif

