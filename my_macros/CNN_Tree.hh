
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
    std::string CNNPrimaryParticleName;
        
    float theta_seed;
    float phi_seed;    
    
    float image_TT[2025];    
    float image_E1[2025];
    float image_E2[2025];
    
//     float image_DRT_S[2025];
//     float image_DRT_C[2025];
    
          
};

void InitCNNTree(TTree* CNNTreeRun, myCNNTreeVars &CNNtreeVars);


#endif

