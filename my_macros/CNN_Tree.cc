#include "CNN_Tree.hh"
#include <vector>
#include <string>

void InitCNNTree (TTree * CNNTreeRun, myCNNTreeVars& CNNtreeVars)
{
        
    CNNTreeRun->Branch("PrimaryParticleEnergy",   &CNNtreeVars.PrimaryParticleEnergy, "PrimaryParticleEnergy/F");        
    CNNTreeRun->Branch("PrimaryParticleMomentum", CNNtreeVars.PrimaryParticleMomentum, "PrimaryParticleMomentum[3]/F");        
    
    CNNTreeRun->Branch("PrimaryParticleName",    &CNNtreeVars.CNNPrimaryParticleName);        
    
    CNNTreeRun->Branch("theta_seed", &CNNtreeVars.theta_seed, "theta_seed/F");        
    CNNTreeRun->Branch("phi_seed",   &CNNtreeVars.phi_seed,   "phi_seed/F");        
    
    CNNTreeRun->Branch("image_TT", CNNtreeVars.image_TT, "image_TT[225]");            
    CNNTreeRun->Branch("image_E1", CNNtreeVars.image_E1, "image_E1[225]");        
    CNNTreeRun->Branch("image_E2", CNNtreeVars.image_E2, "image_E2[225]");        
    
//     CNNTreeRun->Branch("image_TT", CNNtreeVars.image_TT, "image_TT[2025]");            
//     CNNTreeRun->Branch("image_E1", CNNtreeVars.image_E1, "image_E1[2025]");        
//     CNNTreeRun->Branch("image_E2", CNNtreeVars.image_E2, "image_E2[2025]");        
    
//     CNNTreeRun->Branch("image_DRT_S", CNNtreeVars.image_DRT_S, "image_DRT_S[2025]");        
//     CNNTreeRun->Branch("image_DRT_C", CNNtreeVars.image_DRT_C, "image_DRT_C[2025]");        

    
}
