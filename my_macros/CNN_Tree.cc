#include "CNN_Tree.hh"
#include <vector>
#include <string>

void InitCNNTree (TTree * TreeRun, myCNNTreeVars& treeVars)
{
        
    TreeRun->Branch("PrimaryParticleEnergy",   &treeVars.PrimaryParticleEnergy, "PrimaryParticleEnergy/F");        
    TreeRun->Branch("PrimaryParticleMomentum", treeVars.PrimaryParticleMomentum, "PrimaryParticleMomentum[3]/F");        
    
    TreeRun->Branch("PrimaryParticleName",    &treeVars.PrimaryParticleName);        
//     TreeRun->Branch("PrimaryParticleName",     &treeVars.PrimaryParticleName, "PrimaryParticleName/C", 1024);        
    
    TreeRun->Branch("theta_seed", &treeVars.theta_seed, "theta_seed/F");        
    TreeRun->Branch("phi_seed",   &treeVars.phi_seed,   "phi_seed/F");        
    
    TreeRun->Branch("image_TT", treeVars.image_TT, "image_TT[225]/F");        
    
    TreeRun->Branch("image_E1", treeVars.image_E1, "image_E1[225]/F");        
    TreeRun->Branch("image_E2", treeVars.image_E2, "image_E2[225]/F");        
//     TreeRun->Branch("image_ET", treeVars.image_ET, "image_ET[225]/F");        
    
    TreeRun->Branch("image_DRT_S", treeVars.image_DRT_S, "image_DRT_S[225]/F");        
    TreeRun->Branch("image_DRT_C", treeVars.image_DRT_C, "image_DRT_C[225]/F");        
    
//     TreeRun->Branch("image_E1", treeVars.image_E1, "image_E1[30][30]/F");        
//     TreeRun->Branch("image_E2", treeVars.image_E2, "image_E2[30][30]/F");        
//     TreeRun->Branch("image_ET", treeVars.image_ET, "image_ET[30][30]/F");        
    
    
}
