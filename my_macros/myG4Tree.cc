#include "myG4Tree.hh"
#include <vector>

void InitG4Tree (TTree * TreeRun, myG4TreeVars& treeVars)
{
    TreeRun->SetBranchAddress("Energyem",               &treeVars.Energyem);    
    
    TreeRun->SetBranchAddress("EnergyScin",             &treeVars.EnergyScin);    
    TreeRun->SetBranchAddress("EnergyCher",             &treeVars.EnergyCher);    
    TreeRun->SetBranchAddress("NofCherenkovDetected",   &treeVars.NofCherenkovDetected);    
    TreeRun->SetBranchAddress("EnergyTot",              &treeVars.EnergyTot);    
    
    TreeRun->SetBranchAddress("PrimaryParticleEnergy",  &treeVars.PrimaryParticleEnergy);    
    TreeRun->SetBranchAddress("PrimaryParticleName",    &treeVars.PrimaryParticleName);    
    TreeRun->SetBranchAddress("neutrinoleakage",        &treeVars.neutrinoleakage);    
    TreeRun->SetBranchAddress("leakage",                &treeVars.leakage);    
        
    TreeRun->SetBranchAddress("SCEP_EnergyDepF",        &treeVars.SCEP_EnergyDepF);
    TreeRun->SetBranchAddress("SCEP_NCherProdF",        &treeVars.SCEP_NCherProdF);    
    TreeRun->SetBranchAddress("SCEP_EnergyDepR",        &treeVars.SCEP_EnergyDepR);
    TreeRun->SetBranchAddress("SCEP_NCherProdR",        &treeVars.SCEP_NCherProdR);
    
    
    
    treeVars.VecHit_CrystalID       = new std::vector<double>;
    treeVars.VecHit_ScepEneDepF     = new std::vector<double>;
    treeVars.VecHit_ScepEneDepR     = new std::vector<double>;
    treeVars.VecHit_ScepCherF       = new std::vector<double>;
    treeVars.VecHit_ScepCherR       = new std::vector<double>;
    
    treeVars.VecHit_Timing_CrystalID_F = new std::vector<double>;
    treeVars.VecHit_Timing_CrystalID_R = new std::vector<double>;
    treeVars.VecHit_Timing_ScepEneDepF = new std::vector<double>;
    treeVars.VecHit_Timing_ScepEneDepR = new std::vector<double>;
    treeVars.VecHit_Timing_ScepTimeF   = new std::vector<double>;
    treeVars.VecHit_Timing_ScepTimeR   = new std::vector<double>;
    
    treeVars.PrimaryParticleMomentum = new std::vector<double>;
    
    treeVars.VectorL            = new std::vector<double>;
    treeVars.VectorR            = new std::vector<double>;
    treeVars.VectorL_loop       = new std::vector<double>;
    treeVars.VectorR_loop       = new std::vector<double>;
    treeVars.VectorSignalsR     = new std::vector<double>;
    treeVars.VectorSignalsL     = new std::vector<double>;
    treeVars.VectorSignalsCherL = new std::vector<double>;
    treeVars.VectorSignalsCherR = new std::vector<double>;
    
    
    TreeRun->SetBranchAddress("VecHit_CrystalID",       &treeVars.VecHit_CrystalID);    
    TreeRun->SetBranchAddress("VecHit_ScepEneDepF",     &treeVars.VecHit_ScepEneDepF);    
    TreeRun->SetBranchAddress("VecHit_ScepEneDepR",     &treeVars.VecHit_ScepEneDepR);    
    TreeRun->SetBranchAddress("VecHit_ScepCherF",       &treeVars.VecHit_ScepCherF);    
    TreeRun->SetBranchAddress("VecHit_ScepCherR",       &treeVars.VecHit_ScepCherR);    
    
    TreeRun->SetBranchAddress("VecHit_Timing_CrystalID_F",  &treeVars.VecHit_Timing_CrystalID_F);    
    TreeRun->SetBranchAddress("VecHit_Timing_CrystalID_R",  &treeVars.VecHit_Timing_CrystalID_R);    
    TreeRun->SetBranchAddress("VecHit_Timing_ScepEneDepF",  &treeVars.VecHit_Timing_ScepEneDepF);    
    TreeRun->SetBranchAddress("VecHit_Timing_ScepEneDepR",  &treeVars.VecHit_Timing_ScepEneDepR);    
    TreeRun->SetBranchAddress("VecHit_Timing_ScepTimeF",    &treeVars.VecHit_Timing_ScepTimeF);    
    TreeRun->SetBranchAddress("VecHit_Timing_ScepTimeR",    &treeVars.VecHit_Timing_ScepTimeR);    
    
    TreeRun->SetBranchAddress("PrimaryParticleMomentum",&treeVars.PrimaryParticleMomentum);    
    
    
    TreeRun->SetBranchAddress("VectorSignalsR",         &treeVars.VectorSignalsR);        
    TreeRun->SetBranchAddress("VectorSignalsL",         &treeVars.VectorSignalsL);    
    TreeRun->SetBranchAddress("VectorSignalsCherR",     &treeVars.VectorSignalsCherR);    
    TreeRun->SetBranchAddress("VectorSignalsCherL",     &treeVars.VectorSignalsCherL);    
    
    TreeRun->SetBranchAddress("VectorL",                &treeVars.VectorL);        
    TreeRun->SetBranchAddress("VectorR",                &treeVars.VectorR);       
//     TreeRun->SetBranchAddress("VectorL_loop",           &treeVars.VectorL_loop);       
//     TreeRun->SetBranchAddress("VectorR_loop",           &treeVars.VectorR_loop);  
    
//     
    
//     TreeRun->SetBranchStatus("*", 0);
//     
//     TreeRun->SetBranchStatus("SCEP_EnergyDepF", 1);
//     TreeRun->SetBranchStatus("SCEP_NCherProdF", 1);
//     TreeRun->SetBranchStatus("SCEP_EnergyDepR", 1);
//     TreeRun->SetBranchStatus("SCEP_NCherProdR", 1);
//     
//     TreeRun->SetBranchStatus("PrimaryParticleMomentum", 1);
//     
//     
//     TreeRun->SetBranchStatus("EnergyScin", 1);
//     TreeRun->SetBranchStatus("NofCherenkovDetected", 1);
//     
      

    
}
