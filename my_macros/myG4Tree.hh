
#ifndef __myG4Tree__
#define __myG4Tree__

#include <vector>
#include <string>
#include <algorithm>

#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"

struct myG4TreeVars
{

    /// defining tree variables            
                  
    double Energyem;
    double EnergyScin;
    double EnergyCher;
    
    Double_t NofCherenkovDetected;
    double EnergyTot;
    double PrimaryParticleEnergy;
    Char_t PrimaryParticleName;
    double neutrinoleakage;
    double leakage;
    
    double SCEP_EnergyDepF;
    Double_t SCEP_NCherProdF;
    double SCEP_EnergyDepR;
    Double_t SCEP_NCherProdR;
    
    
    std::vector<double> *VecHit_CrystalID;
    std::vector<double> *VecHit_ScepEneDepF;
    std::vector<double> *VecHit_ScepEneDepR;
    std::vector<double> *VecHit_ScepCherF;
    std::vector<double> *VecHit_ScepCherR;
    
    std::vector<double> *VecHit_Timing_CrystalID_F;
    std::vector<double> *VecHit_Timing_CrystalID_R;
    std::vector<double> *VecHit_Timing_ScepEneDepF;
    std::vector<double> *VecHit_Timing_ScepEneDepR;
    std::vector<double> *VecHit_Timing_ScepTimeF;
    std::vector<double> *VecHit_Timing_ScepTimeR;
    
    std::vector<double> *PrimaryParticleMomentum;
    
    std::vector<double> *VectorR_loop;
    std::vector<double> *VectorL_loop;

    std::vector<double> *VectorSignalsR;//Vector filled with scintillating fibers energy deposits
    std::vector<double> *VectorSignalsL;//vector filled for left side
    std::vector<double> *VectorSignalsCherR;//Vector filled with Cherenkov fibers Cherenkov photoelectrons
    std::vector<double> *VectorSignalsCherL;//vector filled for left side
    
    std::vector<double> *VectorR; //vector with energy deposited in towers right
    std::vector<double> *VectorL;   
  
};

void InitG4Tree(TTree* TreeRun, myG4TreeVars &treeVars);


#endif
