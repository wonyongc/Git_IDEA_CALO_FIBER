/*
#ifndef __myTruthTree__
#define __myTruthTree__*/

#include <vector>
#include <string>
#include <algorithm>

#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"

struct myTruthTreeVars
{

    /// defining tree variables            
                  
    UInt_t mcs_n;
    
    std::vector<float> *mcs_E;
    std::vector<float> *mcs_pt;
    std::vector<float> *mcs_m;
    std::vector<float> *mcs_eta;
    std::vector<float> *mcs_phi;    
    std::vector<float> *mcs_charge;
    std::vector<float> *mcs_vx_x;
    std::vector<float> *mcs_vx_y;
    std::vector<float> *mcs_vx_z;        
    
    std::vector<int> *mcs_status;
    std::vector<int> *mcs_barcode;
    std::vector<int> *mcs_pdgId;
  
};

void InitTruthTree(TTree* TreeRun, myTruthTreeVars &treeVars);


// #endif
