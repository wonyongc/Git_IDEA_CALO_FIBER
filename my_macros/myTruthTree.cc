#include "myTruthTree.hh"
#include <vector>
#include "TTree.h"
#include "TChain.h"

void InitTruthTree (TTree * TreeRun, myTruthTreeVars& treeVars)
{
    TreeRun->SetBranchAddress("mcs_n",               &treeVars.mcs_n);    
    
    treeVars.mcs_E       = new std::vector<float>;
    treeVars.mcs_pt      = new std::vector<float>;
    treeVars.mcs_m       = new std::vector<float>;
    treeVars.mcs_eta     = new std::vector<float>;
    treeVars.mcs_phi     = new std::vector<float>;
    treeVars.mcs_charge  = new std::vector<float>;
    treeVars.mcs_vx_x    = new std::vector<float>;
    treeVars.mcs_vx_y    = new std::vector<float>;
    treeVars.mcs_vx_z    = new std::vector<float>;;
    
    treeVars.mcs_status  = new std::vector<int>;
    treeVars.mcs_barcode = new std::vector<int>;
    treeVars.mcs_pdgId   = new std::vector<int>;
    
    
    TreeRun->SetBranchAddress("mcs_E",       &treeVars.mcs_E);
    TreeRun->SetBranchAddress("mcs_pt",      &treeVars.mcs_pt);
    TreeRun->SetBranchAddress("mcs_m",       &treeVars.mcs_m);
    TreeRun->SetBranchAddress("mcs_eta",     &treeVars.mcs_eta);
    TreeRun->SetBranchAddress("mcs_phi",     &treeVars.mcs_phi);
    TreeRun->SetBranchAddress("mcs_charge",  &treeVars.mcs_charge);
    TreeRun->SetBranchAddress("mcs_vx_x",    &treeVars.mcs_vx_x);
    TreeRun->SetBranchAddress("mcs_vx_y",    &treeVars.mcs_vx_y);
    TreeRun->SetBranchAddress("mcs_vx_z",    &treeVars.mcs_vx_z);
    
    TreeRun->SetBranchAddress("mcs_status",  &treeVars.mcs_status);
    TreeRun->SetBranchAddress("mcs_barcode", &treeVars.mcs_barcode);
    TreeRun->SetBranchAddress("mcs_pdgId",   &treeVars.mcs_pdgId);

}
