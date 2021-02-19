// g++ -Wall -o plotMomentumVectorReco plotMomentumVectorReco.C  VectorSmallestInterval.cc myG4Tree.cc myG4Tree.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs`


#include "VectorSmallestInterval.hh"
#include "SCEPCal_GeometryHelper.hh"
#include "myG4Tree.hh"

#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <utility>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TApplication.h"
#include "TFormula.h"
#include "TRandom.h"

#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"

#include "TSpline.h"
#include "TObject.h"
#include "TVector3.h"



int main(int argc, char** argv)
{

  TApplication* theApp = new TApplication("App", &argc, argv);
      
  using namespace std;
  
  gStyle->SetTitleXOffset (1.00) ;                                                                                        
  gStyle->SetTitleYOffset (1.2) ;                                                                                                                                                                                                                 
  gStyle->SetPadLeftMargin (0.13) ;                                                                                       
  gStyle->SetPadBottomMargin (0.13) ;                                                                                                                                                                                                              
  gStyle->SetTitleSize (0.05, "xyz") ;                                                                                    
  gStyle->SetLabelSize (0.035,"xyz") ;  
    
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.035);
  TLegend * leg;
    
  
  
  bool SAVEPLOTS = false;
  
  const int NFILES = 1;
  std::cout << " NFILES = " << NFILES << std::endl;
//   int energies[NFILES] = {50};
  int energies[NFILES] = {100};
  
//   std::string particle_name = "#mu^{-}, 0T";
  std::string particle_name = "#gamma, 0T";
    
  //define histos
  TH1F * hNHits_ECALBarrel[NFILES];
  TH1F * hNHits_ECALEndcap[NFILES];  
  TH1F * hNHits_T1Barrel[NFILES];
  TH1F * hNHits_T2Barrel[NFILES];
  TH1F * hNHits_T1Endcap[NFILES];
  TH1F * hNHits_T2Endcap[NFILES];
  
  TH1F * hDeltaPhi_Tower[NFILES];
  TH1F * hDeltaEta_Tower[NFILES];
  TH1F * hDeltaTheta_Tower[NFILES];
  
  TH1F * hDeltaPhi_Barrel[NFILES];
  TH1F * hDeltaEta_Barrel[NFILES];
  TH1F * hDeltaTheta_Barrel[NFILES];
  TH1F * hDeltaPhi_Endcap[NFILES];
  TH1F * hDeltaEta_Endcap[NFILES];  
  TH1F * hDeltaTheta_Endcap[NFILES];  
  
  TH1F * hDeltaPhi_Barrel_wGravity[NFILES];
  TH1F * hDeltaEta_Barrel_wGravity[NFILES];
  TH1F * hDeltaTheta_Barrel_wGravity[NFILES];
  TH1F * hDeltaPhi_Endcap_wGravity[NFILES];
  TH1F * hDeltaEta_Endcap_wGravity[NFILES];
  TH1F * hDeltaTheta_Endcap_wGravity[NFILES];
  
  TH1F * hDeltaPhi_TimingBarrel_Comb[NFILES];
  TH1F * hDeltaEta_TimingBarrel_Comb[NFILES];
  TH1F * hDeltaTheta_TimingBarrel_Comb[NFILES];
  TH1F * hDeltaPhi_TimingEndcap_Comb[NFILES];
  TH1F * hDeltaEta_TimingEndcap_Comb[NFILES];  
  TH1F * hDeltaTheta_TimingEndcap_Comb[NFILES];  
  
  TH1F * hDeltaPhi_TimingBarrel_F[NFILES];
  TH1F * hDeltaEta_TimingBarrel_F[NFILES];
  TH1F * hDeltaTheta_TimingBarrel_F[NFILES];
  TH1F * hDeltaPhi_TimingEndcap_F[NFILES];
  TH1F * hDeltaEta_TimingEndcap_F[NFILES];
  TH1F * hDeltaTheta_TimingEndcap_F[NFILES];
  
  TH1F * hDeltaPhi_TimingBarrel_R[NFILES];
  TH1F * hDeltaEta_TimingBarrel_R[NFILES];
  TH1F * hDeltaTheta_TimingBarrel_R[NFILES];
  TH1F * hDeltaPhi_TimingEndcap_R[NFILES];
  TH1F * hDeltaEta_TimingEndcap_R[NFILES];    
  TH1F * hDeltaTheta_TimingEndcap_R[NFILES];
    
  TH1F * hDeltaPhi_TimingBarrel_W_F[NFILES];
  TH1F * hDeltaEta_TimingBarrel_W_F[NFILES];
  TH1F * hDeltaTheta_TimingBarrel_W_F[NFILES];
  TH1F * hDeltaPhi_TimingEndcap_W_F[NFILES];
  TH1F * hDeltaEta_TimingEndcap_W_F[NFILES];
  TH1F * hDeltaTheta_TimingEndcap_W_F[NFILES];
  
  TH1F * hDeltaPhi_TimingBarrel_W_R[NFILES];
  TH1F * hDeltaEta_TimingBarrel_W_R[NFILES];
  TH1F * hDeltaTheta_TimingBarrel_W_R[NFILES];
  TH1F * hDeltaPhi_TimingEndcap_W_R[NFILES];
  TH1F * hDeltaEta_TimingEndcap_W_R[NFILES];
  TH1F * hDeltaTheta_TimingEndcap_W_R[NFILES];
  
  TH2F * hScatterTOF_TimingBarrel[NFILES];    
  TH2F * hScatterTOF_TimingEndcap[NFILES];
  
  TH1F * hTOF_TimingBarrel[NFILES];    
  TH1F * hTOF_TimingEndcap[NFILES];
  
//   TH2F * hScatterPosTime = new TH2F("hScatterPosTime", "hScatterPosTime", 1000, -2000, 2000, 1000, -2000, 2000);
  
  
  int NBINS = 2000;
  double scale = 1;
  
  for (int iFile = 0; iFile<NFILES; iFile++)
  {
      
      std::cout << "energies[" << iFile << "] = " << energies[iFile] << std::endl;
      
      hNHits_ECALBarrel[iFile] = new TH1F (Form("hNHits_ECALBarrel_%d", energies[iFile]), Form("hNHits_ECALBarrel_%d", energies[iFile]), 100, -0.5, 100-0.5);
      hNHits_ECALEndcap[iFile] = new TH1F (Form("hNHits_ECALEndcap_%d", energies[iFile]), Form("hNHits_ECALEndcap_%d", energies[iFile]), 100, -0.5, 100-0.5);
      hNHits_T1Barrel[iFile] = new TH1F (Form("hNHits_T1Barrel_%d", energies[iFile]), Form("hNHits_T1Barrel_%d", energies[iFile]), 100, -0.5, 100-0.5);
      hNHits_T2Barrel[iFile] = new TH1F (Form("hNHits_T2Barrel_%d", energies[iFile]), Form("hNHits_T2Barrel_%d", energies[iFile]), 100, -0.5, 100-0.5);
      hNHits_T1Endcap[iFile] = new TH1F (Form("hNHits_T1Endcap_%d", energies[iFile]), Form("hNHits_T1Endcap_%d", energies[iFile]), 100, -0.5, 100-0.5);
      hNHits_T2Endcap[iFile] = new TH1F (Form("hNHits_T2Endcap_%d", energies[iFile]), Form("hNHits_T2Endcap_%d", energies[iFile]), 100, -0.5, 100-0.5);
      
      hDeltaPhi_Tower[iFile] = new TH1F (Form("hDeltaPhi_Tower_%d", energies[iFile]), Form("hDeltaPhi_Tower_%d", energies[iFile]), NBINS, -0.3*scale, 0.3*scale);
      hDeltaEta_Tower[iFile] = new TH1F (Form("hDeltaEta_Tower_%d", energies[iFile]), Form("hDeltaEta_Tower_%d", energies[iFile]), NBINS, -0.3*scale, 0.3*scale);
      hDeltaTheta_Tower[iFile] = new TH1F (Form("hDeltaTheta_Tower_%d", energies[iFile]), Form("hDeltaTheta_Tower_%d", energies[iFile]), NBINS, -0.3*scale, 0.3*scale);
      
      hDeltaPhi_Barrel[iFile] = new TH1F (Form("hDeltaPhi_Barrel_%d", energies[iFile]), Form("hDeltaPhi_Barrel_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_Barrel[iFile] = new TH1F (Form("hDeltaEta_Barrel_%d", energies[iFile]), Form("hDeltaEta_Barrel_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaTheta_Barrel[iFile] = new TH1F (Form("hDeltaTheta_Barrel_%d", energies[iFile]), Form("hDeltaTheta_Barrel_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaPhi_Endcap[iFile] = new TH1F (Form("hDeltaPhi_Endcap_%d", energies[iFile]), Form("hDeltaPhi_Endcap_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_Endcap[iFile] = new TH1F (Form("hDeltaEta_Endcap_%d", energies[iFile]), Form("hDeltaEta_Endcap_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaTheta_Endcap[iFile] = new TH1F (Form("hDeltaTheta_Endcap_%d", energies[iFile]), Form("hDeltaTheta_Endcap_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      
      hDeltaPhi_Barrel_wGravity[iFile] = new TH1F (Form("hDeltaPhi_Barrel_wGravity_%d", energies[iFile]), Form("hDeltaPhi_Barrel_wGravity_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_Barrel_wGravity[iFile] = new TH1F (Form("hDeltaEta_Barrel_wGravity_%d", energies[iFile]), Form("hDeltaEta_Barrel_wGravity_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaTheta_Barrel_wGravity[iFile] = new TH1F (Form("hDeltaTheta_Barrel_wGravity_%d", energies[iFile]), Form("hDeltaTheta_Barrel_wGravity_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaPhi_Endcap_wGravity[iFile] = new TH1F (Form("hDeltaPhi_Endcap_wGravity_%d", energies[iFile]), Form("hDeltaPhi_Endcap_wGravity_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_Endcap_wGravity[iFile] = new TH1F (Form("hDeltaEta_Endcap_wGravity_%d", energies[iFile]), Form("hDeltaEta_Endcap_wGravity_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaTheta_Endcap_wGravity[iFile] = new TH1F (Form("hDeltaTheta_Endcap_wGravity_%d", energies[iFile]), Form("hDeltaTheta_Endcap_wGravity_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      
      hDeltaPhi_TimingBarrel_Comb[iFile] = new TH1F (Form("hDeltaPhi_TimingBarrel_Comb_%d", energies[iFile]), Form("hDeltaPhi_TimingBarrel_Comb_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_TimingBarrel_Comb[iFile] = new TH1F (Form("hDeltaEta_TimingBarrel_Comb_%d", energies[iFile]), Form("hDeltaEta_TimingBarrel_Comb_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingBarrel_Comb[iFile] = new TH1F (Form("hDeltaTheta_TimingBarrel_Comb_%d", energies[iFile]), Form("hDeltaTheta_TimingBarrel_Comb_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaPhi_TimingEndcap_Comb[iFile] = new TH1F (Form("hDeltaPhi_TimingEndcap_Comb_%d", energies[iFile]), Form("hDeltaPhi_TimingEndcap_Comb_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_TimingEndcap_Comb[iFile] = new TH1F (Form("hDeltaEta_TimingEndcap_Comb_%d", energies[iFile]), Form("hDeltaEta_TimingEndcap_Comb_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingEndcap_Comb[iFile] = new TH1F (Form("hDeltaTheta_TimingEndcap_Comb_%d", energies[iFile]), Form("hDeltaTheta_TimingEndcap_Comb_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      
      hDeltaPhi_TimingBarrel_F[iFile] = new TH1F (Form("hDeltaPhi_TimingBarrel_F_%d", energies[iFile]), Form("hDeltaPhi_TimingBarrel_F_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_TimingBarrel_F[iFile] = new TH1F (Form("hDeltaEta_TimingBarrel_F_%d", energies[iFile]), Form("hDeltaEta_TimingBarrel_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingBarrel_F[iFile] = new TH1F (Form("hDeltaTheta_TimingBarrel_F_%d", energies[iFile]), Form("hDeltaTheta_TimingBarrel_F_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaPhi_TimingEndcap_F[iFile] = new TH1F (Form("hDeltaPhi_TimingEndcap_F_%d", energies[iFile]), Form("hDeltaPhi_TimingEndcap_F_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_TimingEndcap_F[iFile] = new TH1F (Form("hDeltaEta_TimingEndcap_F_%d", energies[iFile]), Form("hDeltaEta_TimingEndcap_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingEndcap_F[iFile] = new TH1F (Form("hDeltaTheta_TimingEndcap_F_%d", energies[iFile]), Form("hDeltaTheta_TimingEndcap_F_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      
      hDeltaPhi_TimingBarrel_R[iFile] = new TH1F (Form("hDeltaPhi_TimingBarrel_R_%d", energies[iFile]), Form("hDeltaPhi_TimingBarrel_R_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_TimingBarrel_R[iFile] = new TH1F (Form("hDeltaEta_TimingBarrel_R_%d", energies[iFile]), Form("hDeltaEta_TimingBarrel_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingBarrel_R[iFile] = new TH1F (Form("hDeltaTheta_TimingBarrel_R_%d", energies[iFile]), Form("hDeltaTheta_TimingBarrel_R_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaPhi_TimingEndcap_R[iFile] = new TH1F (Form("hDeltaPhi_TimingEndcap_R_%d", energies[iFile]), Form("hDeltaPhi_TimingEndcap_R_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
      hDeltaEta_TimingEndcap_R[iFile] = new TH1F (Form("hDeltaEta_TimingEndcap_R_%d", energies[iFile]), Form("hDeltaEta_TimingEndcap_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingEndcap_R[iFile] = new TH1F (Form("hDeltaTheta_TimingEndcap_R_%d", energies[iFile]), Form("hDeltaTheta_TimingEndcap_R_%d", energies[iFile]), NBINS, -0.1*scale, 0.1*scale);
                  
      hDeltaPhi_TimingBarrel_W_F[iFile] = new TH1F (Form("hDeltaPhi_TimingBarrel_W_F_%d", energies[iFile]), Form("hDeltaPhi_TimingBarrel_W_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaEta_TimingBarrel_W_F[iFile] = new TH1F (Form("hDeltaEta_TimingBarrel_W_F_%d", energies[iFile]), Form("hDeltaEta_TimingBarrel_W_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingBarrel_W_F[iFile] = new TH1F (Form("hDeltaTheta_TimingBarrel_W_F_%d", energies[iFile]), Form("hDeltaTheta_TimingBarrel_W_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaPhi_TimingEndcap_W_F[iFile] = new TH1F (Form("hDeltaPhi_TimingEndcap_W_F_%d", energies[iFile]), Form("hDeltaPhi_TimingEndcap_W_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaEta_TimingEndcap_W_F[iFile] = new TH1F (Form("hDeltaEta_TimingEndcap_W_F_%d", energies[iFile]), Form("hDeltaEta_TimingEndcap_W_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingEndcap_W_F[iFile] = new TH1F (Form("hDeltaTheta_TimingEndcap_W_F_%d", energies[iFile]), Form("hDeltaTheta_TimingEndcap_W_F_%d", energies[iFile]), NBINS, -0.1, 0.1);
      
      hDeltaPhi_TimingBarrel_W_R[iFile] = new TH1F (Form("hDeltaPhi_TimingBarrel_W_R_%d", energies[iFile]), Form("hDeltaPhi_TimingBarrel_W_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaEta_TimingBarrel_W_R[iFile] = new TH1F (Form("hDeltaEta_TimingBarrel_W_R_%d", energies[iFile]), Form("hDeltaEta_TimingBarrel_W_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingBarrel_W_R[iFile] = new TH1F (Form("hDeltaTheta_TimingBarrel_W_R_%d", energies[iFile]), Form("hDeltaTheta_TimingBarrel_W_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaPhi_TimingEndcap_W_R[iFile] = new TH1F (Form("hDeltaPhi_TimingEndcap_W_R_%d", energies[iFile]), Form("hDeltaPhi_TimingEndcap_W_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaEta_TimingEndcap_W_R[iFile] = new TH1F (Form("hDeltaEta_TimingEndcap_W_R_%d", energies[iFile]), Form("hDeltaEta_TimingEndcap_W_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      hDeltaTheta_TimingEndcap_W_R[iFile] = new TH1F (Form("hDeltaTheta_TimingEndcap_W_R_%d", energies[iFile]), Form("hDeltaTheta_TimingEndcap_W_R_%d", energies[iFile]), NBINS, -0.1, 0.1);
      
      hTOF_TimingBarrel[iFile] = new TH1F (Form("hTOF_TimingBarrel_%d", energies[iFile]), Form("hTOF_TimingBarrel_%d", energies[iFile]), 1000, -100, 100);
      hTOF_TimingEndcap[iFile] = new TH1F (Form("hTOF_TimingEndcap_%d", energies[iFile]), Form("hTOF_TimingEndcap_%d", energies[iFile]), 1000, -100, 100);
      
      hScatterTOF_TimingBarrel[iFile] = new TH2F (Form("hScatterTOF_TimingBarrel_%d", energies[iFile]), Form("hScatterTOF_TimingBarrel_%d", energies[iFile]), 200, 1500, 2500, 1000, 4, 10);
      hScatterTOF_TimingEndcap[iFile] = new TH2F (Form("hScatterTOF_TimingEndcap_%d", energies[iFile]), Form("hScatterTOF_TimingEndcap_%d", energies[iFile]), 200, 1500, 2500, 1000, 4, 10);
  }      
  
  
  
  //run over energy scan
  TFile * RunFile[NFILES];  
  
  SCEPCal_GeometryHelper myGeometry;
  
  int ENE_BINS = 8;  
  float maxEneRange = 100;
  float minEneRange = 5;
  float vMinEneBin [ENE_BINS];
  float vMaxEneBin [ENE_BINS];
  std::vector<double> *vDeltaThetaECAL[ENE_BINS];
  
  for (int iEne = 0; iEne<ENE_BINS; iEne++)
  {
      float stepEne = (maxEneRange-minEneRange)/ENE_BINS;
      float min_ene = minEneRange+stepEne*iEne;
      float max_ene = min_ene+stepEne;
      float ave_ene = (max_ene+min_ene)/2.;
      
      vMinEneBin[iEne]=min_ene;
      vMaxEneBin[iEne]=max_ene;
      
      vDeltaThetaECAL[iEne]     = new std::vector<double>();
  }
  
  int THETA_BINS = 16;  
  float maxThetaRange = 3.14/4*0.9;
  float minThetaRange = 0;
  float vMinThetaBin [THETA_BINS];
  float vMaxThetaBin [THETA_BINS];
  std::vector<double> *vDeltaThetaTiming[THETA_BINS];
  
  for (int iTheta = 0; iTheta<THETA_BINS; iTheta++)
  {
      float stepTheta = (maxThetaRange-minThetaRange)/THETA_BINS;
      float min_theta = minThetaRange+stepTheta*iTheta;
      float max_theta = min_theta+stepTheta;
      float ave_theta = (max_theta+min_theta)/2.;
      
      vMinThetaBin[iTheta]=min_theta;
      vMaxThetaBin[iTheta]=max_theta;
      
      vDeltaThetaTiming[iTheta]     = new std::vector<double>();
  }
  
  
  
  
  for (int iFile = 0; iFile < NFILES; iFile++)
  {
      
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_mu-_50GeV.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_ele_10GeV.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_gamma_10GeV_v4.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_mu_100GeV_timing.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_mu_100GeV_timing_v2.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_mu_100GeV_timing_v3.root","READ"); 
      RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_mu_100GeV_T+E.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/prod/output_SCEPCal_Iso_e-_40GeV.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/prod/output_SCEPCal_fixedPos_gamma_Iso+Uniform1-100_GeV.root","READ");      
//       RunFile[iFile] = new TFile("../root_files/prod/output_SCEPCal_gamma_Iso+Uniform+noMagenticField.root","READ");
      
//       RunFile[iFile] = new TFile("../root_files/prod/output_SCEPCal_fixedPos_e-_Iso+Uniform1-100_GeV.root","READ");      
      
      TTree* TreeRun = (TTree*) RunFile[iFile]->Get("B4");
      myG4TreeVars myTV;
      InitG4Tree (TreeRun, myTV);
      
              
      TreeRun->SetBranchStatus("*", 0);      
      TreeRun->SetBranchStatus("PrimaryParticleEnergy", 1);  
      TreeRun->SetBranchStatus("PrimaryParticleMomentum", 1);  
      
      TreeRun->SetBranchStatus("VecHit_CrystalID", 1);
      TreeRun->SetBranchStatus("VecHit_ScepEneDepF", 1);
      TreeRun->SetBranchStatus("VecHit_ScepEneDepR", 1);
      TreeRun->SetBranchStatus("VecHit_ScepCherF", 1);
      TreeRun->SetBranchStatus("VecHit_ScepCherR", 1);   
      
      
      TreeRun->SetBranchStatus("VecHit_Timing_CrystalID_F", 1); 
      TreeRun->SetBranchStatus("VecHit_Timing_CrystalID_R", 1); 
      TreeRun->SetBranchStatus("VecHit_Timing_ScepEneDepF", 1); 
      TreeRun->SetBranchStatus("VecHit_Timing_ScepEneDepR", 1); 
      TreeRun->SetBranchStatus("VecHit_Timing_ScepTimeF", 1); 
      TreeRun->SetBranchStatus("VecHit_Timing_ScepTimeR", 1); 
 

      ///*******************************************///
      ///		 Run over events	    ///
      ///*******************************************///
    
      int NEVENTS = TreeRun->GetEntries();
      std::cout << "NEVENTS = " << NEVENTS << std::endl;

//       NEVENTS = 10000;
      
      for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
      {
          
          TreeRun->GetEntry(iEvt);
          if (iEvt%1000==0) std::cout << "processing event: " << iEvt << "\r" << std::flush;
        
          
          double px  = myTV.PrimaryParticleMomentum->at(0);
          double py  = myTV.PrimaryParticleMomentum->at(1);
          double pz  = myTV.PrimaryParticleMomentum->at(2);
          double P   = sqrt(px*px+py*py+pz*pz);
          px/= P;
          py/= P;
          pz/= P;
          
//           std::cout << " P = " << P << std::endl;
//           if (P/1000 <70) continue;
          
          double phi   = atan(py/px);
          if (px<0. && py <0.)   {phi = phi - M_PI;}
          if (px<0. && py >0.)   {phi = M_PI + phi;}
          double eta   = -atanh(pz);
          
          double theta = 2*atan(exp(-eta));
          
          //**************************************************************//
          //                             HCAL
          //**************************************************************//
          
          // find hit with max energy (seed)
          int tower_maxHit_index_L;
          int tower_maxHit_index_R;
          double tower_phi_seed = 0.;
          double tower_eta_seed = 0.;
          double tower_theta_seed = 0.;
                    
          std::vector<double>::iterator max_ene_TL;
          std::vector<double>::iterator max_ene_TR;
          
          if (myTV.VectorL->size()>0 || myTV.VectorR->size()>0 )
          {
              float this_ene_L = -999;
              float this_ene_R = -999;
              if (myTV.VectorL->size()>0) 
              {                
                  max_ene_TL = std::max_element(myTV.VectorL->begin(), myTV.VectorL->end());
                  tower_maxHit_index_L = std::distance(myTV.VectorL->begin(), max_ene_TL);
                  this_ene_L = myTV.VectorL->at(tower_maxHit_index_L);  
              }              
              if (myTV.VectorR->size()>0) 
              {
                  max_ene_TR = std::max_element(myTV.VectorR->begin(), myTV.VectorR->end());
                  tower_maxHit_index_R = std::distance(myTV.VectorR->begin(), max_ene_TR);
                  this_ene_R = myTV.VectorR->at(tower_maxHit_index_R);
              }
              
              TVector3 seed_vec;              
              
              
              if (this_ene_L>=this_ene_R) 
              {                                                    
                  seed_vec =  myGeometry.GetTowerVec(tower_maxHit_index_L, 'l');
              }
              else
              {
                  seed_vec =  myGeometry.GetTowerVec(tower_maxHit_index_R, 'r');
              }
                                                                                                                
              tower_phi_seed = seed_vec.Phi();
              tower_eta_seed = seed_vec.Eta();
              tower_theta_seed = seed_vec.Theta();
          }
          
          if (!(tower_phi_seed ==0 && tower_eta_seed ==0 && tower_theta_seed ==0 ))
          {              
              hDeltaPhi_Tower[iFile]->Fill(tower_phi_seed-phi);
              hDeltaEta_Tower[iFile]->Fill(tower_eta_seed-eta);                                                                          
              hDeltaTheta_Tower[iFile]->Fill(tower_theta_seed-theta);            
          }
          
          
          
          //**************************************************************//
          //                             ECAL
          //**************************************************************//
          
          // find hit with max energy (seed)
          int maxHit_index = -1;                
          int seed_crystal_ID = 0;
          double phi_seed = 0.;
          double eta_seed = 0.;
          double theta_seed = 0.;
          
          std::vector<double>::iterator max_ene;
          if (myTV.VecHit_CrystalID->size()>0)
          {
              max_ene = std::max_element(myTV.VecHit_ScepEneDepR->begin(), myTV.VecHit_ScepEneDepR->end());
              maxHit_index = std::distance(myTV.VecHit_ScepEneDepR->begin(), max_ene);               
              seed_crystal_ID = myTV.VecHit_CrystalID->at(maxHit_index);
              TVector3 seed_vec =  myGeometry.GetCrystalVec(seed_crystal_ID);
              phi_seed = seed_vec.Phi();
              eta_seed = seed_vec.Eta();
              theta_seed = seed_vec.Theta();
          }
                                        
          //find center of gravity for phi and eta
          double phi_weighed = 0.;
          double eta_weighed = 0.;
          double theta_weighed = 0.;
          double w_tot = 0.;
          double w_tot_phi = 0.;
          
          
          for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
          {                            
              TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
              double this_phi = this_vec.Phi();
              double this_eta = this_vec.Eta();
              double this_theta = this_vec.Theta();
              double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;              
              double maxHitEne = (myTV.VecHit_ScepEneDepF->at(maxHit_index)+myTV.VecHit_ScepEneDepR->at(maxHit_index))/1000. ;
              
//               if (this_ene < maxHitEne*0.01) continue;
//               if (this_ene < 0.1) continue;
//               double this_ene = (myTV.VecHit_ScepEneDepF->at(i))/1000.;
              
              phi_weighed += this_phi*this_ene;
//               phi_weighed += this_phi*pow(this_ene,2);
              eta_weighed += this_eta*this_ene;
              theta_weighed += this_theta*this_ene;
              w_tot += this_ene;     
              w_tot_phi+=this_ene;
//               w_tot_phi+=pow(this_ene,2);
//               std::cout << "iEvt = " << iEvt << " :: ene_hit[" << i << "] = " << myTV.VecHit_ScepEneDepR->at(i)/1000. << " :: while max_ene is at " << maxHit_index << std::endl;              
          }
          phi_weighed/=w_tot_phi;
          eta_weighed/=w_tot;
          theta_weighed/=w_tot;
          
          if (maxHit_index!= -1)
          {
            double max_hit_ene = (myTV.VecHit_ScepEneDepF->at(maxHit_index)+myTV.VecHit_ScepEneDepR->at(maxHit_index))/1000.;
                                  
            if (fabs(seed_crystal_ID)<1000000  && max_hit_ene>0.2)
            {              
              //                 std::cout << "event (" << iEvt << "), hit in barrel (" << max_hit_ene  << ") --> this_ch_vec[" << seed_crystal_ID << "]: phi = " << phi_seed << " :: theta = " << this_ch_vec.Theta() << " :: Eta = " << eta_seed << std::endl;
              //                 std::cout << " input direction: phi = " << phi << " :: theta = "  << theta << " :: eta = " << eta << std::endl;                                
                hDeltaPhi_Barrel[iFile]->Fill(phi_seed-phi);
                hDeltaEta_Barrel[iFile]->Fill(eta_seed-eta);                                                            
                hDeltaTheta_Barrel[iFile]->Fill(theta_seed-theta);            
                hDeltaPhi_Barrel_wGravity[iFile]->Fill(phi_weighed-phi);
                hDeltaEta_Barrel_wGravity[iFile]->Fill(eta_weighed-eta);     
                hDeltaTheta_Barrel_wGravity[iFile]->Fill(theta_weighed-theta);     
                
//                 if (fabs(theta_seed-theta)>0.2) std::cout << "attention for ECAL barrel! :: theta_seed = " << theta_seed << " :: theta = " << theta 
//                                                           << " :: deltaTheta = " << theta_seed-theta 
//                                                           << " :: seed_crystal_ID = " << seed_crystal_ID
//                                                           << std::endl;
                                
                for (int iEne = 0; iEne<ENE_BINS; iEne++)
                {                
                    if (P/1000>=vMinEneBin[iEne] && P/1000<vMaxEneBin[iEne]) 
                    {
                        vDeltaThetaECAL[iEne]->push_back(theta_weighed-theta);
                    }
                }
            }            
            else if (fabs(seed_crystal_ID)>1000000  && max_hit_ene>0.2)
            {              
//                 std::cout << "event (" << iEvt << "), hit in endcap (" << max_hit_ene  << ") --> this_ch_vec[" << seed_crystal_ID << "]: phi = " << phi_seed << " :: Eta = "<< eta_seed << std::endl;
//                 std::cout << " input direction: phi = " << phi << " :: eta = " << eta << std::endl;                  

                hDeltaPhi_Endcap[iFile]->Fill(phi_seed-phi);
                hDeltaEta_Endcap[iFile]->Fill(eta_seed-eta);
                hDeltaTheta_Endcap[iFile]->Fill(theta_seed-theta);
                hDeltaPhi_Endcap_wGravity[iFile]->Fill(phi_weighed-phi);
                hDeltaEta_Endcap_wGravity[iFile]->Fill(eta_weighed-eta);
                hDeltaTheta_Endcap_wGravity[iFile]->Fill(theta_weighed-theta);
                
                
                if (fabs(phi_seed-phi)>0.2) std::cout << "attention for ECAL endcap! :: phi_seed = " << phi_seed << " :: phi = " << phi 
                                                          << " :: deltaPhi = " << phi_seed-phi 
                                                          << " :: seed_crystal_ID = " << seed_crystal_ID
                                                          << std::endl;
            }
          }
                 
                 
        
          
          //**************************************************************//
          //                            Timing
          //**************************************************************//
          
          
          
          double minEneHitTiming = 1;   //MeV
          double maxEneHitTiming = 15;
          long unsigned int maxTimingHits      = 10;
          
          double c_speed = 1./299792458*1e9; //mm per picosecond
          
          double time_acceptance = 3; //time window to reject out of time hits
          
          //**********//
          //front layer
          int    T1_maxHit_index = -1;                
          int    T1_seed_crystal_ID = 0;
          double T1_phi_seed = 0.;
          double T1_eta_seed = 0.;
          double T1_theta_seed = 0.;
          double T1_maxHit_ene = 0;
          double T1_maxHit_time = 0;
          double T1_distance = 0;
          
//           std::cout << "getting here? with vec size = " << myTV.VecHit_Timing_CrystalID_F->size() << std::endl;
          //get first coordinate (phi?) from front layer
          if (myTV.VecHit_Timing_CrystalID_F->size()>0 && myTV.VecHit_Timing_CrystalID_F->size()<maxTimingHits )
          {
              max_ene         = std::max_element(myTV.VecHit_Timing_ScepEneDepF->begin(), myTV.VecHit_Timing_ScepEneDepF->end());
              T1_maxHit_index = std::distance(myTV.VecHit_Timing_ScepEneDepF->begin(), max_ene);
              T1_seed_crystal_ID = myTV.VecHit_Timing_CrystalID_F->at(T1_maxHit_index);
              T1_maxHit_ene      = myTV.VecHit_Timing_ScepEneDepF->at(T1_maxHit_index);              
              T1_maxHit_time     = myTV.VecHit_Timing_ScepTimeF->at(T1_maxHit_index);              
                                          
              TVector3 seed_vec =  myGeometry.GetCrystalTimingVec(T1_seed_crystal_ID, 1);
              T1_phi_seed = seed_vec.Phi();
              T1_eta_seed = seed_vec.Eta();                            
              T1_theta_seed = seed_vec.Theta();          
              T1_distance =  sqrt(pow(seed_vec.X(),2) + pow(seed_vec.Y(),2) + pow(seed_vec.Z(),2));              
//               std::cout << "phi = " << phi << " :: T1_phi_seed = " << T1_phi_seed << " :: deltaPhi = " << T1_phi_seed-phi << " :: eta = " << eta << " :: T1_eta_seed = " << eta_seed << std::endl;
          }
          
          //find center of gravity for phi and eta
          double T1_phi_weighed = 0.;
          double T1_eta_weighed = 0.;
          double T1_theta_weighed = 0.;
          double T1_w_tot = 0.;
          
          for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_F->size(); i++)
          {                            
              TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_F->at(i),1);
              double this_phi = this_vec.Phi();
              double this_eta = this_vec.Eta();
              double this_theta = this_vec.Theta();
              double this_ene = myTV.VecHit_Timing_ScepEneDepF->at(i);
              
              if (this_ene > minEneHitTiming && this_ene < maxEneHitTiming)
              {
                  T1_phi_weighed += this_phi*this_ene;
                  T1_eta_weighed += this_eta*this_ene;
                  T1_theta_weighed += this_theta*this_ene;
                  T1_w_tot += this_ene;              
              }
          }
          T1_phi_weighed/=T1_w_tot;
          T1_eta_weighed/=T1_w_tot;
          T1_theta_weighed/=T1_w_tot;

          
          
          
          //**********//
          //rear layer
          int    T2_maxHit_index = -1;                
          int    T2_seed_crystal_ID = 0;
          double T2_phi_seed = 0.;
          double T2_eta_seed = 0.;
          double T2_theta_seed = 0.;
//           double T2_maxHit_ene = 0;
//           double T2_maxHit_time = 0;
          
          if (myTV.VecHit_Timing_CrystalID_R->size()>0 && myTV.VecHit_Timing_CrystalID_R->size()<maxTimingHits)
          {
              max_ene         = std::max_element(myTV.VecHit_Timing_ScepEneDepR->begin(), myTV.VecHit_Timing_ScepEneDepR->end());
              T2_maxHit_index = std::distance(myTV.VecHit_Timing_ScepEneDepR->begin(), max_ene); 
              T2_seed_crystal_ID = myTV.VecHit_Timing_CrystalID_R->at(T2_maxHit_index);
//               T2_maxHit_ene      = myTV.VecHit_Timing_ScepEneDepR->at(T2_maxHit_index);
//               T2_maxHit_time     = myTV.VecHit_Timing_ScepTimeR->at(T2_maxHit_index);              
              
              TVector3 seed_vec =  myGeometry.GetCrystalTimingVec(T2_seed_crystal_ID, 2);
              T2_phi_seed = seed_vec.Phi();
              T2_eta_seed = seed_vec.Eta();                                          
              T2_theta_seed = seed_vec.Theta();
//               std::cout << "phi = " << phi << " :: T2_phi_seed = " << T2_phi_seed << " :: deltaPhi = " << T2_phi_seed-phi << " :: eta = " << eta << " :: T2_eta_seed = " << T2_eta_seed << std::endl;
          }
           
          //find center of gravity for phi and eta
          double T2_phi_weighed = 0.;
          double T2_eta_weighed = 0.;
          double T2_theta_weighed = 0.;
          double T2_w_tot = 0.;
          
          for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_R->size(); i++)
          {                            
              TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_R->at(i),2);
              double this_phi = this_vec.Phi();
              double this_eta = this_vec.Eta();
              double this_theta = this_vec.Theta();
              double this_ene = myTV.VecHit_Timing_ScepEneDepR->at(i);
//               std::cout << "this_ene = " << this_ene << std::endl;
              if (this_ene > minEneHitTiming && this_ene < maxEneHitTiming)
              {
                  T2_phi_weighed += this_phi*this_ene;
                  T2_eta_weighed += this_eta*this_ene;
                  T2_theta_weighed += this_theta*this_ene;
                  T2_w_tot += this_ene;              
              }
          }
          T2_phi_weighed/=T2_w_tot;
          T2_eta_weighed/=T2_w_tot;
          T2_theta_weighed/=T2_w_tot;
          
          
          //**********//
          //combined layers                    
          double TT_phi_seed = 0.;
          double TT_eta_seed = 0.;          
          double TT_theta_seed = 0.;      
          double TT_maxHit_time = 0;
          double TT_maxHit_distance = 0;
          
          if (myTV.VecHit_Timing_CrystalID_F->size()>0 && myTV.VecHit_Timing_CrystalID_R->size()>0 && myTV.VecHit_Timing_CrystalID_F->size()<maxTimingHits)
          {
              max_ene         = std::max_element(myTV.VecHit_Timing_ScepEneDepF->begin(), myTV.VecHit_Timing_ScepEneDepF->end());
              int T1_maxHit_index_temp = std::distance(myTV.VecHit_Timing_ScepEneDepF->begin(), max_ene);                             
              int T1_seed_crystal_ID_temp = myTV.VecHit_Timing_CrystalID_F->at(T1_maxHit_index_temp);
              
              max_ene         = std::max_element(myTV.VecHit_Timing_ScepEneDepR->begin(), myTV.VecHit_Timing_ScepEneDepR->end());
              int T2_maxHit_index_temp = std::distance(myTV.VecHit_Timing_ScepEneDepR->begin(), max_ene);               
              int T2_seed_crystal_ID_temp = myTV.VecHit_Timing_CrystalID_R->at(T2_maxHit_index_temp);
              
              TT_maxHit_time     = myTV.VecHit_Timing_ScepTimeF->at(T1_maxHit_index_temp)/myTV.VecHit_Timing_ScepEneDepF->at(T1_maxHit_index_temp);              
              
//               TVector3 seed_vec =  myGeometry.GetCrystalTimingBothVec(T1_seed_crystal_ID_temp, T2_seed_crystal_ID_temp);
              TVector3 seed_vec =  myGeometry.GetCrystalTimingBothVec(T1_seed_crystal_ID_temp, T2_seed_crystal_ID_temp);              
              TT_phi_seed = seed_vec.Phi();
              TT_eta_seed = seed_vec.Eta();
              TT_theta_seed = seed_vec.Theta();
              TT_maxHit_distance =  sqrt(pow(seed_vec.X(),2) + pow(seed_vec.Y(),2) + pow(seed_vec.Z(),2));                                          
                                          
          }
              
          
          if (      true
//                    && myTV.VecHit_Timing_CrystalID_F->size()>0 && myTV.VecHit_Timing_CrystalID_F->size()<maxTimingHits
//                 && myTV.VecHit_Timing_CrystalID_R->size()>0 && myTV.VecHit_Timing_CrystalID_R->size()<maxTimingHits
                && fabs(T1_maxHit_time/T1_maxHit_ene - c_speed*T1_distance ) < time_acceptance
                && T2_theta_seed !=0 && T1_theta_seed !=0
               )
          {
              if (fabs(T1_seed_crystal_ID)<1000000)// && T1_seed_crystal_ID<0)
              {
                  hNHits_T1Barrel[iFile]->Fill(myTV.VecHit_Timing_CrystalID_F->size());
                  hNHits_T2Barrel[iFile]->Fill(myTV.VecHit_Timing_CrystalID_R->size());
                  
                  hDeltaPhi_TimingBarrel_Comb[iFile]->Fill(TT_phi_seed-phi);              
                  hDeltaEta_TimingBarrel_Comb[iFile]->Fill(TT_eta_seed-eta);
                  hDeltaTheta_TimingBarrel_Comb[iFile]->Fill(TT_theta_seed-theta);
                                    
                  if (fabs(T2_phi_seed)>3.140750) std::cout << "T2_phi_seed = " << T2_phi_seed << std::endl;
                  
                  hDeltaPhi_TimingBarrel_F[iFile]->Fill(T1_phi_seed-phi);              
                  hDeltaEta_TimingBarrel_F[iFile]->Fill(T1_eta_seed-eta);
                  hDeltaTheta_TimingBarrel_F[iFile]->Fill(T1_theta_seed-theta);
                  
                  hDeltaPhi_TimingBarrel_W_F[iFile]->Fill(T1_phi_weighed-phi);            
                  hDeltaEta_TimingBarrel_W_F[iFile]->Fill(T1_eta_weighed-eta);          
                  hDeltaTheta_TimingBarrel_W_F[iFile]->Fill(T1_theta_weighed-theta);          
//                   std::cout << " theta = " << theta << std::endl;
                  for (int iTheta = 0; iTheta<THETA_BINS; iTheta++)
                  {                
                    if (fabs(theta-3.14/2)>=vMinThetaBin[iTheta] && fabs(theta-3.14/2)<vMaxThetaBin[iTheta]) 
                    {
                        vDeltaThetaTiming[iTheta]->push_back(T1_theta_weighed-theta);
                    }
                  }
                  
                  hDeltaPhi_TimingBarrel_R[iFile]->Fill(T2_phi_seed-phi);
                  hDeltaEta_TimingBarrel_R[iFile]->Fill(T2_eta_seed-eta);
                  hDeltaTheta_TimingBarrel_R[iFile]->Fill(T2_theta_seed-theta);
                  
                  hDeltaPhi_TimingBarrel_W_R[iFile]->Fill(T2_phi_weighed-phi);
                  hDeltaEta_TimingBarrel_W_R[iFile]->Fill(T2_eta_weighed-eta);                                
                  hDeltaTheta_TimingBarrel_W_R[iFile]->Fill(T2_theta_weighed-theta);                                
                  
                  if (myTV.VecHit_Timing_CrystalID_F->size() == 1) 
                  {
                      hScatterTOF_TimingBarrel[iFile]->Fill(TT_maxHit_distance-1.5, TT_maxHit_time/1000);
                      hTOF_TimingBarrel[iFile]->Fill((TT_maxHit_distance-1.5)*c_speed - TT_maxHit_time);
//                       std::cout << " 
                  }
//                   std::cout << "TT_maxHit_distance = " << TT_maxHit_distance << " :: TT_maxHit_time = " << TT_maxHit_time <<  " :: TT_maxHit_distance/c_speed = " << TT_maxHit_distance*c_speed << std::endl;
                
                if (fabs(T1_eta_seed-eta)<0.003)    
                {
//                     std::cout << "T1_eta_seed-eta = " << T1_eta_seed-eta  << " at module ID = " << myGeometry.GetTimingModuleID(myTV.VecHit_Timing_CrystalID_F->at(T1_maxHit_index)) << " and crystal ID = " << myGeometry.GetTimingCrystalID(myTV.VecHit_Timing_CrystalID_F->at(T1_maxHit_index)) << " --> ene_hit = " << T1_maxHit_ene << " :: time_hit = " << T1_maxHit_time/T1_maxHit_ene << " :: distance = " << T1_distance <<  " :: t_diff = " << T1_maxHit_time/T1_maxHit_ene - c_speed*T1_distance  << std::endl;
                }
                
                if (fabs(T1_eta_seed-eta)>0.003)                    
                {
//                     std::cout << "out of peak! --> T1_eta_seed-eta = " << T1_eta_seed-eta  << " at module ID = " << myGeometry.GetTimingModuleID(myTV.VecHit_Timing_CrystalID_F->at(T1_maxHit_index)) << " and crystal ID = " << myGeometry.GetTimingCrystalID(myTV.VecHit_Timing_CrystalID_F->at(T1_maxHit_index)) << " --> ene_hit = " << T1_maxHit_ene << " :: time_hit = " << T1_maxHit_time/T1_maxHit_ene << " :: distance = " << T1_distance << " :: t_diff = " << T1_maxHit_time/T1_maxHit_ene - c_speed*T1_distance  << std::endl;
                }
            }
            
            else if (fabs(T1_seed_crystal_ID)>=1000000)
            {
                hNHits_T1Endcap[iFile]->Fill(myTV.VecHit_Timing_CrystalID_F->size());
                hNHits_T2Endcap[iFile]->Fill(myTV.VecHit_Timing_CrystalID_R->size());
                
//                 if (fabs(T1_phi_seed)>3.14) std::cout << "T1_phi_seed = " << T1_phi_seed << std::endl;
                
                hDeltaPhi_TimingEndcap_Comb[iFile]->Fill(TT_phi_seed-phi);              
                hDeltaEta_TimingEndcap_Comb[iFile]->Fill(TT_eta_seed-eta);
                hDeltaTheta_TimingEndcap_Comb[iFile]->Fill(TT_theta_seed-theta);
                
                hDeltaPhi_TimingEndcap_F[iFile]->Fill(T1_phi_seed-phi);              
                hDeltaEta_TimingEndcap_F[iFile]->Fill(T1_eta_seed-eta);
                hDeltaTheta_TimingEndcap_F[iFile]->Fill(T1_theta_seed-theta);
                
                hDeltaPhi_TimingEndcap_W_F[iFile]->Fill(T1_phi_weighed-phi);            
                hDeltaEta_TimingEndcap_W_F[iFile]->Fill(T1_eta_weighed-eta);          
                hDeltaTheta_TimingEndcap_W_F[iFile]->Fill(T1_theta_weighed-theta);          
                
                hDeltaPhi_TimingEndcap_R[iFile]->Fill(T2_phi_seed-phi);
                hDeltaEta_TimingEndcap_R[iFile]->Fill(T2_eta_seed-eta);
                hDeltaTheta_TimingEndcap_R[iFile]->Fill(T2_theta_seed-theta);
//                 if (T2_theta_seed-theta<-0.05) std::cout << "attention! T2_theta_seed = " << T2_theta_seed << " :: theta = "  << theta << std::endl;
                
                hDeltaPhi_TimingEndcap_W_R[iFile]->Fill(T2_phi_weighed-phi);
                hDeltaEta_TimingEndcap_W_R[iFile]->Fill(T2_eta_weighed-eta);
                hDeltaTheta_TimingEndcap_W_R[iFile]->Fill(T2_theta_weighed-theta);
                
                if (myTV.VecHit_Timing_CrystalID_F->size() == 1) 
                {
                    hScatterTOF_TimingEndcap[iFile]->Fill(TT_maxHit_distance-3, TT_maxHit_time/1000);
                    hTOF_TimingEndcap[iFile]->Fill((TT_maxHit_distance-3)*c_speed - TT_maxHit_time);
                }
                
//                 std::cout << "TT_maxHit_distance = " << TT_maxHit_distance << " :: TT_maxHit_time = " << TT_maxHit_time <<  " :: TT_maxHit_distance/c_speed = " << TT_maxHit_distance*c_speed << std::endl;
//                 std::cout << "T1_eta_seed-eta = " << T1_eta_seed-eta  << " :: T1_phi_seed - phi = " << T1_phi_seed - phi << std::endl;
            }              
          }          
//           std::cout << "phi = " << phi << " :: T2_phi_weighed = " << T2_phi_weighed << " :: deltaPhi = " << T2_phi_seed-phi << " :: eta = " << eta << " :: T2_eta_weighed = " << T2_eta_weighed << std::endl;          
      }
  }
  
  
  
  
  
  // PLOTTING
  std::cout << " plotting ...  " << std::endl;
  
  float tower_range = 0.3;
  TCanvas * cDeltaPhiTower = new TCanvas ("cDeltaPhiTower", "cDeltaPhiTower", 600, 500);
  cDeltaPhiTower->cd();
  hDeltaPhi_Tower[NFILES-1]->Draw();
  hDeltaPhi_Tower[NFILES-1]->SetStats(0);
  hDeltaPhi_Tower[NFILES-1]->SetTitle(Form("DR Tower: %d GeV %s",energies[0], particle_name.c_str()));
  hDeltaPhi_Tower[NFILES-1]->GetXaxis()->SetRangeUser(-tower_range, tower_range);
  hDeltaPhi_Tower[NFILES-1]->GetYaxis()->SetRangeUser(0,hDeltaPhi_Tower[NFILES-1]->GetMaximum()*1.5);
  hDeltaPhi_Tower[NFILES-1]->GetXaxis()->SetTitle("#phi_{reco} - #phi_{truth@VTX} [rad]");
  hDeltaPhi_Tower[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaPhi_Tower[NFILES-1]->SetLineColor(kBlack);
  
//   leg = new TLegend(0.15,0.68,0.58,0.88,NULL,"brNDC");

//   for (int iFile = NFILES-1; iFile>=0; iFile--)
//   {
            
/*      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -tower_range, tower_range);
      fitGaus->SetLineColor(kBlack);
      hDeltaPhi_Tower[iFile]->Fit(fitGaus, "QR");
      std::cout << "---- BARREL ----" << std::endl;      
      std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      float res = fitGaus->GetParameter(2);
      
      leg->AddEntry(hDeltaPhi_Tower[iFile],          Form("Seed:~#sigma_{#phi}=%.2f mrad",res*1000), "lp");         */ 
//   }
//   leg->Draw();
  if (SAVEPLOTS) cDeltaPhiTower->SaveAs("plots/cDeltaPhi_Tower.png");
  
  
  
  float fit_range = 0.01;
  float phi_res_b, phi_res_e, phi_res_b_cg, phi_res_e_cg;
  float eta_res_b, eta_res_e, eta_res_b_cg, eta_res_e_cg;
  TCanvas * cDeltaPhi = new TCanvas ("cDeltaPhi", "cDeltaPhi", 600, 500);
  cDeltaPhi->cd();
  hDeltaPhi_Barrel[NFILES-1]->Draw();
  hDeltaPhi_Barrel[NFILES-1]->SetStats(0);
  hDeltaPhi_Barrel[NFILES-1]->SetTitle(Form("Barrel: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaPhi_Barrel[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range, fit_range);
  hDeltaPhi_Barrel[NFILES-1]->GetYaxis()->SetRangeUser(1,hDeltaPhi_Barrel_wGravity[NFILES-1]->GetMaximum()*5);
  hDeltaPhi_Barrel[NFILES-1]->GetXaxis()->SetTitle("#phi_{reco} - #phi_{truth@VTX} [rad]");
  hDeltaPhi_Barrel[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaPhi_Barrel[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaPhi_Barrel_wGravity[NFILES-1]->SetLineColor(kGreen+1);
  hDeltaPhi_Barrel_wGravity[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.58,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
            
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range, fit_range);
      fitGaus->SetLineColor(kBlack);
      hDeltaPhi_Barrel[iFile]->Fit(fitGaus, "QR");
      std::cout << "---- BARREL ----" << std::endl;      
      std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_b = fitGaus->GetParameter(2);
//       phi_res_b = fitGaus->GetParameter(hDeltaPhi_Barrel[iFile]->GetRMS());
      
      fitGaus->SetLineColor(kGreen);
      hDeltaPhi_Barrel_wGravity[iFile]->Fit(fitGaus, "QR");
      std::cout << "WEIGHTED: eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_b_cg = fitGaus->GetParameter(2);      
//       phi_res_b_cg = hDeltaPhi_Barrel_wGravity[iFile]->GetRMS();      

//       leg->AddEntry(hDeltaPhi_Barrel[iFile],          Form("Seed:~RMS_{#phi}=%.2f mrad",phi_res_b*1000), "lp");          
//       leg->AddEntry(hDeltaPhi_Barrel_wGravity[iFile], Form("CoG:~~RMS_{#phi}=%.2f mrad",phi_res_b_cg*1000),  "lp");        
      
      leg->AddEntry(hDeltaPhi_Barrel[iFile],          Form("Seed:~#sigma_{#phi}=%.2f mrad",phi_res_b*1000), "lp");          
      leg->AddEntry(hDeltaPhi_Barrel_wGravity[iFile], Form("CoG:~~#sigma_{#phi}=%.2f mrad",phi_res_b_cg*1000),  "lp");        
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaPhi->SaveAs("plots/cDeltaPhi_Barrel.png");
  
  
  
  TCanvas * cDeltaPhi_Endcap = new TCanvas ("cDeltaPhi_Endcap", "cDeltaPhi_Endcap", 600, 500);
  cDeltaPhi_Endcap->cd();
  hDeltaPhi_Endcap[NFILES-1]->Draw();
  hDeltaPhi_Endcap[NFILES-1]->SetStats(0);
  hDeltaPhi_Endcap[NFILES-1]->SetTitle(Form("Endcaps: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaPhi_Endcap[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range, fit_range);
  hDeltaPhi_Endcap[NFILES-1]->GetYaxis()->SetRangeUser(1,hDeltaPhi_Endcap_wGravity[NFILES-1]->GetMaximum()*5);
  hDeltaPhi_Endcap[NFILES-1]->GetXaxis()->SetTitle("#phi_{reco} - #phi_{truth@VTX} [rad]");
  hDeltaPhi_Endcap[NFILES-1]->GetYaxis()->SetTitle("Counts");  
  hDeltaPhi_Endcap[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaPhi_Endcap_wGravity[NFILES-1]->SetLineColor(kGreen+1);
  hDeltaPhi_Endcap_wGravity[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.58,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {

      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range, fit_range);
      std::cout << "---- ENDCAP ----" << std::endl;
      fitGaus->SetLineColor(kBlack);
      hDeltaPhi_Endcap[iFile]->Fit(fitGaus, "QR");
      std::cout << "eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_e = fitGaus->GetParameter(2);
//       phi_res_e = hDeltaPhi_Endcap[iFile]->GetRMS();      
      
      fitGaus->SetLineColor(kGreen);
      hDeltaPhi_Endcap_wGravity[iFile]->Fit(fitGaus, "QR");
      std::cout << "WEIGHTED: eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_e_cg = fitGaus->GetParameter(2);
//       phi_res_e_cg = hDeltaPhi_Endcap_wGravity[iFile]->GetRMS();      
          
//       leg->AddEntry(hDeltaPhi_Endcap[iFile],          Form("Seed:~RMS_{#phi}=%.2f mrad",phi_res_e*1000), "lp");  
//       leg->AddEntry(hDeltaPhi_Endcap_wGravity[iFile], Form("CoG:~~RMS_{#phi}=%.2f mrad",phi_res_e_cg*1000),  "lp");  
            
      leg->AddEntry(hDeltaPhi_Endcap[iFile],          Form("Seed:~#sigma_{#phi}=%.2f mrad",phi_res_e*1000), "lp");  
      leg->AddEntry(hDeltaPhi_Endcap_wGravity[iFile], Form("CoG:~~#sigma_{#phi}=%.2f mrad",phi_res_e_cg*1000),  "lp");  
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaPhi_Endcap->SaveAs("plots/cDeltaPhi_Endcap.png");
  
/*
  TCanvas * cDeltaEta = new TCanvas ("cDeltaEta", "cDeltaEta", 600, 500);
  cDeltaEta->cd();
  hDeltaEta_Barrel[NFILES-1]->Draw();
  hDeltaEta_Barrel[NFILES-1]->SetStats(0);
  hDeltaEta_Barrel[NFILES-1]->SetTitle(Form("Barrel: %d GeV %s",energies[0], particle_name.c_str()));
  hDeltaEta_Barrel[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range, fit_range);
  hDeltaEta_Barrel[NFILES-1]->GetYaxis()->SetRangeUser(0,hDeltaEta_Barrel_wGravity[NFILES-1]->GetMaximum()*1.5);
  hDeltaEta_Barrel[NFILES-1]->GetXaxis()->SetTitle("#eta_{reco} - #eta_{truth@VTX} [rad]");
  hDeltaEta_Barrel[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaEta_Barrel[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaEta_Barrel_wGravity[NFILES-1]->SetLineColor(kGreen+1);
  hDeltaEta_Barrel_wGravity[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.58,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      
      std::cout << "---- BARREL ----" << std::endl;      
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range, fit_range);
      fitGaus->SetLineColor(kBlack);
      hDeltaEta_Barrel[iFile]->Fit(fitGaus, "QR");
      std::cout << "eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_b = fitGaus->GetParameter(2);
//       eta_res_b = hDeltaEta_Barrel[iFile]->GetRMS();      
      
      fitGaus->SetLineColor(kGreen);
      hDeltaEta_Barrel_wGravity[iFile]->Fit(fitGaus, "QR");
      std::cout << "WEIGHTED: eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_b_cg = fitGaus->GetParameter(2);
//       eta_res_b_cg = hDeltaEta_Barrel_wGravity[iFile]->GetRMS();      

//       leg->AddEntry(hDeltaEta_Barrel[iFile],           Form("Seed:~RMS_{#eta}=%.1fx10^{-3}",eta_res_b*1000), "lp");          
//       leg->AddEntry(hDeltaEta_Barrel_wGravity[iFile],  Form("CoG:~~RMS_{#eta}=%.1fx10^{-3}",eta_res_b_cg*1000),  "lp");       
      
      leg->AddEntry(hDeltaEta_Barrel[iFile],           Form("Seed:~#sigma_{#eta}=%.2fx10^{-3}",eta_res_b*1000), "lp");          
      leg->AddEntry(hDeltaEta_Barrel_wGravity[iFile],  Form("CoG:~~#sigma_{#eta}=%.2fx10^{-3}",eta_res_b_cg*1000),  "lp");      
  }
  leg->Draw();
  if (SAVEPLOTS) cDeltaEta->SaveAs("plots/cDeltaEta_Barrel.png");
  */
/*
  TCanvas * cDeltaEta_Endcap = new TCanvas ("cDeltaEta_Endcap", "cDeltaEta_Endcap", 600, 500);
  cDeltaEta_Endcap->cd();
  hDeltaEta_Endcap[NFILES-1]->Draw();
  hDeltaEta_Endcap[NFILES-1]->SetStats(0);
  hDeltaEta_Endcap[NFILES-1]->SetTitle(Form("Endcaps: %d GeV %s",energies[0], particle_name.c_str()));
  hDeltaEta_Endcap[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range, fit_range);
  hDeltaEta_Endcap[NFILES-1]->GetYaxis()->SetRangeUser(0,hDeltaEta_Endcap_wGravity[NFILES-1]->GetMaximum()*1.5);
  hDeltaEta_Endcap[NFILES-1]->GetXaxis()->SetTitle("#eta_{reco} - #eta_{truth@VTX} [rad]");
  hDeltaEta_Endcap[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaEta_Endcap[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaEta_Endcap_wGravity[NFILES-1]->SetLineColor(kGreen+1);
  hDeltaEta_Endcap_wGravity[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.58,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      std::cout << "---- ENDCAP ----" << std::endl;
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range, fit_range);
      fitGaus->SetLineColor(kBlack);
      hDeltaEta_Endcap[iFile]->Fit(fitGaus, "QR");
      std::cout << "eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_e = fitGaus->GetParameter(2);
//       eta_res_e = hDeltaEta_Endcap[iFile]->GetRMS();      
      
      fitGaus->SetLineColor(kGreen);
      hDeltaEta_Endcap_wGravity[iFile]->Fit(fitGaus, "QR");
      std::cout << "WEIGHTED: eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_e_cg = fitGaus->GetParameter(2);
//       eta_res_e_cg = hDeltaEta_Endcap_wGravity[iFile]->GetRMS();      
      
//       leg->AddEntry(hDeltaEta_Endcap[iFile],           Form("Seed:~RMS_{#eta}=%.1fx10^{-3}",eta_res_e*1000), "lp");  
//       leg->AddEntry(hDeltaEta_Endcap_wGravity[iFile],  Form("CoG:~~RMS_{#eta}=%.1fx10^{-3}",eta_res_e_cg*1000),  "lp");     
      
      leg->AddEntry(hDeltaEta_Endcap[iFile],           Form("Seed:~#sigma_{#eta}=%.2fx10^{-3}",eta_res_e*1000), "lp");  
      leg->AddEntry(hDeltaEta_Endcap_wGravity[iFile],  Form("CoG:~~#sigma_{#eta}=%.2fx10^{-3}",eta_res_e_cg*1000),  "lp");     
  }
  leg->Draw();
  if (SAVEPLOTS) cDeltaEta_Endcap->SaveAs("plots/cDeltaEta_Endcap.png");
  
  */
  
  TCanvas * cDeltaTheta = new TCanvas ("cDeltaTheta", "cDeltaTheta", 600, 500);
  cDeltaTheta->cd();
  hDeltaTheta_Barrel[NFILES-1]->Draw();
  hDeltaTheta_Barrel[NFILES-1]->SetStats(0);
  hDeltaTheta_Barrel[NFILES-1]->SetTitle(Form("Barrel: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaTheta_Barrel[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range, fit_range);
  hDeltaTheta_Barrel[NFILES-1]->GetYaxis()->SetRangeUser(1,hDeltaTheta_Barrel_wGravity[NFILES-1]->GetMaximum()*5);
  hDeltaTheta_Barrel[NFILES-1]->GetXaxis()->SetTitle("#theta_{reco} - #theta_{truth@VTX} [rad]");
  hDeltaTheta_Barrel[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaTheta_Barrel[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaTheta_Barrel_wGravity[NFILES-1]->SetLineColor(kGreen+1);
  hDeltaTheta_Barrel_wGravity[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.58,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      
      std::cout << "---- BARREL ----" << std::endl;      
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range, fit_range);
      fitGaus->SetLineColor(kBlack);
      hDeltaTheta_Barrel[iFile]->Fit(fitGaus, "QR");
      std::cout << "eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_b = fitGaus->GetParameter(2);
//       eta_res_b = hDeltaTheta_Barrel[iFile]->GetRMS();      
      
      fitGaus->SetLineColor(kGreen);
      hDeltaTheta_Barrel_wGravity[iFile]->Fit(fitGaus, "QR");
      std::cout << "WEIGHTED: eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_b_cg = fitGaus->GetParameter(2);
//       eta_res_b_cg = hDeltaTheta_Barrel_wGravity[iFile]->GetRMS();      

//       leg->AddEntry(hDeltaTheta_Barrel[iFile],           Form("Seed:~RMS_{#eta}=%.1fx10^{-3}",eta_res_b*1000), "lp");          
//       leg->AddEntry(hDeltaTheta_Barrel_wGravity[iFile],  Form("CoG:~~RMS_{#eta}=%.1fx10^{-3}",eta_res_b_cg*1000),  "lp");       
      
      leg->AddEntry(hDeltaTheta_Barrel[iFile],           Form("Seed:~#sigma_{#theta}=%.2f mrad",eta_res_b*1000), "lp");          
      leg->AddEntry(hDeltaTheta_Barrel_wGravity[iFile],  Form("CoG:~~#sigma_{#theta}=%.2f mrad",eta_res_b_cg*1000),  "lp");      
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaTheta->SaveAs("plots/cDeltaTheta_Barrel.png");
  
    TCanvas * cDeltaTheta_Endcap = new TCanvas ("cDeltaTheta_Endcap", "cDeltaTheta_Endcap", 600, 500);
  cDeltaTheta_Endcap->cd();
  hDeltaTheta_Endcap[NFILES-1]->Draw();
  hDeltaTheta_Endcap[NFILES-1]->SetStats(0);
  hDeltaTheta_Endcap[NFILES-1]->SetTitle(Form("Endcaps: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaTheta_Endcap[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range, fit_range);
  hDeltaTheta_Endcap[NFILES-1]->GetYaxis()->SetRangeUser(1,hDeltaTheta_Endcap_wGravity[NFILES-1]->GetMaximum()*5);
  hDeltaTheta_Endcap[NFILES-1]->GetXaxis()->SetTitle("#theta_{reco} - #theta_{truth@VTX} [rad]");
  hDeltaTheta_Endcap[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaTheta_Endcap[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaTheta_Endcap_wGravity[NFILES-1]->SetLineColor(kGreen+1);
  hDeltaTheta_Endcap_wGravity[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.58,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      std::cout << "---- ENDCAP ----" << std::endl;
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range, fit_range);
      fitGaus->SetLineColor(kBlack);
      hDeltaTheta_Endcap[iFile]->Fit(fitGaus, "QR");
      std::cout << "eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_e = fitGaus->GetParameter(2);
//       eta_res_e = hDeltaTheta_Endcap[iFile]->GetRMS();      
      
      fitGaus->SetLineColor(kGreen);
      hDeltaTheta_Endcap_wGravity[iFile]->Fit(fitGaus, "QR");
      std::cout << "WEIGHTED: eta mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      eta_res_e_cg = fitGaus->GetParameter(2);
//       eta_res_e_cg = hDeltaTheta_Endcap_wGravity[iFile]->GetRMS();      
      
//       leg->AddEntry(hDeltaTheta_Endcap[iFile],           Form("Seed:~RMS_{#eta}=%.1fx10^{-3}",eta_res_e*1000), "lp");  
//       leg->AddEntry(hDeltaTheta_Endcap_wGravity[iFile],  Form("CoG:~~RMS_{#eta}=%.1fx10^{-3}",eta_res_e_cg*1000),  "lp");     
      
      leg->AddEntry(hDeltaTheta_Endcap[iFile],           Form("Seed:~#sigma_{#theta}=%.2f mrad",eta_res_e*1000), "lp");  
      leg->AddEntry(hDeltaTheta_Endcap_wGravity[iFile],  Form("CoG:~~#sigma_{#theta}=%.2f marad",eta_res_e_cg*1000),  "lp");     
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaTheta_Endcap->SaveAs("plots/cDeltaTheta_Endcap.png");
  
  
  //*********************************************************************************//
  //                                Timing layers
  //*********************************************************************************//
  
  fit_range = 0.03;
  
  TCanvas * cDeltaPhi_Timing = new TCanvas ("cDeltaPhi_Timing", "cDeltaPhi_Timing", 600, 500);
  cDeltaPhi_Timing->cd();
  hDeltaPhi_TimingBarrel_F[NFILES-1]->Draw();
  hDeltaPhi_TimingBarrel_F[NFILES-1]->SetStats(0);
  hDeltaPhi_TimingBarrel_F[NFILES-1]->SetTitle(Form("Barrel: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaPhi_TimingBarrel_F[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range*scale, fit_range*scale);
  hDeltaPhi_TimingBarrel_F[NFILES-1]->GetYaxis()->SetRangeUser(1,hDeltaPhi_TimingBarrel_Comb[NFILES-1]->GetMaximum()*5);
  hDeltaPhi_TimingBarrel_F[NFILES-1]->GetXaxis()->SetTitle("#phi_{reco} - #phi_{truth@VTX} [rad]");
  hDeltaPhi_TimingBarrel_F[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaPhi_TimingBarrel_F[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaPhi_TimingBarrel_R[NFILES-1]->SetLineColor(kRed+1);
  hDeltaPhi_TimingBarrel_R[NFILES-1]->Draw("same");
  
//   hDeltaPhi_TimingBarrel_W_R[NFILES-1]->SetLineColor(kBlue+1);
  hDeltaPhi_TimingBarrel_W_R[NFILES-1]->Draw("same");
  
  hDeltaPhi_TimingBarrel_Comb[NFILES-1]->SetLineColor(kGreen+2);
  hDeltaPhi_TimingBarrel_Comb[NFILES-1]->SetFillColor(kGreen+1);
  hDeltaPhi_TimingBarrel_Comb[NFILES-1]->SetFillStyle(3004);
//   hDeltaPhi_TimingBarrel_Comb[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.73,0.5,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
            
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range/10, fit_range/10);
      fitGaus->SetNpx(100);
      fitGaus->SetLineColor(kGreen);
      hDeltaPhi_TimingBarrel_Comb[iFile]->Fit(fitGaus, "0SQR");
      std::cout << "---- BARREL ----" << std::endl;      
      std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_b = fitGaus->GetParameter(2);
      
      leg->AddEntry(hDeltaPhi_TimingBarrel_F[iFile],          "T1 Seed (front)", "lp");          
      leg->AddEntry(hDeltaPhi_TimingBarrel_R[iFile],          "T2 Seed (rear)", "lp");          
      leg->AddEntry(hDeltaTheta_TimingBarrel_W_R[iFile],  Form("CoG: #sigma_{#theta}=%.2f mrad", phi_res_b*1000), "lp");          
//       leg->AddEntry(hDeltaPhi_TimingBarrel_Comb[iFile],  Form("Combined: #sigma_{#phi}=%.2f mrad",phi_res_b*1000), "lp");          
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaPhi_Timing->SaveAs("plots/cDeltaPhi_TimingBarrel.png");
  
  
//   TCanvas * cDeltaEta_Timing = new TCanvas ("cDeltaEta_Timing", "cDeltaEta_Timing", 600, 500);
//   cDeltaEta_Timing->cd();
//   hDeltaEta_TimingBarrel_F[NFILES-1]->Draw();
//   hDeltaEta_TimingBarrel_F[NFILES-1]->SetStats(0);
//   hDeltaEta_TimingBarrel_F[NFILES-1]->SetTitle(Form("Barrel: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaEta_TimingBarrel_F[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range*scale, fit_range*scale);
//   hDeltaEta_TimingBarrel_F[NFILES-1]->GetYaxis()->SetRangeUser(0,hDeltaEta_TimingBarrel_Comb[NFILES-1]->GetMaximum()*1.5);
//   hDeltaEta_TimingBarrel_F[NFILES-1]->GetXaxis()->SetTitle("#eta_{reco} - #eta_{truth@VTX} [rad]");
//   hDeltaEta_TimingBarrel_F[NFILES-1]->GetYaxis()->SetTitle("Counts");
//   hDeltaEta_TimingBarrel_F[NFILES-1]->SetLineColor(kBlack);
//   
//   hDeltaEta_TimingBarrel_R[NFILES-1]->SetLineColor(kRed+1);
//   hDeltaEta_TimingBarrel_R[NFILES-1]->Draw("same");
//   
//   hDeltaEta_TimingBarrel_Comb[NFILES-1]->SetLineColor(kGreen+2);
//   hDeltaEta_TimingBarrel_Comb[NFILES-1]->SetFillColor(kGreen+1);
//   hDeltaEta_TimingBarrel_Comb[NFILES-1]->SetFillStyle(3004);
//   hDeltaEta_TimingBarrel_Comb[NFILES-1]->Draw("same");
// 
//   leg = new TLegend(0.15,0.73,0.5,0.88,NULL,"brNDC");
// 
//   for (int iFile = NFILES-1; iFile>=0; iFile--)
//   {
//             
//       TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range/10, fit_range/10);
//       fitGaus->SetNpx(100);
//       fitGaus->SetLineColor(kGreen);
//       hDeltaEta_TimingBarrel_Comb[iFile]->Fit(fitGaus, "0SQR");
//       std::cout << "---- BARREL ----" << std::endl;      
//       std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
//       phi_res_b = fitGaus->GetParameter(2);
//       
//       leg->AddEntry(hDeltaEta_TimingBarrel_F[iFile],          "T1 Seed (front)", "lp");          
//       leg->AddEntry(hDeltaEta_TimingBarrel_R[iFile],          "T2 Seed (rear)", "lp");          
//       leg->AddEntry(hDeltaEta_TimingBarrel_Comb[iFile],  Form("Combined: #sigma_{#eta}=%.2fx10^{-3}",phi_res_b*1000), "lp");          
//   }
//   leg->Draw();
//   if (SAVEPLOTS) cDeltaEta_Timing->SaveAs("plots/cDeltaEta_TimingBarrel.png");
  
  TCanvas * cDeltaTheta_Timing = new TCanvas ("cDeltaTheta_Timing", "cDeltaTheta_Timing", 600, 600);
  cDeltaTheta_Timing->cd();
  hDeltaTheta_TimingBarrel_F[NFILES-1]->Draw();
  hDeltaTheta_TimingBarrel_F[NFILES-1]->SetStats(0);
  hDeltaTheta_TimingBarrel_F[NFILES-1]->SetTitle(Form("Barrel: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaTheta_TimingBarrel_F[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range*scale, fit_range*scale);
//   hDeltaTheta_TimingBarrel_F[NFILES-1]->GetXaxis()->SetRangeUser(-0.003, 0.003);
  hDeltaTheta_TimingBarrel_F[NFILES-1]->GetYaxis()->SetRangeUser(1 , hDeltaTheta_TimingBarrel_W_F[NFILES-1]->GetMaximum()*5);
  hDeltaTheta_TimingBarrel_F[NFILES-1]->GetXaxis()->SetTitle("#theta_{reco} - #theta_{truth@VTX} [rad]");
  hDeltaTheta_TimingBarrel_F[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaTheta_TimingBarrel_F[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaTheta_TimingBarrel_R[NFILES-1]->SetLineColor(kRed+1);
  hDeltaTheta_TimingBarrel_R[NFILES-1]->Draw("same");
  
  hDeltaTheta_TimingBarrel_W_F[NFILES-1]->SetLineColor(kGreen+1);
  hDeltaTheta_TimingBarrel_W_F[NFILES-1]->Draw("same");
  
  hDeltaTheta_TimingBarrel_Comb[NFILES-1]->SetLineColor(kGreen+2);
  hDeltaTheta_TimingBarrel_Comb[NFILES-1]->SetFillColor(kGreen+1);
  hDeltaTheta_TimingBarrel_Comb[NFILES-1]->SetFillStyle(3004);
//   hDeltaTheta_TimingBarrel_Comb[NFILES-1]->Draw("same");

  leg = new TLegend(0.15,0.73,0.5,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
            
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range/60, fit_range/60);
      fitGaus->SetNpx(100);
      fitGaus->SetLineColor(kGreen+1);
//       hDeltaTheta_TimingBarrel_Comb[iFile]->Fit(fitGaus, "0SQR");
      hDeltaTheta_TimingBarrel_W_F[iFile]->Fit(fitGaus, "QR");
      
      std::cout << "---- BARREL ----" << std::endl;      
      std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_b = fitGaus->GetParameter(2);
      
      leg->AddEntry(hDeltaTheta_TimingBarrel_F[iFile],          "T1 Seed (front)", "lp");          
      leg->AddEntry(hDeltaTheta_TimingBarrel_W_F[iFile],  Form("CoG: #sigma_{#theta}=%.2f mrad", phi_res_b*1000), "lp");          
      leg->AddEntry(hDeltaTheta_TimingBarrel_R[iFile],          "T2 Seed (rear)", "lp");          
//       leg->AddEntry(hDeltaTheta_TimingBarrel_Comb[iFile],  Form("Combined: #sigma_{#theta}=%.2f mrad",phi_res_b*1000), "lp");          
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaTheta_Timing->SaveAs("plots/cDeltaTheta_TimingBarrel.png");
  
  
  
  TCanvas * cDeltaPhi_TimingEndcap = new TCanvas ("cDeltaPhi_TimingEndcap", "cDeltaPhi_TimingEndcap", 600, 500);
  cDeltaPhi_TimingEndcap->cd();
  hDeltaPhi_TimingEndcap_F[NFILES-1]->Draw();
  hDeltaPhi_TimingEndcap_F[NFILES-1]->SetStats(0);
  hDeltaPhi_TimingEndcap_F[NFILES-1]->SetTitle(Form("Endcap: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaPhi_TimingEndcap_F[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range*scale, fit_range*scale);
  hDeltaPhi_TimingEndcap_F[NFILES-1]->GetYaxis()->SetRangeUser(1,hDeltaPhi_TimingEndcap_Comb[NFILES-1]->GetMaximum()*5);
  hDeltaPhi_TimingEndcap_F[NFILES-1]->GetXaxis()->SetTitle("#phi_{reco} - #phi_{truth@VTX} [rad]");
  hDeltaPhi_TimingEndcap_F[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaPhi_TimingEndcap_F[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaPhi_TimingEndcap_R[NFILES-1]->SetLineColor(kRed+1);
  hDeltaPhi_TimingEndcap_R[NFILES-1]->Draw("same");
  
  hDeltaPhi_TimingEndcap_Comb[NFILES-1]->SetLineColor(kGreen+2);
  hDeltaPhi_TimingEndcap_Comb[NFILES-1]->SetFillColor(kGreen+1);
  hDeltaPhi_TimingEndcap_Comb[NFILES-1]->SetFillStyle(3004);
  hDeltaPhi_TimingEndcap_Comb[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.73,0.5,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
            
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range/10, fit_range/10);
      fitGaus->SetNpx(100);
      fitGaus->SetLineColor(kGreen);
      hDeltaPhi_TimingEndcap_Comb[iFile]->Fit(fitGaus, "0SQR");
      std::cout << "---- ENDCAP ----" << std::endl;      
      std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_b = fitGaus->GetParameter(2);
      
      leg->AddEntry(hDeltaPhi_TimingEndcap_F[iFile],          "T1 Seed (front)", "lp");          
      leg->AddEntry(hDeltaPhi_TimingEndcap_R[iFile],          "T2 Seed (rear)", "lp");          
      leg->AddEntry(hDeltaPhi_TimingEndcap_Comb[iFile],  Form("Combined: #sigma_{#phi}=%.2f mrad",phi_res_b*1000), "lp");  
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaPhi_TimingEndcap->SaveAs("plots/cDeltaPhi_TimingEndcap.png");
  
  
//   TCanvas * cDeltaEta_TimingEndcap = new TCanvas ("cDeltaEta_TimingEndcap", "cDeltaEta_TimingEndcap", 600, 500);
//   cDeltaEta_TimingEndcap->cd();
//   hDeltaEta_TimingEndcap_F[NFILES-1]->Draw();
//   hDeltaEta_TimingEndcap_F[NFILES-1]->SetStats(0);
//   hDeltaEta_TimingEndcap_F[NFILES-1]->SetTitle(Form("Endcap: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaEta_TimingEndcap_F[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range*scale, fit_range*scale);
//   hDeltaEta_TimingEndcap_F[NFILES-1]->GetYaxis()->SetRangeUser(0,hDeltaEta_TimingEndcap_Comb[NFILES-1]->GetMaximum()*1.5);
//   hDeltaEta_TimingEndcap_F[NFILES-1]->GetXaxis()->SetTitle("#eta_{reco} - #eta_{truth@VTX} [rad]");
//   hDeltaEta_TimingEndcap_F[NFILES-1]->GetYaxis()->SetTitle("Counts");
//   hDeltaEta_TimingEndcap_F[NFILES-1]->SetLineColor(kBlack);
//   
//   hDeltaEta_TimingEndcap_R[NFILES-1]->SetLineColor(kRed+1);
//   hDeltaEta_TimingEndcap_R[NFILES-1]->Draw("same");
//   
//   hDeltaEta_TimingEndcap_Comb[NFILES-1]->SetLineColor(kGreen+2);
//   hDeltaEta_TimingEndcap_Comb[NFILES-1]->SetFillColor(kGreen+1);
//   hDeltaEta_TimingEndcap_Comb[NFILES-1]->SetFillStyle(3004);
//   hDeltaEta_TimingEndcap_Comb[NFILES-1]->Draw("same");
//   
//   leg = new TLegend(0.15,0.73,0.5,0.88,NULL,"brNDC");
// 
//   for (int iFile = NFILES-1; iFile>=0; iFile--)
//   {
//             
//       TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range/10, fit_range/10);
//       fitGaus->SetNpx(100);
//       fitGaus->SetLineColor(kGreen);
//       hDeltaEta_TimingEndcap_Comb[iFile]->Fit(fitGaus, "0SQR");
//       std::cout << "---- ENDCAP ----" << std::endl;      
//       std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
//       phi_res_b = fitGaus->GetParameter(2);
//       
//       leg->AddEntry(hDeltaEta_TimingEndcap_F[iFile],          "T1 Seed (front)", "lp");          
//       leg->AddEntry(hDeltaEta_TimingEndcap_R[iFile],          "T2 Seed (rear)", "lp");          
//       leg->AddEntry(hDeltaEta_TimingEndcap_Comb[iFile],  Form("Combined: #sigma_{#eta}=%.2fx10^{-3}",phi_res_b*1000), "lp");        
//       
//   }
//   leg->Draw();
//   if (SAVEPLOTS) cDeltaEta_TimingEndcap->SaveAs("plots/cDeltaEta_TimingEndcap.png");
//   
  
  TCanvas * cDeltaTheta_TimingEndcap = new TCanvas ("cDeltaTheta_TimingEndcap", "cDeltaTheta_TimingEndcap", 600, 500);
  cDeltaTheta_TimingEndcap->cd();
  hDeltaTheta_TimingEndcap_F[NFILES-1]->Draw();
  hDeltaTheta_TimingEndcap_F[NFILES-1]->SetStats(0);
  hDeltaTheta_TimingEndcap_F[NFILES-1]->SetTitle(Form("Endcap: %d GeV %s",energies[0], particle_name.c_str()));
//   hDeltaTheta_TimingEndcap_F[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range*scale, fit_range*scale);
  hDeltaTheta_TimingEndcap_F[NFILES-1]->GetYaxis()->SetRangeUser(1,hDeltaTheta_TimingEndcap_Comb[NFILES-1]->GetMaximum()*5);
  hDeltaTheta_TimingEndcap_F[NFILES-1]->GetXaxis()->SetTitle("#theta_{reco} - #theta_{truth@VTX} [rad]");
  hDeltaTheta_TimingEndcap_F[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hDeltaTheta_TimingEndcap_F[NFILES-1]->SetLineColor(kBlack);
  
  hDeltaTheta_TimingEndcap_R[NFILES-1]->SetLineColor(kRed+1);
  hDeltaTheta_TimingEndcap_R[NFILES-1]->Draw("same");
  
  hDeltaTheta_TimingEndcap_Comb[NFILES-1]->SetLineColor(kGreen+2);
  hDeltaTheta_TimingEndcap_Comb[NFILES-1]->SetFillColor(kGreen+1);
  hDeltaTheta_TimingEndcap_Comb[NFILES-1]->SetFillStyle(3004);
  hDeltaTheta_TimingEndcap_Comb[NFILES-1]->Draw("same");
  
  leg = new TLegend(0.15,0.73,0.5,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
            
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -fit_range/10, fit_range/10);
      fitGaus->SetNpx(100);
      fitGaus->SetLineColor(kGreen);
      hDeltaTheta_TimingEndcap_Comb[iFile]->Fit(fitGaus, "0SQR");
      std::cout << "---- ENDCAP ----" << std::endl;      
      std::cout << "phi mean = " << fitGaus->GetParameter(1) << " :: resolution = " << fitGaus->GetParameter(2) << " rad  (" << fitGaus->GetParameter(2)/TMath::DegToRad() << " deg)" << std::endl;
      phi_res_b = fitGaus->GetParameter(2);
      
      leg->AddEntry(hDeltaTheta_TimingEndcap_F[iFile],          "T1 Seed (front)", "lp");          
      leg->AddEntry(hDeltaTheta_TimingEndcap_R[iFile],          "T2 Seed (rear)", "lp");          
      leg->AddEntry(hDeltaTheta_TimingEndcap_Comb[iFile],  Form("Combined: #sigma_{#theta}=%.2f mrad",phi_res_b*1000), "lp");        
      
  }
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cDeltaTheta_TimingEndcap->SaveAs("plots/cDeltaTheta_TimingEndcap.png");
  
  
  
  
  /// other plots
  
  TCanvas * cNHitsTiming = new TCanvas ("cNHitsTiming", "cNHitsTiming", 600, 500);
  cNHitsTiming->cd();
  hNHits_T1Barrel[NFILES-1]->Draw();
  hNHits_T1Barrel[NFILES-1]->SetStats(0);
  hNHits_T1Barrel[NFILES-1]->SetTitle(Form("Barrel: %d GeV %s",energies[0], particle_name.c_str()));
  hNHits_T1Barrel[NFILES-1]->GetXaxis()->SetRangeUser(0, 34);
  hNHits_T1Barrel[NFILES-1]->GetYaxis()->SetRangeUser(0,hNHits_T1Barrel[NFILES-1]->GetMaximum()*1.5);
  hNHits_T1Barrel[NFILES-1]->GetXaxis()->SetTitle("Number of hits / event");
  hNHits_T1Barrel[NFILES-1]->GetYaxis()->SetTitle("Counts");
  hNHits_T1Barrel[NFILES-1]->SetLineColor(kBlack);
  
  hNHits_T2Barrel[NFILES-1]->SetLineColor(kRed);
  hNHits_T2Barrel[NFILES-1]->Draw("same");
  
  TCanvas * cScatterTOF_barrel = new TCanvas("cScatterTOF_barrel", "cScatterTOF_barrel", 600, 500);
  cScatterTOF_barrel->cd();
  hScatterTOF_TimingBarrel[0]->Draw("COLZ");
  hScatterTOF_TimingBarrel[0]->SetStats(0);
  hScatterTOF_TimingBarrel[0]->SetTitle("Barrel TOF");
  hScatterTOF_TimingBarrel[0]->GetXaxis()->SetTitle("Distance from IP [mm]");
  hScatterTOF_TimingBarrel[0]->GetYaxis()->SetTitle("Time of hit at timing layer [ns]");
  hScatterTOF_TimingBarrel[0]->GetXaxis()->SetRangeUser(1760,2440);
  hScatterTOF_TimingBarrel[0]->GetYaxis()->SetRangeUser(5.5,8.5);
  if (SAVEPLOTS)   cScatterTOF_barrel->SaveAs("plots/cScatterTOF_barrel.png");
  
  
    TCanvas * cScatterTOF_endcap = new TCanvas("cScatterTOF_endcap", "cScatterTOF_endcap", 600, 500);
  cScatterTOF_endcap->cd();
  hScatterTOF_TimingEndcap[0]->Draw("COLZ");
  hScatterTOF_TimingEndcap[0]->SetStats(0);
  hScatterTOF_TimingEndcap[0]->SetTitle("Endcap TOF");
  hScatterTOF_TimingEndcap[0]->GetXaxis()->SetTitle("Distance from IP [mm]");
  hScatterTOF_TimingEndcap[0]->GetYaxis()->SetTitle("Time of hit at timing layer [ns]");
  hScatterTOF_TimingEndcap[0]->GetXaxis()->SetRangeUser(1760,2440);
  hScatterTOF_TimingEndcap[0]->GetYaxis()->SetRangeUser(5.5,8.5);
  if (SAVEPLOTS)   cScatterTOF_endcap->SaveAs("plots/cScatterTOF_endcap.png");
  
  
    TCanvas * cTimeRes_TOF_barrel = new TCanvas("cTimeRes_TOF_barrel", "cTimeRes_TOF_barrel", 600, 500);
  cTimeRes_TOF_barrel->cd();
  hTOF_TimingBarrel[0]->Draw("hist");
//   hTOF_TimingBarrel[0]->SetStats(0);
  hTOF_TimingBarrel[0]->SetFillColor(kCyan+1);
  hTOF_TimingBarrel[0]->SetTitle("Barrel time resolution");
  hTOF_TimingBarrel[0]->GetXaxis()->SetTitle("distance/c -TOF [ps]");
  hTOF_TimingBarrel[0]->GetYaxis()->SetTitle("Counts");
  hTOF_TimingBarrel[0]->GetXaxis()->SetRangeUser( -40, 40  );
//   hTOF_TimingBarrel[0]->Fit("gaus");
  if (SAVEPLOTS) cTimeRes_TOF_barrel->SaveAs("plots/cTimeRes_TOF_barrel.png");
    
    
  TCanvas * cTimeRes_TOF_endcap = new TCanvas("cTimeRes_TOF_endcap", "cTimeRes_TOF_endcap", 600, 500);
  cTimeRes_TOF_endcap->cd();
  hTOF_TimingEndcap[0]->Draw("hist");
//   hTOF_TimingEndcap[0]->SetStats(0);
  hTOF_TimingEndcap[0]->SetFillColor(kCyan+1);
  hTOF_TimingEndcap[0]->SetTitle("Endcap time resolution");
  hTOF_TimingEndcap[0]->GetXaxis()->SetTitle("distance/c -TOF [ps]");
  hTOF_TimingEndcap[0]->GetYaxis()->SetTitle("Counts");
  hTOF_TimingEndcap[0]->GetXaxis()->SetRangeUser( -40, 40  );
//   hTOF_TimingEndcap[0]->Fit("gaus");
  if (SAVEPLOTS) cTimeRes_TOF_endcap->SaveAs("plots/cTimeRes_TOF_endcap.png");
  
  
  
  
  TGraphErrors * gEcalThetaRes_vs_ene       = new TGraphErrors();    
  for (int iEne = 0; iEne<ENE_BINS; iEne++)
  {
      float ave_ene = (vMinEneBin[iEne]+vMaxEneBin[iEne])/2.;
      double sigma_eff = FindSmallestInterval(vDeltaThetaECAL[iEne], 0.68, false) /2. ;
      
      float mean = 0;      
      for (long unsigned int it = 0; it < vDeltaThetaECAL[iEne]->size(); it++)
      {
          mean += vDeltaThetaECAL[iEne]->at(it);
      }
      mean/=vDeltaThetaECAL[iEne]->size();
      gEcalThetaRes_vs_ene->SetPoint(iEne, ave_ene, sigma_eff*1000);
      
      
      
  }
  
  TCanvas * cEcalThetaRes_vs_Ene = new TCanvas ("cEcalThetaRes_vs_Ene", "cEcalThetaRes_vs_Ene", 600, 600);
  cEcalThetaRes_vs_Ene->cd();
  gEcalThetaRes_vs_ene->Draw("APE");
  gEcalThetaRes_vs_ene->SetTitle("; Energy [GeV];#sigma_{#theta} [mrad]");
  gEcalThetaRes_vs_ene->SetMarkerStyle(20);
  gEcalThetaRes_vs_ene->SetMaximum(0.55);
  gEcalThetaRes_vs_ene->SetMinimum(0.35);
  gEcalThetaRes_vs_ene->GetXaxis()->SetLimits(0, 100);
  
  TF1* fitResolution = new TF1 ("fitResolution", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 4, 100);    
  fitResolution->SetLineWidth(2);    
//   fitResolution->SetLineStyle(7);    
  fitResolution->SetParameters(1.4, 0.3);
  fitResolution->SetParLimits(1,0.1, 0.5);
  
  
  leg = new TLegend(0.35,0.73,0.88,0.88,NULL,"brNDC");
  fitResolution->SetLineColor(kRed);    
  for (int i = 0; i<30; i++) gEcalThetaRes_vs_ene->Fit(fitResolution, "QRS");
  float stoch_term = fitResolution->GetParameter(0);
  float const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gEcalThetaRes_vs_ene, Form("#sigma_{theta} (mrad) = %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
  fitResolution->Draw("same");
  
  leg->Draw();
  
  
  
  
  
/*  
  TGraphErrors * gTimingThetaRes_vs_theta       = new TGraphErrors();    
  for (int iTheta = 0; iTheta<THETA_BINS; iTheta++)
  {
      
      float ave_theta = (vMinThetaBin[iTheta]+vMaxThetaBin[iTheta])/2.;
      std::cout << "ave_theta = " << ave_theta << std::endl;
      double sigma_eff = 0;
      if (vDeltaThetaTiming[iTheta]->size()>0) sigma_eff = FindSmallestInterval(vDeltaThetaTiming[iTheta], 0.68, false) /2. ;
      std::cout << "sigma_eff = " << sigma_eff << std::endl;
      
      float mean = 0;      
      for (long unsigned int it = 0; it < vDeltaThetaTiming[iTheta]->size(); it++)
      {
          mean += vDeltaThetaTiming[iTheta]->at(it);
      }
      mean/=vDeltaThetaTiming[iTheta]->size();
//       gTimingThetaRes_vs_theta->SetPoint(iTheta, ave_theta/M_PI*180, sigma_eff*1000);
      gTimingThetaRes_vs_theta->SetPoint(iTheta, ave_theta/M_PI*180, sin(sigma_eff)*1785);
      
      
      
  }
  
  TCanvas * cEcalThetaRes_vs_Theta = new TCanvas ("cEcalThetaRes_vs_Theta", "cEcalThetaRes_vs_Theta", 600, 600);
  cEcalThetaRes_vs_Theta->cd();
  gTimingThetaRes_vs_theta->Draw("APE");
//   gTimingThetaRes_vs_theta->SetTitle(";#theta [rad];#sigma_{#theta} [mm]");
  gTimingThetaRes_vs_theta->SetTitle(";#theta [degrees]; #sigma_{z} [mm]");
  gTimingThetaRes_vs_theta->SetMarkerStyle(20);
  gTimingThetaRes_vs_theta->SetMaximum(1.2);
  gTimingThetaRes_vs_theta->SetMinimum(0);
  gTimingThetaRes_vs_theta->GetXaxis()->SetLimits(0, 45);
  
  fitResolution = new TF1 ("fitResolution", "pol2", 0, 45);    
  fitResolution->SetLineWidth(2);    
//   fitResolution->SetLineStyle(7);    
  fitResolution->SetParameters(1., -0.02, 0.0002);
  gTimingThetaRes_vs_theta->Fit(fitResolution, "R");*/
//   fitResolution->SetParLimits(1,0.1, 0.5);
//   
//   
//   leg = new TLegend(0.35,0.73,0.88,0.88,NULL,"brNDC");
//   fitResolution->SetLineColor(kRed);    
//   for (int i = 0; i<30; i++) gTimingThetaRes_vs_theta->Fit(fitResolution, "QRS");
//   stoch_term = fitResolution->GetParameter(0);
//   const_term = fitResolution->GetParameter(1);        
//   leg->AddEntry(gTimingThetaRes_vs_theta, Form("#sigma_{theta} (mrad) = %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
//   fitResolution->Draw("same");
  
//   leg->Draw();
  
   
   theApp->Run();
}


