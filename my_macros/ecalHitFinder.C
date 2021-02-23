// g++ -Wall -o ecalHitFinder ecalHitFinder.C  myG4Tree.cc myG4Tree.hh myTruthTree.cc myTruthTree.hh recoUtils.cc recoUtils.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs`


#include "SCEPCal_GeometryHelper.hh"
#include "myG4Tree.hh"
#include "myTruthTree.hh"
#include "recoUtils.hh"

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
#include "THStack.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"

#include "TSpline.h"
#include "TEllipse.h"
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
//   int energy = 100;
  
  std::string output_tag = "wwln";
  int selEv = 0;
  if (argc > 1) selEv = atoi(argv[1]);
    
  //define histos
  
  int NPHI_TL1   = 186-1;
  int NTHETA_TL1 = 29*16-1;
  int NPHI_TL2   = 186*16-1;
  int NTHETA_TL2 = 29-1;
  
//   int NPHI_EC    = 1130;
//   int NTHETA_EC  = 180+180-1;
  int NPHI_EC    = 1130/4;
  int NTHETA_EC  = (180+180-1)/4;
  int NPHI_DRT   = 36;
  int NTHETA_DRT = 40+40-1;
  
  double minPhi = -M_PI;
  double maxPhi = M_PI;  
  double minTheta = 0;
  double maxTheta = M_PI;
  
  double bin_width_theta_TL1 = (maxTheta-minTheta)/NTHETA_TL1;
  double bin_width_theta_TL2 = (maxTheta-minTheta)/NTHETA_TL2;
  double bin_width_theta_EC  = (maxTheta-minTheta)/NTHETA_EC;
  double bin_width_theta_DRT = (maxTheta-minTheta)/NTHETA_DRT;
  
  double bin_width_phi_TL1 = (maxPhi-minPhi)/NPHI_TL1;
  double bin_width_phi_TL2 = (maxPhi-minPhi)/NPHI_TL2;
  double bin_width_phi_EC  = (maxPhi-minPhi)/NPHI_EC;
  double bin_width_phi_DRT = (maxPhi-minPhi)/NPHI_DRT;
  
  
  TH2F * hGrid_T1 = new TH2F ("hGrid_T1", "hGrid_T1", NTHETA_TL1, minTheta-bin_width_theta_TL1/2, maxTheta-bin_width_theta_TL1/2, NPHI_TL1, minPhi-bin_width_phi_TL1/2, maxPhi-bin_width_phi_TL1/2);
  TH2F * hGrid_T2 = new TH2F ("hGrid_T2", "hGrid_T2", NTHETA_TL2, minTheta-bin_width_theta_TL2/2, maxTheta-bin_width_theta_TL2/2, NPHI_TL2, minPhi-bin_width_phi_TL2/2, maxPhi-bin_width_phi_TL2/2);
  
  TH2F * hGrid_EC_F = new TH2F ("hGrid_EC_F", "hGrid_EC_F", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hGrid_EC_R = new TH2F ("hGrid_EC_R", "hGrid_EC_R", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hGrid_EC_T = new TH2F ("hGrid_EC_T", "hGrid_EC_T", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  
  TH2F * hGrid_DRT_S = new TH2F ("hGrid_DRT_S", "hGrid_DRT_S", NTHETA_DRT, minTheta-bin_width_theta_DRT/2, maxTheta-bin_width_theta_DRT/2, NPHI_DRT, minPhi-bin_width_phi_DRT/2, maxPhi-bin_width_phi_DRT/2);
  TH2F * hGrid_DRT_C = new TH2F ("hGrid_DRT_C", "hGrid_DRT_C", NTHETA_DRT, minTheta-bin_width_theta_DRT/2, maxTheta-bin_width_theta_DRT/2, NPHI_DRT, minPhi-bin_width_phi_DRT/2, maxPhi-bin_width_phi_DRT/2);
  
  TH2F * hTruthFloor = new TH2F ("hTruthFloor", "hTruthFloor", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hTruthChargedEM = new TH2F ("hTruthChargedEM", "hTruthChargedEM", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hTruthNeutralEM = new TH2F ("hTruthNeutralEM", "hTruthNeutralEM", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hTruthChargedHAD = new TH2F ("hTruthChargedHAD", "hTruthChargedHAD", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hTruthNeutralHAD = new TH2F ("hTruthNeutralHAD", "hTruthNeutralHAD", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  
  
  
  SCEPCal_GeometryHelper myGeometry;

  
  //run over energy scan
//   TFile * RunFile = new TFile("../root_files/iso_gun/iso_gun_mu_100GeV_T+E.root","READ"); 
//   TFile * RunFile = new TFile("../root_files/iso_gun/central_kaon0L_10GeV_ALL.root","READ"); 
  
//   TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_gamma_Iso+Uniform1-100_GeV.root","READ");       
//   TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_e-_Iso+Uniform1-100_GeV.root","READ");         
  
  
//   TFile * RecoFile = new TFile("../root_files/hep_outputs/output_hep_test.root","READ");       
  TFile * RecoFile = new TFile("../root_files/hep_outputs/output_SCEPCal_wwlj100k_job_12.root","READ");       
  
  
  TTree* TreeRun = (TTree*) RecoFile->Get("B4");
  myG4TreeVars myTV;
  InitG4Tree (TreeRun, myTV);
      
  
  ///*******************************************///
  ///		 Run over events	    ///
  ///*******************************************///
    
  int NEVENTS = TreeRun->GetEntries();
  std::cout << "NEVENTS = " << NEVENTS << std::endl;
  
  
//       TFile * TruthFile = new TFile("../root_files/hep_outputs/hep_truth.root","READ");
  TFile * TruthFile = new TFile("../../HepMC_Files/wwlj100k_job_12_output_tuple.root","READ");
  TTree* TruthTree = (TTree*) TruthFile->Get("truth");
  myTruthTreeVars myTruthTV;
  InitTruthTree (TruthTree, myTruthTV);
  
  TruthTree->GetEntry(selEv);
  std::cout << "n particles found: " << myTruthTV.mcs_n << std::endl;
  std::cout << "\n*******************************************************\n" << std::endl;

  
  for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
  {
      std::cout << "particle: " << i 
      << ", \npdgId= " <<  myTruthTV.mcs_pdgId->at(i) 
      << ", \nmass= " <<  myTruthTV.mcs_m->at(i) 
      << ", \ncharge= " <<  myTruthTV.mcs_charge->at(i) 
      << ", \nenergy= " << myTruthTV.mcs_E->at(i) << " GeV" << std::endl;

      int    pdgId = myTruthTV.mcs_pdgId->at(i);
      double ene   = myTruthTV.mcs_E->at(i);
      double phi   = myTruthTV.mcs_phi->at(i);
      double eta   = myTruthTV.mcs_eta->at(i);
      int charge = myTruthTV.mcs_charge->at(i);
      double theta = 2*atan(exp(-eta));
      theta = M_PI- theta;

      
      if (charge!=0 && abs(pdgId)!= 11) 
      {
          hTruthChargedHAD ->Fill(theta, phi, ene);
          std::cout << "this is a non electron charged particle!" << std::endl;          
      }
      
      if (charge!=0 && abs(pdgId)== 11) 
      {
          hTruthChargedEM  ->Fill(theta, phi, ene);
          std::cout << "this is an electron!" << std::endl;
      }
      if (charge==0 && abs(pdgId)!= 22) 
      {
          hTruthNeutralHAD ->Fill(theta, phi, ene);
            std::cout << "this is a non photon neutral particle!" << std::endl;
      }
      if (charge==0 && pdgId== 22) 
      {
          std::cout << "this is a photon!" << std::endl;
          hTruthNeutralEM  ->Fill(theta, phi, ene);
      }
        
      std::cout << "******************************************************* " << std::endl;
      
  }
      
      
            

  

  float maxDeltaR = 0.1;
//   NEVENTS = 1000;
  
//   for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
//   {
      int iEvt = selEv;
                                        
     
      TreeRun->GetEntry(iEvt);
      std::cout << "processing event: " << iEvt << "\r" << std::flush;

          
      //**************************************************************//
      //                           DR HCAL
      //**************************************************************//

      float totDRHScint = 0;
      float totDRHCher  = 0;
      
      for (unsigned int i = 0; i<myTV.VectorL->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'l');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorL->at(i)/1000.;      
          double this_scint = myTV.VectorSignalsL->at(i);                
          double this_cher  = myTV.VectorSignalsCherL->at(i);
          hGrid_DRT_S ->Fill(this_theta, this_phi, this_scint);                              
          hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher);    
          totDRHScint+=this_scint;
          totDRHCher+=this_cher;
      }
      for (unsigned int i = 0; i<myTV.VectorR->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'r');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorR->at(i)/1000.;                    
          double this_scint = myTV.VectorSignalsR->at(i);     
          double this_cher  = myTV.VectorSignalsCherR->at(i);      
          hGrid_DRT_S ->Fill(this_theta, this_phi, this_scint);                              
          hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher);                              
          totDRHScint+=this_scint;
          totDRHCher+=this_cher;
      }
      
      std::cout << "Total scint in HCAL: " << totDRHScint << std::endl;
      std::cout << "Total cher in HCAL: " << totDRHCher << std::endl;
      
      //**************************************************************//
      //                             ECAL
      //**************************************************************//

      float ene_EC_th = 0.03;
      float totEcalEne = 0;
      for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
      {                            
              
          TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;                    
          
//           if (fabs(myTV.VecHit_CrystalID->at(i)) <1000000)
          if (this_ene>ene_EC_th)
          {
            hGrid_EC_F ->Fill(this_theta, this_phi, myTV.VecHit_ScepEneDepF->at(i)/1000);          
            hGrid_EC_R ->Fill(this_theta, this_phi, myTV.VecHit_ScepEneDepR->at(i)/1000);          
            hGrid_EC_T ->Fill(this_theta, this_phi, this_ene*100);
            totEcalEne+=this_ene;
          }
      }
      std::cout << "Total energy in ECAL: " << totEcalEne << std::endl;
      
      
      // find hit with energy above seed threshold
      float EC_seed_th = 0.4;
      std::vector<EcalSeed> myEcSeeds;
      for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
      {
              
          TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;
          
//           if (fabs(myTV.VecHit_CrystalID->at(i)) <1000000)
          if (this_ene>EC_seed_th)
          {
              EcalSeed new_seed;
              new_seed.Init(myTV.VecHit_CrystalID->at(i), this_theta, this_phi, this_ene);
              myEcSeeds.push_back(new_seed);
          }
      }
      
      
      std::cout << "Number of ECAL seeds found: " << myEcSeeds.size() << std::endl;
      
      std::cout << "Cleaning up seeds too close to each other" << std::endl;
      std::vector<EcalSeed>  myEcSeedsCleaned;
      for (long unsigned int iseed = 0; iseed < myEcSeeds.size(); iseed++)
      {
          EcalSeed i_seed = myEcSeeds.at(iseed);
          float i_theta = i_seed.GetTheta();
          float i_phi   = i_seed.GetPhi(); 
          bool maxIsolatedHit = true;
                    
          for (long unsigned int jseed = 0; jseed < myEcSeeds.size(); jseed++)
          {
              EcalSeed j_seed = myEcSeeds.at(jseed);
              float j_theta = j_seed.GetTheta();
              float j_phi   = j_seed.GetPhi();
              float dd = sqrt(pow(j_theta-i_theta,2) + pow(j_phi-i_phi,2));
              
              if (dd < maxDeltaR)
              {
                  std::cout << " dd = " << dd << " :: j_theta-i_theta= " << j_theta-i_theta << " :: j_phi-i_phi " << j_phi-i_phi << std::endl;
                  std::cout << "Removing neighboring seed with lowest energy" << std::endl;
                  if (i_seed.GetEne()<j_seed.GetEne()) maxIsolatedHit = false;
              }
          }
          
          if (maxIsolatedHit) myEcSeedsCleaned.push_back(i_seed);
      }
      
    
    
      std::cout << "Matching clusters with gen level" << std::endl;
      for (long unsigned int iseed = 0; iseed < myEcSeedsCleaned.size(); iseed++)
      { 
          EcalSeed this_seed = myEcSeedsCleaned.at(iseed);
          float seed_theta = this_seed.GetTheta();
          float seed_phi   = this_seed.GetPhi();
          
          
          TruthTree->GetEntry(selEv);
          for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
          {
              
              int    pdgId = myTruthTV.mcs_pdgId->at(i);
              double ene   = myTruthTV.mcs_E->at(i);
              if (ene<0.1) continue;
              
              double truth_phi   = myTruthTV.mcs_phi->at(i);
              double eta   = myTruthTV.mcs_eta->at(i);
              int charge = myTruthTV.mcs_charge->at(i);              
              double truth_theta = 2*atan(exp(-eta));
              truth_theta = M_PI- truth_theta;
              
              
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaR)
              {
                  std::cout  << "ECAL cluster " << iseed << "matched to MC truth gen level particle: " << pdgId << std::endl;
                  this_seed.AddGenMatch(pdgId);
              }
            }                    
      }
      
      
      std::cout << "Matching gen level with clusters" << std::endl;
      
      TruthTree->GetEntry(selEv);
      int nGenMatched = 0;
      for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
      {
          int    pdgId = myTruthTV.mcs_pdgId->at(i);
          double ene   = myTruthTV.mcs_E->at(i);
          if (ene<0.1) continue;
          
          int matchedToCluster = 0;
          
          double truth_phi   = myTruthTV.mcs_phi->at(i);
          double eta   = myTruthTV.mcs_eta->at(i);
          int charge = myTruthTV.mcs_charge->at(i);              
          double truth_theta = 2*atan(exp(-eta));
          truth_theta = M_PI- truth_theta;
          
          
          for (long unsigned int iseed = 0; iseed < myEcSeedsCleaned.size(); iseed++)
          { 
              EcalSeed this_seed = myEcSeedsCleaned.at(iseed);
              float seed_theta = this_seed.GetTheta();
              float seed_phi   = this_seed.GetPhi();
            
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaR)
              {                  
                  std::cout << "MC truth gen level particle: " << pdgId << " matched to ECAL cluster " << iseed << std::endl;
                  matchedToCluster ++;
              }
            }
            
            if (matchedToCluster==0) std::cout << "MC truth gen level particle: " << pdgId << " NOT matched to any ECAL cluster " << std::endl;
            else                     nGenMatched++;
      }
      
      std::cout << "Fraction of gen particles matched to at least one ECAL cluster: " << float(nGenMatched)/float(myTruthTV.mcs_E->size()) << std::endl;
          
      
      
      
      //**************************************************************//
      //                            Timing
      //**************************************************************//
      
//       double c_speed = 1./299792458*1e9; //mm per picosecond          
//       double time_acceptance = 3; //time window to reject out of time hits
      float ene_T_th = 1;
      float totTimingEne = 0;
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_F->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_F->at(i),1);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepF->at(i);          
          if (this_ene>ene_T_th) hGrid_T1 ->Fill(this_theta, this_phi,this_ene);                                             
          totTimingEne+= this_ene;
      }
      
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_R->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_R->at(i),2);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepR->at(i);
              
              
          if (this_ene>ene_T_th) hGrid_T2 ->Fill(this_theta, this_phi,this_ene);                                             
          totTimingEne+= this_ene;
      }
      
      std::cout << "Total energy in TIMING: " << totTimingEne << std::endl;

  
  
//   float fit_range = 0.01;
//   float phi_res_b, phi_res_e, phi_res_b_cg, phi_res_e_cg;
//   float eta_res_b, eta_res_e, eta_res_b_cg, eta_res_e_cg;

  /// other plots
  
  TLine * lEndcapMinus = new TLine(M_PI/4*1, minPhi, M_PI/4*1, maxPhi);
  TLine * lEndcapPlus  = new TLine(M_PI/4*3, minPhi, M_PI/4*3, maxPhi);
  lEndcapMinus->SetLineColor(kRed);
  lEndcapMinus->SetLineWidth(2);
  lEndcapPlus->SetLineColor(kRed);
  lEndcapPlus->SetLineWidth(2);
  
//   float phiRange = (maxPhi-minPhi)/2;
//   float thetaRange = M_PI/8;
  
/*  
  TCanvas * cGrid_DRT_S = new TCanvas("cGrid_DRT_S", "cGrid_DRT_S", 900, 600);
  cGrid_DRT_S->cd();
  hGrid_DRT_S->Draw("LEGO2Z");
  hGrid_DRT_S->SetStats(0);
  hGrid_DRT_S->SetTitle("Dual Readout HCAL Tower");
  hGrid_DRT_S->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_DRT_S->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_DRT_S->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_DRT_S->GetYaxis()->SetRangeUser(minPhi, maxPhi);
//   lEndcapMinus->Draw("same");
//   lEndcapPlus->Draw("same");
//   gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_DRT_S->SaveAs(Form("plots/cGrid_DRT_S_%s.png", output_tag.c_str()));
  

  TCanvas * cGrid_DRT_C = new TCanvas("cGrid_DRT_C", "cGrid_DRT_C", 900, 600);
  cGrid_DRT_C->cd();
  hGrid_DRT_C->Draw("LEGO2Z");
  hGrid_DRT_C->SetStats(0);
  hGrid_DRT_C->SetTitle("Dual Readout HCAL Tower");
  hGrid_DRT_C->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_DRT_C->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_DRT_C->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_DRT_C->GetYaxis()->SetRangeUser(minPhi, maxPhi);
//   lEndcapMinus->Draw("same");
//   lEndcapPlus->Draw("same");
//   gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_DRT_C->SaveAs(Form("plots/cGrid_DRT_C_%s.png", output_tag.c_str()));
  
  

  
     TCanvas * cGrid_EC_F = new TCanvas("cGrid_EC_F", "cGrid_EC_F", 900, 600);
  cGrid_EC_F->cd();
  hGrid_EC_F->Draw("LEGO2Z");
  hGrid_EC_F->SetStats(0);
  hGrid_EC_F->SetTitle("E1");
  hGrid_EC_F->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_F->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_F->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_EC_F->GetYaxis()->SetRangeUser(minPhi, maxPhi);
//   lEndcapMinus->Draw("same");
//   lEndcapPlus->Draw("same");
//   gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_EC_F->SaveAs(Form("plots/cGrid_EC_F_%s.png", output_tag.c_str()));
  
  
    TCanvas * cGrid_EC_R = new TCanvas("cGrid_EC_R", "cGrid_EC_R", 900, 600);
  cGrid_EC_R->cd();
  hGrid_EC_R->Draw("LEGO2Z");
  hGrid_EC_R->SetStats(0);
  hGrid_EC_R->SetTitle("E2");
  hGrid_EC_R->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_R->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_R->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_EC_R->GetYaxis()->SetRangeUser(minPhi, maxPhi);
//   lEndcapMinus->Draw("same");
//   lEndcapPlus->Draw("same");
//   gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_EC_R->SaveAs(Form("plots/cGrid_EC_R_%s.png", output_tag.c_str()));
  
  
    TCanvas * cGrid_T1 = new TCanvas("cGrid_T1", "cGrid_T1", 900, 600);
  cGrid_T1->cd();
  hGrid_T1->Draw("LEGO2Z");
  hGrid_T1->SetStats(0);
  hGrid_T1->SetTitle("T1");
  hGrid_T1->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_T1->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_T1->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_T1->GetYaxis()->SetRangeUser(minPhi, maxPhi);
//   lEndcapMinus->Draw("same");
//   lEndcapPlus->Draw("same");
//   gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_T1->SaveAs(Form("plots/cGrid_T1_%s.png", output_tag.c_str()));
  
  
    TCanvas * cGrid_T2 = new TCanvas("cGrid_T2", "cGrid_T2", 900, 600);
  cGrid_T2->cd();
  hGrid_T2->Draw("LEGO2Z");
  hGrid_T2->SetStats(0);
  hGrid_T2->SetTitle("T2");
  hGrid_T2->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_T2->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_T2->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_T2->GetYaxis()->SetRangeUser(minPhi, maxPhi);
//   lEndcapMinus->Draw("same");
//   lEndcapPlus->Draw("same");
//   gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_T2->SaveAs(Form("plots/cGrid_T2_%s.png", output_tag.c_str()));
  */

    TCanvas * cGrid_EC_T = new TCanvas("cGrid_EC_T", "cGrid_EC_T", 900, 1600);
  cGrid_EC_T->cd();
  hGrid_EC_T->Draw("BOX");
  hGrid_EC_T->SetStats(0);
  hGrid_EC_T->SetTitle("E1+E2");
  hGrid_EC_T->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_T->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_T->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_EC_T->GetYaxis()->SetRangeUser(minPhi, maxPhi);
      
  for (long unsigned int iseed = 0; iseed < myEcSeedsCleaned.size(); iseed++)
  {          
      EcalSeed this_seed = myEcSeedsCleaned.at(iseed);
      TEllipse * el1 = new TEllipse(this_seed.GetTheta(),this_seed.GetPhi(), maxDeltaR, maxDeltaR);
      std::cout << "seed: " << iseed << " :: theta = " << this_seed.GetTheta() <<  ", phi = " << this_seed.GetPhi() << std::endl;
      el1->SetLineColor(kRed);
      el1->SetFillStyle(0);
      el1->SetLineWidth(2);
      el1->Draw();
  }
  if (SAVEPLOTS)   cGrid_EC_T->SaveAs(Form("plots/cGrid_EC_T_%s.png", output_tag.c_str()));
  
  
  
  THStack *hStackedTruth = new THStack("hStackedTruth","hStackedTruth");
  
  hTruthFloor->SetFillColor(kWhite);
  hTruthFloor->SetLineColor(kWhite);
  
  hTruthChargedHAD->SetFillColor(kRed+1);
  hTruthChargedHAD->SetLineColor(kRed+1);
  hTruthChargedEM->SetFillColor(kYellow+2);
  hTruthChargedEM->SetLineColor(kYellow+2);
  
  hTruthNeutralHAD->SetFillColor(kBlue);
  hTruthNeutralHAD->SetLineColor(kBlue);
  hTruthNeutralEM->SetFillColor(kGreen+1);
  hTruthNeutralEM->SetLineColor(kGreen+1);
  
  hStackedTruth->Add(hTruthFloor);
  hStackedTruth->Add(hTruthChargedEM);
  hStackedTruth->Add(hTruthNeutralEM);
  hStackedTruth->Add(hTruthChargedHAD);
  hStackedTruth->Add(hTruthNeutralHAD);
     
//   TCanvas * cGrid_Truth = new TCanvas("cGrid_Truth", "cGrid_Truth", 900, 600);
//   cGrid_Truth->cd();
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");
  leg->AddEntry(hTruthChargedEM, "Electrons", "lpf");
  leg->AddEntry(hTruthNeutralEM, "Photons", "lpf");
  leg->AddEntry(hTruthChargedHAD, "Charged (except e^{-})", "lpf");
  leg->AddEntry(hTruthNeutralHAD, "Neutrals (except #gamma)", "lpf");
//   hGrid_EC_R->Draw("LEGO2Z");
//   hStackedTruth->SetStats(0);
//   hStackedTruth->SetTitle("Truth");
//   hStackedTruth->GetXaxis()->SetTitle("#theta [rad]");
//   hStackedTruth->GetYaxis()->SetTitle("#phi [rad]");
//   hStackedTruth->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hStackedTruth->GetYaxis()->SetRangeUser(minPhi, maxPhi);
//   hStackedTruth->Draw();
//   leg->Draw();
  
  
  
  
    
  TCanvas * cOverview = new TCanvas ("cOverview", "cOverview", 700, 1400);
  cOverview->Divide(1,4);
  
  cOverview->cd(1);
  hGrid_DRT_S->Draw("BOX");
  hGrid_DRT_S->SetStats(0);
  hGrid_DRT_S->SetTitle("Dual Readout HCAL Tower");
  hGrid_DRT_S->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_DRT_S->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_DRT_S->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_DRT_S->GetYaxis()->SetRangeUser(minPhi, maxPhi);
  
  cOverview->cd(2);
  hGrid_EC_T->Draw("BOX");
  hGrid_EC_T->SetStats(0);
  hGrid_EC_T->SetTitle("E1+E2 [GeV]");
  hGrid_EC_T->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_T->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_T->GetZaxis()->SetTitle("Energy [GeV]");
  hGrid_EC_T->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_EC_T->GetYaxis()->SetRangeUser(minPhi, maxPhi);
  

        
  
  cOverview->cd(3);
  hGrid_T2->Draw("BOX");
  hGrid_T2->SetStats(0);
  hGrid_T2->SetTitle("T2 [MeV]");
  hGrid_T2->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_T2->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_T2->GetZaxis()->SetTitle("Energy [MeV]");
  hGrid_T2->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_T2->GetYaxis()->SetRangeUser(minPhi, maxPhi);
  
  cOverview->cd(4);
  hStackedTruth->Draw("BOX");
  leg->Draw();
  
//   gPad->SetLogz();
  cOverview->SaveAs(Form("plots/event_display/cOverview_%d.png", selEv));
  
//   hGrid_TT->Reset();
//   hGrid_EC_F->Reset();
//   hGrid_EC_R->Reset();       
//   hGrid_EC_T->Reset();       
//   hGrid_DRT_S->Reset();
//   hGrid_DRT_C->Reset();
//   
//   
//   
  
  
  
  
  theApp->Run();
}







