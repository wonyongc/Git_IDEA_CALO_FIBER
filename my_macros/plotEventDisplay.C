// g++ -Wall -o plotEventDisplay plotEventDisplay.C myG4Tree.cc myG4Tree.hh myTruthTree.cc myTruthTree.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs`


#include "SCEPCal_GeometryHelper.hh"
#include "myG4Tree.hh"
#include "myTruthTree.hh"
// #include "ppfa.hh"

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
#include "TObject.h"
#include "TVector3.h"



int main(int argc, char** argv)
{

//   TApplication* theApp = new TApplication("App", &argc, argv);
      
  using namespace std;
  
  gStyle->SetTitleXOffset (1.00) ;                                                                                        
  gStyle->SetTitleYOffset (1.2) ;                                                                                                                                                                                                                 
  gStyle->SetPadLeftMargin (0.12) ;
  gStyle->SetPadRightMargin (0.16) ;
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
  
  std::string particle_name = "K^{0L}, 0T";
//   std::string output_tag = "kaon_central";
  std::string output_tag = "muon_iso";
  int selEv = 0;
  if (argc > 1) selEv = atoi(argv[1]);
    
  //define histos
  
  int NPHI_TL1   = 186-1;
  int NTHETA_TL1 = 29*16-1;
  int NPHI_TL2   = 186*16-1;
  int NTHETA_TL2 = 29-1;
  
//   int NPHI_EC    = 1130;
//   int NTHETA_EC  = 180+180-1;
  int NPHI_EC    = 1130;
  int NTHETA_EC  = (180+180-1);
  int NPHI_DRT   = 252;
  int NTHETA_DRT = 40+40-1;
  
  double drh_S_norm  = 407;
//   double drh_C_norm  = 103.2;  
  int phiGran = NPHI_DRT;
  
  
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
  TH2F * hTruthMuon = new TH2F ("hTruthMuon", "hTruthMuon", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hTruthNeutralEM = new TH2F ("hTruthNeutralEM", "hTruthNeutralEM", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hTruthChargedHAD = new TH2F ("hTruthChargedHAD", "hTruthChargedHAD", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hTruthNeutralHAD = new TH2F ("hTruthNeutralHAD", "hTruthNeutralHAD", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  
  
  int NRADIUS = 1;
  float minRadiusT1 = 1750;
  float maxRadiusT1 = 1775;
  float minRadiusT2 = 1775;
  float maxRadiusT2 = 1800;
  
  float minRadiusF = 1800;
  float maxRadiusF = 1860;  
  float minRadiusR = 1860;
  float maxRadiusR = 2000;
  
  float minRadius = 1750;
  float maxRadius = 2000;
  
  float minRadius_DRT = 2500;
  float maxRadius_DRT = 4500;
//   float maxRadius_DRT = 2700;
  
  TH2D * hPolar_T1 = new TH2D ("hPolar_T1", "hPolar_T1", NPHI_TL1, minPhi-bin_width_phi_TL1/2, maxPhi-bin_width_phi_TL1/2, NRADIUS, minRadiusT1, maxRadiusT1);
  TH2D * hPolar_T2 = new TH2D ("hPolar_T2", "hPolar_T2", NPHI_TL2, minPhi-bin_width_phi_TL2/2, maxPhi-bin_width_phi_TL2/2, NRADIUS, minRadiusT2, maxRadiusT2);
  
  TH2D * hPolar_EC_F = new TH2D ("hPolar_EC_F", "hPolar_EC_F", NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2, NRADIUS, minRadiusF, maxRadiusF);
  TH2D * hPolar_EC_R = new TH2D ("hPolar_EC_R", "hPolar_EC_R", NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2, NRADIUS, minRadiusR, maxRadiusR);
  TH2D * hPolar_EC_T = new TH2D ("hPolar_EC_T", "hPolar_EC_T", NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2, NRADIUS, minRadius, maxRadius);
  
  TH2D * hPolar_DRT_S = new TH2D ("hPolar_DRT_S", "hPolar_DRT_S", NPHI_DRT, minPhi-bin_width_phi_DRT/2, maxPhi-bin_width_phi_DRT/2, NRADIUS, minRadius_DRT, maxRadius_DRT);
  
  
  SCEPCal_GeometryHelper myGeometry;

  
  //run over energy scan
//   TFile * RunFile = new TFile("../root_files/iso_gun/iso_gun_mu_100GeV_T+E.root","READ"); 
//   TFile * RunFile = new TFile("../root_files/iso_gun/central_kaon0L_10GeV_ALL.root","READ"); 
  
//   TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_gamma_Iso+Uniform1-100_GeV.root","READ");       
//   TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_e-_Iso+Uniform1-100_GeV.root","READ");         
//     TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_kaon0L_Iso+Uniform1-100_GeV.root","READ");       
  
  
//   TFile * RecoFile = new TFile("../root_files/hep_outputs/output_hep_test.root","READ");       
//   TFile * RecoFile = new TFile("../root_files/hep_outputs/output_SCEPCal_B2T_HG_wwlj100k_job_15.root","READ");
//   TFile * RecoFile = new TFile("../root_files/hep_outputs/output_SCEPCal_B0T_zjj_scan_100_job_0.root","READ");       
//   TFile * RecoFile = new TFile("../root_files/hep_outputs/output_SCEPCal_B0T_HG_zjj_scan_90_job_0.root","READ");       
  TFile * RecoFile = new TFile("../root_files/hep_outputs/output_SCEPCal_B2T_HG_zjj_scan_90_job_0.root","READ");       
  
//   TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_Iso+Uniform_B2T__pi-_Iso+Uniform1-100_GeV_job_0.root","READ");       
//   TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_Iso+Uniform_B2T__mu-_Iso+Uniform1-100_GeV_job_0.root","READ");       
//   TFile * RecoFile = new TFile("../root_files/prod/output_SCEPCal_Iso+Uniform_B2T__mu+_Iso+Uniform1-100_GeV_job_0.root","READ");       
  

  
  TTree* TreeRun = (TTree*) RecoFile->Get("B4");
  myG4TreeVars myTV;
  InitG4Tree (TreeRun, myTV);
      
      
  bool isHepMC = true;

  ///*******************************************///
  ///		 Run over events	    ///
  ///*******************************************///
    
  int NEVENTS = TreeRun->GetEntries();
  std::cout << "NEVENTS = " << NEVENTS << std::endl;
  
//   std::vector<TLine*> charged_tracks;
//   std::vector<TLine*> muon_tracks;
  
  std::vector<TGraph*> charged_tracks;
  std::vector<TGraph*> muon_tracks;
  
  std::vector<TLine*> photon_tracks;
  std::vector<TLine*> neutrhad_tracks;
  
  float max_radius = 1740;
  float max_muon_radius = maxRadius_DRT*1.2;
  
  float MC_ene_th = 0.3;
  float EC_hit_th = 0.05;
  float HC_hit_th = 0.1;
  float TT_hit_th = 0.002;
  
  float Bfield = 2.;
  
  if (isHepMC)
  {
//       TFile * TruthFile = new TFile("../root_files/hep_outputs/hep_truth.root","READ");
//       TFile * TruthFile = new TFile("../../HepMC_Files/wwlj100k_job_12_output_tuple.root","READ");
//       TFile * TruthFile = new TFile("../../HepMC_Files/B0T/zjj_scan_90_job_0_output_tuple.root","READ");
      TFile * TruthFile = new TFile("../../HepMC_Files/B2T/zjj_scan_90_job_0_output_tuple.root","READ");
//       TFile * TruthFile = new TFile("../../HepMC_Files/B2T/wwlj100k_job_15_output_tuple.root","READ");
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
        double pT    = myTruthTV.mcs_pt->at(i);
        int charge = myTruthTV.mcs_charge->at(i);
        double theta = 2*atan(exp(-eta));
        theta = M_PI- theta;
        
        double px, py;
        px = pT*cos(phi);
        py = pT*sin(phi);              
        double pz = -pT*sinh(eta);                    
        
        
        if      (charge!=0 && abs(pdgId)!= 11 && abs(pdgId)!= 13) 
        {
            hTruthChargedHAD ->Fill(theta, phi, ene);
            std::cout << "this is a non electron charged particle!" << std::endl;
        }
        if (abs(pdgId)== 11) 
        {
            hTruthChargedEM  ->Fill(theta, phi, ene);
            std::cout << "this is an electron!" << std::endl;
        }
        if (charge==0 && abs(pdgId)!= 22 && abs(pdgId)!= 12 && abs(pdgId)!= 14 && abs(pdgId)!= 16 && abs(pdgId)< 10000) 
        {
            hTruthNeutralHAD ->Fill(theta, phi, ene);
            std::cout << "this is a  neutral hadron!" << std::endl;
        }
        if (charge==0 && pdgId== 22) 
        {
            std::cout << "this is a photon!" << std::endl;
            hTruthNeutralEM  ->Fill(theta, phi, ene);
        }
        
        if (abs(pdgId)== 13) 
        {
            std::cout << "this is a muon!" << std::endl;
            hTruthMuon  ->Fill(theta, phi, ene);
        }
        
        
        if (ene < MC_ene_th) continue;
        
        //define trajectories
        if (abs(pdgId)== 22)
        {
            TLine* new_neutral_line = new TLine (0, 0, cos(phi)*max_radius, sin(phi)*max_radius) ;
            photon_tracks.push_back(new_neutral_line);
        }
        else if (fabs(pdgId) == 130 || fabs(pdgId) == 2112)
        {
            TLine* new_neutral_line = new TLine (0, 0, cos(phi)*max_radius, sin(phi)*max_radius) ;
            neutrhad_tracks.push_back(new_neutral_line);
        }
            
/*        if (Bfield==0)
        {
            if (charge!=0)
            {                
                if (abs(pdgId)!= 13)
                {                            
                    TLine* new_charged_line = new TLine (0, 0, cos(phi)*max_radius, sin(phi)*max_radius) ;
                    charged_tracks.push_back(new_charged_line);
                }
                if (abs(pdgId)== 13)
                {
                    TLine* new_charged_line = new TLine (0, 0, cos(phi)*max_muon_radius, sin(phi)*max_muon_radius) ;
                    muon_tracks.push_back(new_charged_line);
                }
            }            
        }
        else */
//         if (Bfield=0)
        {
            if (charge!=0)
            {                
                if (abs(pdgId)!= 13)
                {
                    charged_tracks.push_back(getTrajectory(Bfield, px, py, pz, charge, max_radius));
                }
                if (abs(pdgId)== 13)
                {
                    muon_tracks.push_back(getTrajectory(Bfield, px, py, pz, charge, max_radius));
//                     muon_tracks.push_back(getTrajectory(0, px, py, pz, charge, max_radius));
                }
            }
            
        }
        
        std::cout << "******************************************************* " << std::endl;
        
      }
      
            

      
  }
  

//   NEVENTS = 1000;
  
//   for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
      int iEvt = selEv;
                                        
     
      TreeRun->GetEntry(iEvt);
      std::cout << "processing event: " << iEvt << "\r" << std::flush;
      if (!isHepMC)
      {  
        double px  = myTV.PrimaryParticleMomentum->at(0)/1000.;
        double py  = myTV.PrimaryParticleMomentum->at(1)/1000.;
        double pz  = myTV.PrimaryParticleMomentum->at(2)/1000.;
        double P   = sqrt(px*px+py*py+pz*pz);
        px/= P;
        py/= P;
        pz/= P;

      
        double phi   = atan(py/px);
        if (px<0. && py <0.)   {phi = phi - M_PI;}
        if (px<0. && py >0.)   {phi = M_PI + phi;}
        double eta   = -atanh(pz);      
        double theta = 2*atan(exp(-eta));
        std::cout << "Single primary particle with momentum and direction:" << std::endl;  
        std::cout << "input P = " << P << " GeV :: phi = " << phi << " :: theta = " << theta << std::endl;
        
        hTruthChargedEM  ->Fill(theta, phi, P);
        
        muon_tracks.push_back(getTrajectory(Bfield, px*P, py*P, pz*P, -1, max_radius));
//         std::cout << "truth theta = " << theta << " :: truth phi = " << phi << std::endl;
      }
          
      //**************************************************************//
      //                           DR HCAL
      //**************************************************************//

      float totDRHScint = 0;
      float totDRHCher  = 0;
      
      for (unsigned int i = 0; i<myTV.VectorL->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'l', phiGran);
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorL->at(i)/1000.;      
          double this_scint = myTV.VectorSignalsL->at(i);                
          double this_cher  = myTV.VectorSignalsCherL->at(i);
          if (this_ene<HC_hit_th) continue;
          hGrid_DRT_S ->Fill(this_theta, this_phi, this_scint);                              
          hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher);    
          hPolar_DRT_S ->Fill(this_phi, minRadius_DRT, this_scint/drh_S_norm);
          totDRHScint+=this_scint;
          totDRHCher+=this_cher;
          
//           std::cout << "hcal reco theta = " << this_theta << " :: hcal reco phi = " << this_phi << std::endl;
      }
      for (unsigned int i = 0; i<myTV.VectorR->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'r', phiGran);
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorR->at(i)/1000.;                    
          double this_scint = myTV.VectorSignalsR->at(i);     
          double this_cher  = myTV.VectorSignalsCherR->at(i);      
          if (this_ene<HC_hit_th) continue;
          
          hGrid_DRT_S ->Fill(this_theta, this_phi, this_scint);                              
          hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher);           
          hPolar_DRT_S ->Fill(this_phi, minRadius_DRT, this_scint/drh_S_norm);
          totDRHScint+=this_scint;
          totDRHCher+=this_cher;
          
//           std::cout << "hcal reco theta = " << this_theta << " :: hcal reco phi = " << this_phi << std::endl;
      }
      
      std::cout << "Total scint in HCAL: " << totDRHScint << std::endl;
      std::cout << "Total cher in HCAL: " << totDRHCher << std::endl;
      
      //**************************************************************//
      //                             ECAL
      //**************************************************************//
      
      float totEcalEne = 0;
      for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
      {                            
              
          TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;                    
          
//           if (fabs(myTV.VecHit_CrystalID->at(i)) <1000000)
          if (this_ene>EC_hit_th)
          {
            hGrid_EC_F ->Fill(this_theta, this_phi, myTV.VecHit_ScepEneDepF->at(i)/1000);
            hGrid_EC_R ->Fill(this_theta, this_phi, myTV.VecHit_ScepEneDepR->at(i)/1000);
            hGrid_EC_T ->Fill(this_theta, this_phi, this_ene);
            
            hPolar_EC_F ->Fill(this_phi, minRadiusF, myTV.VecHit_ScepEneDepF->at(i)/1000);
            hPolar_EC_R ->Fill(this_phi, minRadiusR, myTV.VecHit_ScepEneDepR->at(i)/1000);
            
            
            totEcalEne+=this_ene;
            
//             std::cout << "ecal reco theta = " << this_theta << " :: ecal reco phi = " << this_phi << std::endl;
          }
      }
      
      std::cout << "Total energy in ECAL: " << totEcalEne << std::endl;

      
      
      //**************************************************************//
      //                            Timing
      //**************************************************************//
      
//       double c_speed = 1./299792458*1e9; //mm per picosecond          
//       double time_acceptance = 3; //time window to reject out of time hits
      
      float totTimingEne = 0;
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_F->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_F->at(i),1);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepF->at(i)/1000;          
          if (this_ene<TT_hit_th) continue;
          hGrid_T1 ->Fill(this_theta, this_phi,this_ene);                                             
          totTimingEne+= this_ene;
          hPolar_T1 ->Fill(this_phi, minRadiusT1, this_ene);
        
      }
      
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_R->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_R->at(i),2);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepR->at(i)/1000;
          if (this_ene<TT_hit_th) continue;
          hGrid_T2 ->Fill(this_theta, this_phi,this_ene);                                             
          totTimingEne+= this_ene;
          hPolar_T2 ->Fill(this_phi, minRadiusT2, this_ene);
      }
      
      std::cout << "Total energy in TIMING: " << totTimingEne << std::endl;
      
  }
  
  
//   float fit_range = 0.01;
//   float phi_res_b, phi_res_e, phi_res_b_cg, phi_res_e_cg;
//   float eta_res_b, eta_res_e, eta_res_b_cg, eta_res_e_cg;

/*  
  TCanvas * cDeltaTheta_TimingEndcap = new TCanvas ("cDeltaTheta_TimingEndcap", "cDeltaTheta_TimingEndcap", 600, 500);
  cDeltaTheta_TimingEndcap->cd();
  hDeltaTheta_TimingEndcap_F[NFILES-1]->Draw();
  hDeltaTheta_TimingEndcap_F[NFILES-1]->SetStats(0);
  hDeltaTheta_TimingEndcap_F[NFILES-1]->SetTitle(Form("Endcap: %d GeV %s",energies[0], particle_name.c_str()));
  hDeltaTheta_TimingEndcap_F[NFILES-1]->GetXaxis()->SetRangeUser(-fit_range, fit_range);
  hDeltaTheta_TimingEndcap_F[NFILES-1]->GetYaxis()->SetRangeUser(0,hDeltaTheta_TimingEndcap_Comb[NFILES-1]->GetMaximum()*1.5);
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
  if (SAVEPLOTS) cDeltaTheta_TimingEndcap->SaveAs("plots/cDeltaTheta_TimingEndcap.png");
  
  
  
  */
  /// other plots
  
  TLine * lEndcapMinus = new TLine(M_PI/4*1, minPhi, M_PI/4*1, maxPhi);
  TLine * lEndcapPlus  = new TLine(M_PI/4*3, minPhi, M_PI/4*3, maxPhi);
  lEndcapMinus->SetLineColor(kRed);
  lEndcapMinus->SetLineWidth(2);
  lEndcapPlus->SetLineColor(kRed);
  lEndcapPlus->SetLineWidth(2);
  
//   float phiRange = (maxPhi-minPhi)/2;
//   float thetaRange = M_PI/8;
  
  
//   TCanvas * cGrid_DRT_S = new TCanvas("cGrid_DRT_S", "cGrid_DRT_S", 900, 600);
//   cGrid_DRT_S->cd();
//   hGrid_DRT_S->Draw("LEGO2Z");
//   hGrid_DRT_S->SetStats(0);
//   hGrid_DRT_S->SetTitle("Dual Readout HCAL Tower");
//   hGrid_DRT_S->GetXaxis()->SetTitle("#theta [rad]");
//   hGrid_DRT_S->GetYaxis()->SetTitle("#phi [rad]");
//   hGrid_DRT_S->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hGrid_DRT_S->GetYaxis()->SetRangeUser(minPhi, maxPhi);
// //   lEndcapMinus->Draw("same");
// //   lEndcapPlus->Draw("same");
// //   gPad->SetLogz();
//   if (SAVEPLOTS)   cGrid_DRT_S->SaveAs(Form("plots/cGrid_DRT_S_%s.png", output_tag.c_str()));
//   
// 
//   TCanvas * cGrid_DRT_C = new TCanvas("cGrid_DRT_C", "cGrid_DRT_C", 900, 600);
//   cGrid_DRT_C->cd();
//   hGrid_DRT_C->Draw("LEGO2Z");
//   hGrid_DRT_C->SetStats(0);
//   hGrid_DRT_C->SetTitle("Dual Readout HCAL Tower");
//   hGrid_DRT_C->GetXaxis()->SetTitle("#theta [rad]");
//   hGrid_DRT_C->GetYaxis()->SetTitle("#phi [rad]");
//   hGrid_DRT_C->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hGrid_DRT_C->GetYaxis()->SetRangeUser(minPhi, maxPhi);
// //   lEndcapMinus->Draw("same");
// //   lEndcapPlus->Draw("same");
// //   gPad->SetLogz();
//   if (SAVEPLOTS)   cGrid_DRT_C->SaveAs(Form("plots/cGrid_DRT_C_%s.png", output_tag.c_str()));
//   
//   
//   TCanvas * cGrid_EC_T = new TCanvas("cGrid_EC_T", "cGrid_EC_T", 900, 600);
//   cGrid_EC_T->cd();
//   hGrid_EC_T->Draw("LEGO2Z");
//   hGrid_EC_T->SetStats(0);
//   hGrid_EC_T->SetTitle("E1+E2");
//   hGrid_EC_T->GetXaxis()->SetTitle("#theta [rad]");
//   hGrid_EC_T->GetYaxis()->SetTitle("#phi [rad]");
//   hGrid_EC_T->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hGrid_EC_T->GetYaxis()->SetRangeUser(minPhi, maxPhi);
// //   lEndcapMinus->Draw("same");
// //   lEndcapPlus->Draw("same");
// //   gPad->SetLogz();
//   if (SAVEPLOTS)   cGrid_EC_T->SaveAs(Form("plots/cGrid_EC_T_%s.png", output_tag.c_str()));
//   
  
//      TCanvas * cGrid_EC_F = new TCanvas("cGrid_EC_F", "cGrid_EC_F", 900, 600);
//   cGrid_EC_F->cd();
//   hGrid_EC_F->Draw("LEGO2Z");
//   hGrid_EC_F->SetStats(0);
//   hGrid_EC_F->SetTitle("E1");
//   hGrid_EC_F->GetXaxis()->SetTitle("#theta [rad]");
//   hGrid_EC_F->GetYaxis()->SetTitle("#phi [rad]");
//   hGrid_EC_F->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hGrid_EC_F->GetYaxis()->SetRangeUser(minPhi, maxPhi);
// //   lEndcapMinus->Draw("same");
// //   lEndcapPlus->Draw("same");
// //   gPad->SetLogz();
//   if (SAVEPLOTS)   cGrid_EC_F->SaveAs(Form("plots/cGrid_EC_F_%s.png", output_tag.c_str()));
//   
//   
//     TCanvas * cGrid_EC_R = new TCanvas("cGrid_EC_R", "cGrid_EC_R", 900, 600);
//   cGrid_EC_R->cd();
//   hGrid_EC_R->Draw("LEGO2Z");
//   hGrid_EC_R->SetStats(0);
//   hGrid_EC_R->SetTitle("E2");
//   hGrid_EC_R->GetXaxis()->SetTitle("#theta [rad]");
//   hGrid_EC_R->GetYaxis()->SetTitle("#phi [rad]");
//   hGrid_EC_R->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hGrid_EC_R->GetYaxis()->SetRangeUser(minPhi, maxPhi);
// //   lEndcapMinus->Draw("same");
// //   lEndcapPlus->Draw("same");
// //   gPad->SetLogz();
//   if (SAVEPLOTS)   cGrid_EC_R->SaveAs(Form("plots/cGrid_EC_R_%s.png", output_tag.c_str()));
//   
//   
//     TCanvas * cGrid_T1 = new TCanvas("cGrid_T1", "cGrid_T1", 900, 600);
//   cGrid_T1->cd();
//   hGrid_T1->Draw("LEGO2Z");
//   hGrid_T1->SetStats(0);
//   hGrid_T1->SetTitle("T1");
//   hGrid_T1->GetXaxis()->SetTitle("#theta [rad]");
//   hGrid_T1->GetYaxis()->SetTitle("#phi [rad]");
//   hGrid_T1->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hGrid_T1->GetYaxis()->SetRangeUser(minPhi, maxPhi);
// //   lEndcapMinus->Draw("same");
// //   lEndcapPlus->Draw("same");
// //   gPad->SetLogz();
//   if (SAVEPLOTS)   cGrid_T1->SaveAs(Form("plots/cGrid_T1_%s.png", output_tag.c_str()));
//   
//   
//     TCanvas * cGrid_T2 = new TCanvas("cGrid_T2", "cGrid_T2", 900, 600);
//   cGrid_T2->cd();
//   hGrid_T2->Draw("LEGO2Z");
//   hGrid_T2->SetStats(0);
//   hGrid_T2->SetTitle("T2");
//   hGrid_T2->GetXaxis()->SetTitle("#theta [rad]");
//   hGrid_T2->GetYaxis()->SetTitle("#phi [rad]");
//   hGrid_T2->GetXaxis()->SetRangeUser(minTheta, maxTheta);
//   hGrid_T2->GetYaxis()->SetRangeUser(minPhi, maxPhi);
// //   lEndcapMinus->Draw("same");
// //   lEndcapPlus->Draw("same");
// //   gPad->SetLogz();
//   if (SAVEPLOTS)   cGrid_T2->SaveAs(Form("plots/cGrid_T2_%s.png", output_tag.c_str()));
//   

  
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
  
  hTruthMuon->SetFillColor(kGray);
  hTruthMuon->SetLineColor(kGray);
  
  hStackedTruth->Add(hTruthFloor);
  hStackedTruth->Add(hTruthChargedEM);
  hStackedTruth->Add(hTruthNeutralEM);
  hStackedTruth->Add(hTruthChargedHAD);
  hStackedTruth->Add(hTruthNeutralHAD);
//   hStackedTruth->Add(hTruthMuon);
     
//   TCanvas * cGrid_Truth = new TCanvas("cGrid_Truth", "cGrid_Truth", 900, 600);
//   cGrid_Truth->cd();
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");
  leg->AddEntry(hTruthChargedEM, "Electrons", "lpf");
  leg->AddEntry(hTruthNeutralEM, "Photons", "lpf");
  leg->AddEntry(hTruthChargedHAD, "Charged (except e^{-}, #mu^{-})", "lpf");
  leg->AddEntry(hTruthNeutralHAD, "Neutral hadron", "lpf");
//   leg->AddEntry(hTruthMuon, "Muon", "lpf");
// 
//   hStackedTruth->Draw();
//   leg->Draw();
//   
//   
//   
  
    
  TCanvas * cOverview = new TCanvas ("cOverview", "cOverview", 700, 1400);
  cOverview->Divide(1,4);
  
  cOverview->cd(1);
  hGrid_DRT_S->Draw("LEGO2Z");
  hGrid_DRT_S->SetStats(0);
  hGrid_DRT_S->SetTitle("Dual Readout HCAL Tower");
  hGrid_DRT_S->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_DRT_S->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_DRT_S->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_DRT_S->GetYaxis()->SetRangeUser(minPhi, maxPhi);
  
  cOverview->cd(2);
  hGrid_EC_T->Draw("LEGO2Z");
  hGrid_EC_T->SetStats(0);
  hGrid_EC_T->SetTitle("E1+E2 [GeV]");
  hGrid_EC_T->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_T->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_T->GetZaxis()->SetTitle("Energy [GeV]");
  hGrid_EC_T->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_EC_T->GetYaxis()->SetRangeUser(minPhi, maxPhi);
  
  cOverview->cd(3);
  hGrid_T2->Draw("LEGO2Z");
  hGrid_T2->SetStats(0);
  hGrid_T2->SetTitle("T2 [MeV]");
  hGrid_T2->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_T2->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_T2->GetZaxis()->SetTitle("Energy [MeV]");
  hGrid_T2->GetXaxis()->SetRangeUser(minTheta, maxTheta);
  hGrid_T2->GetYaxis()->SetRangeUser(minPhi, maxPhi);
  
  cOverview->cd(4);
  hStackedTruth->Draw();
  leg->Draw();
  
//   gPad->SetLogz();
//   cOverview->SaveAs(Form("plots/event_display/wwln/cOverview_%d.png", selEv));
  cOverview->SaveAs(Form("plots/event_display/zjj_100/cOverview_%d.png", selEv));
  
  
  
  TCanvas * cPolar = new TCanvas ("cPolar", "cPolar", 5000, 5000);
  cPolar->cd();
  gPad->DrawFrame(-maxRadius_DRT*1.2, -maxRadius_DRT*1.2, maxRadius_DRT*1.2, maxRadius_DRT*1.2);
//   gPad->SetLogz();
  for (long unsigned int iline = 0; iline< charged_tracks.size(); iline++)
  {
      charged_tracks[iline]->SetLineColor(kRed);
      charged_tracks[iline]->SetLineWidth(3);
      charged_tracks[iline]->Draw("same");
  }
  for (long unsigned int iline = 0; iline< muon_tracks.size(); iline++)
  {
      muon_tracks[iline]->SetLineColor(kBlue);
      muon_tracks[iline]->SetLineWidth(3);
      muon_tracks[iline]->Draw("same");
  }
  for (long unsigned int iline = 0; iline< photon_tracks.size(); iline++)
  {
      photon_tracks[iline]->SetLineColor(kGreen+1);
      photon_tracks[iline]->SetLineStyle(7);
      photon_tracks[iline]->SetLineWidth(2);
      photon_tracks[iline]->Draw("same");
  }
  for (long unsigned int iline = 0; iline< neutrhad_tracks.size(); iline++)
  {
      neutrhad_tracks[iline]->SetLineColor(kBlack);
      neutrhad_tracks[iline]->SetLineStyle(7);
      neutrhad_tracks[iline]->SetLineWidth(2);
      neutrhad_tracks[iline]->Draw("same");
  }
  
  
  gStyle->SetPalette(kRainBow);
  hPolar_EC_F->SetTitle(Form("event: %d", selEv));
  hPolar_EC_F->Draw("same colz pol");
  hPolar_EC_R->Draw("same col pol");
  
  hPolar_T1->Draw("same col pol");
  hPolar_T2->Draw("same col pol");
  
//   gStyle->SetPalette(kSolar);
    gStyle->SetPalette(kBlueRedYellow);
//     gStyle->SetPalette(kViridis);
  hPolar_DRT_S->Draw("same col pol");
  gPad->SetLogz();
  gPad->SetGrid();
  
  
  cPolar->SaveAs(Form("plots/event_display/zjj_100/phi_plots/cPolar_%d.png", selEv));
  cPolar->SaveAs(Form("plots/event_display/zjj_100/phi_plots/cPolar_%d.pdf", selEv));
  
  
//   theApp->Run();
}







