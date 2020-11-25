// g++ -Wall -o plotThetaPhiGrids plotThetaPhiGrids.C  VectorSmallestInterval.cc myG4Tree.cc myG4Tree.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs`


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
//   int energy = 100;
  
  std::string particle_name = "K^{0L}, 0T";
//   std::string output_tag = "kaon_central";
  std::string output_tag = "muon_iso";
    
  //define histos
  
  int NPHI_TL1   = 186-1;
  int NTHETA_TL1 = 29*16-1;
  int NPHI_TL2   = 186*16-1;
  int NTHETA_TL2 = 29-1;
  
  int NPHI_EC    = 1130;
  int NTHETA_EC  = 180+180-1;
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
  TH2F * hGrid_TT = new TH2F ("hGrid_TT", "hGrid_TT", NTHETA_TL1, minTheta-bin_width_theta_TL1/2, maxTheta-bin_width_theta_TL1/2, NPHI_TL2, minPhi-bin_width_phi_TL2/2, maxPhi-bin_width_phi_TL2/2);
  
  TH2F * hGrid_EC_F = new TH2F ("hGrid_EC_F", "hGrid_EC_F", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hGrid_EC_R = new TH2F ("hGrid_EC_R", "hGrid_EC_R", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  TH2F * hGrid_EC_T = new TH2F ("hGrid_EC_T", "hGrid_EC_T", NTHETA_EC, minTheta-bin_width_theta_EC/2, maxTheta-bin_width_theta_EC/2, NPHI_EC, minPhi-bin_width_phi_EC/2, maxPhi-bin_width_phi_EC/2);
  
  TH2F * hGrid_DRT_S = new TH2F ("hGrid_DRT_S", "hGrid_DRT_S", NTHETA_DRT, minTheta-bin_width_theta_DRT/2, maxTheta-bin_width_theta_DRT/2, NPHI_DRT, minPhi-bin_width_phi_DRT/2, maxPhi-bin_width_phi_DRT/2);
  TH2F * hGrid_DRT_C = new TH2F ("hGrid_DRT_C", "hGrid_DRT_C", NTHETA_DRT, minTheta-bin_width_theta_DRT/2, maxTheta-bin_width_theta_DRT/2, NPHI_DRT, minPhi-bin_width_phi_DRT/2, maxPhi-bin_width_phi_DRT/2);
  
  
  
  SCEPCal_GeometryHelper myGeometry;
  
  //run over energy scan
  TFile * RunFile = new TFile("../root_files/iso_gun/iso_gun_mu_100GeV_T+E.root","READ"); 
//   TFile * RunFile = new TFile("../root_files/iso_gun/central_kaon0L_10GeV_ALL.root","READ"); 
      
  TTree* TreeRun = (TTree*) RunFile->Get("B4");
  myG4TreeVars myTV;
  InitG4Tree (TreeRun, myTV);
      
      
      

  ///*******************************************///
  ///		 Run over events	    ///
  ///*******************************************///
    
  int NEVENTS = TreeRun->GetEntries();
  std::cout << "NEVENTS = " << NEVENTS << std::endl;

//   NEVENTS = 5000;
     
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
                                        
     
      TreeRun->GetEntry(iEvt);
      std::cout << "processing event: " << iEvt << "\r" << std::flush;
      
      double px  = myTV.PrimaryParticleMomentum->at(0);
      double py  = myTV.PrimaryParticleMomentum->at(1);
      double pz  = myTV.PrimaryParticleMomentum->at(2);
      double P   = sqrt(px*px+py*py+pz*pz);
      px/= P;
      py/= P;
      pz/= P;
      
      double phi   = atan(py/px);
      if (px<0. && py <0.)   {phi = phi - M_PI;}
      if (px<0. && py >0.)   {phi = M_PI + phi;}
      double eta   = -atanh(pz);      
      double theta = 2*atan(exp(-eta));
          
//       std::cout << "phi = " << phi << " :: theta = " << theta << std::endl;
          
      //**************************************************************//
      //                           DR HCAL
      //**************************************************************//

      
      
      for (unsigned int i = 0; i<myTV.VectorL->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'l');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorL->at(i)/1000.;      
          double this_cher  = myTV.VectorSignalsCherL->at(i)/1000.;
          hGrid_DRT_S ->Fill(this_theta, this_phi, this_ene);                              
          hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher);                              
      }
      for (unsigned int i = 0; i<myTV.VectorR->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'r');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorR->at(i)/1000.;                    
          double this_cher  = myTV.VectorSignalsCherR->at(i)/1000.;      
          hGrid_DRT_S ->Fill(this_theta, this_phi, this_ene);                              
          hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher);                              
      }
      
      
      
      //**************************************************************//
      //                             ECAL
      //**************************************************************//

      
      
      for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
      {                            
              
          TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;                    
          
//           if (fabs(myTV.VecHit_CrystalID->at(i)) <1000000)
          if (true)
          {
            hGrid_EC_F ->Fill(this_theta, this_phi, myTV.VecHit_ScepEneDepF->at(i)/1000);          
            hGrid_EC_R ->Fill(this_theta, this_phi, myTV.VecHit_ScepEneDepR->at(i)/1000);          
            hGrid_EC_T ->Fill(this_theta, this_phi, this_ene);                    
          }

      }

      
      
      //**************************************************************//
      //                            Timing
      //**************************************************************//
      
//       double c_speed = 1./299792458*1e9; //mm per picosecond          
//       double time_acceptance = 3; //time window to reject out of time hits
            
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_F->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_F->at(i),1);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepF->at(i);          
          hGrid_T1 ->Fill(this_theta, this_phi,this_ene);                                             
      }
      
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_R->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_R->at(i),2);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepR->at(i);          
          hGrid_T2 ->Fill(this_theta, this_phi,this_ene);                                             
      }
      
       
      std::vector<double>::iterator max_ene;
      if (myTV.VecHit_Timing_CrystalID_F->size()>0 && myTV.VecHit_Timing_CrystalID_R->size()>0) // && myTV.VecHit_Timing_CrystalID_F->size()<maxTimingHits)
      {
          max_ene         = std::max_element(myTV.VecHit_Timing_ScepEneDepF->begin(), myTV.VecHit_Timing_ScepEneDepF->end());
          int T1_maxHit_index_temp = std::distance(myTV.VecHit_Timing_ScepEneDepF->begin(), max_ene);                             
          int T1_seed_crystal_ID_temp = myTV.VecHit_Timing_CrystalID_F->at(T1_maxHit_index_temp);
          
          max_ene         = std::max_element(myTV.VecHit_Timing_ScepEneDepR->begin(), myTV.VecHit_Timing_ScepEneDepR->end());
          int T2_maxHit_index_temp = std::distance(myTV.VecHit_Timing_ScepEneDepR->begin(), max_ene);               
          int T2_seed_crystal_ID_temp = myTV.VecHit_Timing_CrystalID_R->at(T2_maxHit_index_temp);
                        
//           TT_maxHit_time     = myTV.VecHit_Timing_ScepTimeF->at(T1_maxHit_index_temp)/myTV.VecHit_Timing_ScepEneDepF->at(T1_maxHit_index_temp);              
          
          //               TVector3 seed_vec =  myGeometry.GetCrystalTimingBothVec(T1_seed_crystal_ID_temp, T2_seed_crystal_ID_temp);
          TVector3 seed_vec =  myGeometry.GetCrystalTimingBothVec(T1_seed_crystal_ID_temp, T2_seed_crystal_ID_temp);              
          double this_phi = seed_vec.Phi();
          double this_theta = seed_vec.Theta();
//           TT_maxHit_distance =  sqrt(pow(seed_vec.X(),2) + pow(seed_vec.Y(),2) + pow(seed_vec.Z(),2));
          
          hGrid_TT ->Fill(this_theta, this_phi, myTV.VecHit_Timing_ScepEneDepF->at(T1_maxHit_index_temp)+myTV.VecHit_Timing_ScepEneDepR->at(T2_maxHit_index_temp));                                             
                    
       }                  
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
  
  TLine * lEndcapMinus = new TLine(M_PI/4*1, -M_PI, M_PI/4*1, M_PI);
  TLine * lEndcapPlus  = new TLine(M_PI/4*3, -M_PI, M_PI/4*3, M_PI);
  lEndcapMinus->SetLineColor(kRed);
  lEndcapMinus->SetLineWidth(2);
  lEndcapPlus->SetLineColor(kRed);
  lEndcapPlus->SetLineWidth(2);
  
  float phiRange = M_PI/4;
  float thetaRange = M_PI/8;
  
  
  TCanvas * cGrid_DRT_S = new TCanvas("cGrid_DRT_S", "cGrid_DRT_S", 600, 900);
  cGrid_DRT_S->cd();
  hGrid_DRT_S->Draw("COLZ");
  hGrid_DRT_S->SetStats(0);
  hGrid_DRT_S->SetTitle("Dual Readout HCAL Tower");
  hGrid_DRT_S->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_DRT_S->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_DRT_S->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_DRT_S->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_DRT_S->SaveAs(Form("plots/cGrid_DRT_S_%s.png", output_tag.c_str()));
  

  TCanvas * cGrid_DRT_C = new TCanvas("cGrid_DRT_C", "cGrid_DRT_C", 600, 900);
  cGrid_DRT_C->cd();
  hGrid_DRT_C->Draw("COLZ");
  hGrid_DRT_C->SetStats(0);
  hGrid_DRT_C->SetTitle("Dual Readout HCAL Tower");
  hGrid_DRT_C->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_DRT_C->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_DRT_C->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_DRT_C->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_DRT_C->SaveAs(Form("plots/cGrid_DRT_C_%s.png", output_tag.c_str()));
  
  
  TCanvas * cGrid_EC_T = new TCanvas("cGrid_EC_T", "cGrid_EC_T", 600, 900);
  cGrid_EC_T->cd();
  hGrid_EC_T->Draw("COLZ");
  hGrid_EC_T->SetStats(0);
  hGrid_EC_T->SetTitle("E1+E2");
  hGrid_EC_T->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_T->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_T->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_EC_T->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_EC_T->SaveAs(Form("plots/cGrid_EC_T_%s.png", output_tag.c_str()));
  
  
     TCanvas * cGrid_EC_F = new TCanvas("cGrid_EC_F", "cGrid_EC_F", 600, 900);
  cGrid_EC_F->cd();
  hGrid_EC_F->Draw("COLZ");
  hGrid_EC_F->SetStats(0);
  hGrid_EC_F->SetTitle("E1");
  hGrid_EC_F->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_F->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_F->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_EC_F->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_EC_F->SaveAs(Form("plots/cGrid_EC_F_%s.png", output_tag.c_str()));
  
  
    TCanvas * cGrid_EC_R = new TCanvas("cGrid_EC_R", "cGrid_EC_R", 600, 900);
  cGrid_EC_R->cd();
  hGrid_EC_R->Draw("COLZ");
  hGrid_EC_R->SetStats(0);
  hGrid_EC_R->SetTitle("E2");
  hGrid_EC_R->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_R->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_R->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_EC_R->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_EC_R->SaveAs(Form("plots/cGrid_EC_R_%s.png", output_tag.c_str()));
  
  
    TCanvas * cGrid_T1 = new TCanvas("cGrid_T1", "cGrid_T1", 600, 900);
  cGrid_T1->cd();
  hGrid_T1->Draw("COLZ");
  hGrid_T1->SetStats(0);
  hGrid_T1->SetTitle("T1");
  hGrid_T1->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_T1->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_T1->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_T1->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_T1->SaveAs(Form("plots/cGrid_T1_%s.png", output_tag.c_str()));
  
  
    TCanvas * cGrid_T2 = new TCanvas("cGrid_T2", "cGrid_T2", 600, 900);
  cGrid_T2->cd();
  hGrid_T2->Draw("COLZ");
  hGrid_T2->SetStats(0);
  hGrid_T2->SetTitle("T2");
  hGrid_T2->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_T2->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_T2->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_T2->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_T2->SaveAs(Form("plots/cGrid_T2_%s.png", output_tag.c_str()));
  
  
    TCanvas * cGrid_TT = new TCanvas("cGrid_TT", "cGrid_TT", 600, 900);
  cGrid_TT->cd();
  hGrid_TT->Draw("COLZ");
  hGrid_TT->SetStats(0);
  hGrid_TT->SetTitle("MaxHit (T1+T2)");
  hGrid_TT->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_TT->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_TT->GetXaxis()->SetRangeUser(-thetaRange+M_PI/2, thetaRange+M_PI/2);
  hGrid_TT->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_TT->SaveAs(Form("plots/cGrid_TT_%s.png", output_tag.c_str()));

   
   theApp->Run();
}


