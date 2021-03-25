// g++ -Wall -o makeImagesForCNN makeImagesForCNN.C myG4Tree.cc myG4Tree.hh CNN_Tree.cc CNN_Tree.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs`


#include "SCEPCal_GeometryHelper.hh"
#include "myG4Tree.hh"
#include "CNN_Tree.hh"

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
    
  
  
  bool SAVEPLOTS = true;  
  bool WRITEOUTPUT = false;  
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
  
  const int imageSize = 45;
//   float image_E1[imageSize][imageSize];
  
  
  double minPhi   = -M_PI;
  double maxPhi   =  M_PI;  
  double minTheta = -M_PI/2;
  double maxTheta =  M_PI/2;
  
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
  
  //images
  TH2F * hImage_EC_F  = new TH2F ("hImage_EC_F", "hImage_EC_F",     imageSize, 0, imageSize, imageSize, 0, imageSize);
  TH2F * hImage_EC_R  = new TH2F ("hImage_EC_R", "hImage_EC_R",     imageSize, 0, imageSize, imageSize, 0, imageSize);
  TH2F * hImage_EC_T  = new TH2F ("hImage_EC_T", "hImage_EC_T",     imageSize, 0, imageSize, imageSize, 0, imageSize);  
  TH2F * hImage_DRT_S = new TH2F ("hImage_DRT_S", "hImage_DRT_S",   imageSize, 0, imageSize, imageSize, 0, imageSize);
  TH2F * hImage_DRT_C = new TH2F ("hImage_DRT_C", "hImage_DRT_C",   imageSize, 0, imageSize, imageSize, 0, imageSize);  
  TH2F * hImage_TT    = new TH2F ("hImage_TT", "hImage_TT",         imageSize, 0, imageSize, imageSize, 0, imageSize);
      
  
  SCEPCal_GeometryHelper myGeometry;
  
  std::vector<std::string> filenames;
  
//   filenames.push_back("../root_files/prod/output_SCEPCal_fixedPos_e-_Iso+Uniform1-100_GeV.root");
//   filenames.push_back("../root_files/prod/output_SCEPCal_fixedPos_pi-_Iso+Uniform1-100_GeV.root");  
  
//   filenames.push_back("../root_files/prod/output_SCEPCal_fixedPos_gamma_Iso+Uniform1-100_GeV.root");
//   filenames.push_back("../root_files/prod/output_SCEPCal_fixedPos_pi0_Iso+Uniform1-100_GeV.root");
  
  filenames.push_back("../root_files/prod/output_SCEPCal_fixedPos_gamma_Iso+Uniform1-100_GeV.root");
  filenames.push_back("../root_files/prod/output_SCEPCal_fixedPos_kaon0L_Iso+Uniform1-100_GeV.root");
  
  
  TChain * TreeRun = new TChain("B4", "B4");
  for (long unsigned int iFile = 0; iFile < filenames.size(); iFile++)
  {
    TreeRun->Add( filenames.at(iFile).c_str() );
    std::cout << "adding file = " << filenames.at(iFile) << std::endl;
  }

  
  
  //read input tree from Geant4f root file
//   TTree* TreeRun = (TTree*) RunFile->Get("B4");
  myG4TreeVars myTV;
  InitG4Tree (TreeRun, myTV);
  TreeRun->GetEntry(0);
  std::cout << "primary name: " << myTV.PrimaryParticleName << std::endl;  

  //create output tree for CNN
//   TFile* outputFile = new TFile("../CNN_trees/output_forCNN_pi-e-_Iso+Uniform.root","RECREATE"); 
//   TFile* outputFile = new TFile("../CNN_trees/output_forCNN_gamma_pi0_Iso+Uniform.root","RECREATE");
   TFile* outputFile = new TFile("../CNN_trees/output_forCNN_gamma_kaon0L_Iso+Uniform.root","RECREATE");
  
//   TFile* outputFile = new TFile("../CNN_trees/output_temp.root","RECREATE"); 
  
  
  TTree* outputTree = new TTree("CNN", "Tree for CNN");
  myCNNTreeVars cnnTV;
  InitCNNTree (outputTree, cnnTV);

  

  ///*******************************************///
  ///		 Run over events	    ///
  ///*******************************************///
    
  int NEVENTS = TreeRun->GetEntries();
  std::cout << "NEVENTS = " << NEVENTS << std::endl;

  int MAX_EVENTS = 200000;
  if (MAX_EVENTS<NEVENTS) NEVENTS = MAX_EVENTS;
  
  float minEneECAL = 10; // MeV
  float minEneHCAL = 10; // MeV
  float minEneTiming = 0.5; // MeV
   
     
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
                                        
//     std::cout <<"iEvt = " << iEvt << std::endl;
      
      TreeRun->GetEntry(iEvt);
      if (iEvt%100 == 0) std::cout << "processing event: " << iEvt << "\r" << std::flush;
      
      double px  = myTV.PrimaryParticleMomentum->at(0);
      double py  = myTV.PrimaryParticleMomentum->at(1);
      double pz  = myTV.PrimaryParticleMomentum->at(2);
      double P   = sqrt(px*px+py*py+pz*pz);
      px/= P;
      py/= P;
      pz/= P;
      
      
      

//       if (P/1000<8 || P/1000>10) continue;
//       std::cout << "P = " << P/1000 << " GeV" << std::endl;
      
      
      double phi   = atan(py/px);
      if (px<0. && py <0.)   {phi = phi - M_PI;}
      if (px<0. && py >0.)   {phi = M_PI + phi;}
//       double eta   = -atanh(pz);      
//       double theta = 2*atan(exp(-eta));
          
//       std::cout << "phi = " << phi << " :: theta = " << theta << std::endl;

      
      //**************************************************************//
      //                             ECAL
      //**************************************************************//
//       std::cout << "ECAL" << std::endl;
      //use ECAL seed to align pictures of also HCAL and Timing layers
      std::vector<double>::iterator max_ene;

      // find hit with max energy (seed)
      int maxHit_index = -1;
      int seed_crystal_ID = 0;
      double phi_seed = 0.;
      double theta_seed = 0.;
          
      if (myTV.VecHit_CrystalID->size()>0)
      {
          max_ene = std::max_element(myTV.VecHit_ScepEneDepR->begin(), myTV.VecHit_ScepEneDepR->end());
          maxHit_index = std::distance(myTV.VecHit_ScepEneDepR->begin(), max_ene);
          seed_crystal_ID = myTV.VecHit_CrystalID->at(maxHit_index);
          TVector3 seed_vec =  myGeometry.GetCrystalVec(seed_crystal_ID);
          phi_seed = seed_vec.Phi();
          theta_seed = seed_vec.Theta();
      }
      
      for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
      {                            
              
          TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;                    
          
//           if (fabs(myTV.VecHit_CrystalID->at(i)) <1000000)
          if (true)
          {
            if (myTV.VecHit_ScepEneDepF->at(i) >minEneECAL) hGrid_EC_F ->Fill(this_theta-theta_seed, this_phi-phi_seed, myTV.VecHit_ScepEneDepF->at(i)/1000);
            if (myTV.VecHit_ScepEneDepR->at(i) >minEneECAL) hGrid_EC_R ->Fill(this_theta-theta_seed, this_phi-phi_seed, myTV.VecHit_ScepEneDepR->at(i)/1000);
            if (myTV.VecHit_ScepEneDepR->at(i) >minEneECAL) hGrid_EC_T ->Fill(this_theta-theta_seed, this_phi-phi_seed, this_ene);                        
          }
      }

      
      
                
      //**************************************************************//
      //                           DR HCAL
      //**************************************************************//
//       std::cout << "DR HCAL" << std::endl;
      
      
      for (unsigned int i = 0; i<myTV.VectorL->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'l');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
//           double this_ene   = myTV.VectorL->at(i)/1000.;      
          double this_scint = myTV.VectorSignalsL->at(i)/1000.;      
          double this_cher  = myTV.VectorSignalsCherL->at(i)/1000.;
          hGrid_DRT_S ->Fill(this_theta-theta_seed, this_phi-phi_seed, this_scint);                              
          hGrid_DRT_C ->Fill(this_theta-theta_seed, this_phi-phi_seed, this_cher);                              
      }
      for (unsigned int i = 0; i<myTV.VectorR->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'r');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
//           double this_ene   = myTV.VectorR->at(i)/1000.;                    
          double this_scint = myTV.VectorSignalsR->at(i)/1000.;                    
          double this_cher  = myTV.VectorSignalsCherR->at(i)/1000.;      
          hGrid_DRT_S ->Fill(this_theta-theta_seed, this_phi-phi_seed, this_scint);                              
          hGrid_DRT_C ->Fill(this_theta-theta_seed, this_phi-phi_seed, this_cher);                              
      }
      
      
      
      
      //**************************************************************//
      //                            Timing
      //**************************************************************//
//       std::cout << "Timing" << std::endl;
//       double c_speed = 1./299792458*1e9; //mm per picosecond          
//       double time_acceptance = 3; //time window to reject out of time hits
            
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_F->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_F->at(i),1);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepF->at(i);          
          if (this_ene >minEneTiming) hGrid_T1 ->Fill(this_theta-theta_seed, this_phi-phi_seed,this_ene);                                             
      }
      
      for (long unsigned int i = 0; i<myTV.VecHit_Timing_CrystalID_R->size(); i++)
      {                            
          TVector3 this_vec =  myGeometry.GetCrystalTimingVec(myTV.VecHit_Timing_CrystalID_R->at(i),2);
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();                     
          double this_ene = myTV.VecHit_Timing_ScepEneDepR->at(i);          
          if (this_ene >minEneTiming)hGrid_T2 ->Fill(this_theta-theta_seed, this_phi-phi_seed,this_ene);                                             
      }
      
       
//       std::vector<double>::iterator max_ene;
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
          double this_ene = myTV.VecHit_Timing_ScepEneDepF->at(T1_maxHit_index_temp)+myTV.VecHit_Timing_ScepEneDepR->at(T2_maxHit_index_temp);
          
          if (this_ene >minEneTiming) hGrid_TT ->Fill(this_theta-theta_seed, this_phi-phi_seed, this_ene);                                             
                    
       }    
       
       
//        std::cout << "primary name: " << myTV.PrimaryParticleName << std::endl;       
       
       
       //fill CNN tree
//        std::cout << "Filling CNN tree" << std::endl;
       
       cnnTV.PrimaryParticleEnergy      =  myTV.PrimaryParticleEnergy;
       cnnTV.CNNPrimaryParticleName        =  myTV.PrimaryParticleName;
//        std::cout << "primary energy: " << myTV.PrimaryParticleEnergy << std::endl;
       cnnTV.PrimaryParticleMomentum[0] = px;
       cnnTV.PrimaryParticleMomentum[1] = py;
       cnnTV.PrimaryParticleMomentum[2] = pz;              
       
//        cnnTV.PrimaryParticleName        = "test";
       
//        std::cout  << " particle name = " << myTV.PrimaryParticleName << " :: " << cnnTV.PrimaryParticleName << std::endl;
       
       cnnTV.theta_seed                 = theta_seed;
       cnnTV.phi_seed                   = phi_seed;
       
       //cut out images
       //ECAL
//        std::cout << "Filling images ECAL" << std::endl;
       int iX_center_ECAL = NTHETA_EC/2.+1;
       int iY_center_ECAL = NPHI_EC/2.+1;
       for (int iPixelX = 0; iPixelX<imageSize; iPixelX++)
       {
           int iBinX = iX_center_ECAL-imageSize/2 + iPixelX + 1;
           
           for (int iPixelY = 0; iPixelY<imageSize; iPixelY++)
           {           
               int iBinY = iY_center_ECAL-imageSize/2 + iPixelY + 1;
               hImage_EC_F->Fill(iPixelX, iPixelY, hGrid_EC_F->GetBinContent(iBinX, iBinY));
               hImage_EC_R->Fill(iPixelX, iPixelY, hGrid_EC_R->GetBinContent(iBinX, iBinY));
               hImage_EC_T->Fill(iPixelX, iPixelY, hGrid_EC_T->GetBinContent(iBinX, iBinY));               
               cnnTV.image_E1[iPixelX+iPixelY*imageSize] = hGrid_EC_F->GetBinContent(iBinX, iBinY);
               cnnTV.image_E2[iPixelX+iPixelY*imageSize] = hGrid_EC_R->GetBinContent(iBinX, iBinY);
//                cnnTV.image_ET[iPixelX+iPixelY*imageSize] = hGrid_EC_T->GetBinContent(iBinX, iBinY);
               
//                std::cout  << "image_E1 ["<< iPixelX<< "][" << iPixelY << "] : " << hGrid_EC_F->GetBinContent(iBinX, iBinY) << " ::: pixel_id = " << iPixelX+iPixelY*imageSize << std::endl;

           }
       }
       
       //HCAL
//        std::cout << "Filling images DR HCAL" << std::endl;
       int iX_center_HCAL = NTHETA_DRT/2.+1;
       int iY_center_HCAL = NPHI_DRT/2.+1;
       for (int iPixelX = 0; iPixelX<imageSize; iPixelX++)
       {
           int iBinX = iX_center_HCAL-imageSize/2 + iPixelX + 1;
           
           for (int iPixelY = 0; iPixelY<imageSize; iPixelY++)
           {           
               int iBinY = iY_center_HCAL-imageSize/2 + iPixelY + 1;
               hImage_DRT_S->Fill(iPixelX, iPixelY, hGrid_DRT_S->GetBinContent(iBinX, iBinY));
               hImage_DRT_C->Fill(iPixelX, iPixelY, hGrid_DRT_C->GetBinContent(iBinX, iBinY));               
//                cnnTV.image_DRT_S[iPixelX+iPixelY*imageSize] = hGrid_DRT_S->GetBinContent(iBinX, iBinY);
//                cnnTV.image_DRT_C[iPixelX+iPixelY*imageSize] = hGrid_DRT_C->GetBinContent(iBinX, iBinY);                              
           }
       }
       
       //Timing Grid
//        std::cout << "Filling images Timing" << std::endl;
       int iX_center_TT = NTHETA_TL1/2.+1;
       int iY_center_TT = NPHI_TL2/2.+1;
       for (int iPixelX = 0; iPixelX<imageSize; iPixelX++)
       {
           int iBinX = iX_center_TT-imageSize/2 + iPixelX + 1;
           
           for (int iPixelY = 0; iPixelY<imageSize; iPixelY++)
           {           
               int iBinY = iY_center_TT-imageSize/2 + iPixelY + 1;
               hImage_TT->Fill(iPixelX, iPixelY, hGrid_TT->GetBinContent(iBinX, iBinY));                              
               cnnTV.image_TT[iPixelX+iPixelY*imageSize] = hGrid_TT->GetBinContent(iBinX, iBinY);                              
           }
       }
       
       outputTree->Fill();
       
       
       //histo resets
       hGrid_TT->Reset();
       hGrid_EC_F->Reset();
       hGrid_EC_R->Reset();       
       hGrid_EC_T->Reset();       
       hGrid_DRT_S->Reset();
       hGrid_DRT_C->Reset();
       
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
  hGrid_DRT_S->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
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
  hGrid_DRT_C->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
  hGrid_DRT_C->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_DRT_C->SaveAs(Form("plots/cGrid_DRT_C_%s.png", output_tag.c_str()));
  
  TCanvas * cImage_DRT_S = new TCanvas("cImage_DRT_S", "cImage_DRT_S", 600, 600);
  cImage_DRT_S->cd();
  hImage_DRT_S->Draw("COLZ");
  hImage_DRT_S->SetStats(0);
  hImage_DRT_S->SetTitle("DRT_S");
  hImage_DRT_S->GetXaxis()->SetTitle("pixelX");
  hImage_DRT_S->GetYaxis()->SetTitle("pixelY");  
  gPad->SetLogz();
  if (SAVEPLOTS)   cImage_DRT_S->SaveAs(Form("plots/cImage_DRT_S_%s.png", output_tag.c_str()));
  
  TCanvas * cImage_DRT_C = new TCanvas("cImage_DRT_C", "cImage_DRT_C", 600, 600);
  cImage_DRT_C->cd();
  hImage_DRT_C->Draw("COLZ");
  hImage_DRT_C->SetStats(0);
  hImage_DRT_C->SetTitle("DRT_C");
  hImage_DRT_C->GetXaxis()->SetTitle("pixelX");
  hImage_DRT_C->GetYaxis()->SetTitle("pixelY");  
  gPad->SetLogz();
  if (SAVEPLOTS)   cImage_DRT_C->SaveAs(Form("plots/cImage_DRT_C_%s.png", output_tag.c_str()));
  
  
  
  
  TCanvas * cGrid_EC_T = new TCanvas("cGrid_EC_T", "cGrid_EC_T", 600, 900);
  cGrid_EC_T->cd();
  hGrid_EC_T->Draw("COLZ");
  hGrid_EC_T->SetStats(0);
  hGrid_EC_T->SetTitle("E1+E2");
  hGrid_EC_T->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_EC_T->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_EC_T->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
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
  hGrid_EC_F->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
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
  hGrid_EC_R->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
  hGrid_EC_R->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_EC_R->SaveAs(Form("plots/cGrid_EC_R_%s.png", output_tag.c_str()));
  
  
  TCanvas * cImage_EC_F = new TCanvas("cImage_EC_F", "cImage_EC_F", 600, 600);
  cImage_EC_F->cd();
  hImage_EC_F->Draw("COLZ");
  hImage_EC_F->SetStats(0);
  hImage_EC_F->SetTitle("E1");
  hImage_EC_F->GetXaxis()->SetTitle("pixelX");
  hImage_EC_F->GetYaxis()->SetTitle("pixelY");  
  gPad->SetLogz();
  if (SAVEPLOTS)   cImage_EC_F->SaveAs(Form("plots/cImage_EC_F_%s.png", output_tag.c_str()));
  
  TCanvas * cImage_EC_R = new TCanvas("cImage_EC_R", "cImage_EC_R", 600, 600);
  cImage_EC_R->cd();
  hImage_EC_R->Draw("COLZ");
  hImage_EC_R->SetStats(0);
  hImage_EC_R->SetTitle("E2");
  hImage_EC_R->GetXaxis()->SetTitle("pixelX");
  hImage_EC_R->GetYaxis()->SetTitle("pixelY");  
  gPad->SetLogz();
  if (SAVEPLOTS)   cImage_EC_R->SaveAs(Form("plots/cImage_EC_R_%s.png", output_tag.c_str()));
  
  
  
  
    TCanvas * cGrid_T1 = new TCanvas("cGrid_T1", "cGrid_T1", 600, 900);
  cGrid_T1->cd();
  hGrid_T1->Draw("COLZ");
  hGrid_T1->SetStats(0);
  hGrid_T1->SetTitle("T1");
  hGrid_T1->GetXaxis()->SetTitle("#theta [rad]");
  hGrid_T1->GetYaxis()->SetTitle("#phi [rad]");
  hGrid_T1->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
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
  hGrid_T2->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
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
  hGrid_TT->GetXaxis()->SetRangeUser(-thetaRange, thetaRange);
  hGrid_TT->GetYaxis()->SetRangeUser(-phiRange, phiRange);
  lEndcapMinus->Draw("same");
  lEndcapPlus->Draw("same");
  gPad->SetLogz();
  if (SAVEPLOTS)   cGrid_TT->SaveAs(Form("plots/cGrid_TT_%s.png", output_tag.c_str()));

  
  TCanvas * cImage_TT = new TCanvas("cImage_TT", "cImage_TT", 600, 600);
  cImage_TT->cd();
  hImage_TT->Draw("COLZ");
  hImage_TT->SetStats(0);
  hImage_TT->SetTitle("TT");
  hImage_TT->GetXaxis()->SetTitle("pixelX");
  hImage_TT->GetYaxis()->SetTitle("pixelY");  
  gPad->SetLogz();
  if (SAVEPLOTS)   cImage_TT->SaveAs(Form("plots/cImage_TT_%s.png", output_tag.c_str()));
  
  
  if (WRITEOUTPUT)
  {
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
    std::cout << "done writing output file!" << std::endl;
  }
  
  
  theApp->Run();
  
  
  
  
}


