// g++ -Wall -o truthClusterMatchEfficiency truthClusterMatchEfficiency.C  myG4Tree.cc myG4Tree.hh myTruthTree.cc myTruthTree.hh recoUtils.cc recoUtils.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs`


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


using namespace std;

bool FileExists(const char * filename)
{
  bool data = true;
  ifstream file(filename);
  if(file)
  {
      return data;
  }
  else
  {
    data = false;
    return data;
  }
}

bool RootFileExists(const char *filename)
{
  bool data = false;
  TFile *f = TFile::Open(filename);
  if ((!f) || (f->IsZombie()))
  {
      data = false;
      return data;
  }
  else
  {
    data = true;
    return data;
  }
  f->Close();
}

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

  int mycolors[10] = {kBlack, kGreen+1, kBlue, kRed, kYellow+1, kCyan+1, kViolet, kOrange+1, kGray+1, kPink};

  bool SAVEPLOTS = false;
  
  std::string output_tag = "zjj_scan_100";
  int NFILES = 4;
  if (argc>1) output_tag = argv[1];   
  if (argc>2) NFILES = atoi(argv[2]);   
  std::cout << "processing sample of: " << output_tag.c_str() << std::endl;  

  
  double drh_S_norm  = 407;
  float maxDeltaRMatch = 0.1;
  float maxDeltaRSeed  = 0.1;
  float etaAcceptance = 1.4;
  
  float ene_EC_th = 0.01;
  float EC_seed_th = 0.1;
  
  float ene_HC_th   = 0.01;    
  float HC_seed_th = 0.1;
  
  float MC_ene_th = 0.15;
      
  
  
  std::map <int, std::string> myPdgId;
  myPdgId [22] = "#gamma";
  myPdgId [130] = "K^{0,L}";
  myPdgId [2112] = "n";
  
  myPdgId [211] = "#pi^{#pm}";
  myPdgId [321] = "K^{#pm}";
  myPdgId [2212] = "p";
  
  myPdgId [11] = "e^{#pm}";
  myPdgId [13] = "#mu^{#pm}";
  
  
  int NBIN = 300;
  float minEff = 0; float maxEff = 3;
  
  std::map <int, TH1F*> hNTotGen;
  
  std::map <int, TH1F*> hEffGenMatchedToEcalCluster;
  std::map <int, TH1F*> hNEcalClustersMatchedToGen;
  std::map <int, TProfile*> pNEcalClustersMatchedToGen_vsEne;
  
  std::map <int, TH1F*> hEffGenMatchedToHcalCluster;
  std::map <int, TH1F*> hNHcalClustersMatchedToGen;
  std::map <int, TProfile*> pNHcalClustersMatchedToGen_vsEne;
  
  TH1F * hNEcalSeeds = new TH1F ("hNEcalSeeds", "hNEcalSeeds", 100, -0.5, 99.5);
  TH1F * hNHcalSeeds = new TH1F ("hNHcalSeeds", "hNHcalSeeds", 100, -0.5, 99.5);
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      
      hNTotGen[it->first] = new TH1F(Form("hNTotGen_%d", it->first),Form("hNTotGen_%d", it->first), 50, -0.5, 49.5);
      
      hEffGenMatchedToEcalCluster[it->first] = new TH1F(Form("hEffGenMatchedToEcalCluster_%d", it->first), Form("hEffGenMatchedToEcalCluster_%d", it->first), NBIN, minEff, maxEff);
      hNEcalClustersMatchedToGen[it->first]  = new TH1F(Form("hNEcalClustersMatchedToGen_%d", it->first),Form("hNEcalClustersMatchedToGen_%d", it->first), 10, -0.5, 9.5);
      pNEcalClustersMatchedToGen_vsEne[it->first] = new TProfile(Form("pNEcalClustersMatchedToGen_vsEne_%d", it->first),Form("pNEcalClustersMatchedToGen_vsEne_%d", it->first), 50, 0, 100);
      
      hEffGenMatchedToHcalCluster[it->first] = new TH1F(Form("hEffGenMatchedToHcalCluster_%d", it->first), Form("hEffGenMatchedToHcalCluster_%d", it->first), NBIN, minEff, maxEff);
      hNHcalClustersMatchedToGen[it->first]  = new TH1F(Form("hNHcalClustersMatchedToGen_%d", it->first),Form("hNHcalClustersMatchedToGen_%d", it->first), 10, -0.5, 9.5);
      pNHcalClustersMatchedToGen_vsEne[it->first] = new TProfile(Form("pNHcalClustersMatchedToGen_vsEne_%d", it->first),Form("pNHcalClustersMatchedToGen_vsEne_%d", it->first), 50, 0, 100);
      
  }
  
  
  TH1F * hNGenMatchedToCluster = new TH1F ("hNGenMatchedToCluster", "hNGenMatchedToCluster", 20, -0.5, 19.5);
  TH1F * hNGenMatchedToHcalCluster = new TH1F ("hNGenMatchedToHcalCluster", "hNGenMatchedToHcalCluster", 20, -0.5, 19.5);
  
  
  SCEPCal_GeometryHelper myGeometry;
  
  
  TChain * TreeRun = new TChain("B4", "B4");      
  TChain * TruthTree = new TChain("truth", "truth");  
  
  for (int iFile = 0; iFile<NFILES; iFile++)
  {
    std::string fname_reco;
    std::string fname_truth;

    if (output_tag == "hznb" || output_tag == "wwlj" || output_tag == "hzjnbn")
    {
//        fname_reco  = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/reco/output_SCEPCal_B0T_%s100k_job_%d.root", output_tag.c_str(), iFile);
//        fname_truth = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/mc_truth/B0T/%s100k_job_%d_output_tuple.root", output_tag.c_str(), iFile);
    }
    else 
    {
//        fname_reco  = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/reco/output_SCEPCal_B0T_%s_job_%d.root", output_tag.c_str(), iFile);
//        fname_truth = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/mc_truth/B0T/%s_job_%d_output_tuple.root", output_tag.c_str(), iFile);
        fname_reco  = Form("../root_files/hep_outputs/output_SCEPCal_B0T_%s_job_%d.root", output_tag.c_str(), iFile);
        fname_truth = Form("../../HepMC_Files/B0T/%s_job_%d_output_tuple.root", output_tag.c_str(), iFile);

    }

    if (RootFileExists(fname_reco.c_str())  && RootFileExists(fname_truth.c_str()) )
    {
      std::cout << "adding file: " << iFile << std::endl;
      TreeRun->Add(fname_reco.c_str());
      TruthTree->Add(fname_truth.c_str());    
    }
  }

  myG4TreeVars myTV;
  InitG4Tree (TreeRun, myTV);
  
  TreeRun->SetBranchStatus("*", 0);
  TreeRun->SetBranchStatus("VectorSignalsL", 1);
  TreeRun->SetBranchStatus("VectorSignalsR", 1);
  TreeRun->SetBranchStatus("VecHit_CrystalID", 1);
  TreeRun->SetBranchStatus("VecHit_ScepEneDepF", 1);
  TreeRun->SetBranchStatus("VecHit_ScepEneDepR", 1);
  
  
  myTruthTreeVars myTruthTV;
  InitTruthTree (TruthTree, myTruthTV);
  
  
  ///*******************************************///
  ///		 Run over events	        ///
  ///*******************************************///
  
  int maxEVENTS = 100000;
  int NEVENTS = TreeRun->GetEntries();
  
  std::cout << "NEVENTS = " << NEVENTS << std::endl;
  if (NEVENTS>maxEVENTS)  NEVENTS = maxEVENTS;
  std::cout << "... running on " << NEVENTS << " events" << std::endl;  
  
  TCanvas * cPlotImage = new TCanvas ("cPlotImage", "cPlotImage", 1500, 500);
  cPlotImage->Divide(3,1);
  
  TCanvas * cPlotImageSum = new TCanvas ("cPlotImageSum", "cPlotImageSum", 1500, 500);
  cPlotImageSum->Divide(3,1);
    
  int imageSize = 15;
  
  TH2F * hImage_E1  = new TH2F ("hImage_E1", "hImage_E1", imageSize, 0, imageSize, imageSize, 0, imageSize);
  hImage_E1->SetStats(0);
  TH2F * hImage_E1_Sum  = new TH2F ("hImage_E1_Sum", "hImage_E1_Sum", imageSize, 0, imageSize, imageSize, 0, imageSize);
  hImage_E1_Sum->SetStats(0);
  
  TH2F * hImage_E2  = new TH2F ("hImage_E2", "hImage_E2", imageSize, 0, imageSize, imageSize, 0, imageSize);
  hImage_E2->SetStats(0);
  TH2F * hImage_E2_Sum  = new TH2F ("hImage_E2_Sum", "hImage_E2_Sum", imageSize, 0, imageSize, imageSize, 0, imageSize);
  hImage_E2_Sum->SetStats(0);
  
  TH2F * hImage_HC  = new TH2F ("hImage_HC", "hImage_HC", imageSize, 0, imageSize, imageSize, 0, imageSize);
  hImage_HC->SetStats(0);
  TH2F * hImage_HC_Sum  = new TH2F ("hImage_HC_Sum", "hImage_HC_Sum", imageSize, 0, imageSize, imageSize, 0, imageSize);
  hImage_HC_Sum->SetStats(0);
  
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {

      
      TreeRun->GetEntry(iEvt);
      TruthTree->GetEntry(iEvt);
      
      
      std::cout << "processing event: " << iEvt << "\r" << std::flush;

      
                  
      //**************************************************************//
      //                           DR HCAL
      //**************************************************************//
      
      std::vector<CalHit> myHcHits;
      std::vector<CalSeed> myHcSeeds;
      
      for (unsigned int i = 0; i<myTV.VectorSignalsL->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'l');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_scint = myTV.VectorSignalsL->at(i);                
          if (this_scint/drh_S_norm>ene_HC_th)
          {
              CalHit new_hit;
              new_hit.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_hit.SetSide(-1);
              myHcHits.push_back(new_hit);
          }
          if (this_scint/drh_S_norm>HC_seed_th)
          {
              CalSeed new_seed;
              new_seed.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_seed.SetSide(-1);              
              myHcSeeds.push_back(new_seed);
          }
      }
      
      for (unsigned int i = 0; i<myTV.VectorSignalsR->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'r');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_scint = myTV.VectorSignalsR->at(i);
//           std::cout << "this_scint/drh_S_norm = " << this_scint/drh_S_norm << std::endl;
          if (this_scint/drh_S_norm>ene_HC_th)
          {
              CalHit new_hit;
              new_hit.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_hit.SetSide(1);
              myHcHits.push_back(new_hit);
//               std::cout << "passing ene hcal cut" << std::endl;
          }
          if (this_scint/drh_S_norm>HC_seed_th)
          {
              CalSeed new_seed;
              new_seed.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_seed.SetSide(1);
              myHcSeeds.push_back(new_seed);
          }
      }
      
//       std::cout << "Number of HCAL seeds found: " << myHcSeeds.size() << std::endl;
//       std::cout << "Cleaning up HCAL seeds too close to each other" << std::endl;

      std::vector<CalSeed>  myHcSeedsCleaned = CleanSeeds(myHcSeeds, maxDeltaRSeed);                        
      hNHcalSeeds->Fill(myHcSeedsCleaned.size());
      
      
//       std::cout << "Matching HCAL clusters with gen level" << std::endl;
      for (long unsigned int iseed = 0; iseed < myHcSeedsCleaned.size(); iseed++)
      { 
          CalSeed this_seed = myHcSeedsCleaned.at(iseed);
          float seed_theta = this_seed.GetTheta();
          float seed_phi   = this_seed.GetPhi();
          
          int nGenMatchedToCluster = 0;
          for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
          {
              
              int    pdgId = myTruthTV.mcs_pdgId->at(i);
              double ene   = myTruthTV.mcs_E->at(i);
//               if (ene<0.1) continue;
              if (pdgId == 12 || pdgId == 14 || pdgId == 16) continue; //ignore neutrinos
              
              double truth_phi   = myTruthTV.mcs_phi->at(i);
              double eta   = myTruthTV.mcs_eta->at(i);
              double truth_theta = 2*atan(exp(-eta));
              truth_theta = M_PI- truth_theta;
              
              
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatch)
              {
//                   std::cout  << "HCAL cluster (seedEne = "<< this_seed.GetEne() << " GeV) " << iseed << " matched to MC truth gen level particle " << pdgId << " (energy = " << ene << " GeV)" << std::endl;
                  this_seed.AddGenMatch(pdgId);
                  nGenMatchedToCluster++;
              }
            }
            
            hNGenMatchedToHcalCluster->Fill(nGenMatchedToCluster);   
      }
      
      
      //**************************************************************//
      //                             ECAL
      //**************************************************************//


      
      std::vector<CalHit> myEcHits;
      std::vector<CalHit> myEcHitsF;
      std::vector<CalHit> myEcHitsR;
      std::vector<CalSeed> myEcSeeds;
      
      for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
      {                            
              
          TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
          double this_phi = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;                    
          if (this_ene>ene_EC_th)
          {
              CalHit new_hit;
              new_hit.Init(myTV.VecHit_CrystalID->at(i), this_theta, this_phi, this_ene);
              myEcHits.push_back(new_hit);
              
              new_hit.Init(myTV.VecHit_CrystalID->at(i), this_theta, this_phi, myTV.VecHit_ScepEneDepF->at(i));
              myEcHitsF.push_back(new_hit);
              new_hit.Init(myTV.VecHit_CrystalID->at(i), this_theta, this_phi, myTV.VecHit_ScepEneDepR->at(i));
              myEcHitsR.push_back(new_hit);
          }

          // find hit with energy above seed threshold
          if (this_ene>EC_seed_th)
          {
              CalSeed new_seed;
              new_seed.Init(myTV.VecHit_CrystalID->at(i), this_theta, this_phi, this_ene);
              myEcSeeds.push_back(new_seed);
          }
      }
//       std::cout << "Number of ECAL seeds found: " << myEcSeeds.size() << std::endl;      
//       std::cout << "Cleaning up ECAL seeds too close to each other" << std::endl;
      std::vector<CalSeed>  myEcSeedsCleaned = CleanSeeds(myEcSeeds, maxDeltaRSeed);          
      hNEcalSeeds->Fill(myEcSeedsCleaned.size());
      
      
      

      
      //       std::cout << "Matching ECAL clusters with gen level" << std::endl;
      //Creating calo clusters      
      std::vector<CalCluster> myCalClusters;
      
      for (long unsigned int iseed = 0; iseed < myEcSeedsCleaned.size(); iseed++)
      { 
          CalSeed this_seed = myEcSeedsCleaned.at(iseed);
          float seed_theta = this_seed.GetTheta();
          float seed_phi   = this_seed.GetPhi();
          if (abs(this_seed.GetEta())>etaAcceptance) continue;
          int nGenMatchedToCluster = 0;
          
          for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
          {
              
              int    pdgId = myTruthTV.mcs_pdgId->at(i);
//               double ene   = myTruthTV.mcs_E->at(i);
//               if (ene<) continue;
              if (pdgId == 12 || pdgId == 14 || pdgId == 16) continue; //ignore neutrinos

              double mc_ene   = myTruthTV.mcs_E->at(i);
              double truth_phi   = myTruthTV.mcs_phi->at(i);
              double eta   = myTruthTV.mcs_eta->at(i);
              double truth_theta = 2*atan(exp(-eta));
              truth_theta = M_PI- truth_theta;
              
              
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatch)
              {
                  this_seed.AddGenMatch(pdgId);
                  
                  CalCluster thisCluster;
                  thisCluster.Init(this_seed, maxDeltaRSeed, 15);
                  thisCluster.Clusterize(myEcHits, myHcHits, myEcHitsF, myEcHitsR);
                  
                  float *image_E1;
                  float *image_E2;
                  float *image_HC;
                  image_E1 = thisCluster.GetImage("E1");
                  image_E2 = thisCluster.GetImage("E2");
                  image_HC = thisCluster.GetImage("HC");
                  
                  for (int iBinX = 0; iBinX<imageSize; iBinX++)
                  { 
                      for (int iBinY = 0; iBinY<imageSize; iBinY++)
                      {
                          int pixel = iBinX+iBinY*imageSize;
                          hImage_E1->SetBinContent(iBinX+1, iBinY+1, image_E1[pixel]);
                          hImage_E1_Sum->Fill(iBinX+0.5, iBinY+0.5, image_E1[pixel]);
                          
                          hImage_E2->SetBinContent(iBinX+1, iBinY+1, image_E2[pixel]);
                          hImage_E2_Sum->Fill(iBinX+0.5, iBinY+0.5, image_E2[pixel]);
                          
                          hImage_HC->SetBinContent(iBinX+1, iBinY+1, image_HC[pixel]);
                          hImage_HC_Sum->Fill(iBinX+0.5, iBinY+0.5, image_HC[pixel]);
                      }
                  }
                  
                  if (pdgId==22 && mc_ene>5)
//                   if (pdgId==11 && mc_ene>2)
                  {
                    cPlotImage->cd(1);
                    gPad->SetLogz();
                    hImage_E1->Draw("COLZ");
                    cPlotImage->cd(2);
                    gPad->SetLogz();
                    hImage_E2->Draw("COLZ");
                    cPlotImage->cd(3);
                    gPad->SetLogz();
                    hImage_HC->Draw("COLZ");
                    
                    cPlotImage->Update();
                    hImage_E1->Reset();
                    hImage_E2->Reset();
                    hImage_HC->Reset();
                    
                    cPlotImageSum->cd(1);
                    gPad->SetLogz();
                    hImage_E1_Sum->Draw("COLZ");
                    cPlotImageSum->cd(2);
                    gPad->SetLogz();
                    hImage_E2_Sum->Draw("COLZ");
                    cPlotImageSum->cd(3);
                    gPad->SetLogz();
                    hImage_HC_Sum->Draw("COLZ");
                    cPlotImageSum->Update();
                  }
                  
                  nGenMatchedToCluster++;
              }
            }
            
            hNGenMatchedToCluster->Fill(nGenMatchedToCluster);                        
      }
      
      for (long unsigned int iseed = 0; iseed < myHcSeedsCleaned.size(); iseed++)
      { 
          CalSeed this_seed = myHcSeedsCleaned.at(iseed);
          float seed_theta = this_seed.GetTheta();
          float seed_phi   = this_seed.GetPhi();            
          float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
      }
      
      
      //match of truth to some clusters (both ECAL and HCAL)
      
      std::map<int,int> nGenMatchToEcalCluster;
      std::map<int,int> nGenMatchToHcalCluster;
      std::map<int,int> nTotGen;      
      
      for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
      {
          nGenMatchToEcalCluster[it->first] = 0;
          nGenMatchToHcalCluster[it->first] = 0;
          nTotGen[it->first] = 0;
      }
      
      for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
      {
          int    pdgId = myTruthTV.mcs_pdgId->at(i);
          double ene   = myTruthTV.mcs_E->at(i);
          if (ene<MC_ene_th) continue;
          auto it = myPdgId.find(abs(pdgId));
          if (it== myPdgId.end())
          {
              if (abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16) std::cout << "skipping particle not in my pdg id list: " << pdgId << std::endl;
              continue;
          }
          
          double truth_phi   = myTruthTV.mcs_phi->at(i);
          double eta   = myTruthTV.mcs_eta->at(i);
          if (abs(eta)>etaAcceptance) continue;                    
          double truth_theta = 2*atan(exp(-eta));
          truth_theta = M_PI- truth_theta;
          
          int matchedToEcalCluster = 0;
          int matchedToHcalCluster = 0;
          nTotGen[abs(pdgId)]++;
          
          //ECAL match
          for (long unsigned int iseed = 0; iseed < myEcSeedsCleaned.size(); iseed++)
          { 
              CalSeed this_seed = myEcSeedsCleaned.at(iseed);
              float seed_theta = this_seed.GetTheta();
              float seed_phi   = this_seed.GetPhi();            
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatch)   matchedToEcalCluster ++;                            
          }
          hNEcalClustersMatchedToGen[abs(pdgId)]->Fill(matchedToEcalCluster);
          pNEcalClustersMatchedToGen_vsEne[abs(pdgId)]->Fill(ene, matchedToEcalCluster);
          if (matchedToEcalCluster>0) nGenMatchToEcalCluster[abs(pdgId)]++;
          
          //HCAL match
          for (long unsigned int iseed = 0; iseed < myHcSeedsCleaned.size(); iseed++)
          { 
              CalSeed this_seed = myHcSeedsCleaned.at(iseed);
              float seed_theta = this_seed.GetTheta();
              float seed_phi   = this_seed.GetPhi();            
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatch)   matchedToHcalCluster ++;                            
          }                    
          hNHcalClustersMatchedToGen[abs(pdgId)]->Fill(matchedToHcalCluster);
          pNHcalClustersMatchedToGen_vsEne[abs(pdgId)]->Fill(ene, matchedToHcalCluster);
          if (matchedToHcalCluster>0) nGenMatchToHcalCluster[abs(pdgId)]++;
          
      }
                  
      for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
      {
          hNTotGen[it->first]->Fill(nTotGen[it->first]);
          if (nTotGen[it->first]>0)
          {
              float eff = float(nGenMatchToEcalCluster[it->first])/float(nTotGen[it->first]);
              hEffGenMatchedToEcalCluster[it->first]->Fill(eff);
              
              eff = float(nGenMatchToHcalCluster[it->first])/float(nTotGen[it->first]);
              hEffGenMatchedToHcalCluster[it->first]->Fill(eff);
          }
      }
      
      
      
  }
  
        
  TCanvas * cNTotGen = new TCanvas ("cNTotGen", "cNTotGen", 600, 500);
  cNTotGen->cd();
  hNTotGen[22]->SetStats(0);
  hNTotGen[22]->SetTitle(0);
  hNTotGen[22]->Draw();
//   hNTotGen[22]->GetXaxis()->SetRangeUser(0, 1.4);
  hNTotGen[22]->GetYaxis()->SetRangeUser(1, hNTotGen[22]->GetMaximum()*5);
  hNTotGen[22]->GetXaxis()->SetTitle("Total number of gen-match particles");
//   hNTotGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  int color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hNTotGen[it->first]->Draw("same");
      hNTotGen[it->first]->SetLineWidth(2);
      hNTotGen[it->first]->SetLineColor(mycolors[color_it]);
      hNTotGen[it->first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(hNTotGen[it->first], it->second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNTotGen->SaveAs("plots_pfa/cNTotGen.png");
  
  
  
  //gen to ECAL cluster matches
  TCanvas * cNEcalClustersMatchedToGen = new TCanvas ("cNEcalClustersMatchedToGen", "cNEcalClustersMatchedToGen", 600, 500);
  cNEcalClustersMatchedToGen->cd();
  hNEcalClustersMatchedToGen[22]->SetStats(0);
  hNEcalClustersMatchedToGen[22]->SetTitle(0);
  hNEcalClustersMatchedToGen[22]->Draw();
  hNEcalClustersMatchedToGen[22]->GetXaxis()->SetRangeUser(-0.5, 5.5);
//   hNEcalClustersMatchedToGen[22]->GetYaxis()->SetRangeUser(1, hNEcalClustersMatchedToGen[22]->GetMaximum()*5);
  hNEcalClustersMatchedToGen[22]->GetXaxis()->SetTitle("N of ECAL clusters matched to gen particle");
//   hNEcalClustersMatchedToGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hNEcalClustersMatchedToGen[it->first]->Draw("same");
      hNEcalClustersMatchedToGen[it->first]->SetLineWidth(2);
      hNEcalClustersMatchedToGen[it->first]->SetLineColor(mycolors[color_it]);
      hNEcalClustersMatchedToGen[it->first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(hNEcalClustersMatchedToGen[it->first], it->second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNEcalClustersMatchedToGen->SaveAs("plots_pfa/cNEcalClustersMatchedToGen.png");

  
  TCanvas * cEffGenMatchedToEcalCluster = new TCanvas ("cEffGenMatchedToEcalCluster", "cEffGenMatchedToEcalCluster", 600, 500);
  cEffGenMatchedToEcalCluster->cd();
  hEffGenMatchedToEcalCluster[22]->SetStats(0);
  hEffGenMatchedToEcalCluster[22]->SetTitle(0);
  hEffGenMatchedToEcalCluster[22]->Draw();
  hEffGenMatchedToEcalCluster[22]->GetXaxis()->SetRangeUser(0, 1.4);
  hEffGenMatchedToEcalCluster[22]->GetYaxis()->SetRangeUser(1, hEffGenMatchedToEcalCluster[22]->GetMaximum()*5);
  hEffGenMatchedToEcalCluster[22]->GetXaxis()->SetTitle("Fraction of gen-match to at least one ECAL cluster");
//   hEffGenMatchedToEcalCluster[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hEffGenMatchedToEcalCluster[it->first]->Draw("same");
      hEffGenMatchedToEcalCluster[it->first]->SetLineWidth(2);
      hEffGenMatchedToEcalCluster[it->first]->SetLineColor(mycolors[color_it]);
      hEffGenMatchedToEcalCluster[it->first]->SetMarkerColor(mycolors[color_it]);    
      
      color_it++;
      leg->AddEntry(hEffGenMatchedToEcalCluster[it->first], it->second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cEffGenMatchedToEcalCluster->SaveAs("plots_pfa/cEffGenMatchedToEcalCluster.png");
  
  
    TCanvas * cNEcalClustersMatchedToGen_vsEne = new TCanvas ("cNEcalClustersMatchedToGen_vsEne", "cNEcalClustersMatchedToGen_vsEne", 600, 500);
  cNEcalClustersMatchedToGen_vsEne->cd();
  pNEcalClustersMatchedToGen_vsEne[22]->SetStats(0);
  pNEcalClustersMatchedToGen_vsEne[22]->SetTitle(0);
  pNEcalClustersMatchedToGen_vsEne[22]->Draw();
  pNEcalClustersMatchedToGen_vsEne[22]->GetXaxis()->SetRangeUser(0, 100);
  pNEcalClustersMatchedToGen_vsEne[22]->GetYaxis()->SetRangeUser(0, 2);
  pNEcalClustersMatchedToGen_vsEne[22]->GetXaxis()->SetTitle("MC truth particle energy [GeV]");
  pNEcalClustersMatchedToGen_vsEne[22]->GetYaxis()->SetTitle("<N> of HCAL clusters matched to gen particle");
//   pNEcalClustersMatchedToGen_vsEne[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
//   gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      pNEcalClustersMatchedToGen_vsEne[it->first]->Draw("same");
      pNEcalClustersMatchedToGen_vsEne[it->first]->SetLineWidth(2);
      pNEcalClustersMatchedToGen_vsEne[it->first]->SetLineColor(mycolors[color_it]);
      pNEcalClustersMatchedToGen_vsEne[it->first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(pNEcalClustersMatchedToGen_vsEne[it->first], it->second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNEcalClustersMatchedToGen_vsEne->SaveAs("plots_pfa/cNEcalClustersMatchedToGen_vsEne.png");
  

  //gen to HCAL cluster matches
  TCanvas * cNHcalClustersMatchedToGen = new TCanvas ("cNHcalClustersMatchedToGen", "cNHcalClustersMatchedToGen", 600, 500);
  cNHcalClustersMatchedToGen->cd();
  hNHcalClustersMatchedToGen[22]->SetStats(0);
  hNHcalClustersMatchedToGen[22]->SetTitle(0);
  hNHcalClustersMatchedToGen[22]->Draw();
  hNHcalClustersMatchedToGen[22]->GetXaxis()->SetRangeUser(-0.5, 5.5);
//   hNHcalClustersMatchedToGen[22]->GetYaxis()->SetRangeUser(1, hNHcalClustersMatchedToGen[22]->GetMaximum()*5);
  hNHcalClustersMatchedToGen[22]->GetXaxis()->SetTitle("N of HCAL clusters matched to gen particle");
//   hNHcalClustersMatchedToGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hNHcalClustersMatchedToGen[it->first]->Draw("same");
      hNHcalClustersMatchedToGen[it->first]->SetLineWidth(2);
      hNHcalClustersMatchedToGen[it->first]->SetLineColor(mycolors[color_it]);
      hNHcalClustersMatchedToGen[it->first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(hNHcalClustersMatchedToGen[it->first], it->second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNHcalClustersMatchedToGen->SaveAs("plots_pfa/cNHcalClustersMatchedToGen.png");

  
  TCanvas * cEffGenMatchedToHcalCluster = new TCanvas ("cEffGenMatchedToHcalCluster", "cEffGenMatchedToHcalCluster", 600, 500);
  cEffGenMatchedToHcalCluster->cd();
  hEffGenMatchedToHcalCluster[22]->SetStats(0);
  hEffGenMatchedToHcalCluster[22]->SetTitle(0);
  hEffGenMatchedToHcalCluster[22]->Draw();
  hEffGenMatchedToHcalCluster[22]->GetXaxis()->SetRangeUser(0, 1.4);
  hEffGenMatchedToHcalCluster[22]->GetYaxis()->SetRangeUser(1, hEffGenMatchedToHcalCluster[22]->GetMaximum()*5);
  hEffGenMatchedToHcalCluster[22]->GetXaxis()->SetTitle("Fraction of gen-match to at least one HCAL cluster");
//   hEffGenMatchedToHcalCluster[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hEffGenMatchedToHcalCluster[it->first]->Draw("same");
      hEffGenMatchedToHcalCluster[it->first]->SetLineWidth(2);
      hEffGenMatchedToHcalCluster[it->first]->SetLineColor(mycolors[color_it]);
      hEffGenMatchedToHcalCluster[it->first]->SetMarkerColor(mycolors[color_it]);    
      
      color_it++;
      leg->AddEntry(hEffGenMatchedToHcalCluster[it->first], it->second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cEffGenMatchedToHcalCluster->SaveAs("plots_pfa/cEffGenMatchedToHcalCluster.png");
  
  
  TCanvas * cNHcalClustersMatchedToGen_vsEne = new TCanvas ("cNHcalClustersMatchedToGen_vsEne", "cNHcalClustersMatchedToGen_vsEne", 600, 500);
  cNHcalClustersMatchedToGen_vsEne->cd();
  pNHcalClustersMatchedToGen_vsEne[22]->SetStats(0);
  pNHcalClustersMatchedToGen_vsEne[22]->SetTitle(0);
  pNHcalClustersMatchedToGen_vsEne[22]->Draw();
  pNHcalClustersMatchedToGen_vsEne[22]->GetXaxis()->SetRangeUser(0, 100);
  pNHcalClustersMatchedToGen_vsEne[22]->GetYaxis()->SetRangeUser(0, 2);
  pNHcalClustersMatchedToGen_vsEne[22]->GetXaxis()->SetTitle("MC truth particle energy [GeV]");
  pNHcalClustersMatchedToGen_vsEne[22]->GetYaxis()->SetTitle("<N> of HCAL clusters matched to gen particle");
//   pNHcalClustersMatchedToGen_vsEne[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
//   gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      pNHcalClustersMatchedToGen_vsEne[it->first]->Draw("same");
      pNHcalClustersMatchedToGen_vsEne[it->first]->SetLineWidth(2);
      pNHcalClustersMatchedToGen_vsEne[it->first]->SetLineColor(mycolors[color_it]);
      pNHcalClustersMatchedToGen_vsEne[it->first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(pNHcalClustersMatchedToGen_vsEne[it->first], it->second.c_str(), "lp");
  }
    
  leg->Draw();
  if (SAVEPLOTS) cNHcalClustersMatchedToGen_vsEne->SaveAs("plots_pfa/cNHcalClustersMatchedToGen_vsEne.png");
  
  
  
  TCanvas * cNHcalSeeds = new TCanvas ("cNHcalSeeds", "cNHcalSeeds", 600, 500);
  cNHcalSeeds->cd();
  hNHcalSeeds->Draw();
  hNHcalSeeds->GetXaxis()->SetTitle("N of HCAL clusters");  
  if (SAVEPLOTS) cNHcalSeeds->SaveAs("plots_pfa/cNHcalSeeds.png");
  
  TCanvas * cNGenMatchedToHcalCluster = new TCanvas ("cNGenMatchedToHcalCluster", "nGenMatchedToHcalCluster", 600, 500);
  cNGenMatchedToHcalCluster->cd();
  hNGenMatchedToHcalCluster->Draw();
  hNGenMatchedToHcalCluster->GetXaxis()->SetTitle("N of gen particles matched to HCAL cluster");
  if (SAVEPLOTS) cNGenMatchedToHcalCluster->SaveAs("plots_pfa/cNGenMatchedToHcalCluster.png");
  
  TCanvas * cNEcalSeeds = new TCanvas ("cNEcalSeeds", "cNEcalSeeds", 600, 500);
  cNEcalSeeds->cd();
  hNEcalSeeds->Draw();
  hNEcalSeeds->GetXaxis()->SetTitle("N of ECAL clusters");
  if (SAVEPLOTS) cNEcalSeeds->SaveAs("plots_pfa/cNEcalSeeds.png");
  
  TCanvas * cNGenMatchedToCluster = new TCanvas ("cNGenMatchedToCluster", "nGenMatchedToCluster", 600, 500);
  cNGenMatchedToCluster->cd();
  hNGenMatchedToCluster->Draw();
  hNGenMatchedToCluster->GetXaxis()->SetTitle("N of gen particles matched to ECAL cluster");
  if (SAVEPLOTS) cNGenMatchedToCluster->SaveAs("plots_pfa/cNGenMatchedToCluster.png");
  
  
  
  theApp->Run();
}







