// g++ -Wall -o truthClusterMatchEfficiency truthClusterMatchEfficiency.C  myG4Tree.cc myG4Tree.hh myTruthTree.cc myTruthTree.hh recoUtils.cc recoUtils.hh CNN_Tree.cc CNN_Tree.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs`


#include "SCEPCal_GeometryHelper.hh"
#include "myG4Tree.hh"
#include "myTruthTree.hh"
#include "recoUtils.hh"
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

  bool SAVEPLOTS   = true;
  bool WRITEOUTPUT = true;
  bool makePlots   = false;
  
  std::string output_tag = "zjj_scan_100";
  int NFILES = 4;
  if (argc>1) output_tag = argv[1];   
  if (argc>2) NFILES = atoi(argv[2]);   
  std::cout << "processing sample of: " << output_tag.c_str() << std::endl;  
  
  

  double ecal_S_norm = 0.985;
//   double ecal_S_norm = 1.;
  double LO  = 2000;
  double drh_S_norm  = 407;
  
  float maxDeltaRMatchEcal   = 0.01;
  float maxDeltaRSeedEcal    = 0.01;
  float maxDeltaRClusterEcal = 0.01;
  
  float maxDeltaRMatchHcal = 0.1;  
  float maxDeltaRSeedHcal  = 0.1;
  
  float etaAcceptance = 1.4;
  
  float ene_EC_th  = 0.01;
  float EC_seed_th = 0.15;
  
  float ene_HC_th  = 0.01;    
  float HC_seed_th = 0.15;
  
  float MC_ene_th = 0.4;
      
//   std::vector<float> ene_thresholds;
//   ene_thresholds.push_back(0.025);
//   ene_thresholds.push_back(0.050);
//   ene_thresholds.push_back(0.100);
//   ene_thresholds.push_back(0.200);
//   ene_thresholds.push_back(0.300);
//   ene_thresholds.push_back(0.400);
//   ene_thresholds.push_back(0.500);
//   
  
  std::map <int, std::string> myPdgId;
  myPdgId [22]   = "#gamma";
  myPdgId [130]  = "K^{0,L}";
  myPdgId [2112] = "n";
  
  myPdgId [211]  = "#pi^{#pm}";
  myPdgId [321]  = "K^{#pm}";
  myPdgId [2212] = "p";
  
  myPdgId [11] = "e^{#pm}";
  myPdgId [13] = "#mu^{#pm}";
  
  
  std::map <int, std::string> pdgToName;
  pdgToName [22]   = "gamma";
  pdgToName [130]  = "kaon0L";
  pdgToName [2112] = "n";
  
  pdgToName [211]  = "pi+-";
  pdgToName [321]  = "K+-";
  pdgToName [2212] = "p";
  
  pdgToName [11] = "e+-";
  pdgToName [13] = "mu+-";
  
  
  int NBIN = 300;
  float minEff = 0; float maxEff = 3;
  
  std::map <int, TH1F*> hNTotGen;
  std::map <int, TH1F*> hEneGen;
  std::map <int, TH1F*> hEneTotRes;
  std::map <int, TH1F*> hEneTotNarrowRes;
  std::map <int, TH1F*> hEneECALRes;
  
  std::map <int, TH1F*> hEffGenMatchedToEcalCluster;
  std::map <int, TH1F*> hNEcalClustersMatchedToGen;
  std::map <int, TProfile*> pNEcalClustersMatchedToGen_vsEne;
  
  std::map <int, TH1F*> hEffGenMatchedToHcalCluster;
  std::map <int, TH1F*> hNHcalClustersMatchedToGen;
  std::map <int, TProfile*> pNHcalClustersMatchedToGen_vsEne;
  
  std::map <int, TH1F*> hEffGenMatchedToCaloCluster;
  std::map <int, TH1F*> hNCaloClustersMatchedToGen;
  std::map <int, TProfile*> pNCaloClustersMatchedToGen_vsEne;

  std::map <int, TH2F*> h2EneReco_vsEneTruth;
  
  std::map <int, TH2F*> h2NClustersMatchedToGen;
  std::map <int, TH2F*> h2EneFracECAL_HCAL;
  std::map <int, TH2F*> h2EneNarrowFracECAL_HCAL;
  std::map <int, TH2F*> h2EneECALFracECAL_HCAL;
  std::map <int, TH2F*> h2EneECALFracECAL_Seed;
  
  TH1F * hNEcalSeeds = new TH1F ("hNEcalSeeds", "hNEcalSeeds", 100, -0.5, 99.5);
  TH1F * hNHcalSeeds = new TH1F ("hNHcalSeeds", "hNHcalSeeds", 100, -0.5, 99.5);
  TH1F * hNCaloSeeds = new TH1F ("hNCaloSeeds", "hNCaloSeeds", 100, -0.5, 99.5);
  
//   for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  for (auto it : myPdgId)
  {
      
      hNTotGen[it.first]         = new TH1F(Form("hNTotGen_%d", it.first),Form("hNTotGen_%d", it.first),                 50, -0.5, 49.5);
      hEneGen[it.first]          = new TH1F(Form("hEneGen_%d", it.first),Form("hEneGen_%d", it.first),                   100, 0., 50.);
      hEneTotRes[it.first]       = new TH1F(Form("hEneTotRes_%d", it.first),Form("hEneTotRes_%d", it.first),             400, -5., 5.);
      hEneTotNarrowRes[it.first] = new TH1F(Form("hEneTotNarrowRes_%d", it.first),Form("hEneTotNarrowRes_%d", it.first), 400, -5., 5.);
      hEneECALRes[it.first]      = new TH1F(Form("hEneECALRes_%d", it.first),Form("hEneECALRes_%d", it.first),           400, -5., 5.);
      
      hEffGenMatchedToEcalCluster[it.first] = new TH1F(Form("hEffGenMatchedToEcalCluster_%d", it.first), Form("hEffGenMatchedToEcalCluster_%d", it.first), NBIN, minEff, maxEff);
      hNEcalClustersMatchedToGen[it.first]  = new TH1F(Form("hNEcalClustersMatchedToGen_%d", it.first),Form("hNEcalClustersMatchedToGen_%d", it.first), 10, -0.5, 9.5);
      pNEcalClustersMatchedToGen_vsEne[it.first] = new TProfile(Form("pNEcalClustersMatchedToGen_vsEne_%d", it.first),Form("pNEcalClustersMatchedToGen_vsEne_%d", it.first), 50, 0, 100);
      
      hEffGenMatchedToHcalCluster[it.first] = new TH1F(Form("hEffGenMatchedToHcalCluster_%d", it.first), Form("hEffGenMatchedToHcalCluster_%d", it.first), NBIN, minEff, maxEff);
      hNHcalClustersMatchedToGen[it.first]  = new TH1F(Form("hNHcalClustersMatchedToGen_%d", it.first),Form("hNHcalClustersMatchedToGen_%d", it.first), 10, -0.5, 9.5);
      pNHcalClustersMatchedToGen_vsEne[it.first] = new TProfile(Form("pNHcalClustersMatchedToGen_vsEne_%d", it.first),Form("pNHcalClustersMatchedToGen_vsEne_%d", it.first), 50, 0, 100);
      
      hEffGenMatchedToCaloCluster[it.first] = new TH1F(Form("hEffGenMatchedToCaloCluster_%d", it.first), Form("hEffGenMatchedToCaloCluster_%d", it.first), NBIN, minEff, maxEff);
      hNCaloClustersMatchedToGen[it.first]  = new TH1F(Form("hNCaloClustersMatchedToGen_%d", it.first),Form("hNCaloClustersMatchedToGen_%d", it.first), 10, -0.5, 9.5);
      pNCaloClustersMatchedToGen_vsEne[it.first] = new TProfile(Form("pNCaloClustersMatchedToGen_vsEne_%d", it.first),Form("pNCaloClustersMatchedToGen_vsEne_%d", it.first), 50, 0, 100);
      
      h2EneReco_vsEneTruth[it.first] = new TH2F(Form("h2EneReco_vsEneTruth_%d", it.first),Form("h2EneReco_vsEneTruth_%d", it.first), 100, 0, 50, 500, 0, 5);
      
      h2NClustersMatchedToGen[it.first]  = new TH2F(Form("h2NClustersMatchedToGen_%d", it.first),Form("h2NClustersMatchedToGen_%d", it.first), 10, -0.5, 9.5, 10, -0.5, 9.5);
      
      h2EneFracECAL_HCAL[it.first]          = new TH2F(Form("h2EneFracECAL_HCAL_%d", it.first),Form("h2EneFracECAL_HCAL_%d", it.first),             100, 0., 10, 100, 0., 1.2);
      h2EneNarrowFracECAL_HCAL[it.first]    = new TH2F(Form("h2EneNarrowFracECAL_HCAL_%d", it.first),Form("h2EneNarrowFracECAL_HCAL_%d", it.first), 100, 0., 10, 100, 0., 1.2);
      h2EneECALFracECAL_HCAL[it.first]      = new TH2F(Form("h2EneECALFracECAL_HCAL_%d", it.first),Form("h2EneECALFracECAL_HCAL_%d", it.first),     100, 0., 10, 100, 0., 1.2);
      h2EneECALFracECAL_Seed[it.first]      = new TH2F(Form("h2EneECALFracECAL_Seed_%d", it.first),Form("h2EneECALFracECAL_Seed_%d", it.first),     100, 0., 1.2,500, 0., 5);
      
  }
  
  
  TH1F * hNGenMatchedToEcalCluster = new TH1F ("hNGenMatchedToEcalCluster", "hNGenMatchedToEcalCluster", 20, -0.5, 19.5);
  TH1F * hNGenMatchedToHcalCluster = new TH1F ("hNGenMatchedToHcalCluster", "hNGenMatchedToHcalCluster", 20, -0.5, 19.5);
  TH1F * hNGenMatchedToCaloCluster = new TH1F ("hNGenMatchedToCaloCluster", "hNGenMatchedToCaloCluster", 20, -0.5, 19.5);
  TH1F * hFracHcalOnlyClusters = new TH1F ("hFracHcalOnlyClusters", "hFracHcalOnlyClusters", 100, 0, 1);
  
  
  SCEPCal_GeometryHelper myGeometry;
  
  
  TChain * TreeRun = new TChain("B4", "B4");      
  TChain * TruthTree = new TChain("truth", "truth");  
  
  for (int iFile = 0; iFile<NFILES; iFile++)
  {
    std::string fname_reco;
    std::string fname_truth;

    if (output_tag == "hznb" || output_tag == "wwlj" || output_tag == "hzjnbn")
    {
       fname_reco  = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/reco/output_SCEPCal_B0T_%s100k_job_%d.root", output_tag.c_str(), iFile);
       fname_truth = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/mc_truth/B0T/%s100k_job_%d_output_tuple.root", output_tag.c_str(), iFile);
    }
    else 
    {
       fname_reco  = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/reco/output_SCEPCal_B0T_%s_job_%d.root", output_tag.c_str(), iFile);
       fname_truth = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/mc_truth/B0T/%s_job_%d_output_tuple.root", output_tag.c_str(), iFile);
//         fname_reco  = Form("../root_files/hep_outputs/output_SCEPCal_B0T_%s_job_%d.root", output_tag.c_str(), iFile);
//         fname_truth = Form("../../HepMC_Files/B0T/%s_job_%d_output_tuple.root", output_tag.c_str(), iFile);

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
  
  
  TFile* outputFile = new TFile(Form("../CNN_trees/output_JetClusters_forCNN_%s.root", output_tag.c_str()),"RECREATE");
  TTree* outputTree = new TTree("CNN", "Tree for CNN");
  myCNNTreeVars cnnTV;
  InitCNNTree (outputTree, cnnTV);
  
  
  
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
          if (this_scint/drh_S_norm>ene_HC_th/4)
          {
              CalHit new_hit;
              new_hit.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_hit.SetSide(1);
              myHcHits.push_back(new_hit);
//               std::cout << "passing ene hcal cut" << std::endl;
          }
          if (this_scint/drh_S_norm>HC_seed_th/4)
          {
              CalSeed new_seed;
              new_seed.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_seed.SetSide(1);
              myHcSeeds.push_back(new_seed);
          }
      }
      
//       std::cout << "Number of HCAL seeds found: " << myHcSeeds.size() << std::endl;
//       std::cout << "Cleaning up HCAL seeds too close to each other" << std::endl;

      std::vector<CalSeed>  myHcSeedsCleaned = CleanSeeds(myHcSeeds, maxDeltaRSeedHcal);                        
      std::vector<CalSeed>  myHcSuperSeeds   = MakeSuperSeeds(myHcSeedsCleaned, myHcHits, maxDeltaRSeedHcal, HC_seed_th);
      hNHcalSeeds->Fill(myHcSuperSeeds.size());
      
      
//       std::cout << "Matching HCAL clusters with gen level" << std::endl;
      for (long unsigned int iseed = 0; iseed < myHcSuperSeeds.size(); iseed++)
      { 
          CalSeed this_seed = myHcSuperSeeds.at(iseed);
          float seed_theta = this_seed.GetWeighedTheta();
          float seed_phi   = this_seed.GetWeighedPhi();
          
          int nGenMatchedToHcalCluster = 0;
          for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
          {
              
              int    pdgId = myTruthTV.mcs_pdgId->at(i);
              double mc_ene   = myTruthTV.mcs_E->at(i);
              if (mc_ene<MC_ene_th) continue;
              if (pdgId == 12 || pdgId == 14 || pdgId == 16) continue; //ignore neutrinos
              
              double truth_phi   = myTruthTV.mcs_phi->at(i);
              double eta   = myTruthTV.mcs_eta->at(i);
              double truth_theta = 2*atan(exp(-eta));
              truth_theta = M_PI- truth_theta;
              
              
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatchHcal)
              {
//                   std::cout  << "HCAL cluster (seedEne = "<< this_seed.GetEne() << " GeV) " << iseed << " matched to MC truth gen level particle " << pdgId << " (energy = " << ene << " GeV)" << std::endl;
                  this_seed.AddGenMatch(pdgId);
                  nGenMatchedToHcalCluster++;
              }
            }
            
            hNGenMatchedToHcalCluster->Fill(nGenMatchedToHcalCluster);   
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
          double ecal_S = gRandom->Poisson(this_ene*LO)/ecal_S_norm/LO;
          this_ene = ecal_S;
          if (this_ene>ene_EC_th/4)
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
          if (this_ene>EC_seed_th/4)
          {
              CalSeed new_seed;
              new_seed.Init(myTV.VecHit_CrystalID->at(i), this_theta, this_phi, this_ene);
              myEcSeeds.push_back(new_seed);
          }
      }
//       std::cout << "Number of ECAL seeds found: " << myEcSeeds.size() << std::endl;      
//       std::cout << "Cleaning up ECAL seeds too close to each other" << std::endl;
//       int seed_clust_size = 2;
//       std::vector<CalSeed>  myEcSeedsCleaned = FindSeeds(myEcHits, maxDeltaRSeedEcal, seed_clust_size);

      std::vector<CalSeed>  myEcSeedsCleaned = CleanSeeds(myEcSeeds, maxDeltaRSeedEcal);
      std::vector<CalSeed>  myEcSuperSeeds   = MakeSuperSeeds(myEcSeedsCleaned, myEcHits, maxDeltaRSeedEcal, EC_seed_th);
      hNEcalSeeds->Fill(myEcSuperSeeds.size());
      
      
      

      
      //       std::cout << "Matching ECAL clusters with gen level" << std::endl;
      //Creating calo clusters      
      std::vector<CalCluster> myCalClusters;
      
      for (long unsigned int iseed = 0; iseed < myEcSuperSeeds.size(); iseed++)
      { 
          CalSeed this_seed = myEcSuperSeeds.at(iseed);
          float seed_theta = this_seed.GetWeighedTheta();
          float seed_phi   = this_seed.GetWeighedPhi();
//           if (abs(this_seed.GetEta())>etaAcceptance) continue;
          int nGenMatchedToCaloCluster = 0;
          
          for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
          {
              
              int pdgId = myTruthTV.mcs_pdgId->at(i);                          
              if (pdgId == 12 || pdgId == 14 || pdgId == 16) continue; //ignore neutrinos

              double mc_ene   = myTruthTV.mcs_E->at(i);
              if (mc_ene<MC_ene_th) continue;
              double truth_phi   = myTruthTV.mcs_phi->at(i);
              double eta   = myTruthTV.mcs_eta->at(i);
              double truth_theta = 2*atan(exp(-eta));
              truth_theta = M_PI- truth_theta;
              
              
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatchEcal)
              {
                  this_seed.AddGenMatch(pdgId);
                  this_seed.AddGenEne(mc_ene);                  
                  nGenMatchedToCaloCluster++;
              }
          }
        
        
          CalCluster thisCluster;
          thisCluster.Init(this_seed, maxDeltaRClusterEcal, maxDeltaRSeedHcal, 15, "ecal");
          thisCluster.Clusterize(myEcHits, myHcHits, myEcHitsF, myEcHitsR);
          myCalClusters.push_back(thisCluster);
        
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
                  
                  cnnTV.image_E1[pixel]   = image_E1[pixel];
                  cnnTV.image_E2[pixel]   = image_E2[pixel];
              }
          }
        
          CalSeed clust_seed = thisCluster.GetSeed();
          std::vector<int> truth_in_seed =  clust_seed.GetGenMatch();
          std::vector<float> truth_in_seed_ene =  clust_seed.GetGenEne();
          if (truth_in_seed.size()>0)
          {
              if (truth_in_seed.at(0)==22 && makePlots)
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
              
              
                                          
          }
          if (truth_in_seed.size()==1)
          {
              //for all ecal seeds matched to a truth particle fill the CNN tree
              int this_pdgId = truth_in_seed.at(0);
              auto it = myPdgId.find(abs(this_pdgId));
              if (it== myPdgId.end()) continue; 
              
              cnnTV.PrimaryParticleEnergy      =  truth_in_seed_ene.at(0);
              cnnTV.CNNPrimaryParticleName     =  pdgToName[abs(this_pdgId)];
              cnnTV.theta_seed = seed_theta;
              cnnTV.phi_seed   = seed_phi;              
              outputTree->Fill();
          }
          
          hNGenMatchedToEcalCluster->Fill(nGenMatchedToCaloCluster);                        
      }
      
      int nHcalOnlyCluster = 0;
      for (long unsigned int iseed = 0; iseed < myHcSuperSeeds.size(); iseed++)
      { 
          CalSeed this_seed = myHcSuperSeeds.at(iseed);          
          float seed_theta = this_seed.GetWeighedTheta();
          float seed_phi   = this_seed.GetWeighedPhi();
          bool HcalMatchedToEcalCluster = false;
          
          //check if HCAL cluster is matched to ECAL cluster, if so skip HCAL cluster since already present in CalCluster collection
          for (long unsigned int iseed = 0; iseed < myEcSuperSeeds.size(); iseed++)
          { 
            CalSeed ec_this_seed = myEcSuperSeeds.at(iseed);
            float ec_seed_theta = ec_this_seed.GetWeighedTheta();
            float ec_seed_phi   = ec_this_seed.GetWeighedPhi();
            float dd = sqrt(pow(seed_theta-ec_seed_theta,2) + pow(seed_phi-ec_seed_phi,2));              
            if (dd < maxDeltaRMatchHcal)
            {
                HcalMatchedToEcalCluster = true;
                break;
            }
          }
          if (HcalMatchedToEcalCluster) continue;
          
          for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
          {
              int    pdgId = myTruthTV.mcs_pdgId->at(i);
              if (pdgId == 12 || pdgId == 14 || pdgId == 16) continue; //ignore neutrinos
              double mc_ene   = myTruthTV.mcs_E->at(i);
              if (mc_ene<MC_ene_th) continue;
              double truth_phi   = myTruthTV.mcs_phi->at(i);
              double eta   = myTruthTV.mcs_eta->at(i);
              double truth_theta = 2*atan(exp(-eta));
              truth_theta = M_PI- truth_theta;
              
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatchHcal)
              {
                  this_seed.AddGenMatch(pdgId);                  
                  nHcalOnlyCluster++;
              }
          }
          
          CalCluster thisCluster;
          thisCluster.Init(this_seed, maxDeltaRClusterEcal, maxDeltaRSeedHcal, 15, "hcal");
          thisCluster.Clusterize(myEcHits, myHcHits, myEcHitsF, myEcHitsR);
          myCalClusters.push_back(thisCluster);
      }
      
      for (long unsigned int icluster = 0; icluster < myCalClusters.size(); icluster++)
      {
          CalCluster this_cluster = myCalClusters.at(icluster);
          CalSeed clust_seed = this_cluster.GetSeed();
          std::vector<int> truth_in_seed =  clust_seed.GetGenMatch();
          hNGenMatchedToCaloCluster->Fill(truth_in_seed.size());         
//           std::cout << "n matched to calo clust = " << truth_in_seed.size() << std::endl;
      }
      
      float fracHcalOnlyCluster = (float)nHcalOnlyCluster/(float)myHcSuperSeeds.size();
//       std::cout << " fracHcalOnlyCluster = " << fracHcalOnlyCluster << std::endl;
      hFracHcalOnlyClusters->Fill(fracHcalOnlyCluster);
      
      hNCaloSeeds->Fill(myCalClusters.size());
      
      
      //match of truth to some clusters (both ECAL and HCAL)
      
      std::map<int,int> nEcalClusterMatchedToGen;
      std::map<int,int> nHcalClusterMatchedToGen;
      std::map<int,int> nCaloClusterMatchedToGen;
      std::map<int,int> nTotGen;      
      
      for (auto it : myPdgId)
      {
          nEcalClusterMatchedToGen[it.first] = 0;
          nHcalClusterMatchedToGen[it.first] = 0;
          nTotGen[it.first] = 0;
      }
      
      for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
      {
          int    pdgId = myTruthTV.mcs_pdgId->at(i);
          double ene   = myTruthTV.mcs_E->at(i);
          
          if (ene < MC_ene_th) continue;
          auto it = myPdgId.find(abs(pdgId));
          if (it== myPdgId.end())
          {
              if (abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16) std::cout << "skipping particle not in my pdg id list: " << pdgId << std::endl;
              continue;
          }
          hEneGen[abs(pdgId)]->Fill(ene);
          
          double truth_phi   = myTruthTV.mcs_phi->at(i);
          double eta   = myTruthTV.mcs_eta->at(i);
          if (abs(eta)>etaAcceptance) continue;                    
          double truth_theta = 2*atan(exp(-eta));
          truth_theta = M_PI- truth_theta;
          
          int matchedToEcalCluster = 0;
          int matchedToHcalCluster = 0;
          int matchedToCaloCluster = 0;
          nTotGen[abs(pdgId)]++;
          
          //ECAL match
          for (long unsigned int iseed = 0; iseed < myEcSuperSeeds.size(); iseed++)
          { 
              CalSeed this_seed = myEcSuperSeeds.at(iseed);
              float seed_theta = this_seed.GetWeighedTheta();
              float seed_phi   = this_seed.GetWeighedPhi();          
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatchEcal)   matchedToEcalCluster ++;                            
          }
          hNEcalClustersMatchedToGen[abs(pdgId)]->Fill(matchedToEcalCluster);
          pNEcalClustersMatchedToGen_vsEne[abs(pdgId)]->Fill(ene, matchedToEcalCluster);
          if (matchedToEcalCluster>0) nEcalClusterMatchedToGen[abs(pdgId)]++;
          
          //HCAL match
          for (long unsigned int iseed = 0; iseed < myHcSuperSeeds.size(); iseed++)
          { 
              CalSeed this_seed = myHcSuperSeeds.at(iseed);
              float seed_theta = this_seed.GetWeighedTheta();
              float seed_phi   = this_seed.GetWeighedPhi();         
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaRMatchHcal)   matchedToHcalCluster ++;                            
          }                    
          hNHcalClustersMatchedToGen[abs(pdgId)]->Fill(matchedToHcalCluster);
          pNHcalClustersMatchedToGen_vsEne[abs(pdgId)]->Fill(ene, matchedToHcalCluster);
          if (matchedToHcalCluster>0) nHcalClusterMatchedToGen[abs(pdgId)]++;
          
          h2NClustersMatchedToGen[abs(pdgId)]->Fill(matchedToEcalCluster, matchedToHcalCluster);
          
          
          //Calo Cluster Match
          for (long unsigned int icluster = 0; icluster < myCalClusters.size(); icluster++)
          { 
              CalCluster this_cluster = myCalClusters.at(icluster);
              CalSeed this_seed = this_cluster.GetSeed();
              
              float seed_theta = this_seed.GetWeighedTheta();
              float seed_phi   = this_seed.GetWeighedPhi();
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));
              float this_cluster_deltaR;
              if (this_cluster.GetType() == "ecal") this_cluster_deltaR = maxDeltaRMatchEcal;
              if (this_cluster.GetType() == "hcal") this_cluster_deltaR = maxDeltaRMatchHcal;
              if (dd < this_cluster_deltaR)   
              {
//                   if (abs(pdgId)==22 || abs(pdgId)==130) std::cout << "i: " << i << " :: pdgId (" << pdgId << ") : matched to " << this_cluster.GetType() << " cluster with seed ( "<< seed_theta <<", " << seed_phi << ")" << std::endl;
                  if (this_cluster.GetType() == "ecal") 
                  {
                      h2EneFracECAL_HCAL[abs(pdgId)]->Fill(this_cluster.GetTotEne()/ene, 
                                                           this_cluster.GetEcalClusterEne()/this_cluster.GetTotEne());
                      
                      h2EneNarrowFracECAL_HCAL[abs(pdgId)]->Fill(this_cluster.GetTotEneNarrow()/ene, 
                                                                 this_cluster.GetEcalClusterEne()/this_cluster.GetTotEneNarrow());
                      
                      h2EneECALFracECAL_HCAL[abs(pdgId)]->Fill(this_cluster.GetEcalClusterEne()/ene, 
                                                               this_cluster.GetEcalClusterEne()/this_cluster.GetTotEneNarrow());
                      
                      h2EneECALFracECAL_Seed[abs(pdgId)]->Fill(this_seed.GetEne()/this_cluster.GetEcalClusterEne(), 
                                                               this_cluster.GetEcalClusterEneFront()/this_cluster.GetEcalClusterEneRear());
                      
                      
                      hEneTotRes[abs(pdgId)]->Fill((this_cluster.GetTotEne() - ene)/ene);
                      hEneTotNarrowRes[abs(pdgId)]->Fill((this_cluster.GetTotEneNarrow() - ene)/ene);
                      hEneECALRes[abs(pdgId)]->Fill((this_cluster.GetEcalClusterEne() - ene)/ene);
                      h2EneReco_vsEneTruth[abs(pdgId)]->Fill(ene, this_cluster.GetEcalClusterEne()/ene);
                  }
                  matchedToCaloCluster ++;
              }
          }                    
          hNCaloClustersMatchedToGen[abs(pdgId)]->Fill(matchedToCaloCluster);
          pNCaloClustersMatchedToGen_vsEne[abs(pdgId)]->Fill(ene, matchedToCaloCluster);
          if (matchedToCaloCluster>0) nCaloClusterMatchedToGen[abs(pdgId)]++;
          
      }
                  
      for (auto it : myPdgId)
      {
          hNTotGen[it.first]->Fill(nTotGen[it.first]);
          if (nTotGen[it.first]>0)
          {
              float eff = float(nEcalClusterMatchedToGen[it.first])/float(nTotGen[it.first]);
              hEffGenMatchedToEcalCluster[it.first]->Fill(eff);
              
              eff = float(nHcalClusterMatchedToGen[it.first])/float(nTotGen[it.first]);
              hEffGenMatchedToHcalCluster[it.first]->Fill(eff);
              
              eff = float(nCaloClusterMatchedToGen[it.first])/float(nTotGen[it.first]);
              hEffGenMatchedToCaloCluster[it.first]->Fill(eff);
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
  hNTotGen[22]->GetXaxis()->SetTitle("Total number of gen particles");
  hNTotGen[22]->GetYaxis()->SetTitle("Counts/event");
  gPad->SetLogy();
  
  int color_it = 0;
  leg = new TLegend(0.55,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      hNTotGen[it.first]->Draw("same");
      hNTotGen[it.first]->SetLineWidth(2);
      hNTotGen[it.first]->SetLineColor(mycolors[color_it]);
      hNTotGen[it.first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
//       leg->AddEntry(hNTotGen[it.first], it.second.c_str(), "lp");
      leg->AddEntry(hNTotGen[it.first], Form("%s : <N> = %.1f", it.second.c_str(), hNTotGen[it.first]->GetMean()), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNTotGen->SaveAs("plots_pfa/cNTotGen.png");
  
  
  TCanvas * cEneGen = new TCanvas ("cEneGen", "cEneGen", 600, 500);
  cEneGen->cd();
  hEneGen[22]->SetStats(0);
  hEneGen[22]->SetTitle(0);
  hEneGen[22]->Draw();
//   hEneGen[22]->GetXaxis()->SetRangeUser(0, 1.4);
  hEneGen[22]->GetYaxis()->SetRangeUser(1, hEneGen[22]->GetMaximum()*5);
  hEneGen[22]->GetXaxis()->SetTitle("Energy [GeV]");
//   hEneGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.57,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      hEneGen[it.first]->Draw("same");
      hEneGen[it.first]->SetLineWidth(2);
      hEneGen[it.first]->SetLineColor(mycolors[color_it]);
      hEneGen[it.first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(hEneGen[it.first], Form("%s : <E> = %.1f GeV", it.second.c_str(), hEneGen[it.first]->GetMean()), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cEneGen->SaveAs("plots_pfa/cEneGen.png");
  
  //gen to ECAL cluster matches
  TCanvas * cNEcalClustersMatchedToGen = new TCanvas ("cNEcalClustersMatchedToGen", "cNEcalClustersMatchedToGen", 600, 500);
  cNEcalClustersMatchedToGen->cd();
  hNEcalClustersMatchedToGen[22]->SetStats(0);
  hNEcalClustersMatchedToGen[22]->SetTitle(0);  
  hNEcalClustersMatchedToGen[22]->GetXaxis()->SetRangeUser(-0.5, 5.5);
//   hNEcalClustersMatchedToGen[22]->GetYaxis()->SetRangeUser(1, hNEcalClustersMatchedToGen[22]->GetMaximum()*5);
  hNEcalClustersMatchedToGen[22]->GetXaxis()->SetTitle("N of ECAL clusters matched to gen particle");
//   hNEcalClustersMatchedToGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
//   gPad->SetLogy();
  hNEcalClustersMatchedToGen[22]->DrawCopy();
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      
      hNEcalClustersMatchedToGen[it.first]->SetLineWidth(2);
      hNEcalClustersMatchedToGen[it.first]->SetLineColor(mycolors[color_it]);
      hNEcalClustersMatchedToGen[it.first]->SetMarkerColor(mycolors[color_it]);
      hNEcalClustersMatchedToGen[it.first]->DrawCopy("same");
            
      color_it++;
      leg->AddEntry(hNEcalClustersMatchedToGen[it.first], it.second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNEcalClustersMatchedToGen->SaveAs("plots_pfa/cNEcalClustersMatchedToGen.png");

  
  TCanvas * cEffGenMatchedToEcalCluster = new TCanvas ("cEffGenMatchedToEcalCluster", "cEffGenMatchedToEcalCluster", 600, 500);
  cEffGenMatchedToEcalCluster->cd();
  hEffGenMatchedToEcalCluster[22]->SetStats(0);
  hEffGenMatchedToEcalCluster[22]->SetTitle(0);
  hEffGenMatchedToEcalCluster[22]->Draw();
  hEffGenMatchedToEcalCluster[22]->GetXaxis()->SetRangeUser(0, 1.1);
  hEffGenMatchedToEcalCluster[22]->GetYaxis()->SetRangeUser(1, hEffGenMatchedToEcalCluster[22]->GetMaximum()*5);
  hEffGenMatchedToEcalCluster[22]->GetXaxis()->SetTitle("Fraction of gen-match to at least one ECAL cluster");
//   hEffGenMatchedToEcalCluster[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      hEffGenMatchedToEcalCluster[it.first]->Draw("same");
      hEffGenMatchedToEcalCluster[it.first]->SetLineWidth(2);
      hEffGenMatchedToEcalCluster[it.first]->SetLineColor(mycolors[color_it]);
      hEffGenMatchedToEcalCluster[it.first]->SetMarkerColor(mycolors[color_it]);    
      
      color_it++;
      leg->AddEntry(hEffGenMatchedToEcalCluster[it.first], it.second.c_str(), "lp");
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
  
  for (auto it : myPdgId)
  {
      pNEcalClustersMatchedToGen_vsEne[it.first]->Draw("same");
      pNEcalClustersMatchedToGen_vsEne[it.first]->SetLineWidth(2);
      pNEcalClustersMatchedToGen_vsEne[it.first]->SetLineColor(mycolors[color_it]);
      pNEcalClustersMatchedToGen_vsEne[it.first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(pNEcalClustersMatchedToGen_vsEne[it.first], it.second.c_str(), "lp");
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
  
  for (auto it : myPdgId)
  {
      hNHcalClustersMatchedToGen[it.first]->Draw("same");
      hNHcalClustersMatchedToGen[it.first]->SetLineWidth(2);
      hNHcalClustersMatchedToGen[it.first]->SetLineColor(mycolors[color_it]);
      hNHcalClustersMatchedToGen[it.first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(hNHcalClustersMatchedToGen[it.first], it.second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNHcalClustersMatchedToGen->SaveAs("plots_pfa/cNHcalClustersMatchedToGen.png");

  
  TCanvas * cEffGenMatchedToHcalCluster = new TCanvas ("cEffGenMatchedToHcalCluster", "cEffGenMatchedToHcalCluster", 600, 500);
  cEffGenMatchedToHcalCluster->cd();
  hEffGenMatchedToHcalCluster[22]->SetStats(0);
  hEffGenMatchedToHcalCluster[22]->SetTitle(0);
  hEffGenMatchedToHcalCluster[22]->Draw();
  hEffGenMatchedToHcalCluster[22]->GetXaxis()->SetRangeUser(0, 1.1);
  hEffGenMatchedToHcalCluster[22]->GetYaxis()->SetRangeUser(1, hEffGenMatchedToHcalCluster[22]->GetMaximum()*5);
  hEffGenMatchedToHcalCluster[22]->GetXaxis()->SetTitle("Fraction of gen-match to at least one HCAL cluster");
//   hEffGenMatchedToHcalCluster[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      hEffGenMatchedToHcalCluster[it.first]->Draw("same");
      hEffGenMatchedToHcalCluster[it.first]->SetLineWidth(2);
      hEffGenMatchedToHcalCluster[it.first]->SetLineColor(mycolors[color_it]);
      hEffGenMatchedToHcalCluster[it.first]->SetMarkerColor(mycolors[color_it]);    
      
      color_it++;
      leg->AddEntry(hEffGenMatchedToHcalCluster[it.first], it.second.c_str(), "lp");
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
  
  for (auto it : myPdgId)
  {
      pNHcalClustersMatchedToGen_vsEne[it.first]->Draw("same");
      pNHcalClustersMatchedToGen_vsEne[it.first]->SetLineWidth(2);
      pNHcalClustersMatchedToGen_vsEne[it.first]->SetLineColor(mycolors[color_it]);
      pNHcalClustersMatchedToGen_vsEne[it.first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(pNHcalClustersMatchedToGen_vsEne[it.first], it.second.c_str(), "lp");
  }
    
  leg->Draw();
  if (SAVEPLOTS) cNHcalClustersMatchedToGen_vsEne->SaveAs("plots_pfa/cNHcalClustersMatchedToGen_vsEne.png");
  
  
  //matches to Calo Clusters
  
  TCanvas * cNCaloClustersMatchedToGen = new TCanvas ("cNCaloClustersMatchedToGen", "cNCaloClustersMatchedToGen", 600, 500);
  cNCaloClustersMatchedToGen->cd();
  hNCaloClustersMatchedToGen[22]->SetStats(0);
  hNCaloClustersMatchedToGen[22]->SetTitle(0);
  hNCaloClustersMatchedToGen[22]->Draw();
  hNCaloClustersMatchedToGen[22]->GetXaxis()->SetRangeUser(-0.5, 5.5);
//   hNCaloClustersMatchedToGen[22]->GetYaxis()->SetRangeUser(1, hNCaloClustersMatchedToGen[22]->GetMaximum()*5);
  hNCaloClustersMatchedToGen[22]->GetXaxis()->SetTitle("N of CALO clusters matched to gen particle");
//   hNCaloClustersMatchedToGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      hNCaloClustersMatchedToGen[it.first]->Draw("same");
      hNCaloClustersMatchedToGen[it.first]->SetLineWidth(2);
      hNCaloClustersMatchedToGen[it.first]->SetLineColor(mycolors[color_it]);
      hNCaloClustersMatchedToGen[it.first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(hNCaloClustersMatchedToGen[it.first], it.second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cNCaloClustersMatchedToGen->SaveAs("plots_pfa/cNCaloClustersMatchedToGen.png");
  
  TCanvas * cEffGenMatchedToCaloCluster = new TCanvas ("cEffGenMatchedToCaloCluster", "cEffGenMatchedToCaloCluster", 600, 500);
  cEffGenMatchedToCaloCluster->cd();
  hEffGenMatchedToCaloCluster[22]->SetStats(0);
  hEffGenMatchedToCaloCluster[22]->SetTitle(0);
  hEffGenMatchedToCaloCluster[22]->Draw();
  hEffGenMatchedToCaloCluster[22]->GetXaxis()->SetRangeUser(0, 1.1);
  hEffGenMatchedToCaloCluster[22]->GetYaxis()->SetRangeUser(1, hEffGenMatchedToCaloCluster[22]->GetMaximum()*5);
  hEffGenMatchedToCaloCluster[22]->GetXaxis()->SetTitle("Fraction of gen-match to at least one CALO cluster");
//   hEffGenMatchedToCaloCluster[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      hEffGenMatchedToCaloCluster[it.first]->Draw("same");
      hEffGenMatchedToCaloCluster[it.first]->SetLineWidth(2);
      hEffGenMatchedToCaloCluster[it.first]->SetLineColor(mycolors[color_it]);
      hEffGenMatchedToCaloCluster[it.first]->SetMarkerColor(mycolors[color_it]);    
      
      color_it++;
      leg->AddEntry(hEffGenMatchedToCaloCluster[it.first], it.second.c_str(), "lp");
  }    
  leg->Draw();
  if (SAVEPLOTS) cEffGenMatchedToCaloCluster->SaveAs("plots_pfa/cEffGenMatchedToCaloCluster.png");
  
  
  TCanvas * cNCaloClustersMatchedToGen_vsEne = new TCanvas ("cNCaloClustersMatchedToGen_vsEne", "cNCaloClustersMatchedToGen_vsEne", 600, 500);
  cNCaloClustersMatchedToGen_vsEne->cd();
  pNCaloClustersMatchedToGen_vsEne[22]->SetStats(0);
  pNCaloClustersMatchedToGen_vsEne[22]->SetTitle(0);
  pNCaloClustersMatchedToGen_vsEne[22]->Draw();
  pNCaloClustersMatchedToGen_vsEne[22]->GetXaxis()->SetRangeUser(0, 100);
  pNCaloClustersMatchedToGen_vsEne[22]->GetYaxis()->SetRangeUser(0, 2);
  pNCaloClustersMatchedToGen_vsEne[22]->GetXaxis()->SetTitle("MC truth particle energy [GeV]");
  pNCaloClustersMatchedToGen_vsEne[22]->GetYaxis()->SetTitle("<N> of CALO clusters matched to gen particle");
//   pNCaloClustersMatchedToGen_vsEne[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
//   gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it : myPdgId)
  {
      pNCaloClustersMatchedToGen_vsEne[it.first]->Draw("same");
      pNCaloClustersMatchedToGen_vsEne[it.first]->SetLineWidth(2);
      pNCaloClustersMatchedToGen_vsEne[it.first]->SetLineColor(mycolors[color_it]);
      pNCaloClustersMatchedToGen_vsEne[it.first]->SetMarkerColor(mycolors[color_it]);
            
      color_it++;
      leg->AddEntry(pNCaloClustersMatchedToGen_vsEne[it.first], it.second.c_str(), "lp");
  }
    
  leg->Draw();
  if (SAVEPLOTS) cNCaloClustersMatchedToGen_vsEne->SaveAs("plots_pfa/cNCaloClustersMatchedToGen_vsEne.png");
  
  
  
  TCanvas * cNSeeds = new TCanvas ("cNSeeds", "cNSeeds", 600, 500);
  cNSeeds->cd();
  hNEcalSeeds->Draw();
  hNEcalSeeds->SetStats(0);
  hNEcalSeeds->SetLineColor(kGreen+1);
  hNEcalSeeds->SetFillColor(kGreen+1);
  hNEcalSeeds->SetLineWidth(2);
  hNEcalSeeds->GetXaxis()->SetTitle("N of clusters / event");
  hNEcalSeeds->GetYaxis()->SetRangeUser(1, hNEcalSeeds->GetMaximum()*10);
  
  hNHcalSeeds->Draw("same");
  hNHcalSeeds->SetLineColor(kRed+1);
  hNHcalSeeds->SetLineWidth(2);
  hNCaloSeeds->Draw("same");
  hNCaloSeeds->SetLineColor(kBlack);
  hNCaloSeeds->SetLineWidth(2);
  leg = new TLegend(0.65,0.5,0.88,0.88,NULL,"brNDC");
  leg->AddEntry(hNEcalSeeds, "ECAL seeds", "lp");
  leg->AddEntry(hNHcalSeeds, "HCAL seeds", "lp");
  leg->AddEntry(hNCaloSeeds, "CALO seeds", "lp");
  leg->Draw();
  gPad->SetLogy();
  if (SAVEPLOTS) cNSeeds->SaveAs("plots_pfa/cNSeeds.png");
  
  TCanvas * cNGenMatchedToEcalCluster = new TCanvas ("cNGenMatchedToEcalCluster", "cNGenMatchedToEcalCluster", 600, 500);
  hNGenMatchedToEcalCluster->Scale(1./hNGenMatchedToEcalCluster->Integral());
  hNGenMatchedToHcalCluster->Scale(1./hNGenMatchedToHcalCluster->Integral());
  hNGenMatchedToCaloCluster->Scale(1./hNGenMatchedToCaloCluster->Integral());
  cNGenMatchedToEcalCluster->cd();
  hNGenMatchedToEcalCluster->Draw("histo");
  hNGenMatchedToEcalCluster->GetXaxis()->SetTitle("N of gen particles matched to cluster");
  hNGenMatchedToEcalCluster->GetYaxis()->SetTitle("Frequency");
  hNGenMatchedToEcalCluster->SetStats(0);
  hNGenMatchedToEcalCluster->SetLineColor(kGreen+1);
  hNGenMatchedToEcalCluster->SetFillColor(kGreen+1);
  hNGenMatchedToEcalCluster->SetLineWidth(2);  
//   hNGenMatchedToEcalCluster->GetYaxis()->SetRangeUser(1, hNGenMatchedToEcalCluster->GetMaximum()*10);
  hNGenMatchedToEcalCluster->GetYaxis()->SetRangeUser(0, 1.1);
  hNGenMatchedToEcalCluster->GetXaxis()->SetRangeUser(-0.5, 5.5);
  
  hNGenMatchedToHcalCluster->Draw("same histo");
  hNGenMatchedToHcalCluster->SetLineColor(kRed+1);
  hNGenMatchedToHcalCluster->SetLineWidth(2);
  hNGenMatchedToCaloCluster->Draw("same histo");
  hNGenMatchedToCaloCluster->SetLineColor(kBlack);
  hNGenMatchedToCaloCluster->SetLineWidth(2);
  leg = new TLegend(0.65,0.5,0.88,0.88,NULL,"brNDC");
  leg->AddEntry(hNGenMatchedToEcalCluster, "ECAL seeds", "lp");
  leg->AddEntry(hNGenMatchedToHcalCluster, "HCAL seeds", "lp");
  leg->AddEntry(hNGenMatchedToCaloCluster, "CALO clusters", "lp");
  leg->Draw();
//   gPad->SetLogy();
  
  if (SAVEPLOTS) cNGenMatchedToEcalCluster->SaveAs("plots_pfa/cNGenMatchedToEcalCluster.png");
  
  
  std::map<int, TCanvas *> cNClustersMatchedToGenParticle;
  for (auto it : myPdgId)
  {
    cNClustersMatchedToGenParticle[it.first] = new TCanvas (Form("cNClustersMatchedToGen_%d", it.first), Form("cNClustersMatchedToGen_%d", it.first), 500, 500);
    cNClustersMatchedToGenParticle[it.first]->cd();
    
    hNEcalClustersMatchedToGen[it.first]->Sumw2();
    hNHcalClustersMatchedToGen[it.first]->Sumw2();
    hNCaloClustersMatchedToGen[it.first]->Sumw2();
    hNEcalClustersMatchedToGen[it.first]->Scale(1/hNEcalClustersMatchedToGen[it.first]->Integral());
    hNHcalClustersMatchedToGen[it.first]->Scale(1/hNHcalClustersMatchedToGen[it.first]->Integral());
    hNCaloClustersMatchedToGen[it.first]->Scale(1/hNCaloClustersMatchedToGen[it.first]->Integral());     
      
    hNEcalClustersMatchedToGen[it.first]->Draw("histo");
    hNEcalClustersMatchedToGen[it.first]->SetTitle(Form("N of clusters matched to gen %s", it.second.c_str() ) );
    hNEcalClustersMatchedToGen[it.first]->GetXaxis()->SetTitle(Form("N of clusters matched to gen %s", it.second.c_str()));
    hNEcalClustersMatchedToGen[it.first]->GetYaxis()->SetTitle("Frequency");
    hNEcalClustersMatchedToGen[it.first]->SetStats(0);
    hNEcalClustersMatchedToGen[it.first]->SetLineColor(kGreen+1);
    hNEcalClustersMatchedToGen[it.first]->SetFillColor(kGreen+1);
    hNEcalClustersMatchedToGen[it.first]->SetLineWidth(2);  
    hNEcalClustersMatchedToGen[it.first]->GetYaxis()->SetRangeUser(0, 1.1);
    hNEcalClustersMatchedToGen[it.first]->GetXaxis()->SetRangeUser(-0.5, 5.5);
    
    hNHcalClustersMatchedToGen[it.first]->Draw("same histo");
    hNHcalClustersMatchedToGen[it.first]->SetLineColor(kRed+1);
    hNHcalClustersMatchedToGen[it.first]->SetLineWidth(2);
    hNCaloClustersMatchedToGen[it.first]->Draw("same histo");
    hNCaloClustersMatchedToGen[it.first]->SetLineColor(kBlack);
    hNCaloClustersMatchedToGen[it.first]->SetLineWidth(2);
    leg = new TLegend(0.65,0.5,0.88,0.88,NULL,"brNDC");
    leg->AddEntry(hNEcalClustersMatchedToGen[it.first], "ECAL seeds", "lp");
    leg->AddEntry(hNHcalClustersMatchedToGen[it.first], "HCAL seeds", "lp");
    leg->AddEntry(hNCaloClustersMatchedToGen[it.first], "CALO seeds", "lp");
    leg->Draw();
    //   gPad->SetLogy();  
    if (SAVEPLOTS) cNClustersMatchedToGenParticle[it.first]->SaveAs(Form("plots_pfa/cNClustersMatchedToGen_pdgId%d.png", it.first));
  }
  
  
  
  std::map<int, TCanvas *> cCaloClusterEneRes;
  for (auto it : myPdgId)
  {
    cCaloClusterEneRes[it.first] = new TCanvas (Form("cCaloClusterEneRes_%d", it.first), Form("cCaloClusterEneRes_%d", it.first), 600, 500);
    cCaloClusterEneRes[it.first]->cd();
    
    hEneECALRes[it.first]->Sumw2();
    hEneTotRes[it.first]->Sumw2();
    hEneTotNarrowRes[it.first]->Sumw2();
    hEneECALRes[it.first]->Scale(1/hEneECALRes[it.first]->Integral());
    hEneTotRes[it.first]->Scale(1/hEneTotRes[it.first]->Integral());
    hEneTotNarrowRes[it.first]->Scale(1/hEneTotNarrowRes[it.first]->Integral());
    
      
    hEneECALRes[it.first]->Draw("histo");
    hEneECALRes[it.first]->SetTitle(Form("Gen particle %s", it.second.c_str() ) );
    hEneECALRes[it.first]->GetXaxis()->SetTitle("(E_{calo,cluster} - E_{truth}) / E_{truth}");
    hEneECALRes[it.first]->GetYaxis()->SetTitle("Frequency");
//     hEneECALRes[it.first]->SetStats(0);
    hEneECALRes[it.first]->SetLineColor(kGreen+1);
    hEneECALRes[it.first]->SetFillColor(kGreen+1);
    hEneECALRes[it.first]->SetLineWidth(2);  
    hEneECALRes[it.first]->GetYaxis()->SetRangeUser(0.0001, 1.1);
    
    hEneTotNarrowRes[it.first]->Draw("same histo");
    hEneTotNarrowRes[it.first]->SetLineColor(kRed+1);
    hEneTotNarrowRes[it.first]->SetLineWidth(2);
    hEneTotRes[it.first]->Draw("same histo");
    hEneTotRes[it.first]->SetLineColor(kBlack);
    hEneTotRes[it.first]->SetLineWidth(2);
    
    leg = new TLegend(0.65,0.5,0.88,0.88,NULL,"brNDC");
    leg->AddEntry(hEneECALRes[it.first], "ECAL energy", "lp");
    leg->AddEntry(hEneECALRes[it.first], "CALO narrow", "lp");
    leg->AddEntry(hEneECALRes[it.first], "CALO total", "lp");
    
    
    gPad->SetLogy();
    if (SAVEPLOTS) cCaloClusterEneRes[it.first]->SaveAs(Form("plots_pfa/cCaloClusterEneRes_pdgId%d.png", it.first));
  }
  
    
  std::map<int, TCanvas *> cScatterNClustersMatchedToGenParticle;
  for (auto it : myPdgId)
  {
    cScatterNClustersMatchedToGenParticle[it.first] = new TCanvas (Form("cScatterNClustersMatchedToGenParticle_%d", it.first), Form("cScatterNClustersMatchedToGenParticle_%d", it.first), 500, 500);
    cScatterNClustersMatchedToGenParticle[it.first]->cd();        
    
    if (h2NClustersMatchedToGen[it.first]->Integral()>0) h2NClustersMatchedToGen[it.first]->Scale(1./h2NClustersMatchedToGen[it.first]->Integral());
    h2NClustersMatchedToGen[it.first]->Draw("COL TEXT");
    h2NClustersMatchedToGen[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2NClustersMatchedToGen[it.first]->GetXaxis()->SetTitle(Form("N of ECAL clusters matched to gen %s", it.second.c_str()));
    h2NClustersMatchedToGen[it.first]->GetYaxis()->SetTitle(Form("N of HCAL clusters matched to gen %s", it.second.c_str()));
    h2NClustersMatchedToGen[it.first]->SetStats(0);
//     h2NClustersMatchedToGen[it.first]->SetLineColor(kGreen+1);
//     h2NClustersMatchedToGen[it.first]->SetFillColor(kGreen+1);
//     h2NClustersMatchedToGen[it.first]->SetLineWidth(2);  
//     h2NClustersMatchedToGen[it.first]->GetYaxis()->SetRangeUser(0, 1.1);
    h2NClustersMatchedToGen[it.first]->GetXaxis()->SetRangeUser(-0.5, 3.5);
    h2NClustersMatchedToGen[it.first]->GetYaxis()->SetRangeUser(-0.5, 3.5);
    
    //   gPad->SetLogy();  
    if (SAVEPLOTS) cScatterNClustersMatchedToGenParticle[it.first]->SaveAs(Form("plots_pfa/cScatterNClustersMatchedToGenParticle_pdgId%d.png", it.first));
  }
  
  
  
    std::map<int, TCanvas *> cScatterEneFrac_vsCaloFrac;
  for (auto it : myPdgId)
  {
    cScatterEneFrac_vsCaloFrac[it.first] = new TCanvas (Form("cScatterEneFrac_vsCaloFrac_%d", it.first), Form("cScatterEneFrac_vsCaloFrac_%d", it.first), 1600, 400);
    cScatterEneFrac_vsCaloFrac[it.first]->Divide(4,1);
    
    cScatterEneFrac_vsCaloFrac[it.first]->cd(1);
    h2EneFracECAL_HCAL[it.first]->Draw("COLZ");
    h2EneFracECAL_HCAL[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneFracECAL_HCAL[it.first]->GetXaxis()->SetTitle("Ene TOT cluster / Ene MC");
    h2EneFracECAL_HCAL[it.first]->GetYaxis()->SetTitle("Ene ECAL clust / Ene TOT clust");
    h2EneFracECAL_HCAL[it.first]->SetStats(0);
    
    cScatterEneFrac_vsCaloFrac[it.first]->cd(2);
    h2EneNarrowFracECAL_HCAL[it.first]->Draw("COLZ");
    h2EneNarrowFracECAL_HCAL[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneNarrowFracECAL_HCAL[it.first]->GetXaxis()->SetTitle("Ene TOT narrow cluster / Ene MC");
    h2EneNarrowFracECAL_HCAL[it.first]->GetYaxis()->SetTitle("Ene ECAL clust / Ene TOT narrow clust");
    h2EneNarrowFracECAL_HCAL[it.first]->SetStats(0);
    
    cScatterEneFrac_vsCaloFrac[it.first]->cd(3);
    h2EneECALFracECAL_HCAL[it.first]->Draw("COLZ");
    h2EneECALFracECAL_HCAL[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneECALFracECAL_HCAL[it.first]->GetXaxis()->SetTitle("Ene ECAL cluster / Ene MC");
    h2EneECALFracECAL_HCAL[it.first]->GetYaxis()->SetTitle("Ene ECAL clust / Ene TOT narrow clust");
    h2EneECALFracECAL_HCAL[it.first]->SetStats(0);
    
    
    cScatterEneFrac_vsCaloFrac[it.first]->cd(4);
    h2EneECALFracECAL_Seed[it.first]->Draw("COLZ");
    h2EneECALFracECAL_Seed[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneECALFracECAL_Seed[it.first]->GetXaxis()->SetTitle("Ene ECAL seed / Ene ECAL cluster");
    h2EneECALFracECAL_Seed[it.first]->GetYaxis()->SetTitle("Ene ECAL Front / Ene ECAL Rear");
    h2EneECALFracECAL_Seed[it.first]->SetStats(0);
    
    if (SAVEPLOTS) cScatterEneFrac_vsCaloFrac[it.first]->SaveAs(Form("plots_pfa/cScatterEneFrac_vsCaloFrac_pdgId%d.png", it.first));
  }
  
  
  std::map<int, TCanvas *> cScatterEneReco_vs_EneTruth;
  for (auto it : myPdgId)
  {
    cScatterEneReco_vs_EneTruth[it.first] = new TCanvas (Form("cScatterEneReco_vs_EneTruth_%d", it.first), Form("cScatterEneReco_vs_EneTruth_%d", it.first), 500, 500);
//     cScatterEneReco_vs_EneTruth[it.first]->Divide(4,1);
    
//     cScatterEneFrac_vsCaloFrac[it.first]->cd(1);
    h2EneReco_vsEneTruth[it.first]->Draw("COLZ");
    h2EneReco_vsEneTruth[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneReco_vsEneTruth[it.first]->GetXaxis()->SetTitle("Ene MC");
    h2EneReco_vsEneTruth[it.first]->GetYaxis()->SetTitle("Ene ECAL clust / Ene MC");
    h2EneReco_vsEneTruth[it.first]->SetStats(0);
    
    /*
    cScatterEneFrac_vsCaloFrac[it.first]->cd(2);
    h2EneNarrowFracECAL_HCAL[it.first]->Draw("COLZ");
    h2EneNarrowFracECAL_HCAL[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneNarrowFracECAL_HCAL[it.first]->GetXaxis()->SetTitle("Ene TOT narrow cluster / Ene MC");
    h2EneNarrowFracECAL_HCAL[it.first]->GetYaxis()->SetTitle("Ene ECAL clust / Ene TOT narrow clust");
    h2EneNarrowFracECAL_HCAL[it.first]->SetStats(0);
    
    cScatterEneFrac_vsCaloFrac[it.first]->cd(3);
    h2EneECALFracECAL_HCAL[it.first]->Draw("COLZ");
    h2EneECALFracECAL_HCAL[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneECALFracECAL_HCAL[it.first]->GetXaxis()->SetTitle("Ene ECAL cluster / Ene MC");
    h2EneECALFracECAL_HCAL[it.first]->GetYaxis()->SetTitle("Ene ECAL clust / Ene TOT narrow clust");
    h2EneECALFracECAL_HCAL[it.first]->SetStats(0);
    
    
    cScatterEneFrac_vsCaloFrac[it.first]->cd(4);
    h2EneECALFracECAL_Seed[it.first]->Draw("COLZ");
    h2EneECALFracECAL_Seed[it.first]->SetTitle(Form("Gen particle: %s", it.second.c_str() ) );
    h2EneECALFracECAL_Seed[it.first]->GetXaxis()->SetTitle("Ene ECAL seed / Ene ECAL cluster");
    h2EneECALFracECAL_Seed[it.first]->GetYaxis()->SetTitle("Ene ECAL Front / Ene ECAL Rear");
    h2EneECALFracECAL_Seed[it.first]->SetStats(0);*/
    
    if (SAVEPLOTS) cScatterEneReco_vs_EneTruth[it.first]->SaveAs(Form("plots_pfa/cScatterEneReco_vs_EneTruth_pdgId%d.png", it.first));
  }
  
  TCanvas * cFracHcalOnlyClusters = new TCanvas ("cFracHcalOnlyClusters", "cFracHcalOnlyClusters", 600, 500);
  cFracHcalOnlyClusters->cd();
  hFracHcalOnlyClusters->Draw();
  hFracHcalOnlyClusters->GetXaxis()->SetTitle("Fraction of HCAL only clusters (not matched to ECAL)");
  gPad->SetLogy();
  if (SAVEPLOTS) cFracHcalOnlyClusters->SaveAs("plots_pfa/cFracHcalOnlyClusters.png");
  
  
  
  
  if (WRITEOUTPUT)
  {
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
    std::cout << "done writing output file!" << std::endl;
  }
  
  
  theApp->Run();
}







