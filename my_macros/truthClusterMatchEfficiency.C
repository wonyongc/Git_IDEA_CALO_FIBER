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
  
  
//   bool SAVEPLOTS = false;  
//   int energy = 100;
  
  std::string output_tag = "zjj_scan_100";
  int NFILES = 4;
  if (argc>1) output_tag = argv[1];   
  if (argc>2) NFILES = atoi(argv[2]);   
  std::cout << "processing sample of: " << output_tag.c_str() << std::endl;  

    
  //define histos
  
  
  double drh_S_norm  = 407;
  float maxDeltaR = 0.1;
  float etaAcceptance = 1.4;
  
  float ene_EC_th = 0.01;
  float EC_seed_th = 0.1;
  
  float ene_HC_th   = 0.01;    
  float HC_seed_th = 0.1;
      
  
  
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
  
  std::map <int, TH1F*> hEffGenMatchedToCluster;
  std::map <int, TH1F*> hNClustersMatchedToGen;
  std::map <int, TH1F*> hNTotGen;
  
  
  TH1F * hNEcalSeeds = new TH1F ("hNEcalSeeds", "hNEcalSeeds", 100, -0.5, 99.5);
  TH1F * hNHcalSeeds = new TH1F ("hNHcalSeeds", "hNHcalSeeds", 100, -0.5, 99.5);
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hEffGenMatchedToCluster[it->first] = new TH1F(Form("hEffGenMatchedToCluster_%d", it->first), Form("hEffGenMatchedToCluster_%d", it->first), NBIN, minEff, maxEff);
      hNClustersMatchedToGen[it->first] = new TH1F(Form("hNClustersMatchedToGen_%d", it->first),Form("hNClustersMatchedToGen_%d", it->first), 10, -0.5, 9.5);
      hNTotGen[it->first] = new TH1F(Form("hNTotGen_%d", it->first),Form("hNTotGen_%d", it->first), 50, -0.5, 49.5);
  }
  
  
  TH1F * hNGenMatchedToCluster = new TH1F ("hNGenMatchedToCluster", "hNGenMatchedToCluster", 20, -0.5, 19.5);
  
  
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
      
      for (unsigned int i = 0; i<myTV.VectorL->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'l');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorL->at(i)/1000.;      
          double this_scint = myTV.VectorSignalsL->at(i);                
          double this_cher  = myTV.VectorSignalsCherL->at(i);
          if (this_ene>ene_HC_th)
          {
              hGrid_DRT_S ->Fill(this_theta, this_phi, this_scint/drh_S_norm);                              
              hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher/drh_C_norm);    
              CalHit new_hit;
              new_hit.Init(i, this_theta, this_phi, this_ene);
              new_hit.SetSide(-1);
              myHcHits.push_back(new_hit);
          }
          totDRHEne+=this_ene;
          totDRHScint+=this_scint;
          totDRHCher+=this_cher;
          
          if (this_scint/drh_S_norm>HC_seed_th)
          {
              CalSeed new_seed;
              new_seed.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_seed.SetSide(-1);              
              myHcSeeds.push_back(new_seed);
          }
      }
      for (unsigned int i = 0; i<myTV.VectorR->size(); i++)
      {                                        
          TVector3 this_vec = myGeometry.GetTowerVec(i,'r');
          double this_phi   = this_vec.Phi();
          double this_theta = this_vec.Theta();
          double this_ene   = myTV.VectorR->at(i)/1000.;
          double this_scint = myTV.VectorSignalsR->at(i);
          double this_cher  = myTV.VectorSignalsCherR->at(i);
          if (this_ene>ene_HC_th)
          {
              hGrid_DRT_S ->Fill(this_theta, this_phi, this_scint/drh_S_norm);
              hGrid_DRT_C ->Fill(this_theta, this_phi, this_cher/drh_C_norm);
              CalHit new_hit;
              new_hit.Init(i, this_theta, this_phi, this_ene);
              new_hit.SetSide(1);
              myHcHits.push_back(new_hit);
          }
          totDRHEne+=this_ene;
          totDRHScint+=this_scint;
          totDRHCher+=this_cher;
          if (this_scint/drh_S_norm>HC_seed_th)
          {
              CalSeed new_seed;
              new_seed.Init(i, this_theta, this_phi, this_scint/drh_S_norm);
              new_seed.SetSide(1);
              myHcSeeds.push_back(new_seed);
          }
      }
      
      std::cout << "Total ene in scint HCAL fiber: " << totDRHEne << std::endl;
      std::cout << "Total scint in HCAL: " << totDRHScint << std::endl;
      std::cout << "Total cher in HCAL: " << totDRHCher << std::endl;
      
      std::cout << "Number of HCAL seeds found: " << myHcSeeds.size() << std::endl;
      std::cout << "Cleaning up HCAL seeds too close to each other" << std::endl;
      std::vector<CalSeed>  myHcSeedsCleaned = CleanSeeds(myHcSeeds, maxDeltaR);
                        
      std::cout << "Matching HCAL clusters with gen level" << std::endl;
      for (long unsigned int iseed = 0; iseed < myHcSeedsCleaned.size(); iseed++)
      { 
          CalSeed this_seed = myHcSeedsCleaned.at(iseed);
          float seed_theta = this_seed.GetTheta();
          float seed_phi   = this_seed.GetPhi();
          
          
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
                  std::cout  << "HCAL cluster (seedEne = "<< this_seed.GetEne() << " GeV) " << iseed << " matched to MC truth gen level particle " << pdgId << " (energy = " << ene << " GeV)" << std::endl;
                  this_seed.AddGenMatch(pdgId);
              }
            }                    
      }
      
      
      
      //**************************************************************//
      //                             ECAL
      //**************************************************************//


      
      std::vector<CalHit> myEcHits;
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
      std::vector<CalSeed>  myEcSeedsCleaned = CleanSeeds(myEcSeeds, maxDeltaR);
      
    

      hNEcalSeeds->Fill(myEcSeedsCleaned.size());
//       std::cout << "Matching gen level with ECAL clusters" << std::endl;
      
      
      std::map<int,int> nGenMatchToCluster;
      std::map<int,int> nTotGen;
      
      
      for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
      {
          nGenMatchToCluster[it->first] = 0;
          nTotGen[it->first] = 0;
      }
      
      for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
      {
          int    pdgId = myTruthTV.mcs_pdgId->at(i);
          double ene   = myTruthTV.mcs_E->at(i);
          if (ene<EC_seed_th) continue;
          auto it = myPdgId.find(abs(pdgId));
          if (it== myPdgId.end())
          {
//               std::cout << "skipping particle non in my pdg id list: " << pdgId << std::endl;
              continue;
          }
          
          double truth_phi   = myTruthTV.mcs_phi->at(i);
          double eta   = myTruthTV.mcs_eta->at(i);
          if (abs(eta)>etaAcceptance) continue;                    
          double truth_theta = 2*atan(exp(-eta));
          truth_theta = M_PI- truth_theta;
          
          int matchedToCluster = 0;
          nTotGen[abs(pdgId)]++;
          
          
          for (long unsigned int iseed = 0; iseed < myEcSeedsCleaned.size(); iseed++)
          { 
              CalSeed this_seed = myEcSeedsCleaned.at(iseed);
              float seed_theta = this_seed.GetTheta();
              float seed_phi   = this_seed.GetPhi();
            
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaR)
              {                  
//                   std::cout << "MC truth gen level particle " << pdgId << " matched to ECAL cluster " << iseed << std::endl;
                  matchedToCluster ++;
              }
          }
          
          hNClustersMatchedToGen[abs(pdgId)]->Fill(matchedToCluster);
          if (matchedToCluster>=1) nGenMatchToCluster[abs(pdgId)]++;
      }
                  
      for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
      {
          hNTotGen[it->first]->Fill(nTotGen[it->first]);
          if (nTotGen[it->first]>0)
          {
              float eff = float(nGenMatchToCluster[it->first])/float(nTotGen[it->first]);
              hEffGenMatchedToCluster[it->first]->Fill(eff);
          }
      }

      
      //       std::cout << "Matching ECAL clusters with gen level" << std::endl;
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

              double truth_phi   = myTruthTV.mcs_phi->at(i);
              double eta   = myTruthTV.mcs_eta->at(i);              
              double truth_theta = 2*atan(exp(-eta));
              truth_theta = M_PI- truth_theta;
              
              
              float dd = sqrt(pow(seed_theta-truth_theta,2) + pow(seed_phi-truth_phi,2));              
              if (dd < maxDeltaR)
              {
//                   CalCluster thisCluster;
//                   thisCluster.Init(this_seed, maxDeltaR);
//                   thisCluster.Clusterize(myEcHits, myHcHits);
                  
//                   std::cout  << "ECAL cluster (seedEne = "<< this_seed.GetEne() << " GeV, clusterEcalEne = " << thisCluster.GetEcalClusterEne() << " GeV, clusterTotEne = " << thisCluster.GetTotEne() << " GeV) "  << iseed << " matched to MC truth gen level particle " << pdgId << " (energy = " << ene << " GeV)" <<  std::endl;
//                   this_seed.AddGenMatch(pdgId);
                  nGenMatchedToCluster++;
              }
            }
            
            hNGenMatchedToCluster->Fill(nGenMatchedToCluster);                        
      }            
  }
  
  
  
  TCanvas * cNClustersMatchedToGen = new TCanvas ("cNClustersMatchedToGen", "cNClustersMatchedToGen", 600, 500);
  cNClustersMatchedToGen->cd();
  hNClustersMatchedToGen[22]->SetStats(0);
  hNClustersMatchedToGen[22]->SetTitle(0);
  hNClustersMatchedToGen[22]->Draw();
  hNClustersMatchedToGen[22]->GetXaxis()->SetRangeUser(-0.5, 5.5);
//   hNClustersMatchedToGen[22]->GetYaxis()->SetRangeUser(1, hNClustersMatchedToGen[22]->GetMaximum()*5);
  hNClustersMatchedToGen[22]->GetXaxis()->SetTitle("N of ECAL clusters matched to gen particle");
//   hNClustersMatchedToGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  int color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hNClustersMatchedToGen[it->first]->Draw("same");
      hNClustersMatchedToGen[it->first]->SetLineWidth(2);
      hNClustersMatchedToGen[it->first]->SetLineColor(mycolors[color_it]);
      hNClustersMatchedToGen[it->first]->SetMarkerColor(mycolors[color_it]);
      
      
      color_it++;
      leg->AddEntry(hNClustersMatchedToGen[it->first], it->second.c_str(), "lp");
  }
    
  leg->Draw();

  
  TCanvas * cEffGenMatchedToCluster = new TCanvas ("cEffGenMatchedToCluster", "cEffGenMatchedToCluster", 600, 500);
  cEffGenMatchedToCluster->cd();
  hEffGenMatchedToCluster[22]->SetStats(0);
  hEffGenMatchedToCluster[22]->SetTitle(0);
  hEffGenMatchedToCluster[22]->Draw();
  hEffGenMatchedToCluster[22]->GetXaxis()->SetRangeUser(0, 1.4);
  hEffGenMatchedToCluster[22]->GetYaxis()->SetRangeUser(1, hEffGenMatchedToCluster[22]->GetMaximum()*5);
  hEffGenMatchedToCluster[22]->GetXaxis()->SetTitle("Fraction of gen-match to at least one ECAL cluster");
//   hEffGenMatchedToCluster[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
  leg = new TLegend(0.75,0.5,0.88,0.88,NULL,"brNDC");
  
  for (auto it = myPdgId.begin(); it != myPdgId.end(); ++it)
  {
      hEffGenMatchedToCluster[it->first]->Draw("same");
      hEffGenMatchedToCluster[it->first]->SetLineWidth(2);
      hEffGenMatchedToCluster[it->first]->SetLineColor(mycolors[color_it]);
      hEffGenMatchedToCluster[it->first]->SetMarkerColor(mycolors[color_it]);
      
      
      color_it++;
      leg->AddEntry(hEffGenMatchedToCluster[it->first], it->second.c_str(), "lp");
  }
    
  leg->Draw();
  
  
  TCanvas * cNTotGen = new TCanvas ("cNTotGen", "cNTotGen", 600, 500);
  cNTotGen->cd();
  hNTotGen[22]->SetStats(0);
  hNTotGen[22]->SetTitle(0);
  hNTotGen[22]->Draw();
//   hNTotGen[22]->GetXaxis()->SetRangeUser(0, 1.4);
  hNTotGen[22]->GetYaxis()->SetRangeUser(1, hNTotGen[22]->GetMaximum()*5);
  hNTotGen[22]->GetXaxis()->SetTitle("Fraction of gen-match to at least one ECAL cluster");
//   hNTotGen[22]->GetYaxis()->SetTitle("Frequency [a.u.]");
  gPad->SetLogy();
  
  color_it = 0;
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
  
  
  TCanvas * cNEcalSeeds = new TCanvas ("cNEcalSeeds", "nGenMatchedToCluster", 600, 500);
  cNEcalSeeds->cd();
  hNEcalSeeds->Draw();
  hNEcalSeeds->GetXaxis()->SetTitle("N of ECAL clusters");
  
  
  TCanvas * cNGenMatchedToCluster = new TCanvas ("cNGenMatchedToCluster", "nGenMatchedToCluster", 600, 500);
  cNGenMatchedToCluster->cd();
  hNGenMatchedToCluster->Draw();
  hNGenMatchedToCluster->GetXaxis()->SetTitle("N of gen particles matched to ECAL cluster");
  
  
  
  
  
  theApp->Run();
}







