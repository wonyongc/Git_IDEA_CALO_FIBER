// g++ -Wall -o physicsJetAnalysis physicsJetAnalysis.C  myG4Tree.cc myG4Tree.hh myTruthTree.cc myTruthTree.hh recoUtils.cc recoUtils.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs` `~/fastjet-3.3.2-install/bin/fastjet-config --cxxflags --libs --plugins`

// g++ -Wall -o physicsJetAnalysis physicsJetAnalysis.C  myG4Tree.cc myG4Tree.hh myTruthTree.cc myTruthTree.hh recoUtils.cc recoUtils.hh SCEPCal_GeometryHelper.cc SCEPCal_GeometryHelper.hh `root-config --cflags --glibs` `//afs/cern.ch/work/m/mlucchin//fastjet-3.3.2-install/bin/fastjet-config --cxxflags --libs --plugins`

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

#include "fastjet/ClusterSequence.hh"                                                                                                                                                                                   
#include <iostream>                                                                                                                                                                                                     

using namespace std;
using namespace fastjet; 


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
          
  //set Root style
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
    
  
  
  //init  
  bool SAVEPLOTS = false;
  std::string output_tag = "wwlj";
  if (argc>0) output_tag = argv[1];   
  
  double ecal_S_norm = 0.985;
  double ecal_C_norm = 7286;
  float LO  = 2000;
  float CLO = 160;
  double drh_S_norm  = 407;
  double drh_C_norm  = 103.2;  
  
  double x_factor_hcal = 0.43;
  double x_factor_ecal = 0.371;
  
  
// choose a jet definition
  double R = 2*M_PI;
  int nExpJets = 2;
//       JetDefinition jet_def(antikt_algorithm, R);
  JetDefinition jet_def(ee_genkt_algorithm, R, 1);
  std::cout << "Clustering with " << jet_def.description() << std::endl;
  
    // single truth jet
  JetDefinition jet_mc(ee_genkt_algorithm, 4*M_PI, 1);
  
  SCEPCal_GeometryHelper myGeometry;

  
  int NFILES = 100;
  TChain * TreeRun = new TChain("B4", "B4");      
  TChain * TruthTree = new TChain("truth", "truth");  
  
  for (int iFile = 0; iFile<NFILES; iFile++)
  {
    //    std::string fname_reco  = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/reco/output_SCEPCal_B0T_%s100k_job_%d.root", output_tag.c_str(), iFile);
    //    std::string fname_truth = Form("/eos/user/m/mlucchin/WORKAREA/SCEPCal_IDEA_Samples/hep_outputs/mc_truth/B0T/%s100k_job_%d_output_tuple.root", output_tag.c_str(), iFile);

    std::string fname_reco  = Form("../root_files/hep_outputs/output_SCEPCal_B0T_%s100k_job_%d.root", output_tag.c_str(), iFile);
    std::string fname_truth = Form("../../HepMC_Files/B0T/%s100k_job_%d_output_tuple.root", output_tag.c_str(), iFile);

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
  
  
  
  int flag_MCT = 0;
  int flag_JHS = 1;
  int flag_JHC = 2;
  int flag_JES = 3;
  int flag_JEC = 4;
  
  bool debugMode = false;
  
  //define histos
  
  int NBIN = 150;
  int minMass = 0;
  int maxMass = 150;
  
  TH1F * hMCT_MassJJ = new TH1F ("hMCT_MassJJ", "hMCT_MassJJ", NBIN, minMass, maxMass);
  TH1F * hRAW_MassJJ = new TH1F ("hRAW_MassJJ", "hRAW_MassJJ", NBIN, minMass, maxMass);
  TH1F * hDRO_MassJJ = new TH1F ("hDRO_MassJJ", "hDRO_MassJJ", NBIN, minMass, maxMass);
  TH1F * hRAW_MassDiff = new TH1F ("hRAW_MassDiff", "hRAW_MassDiff", NBIN, -1, 1);
  TH1F * hDRO_MassDiff = new TH1F ("hDRO_MassDiff", "hDRO_MassDiff", NBIN, -1, 1);

  TH2F * hRAW_ScatterEne = new TH2F ("hRAW_ScatterEne", "hRAW_ScatterEne", NBIN, -75, 75, NBIN, -75, 75);
  TH2F * hDRO_ScatterEne = new TH2F ("hDRO_ScatterEne", "hDRO_ScatterEne", NBIN, -75, 75, NBIN, -75, 75);
    
  ///*******************************************///
  ///		 Run over events	        ///
  ///*******************************************///
  int maxEVENTS = 10000;
  int NEVENTS = TreeRun->GetEntries();
  
  std::cout << "NEVENTS = " << NEVENTS << std::endl;
  if (NEVENTS>maxEVENTS)  NEVENTS = maxEVENTS;
  std::cout << "... running on " << NEVENTS << " events" << std::endl;
  
  
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
      
    std::cout << "processing event: " << iEvt << "\r" << std::flush;
    std::vector<PseudoJet> allHitsForJet;
      
    bool goodEvent = true;
    int nMuons = 0;
    int nNeutrinos = 0;
    float neutrinoEne = 0;
      
      //filling truth
    TruthTree->GetEntry(iEvt);
//       std::cout << "n MC truth particles found: " << myTruthTV.mcs_n << std::endl;
//       std::cout << "\n*******************************************************\n" << std::endl;
              
    for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
    {          
          int    pdgId = myTruthTV.mcs_pdgId->at(i);
          double ene   = myTruthTV.mcs_E->at(i);
          double phi   = myTruthTV.mcs_phi->at(i);
          double eta   = myTruthTV.mcs_eta->at(i);
          double pT    = myTruthTV.mcs_pt->at(i);
//           int charge = myTruthTV.mcs_charge->at(i);
          double theta = 2*atan(exp(-eta));
          theta = M_PI- theta;
                      
//           if (debugMode) std::cout << "pdgId[" << i << "] = " << pdgId << std::endl;
          
          if (   fabs(pdgId)!=12 && fabs(pdgId)!=14 && fabs(pdgId)!=16 && fabs(pdgId)!=13  // exclude neutrinos and muons
                && fabs(pdgId)<10000 // exclude BSM
              //           && fabs(eta)<etaAcceptance                    // make sure jets are fully contained in calorimeters
          )      
          {
              float px = pT*cos(phi);
              float py = pT*sin(phi);
              float pz = pT*sinh(eta);                    
              
              PseudoJet this_MCT = PseudoJet(px*1e-19, py*1e-19, pz*1e-19, ene*1e-19);        
              this_MCT.set_user_index(flag_MCT);                        
              allHitsForJet.push_back(this_MCT);    
                                          
          }
          
          else
          {
              if (fabs(pdgId)==13) nMuons++;
              if (fabs(pdgId)==12 || fabs(pdgId)==14 || fabs(pdgId)==16) nNeutrinos++;
              if (fabs(pdgId)==14) 
              {
                  neutrinoEne += ene;
		  //                  std::cout << "neutrino Ene = " <<  ene << std::endl;
              }
          }
    }
  
    if (output_tag == "wwln" && (nMuons>1 || nNeutrinos >1))
    {
        goodEvent = false;
        continue;
    }
  
    if (output_tag == "hzjnbn" && (nMuons>0 || nNeutrinos >0))
    {
        goodEvent = false;
        if (debugMode) std::cout << Form(" skipping %s event with %d muons and %d neutrinos ", output_tag.c_str(), nMuons, nNeutrinos) << std::endl;
        continue;
    } 
                                        
     
    //filling reco

    TreeRun->GetEntry(iEvt);

    //**************************************************************//
    //                           DR HCAL
    //**************************************************************//
    
    if (output_tag == "wwln" && (myTV.leakage/1000. - neutrinoEne > 1))
    {
        if (debugMode)        std::cout << "Leakage = " << myTV.leakage/1000. - neutrinoEne << " GeV " << std::endl;
        goodEvent = false;
        continue;
    }
    
    float ene_HC_th   = 0.01;    
    float totS = 0;
    
    
    for (unsigned int i = 0; i<myTV.VectorL->size(); i++)
    {                                        
        TVector3 this_vec = myGeometry.GetTowerVec(i,'l');
//         double this_phi   = this_vec.Phi();
//         double this_theta = this_vec.Theta();
        double this_ene   = myTV.VectorL->at(i)/1000.;      
        double this_scint = myTV.VectorSignalsL->at(i);                
        double this_cher  = myTV.VectorSignalsCherL->at(i);
        double S = this_scint/drh_S_norm;
        double C = this_cher/drh_C_norm;
        
        if (this_ene>ene_HC_th)
        {            
            PseudoJet this_JHS = PseudoJet(this_vec.X()*S, this_vec.Y()*S, this_vec.Z()*S, S);
            this_JHS.set_user_index(flag_JHS);
            allHitsForJet.push_back(this_JHS);
            
            PseudoJet this_JHC = PseudoJet(this_vec.X()*C, this_vec.Y()*C, this_vec.Z()*C, C);
            this_JHC.set_user_index(flag_JHC);
            allHitsForJet.push_back(this_JHC);            
        }          
        totS+=S;        
    }
    
    for (unsigned int i = 0; i<myTV.VectorR->size(); i++)
    {                                        
        TVector3 this_vec = myGeometry.GetTowerVec(i,'r');
//         double this_phi   = this_vec.Phi();
//         double this_theta = this_vec.Theta();
        double this_ene   = myTV.VectorR->at(i)/1000.;
        double this_scint = myTV.VectorSignalsR->at(i);
        double this_cher  = myTV.VectorSignalsCherR->at(i);
        double S = this_scint/drh_S_norm;
        double C = this_cher/drh_C_norm;
        
        if (this_ene>ene_HC_th)
        {
            PseudoJet this_JHS = PseudoJet(this_vec.X()*S, this_vec.Y()*S, this_vec.Z()*S, S);
            this_JHS.set_user_index(flag_JHS);
            allHitsForJet.push_back(this_JHS);
            
            PseudoJet this_JHC = PseudoJet(this_vec.X()*C, this_vec.Y()*C, this_vec.Z()*C, C);
            this_JHC.set_user_index(flag_JHC);
            allHitsForJet.push_back(this_JHC);            
        }                
        totS+=S;        
    }
        
    //    std::cout << "Total S in HCAL: " << totS <<  " GeV " << std::endl;
    
      
      
    
    //**************************************************************//    
    //                             ECAL
    //**************************************************************//

      
    float ene_EC_th = 0.01;    
    float totEcalEne = 0;
            
    for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
    {
     
        TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
//         double this_phi = this_vec.Phi();
//         double this_theta = this_vec.Theta();
        double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;                    
        
        if (this_ene>ene_EC_th)
        {
            double ecal_S = gRandom->Poisson(this_ene*LO)/ecal_S_norm/LO;
            double cherR  = myTV.VecHit_ScepCherR->at(i);
            double ecal_C = gRandom->Poisson(cherR*CLO/ecal_C_norm)/CLO;
            
            PseudoJet this_JES = PseudoJet(this_vec.X()*ecal_S, this_vec.Y()*ecal_S, this_vec.Z()*ecal_S, ecal_S);
            this_JES.set_user_index(flag_JES);
            allHitsForJet.push_back(this_JES);

            PseudoJet this_JEC = PseudoJet(this_vec.X()*ecal_C, this_vec.Y()*ecal_C, this_vec.Z()*ecal_C, ecal_C);
            this_JEC.set_user_index(flag_JEC);
            allHitsForJet.push_back(this_JEC);
        }
        
        totEcalEne+=this_ene;
      }
            
    //      std::cout << "Total energy in ECAL: " << totEcalEne << std::endl;     
//       std::cout << "Running fastjet for clustering calo hits..." << std::endl;            

      // run the clustering, extract the jets
      ClusterSequence cs(allHitsForJet, jet_def);
      
//       std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      std::vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(nExpJets));
      std::vector<PseudoJet> mct_jets;
      std::vector<PseudoJet> raw_jets;
      std::vector<PseudoJet> dro_jets;


      for (unsigned i = 0; i < jets.size(); i++) 
      {
//           cout << "jet " << i << ": "<< jets[i].pt() << " " << jets[i].rap() << " " << jets[i].phi() << endl;
          
          std::vector<PseudoJet> constituents = jets[i].constituents();
          std::vector<PseudoJet> mct_constituents;
          
          float E_MCT = 0;
          float E_JHS = 0;
          float E_JHC = 0;
          float E_JES = 0;
          float E_JEC = 0;
          
          for (unsigned j = 0; j < constituents.size(); j++) 
          {
//               std::cout << "    constituent " << j << "'s pt: "      << constituents[j].pt()  
//                                                    << " :: E = "     << constituents[j].E() 
//                                                    << " :: theta = " << constituents[j].theta()
//                                                    << " :: eta = "   << constituents[j].eta()
//                                                    << " :: phi = "   << constituents[j].phi()
//                                                    << std::endl;
//               
              if      (constituents[j].user_index() == flag_JES) E_JES += constituents[j].E();
              else if (constituents[j].user_index() == flag_JEC) E_JEC += constituents[j].E();
              else if (constituents[j].user_index() == flag_JHS) E_JHS += constituents[j].E();
              else if (constituents[j].user_index() == flag_JHC) E_JHC += constituents[j].E();
              else if (constituents[j].user_index() == flag_MCT) 
              {
                  mct_constituents.push_back(PseudoJet(constituents[j].px()*1e19,constituents[j].py()*1e19,constituents[j].pz()*1e19,constituents[j].E()*1e19));
                  E_MCT += constituents[j].E()*1e19;              
              }
          }
          
          if (debugMode) std::cout << "E_JES = " << E_JES << " :: E_JEC = " << E_JEC << " :: E_JHS = " << E_JHS << " :: E_JHC = " << E_JHC << std::endl;
          if (mct_constituents.size()>0) 
          {
              ClusterSequence csMC(mct_constituents, jet_mc);
              std::vector<PseudoJet> this_mct_jet = sorted_by_pt(csMC.exclusive_jets(int(1) ) );
              if (this_mct_jet.size()>0) mct_jets.push_back(this_mct_jet[0]);
          }
              
//           PseudoJet this_mct_jet = jets[i]*(E_MCT)/jets[i].E();
//           mct_jets.push_back(this_mct_jet);
          
          PseudoJet this_raw_jet = jets[i]*(E_JES+E_JHS)/jets[i].E();
          raw_jets.push_back(this_raw_jet);          
          
          float E_JE   = (E_JES-x_factor_ecal*E_JEC )/(1-x_factor_ecal);
          float E_JH   = (E_JHS-x_factor_hcal*E_JHC )/(1-x_factor_hcal);
          float E_JTot = E_JE+E_JH;
          PseudoJet dro_corr_jet = jets[i]*E_JTot/jets[i].E();
          dro_jets.push_back(dro_corr_jet);
          
      }
      
      if(mct_jets.size()==2)
      {
          //reject jets not fully contained in the calorimeter
          if (fabs(mct_jets[0].eta()) > 2 || fabs(mct_jets[1].eta()) > 2) goodEvent = false;
      }

      float jjMassMCT = 0;
      if (mct_jets.size()==2 && goodEvent)
      {
          float e1 = std::max(mct_jets[0].E(), mct_jets[1].E());
          float e2 = std::min(mct_jets[0].E(), mct_jets[1].E());
          float p1p2Sum = sqrt(pow(mct_jets[0].px()+mct_jets[1].px(),2) + pow(mct_jets[0].py()+mct_jets[1].py(),2) + pow(mct_jets[0].pz()+mct_jets[1].pz(),2) );
          
          jjMassMCT = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hMCT_MassJJ->Fill(jjMassMCT);
	  if (debugMode) std::cout << "MCT: E_j1 + E_j2 = " << e1+e2 << " :: p_j1 + p_j2 = " << p1p2Sum << " :: jjMass = " << jjMassMCT << " GeV" << std::endl;
      }
      float jjMassRAW = 0;
      if (raw_jets.size()==2 && goodEvent)
      {
          
          float e1 = std::max(raw_jets[0].E(), raw_jets[1].E());
          float e2 = std::min(raw_jets[0].E(), raw_jets[1].E());
          float p1p2Sum = sqrt(pow(raw_jets[0].px()+raw_jets[1].px(),2) + pow(raw_jets[0].py()+raw_jets[1].py(),2) + pow(raw_jets[0].pz()+raw_jets[1].pz(),2) );
          
          jjMassRAW = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hRAW_MassJJ->Fill(jjMassRAW);
          hRAW_MassDiff->Fill((jjMassRAW-jjMassMCT)/jjMassMCT);
	  hRAW_ScatterEne->Fill(e1 - std::max(mct_jets[0].E(), mct_jets[1].E()), e2-std::min(mct_jets[0].E(), mct_jets[1].E()) );
	  if (debugMode) std::cout << "RAW: E_j1 + E_j2 = " << e1+e2 << " :: p_j1 + p_j2 = " << p1p2Sum << " :: jjMass = " << jjMassRAW << " GeV" << std::endl;
      }
      float jjMassDRO = 0;
      if (dro_jets.size()==2 && goodEvent)
      {
          float e1 = std::max(dro_jets[0].E(), dro_jets[1].E());
          float e2 = std::min(dro_jets[0].E(), dro_jets[1].E());
          float p1p2Sum = sqrt(pow(dro_jets[0].px()+dro_jets[1].px(),2) + pow(dro_jets[0].py()+dro_jets[1].py(),2) + pow(dro_jets[0].pz()+dro_jets[1].pz(),2) );
                  
          jjMassDRO = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hDRO_MassJJ->Fill(jjMassDRO);
          hDRO_MassDiff->Fill((jjMassDRO-jjMassMCT)/jjMassMCT);
	  hDRO_ScatterEne->Fill(e1 - std::max(mct_jets[0].E(), mct_jets[1].E()), e2-std::min(mct_jets[0].E(), mct_jets[1].E()) );
	  if (debugMode) std::cout << "DRO: E_j1 + E_j2 = " << e1+e2 << " :: p_j1 + p_j2 = " << p1p2Sum << " :: jjMass = " << jjMassDRO << " GeV" << std::endl;
      }
                              
  }
  
  
  
  
  //plotting
  
  TCanvas * cMassJJ = new TCanvas ("cMassJJ", "cMassJJ", 600, 500);
  cMassJJ->cd();
    
  hMCT_MassJJ->Draw();
  hMCT_MassJJ->SetStats(0);
  hMCT_MassJJ->GetXaxis()->SetTitle("M_{jj} [GeV]");
  hMCT_MassJJ->GetYaxis()->SetTitle("Counts");
  hMCT_MassJJ->GetXaxis()->SetRangeUser(0, 140);
  hMCT_MassJJ->SetLineColor(kBlack);
  
  hRAW_MassJJ->Draw("same");
  hRAW_MassJJ->SetLineColor(kRed+1);
  
  hDRO_MassJJ->Draw("same");
  hDRO_MassJJ->SetLineColor(kGreen+1);
  
  
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");
  leg->AddEntry(hMCT_MassJJ, "MC truth", "lpf");
  leg->AddEntry(hRAW_MassJJ, "Raw calo jet", "lpf");
  leg->AddEntry(hDRO_MassJJ, "DRO calo jet", "lpf");

  leg->Draw();
  
  
    
  TCanvas * cMassJJ_Diff = new TCanvas ("cMassJJ_Diff", "cMassJJ_Diff", 600, 500);
  cMassJJ_Diff->cd();
    
  hDRO_MassDiff->Draw();
  hDRO_MassDiff->SetStats(0);
  hDRO_MassDiff->GetXaxis()->SetTitle("M_{jj}^{reco} - M_{jj}^{truth} / M_{jj}^{reco}");
  hDRO_MassDiff->GetYaxis()->SetTitle("Counts");
//   hMCT_MassJJ->GetXaxis()->SetRangeUser(0, 140);
  hDRO_MassDiff->SetLineColor(kGreen+1);
  
  hRAW_MassDiff->Draw("same");
  hRAW_MassDiff->SetLineColor(kRed+1);
  
  
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");
  leg->AddEntry(hRAW_MassDiff, "Raw calo jet", "lpf");
  leg->AddEntry(hDRO_MassDiff, "DRO calo jet", "lpf");

  leg->Draw();
  
  if (SAVEPLOTS) cMassJJ_Diff->SaveAs("plots/cMassJJ_Diff.png");

  TCanvas * cScatterEnergy = new TCanvas ("cScatterEnergy", "cScatterEnergy", 1000, 500);
  cScatterEnergy->Divide(2,1);
  
  cScatterEnergy->cd(1);    
  hRAW_ScatterEne->Draw("COLZ");
  hRAW_ScatterEne->SetStats(0);
  hRAW_ScatterEne->GetXaxis()->SetTitle("E_{j,1}^{reco} - E_{j,1}^{truth}");
  hRAW_ScatterEne->GetYaxis()->SetTitle("E_{j,2}^{reco} - E_{j,2}^{truth}");
  
  cScatterEnergy->cd(2);
  hDRO_ScatterEne->Draw("COLZ");
  hDRO_ScatterEne->SetStats(0);
  hDRO_ScatterEne->GetXaxis()->SetTitle("E_{j,1}^{reco} - E_{j,1}^{truth}");
  hDRO_ScatterEne->GetYaxis()->SetTitle("E_{j,2}^{reco} - E_{j,2}^{truth}");
    
  
  if (SAVEPLOTS) cScatterEnergy->SaveAs("plots/cScatterEnergy.png");
    
  
  TFile * outputFile = new TFile (Form("output_jjMass_%s_10k.root",output_tag.c_str() ) , "RECREATE");
  outputFile->cd();
  hMCT_MassJJ->Write();
  hRAW_MassJJ->Write();
  hDRO_MassJJ->Write();
  hRAW_MassDiff->Write();
  hDRO_MassDiff->Write();
  hRAW_ScatterEne->Write();
  hDRO_ScatterEne->Write();
  outputFile->Write();
  outputFile->Close();
  
  theApp->Run();
}







