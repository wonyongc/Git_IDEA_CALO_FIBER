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
//   std::string output_tag = "wwlj";
  std::string output_tag = "zjj_scan_90";
  
  int NFILES = 100;
  //  double x_factor_hcal = 0.43;
  double x_factor_hcal = 0.445;
  double x_factor_ecal = 0.371;
  
  float maxDeltaRMatchEcal = 0.013;
  float maxDeltaRMatchHcal = 0.1;
  
  
  if (argc>1) output_tag = argv[1];   
  if (argc>2) NFILES = atoi(argv[2]);
  if (argc>3) x_factor_hcal = atof(argv[3]);   
  if (argc>4) x_factor_ecal = atof(argv[4]);   
  if (argc>5) maxDeltaRMatchEcal = atof(argv[5]);   
  if (argc>6) maxDeltaRMatchHcal = atof(argv[6]);   
  
  double thismass = 100;
  if (output_tag == "wwlj") thismass = 80;
  if (output_tag == "hzjnbn") thismass = 90;
  if (output_tag == "hznb") thismass = 125;
  if (output_tag == "zjj_scan_30")  thismass = 30;
  if (output_tag == "zjj_scan_50")  thismass = 50;
  if (output_tag == "zjj_scan_70")  thismass = 70;
  if (output_tag == "zjj_scan_90")  thismass = 90;
  if (output_tag == "zjj_scan_100") thismass = 100;
  if (output_tag == "zjj_scan_150") thismass = 150;
  if (output_tag == "zjj_scan_250") thismass = 250;
  
  
  std::cout << "processing sample of: " << output_tag.c_str() << std::endl;  
  std::cout << "using x_factor_hcal = " << x_factor_hcal << " and x_factor_ecal = " << x_factor_ecal << std::endl;

  double ecal_S_norm = 0.985;
  double ecal_C_norm = 7286;
  double LO  = 2000;
  double CLO = 160;
  double drh_S_norm  = 407;
  double drh_C_norm  = 103.2;
  
  float ene_EC_th  = 0.01;
  float ene_HC_th  = 0.005;
  
  int phiGran = 252; //vs 36 as default
  
// choose a jet definition
  double R = 2*M_PI;
  int nExpJets = 2;
//       JetDefinition jet_def(antikt_algorithm, R);
  JetDefinition jet_def(ee_genkt_algorithm, R, 1);
  std::cout << "Clustering with " << jet_def.description() << std::endl;
  
    // single truth jet
  JetDefinition jet_mc(ee_genkt_algorithm, 4*M_PI, 1);

  double etaAcceptance = 1.4;  //accept only  jet within this eta
//   double phiAcceptance = 2.0;  //accept only  jet within this phi --> neglect phi border reconstruction effects..
//   double deltaR_match = 1.0;
  
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
    }
//     fname_reco  = Form("../root_files/hep_outputs/output_SCEPCal_B0T_%s100k_job_%d.root", output_tag.c_str(), iFile);
//     fname_truth = Form("../../HepMC_Files/B0T/%s100k_job_%d_output_tuple.root", output_tag.c_str(), iFile);
    
    fname_reco  = Form("../root_files/hep_outputs/output_SCEPCal_B0T_HG_%s_job_%d.root", output_tag.c_str(), iFile);
    fname_truth = Form("../../HepMC_Files/B0T/%s_job_%d_output_tuple.root", output_tag.c_str(), iFile);

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
  
  float ecal_stoch = 0.025;
  float ecal_const = 0.01;
  float hcal_stoch = 0.30;
  float hcal_const = 0.023;
  
  TF1 * funcEcalRes = new TF1 ("funcEcalRes", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 0, 300);
  funcEcalRes->SetParameters(ecal_stoch, ecal_const);
  TF1 * funcHcalRes = new TF1 ("funcHcalRes", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 0, 300);
  funcHcalRes->SetParameters(hcal_stoch, hcal_const);
  
  TF1 * funcTrackerRes = new TF1 ("funcTrackerRes", "[0]/sqrt(x)", 0, 300);
  funcTrackerRes->SetParameter(0,0.003);
  
  //define histos
  
  int NBIN = 1200;
  int minMass = 0;
  int maxMass = 300;
  
  TH1F * hMCT_MassJJ        = new TH1F ("hMCT_MassJJ", "hMCT_MassJJ", NBIN, minMass, maxMass);
  TH1F * hRAW_MassJJ        = new TH1F ("hRAW_MassJJ", "hRAW_MassJJ", NBIN, minMass, maxMass);
  TH1F * hDRO_MassJJ        = new TH1F ("hDRO_MassJJ", "hDRO_MassJJ", NBIN, minMass, maxMass);
  TH1F * hMCTFastSim_MassJJ = new TH1F ("hMCTFastSim_MassJJ", "hMCTFastSim_MassJJ", NBIN, minMass, maxMass);
  TH1F * hPFA_MassJJ        = new TH1F ("hPFA_MassJJ", "hPFA_MassJJ", NBIN, minMass, maxMass);
  TH1F * hPFA_RAW_MassJJ    = new TH1F ("hPFA_RAW_MassJJ", "hPFA_RAW_MassJJ", NBIN, minMass, maxMass);
  
  TH1F * hRAW_MassDiff        = new TH1F ("hRAW_MassDiff", "hRAW_MassDiff", NBIN, -1, 1);
  TH1F * hDRO_MassDiff        = new TH1F ("hDRO_MassDiff", "hDRO_MassDiff", NBIN, -1, 1);
  TH1F * hMCTFastSim_MassDiff = new TH1F ("hMCTFastSim_MassDiff", "hMCTFastSim_MassDiff", NBIN, -1, 1);
  TH1F * hPFA_MassDiff        = new TH1F ("hPFA_MassDiff", "hPFA_MassDiff", NBIN, -1, 1);
  TH1F * hPFA_RAW_MassDiff    = new TH1F ("hPFA_RAW_MassDiff", "hPFA_RAW_MassDiff", NBIN, -1, 1);
  
  TH1F * hDRO_Jet1EneDiff = new TH1F ("hDRO_Jet1EneDiff", "hDRO_Jet1EneDiff", NBIN, -1, 1);
  TH1F * hDRO_Jet2EneDiff = new TH1F ("hDRO_Jet2EneDiff", "hDRO_Jet2EneDiff", NBIN, -1, 1);      
  TH1F * hDRO_JetEneDiff  = new TH1F ("hDRO_JetEneDiff", "hDRO_JetEneDiff", NBIN, -1, 1);
  
  TH1F * hFastSim_Jet1EneDiff = new TH1F ("hFastSim_Jet1EneDiff", "hFastSim_Jet1EneDiff", NBIN, -1, 1);
  TH1F * hFastSim_Jet2EneDiff = new TH1F ("hFastSim_Jet2EneDiff", "hFastSim_Jet2EneDiff", NBIN, -1, 1);      
  TH1F * hFastSim_JetEneDiff  = new TH1F ("hFastSim_JetEneDiff", "hFastSim_JetEneDiff", NBIN, -1, 1);
  
  TH1F * hPFA_RAW_Jet1EneDiff = new TH1F ("hPFA_RAW_Jet1EneDiff", "hPFA_RAW_Jet1EneDiff", NBIN, -1, 1);
  TH1F * hPFA_RAW_Jet2EneDiff = new TH1F ("hPFA_RAW_Jet2EneDiff", "hPFA_RAW_Jet2EneDiff", NBIN, -1, 1);      
  TH1F * hPFA_RAW_JetEneDiff  = new TH1F ("hPFA_RAW_JetEneDiff", "hPFA_RAW_JetEneDiff", NBIN, -1, 1);
  
  TH1F * hPFA_Jet1EneDiff = new TH1F ("hPFA_Jet1EneDiff", "hPFA_Jet1EneDiff", NBIN, -1, 1);
  TH1F * hPFA_Jet2EneDiff = new TH1F ("hPFA_Jet2EneDiff", "hPFA_Jet2EneDiff", NBIN, -1, 1);      
  TH1F * hPFA_JetEneDiff  = new TH1F ("hPFA_JetEneDiff", "hPFA_JetEneDiff", NBIN, -1, 1);
  
  TH1F * hMCT_Jet1Ene = new TH1F ("hMCT_Jet1Ene", "hMCT_Jet1Ene", NBIN, 0, maxMass);
  TH1F * hRAW_Jet1Ene = new TH1F ("hRAW_Jet1Ene", "hRAW_Jet1Ene", NBIN, 0, maxMass);
  TH1F * hDRO_Jet1Ene = new TH1F ("hDRO_Jet1Ene", "hDRO_Jet1Ene", NBIN, 0, maxMass);
  TH1F * hMCT_Jet2Ene = new TH1F ("hMCT_Jet2Ene", "hMCT_Jet2Ene", NBIN, 0, maxMass);
  TH1F * hRAW_Jet2Ene = new TH1F ("hRAW_Jet2Ene", "hRAW_Jet2Ene", NBIN, 0, maxMass);
  TH1F * hDRO_Jet2Ene = new TH1F ("hDRO_Jet2Ene", "hDRO_Jet2Ene", NBIN, 0, maxMass);  
  
  TH1F* hGammaEneMC      = new TH1F ("hGammaEneMC", "hGammaEneMC", NBIN/2, 0, 1);
  TH1F* hNeutrHadMC      = new TH1F ("hNeutrHadMC", "hNeutrHadMC", NBIN/2, 0, 1);
  TH1F* hNeutralsMC      = new TH1F ("hNeutralsMC", "hNeutralsMC", NBIN/2, 0, 1);
  TH1F* hNeutralResidual = new TH1F ("hNeutralResidual", "hNeutralResidual", NBIN, -3, 3);
  TH1F* hECALResidual = new TH1F ("hECALResidual", "hECALResidual", NBIN, -3, 3);
  TH1F* hHCALResidual = new TH1F ("hHCALResidual", "hHCALResidual", NBIN, -3, 3);
  TH1F* hNeutralResidualDRO = new TH1F ("hNeutralResidualDRO", "hNeutralResidualDRO", NBIN, -3, 3);


  TH2F * hRAW_ScatterEne = new TH2F ("hRAW_ScatterEne", "hRAW_ScatterEne", NBIN, -75, 75, NBIN, -75, 75);
  TH2F * hDRO_ScatterEne = new TH2F ("hDRO_ScatterEne", "hDRO_ScatterEne", NBIN, -75, 75, NBIN, -75, 75);
  TH2F * hScatterEneVis  = new TH2F ("hScatterEneVis", "hScatterEneVis", 300, 0, 300, 300, 0, 300);
  TH2F * hScatterEneVisEH = new TH2F ("hScatterEneVisEH", "hScatterEneVisEH", 300, 0, 300, 200, 0, 5);

  ///*******************************************///
  ///		 Run over events	        ///
  ///*******************************************///
  int maxEVENTS = 100000;
  int NEVENTS = TreeRun->GetEntries();
  
  std::cout << "NEVENTS = " << NEVENTS << std::endl;
  if (NEVENTS>maxEVENTS)  NEVENTS = maxEVENTS;
  std::cout << "... running on " << NEVENTS << " events" << std::endl;
  int countGoodEvents = 0;
  
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
      
    std::cout << "processing event: " << iEvt << "\r" << std::flush;
    //    std::cout << "processing event: " << iEvt << std::endl;
    std::vector<PseudoJet> allHitsForJet;
    std::vector<PseudoJet> allMCHitsForJet;
    std::vector<PseudoJet> allMCHitsForJetFastSim;
    std::vector<PseudoJet> allHitsForJetPFA;
      
    bool goodEvent = true;
    int nMuons = 0;
    int nNeutrinos = 0;
    double neutrinoEne = 0;
    double muonEne = 0;

    double mc_phi_muon = -999;
    double mc_theta_muon = -999;

      //filling truth
    TruthTree->GetEntry(iEvt);
//       std::cout << "n MC truth particles found: " << myTruthTV.mcs_n << std::endl;
//       std::cout << "\n*******************************************************\n" << std::endl;

    double gamma_ene = 0;
    double neutralhad_ene = 0;
    double neutrals_ene = 0;
    
    
    for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
    {
          int    pdgId = myTruthTV.mcs_pdgId->at(i);
          double ene   = myTruthTV.mcs_E->at(i);
          if (ene<0) continue;
          double phi   = myTruthTV.mcs_phi->at(i);
          double eta   = myTruthTV.mcs_eta->at(i);
          double pT    = myTruthTV.mcs_pt->at(i);
          int charge   = myTruthTV.mcs_charge->at(i);
          double theta = 2*atan(exp(-eta));
          theta = M_PI- theta;
          
//           std::cout << " phi particle = " << phi << std::endl;
                      
//           if (debugMode) std::cout << "pdgId[" << i << "] = " << pdgId << std::endl;
          
          if (   fabs(pdgId)!=12 && fabs(pdgId)!=14 && fabs(pdgId)!=16 && fabs(pdgId)!=13  // exclude neutrinos and muons
                 && fabs(pdgId)<10000 // exclude BSM
//                         && fabs(eta)<etaAcceptance    // exclude particles outside the calorimeter
          )      
          {
              
              double px, py;
              px = pT*cos(phi);
              py = pT*sin(phi);
              
              double pz = -pT*sinh(eta);                    
              
              PseudoJet this_MCT = PseudoJet(px, py, pz, ene);
              allMCHitsForJet.push_back(this_MCT);
              
              
              
//               float smeared_ene = funcTrackerRes->Eval(pT)*pT;
              float smeared_ene = ene;
              if (charge!=0)// || fabs(pdgId) == 130 || fabs(pdgId) == 2112 )
              {                  
                  PseudoJet this_charged_track = PseudoJet(px*smeared_ene/ene, py*smeared_ene/ene, pz*smeared_ene/ene, smeared_ene);
                  this_charged_track.set_user_index(flag_MCT);
                  allHitsForJetPFA.push_back(this_charged_track);
              }
                  
//               if (charge!=0)
//               {
//                   bool isolated = true;
//                   for (unsigned int j = 0; j< myTruthTV.mcs_E->size(); j++)
//                   {
//                       int    j_pdgId = myTruthTV.mcs_pdgId->at(j);
//                       if (fabs(j_pdgId)!= 130 && fabs(j_pdgId)!= 2112) continue;
//                       
//                       double j_phi   = myTruthTV.mcs_phi->at(j);
//                       double j_theta = M_PI- 2*atan(exp(-myTruthTV.mcs_eta->at(j)));
//                       double dd = sqrt(pow(j_phi-phi,2)+pow(j_theta-theta,2) );
//                       
//                       if (dd<maxDeltaRMatchHcal)
//                       {
//                           isolated = false;
//                           break;
//                       }
//                   }
//                   if (isolated)
//                   {
//                       PseudoJet this_charged_track = PseudoJet(px*smeared_ene/ene, py*smeared_ene/ene, pz*smeared_ene/ene, smeared_ene);
//                       this_charged_track.set_user_index(flag_MCT);
//                       allHitsForJetPFA.push_back(this_charged_track);
//                   }
//               }
              
              //EM showers
              if (fabs(pdgId) == 22) gamma_ene+=ene;
              
              
              if (fabs(pdgId) == 22  || fabs(pdgId) == 11)
                  smeared_ene = gRandom->Gaus(ene, funcEcalRes->Eval(ene)*ene);
              //neutral hadrons
              else if (fabs(pdgId) == 130 || fabs(pdgId) == 2112)
              {
                  neutralhad_ene+= ene;
                  smeared_ene = gRandom->Gaus(ene, funcHcalRes->Eval(ene)*ene);
              }
              //charged hadrons
              else if (fabs(pdgId) == 211 || fabs(pdgId) == 321 || fabs(pdgId)==2212)
                  smeared_ene = gRandom->Gaus(ene, funcHcalRes->Eval(ene)*ene);
//               else
//                   std::cout << "particle not smeared (" << pdgId << "): " << ene << std::endl;
              
  
              if ( smeared_ene > 0.01 
//                   && fabs(pdgId) != 130 
//                   && fabs(pdgId) != 2112
            )
              {
                  PseudoJet this_MCT_FS = PseudoJet(px*smeared_ene/ene, py*smeared_ene/ene, pz*smeared_ene/ene, smeared_ene);
                  allMCHitsForJetFastSim.push_back(this_MCT_FS);
              }
              
              
              PseudoJet this_MCT_ghost = PseudoJet(px*1.0e-20, py*1.0e-20, pz*1.0e-20, ene*1.0e-20);        
              this_MCT_ghost.set_user_index(flag_MCT);                        
              allHitsForJet.push_back(this_MCT_ghost);
                                          
          }
          
          else
          {
              if (fabs(pdgId)==13) 
	      {
                  nMuons++;
		  muonEne+= ene;
		  mc_phi_muon = phi;
		  mc_theta_muon = theta;
	      }
              if (fabs(pdgId)==12 || fabs(pdgId)==14 || fabs(pdgId)==16) nNeutrinos++;
              if (fabs(pdgId)==12 || fabs(pdgId)==14 || fabs(pdgId)==16) 
              {
                  neutrinoEne += ene;
		  //                  std::cout << "neutrino Ene = " <<  ene << std::endl;
              }
          }
    }
    
    hGammaEneMC->Fill(gamma_ene/thismass);
    hNeutrHadMC->Fill(neutralhad_ene/thismass);
    neutrals_ene = gamma_ene+neutralhad_ene;
    hNeutralsMC->Fill(neutrals_ene/thismass);
  
    if (output_tag == "wwlj" && (nMuons>1 || nNeutrinos >1 ))
    {
        goodEvent = false;
        if (debugMode) std::cout << Form(" skipping %s event with %d muons and %d neutrinos ", output_tag.c_str(), nMuons, nNeutrinos) << std::endl;
        continue;
    }
  
    if (output_tag == "hzjnbn" && (nMuons>0 || nNeutrinos >0  ))
    {
        goodEvent = false;
        if (debugMode) std::cout << Form(" skipping %s event with %d muons and %d neutrinos ", output_tag.c_str(), nMuons, nNeutrinos) << std::endl;
        continue;
    } 
                                        
    if (output_tag == "hznb" && (nMuons>0 || nNeutrinos >2  ))
    {
        goodEvent = false;
        if (debugMode) std::cout << Form(" skipping %s event with %d muons and %d neutrinos ", output_tag.c_str(), nMuons, nNeutrinos) << std::endl;      
        continue;
    } 
    if (output_tag.find("zjj_scan") != string::npos && (nMuons>0 || nNeutrinos >0  ))
    {
        goodEvent = false;
        if (debugMode) std::cout << Form(" skipping %s event with %d muons and %d neutrinos", output_tag.c_str(), nMuons, nNeutrinos) << std::endl;
        continue;
    } 
     
    //filling reco

    TreeRun->GetEntry(iEvt);

    //**************************************************************//
    //                           DR HCAL
    //**************************************************************//
    
    if (output_tag == "hznb"  && (myTV.leakage/1000. - neutrinoEne > 1))// || myTV.leakage/1000.-muonEne))
    {
      if (debugMode)        std::cout << "Leakage - E_neutrino = " << myTV.leakage/1000. << "  - " << neutrinoEne <<  " = " << myTV.leakage/1000.- neutrinoEne << " GeV :: E_mu = " << muonEne << std::endl;
        goodEvent = false;
        continue;
    }
    if ((output_tag.find("zjj_scan") != string::npos) && (myTV.leakage/1000. > 0.3)) 
    {
      if (debugMode)        std::cout << "Leakage = " << myTV.leakage/1000. << " = " << myTV.leakage/1000. << std::endl;
        goodEvent = false;
        continue;
    }
    if (output_tag == "hzjnbn" && (myTV.leakage/1000. > 1))
    {
        if (debugMode)        std::cout << "Leakage = " << myTV.leakage/1000. - neutrinoEne << " GeV " << std::endl;
        goodEvent = false;
        continue;
    }
    
//     double ene_HC_th   = 0.01;    
    double totS = 0;
    double totEneDRH = 0;
    double edepMuonCalo = 0;
    
    for (unsigned int i = 0; i<myTV.VectorL->size(); i++)
    {                                        
        TVector3 this_vec = myGeometry.GetTowerVec(i,'l', phiGran);
//         double this_phi   = this_vec.Phi();
//         double this_theta = this_vec.Theta();
        double this_ene   = myTV.VectorL->at(i)/1000.;      
        double this_scint = myTV.VectorSignalsL->at(i);                
        double this_cher  = myTV.VectorSignalsCherL->at(i);
        double S = this_scint/drh_S_norm;
        double C = this_cher/drh_C_norm;
	double tower_phi_seed = this_vec.Phi();
// 	double tower_eta_seed = this_vec.Eta();
        double tower_theta_seed = this_vec.Theta();
	double deltaR = sqrt(pow(tower_phi_seed-mc_phi_muon,2)+pow(tower_theta_seed-mc_theta_muon,2) );
        if (deltaR <0.1)
 	{
	    edepMuonCalo+=S;
	    //	    std::cout << "HCAL --> muon ene in cone R=" <<deltaR << " : " << this_ene << " :: totMuonCaloDep = "  << edepMuonCalo << std::endl;
	}
        
        if (S>ene_HC_th)
        {            
            PseudoJet this_JHS = PseudoJet(this_vec.X()*S, this_vec.Y()*S, this_vec.Z()*S, S);
            this_JHS.set_user_index(flag_JHS);
            allHitsForJet.push_back(this_JHS);
            
            PseudoJet this_JHC = PseudoJet(this_vec.X()*C, this_vec.Y()*C, this_vec.Z()*C, C);
            this_JHC.set_user_index(flag_JHC);
            allHitsForJet.push_back(this_JHC);
            
            bool matchedToCharged  = false;
            bool matchedToNeutral  = false;
            bool matchedToNeutrHad = false;
            bool matchedToGamma = false;
            for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
            {
                int i_charge   = myTruthTV.mcs_charge->at(i);

                int    pdgId = myTruthTV.mcs_pdgId->at(i);
//               double ene   = myTruthTV.mcs_E->at(i);
                double i_phi   = myTruthTV.mcs_phi->at(i);
                double i_eta   = myTruthTV.mcs_eta->at(i);                
                double i_theta = 2*atan(exp(-i_eta));
                i_theta = M_PI-i_theta;
                double dd = sqrt(pow(tower_phi_seed-i_phi,2)+pow(tower_theta_seed-i_theta,2) );
                
                if (dd<maxDeltaRMatchHcal)//matched to charge
                {                    
                    if (i_charge!=0) 
                    {
                         matchedToCharged = true;
                    }
                    else if (fabs(pdgId) == 130 || fabs(pdgId) == 2112 )
                    {
                        matchedToNeutrHad = true;
                        matchedToNeutral = true;
                    }
                    else if (fabs(pdgId) == 22)     
                    {
                        matchedToNeutral = true;
                        matchedToGamma   = true;
                        
                    }
//                        std::cout << "HCAL hit matched to charged MC particle " << pdgId << std::endl;
//                        break;
                }
            }
            if (!matchedToCharged )// && !matchedToNeutrHad)// && matchedToNeutral)
            {
//                 std::cout << "HCAL hit NOT matched to charged MC particle " << std::endl;
                allHitsForJetPFA.push_back(this_JHS);
                allHitsForJetPFA.push_back(this_JHC);
            }
            else if (matchedToCharged && matchedToNeutrHad)
            {
//                 allHitsForJetPFA.push_back(this_JHS);
//                 allHitsForJetPFA.push_back(this_JHC);
            }

        }
        totS+=S;        
	totEneDRH+=this_ene;
    }
    
    for (unsigned int i = 0; i<myTV.VectorR->size(); i++)
    {                                        
        TVector3 this_vec = myGeometry.GetTowerVec(i,'r', phiGran);
//         double this_phi   = this_vec.Phi();
//         double this_theta = this_vec.Theta();
        double this_ene   = myTV.VectorR->at(i)/1000.;
        double this_scint = myTV.VectorSignalsR->at(i);
        double this_cher  = myTV.VectorSignalsCherR->at(i);
        double S = this_scint/drh_S_norm;
        double C = this_cher/drh_C_norm;
	double tower_phi_seed = this_vec.Phi();
// 	double tower_eta_seed = this_vec.Eta();
        double tower_theta_seed = this_vec.Theta();
	double deltaR = sqrt(pow(tower_phi_seed-mc_phi_muon,2)+pow(tower_theta_seed-mc_theta_muon,2) );

        if (deltaR <0.1)
 	{
	    edepMuonCalo+=S;
	    //	    std::cout << "HCAL --> muon ene in cone R=" <<deltaR << " : " << this_ene << " :: totMuonCaloDep = "  << edepMuonCalo << std::endl;
	}
        
        if (S>ene_HC_th)
        {
            PseudoJet this_JHS = PseudoJet(this_vec.X()*S, this_vec.Y()*S, this_vec.Z()*S, S);
            this_JHS.set_user_index(flag_JHS);
            allHitsForJet.push_back(this_JHS);
            
            PseudoJet this_JHC = PseudoJet(this_vec.X()*C, this_vec.Y()*C, this_vec.Z()*C, C);
            this_JHC.set_user_index(flag_JHC);
            allHitsForJet.push_back(this_JHC);   
            
            bool matchedToCharged = false;            
            bool matchedToNeutral = false;  
            bool matchedToNeutrHad = false;
            bool matchedToGamma = false;
            for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
            {
                int i_charge   = myTruthTV.mcs_charge->at(i);
                int    pdgId = myTruthTV.mcs_pdgId->at(i);
//                 double ene   = myTruthTV.mcs_E->at(i);
                double i_phi   = myTruthTV.mcs_phi->at(i);
                double i_eta   = myTruthTV.mcs_eta->at(i);                
                double i_theta = 2*atan(exp(-i_eta));
                i_theta = M_PI-i_theta;
                double dd = sqrt(pow(tower_phi_seed-i_phi,2)+pow(tower_theta_seed-i_theta,2) );
                
                if (dd<maxDeltaRMatchHcal)//matched to charge
                {                    
                    if (i_charge!=0) 
                    {
                        matchedToCharged = true;
                    }
                    else if (fabs(pdgId) == 130 || fabs(pdgId) == 2112 )
                    {
                        matchedToNeutrHad = true;
                        matchedToNeutral = true;
                    }
                    else if (fabs(pdgId) == 22) 
                    {
                        matchedToNeutral = true;
                        matchedToGamma   = true;
                    }
                }
            }
            if (!matchedToCharged)// && !matchedToNeutrHad)// && matchedToNeutral)
            {
//                 std::cout << "HCAL hit NOT matched to charged MC particle " << std::endl;
                allHitsForJetPFA.push_back(this_JHS);
                allHitsForJetPFA.push_back(this_JHC);
            }
            else if (matchedToCharged && matchedToNeutrHad)
            {
//                 allHitsForJetPFA.push_back(this_JHS);
//                 allHitsForJetPFA.push_back(this_JHC);
            }
        }
        totS+=S;        
	totEneDRH+=this_ene;
        
    }
        
    //    std::cout << "Total S in HCAL: " << totS <<  " GeV " << std::endl;
    // double thismass = 100;
    // if (output_tag == "wwlj") thismass = 80;
    // if (output_tag == "hzjnbn") thismass = 90;
    // if (output_tag == "hzbn") thismass = 125;
    // if (output_tag == "zjj_scan_100") thismass = 100;
    // if (totS/thismass<0.82)
    // {
    //   std::cout << "Total S in HCAL: " << totS <<  " GeV --> totS/Mass =  " << totS/thismass << std::endl;
    //   goodEvent = false;
    //   continue;
    // }
      
      
    
    //**************************************************************//    
    //                             ECAL
    //**************************************************************//

      
//     double ene_EC_th = 0.01;    
    double totEcalEne = 0;
            
    for (long unsigned int i = 0; i<myTV.VecHit_CrystalID->size(); i++)
    {
     
        TVector3 this_vec =  myGeometry.GetCrystalVec(myTV.VecHit_CrystalID->at(i));
        double this_ene = (myTV.VecHit_ScepEneDepF->at(i)+myTV.VecHit_ScepEneDepR->at(i))/1000.;                    
	double tower_phi_seed = this_vec.Phi();
// 	double tower_eta_seed = this_vec.Eta();
        double tower_theta_seed = this_vec.Theta();
	double deltaR = sqrt(pow(tower_phi_seed-mc_phi_muon,2)+pow(tower_theta_seed-mc_theta_muon,2) );
        // if (this_ene>0.1 && deltaR < 0.1)
	// {
	//   std::cout << "tower_phi = "  << tower_phi_seed << " :: mc_phi_muon =  " << mc_phi_muon << 
	//                " :: tower_eta = "<< tower_eta_seed << " :: mc_eta_muon = " << mc_eta_muon << 
	//                " :: deltaR = " << deltaR << " :: this_ene = " << this_ene << std::endl;
	// }
 
        if (deltaR<0.1)
 	{
	    edepMuonCalo+=this_ene;
	    //	    std::cout << "ECAL --> muon ene in cone R=" <<deltaR << " : " << this_ene << " :: totMuonCaloDep = "  << edepMuonCalo << std::endl;
	}
        
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
            
            bool matchedToCharged = false;
            bool matchedToNeutral = false;
            bool matchedToNeutrHad = false;
            bool matchedToGamma = false;
            for (unsigned int i = 0; i< myTruthTV.mcs_E->size(); i++)
            {
                int i_charge   = myTruthTV.mcs_charge->at(i);
                int    pdgId = myTruthTV.mcs_pdgId->at(i);
//                 double ene   = myTruthTV.mcs_E->at(i);
                double i_phi   = myTruthTV.mcs_phi->at(i);
                double i_eta   = myTruthTV.mcs_eta->at(i);                
                double i_theta = 2*atan(exp(-i_eta));
                i_theta = M_PI-i_theta;
                double dd = sqrt(pow(tower_phi_seed-i_phi,2)+pow(tower_theta_seed-i_theta,2) );
                
                if (dd<maxDeltaRMatchEcal)//matched to charge
                {                    
                    if (i_charge!=0) 
                    {
                        matchedToCharged = true;
                    }
                    else if (fabs(pdgId) == 130 || fabs(pdgId) == 2112 )
                    {
                        matchedToNeutrHad = true;
                        matchedToNeutral = true;
                    }
                    else if (fabs(pdgId) == 22)    
                    {
                        matchedToNeutral = true;
                        matchedToGamma   = true;
                    }

                }
            }
            
            if (!matchedToCharged && matchedToGamma)
            {
//                 std::cout << "HCAL hit NOT matched to charged MC particle " << std::endl;
                allHitsForJetPFA.push_back(this_JES*1.13);
                allHitsForJetPFA.push_back(this_JEC*1.13);
            }
            if (!matchedToCharged && matchedToNeutrHad)
            {
//                 std::cout << "HCAL hit NOT matched to charged MC particle " << std::endl;
                allHitsForJetPFA.push_back(this_JES);
                allHitsForJetPFA.push_back(this_JEC);
            }
            else if (matchedToCharged && matchedToNeutrHad)
            {
//                 allHitsForJetPFA.push_back(this_JES);
//                 allHitsForJetPFA.push_back(this_JEC);
            }
        
        }
        
        totEcalEne+=this_ene;
    }

    if (output_tag == "wwlj" && (myTV.leakage/1000. + edepMuonCalo - neutrinoEne -muonEne > 1))
    {
      if (debugMode)        std::cout << "Leakage + muonCaloDep - E_neutrino - muonEne = " << myTV.leakage/1000. << " + " << edepMuonCalo << "  - " << neutrinoEne << " - " << muonEne <<  " = " << myTV.leakage/1000. + edepMuonCalo- neutrinoEne - muonEne << " GeV" << std::endl;
        goodEvent = false;
        continue;
    }

    
    hScatterEneVis->Fill(totS,totEneDRH);
//     std::cout << "totS = " << totS << " :: totEneDRH = " << totEneDRH << " :: totS/vis = " << totS/totEneDRH << std::endl;
    //    std::cout << "Total S in HCAL: " << totS <<  " GeV :: totEcalEne = " << totEcalEne << " GeV ::  expMass = " << thismass << " :: --> (totS+totEcalEne)/Mass =  " << (totS+totEcalEne)/thismass << std::endl;

    if ((totS+totEcalEne)/thismass<0.8)
    {

      //      std::cout << "Total S in HCAL: " << totS <<  " GeV :: totEcalEne = " << totEcalEne << " GeV ::  --> (totS+totEcalEne)/Mass =  " << (totS+totEcalEne)/thismass << std::endl;
      goodEvent = false;
      continue;
    }

      // run the clustering, extract the jets
      // Monte Carlo truth
      ClusterSequence csMCOnly(allMCHitsForJet, jet_def);
      
      //fast sim
      ClusterSequence csMCFastSim(allMCHitsForJetFastSim, jet_def);
      
      //reco
      ClusterSequence cs(allHitsForJet, jet_def);
      
      //pfa
      ClusterSequence csPFA(allHitsForJetPFA, jet_def);
      
      
//       std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      std::vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(nExpJets));
      std::vector<PseudoJet> mct_ghost_jets;
      std::vector<PseudoJet> raw_jets;
      std::vector<PseudoJet> dro_jets;
                  
      std::vector<PseudoJet> mct_jets     = sorted_by_pt(csMCOnly.exclusive_jets(nExpJets));
      std::vector<PseudoJet> fastSim_jets = sorted_by_pt(csMCFastSim.exclusive_jets(nExpJets));
      
      std::vector<PseudoJet> pfa_jets     = sorted_by_pt(csPFA.exclusive_jets(nExpJets));
      std::vector<PseudoJet> pfa_jets_raw;
      std::vector<PseudoJet> pfa_jets_dro;

      //reco jets
      for (unsigned i = 0; i < jets.size(); i++) 
      {
//           cout << "jet " << i << ": "<< jets[i].pt() << " " << jets[i].rap() << " " << jets[i].phi() << endl;
          
          std::vector<PseudoJet> constituents = jets[i].constituents();
          std::vector<PseudoJet> mct_constituents;
          
          double E_MCT = 0;
          double E_JHS = 0;
          double E_JHC = 0;
          double E_JES = 0;
          double E_JEC = 0;
          
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
                  mct_constituents.push_back(PseudoJet(constituents[j].px()*1.0e20,constituents[j].py()*1.0e20,constituents[j].pz()*1.0e20,constituents[j].E()*1.0e20));
                  E_MCT += constituents[j].E()*1.0e20;              
              }
          }
          
	  //          if (debugMode) std::cout << "E_JES = " << E_JES << " :: E_JEC = " << E_JEC << " :: E_JHS = " << E_JHS << " :: E_JHC = " << E_JHC << std::endl;
          if (mct_constituents.size()>0) 
          {
              ClusterSequence csMC(mct_constituents, jet_mc);
              std::vector<PseudoJet> this_mct_jet = sorted_by_pt(csMC.exclusive_jets(int(1) ) );
              if (this_mct_jet.size()>0) mct_ghost_jets.push_back(this_mct_jet[0]);
          }
              
//           PseudoJet this_mct_jet = jets[i]*(E_MCT)/jets[i].E();
//           mct_ghost_jets.push_back(this_mct_jet);
          
          PseudoJet this_raw_jet = jets[i]*(E_JES+E_JHS)/jets[i].E();
          raw_jets.push_back(this_raw_jet);          
          
          double E_JE   = (E_JES-x_factor_ecal*E_JEC )/(1-x_factor_ecal);
          double E_JH   = (E_JHS-x_factor_hcal*E_JHC )/(1-x_factor_hcal);
          double E_JTot = E_JE+E_JH;
          PseudoJet dro_corr_jet = jets[i]*E_JTot/jets[i].E();
          dro_jets.push_back(dro_corr_jet);
          
      }
      
      double neutralhad_ene_reco = 0;
      double neutralhad_ene_reco_ecal = 0;
      double neutralhad_ene_reco_hcal = 0;
      double neutralhad_ene_reco_dro = 0;
      
      //pfa jets
      for (unsigned i = 0; i < pfa_jets.size(); i++) 
      {
          std::vector<PseudoJet> constituents = pfa_jets[i].constituents();
          double E_MCT = 0;
          double E_JHS = 0;
          double E_JHC = 0;
          double E_JES = 0;
          double E_JEC = 0;
          
          for (unsigned j = 0; j < constituents.size(); j++) 
          {      
              if      (constituents[j].user_index() == flag_JES) E_JES += constituents[j].E();
              else if (constituents[j].user_index() == flag_JEC) E_JEC += constituents[j].E();
              else if (constituents[j].user_index() == flag_JHS) E_JHS += constituents[j].E();
              else if (constituents[j].user_index() == flag_JHC) E_JHC += constituents[j].E();
              else if (constituents[j].user_index() == flag_MCT) E_MCT += constituents[j].E();
          }          	            
          
          PseudoJet this_raw_jet = pfa_jets[i]*(E_JES+E_JHS+E_MCT)/pfa_jets[i].E();
          pfa_jets_raw.push_back(this_raw_jet);
          
          double E_JE   = (E_JES-x_factor_ecal*E_JEC )/(1-x_factor_ecal);
          double E_JH   = (E_JHS-x_factor_hcal*E_JHC )/(1-x_factor_hcal);
          double E_JTot = E_JE+E_JH+E_MCT;
          PseudoJet dro_corr_jet = pfa_jets[i]*E_JTot/pfa_jets[i].E();
          pfa_jets_dro.push_back(dro_corr_jet);
          
          neutralhad_ene_reco     += E_JES+E_JHS;
          neutralhad_ene_reco_dro += E_JE+E_JH;
          
          neutralhad_ene_reco_ecal += E_JES;
          neutralhad_ene_reco_hcal += E_JHS;
          
//           neutralhad_ene_reco_dro += E_JE;
//           neutralhad_ene_reco_dro += E_JH;
          
      }
      
      
      hECALResidual->Fill((neutralhad_ene_reco_ecal-gamma_ene)/gamma_ene);
      hHCALResidual->Fill((neutralhad_ene_reco_hcal-neutralhad_ene)/neutralhad_ene);
      
      hNeutralResidual->Fill((neutralhad_ene_reco-neutrals_ene)/neutrals_ene);
      hNeutralResidualDRO->Fill((neutralhad_ene_reco_dro-neutrals_ene)/neutrals_ene);
      

//       mct_jets = mct_ghost_jets;
      
      //MCT jets
      if(mct_jets.size()==2)
      {
          //reject jets not fully contained in the calorimeter
          if ( fabs(mct_jets[0].eta()) > etaAcceptance || fabs(mct_jets[1].eta()) > etaAcceptance 
//             || fabs(mct_jets[0].phi()) > phiAcceptance || fabs(mct_jets[1].phi()) > phiAcceptance
        )
              
	  {
	    goodEvent = false;
  	    if (debugMode) std::cout << "skipping event with jets eta1 = " << fabs(mct_jets[0].eta()) << " :: eta2 = " << fabs(mct_jets[1].eta()) <<  " :: phi1 = " << mct_jets[0].phi() << " :: phi2 = " << mct_jets[1].phi() << std::endl;
	    continue;
	  }
      }
      
      //Monte Carlo truth
      double jjMassMCT = 0;
      if (mct_jets.size()==2 && goodEvent)
      {
          double e1 = mct_jets[0].E();
          double e2 = mct_jets[1].E();
          double p1p2Sum = sqrt(pow(mct_jets[0].px()+mct_jets[1].px(),2) + pow(mct_jets[0].py()+mct_jets[1].py(),2) + pow(mct_jets[0].pz()+mct_jets[1].pz(),2) );
          
          jjMassMCT = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hMCT_MassJJ->Fill(jjMassMCT);
          hMCT_Jet1Ene->Fill(e1);
          hMCT_Jet2Ene->Fill(e2);
	  //	  if (debugMode) std::cout << "MCT: E_j1 + E_j2 = " << e1+e2 << " :: p_j1 + p_j2 = " << p1p2Sum << " :: jjMass = " << jjMassMCT << " GeV" << std::endl;
      }
      
      //Fast Sim
      double jjMassMCT_FS = 0;
      if (fastSim_jets.size()==2 && goodEvent)
      {
          double e1 = fastSim_jets[0].E();
          double e2 = fastSim_jets[1].E();
          double p1p2Sum = sqrt(pow(fastSim_jets[0].px()+fastSim_jets[1].px(),2) + pow(fastSim_jets[0].py()+fastSim_jets[1].py(),2) + pow(fastSim_jets[0].pz()+fastSim_jets[1].pz(),2) );
          
          jjMassMCT_FS = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hMCTFastSim_MassJJ  ->Fill(jjMassMCT_FS);
          hMCTFastSim_MassDiff->Fill((jjMassMCT_FS-jjMassMCT)/jjMassMCT);
          
          if (sqrt(pow(fastSim_jets[0].phi() - mct_jets[0].phi(),2) + pow(fastSim_jets[0].theta() - mct_jets[0].theta(),2)) < sqrt(pow(fastSim_jets[0].phi() - mct_jets[1].phi(),2) + pow(fastSim_jets[0].theta() - mct_jets[1].theta(),2))) 
          {                             
            hFastSim_Jet1EneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hFastSim_Jet2EneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
            
            hFastSim_JetEneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hFastSim_JetEneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
          }
          else 
          {
            hFastSim_Jet1EneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hFastSim_Jet2EneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
            
            hFastSim_JetEneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hFastSim_JetEneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
          }          
      }
      
      //PFA Sim raw
      double jjMassMCT_PFA_RAW = 0;
      if (pfa_jets_raw.size()==2 && goodEvent)
      {
          double e1 = pfa_jets_raw[0].E();
          double e2 = pfa_jets_raw[1].E();
          double p1p2Sum = sqrt(pow(pfa_jets_raw[0].px()+pfa_jets_raw[1].px(),2) + pow(pfa_jets_raw[0].py()+pfa_jets_raw[1].py(),2) + pow(pfa_jets_raw[0].pz()+pfa_jets_raw[1].pz(),2) );
          
          jjMassMCT_PFA_RAW = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hPFA_RAW_MassJJ  ->Fill(jjMassMCT_PFA_RAW);
          hPFA_RAW_MassDiff->Fill((jjMassMCT_PFA_RAW-jjMassMCT)/jjMassMCT);
          
          if (sqrt(pow(pfa_jets_raw[0].phi() - mct_jets[0].phi(),2) + pow(pfa_jets_raw[0].theta() - mct_jets[0].theta(),2)) < sqrt(pow(pfa_jets_raw[0].phi() - mct_jets[1].phi(),2) + pow(pfa_jets_raw[0].theta() - mct_jets[1].theta(),2))) 
          {                             
            hPFA_RAW_Jet1EneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hPFA_RAW_Jet2EneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
            
            hPFA_RAW_JetEneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hPFA_RAW_JetEneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
          }
          else 
          {
            hPFA_RAW_Jet1EneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hPFA_RAW_Jet2EneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
            
            hPFA_RAW_JetEneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hPFA_RAW_JetEneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
          } 
      }
      //PFA Sim dro
      double jjMassMCT_PFA = 0;
      if (pfa_jets_dro.size()==2 && goodEvent)
      {
          double e1 = pfa_jets_dro[0].E();
          double e2 = pfa_jets_dro[1].E();
          double p1p2Sum = sqrt(pow(pfa_jets_dro[0].px()+pfa_jets_dro[1].px(),2) + pow(pfa_jets_dro[0].py()+pfa_jets_dro[1].py(),2) + pow(pfa_jets_dro[0].pz()+pfa_jets_dro[1].pz(),2) );
          
          jjMassMCT_PFA = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hPFA_MassJJ  ->Fill(jjMassMCT_PFA);
          hPFA_MassDiff->Fill((jjMassMCT_PFA-jjMassMCT)/jjMassMCT);
          
          if (sqrt(pow(pfa_jets_dro[0].phi() - mct_jets[0].phi(),2) + pow(pfa_jets_dro[0].theta() - mct_jets[0].theta(),2)) < sqrt(pow(pfa_jets_dro[0].phi() - mct_jets[1].phi(),2) + pow(pfa_jets_dro[0].theta() - mct_jets[1].theta(),2))) 
          {                             
            hPFA_Jet1EneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hPFA_Jet2EneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
            
            hPFA_JetEneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hPFA_JetEneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
          }
          else 
          {
            hPFA_Jet1EneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hPFA_Jet2EneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
            
            hPFA_JetEneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hPFA_JetEneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
          }          
      }

      //Reco raw
      double jjMassRAW = 0;
      if (raw_jets.size()==2 && goodEvent)
      {
          
          double e1 = raw_jets[0].E();
          double e2 = raw_jets[1].E();
          double p1p2Sum = sqrt(pow(raw_jets[0].px()+raw_jets[1].px(),2) + pow(raw_jets[0].py()+raw_jets[1].py(),2) + pow(raw_jets[0].pz()+raw_jets[1].pz(),2) );
          
          jjMassRAW = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hRAW_MassJJ->Fill(jjMassRAW);
          hRAW_MassDiff->Fill((jjMassRAW-jjMassMCT)/jjMassMCT);
          hRAW_Jet1Ene->Fill(e1);
          hRAW_Jet2Ene->Fill(e2);
          
          
          
          if (sqrt(pow(raw_jets[0].phi() - mct_jets[0].phi(),2) + pow(raw_jets[0].theta() - mct_jets[0].theta(),2)) < sqrt(pow(raw_jets[0].phi() - mct_jets[1].phi(),2) + pow(raw_jets[0].theta() - mct_jets[1].theta(),2))) 
          {                 
            hRAW_ScatterEne->Fill(e1 - mct_jets[0].E(), e2-mct_jets[1].E() );
          }
          else
          {
            hRAW_ScatterEne->Fill(e1 - mct_jets[1].E(), e2-mct_jets[0].E() );
          }
	  //	  if (debugMode) std::cout << "RAW: E_j1 + E_j2 = " << e1+e2 << " :: p_j1 + p_j2 = " << p1p2Sum << " :: jjMass = " << jjMassRAW << " GeV" << std::endl;
      }
      
      //Reco DRO
      double jjMassDRO = 0;
      if (dro_jets.size()==2 && goodEvent)
      {
          double e1 = dro_jets[0].E();
          double e2 = dro_jets[1].E();
          double p1p2Sum = sqrt(pow(dro_jets[0].px()+dro_jets[1].px(),2) + pow(dro_jets[0].py()+dro_jets[1].py(),2) + pow(dro_jets[0].pz()+dro_jets[1].pz(),2) );
                  
          jjMassDRO = sqrt(pow(e1+e2,2) - pow(p1p2Sum,2) );
          hDRO_MassJJ->Fill(jjMassDRO);
          hDRO_MassDiff->Fill((jjMassDRO-jjMassMCT)/jjMassMCT);
          hDRO_Jet1Ene->Fill(e1);
          hDRO_Jet2Ene->Fill(e2);
          
          
          if (sqrt(pow(dro_jets[0].phi() - mct_jets[0].phi(),2) + pow(dro_jets[0].theta() - mct_jets[0].theta(),2)) < sqrt(pow(dro_jets[0].phi() - mct_jets[1].phi(),2) + pow(dro_jets[0].theta() - mct_jets[1].theta(),2))) 
          {                             
            hDRO_ScatterEne  -> Fill(e1 - mct_jets[0].E(), e2 - mct_jets[1].E() );
            hDRO_Jet1EneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hDRO_Jet2EneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
            
            hDRO_JetEneDiff -> Fill((e1-mct_jets[0].E())/mct_jets[0].E());
            hDRO_JetEneDiff -> Fill((e2-mct_jets[1].E())/mct_jets[1].E());
          }
          else 
          {
            hDRO_ScatterEne->Fill(e1 - mct_jets[1].E(), e2 - mct_jets[0].E() );
            hDRO_Jet1EneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hDRO_Jet2EneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
            
            hDRO_JetEneDiff -> Fill((e1-mct_jets[1].E())/mct_jets[1].E());
            hDRO_JetEneDiff -> Fill((e2-mct_jets[0].E())/mct_jets[0].E());
          }
          
          
          
	  
          
//           if (fabs(e1 -mct_jets[0].E()>10)
//               std::cout << "max dro jet ene = " << e1 << " :: min dro jet ene = " << e2 << " :: max mct jet ene = " << std::max(mct_jets[0].E(), mct_jets[1].E()) << " :: min mct jet ene = " << std::min(mct_jets[0].E(), mct_jets[1].E()) << std::endl;
//           float MC_phi_jet1, DRO_phi_jet1;  
//           std::cout  << "MC_phi_jet1 =  " << MC_phi_jet1 << ":: DRO_phi_jet1 =  " << DRO_phi_jet1 << std::endl;
	  //	  if (debugMode) std::cout << "DRO: E_j1 + E_j2 = " << e1+e2 << " :: p_j1 + p_j2 = " << p1p2Sum << " :: jjMass = " << jjMassDRO << " GeV" << std::endl;
	  hScatterEneVisEH->Fill(jjMassDRO, totEcalEne/totEneDRH);
      }
                     
      countGoodEvents++;
  }

  
  
  std::cout << "selection efficiency: " << double(countGoodEvents)/double(NEVENTS) << std::endl;
  
  //plotting
  

  
    
  TCanvas * cMassJJ_Diff = new TCanvas ("cMassJJ_Diff", "cMassJJ_Diff", 600, 500);
  cMassJJ_Diff->cd();
    
  hDRO_MassDiff->Draw();
  //  hDRO_MassDiff->SetStats(0);
  hDRO_MassDiff->GetXaxis()->SetTitle("M_{jj}^{reco} - M_{jj}^{truth} / M_{jj}^{reco}");
  hDRO_MassDiff->GetYaxis()->SetTitle("Counts");
//   hMCT_MassJJ->GetXaxis()->SetRangeUser(0, 140);
  hDRO_MassDiff->SetLineColor(kGreen+1);
  
  hRAW_MassDiff->Draw("same");
  hRAW_MassDiff->SetLineColor(kRed+1);
  
  hMCTFastSim_MassDiff->Draw("same");
  hMCTFastSim_MassDiff->SetLineColor(kBlue);
  
  
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");
  leg->AddEntry(hMCTFastSim_MassDiff, "Fast sim jet", "lpf");
  leg->AddEntry(hRAW_MassDiff, "Raw calo jet", "lpf");
  leg->AddEntry(hDRO_MassDiff, "DRO calo jet", "lpf");

  leg->Draw();
  
  if (SAVEPLOTS) cMassJJ_Diff->SaveAs("plots/cMassJJ_Diff.png");

  TCanvas * cScatterEnergy = new TCanvas ("cScatterEnergy", "cScatterEnergy", 1000, 500);
  cScatterEnergy->Divide(2,1);
  
  cScatterEnergy->cd(1);    
  hRAW_ScatterEne->Draw("COLZ");
  //  hRAW_ScatterEne->SetStats(0);
  hRAW_ScatterEne->GetXaxis()->SetTitle("E_{j,1}^{reco} - E_{j,1}^{truth}");
  hRAW_ScatterEne->GetYaxis()->SetTitle("E_{j,2}^{reco} - E_{j,2}^{truth}");
  
  cScatterEnergy->cd(2);
  hDRO_ScatterEne->Draw("COLZ");
  //  hDRO_ScatterEne->SetStats(0);
  hDRO_ScatterEne->GetXaxis()->SetTitle("E_{j,1}^{reco} - E_{j,1}^{truth}");
  hDRO_ScatterEne->GetYaxis()->SetTitle("E_{j,2}^{reco} - E_{j,2}^{truth}");
    
  
  if (SAVEPLOTS) cScatterEnergy->SaveAs("plots/cScatterEnergy.png");

  TCanvas * cScatterEnergyVisible = new TCanvas ("cScatterEnergyVisible", "cScatterEnergyVisible", 1000, 500);
  cScatterEnergyVisible->Divide(2,1);    
  cScatterEnergyVisible->cd(1);    
  hScatterEneVis->Draw("COLZ");
  hScatterEneVis->GetXaxis()->SetTitle("E_{visible}^{HCAL}");
  hScatterEneVis->GetYaxis()->SetTitle("S_{tot}^{HCAL}");

  cScatterEnergyVisible->cd(2);    
  hScatterEneVisEH->Draw("COLZ");
  hScatterEneVisEH->GetXaxis()->SetTitle("M_{jj}^{HCAL}");
  hScatterEneVisEH->GetYaxis()->SetTitle("S_{ECAL}/S_{HCAL}");
  
  
  TCanvas * cPFA_Checks = new TCanvas ("cPFA_Checks", "cPFA_Checks", 1000, 500);
  cPFA_Checks->Divide(2,1);
  cPFA_Checks->cd(1);
    
  hGammaEneMC->Draw();
  hGammaEneMC->SetStats(0);
  hGammaEneMC->GetXaxis()->SetTitle("E_{tot} [GeV]");
  hGammaEneMC->GetYaxis()->SetTitle("Counts");
  hGammaEneMC->SetLineColor(kGreen+1);
  
  hNeutrHadMC->Draw("same");
  hNeutrHadMC->SetLineColor(kRed+1);
  
  hNeutralsMC->Draw("same");
  hNeutralsMC->SetLineColor(kBlack);
    
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");  
  leg->AddEntry(hGammaEneMC, "#gamma", "lpf");
  leg->AddEntry(hNeutrHadMC, "K^{0,L}, neutrons", "lpf");
  leg->AddEntry(hNeutralsMC, "All neutrals", "lpf");
  leg->Draw();
  
  
  
  cPFA_Checks->cd(2);
//   hNeutralResidual->Rebin(4);  
  hNeutralResidual->SetLineColor(kBlack);
  hNeutralResidual->Draw();
  hNeutralResidual->SetStats(0);
  hNeutralResidual->GetXaxis()->SetTitle("(E_{reco, neutr} - E_{MC, neutr}) / E_{MC,neutr}");
  hNeutralResidual->GetYaxis()->SetTitle("Counts");
  
//   hNeutralResidualDRO->Rebin(4);
  hNeutralResidualDRO->SetLineColor(kBlue+1);
  hNeutralResidualDRO->Draw("same");
  
//   hECALResidual->Rebin(4);
  hECALResidual->SetLineColor(kGreen+1);
  hECALResidual->Draw("same");
//   
//   hHCALResidual->Rebin(4);
  hHCALResidual->SetLineColor(kRed+1);
  hHCALResidual->Draw("same");
  
  
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");  
  leg->AddEntry(hNeutralResidual, "raw", "lpf");
  leg->AddEntry(hNeutralResidualDRO, "dro", "lpf");
  leg->AddEntry(hECALResidual, "ecal_reco - mc_gamma", "lpf");
  leg->AddEntry(hHCALResidual, "hcal_reco - mc_nhadrons", "lpf");
  leg->Draw();
    
  
  TFile * outputFile = new TFile (Form("output_jjMass_%s_xh%.3f_xe%.3f_dre%.3f_drh%.3f.root",output_tag.c_str(), x_factor_hcal, x_factor_ecal, maxDeltaRMatchEcal, maxDeltaRMatchHcal ) , "RECREATE");
  outputFile->cd();
  hMCT_MassJJ->Write();
  hMCTFastSim_MassJJ->Write();
  hPFA_MassJJ->Write();  
  hPFA_RAW_MassJJ->Write();  
  hRAW_MassJJ->Write();
  hDRO_MassJJ->Write();
  
  hPFA_MassDiff->Write();
  hPFA_RAW_MassDiff->Write();
  hMCTFastSim_MassDiff->Write();
  hRAW_MassDiff->Write();
  hDRO_MassDiff->Write();
  
  hRAW_ScatterEne->Write();
  hDRO_ScatterEne->Write();
  
  hPFA_Jet1EneDiff->Write();
  hPFA_Jet2EneDiff->Write();
  hPFA_JetEneDiff->Write();
  
  hPFA_RAW_Jet1EneDiff->Write();
  hPFA_RAW_Jet2EneDiff->Write();
  hPFA_RAW_JetEneDiff->Write();
  
  hFastSim_Jet1EneDiff->Write();
  hFastSim_Jet2EneDiff->Write();
  hFastSim_JetEneDiff->Write();
  
  hDRO_Jet1EneDiff->Write();
  hDRO_Jet2EneDiff->Write();
  hDRO_JetEneDiff->Write();
  
  hMCT_Jet1Ene->Write();
  hRAW_Jet1Ene->Write();
  hDRO_Jet1Ene->Write();
  hMCT_Jet2Ene->Write();
  hRAW_Jet2Ene->Write();
  hDRO_Jet2Ene->Write();
  
  hScatterEneVis->Write();
  hScatterEneVisEH->Write();
    
  hGammaEneMC->Write();
  hNeutrHadMC->Write();
  hNeutralsMC->Write();
  
  hNeutralResidual->Write();
  hNeutralResidualDRO->Write();
  hECALResidual->Write();
  hHCALResidual->Write();
  
  outputFile->Write();
  outputFile->Close();
  
  
  hECALResidual->Rebin(4);
  hHCALResidual->Rebin(4);
  hNeutralResidual->Rebin(4);
  hNeutralResidualDRO->Rebin(4);
  
  hMCT_MassJJ->Rebin(4);
  hMCTFastSim_MassJJ->Rebin(4);
  hRAW_MassJJ->Rebin(4);
  hDRO_MassJJ->Rebin(4);
  hPFA_MassJJ->Rebin(4);
  hPFA_RAW_MassJJ->Rebin(4);
  
  std::cout << "************************************************" << std::endl;
  std::cout << " ecal_reco - mc_gamma --> " << hECALResidual->GetMean() << " +/- " << hECALResidual->GetRMS() << std::endl;
  std::cout << " hcal_reco - mc_nhadrons  --> " << hHCALResidual->GetMean() << " +/- " << hHCALResidual->GetRMS() << std::endl;
  std::cout << " neutral residual  --> " << hNeutralResidual->GetMean() << " +/- " << hNeutralResidual->GetRMS() << std::endl;
  std::cout << " neutral residual DRO --> " << hNeutralResidualDRO->GetMean() << " +/- " << hNeutralResidualDRO->GetRMS() << std::endl;
  


std::cout << "************************************************" << std::endl;
  TCanvas * cMassJJ = new TCanvas ("cMassJJ", "cMassJJ", 600, 500);
  cMassJJ->cd();
    
  hMCT_MassJJ->Draw();
  //  hMCT_MassJJ->SetStats(0);
  hMCT_MassJJ->GetXaxis()->SetTitle("M_{jj} [GeV]");
  hMCT_MassJJ->GetYaxis()->SetTitle("Counts");
  hMCT_MassJJ->GetXaxis()->SetRangeUser(50, 140);
  hMCT_MassJJ->GetYaxis()->SetRangeUser(1, 1000);
  hMCT_MassJJ->SetLineColor(kBlack);
  hMCT_MassJJ->SetStats(0);
  
  hRAW_MassJJ->Draw("same");
  hRAW_MassJJ->SetLineColor(kRed+1);
  
  TF1 * fitGaus = new TF1 ("fitGaus", "gaus", 0, 140);
  hMCTFastSim_MassJJ->Draw("same");
  hMCTFastSim_MassJJ->SetLineColor(kBlue);
  fitGaus->SetLineColor(kBlue);
  hMCTFastSim_MassJJ->Fit(fitGaus, "QR");
  std::cout << "fast sim mjj resolution = " << fitGaus->GetParameter(2) << " / " << fitGaus->GetParameter(1) << " = " << fitGaus->GetParameter(2)/fitGaus->GetParameter(1) <<std::endl;
  std::cout << "fast sim mjj RMS/mean = " << hMCTFastSim_MassJJ->GetMean() << " / " << hMCTFastSim_MassJJ->GetRMS() << " = " << hMCTFastSim_MassJJ->GetRMS()/hMCTFastSim_MassJJ->GetMean() <<std::endl;
  
  
  fitGaus->SetLineColor(kRed);
  hRAW_MassJJ->Fit(fitGaus, "QR");
  std::cout << "raw mjj resolution = " << fitGaus->GetParameter(2) << " / " << fitGaus->GetParameter(1) << " = " << fitGaus->GetParameter(2)/fitGaus->GetParameter(1) <<std::endl;
  std::cout << "raw mjj RMS/mean = " << hRAW_MassJJ->GetMean() << " / " << hRAW_MassJJ->GetRMS() << " = " << hRAW_MassJJ->GetRMS()/hRAW_MassJJ->GetMean() <<std::endl;
  
  
  hDRO_MassJJ->Draw("same");
  hDRO_MassJJ->SetLineColor(kGreen+1);
  hDRO_MassJJ->SetLineWidth(2);
  fitGaus->SetLineColor(kGreen);
  hDRO_MassJJ->Fit(fitGaus, "QR");
  std::cout << "dro mjj resolution = " << fitGaus->GetParameter(2) << " / " << fitGaus->GetParameter(1) << " = " << fitGaus->GetParameter(2)/fitGaus->GetParameter(1) <<std::endl;
  std::cout << "dro mjj RMS/mean = " << hDRO_MassJJ->GetMean() << " / " << hDRO_MassJJ->GetRMS() << " = " << hDRO_MassJJ->GetRMS()/hDRO_MassJJ->GetMean() <<std::endl;
  
  
//   hPFA_RAW_MassJJ->Draw("same");
//   hPFA_RAW_MassJJ->SetLineColor(kYellow+1);
//   hPFA_RAW_MassJJ->SetLineWidth(2);
// //   hPFA_RAW_MassJJ->SetFillColor(kCyan);
//   fitGaus->SetLineColor(kYellow+1);
//   hPFA_RAW_MassJJ->Fit(fitGaus, "QR");
//   std::cout << "PFA RAW mjj resolution = " << fitGaus->GetParameter(2) << " / " << fitGaus->GetParameter(1) << " = " << fitGaus->GetParameter(2)/fitGaus->GetParameter(1) <<std::endl;
//   std::cout << "PFA RAW mjj RMS/mean = " << hPFA_RAW_MassJJ->GetMean() << " / " << hPFA_RAW_MassJJ->GetRMS() << " = " << hPFA_RAW_MassJJ->GetRMS()/hPFA_RAW_MassJJ->GetMean() <<std::endl;
//   
//   
//   hPFA_MassJJ->Draw("same");
//   hPFA_MassJJ->SetLineColor(kViolet);
// //   hPFA_MassJJ->SetFillColor(kViolet);
//   hPFA_MassJJ->SetLineWidth(2);
//   fitGaus->SetLineColor(kViolet);
//   hPFA_MassJJ->Fit(fitGaus, "QR");
//   std::cout << "PFA DRO mjj resolution = " << fitGaus->GetParameter(2) << " / " << fitGaus->GetParameter(1) << " = " << fitGaus->GetParameter(2)/fitGaus->GetParameter(1) <<std::endl;
//   std::cout << "PFA DRO mjj RMS/mean = " << hPFA_MassJJ->GetMean() << " / " << hPFA_MassJJ->GetRMS() << " = " << hPFA_MassJJ->GetRMS()/hPFA_MassJJ->GetMean() <<std::endl;
//   
  
  leg = new TLegend(0.75,0.75,0.95,0.95,NULL,"brNDC");
  leg->AddEntry(hMCT_MassJJ, "MC truth", "lpf");
  leg->AddEntry(hMCTFastSim_MassJJ, "Fast Sim", "lpf");
  leg->AddEntry(hRAW_MassJJ, "Raw calo jet", "lpf");
  leg->AddEntry(hDRO_MassJJ, "DRO calo jet", "lpf");
  leg->AddEntry(hPFA_MassJJ, "Proto PFA", "lpf");
  leg->AddEntry(hPFA_RAW_MassJJ, "Proto PFA raw", "lpf");

  leg->Draw();
  gPad->SetLogy();
    
  
  theApp->Run();
  
  
  
}







