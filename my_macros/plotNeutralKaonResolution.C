// g++ -Wall -o plotNeutralKaonResolution plotNeutralKaonResolution.C  VectorSmallestInterval.cc myG4Tree.cc myG4Tree.hh `root-config --cflags --glibs`


#include "VectorSmallestInterval.hh"
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

#include "TString.h"

#include "TSpline.h"
#include "TObject.h"
#include "myG4Tree.hh"
#include "CrystalBall.h"


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
  gStyle->SetLegendTextSize(0.03);
  TLegend * leg;
  

  
  std::string particle_name = "#gamma";

  bool SAVEPLOTS = false;
  
  int ENE_BINS = 10;  
//   int ENE_BINS = 14;  
  float maxEneRange = 100;
  float minEneRange = 1;
  float vMinEneBin [ENE_BINS];
  float vMaxEneBin [ENE_BINS];
  std::vector<double> *vTotEnergy[ENE_BINS];
  std::vector<double> *vTotEnergy_C[ENE_BINS];
  
  std::vector<double> *vTotEnergyCorr[ENE_BINS];
  std::vector<double> *vTotEnergyCorr_C[ENE_BINS];
  
  std::vector<double> *vHCALOnly[ENE_BINS];
  std::vector<double> *vHCALCherOnly[ENE_BINS];
  std::vector<double> *vHCALOnlyCorr[ENE_BINS];
  
  std::vector<double> *vECALwShowerEnergy[ENE_BINS];
  std::vector<double> *vECALCountNoMips[ENE_BINS];
  std::vector<double> *vECALCountAll[ENE_BINS];
  std::vector<double> *vInputEnergy[ENE_BINS];
  
  float binBoundary[ENE_BINS] = {0,1.5, 2.5, 4.5, 5.5, 8.5, 10.5, 15, 30, 40, 60, 80, 149, 195};
  
  for (int iEne = 0; iEne<ENE_BINS; iEne++)
  {
      float stepEne = (maxEneRange-minEneRange)/ENE_BINS;
      float min_ene = minEneRange+stepEne*iEne;
      float max_ene = min_ene+stepEne;
      
      //for uniform energy bins
      vMinEneBin[iEne]=min_ene;
      vMaxEneBin[iEne]=max_ene;
      
      //to use energy bins defined above
//       vMinEneBin[iEne]=binBoundary[iEne];
//       if (iEne<ENE_BINS-1) vMaxEneBin[iEne]=binBoundary[iEne+1];
//       else                 vMaxEneBin[iEne]=210;
      
      vECALwShowerEnergy[iEne]     = new std::vector<double>();
      vECALCountNoMips[iEne]     = new std::vector<double>();
      vECALCountAll[iEne]     = new std::vector<double>();
      
      vTotEnergy[iEne]     = new std::vector<double>();
      vTotEnergy_C[iEne]     = new std::vector<double>();
      
      vTotEnergyCorr[iEne] = new std::vector<double>();
      vTotEnergyCorr_C[iEne] = new std::vector<double>();
      
      vHCALOnly[iEne]     = new std::vector<double>();
      vHCALCherOnly[iEne] = new std::vector<double>();
      vHCALOnlyCorr[iEne] = new std::vector<double>();
      vInputEnergy[iEne]     = new std::vector<double>();
  }
    
  float mipInPWO = 0.220;
    
  double ecal_S_norm = 0.985;
//   double ecal_C_norm = 5460;
  double ecal_C_norm = 7286;
  float LO  = 2000;
//   float LO  = 20000;
  float CLO = 160;
//   float CLO = 16000;
  float LCE_C = CLO/ecal_C_norm;
//   double ecal_CS_norm = 1.4;
  
//   double drh_S_norm  = 375;
//   double drh_C_norm  = 375;
  
//   double drh_S_norm  = 410;
//   double drh_C_norm  = 105.6;
  double drh_S_norm  = 407;
  double drh_C_norm  = 103.2;  
  
  double x_factor_hcal = 0.43;
  double x_factor_ecal = 0.371;
//   double drh_S_norm  = 410*0.9;
//   double drh_C_norm  = 105.6*0.9;
  

  
  //define histos
  int NBIN_ENE = 400;
  float minEne = 0;
  float maxEne = 2;
  
  
  double maxEneECAL  = 1.2;
  double maxEneDRH   = drh_S_norm*1.2;
  double maxCherDRH  = drh_C_norm*1.2;
    
  int NBINS = 150;        
  float maxCS_norm = 2;    
  
  
  TH1F * hECAL_Scint = new TH1F ("hECAL_Scint", "hECAL_Scint", NBIN_ENE, 0, maxEneECAL*LO);
  TH1F * hECAL_Cher  = new TH1F ("hECAL_Cher", "hECAL_Cher", NBIN_ENE, 0, maxEneECAL*CLO);
  
  TH1F * hDRH_Scint = new TH1F ("hDRH_Scint", "hDRH_Scint", NBIN_ENE, minEne, maxEneDRH);
  TH1F * hDRH_Cher  = new TH1F ("hDRH_Cher", "hDRH_Cher", NBIN_ENE, minEne, maxCherDRH);
  
  
  TH1F * hECAL_S = new TH1F ("hECAL_S", "hECAL_S", NBIN_ENE, 0, maxEne);
  TH1F * hECAL_C = new TH1F ("hECAL_C", "hECAL_C", NBIN_ENE, 0, maxEne);
  
  TH2F * hECAL_CE_vs_SE = new TH2F ("hECAL_CE_vs_SE", "hECAL_CE_vs_SE", NBIN_ENE/4, 0, maxEne, NBIN_ENE/4, 0, maxEne);
  TH1F * hECAL_XFactor  = new TH1F ("hECAL_XFactor", "hECAL_XFactor", 200, 0, 2);
  
  TH1F * hHCAL_S = new TH1F ("hHCAL_S", "hHCAL_S", NBIN_ENE, 0, maxEne);
  TH1F * hHCAL_C = new TH1F ("hHCAL_C", "hHCAL_C", NBIN_ENE, 0, maxEne);
  
  TH2F * hHCAL_CE_vs_SE = new TH2F ("hHCAL_CE_vs_SE", "hHCAL_CE_vs_SE", NBIN_ENE/4, 0, maxEne, NBIN_ENE/4, 0, maxEne);
  TH1F * hHCAL_XFactor  = new TH1F ("hHCAL_XFactor", "hHCAL_XFactor", 200, 0, 2);
  
  TH1F * hHCAL_S_corr = new TH1F ("hHCAL_S_corr", "hHCAL_S_corr", NBIN_ENE, 0, maxEne);
  TH1F * hHCAL_C_corr = new TH1F ("hHCAL_C_corr", "hHCAL_C_corr", NBIN_ENE, 0, maxEne);
  
  TH1F * hHCALOnly_S = new TH1F ("hHCALOnly_S", "hHCALOnly_S", NBIN_ENE, 0, maxEne);
  TH1F * hHCALOnly_C = new TH1F ("hHCALOnly_C", "hHCALOnly_C", NBIN_ENE, 0, maxEne);
  
  TH1F * hHCALOnly_S_corr = new TH1F ("hHCALOnly_S_corr", "hHCALOnly_S_corr", NBIN_ENE, 0, maxEne);
  TH1F * hHCALOnly_C_corr = new TH1F ("hHCALOnly_C_corr", "hHCALOnly_C_corr", NBIN_ENE, 0, maxEne);


  
  TH1F * hEneTotal_S = new TH1F ("hEneTotal_S", "hEneTotal_S", NBIN_ENE, 0, maxEne);
  TH1F * hEneTotal_C = new TH1F ("hEneTotal_C", "hEneTotal_C", NBIN_ENE, 0, maxEne);
  
  TH1F * hEneTotal_S_corr = new TH1F ("hEneTotal_S_corr", "hEneTotal_S_corr", NBIN_ENE, 0, maxEne);
  TH1F * hEneTotal_C_corr = new TH1F ("hEneTotal_C_corr", "hEneTotal_C_corr", NBIN_ENE, 0, maxEne);
  
  
  TH2F * hScatterDRO_ECAL_S = new TH2F ("hScatterDRO_ECAL_S", "hScatterDRO_ECAL_S",   NBINS, 0, maxCS_norm, NBINS, 0, maxCS_norm);
  TH2F * hScatterDRO_ECAL_C = new TH2F ("hScatterDRO_ECAL_C", "hScatterDRO_ECAL_C",   NBINS, 0, maxCS_norm, NBINS, 0, maxCS_norm);
  
  TH2F * hScatterDRO_HCAL_S = new TH2F ("hScatterDRO_HCAL_S", "hScatterDRO_HCAL_S",   NBINS, 0, maxCS_norm, NBINS, 0, 1.5);
  TH2F * hScatterDRO_HCAL_C = new TH2F ("hScatterDRO_HCAL_C", "hScatterDRO_HCAL_C",   NBINS, 0, maxCS_norm, NBINS, 0, 1.5);
  
//   TProfile2D * hScatterECAL_frac = new TProfile2D ("hScatterECAL_frac", "hScatterECAL_frac",   NBINS/6, 0, 100, NBINS/6, 0, 1.2);
  TH2D * hScatterECAL_frac = new TH2D ("hScatterECAL_frac", "hScatterECAL_frac",   NBINS/6, 0, 100, NBINS/6, 0, 1.2);
  

  
  
  //run over energy scan
  
//   TFile * RunFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_kaon0L_Iso+Uniform1-100_GeV.root","READ");
  
//   TFile * RunFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_pi-_Iso+Uniform1-100_GeV.root","READ");            
//   TTree* TreeRun = (TTree*) RunFile->Get("B4");    
      
    
  TChain * TreeRun = new TChain("B4", "B4");
//   TreeRun->Add("../root_files/prod/output_SCEPCal_fixedPos_kaon0L_Iso+Uniform1-100_GeV.root");
  
  TreeRun->Add("../root_files/prod/output_SCEPCal_fixedPos_pi-_Iso+Uniform1-100_GeV.root");
//   TreeRun->Add("../root_files/prod/output_SCEPCal_fixedEne_pi-_IsoSupplement.root");
  
    
  myG4TreeVars myTV;
  InitG4Tree (TreeRun, myTV);
  
  TreeRun->SetBranchStatus("*", 0);
  
  TreeRun->SetBranchStatus("PrimaryParticleEnergy", 1);  
  TreeRun->SetBranchStatus("PrimaryParticleMomentum", 1);  
  TreeRun->SetBranchStatus("SCEP_EnergyDepF", 1);
  TreeRun->SetBranchStatus("SCEP_EnergyDepR", 1);
  TreeRun->SetBranchStatus("SCEP_NCherProdF", 1);
  TreeRun->SetBranchStatus("SCEP_NCherProdR", 1);    
  
  TreeRun->SetBranchStatus("Energyem", 1);
  TreeRun->SetBranchStatus("leakage", 1);
  TreeRun->SetBranchStatus("EnergyScin", 1);
  TreeRun->SetBranchStatus("EnergyCher", 1);
  TreeRun->SetBranchStatus("NofCherenkovDetected", 1);
  TreeRun->SetBranchStatus("VectorSignalsL", 1);
  TreeRun->SetBranchStatus("VectorSignalsR", 1);
  TreeRun->SetBranchStatus("VectorSignalsCherL", 1);
  TreeRun->SetBranchStatus("VectorSignalsCherR", 1);
  
  
  
      
            
  
  ///*******************************************///
  ///		 Run over events	    ///
  ///*******************************************///
    
  
  int NEVENTS = TreeRun->GetEntries();
  std::cout << "NEVENTS = " << NEVENTS << std::endl;
  
//   NEVENTS = 30000;
      
  //first loop to calculate hcal dro correction
  std::cout << "first loop to calculate hcal dro correction..." << std::endl;
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
      
      TreeRun->GetEntry(iEvt);                  
      //           std::cout << "iEvent = " << iEvt << std::endl;
      
      double this_ene = myTV.PrimaryParticleEnergy/1000; // in GeV
      double px  = myTV.PrimaryParticleMomentum->at(0);
      double py  = myTV.PrimaryParticleMomentum->at(1);
      double pz  = myTV.PrimaryParticleMomentum->at(2);
      double P   = sqrt(px*px+py*py+pz*pz);
      px/= P;
      py/= P;
      pz/= P;
          
      
//       double phi   = atan(py/px);
      double eta   = atanh(pz);
      double theta = 2*atan(exp(-eta));
//       std::cout << " theta = " << theta << std::endl;
      if (theta <0.5 || theta>2.5) continue;
      
//       std::cout << "leakage  =" << myTV.leakage/1000 << " GeV" << std::endl;
      if      (this_ene < 12 && myTV.leakage/1000. > 0.5) continue;
      else if (this_ene < 50 && myTV.leakage/1000. > 1.0) continue;
      else if (this_ene > 50 && myTV.leakage/1000. > 3.0) continue;

      
//       if (this_ene>50) std::cout << "this_ene  = " << this_ene << std::endl;
      
      double eneF   = myTV.SCEP_EnergyDepF/1000.;      
      double eneR   = myTV.SCEP_EnergyDepR/1000.;      
//       double cherF  = myTV.SCEP_NCherProdF;      
      double cherR  = myTV.SCEP_NCherProdR;            
      
      //ecal
      double ecal_ene_tot  = (eneR+eneF) / this_ene;                          
//       double ecal_cher_tot = (cherR) / this_ene;
      
      double ecal_S_F = gRandom->Poisson(eneR*LO)/this_ene/ecal_S_norm/LO;
      double ecal_S_R = gRandom->Poisson(eneF*LO)/this_ene/ecal_S_norm/LO;
      double ecal_S   = ecal_S_F+ecal_S_R;
      double ecal_C   = gRandom->Poisson(cherR*LCE_C)/this_ene/CLO;
//       double ecal_C_S = ecal_C / ecal_S;   

      //hcal
      float sumScintWBirks = 0;
      float sumCherTowers  = 0;
      for (long unsigned int iT = 0; iT < myTV.VectorSignalsL->size(); iT++)
      {
          sumScintWBirks+= myTV.VectorSignalsL->at(iT) + myTV.VectorSignalsR->at(iT);
          sumCherTowers += myTV.VectorSignalsCherL->at(iT) + myTV.VectorSignalsCherR->at(iT);
      }
      double drh_scint = sumScintWBirks/this_ene;
//       double drh_scint = myTV.EnergyScin/this_ene;
      double drh_cher  = sumCherTowers/this_ene;
//       double drh_cher  = myTV.NofCherenkovDetected/this_ene;
      double drh_S     = drh_scint/drh_S_norm;
      double drh_C     = drh_cher/drh_C_norm;
      double hcal_C_S  = drh_C / drh_S;
      
                          
      if (drh_S<0.01 || drh_C<0.01 || ecal_S<0.001 ) continue;
      
      hECAL_Scint->Fill(ecal_S*LO);
      hECAL_Cher->Fill(ecal_C*CLO);
      
      hECAL_S->Fill(ecal_S);
      hECAL_C->Fill(ecal_C);
      
//       std::cout << "drh_S = " << drh_S << " :: drh_C = " << drh_C << " :: hcal_C_S = " << hcal_C_S << std::endl;
      
      hDRH_Scint->Fill(drh_scint);
      hDRH_Cher->Fill(drh_cher);
      
      hHCAL_S->Fill(drh_S);
      hHCAL_C->Fill(drh_C);
//       std::cout << " ecal_ene_tot*this_ene = " << ecal_ene_tot*this_ene << std::endl;
      if (ecal_ene_tot*this_ene<mipInPWO*2 && ecal_ene_tot <0.03)
      {
        hScatterDRO_HCAL_S -> Fill(hcal_C_S, drh_S);
        hScatterDRO_HCAL_C -> Fill(hcal_C_S, drh_C);
        
//         if (this_ene>40) 
        {
            hHCAL_CE_vs_SE -> Fill(drh_S, drh_C);
            hHCAL_XFactor->Fill((drh_S-1)/(drh_C-1));
            
        }
        
        hHCALOnly_S->Fill(drh_S);
        hHCALOnly_C->Fill(drh_C);
      }                                    
  }
  


  //ECAL segment
  
  TCanvas * cECAL_Scint = new TCanvas ("cECAL_Scint", "cECAL_Scint", 600, 500);
  cECAL_Scint->cd();
  hECAL_Scint->Draw();
  hECAL_Scint->GetXaxis()->SetTitle("S/E [phe/GeV]");
  gPad->SetLogy();
  
  TCanvas * cECAL_Cher = new TCanvas ("cECAL_Cher", "cECAL_Cher", 600, 500);
  cECAL_Cher->cd();
  hECAL_Cher->Draw();
  hECAL_Cher->GetXaxis()->SetTitle("C/E [phe/GeV]");
  gPad->SetLogy();
  
  TCanvas * cECAL_S = new TCanvas ("cECAL_S", "cECAL_S", 600, 500);
  cECAL_S->cd();
  hECAL_S->Draw();
  hECAL_S->GetXaxis()->SetTitle("S_{norm}/E");
  gPad->SetLogy();
  
  TCanvas * cECAL_C = new TCanvas ("cECAL_C", "cECAL_C", 600, 500);
  cECAL_C->cd();
  hECAL_C->Draw();
  hECAL_C->GetXaxis()->SetTitle("C_{norm}/E");
  gPad->SetLogy();
  
  
  //HCAL segment
  
  TCanvas * cDRH_Scint = new TCanvas ("cDRH_Scint", "cDRH_Scint", 600, 500);
  cDRH_Scint->cd();
  hDRH_Scint->Draw();
  gPad->SetLogy();
  
  TCanvas * cDRH_Cher = new TCanvas ("cDRH_Cher", "cDRH_Cher", 600, 500);
  cDRH_Cher->cd();
  hDRH_Cher->Draw();
  gPad->SetLogy();
  
  TCanvas * cHCAL_S = new TCanvas ("cHCAL_S", "cHCAL_S", 600, 500);
  cHCAL_S->cd();
  hHCAL_S->SetStats(0);
  hHCAL_S->GetXaxis()->SetTitle("S_{HCAL}");
  hHCAL_S->SetLineColor(kBlack);
  hHCAL_S->Draw();
  hHCALOnly_S->SetLineColor(kRed);
  hHCALOnly_S->Draw("same");
  TLegend* legHCAL_s = new TLegend(0.15,0.15,0.75,0.35,NULL,"brNDC");
  legHCAL_s->AddEntry(hHCAL_S, "All events");
  legHCAL_s->AddEntry(hHCALOnly_S, "no shower in ECAL");  
  gPad->SetLogy();
  
  
  TCanvas * cHCAL_C = new TCanvas ("cHCAL_C", "cHCAL_C", 600, 500);
  cHCAL_C->cd();
  hHCAL_C->SetLineColor(kBlack);
  hHCAL_C->Draw();
  hHCALOnly_C->SetLineColor(kRed);
  hHCALOnly_C->Draw("same");
  gPad->SetLogy();
  
  
  TCanvas * cScatterDRO_HCAL_S = new TCanvas ("cScatterDRO_HCAL_S", "cScatterDRO_HCAL_S", 600, 500);
  cScatterDRO_HCAL_S->cd();
  hScatterDRO_HCAL_S->SetStats(0);
  hScatterDRO_HCAL_S->Draw("COLZ");
  hScatterDRO_HCAL_S->GetXaxis()->SetTitle("C/S_{DRH}");
  hScatterDRO_HCAL_S->GetYaxis()->SetTitle("S_{DRH}/E [ph/GeV?]");
  
  TProfile* pDRO_corr_HCAL_S = hScatterDRO_HCAL_S->ProfileX();
  pDRO_corr_HCAL_S->SetStats(0);
  pDRO_corr_HCAL_S->SetLineColor(kRed);
  pDRO_corr_HCAL_S->SetMarkerColor(kRed);
  pDRO_corr_HCAL_S->SetMarkerStyle(20);
  pDRO_corr_HCAL_S->Draw("same");
  hScatterDRO_HCAL_S->GetXaxis()->SetRangeUser(0., 1.2);
  hScatterDRO_HCAL_S->GetYaxis()->SetRangeUser(0., 1.2);
  
  TF1 * funcHCAL_DRO_S = new TF1 ("funcHCAL_DRO_S", "pol2", 0.5, 1.);
  //[0]*log([1]*x)+[2]
  funcHCAL_DRO_S->SetParameters(1,1,1);
  
  funcHCAL_DRO_S->SetLineColor(kRed);
  for (int i = 0; i < 10; i++) pDRO_corr_HCAL_S->Fit(funcHCAL_DRO_S, "QR");
  funcHCAL_DRO_S->Draw("same");                
  gPad->SetGrid();
  
  
  TCanvas * cScatterDRO_HCAL_C = new TCanvas ("cScatterDRO_HCAL_C", "cScatterDRO_HCAL_C", 600, 500);
  cScatterDRO_HCAL_C->cd();
  hScatterDRO_HCAL_C->SetStats(0);
  hScatterDRO_HCAL_C->Draw("COLZ");
  hScatterDRO_HCAL_C->GetXaxis()->SetTitle("C/S_{DRH}");
  hScatterDRO_HCAL_C->GetYaxis()->SetTitle("C_{DRH}/E [ph/GeV?]");
  
  TProfile* pDRO_corr_HCAL_C = hScatterDRO_HCAL_C->ProfileX();
  pDRO_corr_HCAL_C->SetStats(0);
  pDRO_corr_HCAL_C->SetLineColor(kRed);
  pDRO_corr_HCAL_C->SetMarkerColor(kRed);
  pDRO_corr_HCAL_C->SetMarkerStyle(20);
  pDRO_corr_HCAL_C->Draw("same");
        
  hScatterDRO_HCAL_C->GetXaxis()->SetRangeUser(0., 1.2);
  hScatterDRO_HCAL_C->GetYaxis()->SetRangeUser(0., 1.2);
  
  TF1 * funcHCAL_DRO_C = new TF1 ("funcHCAL_DRO_C", "pol2", 0.2, 1.);
  //[0]*log([1]*x)+[2]
  funcHCAL_DRO_C->SetParameters(1,1,1);
  
  funcHCAL_DRO_C->SetLineColor(kRed);
  for (int i = 0; i < 10; i++) pDRO_corr_HCAL_C->Fit(funcHCAL_DRO_C, "QR");
  funcHCAL_DRO_C->Draw("same");  
  gPad->SetGrid();
  
  
  
//   TGraphErrors* gMPV_CE_vs_SE = new TGraphErrors();
//   TCanvas * cSlices = new TCanvas ("cSlices", "cSlices", 600, 500);
//   cSlices->cd();
//   TH1F* hSlice [NBIN_ENE/4];
//   int ipoint = 0;
//   for (int islice = 0; islice< NBIN_ENE/4; islice++)
//   {
//       hSlice[islice] = new TH1F (Form("hSlice%d", islice), Form("hSlice%d", islice), NBIN_ENE/4, 0, maxEne);
//       for (int ibin = 0; ibin< NBIN_ENE/4; ibin++)
//       {
//           hSlice[islice]->SetBinContent(ibin, hHCAL_CE_vs_SE->GetBinContent(islice,ibin));
// //           std::cout << " islice = " << islice << " :: bincontent = " << hHCAL_CE_vs_SE->GetBinContent(islice,ibin) << std::endl;
//       }
//       if (islice == 0) 
//       {
//           hSlice[0]->Draw("h");
//           hSlice[0]->GetXaxis()->SetRangeUser(0, 1.2);
//           hSlice[0]->GetYaxis()->SetRangeUser(0, 100);
//       }
//       else 
//       {
// //           hSlice[islice]->SetLineColor(kRainBow+islice);
//           hSlice[islice]->Draw("same h");
//       }
//       
//       float sliceMax = hSlice[islice]->GetBinCenter(hSlice[islice]->GetMaximumBin());
// //       TF1 * fitSlice = new TF1 ("fitSlice", "gaus", sliceMax*0.7, sliceMax*1.3);
// 
//       if (sliceMax>0.1) 
//       {
//       
//           TF1 * fitSlice = new TF1 ("fitSlice",crystalBallDouble, sliceMax*0.8, sliceMax*1.2, 5);
//           fitSlice->SetParameters(1, sliceMax, 0.2, 0.2, 0.2);
//       
//           hSlice[islice]->Fit(fitSlice, "QR");
//           float binCenter = maxEne/NBIN_ENE*4*(islice+0.5);
//           float mpv = fitSlice->GetParameter(1);
//           std::cout << " islice = " << islice << " sliceMax = " << sliceMax <<  " :: binCenter = " << binCenter << " :: mpv = " << mpv << std::endl; 
//           gMPV_CE_vs_SE->SetPoint(ipoint, binCenter , mpv);
//           ipoint++;
//       }
//   }
//   
  
  
  TCanvas * cScatterDRO_HCAL_CE_SE = new TCanvas ("cScatterDRO_HCAL_CE_SE", "cScatterDRO_HCAL_CE_SE", 600, 600);
  cScatterDRO_HCAL_CE_SE->cd();
  hHCAL_CE_vs_SE->SetTitle("");
  hHCAL_CE_vs_SE->SetStats(0);
  hHCAL_CE_vs_SE->Draw("COLZ");  
  hHCAL_CE_vs_SE->GetXaxis()->SetTitle("S_{DRH}/E");
  hHCAL_CE_vs_SE->GetYaxis()->SetTitle("C_{DRH}/E ");
  
  
  TProfile* pHCAL_CE_vs_SE = hHCAL_CE_vs_SE->ProfileX();
  pHCAL_CE_vs_SE->SetStats(0);
  pHCAL_CE_vs_SE->SetLineColor(kRed);
  pHCAL_CE_vs_SE->SetMarkerColor(kRed);
  pHCAL_CE_vs_SE->SetMarkerStyle(20);
  pHCAL_CE_vs_SE->SetMarkerSize(0.8);
  pHCAL_CE_vs_SE->Draw("same");
  hHCAL_CE_vs_SE->GetXaxis()->SetRangeUser(0., 1.2);
  hHCAL_CE_vs_SE->GetYaxis()->SetRangeUser(0., 1.2);
  
  TF1 * funcHCAL_CE_vs_SE = new TF1 ("funcHCAL_CE_vs_SE", "1-[0]+[0]*x", 0.5, 1.);
  //[0]*log([1]*x)+[2]
  funcHCAL_CE_vs_SE->SetParameters(1,1,1);
  
  funcHCAL_CE_vs_SE->SetLineColor(kBlack);
  funcHCAL_CE_vs_SE->SetLineWidth(3);
  for (int i = 0; i < 10; i++) pHCAL_CE_vs_SE->Fit(funcHCAL_CE_vs_SE, "QR");
  funcHCAL_CE_vs_SE->Draw("same");        
  
  TF1 * CequalS = new TF1("CequalS", "x", 0, 1.2);
  CequalS->SetLineStyle(7);
  CequalS->SetLineWidth(3);
  CequalS->SetLineColor(kViolet);
  CequalS->Draw("same");
  gPad->SetGrid();
  
//   gMPV_CE_vs_SE->Draw("same PE");
//   gMPV_CE_vs_SE->SetMarkerStyle(20);
//   gMPV_CE_vs_SE->SetMarkerColor(kViolet);
  if (SAVEPLOTS) cScatterDRO_HCAL_CE_SE->SaveAs("plots/cScatterDRO_HCAL_CE_SE.pdf");
  
  
  
  TCanvas * cHCAL_XFactor = new TCanvas ("cHCAL_XFactor", "cHCAL_XFactor", 600, 600);
  cHCAL_XFactor->cd();
  hHCAL_XFactor->SetStats(0);
  hHCAL_XFactor->SetTitle("");
  hHCAL_XFactor->GetXaxis()->SetTitle("#Chi_{HCAL}");
  hHCAL_XFactor->GetYaxis()->SetTitle("Counts");
  hHCAL_XFactor->SetLineColor(kBlack);
  hHCAL_XFactor->Draw();
//   gPad->SetLogy();
  if (SAVEPLOTS) cHCAL_XFactor->SaveAs("plots/cHCAL_XFactor.pdf");
  
  
  
  //second loop to calculate ecal dro correction
  std::cout << "second loop to calculate ecal dro correction..." << std::endl;
  
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
      
      TreeRun->GetEntry(iEvt);                        
      
      
      double this_ene = myTV.PrimaryParticleEnergy/1000; // in GeV
      double px  = myTV.PrimaryParticleMomentum->at(0);
      double py  = myTV.PrimaryParticleMomentum->at(1);
      double pz  = myTV.PrimaryParticleMomentum->at(2);
      double P   = sqrt(px*px+py*py+pz*pz);
      px/= P;
      py/= P;
      pz/= P;
          
      double eta   = atanh(pz);
      double theta = 2*atan(exp(-eta));//       double phi   = atan(py/px);
//       double eta   = atanh(pz);
      if (theta <0.5 || theta>2.5) continue;
      if      (this_ene < 12 && myTV.leakage/1000. > 0.5) continue;
      else if (this_ene < 50 && myTV.leakage/1000. > 1.0) continue;
      else if (this_ene > 50 && myTV.leakage/1000. > 3.0) continue;
//       if (this_ene<5) continue;
      
      double eneF   = myTV.SCEP_EnergyDepF/1000.;      
      double eneR   = myTV.SCEP_EnergyDepR/1000.;      
//       double cherF  = myTV.SCEP_NCherProdF;      
      double cherR  = myTV.SCEP_NCherProdR;            
      
      //ecal
      double ecal_ene_tot  = (eneR+eneF) / this_ene;                          
//       double ecal_cher_tot = (cherR) / this_ene;          
//       double ecal_S        = ecal_ene_tot/ecal_S_norm;
//       double ecal_C        = ecal_cher_tot/ecal_C_norm;
//       double ecal_C_S      = ecal_C / ecal_S;      
      
      double ecal_S_F = gRandom->Poisson(eneR*LO)/this_ene/ecal_S_norm/LO;
      double ecal_S_R = gRandom->Poisson(eneF*LO)/this_ene/ecal_S_norm/LO;
      double ecal_S   = ecal_S_F+ecal_S_R;
      double ecal_C   = gRandom->Poisson(cherR*LCE_C)/this_ene/CLO;
      double ecal_C_S = ecal_C / ecal_S;   

      //hcal
      float sumScintWBirks = 0;
      float sumCherTowers  = 0;
      for (long unsigned int iT = 0; iT < myTV.VectorSignalsL->size(); iT++)
      {
          sumScintWBirks+= myTV.VectorSignalsL->at(iT) + myTV.VectorSignalsR->at(iT);
          sumCherTowers += myTV.VectorSignalsCherL->at(iT) + myTV.VectorSignalsCherR->at(iT);
      }
      double drh_scint = sumScintWBirks/this_ene;
//       double drh_scint = myTV.EnergyScin/this_ene;
      double drh_cher  = sumCherTowers/this_ene;
//       double drh_cher  = myTV.NofCherenkovDetected/this_ene;
      double drh_S     = drh_scint/drh_S_norm;
      double drh_C     = drh_cher/drh_C_norm;
      double hcal_C_S  = drh_C / drh_S;
      

      //apply dro correction o HCAL energy deposit      
      double drh_S_corr = (drh_S-x_factor_hcal*drh_C )/(1-x_factor_hcal);
//       double drh_S_corr = drh_S *(funcHCAL_DRO_S->Eval(1)/funcHCAL_DRO_S->Eval(hcal_C_S));
      double drh_C_corr = drh_C *(funcHCAL_DRO_C->Eval(1)/funcHCAL_DRO_C->Eval(hcal_C_S));      
      
//       if (drh_S_corr< 0.01 ) std::cout << "drh_S_corr = " << drh_S_corr << " :: drh_C_corr = " << drh_C_corr << " :: hcal_C_S = " << hcal_C_S << std::endl;
//       if (drh_S<0.01 || drh_C<0.01 || ecal_S<0.001 ) continue;
      
      
      
      if (ecal_ene_tot*this_ene>mipInPWO*2 && ecal_ene_tot >0.03//more than 3 mips ene in ECAL
          && hcal_C_S >0
        )
      {

      
//         double hcal_ene_frac = drh_S_corr/funcHCAL_DRO_S->Eval(1);
        double hcal_ene_frac = drh_S_corr;
        float ecal_ene_frac = ecal_S/(1 - hcal_ene_frac ) ;
        if (ecal_ene_frac>0.1) hScatterDRO_ECAL_S -> Fill(ecal_C_S , ecal_ene_frac);
        //         std::cout << "ecal_S = " << ecal_S << " :: hcal_ene_frac = " << hcal_ene_frac << " :: ecal_ene_frac = " << ecal_ene_frac << std::endl;                
        
        hEneTotal_S->Fill(ecal_S+hcal_ene_frac);
        
        double hcal_cher_frac = drh_C_corr/funcHCAL_DRO_C->Eval(1);
        float ecal_cher_frac = ecal_C/(1 - hcal_cher_frac ) ;
        if (ecal_cher_frac>0.1) hScatterDRO_ECAL_C -> Fill(ecal_C_S , ecal_cher_frac);
//         std::cout << "ecal_C = " << ecal_C << " :: hcal_ene_frac = " << hcal_ene_frac << " :: ecal_cher_frac = " << ecal_cher_frac << std::endl;                
        hEneTotal_C->Fill(ecal_C+hcal_cher_frac);
        
        if (this_ene>40)
        {
            hECAL_CE_vs_SE -> Fill(ecal_S/(1 - hcal_ene_frac ), ecal_C/(1 - hcal_ene_frac ));
            hECAL_XFactor->Fill((ecal_S/(1 - hcal_ene_frac )-1)/(ecal_C/(1 - hcal_ene_frac )-1));
        }
        hScatterECAL_frac->Fill(this_ene, ecal_ene_tot);

        for (int iEne = 0; iEne<ENE_BINS; iEne++)
        {
            if (this_ene>=vMinEneBin[iEne] && this_ene<vMaxEneBin[iEne]) 
            {
                vTotEnergy[iEne]->push_back(ecal_S+hcal_ene_frac);
                vTotEnergy_C[iEne]->push_back(ecal_C+hcal_cher_frac);
                
                if (ecal_ene_tot>0.02) 
                {
                    vECALCountNoMips[iEne]->push_back(1);
                    vECALwShowerEnergy[iEne]->push_back(ecal_S);
                }
                
            }
        }        
      }
      
      for (int iEne = 0; iEne<ENE_BINS; iEne++)
      {
         if (this_ene>=vMinEneBin[iEne] && this_ene<vMaxEneBin[iEne]) 
         {
             vECALCountAll[iEne]->push_back(1);
             vInputEnergy[iEne]->push_back(this_ene);
         }
      }
            
      hHCAL_S_corr->Fill(drh_S_corr);
      hHCAL_C_corr->Fill(drh_C_corr);
      
      if (ecal_ene_tot*this_ene<mipInPWO*3 && ecal_ene_tot <0.05)
      {
          hHCALOnly_S_corr->Fill(drh_S_corr);
          hHCALOnly_C_corr->Fill(drh_C_corr);
          
          for (int iEne = 0; iEne<ENE_BINS; iEne++)
          {
              if (this_ene>=vMinEneBin[iEne] && this_ene<vMaxEneBin[iEne]) 
              {
                  vHCALOnly[iEne]->push_back(drh_S);              
                  vHCALCherOnly[iEne]->push_back(drh_C);
                  vHCALOnlyCorr[iEne]->push_back(drh_S_corr);              
              }
          }
      }
      
  }
  
  
  
  
  
  cHCAL_S->cd();
  hHCAL_S_corr->SetLineColor(kGreen+1);
  hHCAL_S_corr->Draw("same");
  hHCALOnly_S_corr->SetLineColor(kBlue+1);
  hHCALOnly_S_corr->Draw("same");
  legHCAL_s->AddEntry(hHCAL_S_corr, "All events (with DRO corr)");
  legHCAL_s->AddEntry(hHCALOnly_S_corr, "no shower in ECAL (with DRO corr)");
  legHCAL_s->Draw();
      
  
  cHCAL_C->cd();
  hHCAL_C_corr->SetLineColor(kGreen+1);
  hHCAL_C_corr->Draw("same");
  hHCALOnly_C_corr->SetLineColor(kBlue+1);
  hHCALOnly_C_corr->Draw("same");
  
    
  TCanvas * cScatterDRO_ECAL_S = new TCanvas ("cScatterDRO_ECAL_S", "cScatterDRO_ECAL_S", 600, 500);
  cScatterDRO_ECAL_S->cd();
  hScatterDRO_ECAL_S->SetStats(0);
  hScatterDRO_ECAL_S->Draw("COLZ");  
  hScatterDRO_ECAL_S->GetXaxis()->SetTitle("C/S_{ECAL}");
  hScatterDRO_ECAL_S->GetYaxis()->SetTitle("S_{ECAL}/ (E_{0} - E_{HCAL})");
  
  TProfile* pDRO_corr_ECAL_S = hScatterDRO_ECAL_S->ProfileX();
  pDRO_corr_ECAL_S->SetStats(0);
  pDRO_corr_ECAL_S->SetLineColor(kRed);
  pDRO_corr_ECAL_S->SetMarkerColor(kRed);
  pDRO_corr_ECAL_S->SetMarkerStyle(20);
  pDRO_corr_ECAL_S->Draw("same");
  hScatterDRO_ECAL_S->GetXaxis()->SetRangeUser(0., 1.2);
  hScatterDRO_ECAL_S->GetYaxis()->SetRangeUser(0., 1.2);
  
  TF1 * funcECAL_DRO_S = new TF1 ("funcECAL_DRO_S", "pol2", 0.4, 1.2);
  funcECAL_DRO_S->SetParameters(1, 1,1);//,1,0.3);    
  funcECAL_DRO_S->SetLineColor(kRed);
  for (int i = 0; i < 10; i++) pDRO_corr_ECAL_S->Fit(funcECAL_DRO_S, "QR");
  funcECAL_DRO_S->Draw("same");      
  gPad->SetGrid();
  
  
  TCanvas * cScatterDRO_ECAL_C = new TCanvas ("cScatterDRO_ECAL_C", "cScatterDRO_ECAL_C", 600, 500);
  cScatterDRO_ECAL_C->cd();
  hScatterDRO_ECAL_C->SetStats(0);
  hScatterDRO_ECAL_C->Draw("COLZ");  
  hScatterDRO_ECAL_C->GetXaxis()->SetTitle("C/S_{ECAL}");
  hScatterDRO_ECAL_C->GetYaxis()->SetTitle("C_{ECAL}/ (E_{0} - E_{HCAL})");
  
  TProfile* pDRO_corr_ECAL_C = hScatterDRO_ECAL_C->ProfileX();
  pDRO_corr_ECAL_C->SetStats(0);
  pDRO_corr_ECAL_C->SetLineColor(kRed);
  pDRO_corr_ECAL_C->SetMarkerColor(kRed);
  pDRO_corr_ECAL_C->SetMarkerStyle(20);
  pDRO_corr_ECAL_C->Draw("same");
  
  TF1 * funcECAL_DRO_C = new TF1 ("funcECAL_DRO_C", "pol2", 0.4, 1.2);
  funcECAL_DRO_C->SetParameters(1, 1,1);//,1,0.3);  
  funcECAL_DRO_C->SetLineColor(kRed);
  for (int i = 0; i < 10; i++) pDRO_corr_ECAL_C->Fit(funcECAL_DRO_C, "QR");
  funcECAL_DRO_C->Draw("same");      
  hScatterDRO_ECAL_C->GetXaxis()->SetRangeUser(0., 1.2);
  hScatterDRO_ECAL_C->GetYaxis()->SetRangeUser(0., 1.2);
  gPad->SetGrid();
  
  
  
  
  TCanvas * cScatterECAL_frac = new TCanvas ("cScatterECAL_frac", "cScatterECAL_frac", 600, 500);
  cScatterECAL_frac->cd();
  hScatterECAL_frac->SetStats(0);
  hScatterECAL_frac->Draw("COLZ");  
  hScatterECAL_frac->GetXaxis()->SetTitle("E [GeV]");
  hScatterECAL_frac->GetYaxis()->SetTitle("E_{ECAL} / E_{0}");  
  
  TProfile* pECAL_frac = hScatterECAL_frac->ProfileX();
  pECAL_frac->SetStats(0);
  pECAL_frac->SetLineColor(kRed);
  pECAL_frac->SetMarkerColor(kRed);
  pECAL_frac->SetMarkerStyle(20);
  pECAL_frac->Draw("same");
  hScatterECAL_frac->GetXaxis()->SetRangeUser(0., 100);
  hScatterECAL_frac->GetYaxis()->SetRangeUser(0., 1.2);
  
  
  
/*  
  
  
  TGraphErrors* gMPV_CE_vs_SE_ECAL = new TGraphErrors();
  TCanvas * cSlices_ECAL = new TCanvas ("cSlices_ECAL", "cSlices_ECAL", 600, 500);
  cSlices_ECAL->cd();
  TH1F* hSlice_ECAL [NBIN_ENE/4];
  int ipoint_ECAL = 0;
  for (int islice = 0; islice< NBIN_ENE/4; islice++)
  {
      hSlice_ECAL[islice] = new TH1F (Form("hSlice_ECAL%d", islice), Form("hSlice_ECAL%d", islice), NBIN_ENE/4, 0, maxEne);
      for (int ibin = 0; ibin< NBIN_ENE/4; ibin++)
      {
          hSlice_ECAL[islice]->SetBinContent(ibin, hECAL_CE_vs_SE->GetBinContent(islice,ibin));
//           std::cout << " islice = " << islice << " :: bincontent = " << hHCAL_CE_vs_SE->GetBinContent(islice,ibin) << std::endl;
      }
      if (islice == 0) 
      {
          hSlice_ECAL[0]->Draw("h");
          hSlice_ECAL[0]->GetXaxis()->SetRangeUser(0, 1.2);
          hSlice_ECAL[0]->GetYaxis()->SetRangeUser(0, 100);
      }
      else 
      {
//           hSlice[islice]->SetLineColor(kRainBow+islice);
          hSlice_ECAL[islice]->Draw("same h");
      }
      
      float sliceMax = hSlice_ECAL[islice]->GetBinCenter(hSlice_ECAL[islice]->GetMaximumBin());
//       TF1 * fitSlice = new TF1 ("fitSlice", "gaus", sliceMax*0.7, sliceMax*1.3);

      if (sliceMax>0.1) 
      {
      
          TF1 * fitSlice = new TF1 ("fitSlice",crystalBallDouble, sliceMax*0.8, sliceMax*1.2, 5);
          fitSlice->SetParameters(1, sliceMax, 0.2, 0.2, 0.2);
      
          hSlice_ECAL[islice]->Fit(fitSlice, "QR");
          float binCenter = maxEne/NBIN_ENE*4*(islice+0.5);
          float mpv = fitSlice->GetParameter(1);
          std::cout << " islice = " << islice << " sliceMax = " << sliceMax <<  " :: binCenter = " << binCenter << " :: mpv = " << mpv << std::endl; 
          gMPV_CE_vs_SE_ECAL->SetPoint(ipoint, binCenter , mpv);
          ipoint++;
      }
  }*/
  
  
  TCanvas * cScatterDRO_ECAL_CE_SE = new TCanvas ("cScatterDRO_ECAL_CE_SE", "cScatterDRO_ECAL_CE_SE", 600, 600);
  cScatterDRO_ECAL_CE_SE->cd();
  hECAL_CE_vs_SE->SetStats(0);
  hECAL_CE_vs_SE->SetTitle("");
  hECAL_CE_vs_SE->Draw("COLZ");  
  hECAL_CE_vs_SE->GetXaxis()->SetTitle("S_{ECAL}/E");
  hECAL_CE_vs_SE->GetYaxis()->SetTitle("C_{ECAL}/E ");
  
  TProfile* pECAL_CE_vs_SE = hECAL_CE_vs_SE->ProfileX();
  pECAL_CE_vs_SE->SetStats(0);
  pECAL_CE_vs_SE->SetLineColor(kRed);
  pECAL_CE_vs_SE->SetMarkerColor(kRed);
  pECAL_CE_vs_SE->SetMarkerStyle(20);
  pECAL_CE_vs_SE->SetMarkerSize(0.8);
  pECAL_CE_vs_SE->Draw("same");
  hECAL_CE_vs_SE->GetXaxis()->SetRangeUser(0., 1.2);
  hECAL_CE_vs_SE->GetYaxis()->SetRangeUser(0., 1.2);
  
  TF1 * funcECAL_CE_vs_SE = new TF1 ("funcECAL_CE_vs_SE", "1-[0]+[0]*x", 0.7, 1.);
  //[0]*log([1]*x)+[2]
  funcECAL_CE_vs_SE->SetParameters(1,1,1);
  funcECAL_CE_vs_SE->SetLineColor(kBlack);
  funcHCAL_CE_vs_SE->SetLineWidth(3);
//   funcECAL_CE_vs_SE->SetLineColor(kRed);
  for (int i = 0; i < 10; i++) pECAL_CE_vs_SE->Fit(funcECAL_CE_vs_SE, "QR");
  funcECAL_CE_vs_SE->SetRange(0.5,1);
  funcECAL_CE_vs_SE->Draw("same");        
  
//   TF1 * CequalS = new TF1("CequalS", "x", 0, 1.1);
//   CequalS->SetLineStyle(7);
//   CequalS->SetLineWidth(2);
//   CequalS->SetLineColor(kBlue+1);
  CequalS->Draw("same");
  gPad->SetGrid();
  
    
//   gMPV_CE_vs_SE_ECAL->Draw("same PE");
//   gMPV_CE_vs_SE_ECAL->SetMarkerStyle(20);
//   gMPV_CE_vs_SE_ECAL->SetMarkerColor(kViolet);
  if (SAVEPLOTS) cScatterDRO_ECAL_CE_SE->SaveAs("plots/cScatterDRO_ECAL_CE_SE.pdf");
  
  
  TCanvas * cECAL_XFactor = new TCanvas ("cECAL_XFactor", "cECAL_XFactor", 600, 600);
  cECAL_XFactor->cd();
  hECAL_XFactor->SetStats(0);
  hECAL_XFactor->SetTitle("");
  hECAL_XFactor->GetXaxis()->SetTitle("#Chi_{ECAL}");
  hECAL_XFactor->GetYaxis()->SetTitle("Counts");
  hECAL_XFactor->SetLineColor(kBlack);
  hECAL_XFactor->Draw();
  if (SAVEPLOTS) cECAL_XFactor->SaveAs("plots/cECAL_XFactor.pdf");
  
  
  
  TCanvas * cEneTotal_S = new TCanvas ("cEneTotal_S", "cEneTotal_S", 600, 500);
  cEneTotal_S->cd();
  hEneTotal_S->Draw();
  gPad->SetLogy();
  

  
  TCanvas * cEneTotal_C = new TCanvas ("cEneTotal_C", "cEneTotal_C", 600, 500);
  cEneTotal_C->cd();
  hEneTotal_C->Draw();
  gPad->SetLogy();
  
  
  
  
  
  
  
  
  //third loop to calculate ecal dro correction
  std::cout << "third loop to calculate combine energy..." << std::endl;
  
  for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
  {
      
      TreeRun->GetEntry(iEvt);                        
      
      double this_ene = myTV.PrimaryParticleEnergy/1000; // in GeV
      double px  = myTV.PrimaryParticleMomentum->at(0);
      double py  = myTV.PrimaryParticleMomentum->at(1);
      double pz  = myTV.PrimaryParticleMomentum->at(2);
      double P   = sqrt(px*px+py*py+pz*pz);
      px/= P;
      py/= P;
      pz/= P;
          
      double eta   = atanh(pz);
      double theta = 2*atan(exp(-eta));
//       double phi   = atan(py/px);
//       double eta   = atanh(pz);
      if (theta <0.5 || theta>2.5) continue;
      if      (this_ene < 12 && myTV.leakage/1000. > 0.5) continue;
      else if (this_ene < 50 && myTV.leakage/1000. > 1.0) continue;
      else if (this_ene > 50 && myTV.leakage/1000. > 3.0) continue;
//       if (this_ene<5) continue;
      
      double eneF   = myTV.SCEP_EnergyDepF/1000.;      
      double eneR   = myTV.SCEP_EnergyDepR/1000.;      
      double cherF  = myTV.SCEP_NCherProdF;      
      double cherR  = myTV.SCEP_NCherProdR;            
      
      //ecal
      double ecal_ene_tot  = (eneR+eneF) / this_ene;                          
      double ecal_cher_tot = (cherR) / this_ene;          
//       double ecal_S        = ecal_ene_tot/ecal_S_norm;
//       double ecal_C        = ecal_cher_tot/ecal_C_norm;
                  
      double ecal_S_F = gRandom->Poisson(eneR*LO)/this_ene/ecal_S_norm/LO;
      double ecal_S_R = gRandom->Poisson(eneF*LO)/this_ene/ecal_S_norm/LO;
      double ecal_S   = ecal_S_F+ecal_S_R;
      double ecal_C   = gRandom->Poisson(cherR*LCE_C)/this_ene/CLO;
      double ecal_C_S = ecal_C / ecal_S;    

      //hcal
      float sumScintWBirks = 0;
      float sumCherTowers  = 0;
      for (long unsigned int iT = 0; iT < myTV.VectorSignalsL->size(); iT++)
      {
          sumScintWBirks+= myTV.VectorSignalsL->at(iT) + myTV.VectorSignalsR->at(iT);
          sumCherTowers += myTV.VectorSignalsCherL->at(iT) + myTV.VectorSignalsCherR->at(iT);
      }
      double drh_scint = sumScintWBirks/this_ene;
//       double drh_scint = myTV.EnergyScin/this_ene;
      double drh_cher  = sumCherTowers/this_ene;
//       double drh_cher  = myTV.NofCherenkovDetected/this_ene;
      double drh_S     = drh_scint/drh_S_norm;
      double drh_C     = drh_cher/drh_C_norm;
      double hcal_C_S  = drh_C / drh_S;
      

      //apply dro correction o HCAL energy deposit
      double drh_S_corr = (drh_S-x_factor_hcal*drh_C )/(1-x_factor_hcal);
//       double drh_S_corr = drh_S *(funcHCAL_DRO_S->Eval(1)/funcHCAL_DRO_S->Eval(hcal_C_S));
      double drh_C_corr = drh_C *(funcHCAL_DRO_C->Eval(1)/funcHCAL_DRO_C->Eval(hcal_C_S));      
      
//       std::cout << "drh_S_corr = " << drh_S_corr << " :: drh_C_corr = " << drh_C_corr << " :: hcal_C_S = " << hcal_C_S << std::endl;
      if (drh_S<0.01 || drh_C<0.01 || ecal_S<0.001 ) continue;
      
                  
//       if (ecal_ene_tot*this_ene>mipInPWO*2 && ecal_ene_tot <0.03
//           && hcal_C_S >0
//         )
      {

      
//         double hcal_ene_frac = drh_S_corr/funcHCAL_DRO_S->Eval(1);
        double hcal_ene_frac = drh_S_corr;
        float ecal_ene_frac = ecal_S/(1 - hcal_ene_frac ) ;        
        
                
        double hcal_cher_frac = drh_C_corr/funcHCAL_DRO_C->Eval(1);
        float ecal_cher_frac = ecal_C/(1 - hcal_ene_frac ) ;
        
        
        double ecal_ene_frac_corr = (ecal_S-x_factor_ecal*ecal_C )/(1-x_factor_ecal);
//         double ecal_ene_frac_corr = ecal_S *(funcECAL_DRO_S->Eval(1)/funcECAL_DRO_S->Eval(ecal_C_S));        
        double tot_Edep_S_corr    = hcal_ene_frac + ecal_ene_frac_corr;
        
//         std::cout << "ecal_S = " << ecal_S << " :: hcal_ene_frac = " << hcal_ene_frac << " :: ecal_ene_frac = " << ecal_ene_frac << " :: ecal_ene_frac_corr = " << ecal_ene_frac_corr <<  std::endl;                
        

        hEneTotal_S_corr->Fill(tot_Edep_S_corr);

        
        
        double ecal_cher_frac_corr = ecal_C *(funcECAL_DRO_C->Eval(1)/funcECAL_DRO_C->Eval(ecal_C_S));        
        double tot_Edep_C_corr    = hcal_cher_frac + ecal_cher_frac_corr;
        hEneTotal_C_corr->Fill(tot_Edep_C_corr);
        
        
        for (int iEne = 0; iEne<ENE_BINS; iEne++)
        {
            if (this_ene>=vMinEneBin[iEne] && this_ene<vMaxEneBin[iEne]) 
            {
                vTotEnergyCorr[iEne]->push_back(tot_Edep_S_corr);
                vTotEnergyCorr_C[iEne]->push_back(tot_Edep_C_corr);
            }
        }
      
      }            
  }
  
  
  
  
  cEneTotal_S->cd();
  hEneTotal_S_corr->SetLineColor(kGreen+1);
  hEneTotal_S_corr->Draw("same");
  
  
  cEneTotal_C->cd();
  hEneTotal_C_corr->SetLineColor(kGreen+1);
  hEneTotal_C_corr->Draw("same");
  
  
  
  
  
  TGraphErrors * gEres_vs_ene       = new TGraphErrors();
  TGraphErrors * gEresCorr_vs_ene   = new TGraphErrors();
  TGraphErrors * gEres_vs_ene_C       = new TGraphErrors();
  TGraphErrors * gEresCorr_vs_ene_C   = new TGraphErrors();
  TGraphErrors * gEres_vs_ene_hcalOnly   = new TGraphErrors();
  TGraphErrors * gEres_vs_ene_hcalOnlyCher   = new TGraphErrors();
  TGraphErrors * gEresCorr_vs_ene_hcalOnly   = new TGraphErrors();
  
  
  TGraphErrors * gElin_vs_ene       = new TGraphErrors();
  TGraphErrors * gElinCorr_vs_ene   = new TGraphErrors();
  TGraphErrors * gElin_vs_ene_C       = new TGraphErrors();
  TGraphErrors * gElinCorr_vs_ene_C   = new TGraphErrors();
  TGraphErrors * gElin_vs_ene_hcalOnly   = new TGraphErrors();
  TGraphErrors * gElin_vs_ene_hcalOnlyCher   = new TGraphErrors();
  TGraphErrors * gElinCorr_vs_ene_hcalOnly   = new TGraphErrors();
  
  TGraphErrors * gFracShowerInECAL = new TGraphErrors();

  
  for (int iEne = 0; iEne<ENE_BINS; iEne++)
  {
      float ave_ene = 0;//(vMinEneBin[iEne]+vMaxEneBin[iEne])/2.;            
      for (long unsigned int it = 0; it < vTotEnergy[iEne]->size(); it++)
      {
          ave_ene += vInputEnergy[iEne]->at(it);
      }
      ave_ene/=vInputEnergy[iEne]->size();
      
      
      double sigma_eff = FindSmallestInterval(vTotEnergy[iEne], 0.68, false) /2. ;
      
      float mean = 0;      
      for (long unsigned int it = 0; it < vTotEnergy[iEne]->size(); it++)
      {
          mean += vTotEnergy[iEne]->at(it);
      }
      mean/=vTotEnergy[iEne]->size();
      gEres_vs_ene->SetPoint(iEne, ave_ene, sigma_eff/mean);
      gElin_vs_ene->SetPoint(iEne, ave_ene, mean);
//       std::cout << " vTotEnergy[" << iEne << "]: counts = " << vTotEnergy[iEne]->size() << " :: mean = " << mean << " :: sigma = "  << sigma_eff << std::endl;
      
      double sigma_eff_corr = FindSmallestInterval(vTotEnergyCorr[iEne], 0.68, false) /2. ;      
      mean = 0;      
      for (long unsigned int it = 0; it < vTotEnergyCorr[iEne]->size(); it++)
      {
          mean += vTotEnergyCorr[iEne]->at(it);
      }
      mean/=vTotEnergyCorr[iEne]->size();
      gEresCorr_vs_ene->SetPoint(iEne, ave_ene, sigma_eff_corr/mean);
      gElinCorr_vs_ene->SetPoint(iEne, ave_ene, mean);
      
      double sigma_eff_C = FindSmallestInterval(vTotEnergy_C[iEne], 0.68, false) /2. ;      
      mean = 0;      
      for (long unsigned int it = 0; it < vTotEnergy_C[iEne]->size(); it++)
      {
          mean += vTotEnergy_C[iEne]->at(it);
      }
      mean/=vTotEnergy_C[iEne]->size();
      gEres_vs_ene_C->SetPoint(iEne, ave_ene, sigma_eff_C/mean);
      gElin_vs_ene_C->SetPoint(iEne, ave_ene, mean);
      
      double sigma_eff_corr_C = FindSmallestInterval(vTotEnergyCorr_C[iEne], 0.68, false) /2. ;      
      mean = 0;      
      for (long unsigned int it = 0; it < vTotEnergyCorr_C[iEne]->size(); it++)
      {
          mean += vTotEnergyCorr_C[iEne]->at(it);
      }
      mean/=vTotEnergyCorr_C[iEne]->size();
      gEresCorr_vs_ene_C->SetPoint(iEne, ave_ene, sigma_eff_corr_C/mean);    
      gElinCorr_vs_ene_C->SetPoint(iEne, ave_ene, mean);
      
      double sigma_eff_HOnly = FindSmallestInterval(vHCALOnly[iEne], 0.68, false) /2. ;      
      mean = 0;      
      for (long unsigned int it = 0; it < vHCALOnly[iEne]->size(); it++)
      {
          mean += vHCALOnly[iEne]->at(it);
      }
      mean/=vHCALOnly[iEne]->size();
      gEres_vs_ene_hcalOnly->SetPoint(iEne, ave_ene, sigma_eff_HOnly/mean);
      gElin_vs_ene_hcalOnly->SetPoint(iEne, ave_ene, mean);
      
      double sigma_eff_HOnlyCher = FindSmallestInterval(vHCALCherOnly[iEne], 0.68, false) /2. ;      
      mean = 0;      
      for (long unsigned int it = 0; it < vHCALCherOnly[iEne]->size(); it++)
      {
          mean += vHCALCherOnly[iEne]->at(it);
      }
      mean/=vHCALCherOnly[iEne]->size();
      gEres_vs_ene_hcalOnlyCher->SetPoint(iEne, ave_ene, sigma_eff_HOnlyCher/mean);
      gElin_vs_ene_hcalOnlyCher->SetPoint(iEne, ave_ene, mean);
      
      double sigma_eff_HOnlyCorr = FindSmallestInterval(vHCALOnlyCorr[iEne], 0.68, false) /2. ;            
      mean = 0;      
      for (long unsigned int it = 0; it < vHCALOnlyCorr[iEne]->size(); it++)
      {
          mean += vHCALOnlyCorr[iEne]->at(it);
      }
      mean/=vHCALOnlyCorr[iEne]->size();
      gEresCorr_vs_ene_hcalOnly->SetPoint(iEne, ave_ene, sigma_eff_HOnlyCorr/mean);
      gElinCorr_vs_ene_hcalOnly->SetPoint(iEne, ave_ene, mean);
      std::cout << " vHCALOnlyCorr[" << iEne << "]: counts = " << vHCALOnlyCorr[iEne]->size() << " :: mean = " << mean << " :: sigma = "  << sigma_eff_HOnlyCorr << std::endl;
      
      float frac = float(vECALCountNoMips[iEne]->size())/float(vECALCountAll[iEne]->size());
      std::cout << "frac in ECAL = "  << frac << std::endl;
      gFracShowerInECAL->SetPoint(iEne, ave_ene, frac);
      
  }
  
  
  TCanvas * cFracShowerInECAL = new TCanvas ("cFracShowerInECAL", "cFracShowerInECAL", 600, 600);
  cFracShowerInECAL->cd();
  gFracShowerInECAL->Sort();
  gFracShowerInECAL->Draw("ALPE");
  gFracShowerInECAL->SetTitle(";E_{0} [GeV]; Fraction of hadron showers starting in ECAL");

  

  
  TCanvas * cEnergyResolutionHCAL_Only = new TCanvas ("cEnergyResolutionHCAL_Only", "cEnergyResolutionHCAL_Only", 600, 600);
  cEnergyResolutionHCAL_Only->cd();
  gEres_vs_ene_hcalOnly->Sort();
  gEres_vs_ene_hcalOnly->SetLineColor(kGray+1);
  gEres_vs_ene_hcalOnly->SetMarkerColor(kGray+1);
  gEres_vs_ene_hcalOnly->SetMarkerStyle(21);
  gEres_vs_ene_hcalOnly->Draw("ALPE");
  gEres_vs_ene_hcalOnly->SetTitle("; Particle energy [GeV] ; #sigma_{E}/E");
  gEres_vs_ene_hcalOnly->SetMaximum(0.4);
  gEres_vs_ene_hcalOnly->SetMinimum(0.01);
  gEres_vs_ene_hcalOnly->GetXaxis()->SetLimits(8, 100);
  
    
  gEres_vs_ene_hcalOnlyCher->SetLineColor(kCyan+1);
  gEres_vs_ene_hcalOnlyCher->SetMarkerColor(kCyan+1);
  gEres_vs_ene_hcalOnlyCher->SetMarkerStyle(22);
  gEres_vs_ene_hcalOnlyCher->Draw("same LPE");
  
  gEresCorr_vs_ene_hcalOnly->SetLineColor(kBlack);
  gEresCorr_vs_ene_hcalOnly->SetMarkerStyle(20);
  gEresCorr_vs_ene_hcalOnly->Draw("same lPE");
  
  
  leg = new TLegend(0.15,0.15,0.75,0.35,NULL,"brNDC");
  
  TF1 * fitResolution = new TF1 ("fitResolution", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 4, 100);    
  fitResolution->SetLineWidth(2);    
  fitResolution->SetLineStyle(7);    
  fitResolution->SetParameters(0.3, 0.01);
  fitResolution->SetParLimits(1,0.014, 0.08);
    
  
  fitResolution->SetLineColor(kGray+1);    
  gEres_vs_ene_hcalOnly->Fit(fitResolution, "QRS");
  float stoch_term = fitResolution->GetParameter(0);
  float const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gEres_vs_ene_hcalOnly, Form("HCAL S only:  %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
  
  
  fitResolution->SetLineColor(kCyan);    
  for (int i = 0; i<5; i++) gEres_vs_ene_hcalOnlyCher->Fit(fitResolution, "QRS");
  stoch_term = fitResolution->GetParameter(0);
  const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gEres_vs_ene_hcalOnlyCher, Form("HCAL C only:  %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");  
  
  fitResolution->SetLineColor(kBlack);    
  for (int i = 0; i<5; i++) gEresCorr_vs_ene_hcalOnly->Fit(fitResolution, "QRS");
  stoch_term = fitResolution->GetParameter(0);
  const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gEresCorr_vs_ene_hcalOnly, Form("HCAL only (DRO corr):  %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
            
  
  leg->Draw();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  
  
  
  
  TCanvas * cEnergyResolution = new TCanvas ("cEnergyResolution", "cEnergyResolution", 600, 600);
  cEnergyResolution->cd();
  gEres_vs_ene->Sort();
  gEres_vs_ene->SetLineColor(kRed+1);
  gEres_vs_ene->SetMarkerColor(kRed+1);
  gEres_vs_ene->SetMarkerStyle(21);
  gEres_vs_ene->Draw("APE");
  gEres_vs_ene->SetTitle("; Particle energy [GeV] ; #sigma_{E}/E");
  gEres_vs_ene->SetMaximum(0.2);
  gEres_vs_ene->SetMinimum(0.01);
  gEres_vs_ene->GetXaxis()->SetLimits(8, 100);
  
  gEresCorr_vs_ene->SetLineColor(kGreen+1);
  gEresCorr_vs_ene->SetMarkerColor(kGreen+1);
  gEresCorr_vs_ene->SetMarkerStyle(20);
  gEresCorr_vs_ene->Draw("same PE");

  gEresCorr_vs_ene_hcalOnly->SetLineColor(kBlack);
  gEresCorr_vs_ene_hcalOnly->SetMarkerStyle(20);
  gEresCorr_vs_ene_hcalOnly->Draw("same PE");
    
  
  leg = new TLegend(0.15,0.15,0.75,0.35,NULL,"brNDC");
  
  fitResolution = new TF1 ("fitResolution", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 4, 100);    
  fitResolution->SetLineWidth(2);    
  fitResolution->SetLineStyle(7);    
  fitResolution->SetParameters(0.3, 0.01);
  fitResolution->SetParLimits(1,0.014, 0.08);
  
  fitResolution->SetLineColor(kRed);    
  for (int i = 0; i<5; i++) gEres_vs_ene->Fit(fitResolution, "QRS");
  stoch_term = fitResolution->GetParameter(0);
  const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gEres_vs_ene, Form("ECAL+HCAL (no DRO):  %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
  fitResolution->Draw("same");
  
  
  fitResolution->SetLineColor(kGreen+1);    
  for (int i = 0; i<5; i++) gEresCorr_vs_ene->Fit(fitResolution, "QRS");
  stoch_term = fitResolution->GetParameter(0);
  const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gEresCorr_vs_ene, Form("ECAL+HCAL (with DRO):  %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
  
  
  fitResolution->SetLineColor(kBlack);    
  fitResolution->SetLineStyle(7);    
  for (int i = 0; i<5; i++) gEresCorr_vs_ene_hcalOnly->Fit(fitResolution, "QRS");
  stoch_term = fitResolution->GetParameter(0);
  const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gEresCorr_vs_ene_hcalOnly, Form("HCAL only (with DRO):  %.2f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
            
  
  leg->Draw();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  if (SAVEPLOTS) cEnergyResolution->SaveAs("plots/cHadronEnergyResolution_wCrystals.pdf");
  
  
  
  
  
  TCanvas * cEnergyLinearityHCAL_Only = new TCanvas ("cEnergyLinearityHCAL_Only", "cEnergyLinearityHCAL_Only", 600, 600);
  cEnergyLinearityHCAL_Only->cd();
  gElin_vs_ene_hcalOnly->Sort();
  gElin_vs_ene_hcalOnly->SetLineColor(kGray+1);
  gElin_vs_ene_hcalOnly->SetMarkerColor(kGray+1);
  gElin_vs_ene_hcalOnly->SetMarkerStyle(21);
  gElin_vs_ene_hcalOnly->Draw("AlPE");
  gElin_vs_ene_hcalOnly->SetTitle("; Particle energy [GeV] ; E_{reco}/E_{truth}");
  gElin_vs_ene_hcalOnly->SetMaximum(1.1);
  gElin_vs_ene_hcalOnly->SetMinimum(0.4);
  gElin_vs_ene_hcalOnly->GetXaxis()->SetLimits(0, 100);
  
  gElin_vs_ene_hcalOnlyCher->SetLineColor(kCyan+1);
  gElin_vs_ene_hcalOnlyCher->SetMarkerColor(kCyan+1);
  gElin_vs_ene_hcalOnlyCher->SetMarkerStyle(22);
  gElin_vs_ene_hcalOnlyCher->Draw("same lPE");
  
  gElinCorr_vs_ene_hcalOnly->SetLineColor(kBlack);
  gElinCorr_vs_ene_hcalOnly->SetMarkerStyle(20);
  gElinCorr_vs_ene_hcalOnly->Draw("same lPE");
  
  
  leg = new TLegend(0.45,0.15,0.88,0.35,NULL,"brNDC");  
  leg->AddEntry(gElin_vs_ene_hcalOnly, "HCAL S only", "lp");    
  leg->AddEntry(gElin_vs_ene_hcalOnlyCher, "HCAL C only", "lp");    
  leg->AddEntry(gElinCorr_vs_ene_hcalOnly, "HCAL only (DRO corr)", "lp");    
  leg->Draw();
  gPad->SetGridy();
  
  
  
  TCanvas * cEnergyLinearity = new TCanvas ("cEnergyLinearity", "cEnergyLinearity", 600, 600);
  cEnergyLinearity->cd();  
  gElin_vs_ene->Sort();
  gElin_vs_ene->SetLineColor(kRed+1);
  gElin_vs_ene->SetMarkerColor(kRed+1);
  gElin_vs_ene->SetMarkerStyle(21);
  gElin_vs_ene->Draw("ALPE");
  gElin_vs_ene->SetTitle("; Particle energy [GeV] ; E_{reco}/E_{truth}");
  gElin_vs_ene->SetMaximum(1.1);
  gElin_vs_ene->SetMinimum(0.7);
  gElin_vs_ene->GetXaxis()->SetLimits(0, 100);
  
  gElinCorr_vs_ene->SetLineColor(kGreen+1);
  gElinCorr_vs_ene->SetMarkerColor(kGreen+1);
  gElinCorr_vs_ene->SetMarkerStyle(20);
  gElinCorr_vs_ene->Draw("same lPE");

  gElinCorr_vs_ene_hcalOnly->SetLineColor(kBlack);
  gElinCorr_vs_ene_hcalOnly->SetMarkerStyle(20);
  gElinCorr_vs_ene_hcalOnly->Draw("same lPE");
  
  
  leg = new TLegend(0.45,0.15,0.88,0.35,NULL,"brNDC");  
  leg->AddEntry(gElin_vs_ene, "ECAL+HCAL (w/o DRO)", "lp");    
  leg->AddEntry(gElinCorr_vs_ene, "ECAL+HCAL (w/ DRO)", "lp");    
  leg->AddEntry(gElinCorr_vs_ene_hcalOnly, "HCAL only (DRO corr)", "lp");    
  leg->Draw();
  gPad->SetGridy();
  
  
//   TF1 * fitLinearity = new TF1 ("fitLinearity", "[0]", 0, 100);
//   fitLinearity->SetParameter(0,1);
//   fitLinearity->SetLineColor(kBlack);
// //   gSlin_vs_ene->Fit(fitLinearity, "Q");
//   
//   TLine * upLim  = new TLine(0, fitLinearity->GetParameter(0)+0.01, 100, fitLinearity->GetParameter(0)+0.01);
//   TLine * lowLim = new TLine(0, fitLinearity->GetParameter(0)-0.01, 100, fitLinearity->GetParameter(0)-0.01);
//   
//   upLim->SetLineStyle(7);
//   upLim->SetLineColor(kGray+1);
//   upLim->SetLineWidth(3);
//   
//   lowLim->SetLineStyle(7);
//   lowLim->SetLineColor(kGray+1);
//   lowLim->SetLineWidth(3);
//   
//   upLim->Draw();
//   lowLim->Draw();
//   
  
  
  if (SAVEPLOTS) cEnergyLinearity->SaveAs("plots/cHadronEnergyLinearity_wCrystals.pdf");
  
  
  
  
    
   
   theApp->Run();
}

