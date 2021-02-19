// g++ -Wall -o plotEMResolution plotEMResolution.C  VectorSmallestInterval.cc myG4Tree.cc myG4Tree.hh `root-config --cflags --glibs`


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
  

  
  std::string particle_name = "e-";

  
  
  int ENE_BINS = 14;  
  float maxEneRange = 205;
  float minEneRange = 0;
  float vMinEneBin [ENE_BINS];
  float vMaxEneBin [ENE_BINS];
  std::vector<double> *vTotEnergy[ENE_BINS];
  std::vector<double> *vTotS[ENE_BINS];
  std::vector<double> *vInputEnergy[ENE_BINS];
  
  float binBoundary[ENE_BINS] = {0,1.5, 2.5, 4.5, 5.5, 8.5, 10.5, 15, 30, 40, 60, 80, 149, 195};
    
  for (int iEne = 0; iEne<ENE_BINS; iEne++)
  {
//       float stepEne = (maxEneRange-minEneRange)/ENE_BINS;
//       float min_ene = minEneRange+stepEne*iEne;
//       float max_ene = min_ene+stepEne;
      
      vMinEneBin[iEne]=binBoundary[iEne];
      if (iEne<ENE_BINS-1) vMaxEneBin[iEne]=binBoundary[iEne+1];
      else                 vMaxEneBin[iEne]=210;
      
      vTotEnergy[iEne] = new std::vector<double>();
      vTotS[iEne]      = new std::vector<double>();
      vInputEnergy[iEne]     = new std::vector<double>();
      
  }

  double ecal_S_norm = 0.985;
  double ecal_C_norm = 5460;
  float LO = 2000;
  float LCE_C = 160/ecal_C_norm;
  
  //define histos
  int NBIN_ENE = 500;
  
  double maxEneECAL  = 1.2;
  double maxCherECAL = 10000;
  
  float minS = 0.;
  float maxS = 2;
  float maxC = 2;
  
  TH1F * hECAL_Scint = new TH1F ("hECAL_Scint", "hECAL_Scint", NBIN_ENE, 0, maxEneECAL);
  TH1F * hECAL_Cher  = new TH1F ("hECAL_Cher", "hECAL_Cher", NBIN_ENE, 0, maxCherECAL);
  
  TH1F * hECAL_S = new TH1F ("hECAL_S", "hECAL_S", NBIN_ENE, 0, maxS);
  TH1F * hECAL_C = new TH1F ("hECAL_C", "hECAL_C", NBIN_ENE, 0, maxC);
  
  TH2F * hScatterCS = new TH2F ("hScatterCS", "hScatterCS", NBIN_ENE, 0, maxS, NBIN_ENE, 0, maxC);
  

  
  //run over energy scan

//   TFile * RunFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_e-_Iso+Uniform1-100_GeV.root","READ");      
//   TFile * RunFile = new TFile("../root_files/prod/output_SCEPCal_fixedPos_gamma_Iso+Uniform1-100_GeV.root","READ");      
          
//   TTree* TreeRun = (TTree*) RunFile->Get("B4");
  
  TChain * TreeRun = new TChain("B4", "B4");
  TreeRun->Add("../root_files/prod/output_SCEPCal_fixedPos_e-_Iso+Uniform1-100_GeV.root");
  TreeRun->Add("../root_files/prod/output_SCEPCal_fixedEne_e-_IsoSupplement.root");

  
  
  myG4TreeVars myTV;
  InitG4Tree (TreeRun, myTV);
  
  TreeRun->SetBranchStatus("*", 0);
  
  TreeRun->SetBranchStatus("PrimaryParticleEnergy", 1);  
  TreeRun->SetBranchStatus("PrimaryParticleMomentum", 1);  
  TreeRun->SetBranchStatus("SCEP_EnergyDepF", 1);
  TreeRun->SetBranchStatus("SCEP_EnergyDepR", 1);
  TreeRun->SetBranchStatus("SCEP_NCherProdF", 1);
  TreeRun->SetBranchStatus("SCEP_NCherProdR", 1);
            
  
  ///*******************************************///
  ///		 Run over events	    ///
  ///*******************************************///
    
  
  int NEVENTS = TreeRun->GetEntries();
  std::cout << "NEVENTS = " << NEVENTS << std::endl;
  
//   NEVENTS = 50000;
      
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

      
//       if (this_ene<5) continue;
      
      double eneF   = myTV.SCEP_EnergyDepF/1000.;      
      double eneR   = myTV.SCEP_EnergyDepR/1000.;      
      double cherF  = myTV.SCEP_NCherProdF;      
      double cherR  = myTV.SCEP_NCherProdR;            
      
      //ecal
      double ecal_ene_tot  = (eneR+eneF) / this_ene;                          
      double ecal_cher_tot = (cherR+cherF) / this_ene;    
            
      
      hECAL_Scint->Fill(ecal_ene_tot);
      hECAL_Cher->Fill(ecal_cher_tot);
      
      double ecal_S_F = gRandom->Poisson(eneR*LO)/this_ene;
      double ecal_S_R = gRandom->Poisson(eneF*LO)/this_ene;
      double ecal_S   = ecal_S_F+ecal_S_R;
      double ecal_C   = gRandom->Poisson(cherR*LCE_C)/this_ene;
      
      
      hECAL_S->Fill(ecal_S/LO/ecal_S_norm);
      hECAL_C->Fill(ecal_C/LCE_C/ecal_C_norm);
      
      hScatterCS->Fill(ecal_S_R/LO/ecal_S_norm, ecal_C/LCE_C/ecal_C_norm);
      
         
      for (int iEne = 0; iEne<ENE_BINS; iEne++)
      {
          
        if (this_ene>=vMinEneBin[iEne] && this_ene<vMaxEneBin[iEne])              
        {
            vTotEnergy[iEne]->push_back(ecal_ene_tot);            
            vTotS[iEne]->push_back(ecal_S/LO/ecal_S_norm);
            vInputEnergy[iEne]->push_back(this_ene);
//             vTotEnergy_C[iEne]->push_back(ecal_C+hcal_cher_frac);                        
        }
      }
      
      
                                
  }
  


  //ECAL segment
  
  TCanvas * cECAL_Scint = new TCanvas ("cECAL_Scint", "cECAL_Scint", 600, 500);
  cECAL_Scint->cd();
  hECAL_Scint->Draw();
  hECAL_Scint->GetXaxis()->SetTitle("E_{dep} / E_{truth}");
  gPad->SetLogy();
  
  TCanvas * cECAL_Cher = new TCanvas ("cECAL_Cher", "cECAL_Cher", 600, 500);
  cECAL_Cher->cd();
  hECAL_Cher->Draw();
  hECAL_Cher->GetXaxis()->SetTitle("C_{prod} / E_{truth}");
  gPad->SetLogy();
  
  TCanvas * cECAL_S = new TCanvas ("cECAL_S", "cECAL_S", 600, 500);
  cECAL_S->cd();
  hECAL_S->Draw();
  hECAL_S->GetXaxis()->SetTitle("S / E_{dep}");
  gPad->SetLogy();
  
  TCanvas * cECAL_C = new TCanvas ("cECAL_C", "cECAL_C", 600, 500);
  cECAL_C->cd();
  hECAL_C->Draw();
  hECAL_C->GetXaxis()->SetTitle("C / E_{dep}");
  gPad->SetLogy();
  
  TCanvas * cScatterCS = new TCanvas ("cScatterCS", "cScatterCS", 600, 600);
  cScatterCS->cd();
  hScatterCS->Draw("COLZ");
  

  
  TGraphErrors * gEres_vs_ene       = new TGraphErrors();    
  TGraphErrors * gElin_vs_ene       = new TGraphErrors();
  TGraphErrors * gSres_vs_ene       = new TGraphErrors();    
  TGraphErrors * gSlin_vs_ene       = new TGraphErrors();
  
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
      

      sigma_eff = FindSmallestInterval(vTotS[iEne], 0.68, false) /2. ;
      
      mean = 0;      
      for (long unsigned int it = 0; it < vTotS[iEne]->size(); it++)
      {
          mean += vTotS[iEne]->at(it);
      }
      mean/=vTotS[iEne]->size();
//       std::cout << "meanS = " << mean << std::endl;
      gSres_vs_ene->SetPoint(iEne, ave_ene, sigma_eff/mean);
      gSlin_vs_ene->SetPoint(iEne, ave_ene, mean);
      
  }
  
  
  
  
  

  
  TCanvas * cEnergyResolution = new TCanvas ("cEnergyResolution", "cEnergyResolution", 600, 600);
  cEnergyResolution->cd();
  gEres_vs_ene->SetLineColor(kRed+1);
  gEres_vs_ene->SetMarkerColor(kRed+1);
  gEres_vs_ene->SetMarkerStyle(21);
  gEres_vs_ene->Draw("APE");
  gEres_vs_ene->SetTitle("; Electron energy [GeV] ; #sigma_{E}/E");
  gEres_vs_ene->SetMaximum(0.03);
  gEres_vs_ene->SetMinimum(0.00);
  gEres_vs_ene->GetXaxis()->SetLimits(1, 210);

  gSres_vs_ene->Draw("same LPE"); 
  gSres_vs_ene->SetLineColor(kBlack);
  gSres_vs_ene->SetMarkerColor(kBlack);
  gSres_vs_ene->SetMarkerStyle(20);
  
  leg = new TLegend(0.15,0.7,0.88,0.88,NULL,"brNDC");
  
  TF1 * fitResolution = new TF1 ("fitResolution", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 1, 200);
  fitResolution->SetLineWidth(2);    
  fitResolution->SetLineStyle(7);    
  fitResolution->SetParameters(0.03, 0.001);
  
  fitResolution->SetLineColor(kBlack);    
  gSres_vs_ene->Fit(fitResolution, "QRS");
  float stoch_term = fitResolution->GetParameter(0);
  float const_term = fitResolution->GetParameter(1);        
  leg->AddEntry(gSres_vs_ene, Form("Total resolution:  %.3f /#sqrt{E} #oplus %.3f", stoch_term, const_term), "lp");    
  leg->AddEntry(gEres_vs_ene, "Contribution from shower containment", "lp");    
  fitResolution->Draw("same");
  
  
//   fitResolution = new TF1 ("fitResolution", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 1, 100);    
  
  leg->Draw();
  gPad->SetGridy();
//   gPad->SetLogy();
  gPad->SetLogx();
  cEnergyResolution->SaveAs("plots/SCEPCal_EM_eres.pdf");

  
  TCanvas * cEnergyLinearity = new TCanvas ("cEnergyLinearity", "cEnergyLinearity", 600, 600);
  cEnergyLinearity->cd();  
  gSlin_vs_ene->SetLineColor(kBlack);
  gSlin_vs_ene->SetMarkerColor(kBlack);
  gSlin_vs_ene->SetMarkerStyle(21);
  gSlin_vs_ene->Draw("APE");
  gSlin_vs_ene->SetTitle("; Electron energy [GeV] ; E_{reco}/E_{truth}");
  gSlin_vs_ene->SetMaximum(1.1);
  gSlin_vs_ene->SetMinimum(0.9);
  gSlin_vs_ene->GetXaxis()->SetLimits(0, 210);
  
  gPad->SetGridy();
  
  TF1 * fitLinearity = new TF1 ("fitLinearity", "pol0", 0, 200);
  fitLinearity->SetLineColor(kBlack);
  gSlin_vs_ene->Fit(fitLinearity, "Q");
  
  TLine * upLim  = new TLine(0, fitLinearity->GetParameter(0)+0.01, 200, fitLinearity->GetParameter(0)+0.01);
  TLine * lowLim = new TLine(0, fitLinearity->GetParameter(0)-0.01, 200, fitLinearity->GetParameter(0)-0.01);
  
  upLim->SetLineStyle(7);
  upLim->SetLineColor(kGray+1);
  upLim->SetLineWidth(3);
  
  lowLim->SetLineStyle(7);
  lowLim->SetLineColor(kGray+1);
  lowLim->SetLineWidth(3);
  
  upLim->Draw();
  lowLim->Draw();
  
  
  cEnergyLinearity->SaveAs("plots/SCEPCal_EM_elin.pdf");
  
  
  
  
  theApp->Run();
}

