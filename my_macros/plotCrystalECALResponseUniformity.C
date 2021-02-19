// g++ -Wall -o plotCrystalECALResponseUniformity plotCrystalECALResponseUniformity.C  VectorSmallestInterval.cc myG4Tree.cc myG4Tree.hh `root-config --cflags --glibs`


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
#include "TTree.h"
#include "TBranch.h"

#include "TSpline.h"
#include "TCanvas.h"
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
  gStyle->SetLegendTextSize(0.035);
  TLegend * leg;
  
  
  
  
  const int NFILES = 1;  
  int energies[NFILES] = {40};
  
  
  std::cout << " NFILES = " << NFILES << std::endl;
  
  std::string particle_name = "e-";
  //define vector energy
  std::vector<double> *vEcalEnergy[NFILES];
  float nCherPerGeV = 7300*4;
  
  float LO       = 2000;
  float LCE_C    = 160./7300.;
  
  
  //define histos
  TH1F * hSCEP_EneF[NFILES];
  TH1F * hSCEP_EneR[NFILES];
  TH1F * hSCEP_EneT[NFILES];
  
  TH1F * hSCEP_CherF[NFILES];
  TH1F * hSCEP_CherR[NFILES];
  TH1F * hSCEP_CherT[NFILES];
  
  TH2F * hScatterCS[NFILES];
  
  TProfile* hEne_vs_Phi[NFILES];
  TProfile* hEne_vs_Theta[NFILES];
  
  TProfile* hSF_vs_Phi[NFILES];
  TProfile* hSF_vs_Theta[NFILES];
  
  TProfile* hSR_vs_Phi[NFILES];
  TProfile* hSR_vs_Theta[NFILES];
  
  TProfile* hS_vs_Phi[NFILES];
  TProfile* hS_vs_Theta[NFILES];
  
  TProfile* hC_vs_Phi[NFILES];
  TProfile* hC_vs_Theta[NFILES];
  
  
  int NBIN_ENE = 200;
  double maxEne = 45;
  
  for (int iFile = 0; iFile<NFILES; iFile++)
  {
      
      std::cout << "energies[" << iFile << "] = " << energies[iFile] << std::endl;
      vEcalEnergy[iFile]  = new std::vector<double>();
      
      double minEne = 0;
//       minEne = energies[iFile]*0.7;
//       maxEne = energies[iFile]*1.2;
      
      
      hSCEP_EneF[iFile] = new TH1F (Form("hSCEP_EneF_%d", energies[iFile]), Form("hSCEP_EneF_%d", energies[iFile]), NBIN_ENE, minEne, maxEne);
      hSCEP_EneR[iFile] = new TH1F (Form("hSCEP_EneR_%d", energies[iFile]), Form("hSCEP_EneR_%d", energies[iFile]), NBIN_ENE, minEne, maxEne);
      hSCEP_EneT[iFile] = new TH1F (Form("hSCEP_EneT_%d", energies[iFile]), Form("hSCEP_EneT_%d", energies[iFile]), NBIN_ENE, minEne, maxEne);
      
      hSCEP_CherF[iFile] = new TH1F (Form("hSCEP_CherF_%d", energies[iFile]), Form("hSCEP_CherF_%d", energies[iFile]), NBIN_ENE, 0, energies[iFile]*1.5*nCherPerGeV);
      hSCEP_CherR[iFile] = new TH1F (Form("hSCEP_CherR_%d", energies[iFile]), Form("hSCEP_CherR_%d", energies[iFile]), NBIN_ENE, 0, energies[iFile]*1.5*nCherPerGeV);
      hSCEP_CherT[iFile] = new TH1F (Form("hSCEP_CherT_%d", energies[iFile]), Form("hSCEP_CherT_%d", energies[iFile]), NBIN_ENE, 0, energies[iFile]*1.5*nCherPerGeV);
      
      hScatterCS[iFile] = new TH2F (Form("hScatterCS_%d", energies[iFile]), Form("hScatterCS_%d", energies[iFile]), NBIN_ENE, 0, 3000, NBIN_ENE, 0, 260);
      
      hEne_vs_Phi[iFile] = new TProfile (Form("hEne_vs_Phi_%d", energies[iFile]), Form("hEne_vs_Phi_%d", energies[iFile]), 50, -3.14/2, -3.14/2);
      hEne_vs_Theta[iFile] = new TProfile (Form("hEne_vs_Theta_%d", energies[iFile]), Form("hEne_vs_Theta_%d", energies[iFile]), 50, 0, 3.14);
      
      
      hSF_vs_Phi[iFile] = new TProfile (Form("hSF_vs_Phi_%d", energies[iFile]), Form("hSF_vs_Phi_%d", energies[iFile]), 50, -3.14/2, -3.14/2);
      hSF_vs_Theta[iFile] = new TProfile (Form("hSF_vs_Theta_%d", energies[iFile]), Form("hSF_vs_Theta_%d", energies[iFile]), 50, 0, 3.14);
      
      
      hSR_vs_Phi[iFile] = new TProfile (Form("hSR_vs_Phi_%d", energies[iFile]), Form("hSR_vs_Phi_%d", energies[iFile]), 50, -3.14/2, -3.14/2);
      hSR_vs_Theta[iFile] = new TProfile (Form("hSR_vs_Theta_%d", energies[iFile]), Form("hSR_vs_Theta_%d", energies[iFile]), 50, 0, 3.14);
      
      
      hS_vs_Phi[iFile] = new TProfile (Form("hS_vs_Phi_%d", energies[iFile]), Form("hS_vs_Phi_%d", energies[iFile]), 50, -3.14/2, -3.14/2);
      hS_vs_Theta[iFile] = new TProfile (Form("hS_vs_Theta_%d", energies[iFile]), Form("hS_vs_Theta_%d", energies[iFile]), 50, 0, 3.14);
      
      
      hC_vs_Phi[iFile] = new TProfile (Form("hC_vs_Phi_%d", energies[iFile]), Form("hC_vs_Phi_%d", energies[iFile]), 50, -3.14/2, -3.14/2);
      hC_vs_Theta[iFile] = new TProfile (Form("hC_vs_Theta_%d", energies[iFile]), Form("hC_vs_Theta_%d", energies[iFile]), 50, 0, 3.14);
      
  }
  
  
        
  //run over energy scan
  TFile * RunFile[NFILES];  
  
  for (int iFile = 0; iFile < NFILES; iFile++)
  {
//       RunFile[iFile] = new TFile(Form("../root_files/merged/ele_%dGeV.root", energies[iFile]),"READ"); 
//       RunFile[iFile] = new TFile("../root_files/B4_t0.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_mu-_50GeV.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_ele_10GeV.root","READ"); 
//       RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_gamma_10GeV_job0.root","READ"); 
//         RunFile[iFile] = new TFile("../root_files/iso_gun/iso_gun_ele_10GeV_cher.root","READ"); 
        RunFile[iFile] = new TFile("../root_files/prod/output_SCEPCal_Iso_e-_40GeV.root","READ"); 
      
      
      TTree* TreeRun = (TTree*) RunFile[iFile]->Get("B4");
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

//       NEVENTS = 10000;
      
      for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
      {
          
          
          
          
          TreeRun->GetEntry(iEvt);                  
//           std::cout << "iEvent = " << iEvt << std::endl;
                //ecal
          double this_ene = myTV.PrimaryParticleEnergy/1000; // in GeV          
          double px  = myTV.PrimaryParticleMomentum->at(0);
          double py  = myTV.PrimaryParticleMomentum->at(1);
          double pz  = myTV.PrimaryParticleMomentum->at(2);
          double P   = sqrt(px*px+py*py+pz*pz);
          px/= P;
          py/= P;
          pz/= P;
          
          if (this_ene < 4) std::cout << "this_ene = " << this_ene<< std::endl;
          double phi   = atan(py/px);
          double eta   = atanh(pz);
          double theta = 2*atan(exp(-eta));
          //       std::cout << " theta = " << theta << std::endl;
//           if (theta <0.5 || theta>2.5) continue;
      
      
          
          double eneF   = myTV.SCEP_EnergyDepF/1000.;
//           std::cout << "eneF = " << eneF << std::endl;
          double eneR   = myTV.SCEP_EnergyDepR/1000.;
//           std::cout << "eneR = " << eneR << std::endl;
          double cherF  = myTV.SCEP_NCherProdF;
//           std::cout << "cherF = " << cherF << std::endl;
          double cherR  = myTV.SCEP_NCherProdR;
//           std::cout << "cherR = " << cherR << std::endl;          
          
          double ene_tot = eneR+eneF;
//           hSCEP_EneF[iFile] ->Fill(eneF/energies[iFile]);
//           hSCEP_EneR[iFile] ->Fill(eneR/energies[iFile]);
//           hSCEP_EneT[iFile] ->Fill(ene_tot/energies[iFile]);
          
          hSCEP_EneF[iFile] ->Fill(eneF);
          hSCEP_EneR[iFile] ->Fill(eneR);
          hSCEP_EneT[iFile] ->Fill(ene_tot);
          vEcalEnergy[iFile]->push_back(ene_tot/energies[iFile]);
//           
// //           std::cout << "filled ene hisots... " << std::endl;
          
          float S_F = gRandom->Poisson(eneF*LO)/this_ene;
          float S_R = gRandom->Poisson(eneR*LO)/this_ene;          
          float S   = S_F+S_R;
          
          
          
          hSCEP_CherF[iFile]->Fill(cherF);
          hSCEP_CherR[iFile]->Fill(cherR);
          hSCEP_CherT[iFile]->Fill(cherR+cherF);
          
          
          hEne_vs_Phi[iFile]->Fill(phi, ene_tot);
          hEne_vs_Theta[iFile]->Fill(theta, ene_tot);
          hSF_vs_Phi[iFile]->Fill(phi, S_F);
          hSF_vs_Theta[iFile]->Fill(theta, S_F);
          
          hSR_vs_Phi[iFile]->Fill(phi, S_R);
          hSR_vs_Theta[iFile]->Fill(theta, S_R);
          hS_vs_Phi[iFile]->Fill(phi, S);
          hS_vs_Theta[iFile]->Fill(theta, S);

          if (eneR>0) 
          {
            float C   = gRandom->Poisson(cherR*LCE_C)/eneR;          
            hC_vs_Phi[iFile]->Fill(phi, C);
            hC_vs_Theta[iFile]->Fill(theta, C);
            hScatterCS[iFile]->Fill(S_R , C);
          }
          
          
      }
      
      
      
  }
  
  
  int selEne = 0;
  TCanvas * cSCEP_EneSharing = new TCanvas ("cSCEP_EneSharing", "cSCEP_EneSharing", 600, 600);
  cSCEP_EneSharing->cd();
  hSCEP_EneT[selEne]->Draw();
  hSCEP_EneT[selEne]->SetStats(0);
  hSCEP_EneT[selEne]->SetTitle(0);  
  hSCEP_EneT[selEne]->SetLineWidth(2);
  hSCEP_EneT[selEne]->GetXaxis()->SetTitle("Energy deposited in ECAL [GeV]");
  hSCEP_EneT[selEne]->GetYaxis()->SetTitle("Counts");
  hSCEP_EneT[selEne]->GetYaxis()->SetRangeUser(10, hSCEP_EneT[selEne]->GetMaximum()*5);
  hSCEP_EneT[selEne]->SetLineColor(kBlack);
  
  hSCEP_EneR[selEne]->SetLineColor(kGreen+1);
  hSCEP_EneR[selEne]->SetLineWidth(2);
  hSCEP_EneR[selEne]->Draw("same");
  
  hSCEP_EneF[selEne]->SetLineColor(kBlue);
  hSCEP_EneF[selEne]->SetLineWidth(2);
  hSCEP_EneF[selEne]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.45,0.88,NULL,"brNDC");  
  leg->AddEntry(hSCEP_EneF[selEne], "Front crystal", "lp");    
  leg->AddEntry(hSCEP_EneR[selEne], "Rear crystal", "lp");    
  leg->AddEntry(hSCEP_EneT[selEne], "Total", "lp");    
  leg->Draw();
  gPad->SetLogy();
  cSCEP_EneSharing->SaveAs("plots/cSCEP_EneSharing.pdf");
  
  
  
  TCanvas * cSCEP_CherSharing = new TCanvas ("cSCEP_CherSharing", "cSCEP_CherSharing", 600, 500);
  cSCEP_CherSharing->cd();
  hSCEP_CherT[selEne]->Draw();
  hSCEP_CherT[selEne]->SetStats(0);
  hSCEP_CherT[selEne]->SetTitle(0);  
  hSCEP_CherT[selEne]->GetXaxis()->SetTitle("Cherenkov signal [nb. of photons]");
  hSCEP_CherT[selEne]->GetYaxis()->SetTitle("Counts");
  hSCEP_CherT[selEne]->SetLineWidth(2);
  hSCEP_CherT[selEne]->SetLineColor(kBlack);
  
  hSCEP_CherR[selEne]->SetLineColor(kGreen+1);
  hSCEP_CherR[selEne]->SetLineWidth(2);
  hSCEP_CherR[selEne]->Draw("same");
  hSCEP_CherF[selEne]->SetLineWidth(2);
  hSCEP_CherF[selEne]->SetLineColor(kBlue);
  hSCEP_CherF[selEne]->Draw("same");
  
  leg = new TLegend(0.15,0.68,0.45,0.88,NULL,"brNDC");  
  leg->AddEntry(hSCEP_CherF[selEne], "Front crystal", "lp");    
  leg->AddEntry(hSCEP_CherR[selEne], "Rear crystal", "lp");    
  leg->AddEntry(hSCEP_CherT[selEne], "Total in SCEPCal", "lp");    
  leg->Draw();
  
  
  
  TCanvas * cSCEP_EneT = new TCanvas ("cSCEP_EneT", "cSCEP_EneT", 600, 500);
  cSCEP_EneT->cd();
  hSCEP_EneT[NFILES-1]->Draw();
  hSCEP_EneT[NFILES-1]->SetStats(0);
  hSCEP_EneT[NFILES-1]->SetTitle(0);  
  hSCEP_EneT[NFILES-1]->GetXaxis()->SetTitle("Energy deposited in ECAL [GeV]");
  hSCEP_EneT[NFILES-1]->GetYaxis()->SetTitle("Counts");
  
  leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      hSCEP_EneT[iFile]->SetLineColor(NFILES-iFile);
      hSCEP_EneT[iFile]->Draw("same");
  }
  leg->Draw();
  cSCEP_EneT->SaveAs(Form("plots/cSCEP_EneT_%s.png", particle_name.c_str()));
  

  
  TCanvas * cScatterCS = new TCanvas ("cScatterCS", "cScatterCS", 600, 600);
  cScatterCS->cd();
  hScatterCS[0]->SetStats(0);
  hScatterCS[0]->SetTitle(0);
  hScatterCS[0]->Draw("COLZ");
  hScatterCS[0]->GetXaxis()->SetTitle("S_{R} [phe/GeV]");
  hScatterCS[0]->GetXaxis()->SetRangeUser(0,25000);
  hScatterCS[0]->GetYaxis()->SetTitle("C_{R} [phe/GeV]");
  hScatterCS[0]->GetYaxis()->SetRangeUser(120,200);
  cScatterCS->SaveAs("plots/cSCEPCcalScatterCS.pdf");
  
  
//   TCanvas * cEne_vs_Phi = new TCanvas ("cEne_vs_Phi", "cEne_vs_Phi", 600, 500);
//   cEne_vs_Phi->cd();
//   hEne_vs_Phi[NFILES-1]->Draw();
//   hEne_vs_Phi[NFILES-1]->SetStats(0);
//   hEne_vs_Phi[NFILES-1]->SetTitle(0);
//   hEne_vs_Phi[NFILES-1]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
//   hEne_vs_Phi[NFILES-1]->GetYaxis()->SetRangeUser(0,maxEne*1.5);
//   hEne_vs_Phi[NFILES-1]->GetXaxis()->SetTitle("#phi [rad]");
//   hEne_vs_Phi[NFILES-1]->GetYaxis()->SetTitle("Energy deposited in ECAL [GeV]");
//   
//   leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");
// 
//   for (int iFile = NFILES-1; iFile>=0; iFile--)
//   {
//       hEne_vs_Phi[iFile]->SetLineColor(NFILES-iFile);
//       hEne_vs_Phi[iFile]->SetMarkerColor(NFILES-iFile);
//       hEne_vs_Phi[iFile]->SetMarkerStyle(20);
//       hEne_vs_Phi[iFile]->Draw("same");
//             
//       TF1 * fitLin = new TF1 ("fitLin", "pol1", -2, 2);
//       hEne_vs_Phi[iFile]->Fit(fitLin, "QR");
// 
//       leg->AddEntry(hEne_vs_Phi[iFile], Form("%d GeV %s", energies[iFile], particle_name.c_str()), "lp");          
//   }
//   leg->Draw();
//   cEne_vs_Phi->SaveAs(Form("plots/cEne_vs_Phi_%s.png", particle_name.c_str()));
//   
//   TCanvas * cEne_vs_Theta = new TCanvas ("cEne_vs_Theta", "cEne_vs_Theta", 600, 500);
//   cEne_vs_Theta->cd();
//   hEne_vs_Theta[NFILES-1]->Draw();
//   hEne_vs_Theta[NFILES-1]->SetStats(0);
//   hEne_vs_Theta[NFILES-1]->SetTitle(0);
//   hEne_vs_Theta[NFILES-1]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
//   hEne_vs_Theta[NFILES-1]->GetYaxis()->SetRangeUser(0,maxEne*1.5);
//   hEne_vs_Theta[NFILES-1]->GetXaxis()->SetTitle("#theta [rad]");
//   hEne_vs_Theta[NFILES-1]->GetYaxis()->SetTitle("Energy deposited in ECAL [GeV]");
//   
//   leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");
// 
//   for (int iFile = NFILES-1; iFile>=0; iFile--)
//   {
//       hEne_vs_Theta[iFile]->SetLineColor(NFILES-iFile);
//       hEne_vs_Theta[iFile]->Draw("same");
//       hEne_vs_Theta[iFile]->SetMarkerColor(NFILES-iFile);
//       hEne_vs_Theta[iFile]->SetMarkerStyle(20);
//             
//       TF1 * fitLin = new TF1 ("fitLin", "pol1",-2, 2);
//       hEne_vs_Theta[iFile]->Fit(fitLin, "QR");
// 
//       leg->AddEntry(hEne_vs_Theta[iFile], Form("%d GeV %s", energies[iFile],particle_name.c_str()), "lp");          
//   }
//   leg->Draw();
//   cEne_vs_Theta->SaveAs(Form("plots/cEne_vs_Theta_%s.png", particle_name.c_str()));


  TCanvas * cEne_vs_Phi = new TCanvas ("cEne_vs_Phi", "cEne_vs_Phi", 600, 500);
  cEne_vs_Phi->cd();
  hEne_vs_Phi[0]->Draw();
  hEne_vs_Phi[0]->SetStats(0);
  hEne_vs_Phi[0]->SetTitle(0);
  hEne_vs_Phi[0]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
//   hEne_vs_Phi[0]->GetYaxis()->SetRangeUser(0,maxEne*1.5);
  hEne_vs_Phi[0]->GetXaxis()->SetTitle("#phi [rad]");
  hEne_vs_Phi[0]->GetYaxis()->SetTitle("Energy deposited in crystals [GeV]");
  
  TCanvas * cEne_vs_Theta = new TCanvas ("cEne_vs_Theta", "cEne_vs_Theta", 600, 500);
  cEne_vs_Theta->cd();
  hEne_vs_Theta[0]->Draw();
  hEne_vs_Theta[0]->SetStats(0);
  hEne_vs_Theta[0]->SetTitle(0);
  hEne_vs_Theta[0]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
//   hEne_vs_Theta[0]->GetYaxis()->SetRangeUser(0,maxEne*1.5);
  hEne_vs_Theta[0]->GetXaxis()->SetTitle("#theta [rad]");
  hEne_vs_Theta[0]->GetYaxis()->SetTitle("Energy deposited in crystals [GeV]");
  
    
  
  TCanvas * cS_vs_Phi = new TCanvas ("cS_vs_Phi", "cS_vs_Phi", 600, 600);
  cS_vs_Phi->cd();
  hS_vs_Phi[0]->Draw("P");
  hS_vs_Phi[0]->SetStats(0);
  hS_vs_Phi[0]->SetTitle(0);
  hS_vs_Phi[0]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
  hS_vs_Phi[0]->GetYaxis()->SetRangeUser(0,3000);
  hS_vs_Phi[0]->GetXaxis()->SetTitle("#phi [rad]");
  hS_vs_Phi[0]->GetYaxis()->SetTitle("S [phe/GeV]");
  
  hS_vs_Phi[0]->SetLineColor(kBlack);
  hS_vs_Phi[0]->SetMarkerColor(kBlack);
  hS_vs_Phi[0]->SetMarkerStyle(20);
  hS_vs_Phi[0]->SetLineWidth(2);
  
    
  hSF_vs_Phi[0]->Draw("same");
  hSF_vs_Phi[0]->SetLineColor(kBlue);
  hSF_vs_Phi[0]->SetMarkerColor(kBlue);
  hSF_vs_Phi[0]->SetMarkerStyle(20);
  hSF_vs_Phi[0]->SetLineWidth(2);
  hSR_vs_Phi[0]->Draw("same");
  hSR_vs_Phi[0]->SetLineColor(kGreen+1);  
  hSR_vs_Phi[0]->SetMarkerColor(kGreen+1);  
  hSR_vs_Phi[0]->SetMarkerStyle(20);
  hSR_vs_Phi[0]->SetLineWidth(2);
   
  leg = new TLegend(0.15,0.68,0.55,0.88,NULL,"brNDC");
  leg->AddEntry(hS_vs_Phi[0], "S (total)", "lp");          
  leg->AddEntry(hSF_vs_Phi[0], "S_{F} (front segment)", "lp");          
  leg->AddEntry(hSR_vs_Phi[0], "S_{R} (rear segment)", "lp");          
  leg->Draw();
  
//   TF1 * fitLinearity = new TF1 ("fitLinearity", "pol0");
//   fitLinearity->SetLineColor(kBlack);
//   hS_vs_Phi[0]->Fit(fitLinearity, "Q");
//   
//   TLine * upLim  = new TLine(-3.14/2, fitLinearity->GetParameter(0)*1.01, 3.14/2, fitLinearity->GetParameter(0)*1.01);
//   TLine * lowLim = new TLine(-3.14/2, fitLinearity->GetParameter(0)*0.99, 3.14/2, fitLinearity->GetParameter(0)*0.99);
//   
//   upLim->SetLineStyle(7);
//   upLim->SetLineColor(kGray+2);
//   upLim->SetLineWidth(2);
//   
//   lowLim->SetLineStyle(7);
//   lowLim->SetLineColor(kGray+1);
//   lowLim->SetLineWidth(2);
//   
//   upLim->Draw();
//   lowLim->Draw();
  
  cS_vs_Phi->SaveAs("plots/cS_vs_Phi.pdf");
  
   
  TCanvas * cS_vs_Theta = new TCanvas ("cS_vs_Theta", "cS_vs_Theta", 600, 600);
  cS_vs_Theta->cd();
  hS_vs_Theta[0]->Draw();
  hS_vs_Theta[0]->SetStats(0);
  hS_vs_Theta[0]->SetTitle(0);
//   hS_vs_Theta[0]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
  hS_vs_Theta[0]->GetYaxis()->SetRangeUser(0,3000);
  hS_vs_Theta[0]->GetXaxis()->SetTitle("#theta [rad]");
  hS_vs_Theta[0]->GetYaxis()->SetTitle("S [phe/GeV]");
  
  hS_vs_Theta[0]->SetLineColor(kBlack);
  hS_vs_Theta[0]->SetMarkerColor(kBlack);
  hS_vs_Theta[0]->SetMarkerStyle(20);
  hS_vs_Theta[0]->SetLineWidth(2);
  
    
  hSF_vs_Theta[0]->Draw("same");
  hSF_vs_Theta[0]->SetLineColor(kBlue);
  hSF_vs_Theta[0]->SetMarkerColor(kBlue);
  hSF_vs_Theta[0]->SetMarkerStyle(20);
  hSF_vs_Theta[0]->SetLineWidth(2);
  hSR_vs_Theta[0]->Draw("same");
  hSR_vs_Theta[0]->SetLineColor(kGreen+1);  
  hSR_vs_Theta[0]->SetMarkerColor(kGreen+1);  
  hSR_vs_Theta[0]->SetMarkerStyle(20);
  hSR_vs_Theta[0]->SetLineWidth(2);
   
  leg = new TLegend(0.15,0.68,0.55,0.88,NULL,"brNDC");
  leg->AddEntry(hS_vs_Theta[0], "S (total)", "lp");          
  leg->AddEntry(hSF_vs_Theta[0], "S_{F} (front segment)", "lp");          
  leg->AddEntry(hSR_vs_Theta[0], "S_{R} (rear segment)", "lp");          
  leg->Draw();
  cS_vs_Theta->SaveAs("plots/cS_vs_Theta.pdf");
  
  
  
    TCanvas * cC_vs_Phi = new TCanvas ("cC_vs_Phi", "cC_vs_Phi", 600, 600);
  cC_vs_Phi->cd();
  hC_vs_Phi[0]->Draw();
  hC_vs_Phi[0]->SetStats(0);
  hC_vs_Phi[0]->SetTitle(0);
  hC_vs_Phi[0]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
  hC_vs_Phi[0]->GetYaxis()->SetRangeUser(0,260);
  hC_vs_Phi[0]->GetXaxis()->SetTitle("#phi [rad]");
  hC_vs_Phi[0]->GetYaxis()->SetTitle("C [phe/GeV]");
  
  hC_vs_Phi[0]->SetLineColor(kBlack);
  hC_vs_Phi[0]->SetMarkerColor(kBlack);
  hC_vs_Phi[0]->SetMarkerStyle(20);
  hC_vs_Phi[0]->SetLineWidth(2);
  
  
  cC_vs_Phi->SaveAs("plots/cC_vs_Phi.pdf");
  
   
  TCanvas * cC_vs_Theta = new TCanvas ("cC_vs_Theta", "cC_vs_Theta", 600, 600);
  cC_vs_Theta->cd();
  hC_vs_Theta[0]->Draw();
  hC_vs_Theta[0]->SetStats(0);
  hC_vs_Theta[0]->SetTitle(0);
//   hC_vs_Theta[0]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
  hC_vs_Theta[0]->GetYaxis()->SetRangeUser(0,260);
  hC_vs_Theta[0]->GetXaxis()->SetTitle("#theta [rad]");
  hC_vs_Theta[0]->GetYaxis()->SetTitle("C [phe/GeV]");
  
  hC_vs_Theta[0]->SetLineColor(kBlack);
  hC_vs_Theta[0]->SetMarkerColor(kBlack);
  hC_vs_Theta[0]->SetMarkerStyle(20);
  hC_vs_Theta[0]->SetLineWidth(2);
  
  cC_vs_Theta->SaveAs("plots/cC_vs_Theta.pdf");
  
  theApp->Run();
}

