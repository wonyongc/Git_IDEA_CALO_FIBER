// g++ -Wall -o plotDualReadout4Pi plotDualReadout4Pi.C  VectorSmallestInterval.cc myG4Tree.cc myG4Tree.hh `root-config --cflags --glibs`


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
  int energies[NFILES] = {10};
  
//   const int NFILES = 6;
//   int energies[NFILES] = {1, 5, 10, 30, 60, 120};
  
//   int energies[NFILES] = {5, 10, 15, 30, 60, 120, 300};//, 300};  
  
  std::cout << " NFILES = " << NFILES << std::endl;
  
  std::string particle_name = "#gamma";
  //define vector energy
  std::vector<double> *vEcalEnergy[NFILES];
  float nCherPerGeV = 3000*4;
  
  float LO       = 2000;
  float DCR      = 1; //GHz -- dark counts per nanosecond - -- for 10 um , 0.1 MHz/mm², assume a factor 100 because of radiation and a surface of 100 mm² --> 1000 MHz --> 1GHz
  float int_gate = 100; //integration gate in ns 
  float noise    = sqrt(DCR*int_gate); //phe noise (fluctuations)
  
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
  TProfile* hEne_vs_Eta[NFILES];
  
  int NBIN_ENE = 1000;
  double maxEne = 120.5;  //for MIPs
  
  for (int iFile = 0; iFile<NFILES; iFile++)
  {
      
      std::cout << "energies[" << iFile << "] = " << energies[iFile] << std::endl;
      vEcalEnergy[iFile]  = new std::vector<double>();
      
      double minEne = 0;
      minEne = energies[iFile]*0.7;
      maxEne = energies[iFile]*1.2;
      
      
      hSCEP_EneF[iFile] = new TH1F (Form("hSCEP_EneF_%d", energies[iFile]), Form("hSCEP_EneF_%d", energies[iFile]), NBIN_ENE, minEne, maxEne);
      hSCEP_EneR[iFile] = new TH1F (Form("hSCEP_EneR_%d", energies[iFile]), Form("hSCEP_EneR_%d", energies[iFile]), NBIN_ENE, minEne, maxEne);
      hSCEP_EneT[iFile] = new TH1F (Form("hSCEP_EneT_%d", energies[iFile]), Form("hSCEP_EneT_%d", energies[iFile]), NBIN_ENE, minEne, maxEne);
      
      hSCEP_CherF[iFile] = new TH1F (Form("hSCEP_CherF_%d", energies[iFile]), Form("hSCEP_CherF_%d", energies[iFile]), NBIN_ENE, 0, energies[iFile]*1.5*nCherPerGeV);
      hSCEP_CherR[iFile] = new TH1F (Form("hSCEP_CherR_%d", energies[iFile]), Form("hSCEP_CherR_%d", energies[iFile]), NBIN_ENE, 0, energies[iFile]*1.5*nCherPerGeV);
      hSCEP_CherT[iFile] = new TH1F (Form("hSCEP_CherT_%d", energies[iFile]), Form("hSCEP_CherT_%d", energies[iFile]), NBIN_ENE, 0, energies[iFile]*1.5*nCherPerGeV);
      
      hScatterCS[iFile] = new TH2F (Form("hScatterCS_%d", energies[iFile]), Form("hScatterCS_%d", energies[iFile]), NBIN_ENE, minEne, maxEne, NBIN_ENE, 0, energies[iFile]*1.5*nCherPerGeV);
      
      hEne_vs_Phi[iFile] = new TProfile (Form("hEne_vs_Phi_%d", energies[iFile]), Form("hEne_vs_Phi_%d", energies[iFile]), 50, -2, 2);
      hEne_vs_Theta[iFile] = new TProfile (Form("hEne_vs_Theta_%d", energies[iFile]), Form("hEne_vs_Theta_%d", energies[iFile]), 50, -2, 2);
      hEne_vs_Eta[iFile] = new TProfile (Form("hEne_vs_Eta_%d", energies[iFile]), Form("hEne_vs_Eta_%d", energies[iFile]), 50, -4, 4);
  }
  
  
  
  TGraphErrors * gSCEP_EneRes      = new TGraphErrors();
  TGraphErrors * gSCEP_EneResPhot  = new TGraphErrors();
  TGraphErrors * gSCEP_EneResNoise = new TGraphErrors();
  TGraphErrors * gSCEP_EneResTot   = new TGraphErrors();
  TGraphErrors * gSCEP_EneLin      = new TGraphErrors();
  
  
  
  
  
  
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
        RunFile[iFile] = new TFile("../root_files/merged/ele_10GeV.root","READ"); 
      
      
      TTree* TreeRun = (TTree*) RunFile[iFile]->Get("B4");
      myG4TreeVars myTV;
      InitG4Tree (TreeRun, myTV);
      
      
      

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
          
          hSCEP_CherF[iFile]->Fill(cherF);
          hSCEP_CherR[iFile]->Fill(cherR);
          hSCEP_CherT[iFile]->Fill(cherR+cherF);
          hScatterCS[iFile]->Fill(ene_tot , cherF+cherR);
          
//           double px  = myTV.PrimaryParticleMomentum->at(0);
//           double py  = myTV.PrimaryParticleMomentum->at(1);
//           double pz  = myTV.PrimaryParticleMomentum->at(2);
//           double P   = sqrt(px*px+py*py+pz*pz);
//           px/= P;
//           py/= P;
//           pz/= P;
//           
//           double theta = atan(px/pz);
//           double phi   = atan(py/px);
//           double eta   = atanh(pz);//-log(tan(theta/2));
//           
// //           std::cout << " mom(" << px << ", " << py << ", " << pz << ") :: P = " << P << " :: theta = " << theta << " :: phi = " << phi << " :: eta = " << eta << std::endl;
//           
//           //if in barrel
// //           if (abs(theta)<3.14/2/2)             
//           if (ene_tot > energies[iFile]*0.5 )
//           {
//             hEne_vs_Phi[iFile]->Fill(phi, ene_tot);
//             hEne_vs_Theta[iFile]->Fill(theta, ene_tot);
//           
//             //barrel
//             hEne_vs_Eta[iFile]->Fill(eta, ene_tot);
//           }
//           //endcap
//           
          
      }
      
      
      
  }
  
  
  int selEne = 0;
  TCanvas * cSCEP_EneSharing = new TCanvas ("cSCEP_EneSharing", "cSCEP_EneSharing", 600, 500);
  cSCEP_EneSharing->cd();
  hSCEP_EneT[selEne]->Draw();
  hSCEP_EneT[selEne]->SetStats(0);
  hSCEP_EneT[selEne]->SetTitle(0);  
  hSCEP_EneT[selEne]->GetXaxis()->SetTitle("Energy deposited in ECAL [GeV]");
  hSCEP_EneT[selEne]->GetYaxis()->SetTitle("Counts");
  
  hSCEP_EneR[selEne]->SetLineColor(kBlue+1);
  hSCEP_EneR[selEne]->Draw("same");
  
  hSCEP_EneF[selEne]->SetLineColor(kGreen+1);
  hSCEP_EneF[selEne]->Draw("same");
  
  leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");  
  leg->AddEntry(hSCEP_EneF[selEne], "Front crystal", "lp");    
  leg->AddEntry(hSCEP_EneR[selEne], "Rear crystal", "lp");    
  leg->AddEntry(hSCEP_EneT[selEne], "Total in SCEPCal", "lp");    
  leg->Draw();
  
  
  
  TCanvas * cSCEP_CherSharing = new TCanvas ("cSCEP_CherSharing", "cSCEP_CherSharing", 600, 500);
  cSCEP_CherSharing->cd();
  hSCEP_CherT[selEne]->Draw();
  hSCEP_CherT[selEne]->SetStats(0);
  hSCEP_CherT[selEne]->SetTitle(0);  
  hSCEP_CherT[selEne]->GetXaxis()->SetTitle("Cherenkov signal [nb. of photons]");
  hSCEP_CherT[selEne]->GetYaxis()->SetTitle("Counts");
  
  hSCEP_CherR[selEne]->SetLineColor(kBlue+1);
  hSCEP_CherR[selEne]->Draw("same");
  
  hSCEP_CherF[selEne]->SetLineColor(kGreen+1);
  hSCEP_CherF[selEne]->Draw("same");
  
  leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");  
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
      
      
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", 0.98, 1.05);//maxEnHCAL[iFile]);
      fitGaus->SetParameters(hSCEP_EneT[iFile]->GetMaximum(), hSCEP_EneT[iFile]->GetMean(), hSCEP_EneT[iFile]->GetRMS());
      fitGaus->SetLineColor(NFILES-iFile);
//       TF1 * fitGaus = new TF1 ("fitGaus", "landau", 0, (maxEnECAL+maxEnHCAL)*2);
      hSCEP_EneT[iFile]->Fit(fitGaus, "QR");

      float sigma_eff = FindSmallestInterval(vEcalEnergy[iFile], 0.68, false) /2. ;
//       float mean      = hSCEP_EneT[iFile]->GetMean();
      float mean      = hSCEP_EneT[iFile]->GetBinCenter(hSCEP_EneT[iFile]->GetMaximumBin());
      gSCEP_EneRes     ->SetPoint (iFile, energies[iFile], sigma_eff/mean *100);      
      
              
      float phot_stoch = 1./sqrt(LO*mean)*100;
      gSCEP_EneResPhot->SetPoint(iFile, energies[iFile], phot_stoch);    
      
      float noise_term = noise/(LO*mean )*100;
      gSCEP_EneResNoise->SetPoint(iFile, energies[iFile], noise_term ) ;   
      
      gSCEP_EneResTot->SetPoint(iFile, energies[iFile], sqrt(pow(sigma_eff/mean*100,2) + pow(phot_stoch ,2) + pow(noise_term ,2))  ) ;   
        
      gSCEP_EneLin->SetPoint (iFile, energies[iFile], hSCEP_EneT[iFile]->GetMean());
      
      leg->AddEntry(hSCEP_EneT[iFile], Form("%d GeV %s", energies[iFile], particle_name.c_str()), "lp");    
      
  }
  leg->Draw();
  cSCEP_EneT->SaveAs(Form("plots/cSCEP_EneT_%s.png", particle_name.c_str()));
  
  TCanvas * cSCEP_Res_vsEne = new TCanvas ("cSCEP_Res_vsEne", "cSCEP_Res_vsEne", 600, 500);
  cSCEP_Res_vsEne->cd();
  gPad->SetLogx();
  gPad->SetLogy();
  gSCEP_EneResTot->SetMarkerStyle(20);
  gSCEP_EneResTot->GetXaxis()->SetLimits(0.5, 130);
  gSCEP_EneResTot->SetMinimum(1e-2);
  gSCEP_EneResTot->SetMaximum(5);
  gSCEP_EneResTot->GetXaxis()->SetTitle("Energy [GeV]");
  gSCEP_EneResTot->GetYaxis()->SetTitle("Energy resolution [%]");
  gSCEP_EneResTot->Draw("ALPE");
  
  gSCEP_EneRes->SetLineWidth(2);
  gSCEP_EneRes->SetLineColor(kBlue+1);
  gSCEP_EneRes->SetMarkerColor(kBlue+1);
  gSCEP_EneRes->Draw("same lpe");
  
  gSCEP_EneResPhot->SetLineWidth(2);
  gSCEP_EneResPhot->SetLineColor(kRed+1);
  gSCEP_EneResPhot->SetMarkerColor(kRed+1);
  gSCEP_EneResPhot->Draw("same lpe");
  
  gSCEP_EneResNoise->SetLineWidth(2);
  gSCEP_EneResNoise->SetLineColor(kGreen+1);
  gSCEP_EneResNoise->SetMarkerColor(kGreen+1);
  gSCEP_EneResNoise->Draw("same lpe");
  
  leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");
  leg->AddEntry(gSCEP_EneResTot, "total", "lp");          
  leg->AddEntry(gSCEP_EneResPhot, "photostatistics", "lp");          
  leg->AddEntry(gSCEP_EneResNoise,"noise", "lp");          
  leg->AddEntry(gSCEP_EneRes,     "shower fluctuations", "lp");          
  
  leg->Draw();
  
  TCanvas *  cSCEP_Lin_vsEne = new TCanvas ("cSCEP_Lin_vsEne", "cSCEP_Lin_vsEne", 600, 500);
  cSCEP_Lin_vsEne->cd();
  gSCEP_EneLin->SetMarkerStyle(20);
  gSCEP_EneLin->GetXaxis()->SetLimits(0, 100);
  gSCEP_EneLin->SetMinimum(0.95);
  gSCEP_EneLin->SetMaximum(1.05);
  gSCEP_EneLin->GetXaxis()->SetTitle("Energy [GeV]");
  gSCEP_EneLin->GetYaxis()->SetTitle("Energy linearity");
  gSCEP_EneLin->Draw("ALPE");
  
  
  TCanvas * cScatterCS = new TCanvas ("cScatterCS", "cScatterCS", 600, 500);
  cScatterCS->cd();
  hScatterCS[0]->Draw("COLZ");
  hScatterCS[0]->GetXaxis()->SetTitle("Scintillation signal [GeV]");
  hScatterCS[0]->GetYaxis()->SetTitle("Cherenkov signal [nb. of photons]");
  
  
  
  TCanvas * cEne_vs_Phi = new TCanvas ("cEne_vs_Phi", "cEne_vs_Phi", 600, 500);
  cEne_vs_Phi->cd();
  hEne_vs_Phi[NFILES-1]->Draw();
  hEne_vs_Phi[NFILES-1]->SetStats(0);
  hEne_vs_Phi[NFILES-1]->SetTitle(0);
  hEne_vs_Phi[NFILES-1]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
  hEne_vs_Phi[NFILES-1]->GetYaxis()->SetRangeUser(0,maxEne*1.5);
  hEne_vs_Phi[NFILES-1]->GetXaxis()->SetTitle("#phi [rad]");
  hEne_vs_Phi[NFILES-1]->GetYaxis()->SetTitle("Energy deposited in ECAL [GeV]");
  
  leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      hEne_vs_Phi[iFile]->SetLineColor(NFILES-iFile);
      hEne_vs_Phi[iFile]->SetMarkerColor(NFILES-iFile);
      hEne_vs_Phi[iFile]->SetMarkerStyle(20);
      hEne_vs_Phi[iFile]->Draw("same");
            
      TF1 * fitLin = new TF1 ("fitLin", "pol1", -2, 2);
      hEne_vs_Phi[iFile]->Fit(fitLin, "QR");

      leg->AddEntry(hEne_vs_Phi[iFile], Form("%d GeV %s", energies[iFile], particle_name.c_str()), "lp");          
  }
  leg->Draw();
  cEne_vs_Phi->SaveAs(Form("plots/cEne_vs_Phi_%s.png", particle_name.c_str()));
  
  TCanvas * cEne_vs_Theta = new TCanvas ("cEne_vs_Theta", "cEne_vs_Theta", 600, 500);
  cEne_vs_Theta->cd();
  hEne_vs_Theta[NFILES-1]->Draw();
  hEne_vs_Theta[NFILES-1]->SetStats(0);
  hEne_vs_Theta[NFILES-1]->SetTitle(0);
  hEne_vs_Theta[NFILES-1]->GetXaxis()->SetRangeUser(-3.14/2,3.14/2);
  hEne_vs_Theta[NFILES-1]->GetYaxis()->SetRangeUser(0,maxEne*1.5);
  hEne_vs_Theta[NFILES-1]->GetXaxis()->SetTitle("#theta [rad]");
  hEne_vs_Theta[NFILES-1]->GetYaxis()->SetTitle("Energy deposited in ECAL [GeV]");
  
  leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      hEne_vs_Theta[iFile]->SetLineColor(NFILES-iFile);
      hEne_vs_Theta[iFile]->Draw("same");
      hEne_vs_Theta[iFile]->SetMarkerColor(NFILES-iFile);
      hEne_vs_Theta[iFile]->SetMarkerStyle(20);
            
      TF1 * fitLin = new TF1 ("fitLin", "pol1",-2, 2);
      hEne_vs_Theta[iFile]->Fit(fitLin, "QR");

      leg->AddEntry(hEne_vs_Theta[iFile], Form("%d GeV %s", energies[iFile],particle_name.c_str()), "lp");          
  }
  leg->Draw();
  cEne_vs_Theta->SaveAs(Form("plots/cEne_vs_Theta_%s.png", particle_name.c_str()));
  
  TCanvas * cEne_vs_Eta = new TCanvas ("cEne_vs_Eta", "cEne_vs_Eta", 600, 500);
  cEne_vs_Eta->cd();
  hEne_vs_Eta[NFILES-1]->Draw();
  hEne_vs_Eta[NFILES-1]->SetStats(0);
  hEne_vs_Eta[NFILES-1]->SetTitle(0);
  hEne_vs_Eta[NFILES-1]->GetXaxis()->SetRangeUser(-3.5,3.5);
  hEne_vs_Eta[NFILES-1]->GetYaxis()->SetRangeUser(0,maxEne*1.5);
  hEne_vs_Eta[NFILES-1]->GetXaxis()->SetTitle("#eta [rad]");
  hEne_vs_Eta[NFILES-1]->GetYaxis()->SetTitle("Energy deposited in ECAL [GeV]");
  
  leg = new TLegend(0.7,0.68,0.88,0.88,NULL,"brNDC");

  for (int iFile = NFILES-1; iFile>=0; iFile--)
  {
      hEne_vs_Eta[iFile]->SetLineColor(NFILES-iFile);
      hEne_vs_Eta[iFile]->Draw("same");
      hEne_vs_Eta[iFile]->SetMarkerColor(NFILES-iFile);
      hEne_vs_Eta[iFile]->SetMarkerStyle(20);
            
      TF1 * fitLin = new TF1 ("fitLin", "pol0", -3, 3);
      hEne_vs_Eta[iFile]->Fit(fitLin, "QR");

      leg->AddEntry(hEne_vs_Eta[iFile], Form("%d GeV %s", energies[iFile], particle_name.c_str()), "lp");          
  }
  leg->Draw();
  cEne_vs_Eta->SaveAs(Form("plots/cEne_vs_Eta_%s.png", particle_name.c_str()));
  
  
//     TFile * outFile = new TFile (Form("./outputs/wScepcal_btod_%s_%s_%s.root", config.c_str(), particle_name.c_str(), btod.c_str()),"RECREATE");
//     outFile->cd();
//     
//     gHCALRes_vs_energy->SetName("gHCALRes_vs_energy");
//     gHCALRes_vs_energy->Write();
//     gHCALCherRes_vs_energy->SetName("gHCALCherRes_vs_energy");
//     gHCALCherRes_vs_energy->Write();
//     gHCALRes_vs_energy_DRO->SetName("gHCALRes_vs_energy_DRO");
//     gHCALRes_vs_energy_DRO->Write();    
//     gHCALCherRes_vs_energy_DRO->SetName("gHCALCherRes_vs_energy_DRO");
//     gHCALCherRes_vs_energy_DRO->Write();
//     
//     gTotRes_vs_energy->SetName("gTotRes_vs_energy");
//     gTotRes_vs_energy->Write();    
//     gTotResCorr_vs_energy->SetName("gTotResCorr_vs_energy");
//     gTotResCorr_vs_energy->Write();
//     
//     gTotRes_dro_recal_only_vs_energy->SetName("gTotRes_dro_recal_only_vs_energy");
//     gTotRes_dro_recal_only_vs_energy->Write();    
//     gTotRes_dro_hcal_only_vs_energy->SetName("gTotRes_dro_hcal_only_vs_energy");
//     gTotRes_dro_hcal_only_vs_energy->Write();    
//     
//     hScatterDRO_HCAL_MasterC->Write();    
//     hScatterDRO_HCAL_Master->Write();
//     
//     hScatterDRO_ECAL_Master->Write();
//     hScatterDRO_ECAL_MasterCher->Write();
//     hScatterDRO_ECAL_MasterCorr->Write();
//     hScatterDRO_ECAL_MasterCherCorr->Write();
//     
//     pDRO_corr_ECAL_Master->Write();
//     pDRO_corr_ECAL_MasterCher->Write();
//     pDRO_corr_ECAL_MasterCorr->Write();
//     pDRO_corr_ECAL_MasterCherCorr->Write();
//         
//     
//     hTestScintSmearing->Write();
//     hTestCherSmearing->Write();
//     
//     hTestSinC_frac->Write();
//     hTestCinS_frac->Write();
//                 
//     for (int iFile = 0; iFile<NFILES; iFile++)
//     {
//       hRatioContPureS[iFile]->Write();
//       hRatioContPureC[iFile]->Write();
//     }
//   
//     
//     outFile->Write();
//     outFile->Close();
// 
//     
    
   
   theApp->Run();
}

