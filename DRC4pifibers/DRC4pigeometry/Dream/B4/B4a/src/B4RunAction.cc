//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4RunAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include "stdlib.h"
using namespace std;
#include "B4aEventAction.hh"
#include <sstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
 /* G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  const B4aEventAction* constEventAction = static_cast<const B4aEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
  B4aEventAction* eventAction = const_cast<B4aEventAction*>(constEventAction);


  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetFirstHistoId(0);

  // Book histograms, ntuple
  //
  
  // Creating histograms
  //analysisManager->CreateH1("1","Edep in module", 100, 0., 1000*MeV);
  //analysisManager->CreateH1("2","trackL in module", 100, 0., 1*m);
  //analysisManager->CreateH1("3", "Edep in scintillating fibers", 100, 0., 1000*MeV);
  // Creating ntuple

  //G4double vector[10]={0};
  //vector<int>* interactions = new vector<int>();
  //vector<G4double> vector[10]={0.};
  //
  analysisManager->CreateNtuple("B4", "Edep");
  analysisManager->CreateNtupleDColumn("Eem");
  analysisManager->CreateNtupleDColumn("EScin");
  analysisManager->CreateNtupleDColumn("ECher");
  analysisManager->CreateNtupleDColumn("Cherenkovintheglass");
  analysisManager->CreateNtupleDColumn("Scinintheglass");
  analysisManager->CreateNtupleDColumn("EnergyTot");
  analysisManager->CreateNtupleDColumn("EnergyInside2mm");
  analysisManager->CreateNtupleIColumn("prova",eventAction->GetVectorSignals());
  analysisManager->FinishNtuple();*/


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void B4RunAction::BeginOfRunAction(const G4Run* run, G4string outputFileName="B4")
void B4RunAction::BeginOfRunAction(const G4Run* run)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  G4cout<< "run " <<run->GetRunID()<<G4endl;

  const B4aEventAction* constEventAction = static_cast<const B4aEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
  B4aEventAction* eventAction = const_cast<B4aEventAction*>(constEventAction);
            
  std::stringstream ss;
  std::string myrun;
  ss<<run->GetRunID();
  ss>>myrun;
  
  
  // Creating ntuple
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;
  
//   Open an output file
//   G4String fileName = outputFileName;    
  G4String fileName = "out.root";    
  analysisManager->OpenFile(fileName);    
  
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetFirstHistoId(0);
  
  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");              
  // Creating histograms
  //analysisManager->CreateH1("1","Edep in module", 100, 0., 1000*MeV);       
  
  analysisManager->CreateNtuple("B4", "edep");
  analysisManager->CreateNtupleDColumn("Energyem");
  analysisManager->CreateNtupleDColumn("EnergyScin");
  analysisManager->CreateNtupleDColumn("EnergyCher");
  analysisManager->CreateNtupleDColumn("NofCherenkovDetected"); 
  analysisManager->CreateNtupleDColumn("EnergyTot");
  analysisManager->CreateNtupleDColumn("PrimaryParticleEnergy");
  analysisManager->CreateNtupleSColumn("PrimaryParticleName");  
  analysisManager->CreateNtupleDColumn("neutrinoleakage");
  analysisManager->CreateNtupleDColumn("leakage");
   
  analysisManager->CreateNtupleDColumn("SCEP_EnergyDepF");
  analysisManager->CreateNtupleDColumn("SCEP_NCherProdF");
  analysisManager->CreateNtupleDColumn("SCEP_EnergyDepR");
  analysisManager->CreateNtupleDColumn("SCEP_NCherProdR");
   
  analysisManager->CreateNtupleDColumn("SCEP_Timing_EnergyDepF");
  analysisManager->CreateNtupleDColumn("SCEP_Timing_EnergyDepR");
  analysisManager->CreateNtupleDColumn("SCEP_Timing_TimeF");
  analysisManager->CreateNtupleDColumn("SCEP_Timing_TimeR");
  
  analysisManager->CreateNtupleDColumn("VecHit_CrystalID",  eventAction->GetVecScep_CrystalID());
  analysisManager->CreateNtupleDColumn("VecHit_ScepEneDepF",eventAction->GetVecScep_ScepEneDepF());
  analysisManager->CreateNtupleDColumn("VecHit_ScepEneDepR",eventAction->GetVecScep_ScepEneDepR());
  analysisManager->CreateNtupleDColumn("VecHit_ScepCherF",  eventAction->GetVecScep_ScepCherF());
  analysisManager->CreateNtupleDColumn("VecHit_ScepCherR",  eventAction->GetVecScep_ScepCherR());
  
  analysisManager->CreateNtupleDColumn("VecHit_Timing_CrystalID_F",  eventAction->GetVecScep_Timing_CrystalID_F());
  analysisManager->CreateNtupleDColumn("VecHit_Timing_CrystalID_R",  eventAction->GetVecScep_Timing_CrystalID_R());
  analysisManager->CreateNtupleDColumn("VecHit_Timing_ScepEneDepF",eventAction->GetVecScep_Timing_ScepEneDepF());
  analysisManager->CreateNtupleDColumn("VecHit_Timing_ScepEneDepR",eventAction->GetVecScep_Timing_ScepEneDepR());
  analysisManager->CreateNtupleDColumn("VecHit_Timing_ScepTimeF",  eventAction->GetVecScep_Timing_ScepTimeF());
  analysisManager->CreateNtupleDColumn("VecHit_Timing_ScepTimeR",  eventAction->GetVecScep_Timing_ScepTimeR());
  
  analysisManager->CreateNtupleDColumn("PrimaryParticleMomentum",  eventAction->GetPrimaryParticleMomentum());
    
  analysisManager->CreateNtupleDColumn("VectorSignalsR",eventAction->GetVectorSignalsR());
  analysisManager->CreateNtupleDColumn("VectorSignalsL",eventAction->GetVectorSignalsL());
  analysisManager->CreateNtupleDColumn("VectorSignalsCherR",eventAction->GetVectorSignalsCherR());
  analysisManager->CreateNtupleDColumn("VectorSignalsCherL",eventAction->GetVectorSignalsCherL());
  analysisManager->CreateNtupleDColumn("VectorL", eventAction->GetVectorL());
  analysisManager->CreateNtupleDColumn("VectorR", eventAction->GetVectorR());
  analysisManager->CreateNtupleDColumn("VectorL_loop", eventAction->GetVectorL_loop());
  analysisManager->CreateNtupleDColumn("VectorR_loop", eventAction->GetVectorR_loop()); 
  analysisManager->FinishNtuple();                
    
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
