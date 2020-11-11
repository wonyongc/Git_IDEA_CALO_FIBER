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
// $Id: B4aEventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file B4aEventAction.cc
/// \brief Implementation of the B4aEventAction class

#include "B4aEventAction.hh"
#include "B4RunAction.hh"
#include "B4Analysis.hh"
#include "B4MyMaterials.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"


#include "Randomize.hh"
#include <iomanip>
#include <vector>
#include <map>

#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction()
 : G4UserEventAction(),
   Energyem(0.),
   neutrinoleakage(0.),
   leakage(0.),
   EnergyScin(0.),
   EnergyCher(0.),
   NofCherenkovDetected(0),
   //NofScintillationDetected(0),
   EnergyTot(0.),
   
   SCEP_EnergyDepF(0.),
   SCEP_NCherProdF(0.),
   SCEP_EnergyDepR(0.),
   SCEP_NCherProdR(0.),
   
   VecHit_CrystalID(0.),
   VecHit_ScepEneDepF(0.),
   VecHit_ScepEneDepR(0.),
   VecHit_ScepCherF(0.),
   VecHit_ScepCherR(0.),


   PrimaryParticleEnergy(0.),
   PrimaryParticleMomentum(0.),
   VectorSignalsR(0.),
   VectorSignalsL(0.),
   VectorSignalsCherR(0.),
   VectorSignalsCherL(0.),
   VectorR(0.),
   VectorL(0.),
   VectorR_loop(0.),
   VectorL_loop(0.),
   Fiber_Hits{0.},
   Tracking_Hits{0.}
   
   
   
   
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* event)
{  
  //Time_distribution event
  std::ofstream TimeFile;
  TimeFile.open("Time.txt", std::ios_base::app);
  TimeFile<<"Event "<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()<<" % % %"<<G4endl;
  TimeFile.close();
  
 
  
  
	
  // initialisation per event
  Energyem = 0.;
  EnergyScin = 0.;
  EnergyCher = 0.;
  NofCherenkovDetected = 0;
  //NofScintillationDetected = 0;
  EnergyTot = 0.;
  neutrinoleakage = 0.;
  leakage = 0.;
  
  SCEP_EnergyDepF = 0.;
  SCEP_NCherProdF = 0.;
  SCEP_EnergyDepR = 0.;
  SCEP_NCherProdR = 0.;
  
  VecHit_CrystalID.clear();
  VecHit_ScepEneDepF.clear();
  VecHit_ScepEneDepR.clear();
  VecHit_ScepCherF.clear();
  VecHit_ScepCherR.clear();
  PrimaryParticleMomentum.clear();
  
  
  int fNbOfBarrel = 40;
  int fNbOfEndcap = 35;
  int fNbOfZRot = 36;
  /*for(int i=0;i<64;i++){
   Signalfibre[i]=0;
  }*///only if you want to use SignalFibre[64]
	
    for (int i=0;i<10000;i++)
    Fiber_Hits[i]={0.};
	
    for (int i=0;i<200;i++)
    Tracking_Hits[i]={0.};

    for (int i=0;i<VectorR.size();i++){
    VectorR.at(i)=0.;
    }
    for (int i=0;i<VectorL.size();i++){
    VectorL.at(i)=0.;
    }
	
	for (int i=0;i<VectorR_loop.size();i++){
    VectorR_loop.at(i)=0.;
    }
    for (int i=0;i<VectorL_loop.size();i++){
    VectorL_loop.at(i)=0.;
    }
       
    for (int i=0;i<VectorSignalsR.size();i++){
    VectorSignalsR.at(i)=0.;
    }
    for (int i=0;i<VectorSignalsL.size();i++){
    VectorSignalsL.at(i)=0.;
    }
    for (int i=0;i<VectorSignalsCherR.size();i++){
    VectorSignalsCherR.at(i)=0.;
    }
    for (int i=0;i<VectorSignalsCherL.size();i++){
    VectorSignalsCherL.at(i)=0.;
    }
    PrimaryParticleEnergy = 0;  
    
    
    for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
        if(VectorR.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorR.push_back(0.);}}
    for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
        if(VectorL.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorL.push_back(0.);}}
    
    for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
        if(VectorR_loop.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorR_loop.push_back(0.);}}
    for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
        if(VectorL_loop.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorL_loop.push_back(0.);}}

    for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
        if(VectorSignalsR.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorSignalsR.push_back(0.);}}
    for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
        if(VectorSignalsL.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorSignalsL.push_back(0.);}}
    //VectorSignals.at(i)=0;}
    for(int k=0;k<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);k++){
        if(VectorSignalsCherR.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorSignalsCherR.push_back(0.);}}
    for(int k=0;k<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);k++){
        if(VectorSignalsCherL.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorSignalsCherL.push_back(0.);}}
    
    
  
  // save info on primary particle  
  G4PrimaryVertex* vertex = event->GetPrimaryVertex();
  G4double x = vertex -> GetX0();
  G4double y = vertex -> GetY0();
  G4double z = vertex -> GetZ0();
  
  G4PrimaryParticle * particle = vertex -> GetPrimary();
  G4double InitEnergy = particle -> GetTotalEnergy();
  G4double px = particle -> GetPx();
  G4double py = particle -> GetPy();
  G4double pz = particle -> GetPz();
  
  SavePrimaryParticle(particle->GetParticleDefinition()->GetParticleName());
  SavePrimaryEnergy(InitEnergy);
  SavePrimaryMomentum(px, py, pz);
  
//   std::cout << "ene = " << InitEnergy << ", momentum = (" << px << ", " << py << ", " << pz << std::endl;
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
	
    std::ofstream eventFile;
    eventFile.open("Event.txt", std::ios_base::app);
    /*G4cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
    G4cout<<"\t ID \t Energy(MeV) \t S/C \t Position \t slice \t tower"<<std::endl;
    G4cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;*/
    int v=0;
    G4double E=0.;
    if(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()==0)	eventFile<<"EvtID\tFiberID\tEt\tXt\tYt\tZt\tFlagt\tslicet\ttowert"<<std::endl;
    while(Fiber_Hits[v].F_ID!=0)
    {
    //G4cout<<Fiber_Hits[v].F_ID<<"\t"<<Fiber_Hits[v].F_E<<"\t"<<Fiber_Hits[v].F_Type<<"\t"<<Fiber_Hits[v].F_X<<" "<<Fiber_Hits[v].F_Y<<" "<<Fiber_Hits[v].F_Z<<std::endl;
        E = E+Fiber_Hits[v].F_E;
        eventFile<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()<<"\t"<<std::fixed << std::setprecision(3) <<Fiber_Hits[v].F_ID<<"\t"<<Fiber_Hits[v].F_E<<"\t"<<Fiber_Hits[v].F_X<<"\t"<<Fiber_Hits[v].F_Y<<"\t"<<Fiber_Hits[v].F_Z<<"\t"<<Fiber_Hits[v].F_Type<<"\t"<<Fiber_Hits[v].F_slice<<"\t"<<Fiber_Hits[v].F_tower<<std::endl;
        v++;        
    }
    //eventFile<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
    eventFile.close();
	
    std::ofstream eventFile1;
    eventFile1.open("Event_Track.txt", std::ios_base::app);
    if(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()==0)	eventFile1<<"EvtIDtrack\tTrackID\tXtrackt\tYtrackt\tZtrackt\tparticlesnamet"<<std::endl;
    for (v=0;v<200;v++)
    {
        if(Tracking_Hits[v].T_ID!=0)
        {
            eventFile1<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()<<"\t"<<Tracking_Hits[v].T_ID<<"\t"<<Tracking_Hits[v].T_X<<"\t"<<Tracking_Hits[v].T_Y<<"\t"<<Tracking_Hits[v].T_Z<<"\t"<<Tracking_Hits[v].T_Name<<std::endl;
	}        
    }
    eventFile1.close();
	//eventFile<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
	
	/*G4cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
	G4cout<<"N fired: "<<v<<"\t tot E: "<<E<<"\t Lor: "<<NofCherenkovDetected<<std::endl;*/
	
	
  // fill histograms
  //analysisManager->FillH1(1, Energymodule);
  //analysisManager->FillH1(2, TrackLmodule);
  //analysisManager->FillH1(3, EnergyScin);
  
  std::cout << "filling ntuple with: " << SCEP_EnergyDepF  << " MeV in SCEPCal front and " << SCEP_EnergyDepR  << " MeV in SCEPCal rear" << std::endl;
  // fill ntuple event by event
  analysisManager->FillNtupleDColumn(0, Energyem);
  analysisManager->FillNtupleDColumn(1, EnergyScin);
  analysisManager->FillNtupleDColumn(2, EnergyCher);
  analysisManager->FillNtupleDColumn(3, NofCherenkovDetected);
  analysisManager->FillNtupleDColumn(4, EnergyTot);
  analysisManager->FillNtupleDColumn(5, PrimaryParticleEnergy);
  analysisManager->FillNtupleSColumn(6, PrimaryParticleName);
  analysisManager->FillNtupleDColumn(7, neutrinoleakage);
  analysisManager->FillNtupleDColumn(8, leakage);
  
  analysisManager->FillNtupleDColumn(9, SCEP_EnergyDepF);
  analysisManager->FillNtupleDColumn(10,SCEP_NCherProdF);
  analysisManager->FillNtupleDColumn(11,SCEP_EnergyDepR);
  analysisManager->FillNtupleDColumn(12,SCEP_NCherProdR);
  
  analysisManager->AddNtupleRow();//columns with vector are automatically filled with this function
  
  

  //print here if you need event by event some information of the screen
  //G4cout<<EnergyTot<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
