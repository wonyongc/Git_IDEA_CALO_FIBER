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
// $Id: B4PrimaryGeneratorAction.cc 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"

#include "H02PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fGeneralParticleSource(0)
   //fParticleGun(0)
{
  //G4int nofParticles = 1;
  fGeneralParticleSource = new G4GeneralParticleSource();
  //fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  G4ParticleDefinition* particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fGeneralParticleSource->SetParticleDefinition(particleDefinition);
  //fGeneralParticleSource->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //fGeneralParticleSource->SetParticleEnergy(50.*MeV);
  //fParticleGun->SetParticleDefinition(particleDefinition);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //fParticleGun->SetParticleEnergy(50.*MeV);
}
*/
B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
{
   // default generator is particle gun.
  currentGenerator= particleGun= new G4GeneralParticleSource();
  currentGeneratorName= "gps";
  hepmcAscii= new HepMCG4AsciiReader();
#ifdef G4LIB_USE_PYTHIA
  pythiaGen= new HepMCG4PythiaInterface();
#else
  pythiaGen= 0;
#endif

  gentypeMap["gps"]= particleGun;
  gentypeMap["hepmcAscii"]= hepmcAscii;
  gentypeMap["pythia"]= pythiaGen;

  messenger= new H02PrimaryGeneratorMessenger(this); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fGeneralParticleSource;
  //delete fParticleGun;
}
*/
B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{ 
  delete messenger;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0;
  G4LogicalVolume* worlLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = 0;
  if ( worlLV) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid()); 
  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
  
  // Set gun position
  //fGeneralParticleSource
  //  ->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));
  //fParticleGun
  //  ->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));

  fGeneralParticleSource->GeneratePrimaryVertex(anEvent);
  //fParticleGun->GeneratePrimaryVertex(anEvent);
}
*/
void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(currentGenerator)
    currentGenerator-> GeneratePrimaryVertex(anEvent);
  else
    G4Exception("H02PrimaryGeneratorAction::GeneratePrimaries",
                "InvalidSetup", FatalException,
                "Generator is not instanciated.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

