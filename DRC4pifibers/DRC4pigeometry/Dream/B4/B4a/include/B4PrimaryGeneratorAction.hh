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
// $Id: B4PrimaryGeneratorAction.hh 95508 2016-02-12 13:52:06Z gcosmo $
// 
/// \file B4PrimaryGeneratorAction.hh
/// \brief Definition of the B4PrimaryGeneratorAction class

#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <map>

class G4GeneralParticleSource;
//class G4ParticleGun;
class G4Event;

class G4VPrimaryGenerator;
class H02PrimaryGeneratorMessenger;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class B4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
//   B4PrimaryGeneratorAction(char*);    
  B4PrimaryGeneratorAction();    
  virtual ~B4PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);
  
  // set methods
  void SetRandomFlag(G4bool value);

  
  void SetGenerator(G4VPrimaryGenerator* gen);
  void SetGenerator(G4String genname);

  G4VPrimaryGenerator* GetGenerator() const;
  G4String GetGeneratorName() const;

private:
   G4GeneralParticleSource* fGeneralParticleSource;
  //G4ParticleGun*  fParticleGun; // G4 particle gun

  G4VPrimaryGenerator* particleGun;
  G4VPrimaryGenerator* hepmcAscii;
  G4VPrimaryGenerator* pythiaGen;

  G4VPrimaryGenerator* currentGenerator;
  G4String currentGeneratorName;
  std::map<G4String, G4VPrimaryGenerator*> gentypeMap;

  H02PrimaryGeneratorMessenger* messenger;
//   char* m_outputFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void B4PrimaryGeneratorAction::SetGenerator(G4VPrimaryGenerator* gen)
{
  currentGenerator= gen;
}

inline void B4PrimaryGeneratorAction::SetGenerator(G4String genname)
{
  std::map<G4String, G4VPrimaryGenerator*>::iterator
       pos = gentypeMap.find(genname);
  if(pos != gentypeMap.end()) {
    currentGenerator= pos->second;
    currentGeneratorName= genname;
  }
}

inline G4VPrimaryGenerator* B4PrimaryGeneratorAction::GetGenerator() const
{
  return currentGenerator;
}

inline G4String B4PrimaryGeneratorAction::GetGeneratorName() const
{
  return currentGeneratorName;
}

#endif
