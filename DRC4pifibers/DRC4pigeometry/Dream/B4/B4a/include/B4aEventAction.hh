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
// $Id: B4aEventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4aEventAction.hh
/// \brief Definition of the B4aEventAction class

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>
#include <map>
#include "G4ThreeVector.hh"


/// Event action class 

class B4aEventAction : public G4UserEventAction
{
  public:
    B4aEventAction();
    virtual ~B4aEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    
    void Addneutrinoleakage(G4double de); //add energy of neutrinos in the ball containing the calorimeter
    void Addleakage(G4double de); //add energy of all particles that are not neutrinos (or anti_neutrinos) in the ball containing the calorimeter
    void Addem(G4double de);  //Add em component
    void AddScin(G4double de);//Add energy in scintillating fibers
    void AddCher(G4double de);//Add energy in Cherenkov fibers
    void AddCherenkov();//Add cherenkov photoelectron
    
    void AddScepEneF(G4double de);   //Add hit for given SCEPCal crystals
    void AddScepEneR(G4double de);   //Add hit for given SCEPCal crystals
    void AddScepCherF();             //Add hit for given SCEPCal crystals
    void AddScepCherR();             //Add hit for given SCEPCal crystals
    
    void AddScepTimingEneTime(G4double de, G4double time_hit, G4double layerID);   //Add hit for given SCEPCal Timing crystals
    
    
    
    //void AddScintillation();
    void Addenergy(G4double de);//Add all energy deposited
    //void AddEnergyfibre(G4double de, G4int number);//Add energy in copy number fiber
    //void AddSignalfibre(G4int number);
    void SavePrimaryParticle(G4String name);
    void SavePrimaryEnergy(G4double primaryparticleenergy);
    void SavePrimaryMomentum(G4double primaryparticle_px, G4double primaryparticle_py, G4double primaryparticle_pz);

    //to save vectors in ntuple
    std::vector<G4double>& GetVectorSignalsR() {return VectorSignalsR;}
    std::vector<G4double>& GetVectorSignalsL() {return VectorSignalsL;} 
    std::vector<G4double>& GetVectorSignalsCherR() {return VectorSignalsCherR;}
    std::vector<G4double>& GetVectorSignalsCherL() {return VectorSignalsCherL;}
    std::vector<G4double>& GetVectorR() {return VectorR;}
    std::vector<G4double>& GetVectorL() {return VectorL;}
    std::vector<G4double>& GetVectorR_loop() {return VectorR_loop;}
    std::vector<G4double>& GetVectorL_loop() {return VectorL_loop;}
    
    std::vector<G4double>& GetPrimaryParticleMomentum()    {return PrimaryParticleMomentum;}
    std::vector<G4double>& GetVecScep_CrystalID()   {return VecHit_CrystalID;}
    std::vector<G4double>& GetVecScep_ScepEneDepF() {return VecHit_ScepEneDepF;}
    std::vector<G4double>& GetVecScep_ScepEneDepR() {return VecHit_ScepEneDepR;}
    std::vector<G4double>& GetVecScep_ScepCherF()   {return VecHit_ScepCherF;}
    std::vector<G4double>& GetVecScep_ScepCherR()   {return VecHit_ScepCherR;}
    
    std::vector<G4double>& GetVecScep_Timing_CrystalID_F() {return VecHit_Timing_CrystalID_F;}
    std::vector<G4double>& GetVecScep_Timing_CrystalID_R() {return VecHit_Timing_CrystalID_R;}
    std::vector<G4double>& GetVecScep_Timing_ScepEneDepF() {return VecHit_Timing_ScepEneDepF;}
    std::vector<G4double>& GetVecScep_Timing_ScepEneDepR() {return VecHit_Timing_ScepEneDepR;}
    std::vector<G4double>& GetVecScep_Timing_ScepTimeF()   {return VecHit_Timing_ScepTimeF;}
    std::vector<G4double>& GetVecScep_Timing_ScepTimeR()   {return VecHit_Timing_ScepTimeR;}
    
    
    

    //to fill vectors
    void AddVectorScinEnergyR(G4double de, G4int tower, G4int slice);   //fill vector of scintillating fibers with energy deposition
    void AddVectorScinEnergyL(G4double de, G4int tower, G4int slice);   //fill vector left side
    void AddVectorCherPER(G4int c_signal, G4int tower, G4int slice);    //fill vector of cherenkov fibers with chernekov photoelectrons
    void AddVectorCherPEL(G4int c_signal, G4int tower, G4int slice);
    void AddVectorR(G4double de, G4int tower, G4int slice);
    void AddVectorL(G4double de, G4int tower, G4int slice);
    void AddVectorR_loop(G4double de, G4int tower, G4int slice);
    void AddVectorL_loop(G4double de, G4int tower, G4int slice);
    
    void AddScepHit(G4double chId, G4double de, G4String hitType);
    void AddScepTimingHit(G4double chId, G4double de, G4double time_hit, G4double layerID);
    
    
    typedef struct FiberInfo {
        G4double F_ID, F_E, F_X, F_Y, F_Z; //fiber saturated energy
        G4int F_Type, F_slice, F_tower; //C==0 S ==1;
    } Fiber_Info;
    
    void WriteFiber_Info(G4double FID, G4double FE, G4int FType, G4ThreeVector Fpos, G4int slice, G4int tower);
    
    typedef struct TrackingInfo {
        G4double T_ID, T_X, T_Y, T_Z, T_Ek;
        G4String T_Name;
    } Tracking_Info;
    
    void WriteTracking_Info(G4double T_ID, G4ThreeVector Tpos, G4String Name, G4double Ek);
    
  private:
    G4int cont;
    Fiber_Info Fiber_Hits[1000000];
    Tracking_Info Tracking_Hits[200];
    G4double  Energyem; //Energy of em component
    G4double  EnergyScin; //Energy in scintillating fibers
    G4double  EnergyCher; //Energy in Cherenkov fibers
    G4int     NofCherenkovDetected; //Number of Cherenkov photons detected (in cherenkov fibers)
    //G4int     NofScintillationDetected;//Number of Scintillating photons detected (in scintillating fibers)
    G4double  EnergyTot;//Total energy deposited (does not count invisibile energy)
    //G4double  Signalfibre[64];//Signal in 64 single module fibers, to be used with AddEnergyfibre
    G4String PrimaryParticleName; //Name of primary particle
    G4double PrimaryParticleEnergy;//Primary particle energy
    std::vector<G4double>  PrimaryParticleMomentum;//Primary particle energy
    G4double neutrinoleakage; //leakage neutrino
    G4double leakage; //leakage non neutrino
    
    G4double  SCEP_EnergyDepF; //total energy deposited in the crystal SCEPCal volume due to ionization
    G4int     SCEP_NCherProdF; //total number of cherenkov photons produced in the SCEPCal volume
    G4double  SCEP_EnergyDepR; //total energy deposited in the crystal SCEPCal volume due to ionization
    G4int     SCEP_NCherProdR; //total number of cherenkov photons produced in the SCEPCal volume
    
    G4double  SCEP_Timing_EnergyDepF; //total energy deposited in the crystal SCEPCal timing volume due to ionization
    G4double  SCEP_Timing_EnergyDepR; //total energy deposited in the crystal SCEPCal timing volume due to ionization
    G4double  SCEP_Timing_TimeF; //time stamp of the hit in the crystal SCEPCal timing volume due to ionization
    G4double  SCEP_Timing_TimeR; //time stamp of the hit in the crystal SCEPCal timing volume due to ionization
    
    std::vector<G4double> VecHit_CrystalID;   //for each crystal with hits --> total number of cherenkov photons produced in the SCEPCal volume
    std::vector<G4double> VecHit_ScepEneDepF; //for each crystal with hits --> total energy deposited in the crystal SCEPCal volume due to ionization
    std::vector<G4double> VecHit_ScepEneDepR; //for each crystal with hits --> total number of cherenkov photons produced in the SCEPCal volume
    std::vector<G4double> VecHit_ScepCherF;   //for each crystal with hits --> total energy deposited in the crystal SCEPCal volume due to ionization
    std::vector<G4double> VecHit_ScepCherR;   //for each crystal with hits --> total number of cherenkov photons produced in the SCEPCal volume
    
    std::vector<G4double> VecHit_Timing_CrystalID_F;   //for each crystal with hits --> total number of cherenkov photons produced in the SCEPCal volume
    std::vector<G4double> VecHit_Timing_CrystalID_R;   //for each crystal with hits --> total number of cherenkov photons produced in the SCEPCal volume
    std::vector<G4double> VecHit_Timing_ScepEneDepF; //for each crystal with hits --> total energy deposited in the crystal SCEPCal volume due to ionization
    std::vector<G4double> VecHit_Timing_ScepEneDepR; //for each crystal with hits --> total number of cherenkov photons produced in the SCEPCal volume
    std::vector<G4double> VecHit_Timing_ScepTimeF;   //for each crystal with hits --> total energy deposited in the crystal SCEPCal volume due to ionization
    std::vector<G4double> VecHit_Timing_ScepTimeR;   //for each crystal with hits --> total number of cherenkov photons produced in the SCEPCal volume
    
    std::vector<G4double> VectorR_loop;
    std::vector<G4double> VectorL_loop;

    std::vector<G4double> VectorSignalsR;//Vector filled with scintillating fibers energy deposits
    std::vector<G4double> VectorSignalsL;//vector filled for left side
    std::vector<G4double> VectorSignalsCherR;//Vector filled with Cherenkov fibers Cherenkov photoelectrons
    std::vector<G4double> VectorSignalsCherL;//vector filled for left side
    
    std::vector<G4double> VectorR; //vector with energy deposited in towers right
    std::vector<G4double> VectorL;
};

// inline functions
inline void B4aEventAction::Addneutrinoleakage(G4double de){
    neutrinoleakage += de;
}

inline void B4aEventAction::Addleakage(G4double de){
    leakage += de;
}

inline void B4aEventAction::AddVectorR(G4double de, G4int tower, G4int slice){
	VectorR.at(tower+(slice*75)) += de;	
}

inline void B4aEventAction::AddVectorL(G4double de, G4int tower, G4int slice){
	tower = -1*tower;
	VectorL.at(tower+(slice*75)) += de;
}

inline void B4aEventAction::AddVectorR_loop(G4double de, G4int tower, G4int slice){
    VectorR_loop.at(tower+(slice*75)) = de; 
}

inline void B4aEventAction::AddVectorL_loop(G4double de, G4int tower, G4int slice){
    tower = -1*tower;
    VectorL_loop.at(tower+(slice*75)) = de;
}

inline void B4aEventAction::WriteFiber_Info(G4double FID, G4double FE, G4int FType, G4ThreeVector Fpos, G4int slice, G4int tower){
    int k=0;
    while (Fiber_Hits[k].F_ID!=0 && Fiber_Hits[k].F_ID!=FID){
        k++;}
    Fiber_Hits[k].F_ID = FID;
    Fiber_Hits[k].F_E += FE;
    Fiber_Hits[k].F_Type = FType;
    Fiber_Hits[k].F_X = Fpos[0];
    Fiber_Hits[k].F_Y = Fpos[1];
    Fiber_Hits[k].F_Z = Fpos[2];
    Fiber_Hits[k].F_slice = slice;
    Fiber_Hits[k].F_tower = tower;
}

inline void B4aEventAction::WriteTracking_Info(G4double TID, G4ThreeVector Tpos, G4String Name, G4double Ek){
    int k=0;
    while (Tracking_Hits[k].T_ID!=0 && Tracking_Hits[k].T_ID!=TID){
        k++;}
    Tracking_Hits[k].T_ID = TID;
    Tracking_Hits[k].T_X = Tpos[0];
    Tracking_Hits[k].T_Y = Tpos[1];
    Tracking_Hits[k].T_Z = Tpos[2];
    Tracking_Hits[k].T_Name = Name;
    if (Tracking_Hits[k].T_Ek<=0.) {Tracking_Hits[k].T_Ek = Ek;}
}


inline void B4aEventAction::SavePrimaryParticle(G4String name){
  PrimaryParticleName = name;
}

inline void B4aEventAction::SavePrimaryEnergy(G4double primaryparticleenergy){
  PrimaryParticleEnergy = primaryparticleenergy;
}

inline void B4aEventAction::SavePrimaryMomentum(G4double px, G4double py, G4double pz){
//   std::cout << " momentum = (" << px << ", " << py << ", " << pz << ")" << std::endl;
  PrimaryParticleMomentum.push_back(px);
  PrimaryParticleMomentum.push_back(py);
  PrimaryParticleMomentum.push_back(pz);
}

inline void B4aEventAction::AddVectorScinEnergyR(G4double de, G4int tower, G4int slice) {
    VectorSignalsR.at(tower+(slice*75)) += de;
}

inline void B4aEventAction::AddVectorScinEnergyL(G4double de, G4int tower, G4int slice) {
    tower = -1*tower;
    VectorSignalsL.at(tower+(slice*75)) += de;
}

inline void B4aEventAction::AddVectorCherPEL(G4int c_signal, G4int tower, G4int slice) {
	tower = -1*tower;
    VectorSignalsCherL.at(tower+(slice*75)) = VectorSignalsCherL.at(tower+(slice*75))+c_signal;
}

inline void B4aEventAction::AddVectorCherPER(G4int c_signal, G4int tower, G4int slice) {
    VectorSignalsCherR.at(tower+(slice*75)) = VectorSignalsCherR.at(tower+(slice*75))+c_signal;
}

inline void B4aEventAction::Addem(G4double de) {
  Energyem += de; 
}

inline void B4aEventAction::AddScin(G4double de){
  EnergyScin += de;
}

inline void B4aEventAction::AddCher(G4double de){
  EnergyCher += de;
}

inline void B4aEventAction::AddCherenkov(){
  NofCherenkovDetected = NofCherenkovDetected + 1;
}


inline void B4aEventAction::AddScepEneF(G4double de)   {SCEP_EnergyDepF += de;}
inline void B4aEventAction::AddScepCherF()             {SCEP_NCherProdF = SCEP_NCherProdF + 1;}
inline void B4aEventAction::AddScepEneR(G4double de)   {SCEP_EnergyDepR += de;}
inline void B4aEventAction::AddScepCherR()             {SCEP_NCherProdR = SCEP_NCherProdR + 1;}


inline void B4aEventAction::AddScepTimingEneTime(G4double de, G4double time_hit, G4double layerID)
{ 
    if (layerID == 1) 
    {
        SCEP_Timing_EnergyDepF += de;
        SCEP_Timing_TimeF += de*time_hit;
    }
    else //if (layerID == 2) 
    {
        SCEP_Timing_EnergyDepR += de;
        SCEP_Timing_TimeR += de*time_hit;
    }
}



inline void B4aEventAction::AddScepHit(G4double chId, G4double de, G4String hitType) 
{
    bool found = false;
    for (long unsigned int it = 0; it<VecHit_CrystalID.size(); it++)
    {
    
        if (VecHit_CrystalID.at(it) == chId)
        {
            found = true;
        
            if      (hitType == "FrontEne")  VecHit_ScepEneDepF.at(it) += de;
            else if (hitType == "RearEne")   VecHit_ScepEneDepR.at(it) += de;
            else if (hitType == "FrontCher") VecHit_ScepCherF.at(it)   += de;
            else if (hitType == "RearCher")  VecHit_ScepCherR.at(it)   += de;
            break;
        }
    }    
    if (!found)  //first we get a hit in this detector location (either front or rear) --> add de to hit segment and set entry initialized to 0 also on all the other vectors
    {
        VecHit_CrystalID.push_back(chId);    
        
        if      (hitType == "FrontEne")  
        {
            VecHit_ScepEneDepF.push_back(de);    
            VecHit_ScepEneDepR.push_back(0.);    
            VecHit_ScepCherF.push_back(0.);    
            VecHit_ScepCherR.push_back(0.);    
        }
        else if (hitType == "RearEne")   
        {
            VecHit_ScepEneDepF.push_back(0.);    
            VecHit_ScepEneDepR.push_back(de);    
            VecHit_ScepCherF.push_back(0.);    
            VecHit_ScepCherR.push_back(0.);    
        }
        else if (hitType == "FrontCher") 
        {
            VecHit_ScepEneDepF.push_back(0.);    
            VecHit_ScepEneDepR.push_back(0.);    
            VecHit_ScepCherF.push_back(de);    
            VecHit_ScepCherR.push_back(0.);    
        }
        else if (hitType == "RearCher")  
        {
            VecHit_ScepEneDepF.push_back(0.);    
            VecHit_ScepEneDepR.push_back(0.);    
            VecHit_ScepCherF.push_back(0.);    
            VecHit_ScepCherR.push_back(de);    
        }
    }    
//     std::cout << "added hit @ " << chId << " with Ene = " << de << std::endl;
}


inline void B4aEventAction::AddScepTimingHit(G4double chId, G4double de, G4double time_hit, G4double layerID) 
{    
//     std::cout << "chId  = " << chId << " :: de = " << de << " :: time_hit = " << time_hit << " :: layerID = " << layerID << std::endl;    
    bool found = false;
    if (layerID == 1)
    {
        for (long unsigned int it = 0; it<VecHit_Timing_CrystalID_F.size(); it++)
        {
            if (VecHit_Timing_CrystalID_F.at(it) == chId)
            {
                found = true;                    
                VecHit_Timing_ScepEneDepF.at(it) += de;
                VecHit_Timing_ScepTimeF.at(it) += time_hit*de;                
                break;
            }
        }    
        if (!found)
        {
            VecHit_Timing_CrystalID_F.push_back(chId);    
            VecHit_Timing_ScepEneDepF.push_back(de);
            VecHit_Timing_ScepTimeF.push_back(time_hit*de);            
        }
    }
    else
    {
        for (long unsigned int it = 0; it<VecHit_Timing_CrystalID_R.size(); it++)
        {
            if (VecHit_Timing_CrystalID_R.at(it) == chId)
            {
                found = true;                    
                VecHit_Timing_ScepEneDepR.at(it) += de;
                VecHit_Timing_ScepTimeR.at(it) += time_hit*de;                
                break;
            }
        }    
        if (!found)
        {
            VecHit_Timing_CrystalID_R.push_back(chId);    
            VecHit_Timing_ScepEneDepR.push_back(de);
            VecHit_Timing_ScepTimeR.push_back(time_hit*de);            
        }        
    }
}

/*inline void B4aEventAction::AddScintillation(){
  NofScintillationDetected = NofScintillationDetected +1;
}*/

inline void B4aEventAction::Addenergy(G4double de){
  EnergyTot += de;
}

/*inline void B4aEventAction::AddEnergyfibre(G4double de, G4int number){
    Signalfibre[number] += de;
}*/

/*inline void B4aEventAction::AddSignalfibre(G4int number){
    Signalfibre[number] = Signalfibre[number] + 1;
}*/
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
